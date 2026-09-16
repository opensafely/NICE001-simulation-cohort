/* Purpose: Estimate entropy-balancing weights that align women in the OpenSAFELY
DXA-eligible cohort with women treatment-eligible cohort.*/

version 16.1
clear all
set more off

/* --------------------------- User settings ---------------------------- */
capture adopath ++ "analysis/extra_ados"
capture which ebalance
capture which stddiffi

* Solver controls. 
local ebal_maxiter 500
local ebal_tolerance 0.000001
local verify_tolerance 0.00001

use "output/dxa_eligible_population.dta", clear
keep if dxa_eligible == 1
keep if lower(strtrim(sex)) == "female"
local source_n = r(N)

/* ---------------- aggregate targets ---------------------------------- */
local target_total 16592

* Age bands: 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84, 85-89, 90+.
matrix N_age      = (674, 1080, 1984, 2232, 2691, 2637, 2436, 1762, 1096)
matrix C_lowbmi   = ( 14,   20,   31,   32,   39,   21,   29,   34,   32)
matrix C_prevfx   = (622,  944, 1472, 1379, 1555, 1458, 1353,  974,  606)
matrix C_smoke    = (197,  237,  376,  335,  352,  281,  199,  126,   63)
matrix C_alcohol  = (105,  146,  248,  259,  287,  268,  233,  165,  103)
matrix C_parenthf = ( 96,  109,  158,  163,  157,  195,  167,  110,   95)
matrix C_ra       = ( 70,  122,  167,  178,  193,  172,  159,  110,   60)
matrix C_steroid  = (142,  203,  462,  466,  426,  423,  365,  250,  160)

matrix M_bmi = ( ///
    26.95727002967358, 27.32500000000003, 27.34077620967744, ///
    27.24439964157706, 27.62772203641772, 28.23147516116799, ///
    27.94831691297208, 27.07474460839949, 26.37691605839416)
	
matrix SD_bmi = ( ///
    5.134356696583268, 5.330737368292470, 4.950353176611750, ///
    4.847662751582217, 5.145261361133998, 5.269837424449750, ///
    5.275743691528498, 4.965330463029384, 4.560051467333120)

scalar target_bmi_mean_all = 27.49658269045323
scalar target_bmi_sd_all   =  5.102052591344626

* Verify that the age counts reproduce stated total.
scalar __target_n_check = 0
forvalues g = 1/9 {
    scalar __target_n_check = scalar(__target_n_check) + N_age[1, `g']
}
assert scalar(__target_n_check) == `target_total'

* Convert target counts into within-age proportions.
foreach characteristic in lowbmi prevfx smoke alcohol parenthf ra steroid {
    matrix P_`characteristic' = J(1, 9, .)
    forvalues g = 1/9 {
        matrix P_`characteristic'[1, `g'] = C_`characteristic'[1, `g'] / N_age[1, `g']
    }
}

/* ---------------- Harmonised OpenSAFELY variables --------------------- */
generate byte cal_ageband = .
replace cal_ageband = 1 if inrange(age, 50, 54)
replace cal_ageband = 2 if inrange(age, 55, 59)
replace cal_ageband = 3 if inrange(age, 60, 64)
replace cal_ageband = 4 if inrange(age, 65, 69)
replace cal_ageband = 5 if inrange(age, 70, 74)
replace cal_ageband = 6 if inrange(age, 75, 79)
replace cal_ageband = 7 if inrange(age, 80, 84)
replace cal_ageband = 8 if inrange(age, 85, 89)
replace cal_ageband = 9 if inrange(age, 90, 99)

assert !missing(cal_ageband)
label define cal_ageband_label 1 "50-54" 2 "55-59" 3 "60-64" 4 "65-69" 5 "70-74" 6 "75-79" 7 "80-84" 8 "85-89" 9 "90+"
label values cal_ageband cal_ageband_label

generate double cal_bmi = bmi_raw if inrange(bmi_raw, 10, 80)
generate byte cal_lowbmi = cal_bmi <= 19 if !missing(cal_bmi)

egen byte __previous_fracture_max = rowmax(dx_hip_fracture_n dx_vertebral_fracture_n dx_wrist_fracture_n dx_proximal_humerus_fracture_n)
generate byte cal_prevfx = __previous_fracture_max > 0 if !missing(__previous_fracture_max)
drop __previous_fracture_max

generate byte cal_smoke = inrange(qf_smoke_cat, 2, 4) if !missing(qf_smoke_cat)
generate byte cal_alcohol = inrange(qf_alcohol_cat6, 3, 5) if !missing(qf_alcohol_cat6)
generate byte cal_parenthf = fh_parental_hip_fracture if !missing(fh_parental_hip_fracture)
generate byte cal_ra = b_ra if !missing(b_ra)
generate byte cal_steroid = b_corticosteroids if !missing(b_corticosteroids)

label variable cal_bmi      "Uncapped BMI used for calibration"
label variable cal_lowbmi   "BMI <=19 kg/m2"
label variable cal_prevfx   "Any previous hip/vertebral/wrist/proximal-humerus fracture"
label variable cal_smoke    "Current smoking"
label variable cal_alcohol  "Alcohol >=3 units/day"
label variable cal_parenthf "Parental hip fracture"
label variable cal_ra       "Rheumatoid arthritis"
label variable cal_steroid  ">=2 oral corticosteroid prescriptions in prior 6 months"

local binary_variables ///
    cal_lowbmi cal_prevfx cal_steroid cal_ra ///
    cal_smoke cal_alcohol cal_parenthf
local binary_target_matrices P_lowbmi P_prevfx P_steroid P_ra P_smoke P_alcohol P_parenthf

foreach variable of local binary_variables {
    assert inlist(`variable', 0, 1, .)
}

/* ---------------- Missingness report --------------------------------- */
tempname missing_post
tempfile missing_results
postfile `missing_post' ///
    str40 characteristic str12 age_band ///
    double n_total n_observed n_missing pct_missing ///
    using "`missing_results'", replace

local missing_variables ///
    cal_bmi cal_lowbmi cal_prevfx cal_smoke ///
    cal_alcohol cal_parenthf cal_ra cal_steroid

foreach variable of local missing_variables {
    local characteristic : variable label `variable'

    quietly count
    local ntotal = r(N)
    quietly count if !missing(`variable')
    local nobserved = r(N)
    local nmissing = `ntotal' - `nobserved'
    local pmissing = 100 * `nmissing' / `ntotal'
    post `missing_post' ("`characteristic'") ("Overall") (`ntotal') (`nobserved') (`nmissing') (`pmissing')

    forvalues g = 1/9 {
        local group_label : label cal_ageband_label `g'
        quietly count if cal_ageband == `g'
        local ntotal = r(N)
        quietly count if cal_ageband == `g' & !missing(`variable')
        local nobserved = r(N)
        local nmissing = `ntotal' - `nobserved'
        local pmissing = 100 * `nmissing' / `ntotal'
        post `missing_post' ("`characteristic'") ("`group_label'") (`ntotal') (`nobserved') (`nmissing') (`pmissing')
    }
}
postclose `missing_post'

preserve
use "`missing_results'", clear
sort characteristic age_band
save "output/women_calibration_missingness.dta", replace
export delimited using "output/women_calibration_missingness.csv", replace
restore

/* ---------------- Feasibility and overlap checks --------------------- */
forvalues g = 1/9 {
    quietly count if cal_ageband == `g'
    if r(N) == 0 {
        local group_label : label cal_ageband_label `g'
        display as error "No OpenSAFELY women in age band `group_label'."
        exit 459
    }
}

local n_binary : word count `binary_variables'
forvalues j = 1/`n_binary' {
    local variable : word `j' of `binary_variables'
    local target_matrix : word `j' of `binary_target_matrices'

    forvalues g = 1/9 {
        scalar __target_p = `target_matrix'[1, `g']
        assert scalar(__target_p) > 0 & scalar(__target_p) < 1

        quietly count if cal_ageband == `g' & `variable' == 0
        local n0 = r(N)
        quietly count if cal_ageband == `g' & `variable' == 1
        local n1 = r(N)

        if `n0' == 0 | `n1' == 0 {
            local group_label : label cal_ageband_label `g'
            display as error ///
                "Insufficient overlap: `variable' has no observed 0 or 1 in age band `group_label'."
            exit 459
        }
    }
}

forvalues g = 1/9 {
    quietly summarize cal_bmi if cal_ageband == `g'
    if r(N) < 2 | r(sd) == 0 | missing(r(sd)) {
        local group_label : label cal_ageband_label `g'
        display as error "Insufficient BMI variation in age band `group_label'."
        exit 459
    }
}

/* ---------------- Construct centred calibration moments -------------- */
local f_age ""
forvalues g = 1/9 {
    generate double eb_a`g' = (cal_ageband == `g')
    * ebalance separately normalises the weights to sum to the sample size.
    * Therefore only eight age indicators are supplied to the solver. The
    * ninth age proportion follows automatically because all nine sum to one.
    if `g' <= 8 {
        local f_age "`f_age' eb_a`g'"
    }
}

* The first four binary variables form the core stage. The final three are
* added in the full shared-risk-factor stage.
local f_core_binary ""
local f_full_binary ""

forvalues j = 1/`n_binary' {
    local variable : word `j' of `binary_variables'
    local target_matrix : word `j' of `binary_target_matrices'

    forvalues g = 1/9 {
        scalar __target_p = `target_matrix'[1, `g']
        generate double eb_b`j'_`g' = 0
        replace eb_b`j'_`g' = `variable' - scalar(__target_p) if cal_ageband == `g' & !missing(`variable')
        local f_full_binary "`f_full_binary' eb_b`j'_`g'"
        if `j' <= 4 {
            local f_core_binary "`f_core_binary' eb_b`j'_`g'"
        }
    }
}

local f_bmi ""
forvalues g = 1/9 {
    scalar __target_mean = M_bmi[1, `g']
    * Excel PivotTable StdDev is a sample SD. Convert it to the corresponding
    * population central moment before imposing the weighted moment constraint.
    scalar __target_var_pop = ((N_age[1, `g'] - 1) / N_age[1, `g']) * SD_bmi[1, `g']^2
    generate double eb_m`g' = 0
    replace eb_m`g' = (cal_bmi - scalar(__target_mean)) / 10 if cal_ageband == `g' & !missing(cal_bmi)
    generate double eb_v`g' = 0
    replace eb_v`g' = ((cal_bmi - scalar(__target_mean))^2 - scalar(__target_var_pop)) / 100 if cal_ageband == `g' & !missing(cal_bmi)
    local f_bmi "`f_bmi' eb_m`g' eb_v`g'"
}

local f_core     "`f_age' `f_core_binary'"
local f_full     "`f_age' `f_full_binary'"
local f_full_bmi "`f_age' `f_full_binary' `f_bmi'"

* manualtargets() expects one target mean for every supplied feature, in the
* same order as the feature list. Age indicators target FRAX's proportions;
* every centred binary/BMI calibration feature has a target mean of zero.
local targets_age ""
forvalues g = 1/8 {
    scalar __target_age_p = N_age[1, `g'] / `target_total'
    local target_value : display %21.15g scalar(__target_age_p)
    local targets_age "`targets_age' `target_value'"
}

local targets_core "`targets_age'"
local n_core_binary : word count `f_core_binary'
forvalues j = 1/`n_core_binary' {
    local targets_core "`targets_core' 0"
}

local targets_full "`targets_age'"
local n_full_binary : word count `f_full_binary'
forvalues j = 1/`n_full_binary' {
    local targets_full "`targets_full' 0"
}

local targets_full_bmi "`targets_full'"
local n_bmi_features : word count `f_bmi'
forvalues j = 1/`n_bmi_features' {
    local targets_full_bmi "`targets_full_bmi' 0"
}

foreach stage in age core full full_bmi {
    local features "`f_`stage''"
    local targets "`targets_`stage''"
    local n_features : word count `features'
    local n_targets : word count `targets'

    if `n_features' != `n_targets' {
        display as error "Internal error in `stage': `n_features' features but `n_targets' targets."
        exit 198
    }
}

/* ---------------- Entropy balancing with ebalance -------------------- */
foreach stage in age core full full_bmi {
    display as text "Estimating `stage' calibration weights..."
    local features "`f_`stage''"
    local targets "`targets_`stage''"
    local n_features : word count `features'

    capture drop cw_`stage'
    capture noisily ebalance `features', ///
        manualtargets(`targets') ///
        generate(cw_`stage') ///
        maxiter(`ebal_maxiter') ///
        tolerance(`ebal_tolerance')

    local ebal_rc = _rc
    if `ebal_rc' {
        display as error "ebalance failed at the `stage' stage. Review the output above and the overlap checks."
        exit `ebal_rc'
    }

    confirm variable cw_`stage'
    assert !missing(cw_`stage')
    assert cw_`stage' > 0

    * Normalise each weight set to mean one. This does not alter any weighted
    * means and makes the stage-specific weights directly comparable.
    quietly summarize cw_`stage', meanonly
    replace cw_`stage' = cw_`stage' / r(mean)

    * Independently verify every target rather than relying only on the
    * optimiser's printed convergence message.
    scalar __max_target_deviation = 0
    forvalues j = 1/`n_features' {
        local feature : word `j' of `features'
        local target : word `j' of `targets'
        quietly summarize `feature' [aw=cw_`stage'], meanonly
        scalar __target_deviation = abs(r(mean) - (`target'))
        if scalar(__target_deviation) > scalar(__max_target_deviation) {
            scalar __max_target_deviation = scalar(__target_deviation)
        }
    }

    if scalar(__max_target_deviation) > `verify_tolerance' {
        display as error "The `stage' stage did not meet the verification tolerance."
        display as error "Maximum target deviation = " scalar(__max_target_deviation)
        exit 430
    }

    quietly summarize cw_`stage', meanonly
    matrix R_`stage' = (1, scalar(__max_target_deviation), r(min), r(max))
    matrix colnames R_`stage' = converged max_target_deviation min_weight max_weight
}

generate double calibration_weight = cw_full_bmi
label variable calibration_weight "Women calibration weight: age + shared risks + BMI moments"

* Sensitivity weight truncated at the final weight's 99th percentile.
quietly summarize calibration_weight, detail
scalar __weight_p99 = r(p99)
generate double calibration_weight_p99 = min(calibration_weight, scalar(__weight_p99))
quietly summarize calibration_weight_p99, meanonly
replace calibration_weight_p99 = calibration_weight_p99 / r(mean)
label variable calibration_weight_p99 "Calibration weight truncated at p99 and renormalised"

/* ---------------- Weight diagnostics --------------------------------- */
tempname weight_post
tempfile weight_results
postfile `weight_post' ///
    str16 stage double converged max_target_deviation ///
    mean sd min p1 p50 p99 max effective_sample_size ///
    using "`weight_results'", replace

foreach stage in age core full full_bmi {
    quietly summarize cw_`stage', detail
    local wmean = r(mean)
    local wsd = r(sd)
    local wmin = r(min)
    local wp1 = r(p1)
    local wp50 = r(p50)
    local wp99 = r(p99)
    local wmax = r(max)

    tempvar weight_squared
    generate double `weight_squared' = cw_`stage'^2
    quietly summarize cw_`stage', meanonly
    local sum_weight = r(sum)
    quietly summarize `weight_squared', meanonly
    local sum_weight_squared = r(sum)
    local ess = (`sum_weight'^2) / `sum_weight_squared'
    drop `weight_squared'

    local result_matrix R_`stage'
    local converged = `result_matrix'[1, 1]
    local max_error = `result_matrix'[1, 2]

    post `weight_post' ("`stage'") ///
        (`converged') (`max_error') ///
        (`wmean') (`wsd') (`wmin') (`wp1') (`wp50') ///
        (`wp99') (`wmax') (`ess')
}

quietly summarize calibration_weight_p99, detail
local wmean = r(mean)
local wsd = r(sd)
local wmin = r(min)
local wp1 = r(p1)
local wp50 = r(p50)
local wp99 = r(p99)
local wmax = r(max)
tempvar weight_squared
generate double `weight_squared' = calibration_weight_p99^2
quietly summarize calibration_weight_p99, meanonly
local sum_weight = r(sum)
quietly summarize `weight_squared', meanonly
local sum_weight_squared = r(sum)
local ess = (`sum_weight'^2) / `sum_weight_squared'
drop `weight_squared'
post `weight_post' ("p99_truncated") ///
    (.) (.) (`wmean') (`wsd') (`wmin') (`wp1') (`wp50') ///
    (`wp99') (`wmax') (`ess')
postclose `weight_post'

preserve
use "`weight_results'", clear
save "output/women_calibration_weight_diagnostics.dta", replace
export delimited using "output/women_calibration_weight_diagnostics.csv", replace
restore

/* ---------------- Balance table before and after weighting ------------ */
tempname balance_post
tempfile balance_results
postfile `balance_post' ///
    str55 characteristic str12 age_band str12 measure ///
    double target unweighted weighted smd_unweighted smd_weighted ///
    using "`balance_results'", replace

* Age-band proportions.
forvalues g = 1/9 {
    local group_label : label cal_ageband_label `g'
    scalar __target = N_age[1, `g'] / `target_total'

    quietly summarize eb_a`g', meanonly
    scalar __unweighted = r(mean)
    quietly summarize eb_a`g' [aw=calibration_weight], meanonly
    scalar __weighted = r(mean)

    * Treat each indicator as a 0/1 continuous variable. Its SD is
    * sqrt(p(1-p)); stddiffi then applies the pooled-SD SMD formula.
    quietly stddiffi `=scalar(__unweighted)' `=sqrt(scalar(__unweighted) * (1 - scalar(__unweighted)))' ///
				     `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_u = r(std_diff)
    quietly stddiffi `=scalar(__weighted)' `=sqrt(scalar(__weighted) * (1 - scalar(__weighted)))' ///
					 `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_w = r(std_diff)
    post `balance_post' ("Age `group_label'") ("Overall") ("proportion") ///
        (scalar(__target)) (scalar(__unweighted)) (scalar(__weighted)) (scalar(__smd_u)) (scalar(__smd_w))
}

local binary_labels low_bmi previous_fracture corticosteroids rheumatoid_arthritis ///
    current_smoking alcohol_3plus parental_hip_fracture
forvalues j = 1/`n_binary' {
    local variable : word `j' of `binary_variables'
    local target_matrix : word `j' of `binary_target_matrices'
    local characteristic : word `j' of `binary_labels'

    * Overall target is the age-count-weighted target prevalence.
    scalar __target_count = 0
    forvalues g = 1/9 {
        scalar __target_count = scalar(__target_count) + `target_matrix'[1, `g'] * N_age[1, `g']
    }
    scalar __target = scalar(__target_count) / `target_total'

    quietly summarize `variable', meanonly
    scalar __unweighted = r(mean)
    quietly summarize `variable' [aw=calibration_weight], meanonly
    scalar __weighted = r(mean)

    quietly stddiffi `=scalar(__unweighted)' `=sqrt(scalar(__unweighted) * (1 - scalar(__unweighted)))' ///
                     `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_u = r(std_diff)
    quietly stddiffi `=scalar(__weighted)' `=sqrt(scalar(__weighted) * (1 - scalar(__weighted)))' ///
                     `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_w = r(std_diff)
    post `balance_post' ("`characteristic'") ("Overall") ("proportion") ///
        (scalar(__target)) (scalar(__unweighted)) (scalar(__weighted)) (scalar(__smd_u)) (scalar(__smd_w))

    forvalues g = 1/9 {
        local group_label : label cal_ageband_label `g'
        scalar __target = `target_matrix'[1, `g']

        quietly summarize `variable' if cal_ageband == `g', meanonly
        scalar __unweighted = r(mean)
        quietly summarize `variable' if cal_ageband == `g' [aw=calibration_weight], meanonly
        scalar __weighted = r(mean)

        quietly stddiffi `=scalar(__unweighted)' `=sqrt(scalar(__unweighted) * (1 - scalar(__unweighted)))' ///
						 `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
        scalar __smd_u = r(std_diff)
        quietly stddiffi `=scalar(__weighted)' `=sqrt(scalar(__weighted) * (1 - scalar(__weighted)))' ///
						 `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
        scalar __smd_w = r(std_diff)
        post `balance_post' ("`characteristic'") ("`group_label'") ("proportion") (scalar(__target)) ///
            (scalar(__unweighted)) (scalar(__weighted)) (scalar(__smd_u)) (scalar(__smd_w))
    }
}

* BMI mean and SD, overall and within age band.
forvalues g = 0/9 {
    if `g' == 0 {
        local group_label "Overall"
        local restriction ""
        scalar __target_mean = scalar(target_bmi_mean_all)
        scalar __target_sd = scalar(target_bmi_sd_all)
    }
    else {
        local group_label : label cal_ageband_label `g'
        local restriction "if cal_ageband == `g'"
        scalar __target_mean = M_bmi[1, `g']
        scalar __target_sd = SD_bmi[1, `g']
    }

    quietly summarize cal_bmi `restriction'
    scalar __unweighted_mean = r(mean)
    scalar __unweighted_sd = r(sd)
    quietly summarize cal_bmi `restriction' [aw=calibration_weight]
    scalar __weighted_mean = r(mean)
    scalar __weighted_sd = r(sd)

    quietly stddiffi `=scalar(__unweighted_mean)' `=scalar(__unweighted_sd)' ///
				     `=scalar(__target_mean)' `=scalar(__target_sd)'
    scalar __smd_u = r(std_diff)
    quietly stddiffi `=scalar(__weighted_mean)' `=scalar(__weighted_sd)' ///
                     `=scalar(__target_mean)' `=scalar(__target_sd)'
    scalar __smd_w = r(std_diff)
    post `balance_post' ("BMI") ("`group_label'") ("mean") (scalar(__target_mean)) (scalar(__unweighted_mean)) ///
        (scalar(__weighted_mean)) (scalar(__smd_u)) (scalar(__smd_w))
    post `balance_post' ("BMI") ("`group_label'") ("SD") (scalar(__target_sd)) (scalar(__unweighted_sd)) ///
        (scalar(__weighted_sd)) (.) (.)
}

postclose `balance_post'
preserve
use "`balance_results'", clear
generate double abs_smd_weighted = abs(smd_weighted)
sort characteristic age_band measure
save "output/women_calibration_balance.dta", replace
export delimited using "output/women_calibration_balance.csv", replace
restore

/* ---------------- Save patient-level analysis dataset ---------------- */
drop eb_a* eb_b* eb_m* eb_v*
compress
save "output/dxa_eligible_women_calibration_weights.dta", replace

display as result "Calibration weighting completed for `source_n' DXA-eligible women."
display as result "Patient-level weights: output/dxa_eligible_women_calibration_weights.dta"
display as result "Balance table: output/women_calibration_balance.csv"
display as result "Weight diagnostics: output/women_calibration_weight_diagnostics.csv"
display as result "Missingness table: output/women_calibration_missingness.csv"
