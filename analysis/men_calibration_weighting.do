/* Purpose: Estimate entropy-balancing weights that align men in the OpenSAFELY
DXA-eligible cohort with men treatment-eligible cohort.*/

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
keep if dxa_eligible == 1 & qfracture_calculable == 1
keep if lower(strtrim(sex)) == "male"
quietly count
local source_n = r(N)

/* ---------------- Male aggregate targets ----------------------------- */
local target_total 7795

* Age bands: 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84, 85-89, 90+.
matrix N_age       = (983, 729, 1041, 1052, 1072, 1038, 981, 622, 277)
matrix C_lowbmi    = (  4,   5,    4,    7,    9,    8,   5,   4,   3)
matrix C_prevfx    = (971, 700,  940,  880,  800,  609, 443, 211,  73)
matrix C_smoke     = (344, 253,  328,  322,  283,  215, 200, 102,  46)
matrix C_alcohol   = (373, 273,  410,  407,  425,  365, 328, 209,  85)
matrix C_ra        = ( 32,  23,   44,   63,   51,   55,  59,  47,  30)
matrix C_steroid   = ( 65,  52,   82,   82,   87,   79,  61,  41,   8)
matrix C_secondary = (200, 190,  333,  343,  395,  444, 497, 315, 143)

matrix M_bmi = ( ///
	28.33825025432353, 28.15349794238686, 27.28424591738713, ///
	26.83146387832700, 26.61016791044777, 26.70385356454722, ///
	26.48358817533126, 26.35691318327976, 25.87256317689531)
	
matrix SD_bmi = ( ///
	5.299982879582247, 4.831162235352270, 4.312783060169213, ///
	4.031426091677906, 4.023226122300709, 3.766794547452983, ///
	3.437824681456744, 3.009535375707634, 3.401955682293866)

scalar target_bmi_mean_all = 27.04243745991023
scalar target_bmi_sd_all   =  4.198788601895486

* Fail if any required target is missing.
assert !missing(`target_total') & `target_total' > 0
assert !missing(scalar(target_bmi_mean_all))
assert !missing(scalar(target_bmi_sd_all)) & scalar(target_bmi_sd_all) > 0
foreach target_matrix in N_age C_lowbmi C_prevfx C_smoke C_alcohol C_ra C_steroid C_secondary M_bmi SD_bmi {
    forvalues g = 1/9 {
        assert !missing(`target_matrix'[1, `g'])
    }
}

* Verify that the age counts reproduce stated total.
scalar __target_n_check = 0
forvalues g = 1/9 {
    scalar __target_n_check = scalar(__target_n_check) + N_age[1, `g']
}
assert scalar(__target_n_check) == `target_total'

* Convert target counts into within-age proportions.
foreach characteristic in lowbmi prevfx smoke alcohol ra steroid secondary {
    matrix P_`characteristic' = J(1, 9, .)
    forvalues g = 1/9 {
        matrix P_`characteristic'[1, `g'] = C_`characteristic'[1, `g'] / N_age[1, `g']
    }
}

* Calibrate binary risks and BMI within four broader age groups. The original
* five-year targets are retained below for detailed balance diagnostics.
* Broad groups: 50-59, 60-69, 70-79, 80+.
matrix B_start = (1, 3, 5, 7)
matrix B_end   = (2, 4, 6, 9)
matrix N_broad = J(1, 4, .)
matrix M_bmi_broad = J(1, 4, .)
matrix VAR_bmi_broad = J(1, 4, .)
matrix SD_bmi_broad = J(1, 4, .)

foreach characteristic in lowbmi prevfx smoke alcohol ra steroid secondary {
    matrix C_`characteristic'_broad = J(1, 4, 0)
    matrix P_`characteristic'_broad = J(1, 4, .)
}

forvalues h = 1/4 {
    local first = B_start[1, `h']
    local last  = B_end[1, `h']
    scalar __broad_n = 0
    scalar __broad_bmi_sum = 0
    foreach characteristic in lowbmi prevfx smoke alcohol ra steroid secondary {
        scalar __broad_count_`characteristic' = 0
    }

    forvalues g = `first'/`last' {
        scalar __broad_n = scalar(__broad_n) + N_age[1, `g']
        scalar __broad_bmi_sum = scalar(__broad_bmi_sum) + N_age[1, `g'] * M_bmi[1, `g']
        foreach characteristic in lowbmi prevfx smoke alcohol ra steroid secondary {
            scalar __broad_count_`characteristic' = scalar(__broad_count_`characteristic') + C_`characteristic'[1, `g']
        }
    }

    matrix N_broad[1, `h'] = scalar(__broad_n)
    matrix M_bmi_broad[1, `h'] = scalar(__broad_bmi_sum) / scalar(__broad_n)

    foreach characteristic in lowbmi prevfx smoke alcohol ra steroid secondary {
        matrix C_`characteristic'_broad[1, `h'] = scalar(__broad_count_`characteristic')
        matrix P_`characteristic'_broad[1, `h'] = scalar(__broad_count_`characteristic') / scalar(__broad_n)
    }

    * Pool the five-year sample SDs into a broad-age population variance.
    scalar __broad_ss = 0
    forvalues g = `first'/`last' {
        scalar __broad_ss = scalar(__broad_ss) + (N_age[1, `g'] - 1) * SD_bmi[1, `g']^2 + ///
            N_age[1, `g'] * (M_bmi[1, `g'] - M_bmi_broad[1, `h'])^2
    }
    matrix VAR_bmi_broad[1, `h'] = scalar(__broad_ss) / scalar(__broad_n)
    matrix SD_bmi_broad[1, `h'] = sqrt(scalar(__broad_ss) / (scalar(__broad_n) - 1))
}

scalar __broad_n_check = 0
forvalues h = 1/4 {
    scalar __broad_n_check = scalar(__broad_n_check) + N_broad[1, `h']
}
assert scalar(__broad_n_check) == `target_total'

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

generate byte cal_agebroad = .
replace cal_agebroad = 1 if inrange(cal_ageband, 1, 2)
replace cal_agebroad = 2 if inrange(cal_ageband, 3, 4)
replace cal_agebroad = 3 if inrange(cal_ageband, 5, 6)
replace cal_agebroad = 4 if inrange(cal_ageband, 7, 9)
assert !missing(cal_agebroad)
label define cal_agebroad_label 1 "50-59" 2 "60-69" 3 "70-79" 4 "80+"
label values cal_agebroad cal_agebroad_label

generate double cal_bmi = bmi_raw if inrange(bmi_raw, 10, 80)
generate byte cal_lowbmi = cal_bmi <= 19 if !missing(cal_bmi)

egen byte __previous_fracture_max = rowmax(dx_hip_fracture_n dx_vertebral_fracture_n dx_wrist_fracture_n dx_proximal_humerus_fracture_n)
generate byte cal_prevfx = __previous_fracture_max > 0 if !missing(__previous_fracture_max)
drop __previous_fracture_max


generate byte cal_smoke = .
replace cal_smoke = 0 if inlist(qf_smoke_cat, 0, 1)
replace cal_smoke = 1 if inrange(qf_smoke_cat, 2, 4) | smoking_current_unknown == 1
assert !missing(cal_smoke) if smoking_no_record == 0
generate byte cal_alcohol = inrange(qf_alcohol_cat6, 3, 5) if !missing(qf_alcohol_cat6)
generate byte cal_ra = b_ra if !missing(b_ra)
generate byte cal_steroid = b_corticosteroids if !missing(b_corticosteroids)
generate byte cal_secondary = secondary_osteoporosis if !missing(secondary_osteoporosis)

label variable cal_bmi      "Uncapped BMI used for calibration"
label variable cal_lowbmi   "BMI <=19 kg/m2"
label variable cal_prevfx   "Any previous hip/vertebral/wrist/proximal-humerus fracture"
label variable cal_smoke    "Current smoking"
label variable cal_alcohol  "Alcohol >=3 units/day"
label variable cal_ra       "Rheumatoid arthritis"
label variable cal_steroid  ">=2 oral corticosteroid prescriptions in prior 6 months"
label variable cal_secondary "FRAX secondary osteoporosis proxy"

local core_binary_variables cal_prevfx cal_steroid cal_ra cal_smoke
local core_binary_target_matrices P_prevfx_broad P_steroid_broad P_ra_broad P_smoke_broad

local full_binary_variables `core_binary_variables' cal_secondary cal_alcohol
local full_binary_target_matrices `core_binary_target_matrices' P_secondary_broad P_alcohol_broad

local bmi_binary_variables cal_lowbmi
local bmi_binary_target_matrices P_lowbmi_broad

* Variables and five-year target matrices used for the balance table.
local binary_variables `full_binary_variables' `bmi_binary_variables'
local binary_target_matrices P_prevfx P_steroid P_ra P_smoke ///
    P_secondary P_alcohol P_lowbmi

* Variables and broad-age target matrices supplied to the calibration solver.
local calibration_binary_variables `full_binary_variables' `bmi_binary_variables'
local cali_binary_target_matrices ///
    `full_binary_target_matrices' `bmi_binary_target_matrices'

foreach variable of local binary_variables {
    assert inlist(`variable', 0, 1, .)
}

/* ---------------- Missingness report --------------------------------- */
tempname missing_post
tempfile missing_results
postfile `missing_post' ///
    str60 characteristic str12 age_band ///
    double n_total n_observed n_missing pct_missing ///
    using "`missing_results'", replace

local missing_variables ///
    cal_bmi cal_lowbmi cal_prevfx cal_smoke ///
    cal_alcohol cal_secondary cal_ra cal_steroid

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
save "output/men_calibration_missingness.dta", replace
export delimited using "output/men_calibration_missingness.csv", replace
restore

/* ---------------- Feasibility and overlap checks --------------------- */
tempname feasibility_post
tempfile feasibility_results

postfile `feasibility_post' ///
    str20 check_type str80 characteristic ///
    str12 age_band double target_value source_value ///
    source_min source_max source_sd ///
    n_total n_observed n_missing n_zero n_one ///
    str8 status str100 reason ///
    using "`feasibility_results'", replace

/* Check that every five-year age band is represented. */
forvalues g = 1/9 {
    local group_label : label cal_ageband_label `g'
    quietly count if cal_ageband == `g'
    local ntotal = r(N)
    local target_value = N_age[1, `g'] / `target_total'
    local source_value = `ntotal' / `source_n'
    local status "PASS"
    local reason ""
    if `ntotal' == 0 {
        local status "FAIL"
        local reason "No source observations in this age band"
    }
    post `feasibility_post' ///
        ("age_support") ///
        ("Age distribution") ///
        ("`group_label'") ///
        (`target_value') (`source_value') ///
        (.) (.) (.) ///
        (`ntotal') (`ntotal') (0) (.) (.) ///
        ("`status'") ("`reason'")
}

/* Check binary-variable overlap within broad age groups. */
local n_calibration_binary : word count `calibration_binary_variables'
forvalues j = 1/`n_calibration_binary' {
    local variable : word `j' of `calibration_binary_variables'
    local target_matrix : word `j' of `cali_binary_target_matrices'
    local characteristic : variable label `variable'
    forvalues h = 1/4 {
        local group_label : label cal_agebroad_label `h'
        local target_value = `target_matrix'[1, `h']
        quietly count if cal_agebroad == `h'
        local ntotal = r(N)
        quietly count if cal_agebroad == `h' & `variable' == 0
        local n0 = r(N)
        quietly count if cal_agebroad == `h' & `variable' == 1
        local n1 = r(N)
        local nobserved = `n0' + `n1'
        local nmissing = `ntotal' - `nobserved'
        local source_value = .
        if `nobserved' > 0 {
            local source_value = `n1' / `nobserved'
        }
        local status "PASS"
        local reason ""
        if missing(`target_value') | ///
            `target_value' <= 0 | `target_value' >= 1 {
            local status "FAIL"
            local reason "Target proportion is not between 0 and 1"
        }
        if `n0' == 0 | `n1' == 0 {
            local status "FAIL"
            if "`reason'" != "" {
                local reason "`reason'; no observed 0 or 1"
            }
            else {
                local reason "No observed 0 or 1"
            }
        }
        post `feasibility_post' ///
            ("binary_overlap") ///
            ("`characteristic'") ///
            ("`group_label'") ///
            (`target_value') (`source_value') ///
            (0) (1) (.) ///
            (`ntotal') (`nobserved') (`nmissing') (`n0') (`n1') ///
            ("`status'") ("`reason'")
    }
}

/* Check BMI variation and whether the target mean is within source range. */
forvalues h = 1/4 {
    local group_label : label cal_agebroad_label `h'
    local target_value = M_bmi_broad[1, `h']
    quietly count if cal_agebroad == `h'
    local ntotal = r(N)
    quietly summarize cal_bmi if cal_agebroad == `h'
    local nobserved = r(N)
    local source_value = r(mean)
    local source_min = r(min)
    local source_max = r(max)
    local source_sd = r(sd)
    local nmissing = `ntotal' - `nobserved'
    local status "PASS"
    local reason ""
    if `nobserved' < 2 | missing(`source_sd') | `source_sd' == 0 {
        local status "FAIL"
        local reason "Fewer than two BMI values or no BMI variation"
    }
    else if `target_value' < `source_min' | ///
            `target_value' > `source_max' {
        local status "FAIL"
        local reason "Target BMI mean is outside the source BMI range"
    }
    post `feasibility_post' ///
        ("bmi_overlap") ///
        ("BMI") ///
        ("`group_label'") ///
        (`target_value') (`source_value') ///
        (`source_min') (`source_max') (`source_sd') ///
        (`ntotal') (`nobserved') (`nmissing') (.) (.) ///
        ("`status'") ("`reason'")
}
postclose `feasibility_post'

/* Save the report before stopping for any failure. */
preserve
use "`feasibility_results'", clear
sort check_type characteristic age_band
quietly count if status == "FAIL"
local n_feasibility_failures = r(N)
save "output/men_calibration_feasibility.dta", replace
export delimited using "output/men_calibration_feasibility.csv", replace
restore

if `n_feasibility_failures' > 0 {
    display as error "`n_feasibility_failures' feasibility checks failed."
    display as error "See output/men_calibration_feasibility.csv."
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

local f_core_binary ""
local f_full_binary ""
local f_bmi_low_terms ""

forvalues j = 1/`n_calibration_binary' {
    local variable : word `j' of `calibration_binary_variables'
    local target_matrix : word `j' of `cali_binary_target_matrices'
    local is_core : list variable in core_binary_variables
    local is_full : list variable in full_binary_variables
    local is_bmi : list variable in bmi_binary_variables

    forvalues h = 1/4 {
        scalar __target_p = `target_matrix'[1, `h']
        generate double eb_b`j'_`h' = 0
        replace eb_b`j'_`h' = `variable' - scalar(__target_p) if cal_agebroad == `h' & !missing(`variable')
        if `is_full' {
            local f_full_binary "`f_full_binary' eb_b`j'_`h'"
        }
        if `is_core' {
            local f_core_binary "`f_core_binary' eb_b`j'_`h'"
        }
        if `is_bmi' {
            local f_bmi_low_terms "`f_bmi_low_terms' eb_b`j'_`h'"
        }
    }
}

local f_bmi_mean_terms ""
local f_bmi_var_terms ""
forvalues h = 1/4 {
    scalar __target_mean = M_bmi_broad[1, `h']
    scalar __target_var_pop = VAR_bmi_broad[1, `h']
    generate double eb_m`h' = 0
    replace eb_m`h' = (cal_bmi - scalar(__target_mean)) / 10 if cal_agebroad == `h' & !missing(cal_bmi)
    generate double eb_v`h' = 0
    replace eb_v`h' = ((cal_bmi - scalar(__target_mean))^2 - scalar(__target_var_pop)) / 100 ///
        if cal_agebroad == `h' & !missing(cal_bmi)
    local f_bmi_mean_terms "`f_bmi_mean_terms' eb_m`h'"
    local f_bmi_var_terms "`f_bmi_var_terms' eb_v`h'"
}

local f_core        "`f_age' `f_core_binary'"
local f_full        "`f_age' `f_full_binary'"
local f_bmi_mean    "`f_full' `f_bmi_mean_terms'"
local f_bmi_mean_low "`f_bmi_mean' `f_bmi_low_terms'"
local f_bmi_all     "`f_bmi_mean_low' `f_bmi_var_terms'"

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

local targets_bmi_mean "`targets_full'"
local n_bmi_mean_terms : word count `f_bmi_mean_terms'
forvalues j = 1/`n_bmi_mean_terms' {
    local targets_bmi_mean "`targets_bmi_mean' 0"
}

local targets_bmi_mean_low "`targets_bmi_mean'"
local n_bmi_low_terms : word count `f_bmi_low_terms'
forvalues j = 1/`n_bmi_low_terms' {
    local targets_bmi_mean_low "`targets_bmi_mean_low' 0"
}

local targets_bmi_all "`targets_bmi_mean_low'"
local n_bmi_var_terms : word count `f_bmi_var_terms'
forvalues j = 1/`n_bmi_var_terms' {
    local targets_bmi_all "`targets_bmi_all' 0"
}

foreach stage in age core full bmi_mean bmi_mean_low bmi_all {
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
foreach stage in age core full bmi_mean bmi_mean_low bmi_all {
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

* Primary specification: exact five-year age balance, broad-age binary-risk
* balance, broad-age BMI means and broad-age low-BMI prevalence. The bmi_all
* stage additionally calibrates BMI variance and is retained for sensitivity.
local primary_stage bmi_mean_low
generate double calibration_weight = cw_`primary_stage'
label variable calibration_weight "Men calibration weight: exact age + broad risks + BMI mean/low BMI"
label variable cw_bmi_all "Sensitivity weight: primary specification plus broad-age BMI SD"

* Sensitivity weight truncated at the primary weight's 99th percentile.
quietly summarize calibration_weight, detail
scalar __weight_p99 = r(p99)
generate double calibration_weight_p99 = min(calibration_weight, scalar(__weight_p99))
quietly summarize calibration_weight_p99, meanonly
replace calibration_weight_p99 = calibration_weight_p99 / r(mean)
label variable calibration_weight_p99 "Calibration weight truncated at p99 and renormalised"

* Quantify the loss of exact calibration after truncation against every
* constraint in the primary specification.
scalar __p99_max_target_deviation = 0
local final_features "`f_bmi_mean_low'"
local final_targets "`targets_bmi_mean_low'"
local n_final_features : word count `final_features'
forvalues j = 1/`n_final_features' {
    local feature : word `j' of `final_features'
    local target : word `j' of `final_targets'
    quietly summarize `feature' [aw=calibration_weight_p99], meanonly
    scalar __target_deviation = abs(r(mean) - (`target'))
    if scalar(__target_deviation) > scalar(__p99_max_target_deviation) {
        scalar __p99_max_target_deviation = scalar(__target_deviation)
    }
}

/* ---------------- Weight diagnostics --------------------------------- */
tempname weight_post
tempfile weight_results
postfile `weight_post' ///
    str24 stage double converged max_target_deviation ///
    mean sd min p1 p50 p99 max effective_sample_size effective_sample_size_pct ///
    using "`weight_results'", replace

foreach stage in age core full bmi_mean bmi_mean_low bmi_all {
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
    local ess_pct = 100 * `ess' / `source_n'
    drop `weight_squared'

    local result_matrix R_`stage'
    local converged = `result_matrix'[1, 1]
    local max_error = `result_matrix'[1, 2]

    post `weight_post' ("`stage'") ///
        (`converged') (`max_error') ///
        (`wmean') (`wsd') (`wmin') (`wp1') (`wp50') ///
        (`wp99') (`wmax') (`ess') (`ess_pct')
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
local ess_pct = 100 * `ess' / `source_n'
drop `weight_squared'
post `weight_post' ("p99_truncated") ///
    (.) (scalar(__p99_max_target_deviation)) ///
    (`wmean') (`wsd') (`wmin') (`wp1') (`wp50') ///
    (`wp99') (`wmax') (`ess') (`ess_pct')
postclose `weight_post'

preserve
use "`weight_results'", clear
save "output/men_calibration_weight_diagnostics.dta", replace
export delimited using "output/men_calibration_weight_diagnostics.csv", replace
restore

/* ---------------- Missingness after weighting ------------------------ */
tempname weighted_missing_post
tempfile weighted_missing_results
postfile `weighted_missing_post' ///
    str60 characteristic str12 age_band ///
    double n_total n_missing pct_missing_unweighted ///
    pct_missing_weighted pct_missing_weighted_p99 ///
    using "`weighted_missing_results'", replace

foreach variable of local missing_variables {
    local characteristic : variable label `variable'
    tempvar missing_indicator
    generate byte `missing_indicator' = missing(`variable')

    quietly count
    local ntotal = r(N)
    quietly count if `missing_indicator' == 1
    local nmissing = r(N)
    quietly summarize `missing_indicator', meanonly
    local pmissing_unweighted = 100 * r(mean)
    quietly summarize `missing_indicator' [aw=calibration_weight], meanonly
    local pmissing_weighted = 100 * r(mean)
    quietly summarize `missing_indicator' [aw=calibration_weight_p99], meanonly
    local pmissing_weighted_p99 = 100 * r(mean)
    post `weighted_missing_post' ("`characteristic'") ("Overall") ///
        (`ntotal') (`nmissing') (`pmissing_unweighted') ///
        (`pmissing_weighted') (`pmissing_weighted_p99')

    forvalues g = 1/9 {
        local group_label : label cal_ageband_label `g'
        quietly count if cal_ageband == `g'
        local ntotal = r(N)
        quietly count if cal_ageband == `g' & `missing_indicator' == 1
        local nmissing = r(N)
        quietly summarize `missing_indicator' if cal_ageband == `g', meanonly
        local pmissing_unweighted = 100 * r(mean)
        quietly summarize `missing_indicator' if cal_ageband == `g' [aw=calibration_weight], meanonly
        local pmissing_weighted = 100 * r(mean)
        quietly summarize `missing_indicator' if cal_ageband == `g' [aw=calibration_weight_p99], meanonly
        local pmissing_weighted_p99 = 100 * r(mean)
        post `weighted_missing_post' ("`characteristic'") ("`group_label'") ///
            (`ntotal') (`nmissing') (`pmissing_unweighted') ///
            (`pmissing_weighted') (`pmissing_weighted_p99')
    }
    drop `missing_indicator'
}
postclose `weighted_missing_post'

preserve
use "`weighted_missing_results'", clear
sort characteristic age_band
save "output/men_calibration_weighted_missingness.dta", replace
export delimited using "output/men_calibration_weighted_missingness.csv", replace
restore

/* ---------------- Balance table before and after weighting ------------ */
tempname balance_post
tempfile balance_results
postfile `balance_post' ///
    str55 characteristic str12 age_band str12 measure ///
    double target unweighted weighted weighted_p99 ///
    smd_unweighted smd_weighted smd_weighted_p99 ///
    using "`balance_results'", replace

* Age-band proportions.
forvalues g = 1/9 {
    local group_label : label cal_ageband_label `g'
    scalar __target = N_age[1, `g'] / `target_total'

    quietly summarize eb_a`g', meanonly
    scalar __unweighted = r(mean)
    quietly summarize eb_a`g' [aw=calibration_weight], meanonly
    scalar __weighted = r(mean)
    quietly summarize eb_a`g' [aw=calibration_weight_p99], meanonly
    scalar __weighted_p99 = r(mean)

    * Treat each indicator as a 0/1 continuous variable. Its SD is
    * sqrt(p(1-p)); stddiffi then applies the pooled-SD SMD formula.
    quietly stddiffi `=scalar(__unweighted)' `=sqrt(scalar(__unweighted) * (1 - scalar(__unweighted)))' ///
				     `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_u = r(std_diff)
    quietly stddiffi `=scalar(__weighted)' `=sqrt(scalar(__weighted) * (1 - scalar(__weighted)))' ///
					 `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_w = r(std_diff)
    quietly stddiffi `=scalar(__weighted_p99)' `=sqrt(scalar(__weighted_p99) * (1 - scalar(__weighted_p99)))' ///
                     `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_w_p99 = r(std_diff)
    post `balance_post' ("Age `group_label'") ("Overall") ("proportion") ///
        (scalar(__target)) (scalar(__unweighted)) (scalar(__weighted)) ///
        (scalar(__weighted_p99)) (scalar(__smd_u)) (scalar(__smd_w)) ///
        (scalar(__smd_w_p99))
}

local binary_labels previous_fracture corticosteroids rheumatoid_arthritis current_smoking secondary_osteoporosis alcohol_3plus low_bmi
local n_binary : word count `binary_variables'    
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
    quietly summarize `variable' [aw=calibration_weight_p99], meanonly
    scalar __weighted_p99 = r(mean)

    quietly stddiffi `=scalar(__unweighted)' `=sqrt(scalar(__unweighted) * (1 - scalar(__unweighted)))' ///
                     `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_u = r(std_diff)
    quietly stddiffi `=scalar(__weighted)' `=sqrt(scalar(__weighted) * (1 - scalar(__weighted)))' ///
                     `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_w = r(std_diff)
    quietly stddiffi `=scalar(__weighted_p99)' `=sqrt(scalar(__weighted_p99) * (1 - scalar(__weighted_p99)))' ///
                     `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
    scalar __smd_w_p99 = r(std_diff)
    post `balance_post' ("`characteristic'") ("Overall") ("proportion") ///
        (scalar(__target)) (scalar(__unweighted)) (scalar(__weighted)) ///
        (scalar(__weighted_p99)) (scalar(__smd_u)) (scalar(__smd_w)) ///
        (scalar(__smd_w_p99))

    forvalues g = 1/9 {
        local group_label : label cal_ageband_label `g'
        scalar __target = `target_matrix'[1, `g']

        quietly summarize `variable' if cal_ageband == `g', meanonly
        scalar __unweighted = r(mean)
        quietly summarize `variable' if cal_ageband == `g' [aw=calibration_weight], meanonly
        scalar __weighted = r(mean)
        quietly summarize `variable' if cal_ageband == `g' [aw=calibration_weight_p99], meanonly
        scalar __weighted_p99 = r(mean)

        quietly stddiffi `=scalar(__unweighted)' `=sqrt(scalar(__unweighted) * (1 - scalar(__unweighted)))' ///
						 `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
        scalar __smd_u = r(std_diff)
        quietly stddiffi `=scalar(__weighted)' `=sqrt(scalar(__weighted) * (1 - scalar(__weighted)))' ///
						 `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
        scalar __smd_w = r(std_diff)
        quietly stddiffi `=scalar(__weighted_p99)' `=sqrt(scalar(__weighted_p99) * (1 - scalar(__weighted_p99)))' ///
                         `=scalar(__target)' `=sqrt(scalar(__target) * (1 - scalar(__target)))'
        scalar __smd_w_p99 = r(std_diff)
        post `balance_post' ("`characteristic'") ("`group_label'") ("proportion") (scalar(__target)) ///
            (scalar(__unweighted)) (scalar(__weighted)) (scalar(__weighted_p99)) ///
            (scalar(__smd_u)) (scalar(__smd_w)) (scalar(__smd_w_p99))
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
    quietly summarize cal_bmi `restriction' [aw=calibration_weight_p99]
    scalar __weighted_p99_mean = r(mean)
    scalar __weighted_p99_sd = r(sd)

    quietly stddiffi `=scalar(__unweighted_mean)' `=scalar(__unweighted_sd)' ///
				     `=scalar(__target_mean)' `=scalar(__target_sd)'
    scalar __smd_u = r(std_diff)
    quietly stddiffi `=scalar(__weighted_mean)' `=scalar(__weighted_sd)' ///
                     `=scalar(__target_mean)' `=scalar(__target_sd)'
    scalar __smd_w = r(std_diff)
    quietly stddiffi `=scalar(__weighted_p99_mean)' `=scalar(__weighted_p99_sd)' ///
                     `=scalar(__target_mean)' `=scalar(__target_sd)'
    scalar __smd_w_p99 = r(std_diff)
    post `balance_post' ("BMI") ("`group_label'") ("mean") (scalar(__target_mean)) (scalar(__unweighted_mean)) ///
        (scalar(__weighted_mean)) (scalar(__weighted_p99_mean)) ///
        (scalar(__smd_u)) (scalar(__smd_w)) (scalar(__smd_w_p99))
    post `balance_post' ("BMI") ("`group_label'") ("SD") (scalar(__target_sd)) (scalar(__unweighted_sd)) ///
        (scalar(__weighted_sd)) (scalar(__weighted_p99_sd)) (.) (.) (.)
}

postclose `balance_post'
preserve
use "`balance_results'", clear
generate double abs_smd_weighted = abs(smd_weighted)
generate double abs_smd_weighted_p99 = abs(smd_weighted_p99)
generate double abs_diff_weighted = abs(weighted - target)
generate double abs_diff_weighted_p99 = abs(weighted_p99 - target)
sort characteristic age_band measure
save "output/men_calibration_balance.dta", replace
export delimited using "output/men_calibration_balance.csv", replace
restore

/* ---------------- Save patient-level analysis dataset ---------------- */
drop eb_a* eb_b* eb_m* eb_v*
compress
save "output/dxa_eligible_men_calibration_weights.dta", replace

display as result "Calibration weighting completed for `source_n' DXA-eligible men."
display as result "Patient-level weights: output/dxa_eligible_men_calibration_weights.dta"
display as result "Balance table: output/men_calibration_balance.csv"
display as result "Weight diagnostics: output/men_calibration_weight_diagnostics.csv"
display as result "Missingness table: output/men_calibration_missingness.csv"
display as result "Weighted missingness table: output/men_calibration_weighted_missingness.csv"
