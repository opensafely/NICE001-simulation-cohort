/*Sex-specific weighted Table 1: simple CSV output for women and men.*/

version 16.1
clear all
set more off

/* ----------------------------- Settings ----------------------------- */
local women_data "output/dxa_eligible_women_calibration_weights.dta"
local men_data   "output/dxa_eligible_men_calibration_weights.dta"
local outdir     "output/release"
local internal_dir "output/internal"
local sdc_threshold 7
local minimum_ess 8
local round_base 5
local write_internal_audit 1
local weight_variable calibration_weight
local cohort_scope "model_complete"  

capture mkdir "output"
capture mkdir "`outdir'"
capture mkdir "`internal_dir'"

tempname table_post
tempfile table_results
postfile `table_post' ///
    int row_order int category_order str8 sex str16 cohort_scope ///
    str32 variable str100 characteristic str80 category str16 statistic ///
    double n_cohort n_observed n_missing n_denominator n_category ///
    weighted_pct weighted_mean weighted_sd byte cell_unsafe double minimum_support_ess ///
    using "`table_results'", replace

/* Define t1_support. */
capture program drop t1_support
program define t1_support, rclass
    version 16.1
    syntax varname(numeric), WEIGHT(varname) THRESHOLD(integer) Miness(real)
    quietly count if `varlist' == 1
    local n = r(N)
    local sumw = 0
    local ess = .
    if `n' > 0 {
        quietly summarize `weight' if `varlist' == 1, meanonly
        local sumw = r(sum)
        tempvar wsq
        quietly generate double `wsq' = `weight'^2 if `varlist' == 1
        quietly summarize `wsq', meanonly
        local ess = (`sumw'^2)/r(sum)
    }
    local unsafe = inrange(`n', 1, `threshold')
    if `n' > 0 & (missing(`ess') | `ess' < `miness' - 1e-8) local unsafe = 1
    return scalar n = `n'
    return scalar sumw = `sumw'
    return scalar ess = `ess'
    return scalar unsafe = `unsafe'
end

/* Define program that calculates each characteristic’s Table 1 rows. */
capture program drop t1_weighted_variable
program define t1_weighted_variable
    version 16.1
    syntax varname(numeric), HANDLE(name) SEX(string) WEIGHT(varname) ///
        SCOPE(string) ORDER(integer) THRESHOLD(integer) Miness(real) [LEVELS(numlist)]
    local characteristic : variable label `varlist'
    if `"`characteristic'"' == "" local characteristic "`varlist'"
    quietly count
    local n_cohort = r(N)
    quietly count if !missing(`varlist')
    local n_observed = r(N)
    local n_missing = `n_cohort' - `n_observed'
    tempvar observed missing_flag category_flag complement
    quietly generate byte `observed' = !missing(`varlist')
    quietly generate byte `missing_flag' = missing(`varlist')
    quietly generate byte `category_flag' = 0
    quietly generate byte `complement' = 0
    t1_support `observed', weight(`weight') threshold(`threshold') miness(`miness')
    local unsafe = r(unsafe)
    local ess_min = r(ess)
    t1_support `missing_flag', weight(`weight') threshold(`threshold') miness(`miness')
    local unsafe = max(`unsafe', r(unsafe))
    local ess_min = min(`ess_min', r(ess))
    quietly summarize `missing_flag' [aw=`weight'], meanonly
    local missing_pct = 100 * r(mean)

    if "`levels'" == "" {
        local weighted_mean = .
        local weighted_sd = .
        if `n_observed' > 0 {
            quietly summarize `varlist' [aw=`weight']
            local weighted_mean = r(mean)
            local weighted_sd = r(sd)
        }
        post `handle' (`order') (1) ("`sex'") ("`scope'") ///
            ("`varlist'") (`"`characteristic'"') ("") ("mean_sd") ///
            (`n_cohort') (`n_observed') (`n_missing') (`n_observed') (.) ///
            (.) (`weighted_mean') (`weighted_sd') (`unsafe') (`ess_min')
    }
    else {
        /* levels(1) is shorthand for a binary Yes row; validate both 0 and 1. */
        local checklevels "`levels'"
        if "`levels'" == "1" local checklevels "0 1"
        local allowed "0"
        foreach value of numlist `checklevels' {
            local allowed "`allowed' | `varlist' == `value'"
        }
        quietly assert (`allowed') if !missing(`varlist')
        /* Check all categories and their observed complements before posting. */
        foreach value of numlist `checklevels' {
            quietly replace `category_flag' = (`varlist' == `value')
            quietly replace `complement' = !missing(`varlist') & `varlist' != `value'
            t1_support `category_flag', weight(`weight') threshold(`threshold') miness(`miness')
            local unsafe = max(`unsafe', r(unsafe))
            local ess_min = min(`ess_min', r(ess))
            t1_support `complement', weight(`weight') threshold(`threshold') miness(`miness')
            local unsafe = max(`unsafe', r(unsafe))
            local ess_min = min(`ess_min', r(ess))
        }
        local category_order = 0
        foreach value of numlist `levels' {
            local category_order = `category_order' + 1
            local category : label (`varlist') `value'
            if `"`category'"' == "" local category "`value'"
            quietly replace `category_flag' = (`varlist' == `value')
            quietly count if `category_flag' == 1
            local n_category = r(N)
            quietly summarize `category_flag' [aw=`weight'], meanonly
            local weighted_pct = 100 * r(mean)
            post `handle' (`order') (`category_order') ("`sex'") ("`scope'") ///
                ("`varlist'") (`"`characteristic'"') (`"`category'"') ("percent") ///
                (`n_cohort') (`n_observed') (`n_missing') (`n_cohort') (`n_category') ///
                (`weighted_pct') (.) (.) (`unsafe') (`ess_min')
        }
    }
    post `handle' (`order') (99) ("`sex'") ("`scope'") ///
        ("`varlist'") (`"`characteristic'"') ("Missing") ("percent") ///
        (`n_cohort') (`n_observed') (`n_missing') (`n_cohort') (`n_missing') ///
        (`missing_pct') (.) (.) (`unsafe') (`ess_min')
end

/* -------------------------- Women and men --------------------------- */
foreach sex in women men {
    if "`sex'" == "women" local input_data "`women_data'"
    if "`sex'" == "men"   local input_data "`men_data'"
    confirm file "`input_data'"
    use "`input_data'", clear

    local core_variables ///
        age cal_bmi cal_prevfx cal_smoke cal_alcohol cal_steroid cal_ra cal_secondary
    local additional_variables ///
        qf_ethrisk qf_smoke_cat qf_alcohol_cat6 ///
        fh_parental_hip_fracture fh_parental_osteoporosis ///
        b_type1 b_liver b_malabsorption b_renal b_type2 b_asthmacopd ///
        b_anycancer b_cvd b_epilepsy2 b_parkinsons b_dementia b_sle ///
        b_antidepressant b_carehome b_falls ///
        dx_hip_fracture_n dx_vertebral_fracture_n ///
        dx_wrist_fracture_n dx_proximal_humerus_fracture_n
    foreach variable in `weight_variable' `core_variables' `additional_variables' {
        confirm numeric variable `variable'
    }
    if "`sex'" == "women" {
        confirm numeric variable b_endocrine
        confirm numeric variable b_hrt_oest
    }

    /* Check that a sex-specific input has not been mixed up. */
    confirm string variable sex
    if "`sex'" == "women" assert lower(strtrim(sex)) == "female"
    if "`sex'" == "men"   assert lower(strtrim(sex)) == "male"

    quietly count
    local n_input = r(N)
    keep if !missing(`weight_variable') & `weight_variable' > 0
    if "`cohort_scope'" == "model_complete" {
        keep if !missing(age, cal_bmi, cal_prevfx, cal_smoke, cal_alcohol, ///
            cal_steroid, cal_ra, cal_secondary, `weight_variable')
    }
    quietly count
    local n_cohort = r(N)
    if `n_cohort' == 0 {
        display as error "No `sex' remain in the requested Table 1 population."
        exit 2000
    }
    if `n_cohort' <= `sdc_threshold' {
        display as error "Insufficient patients for sex-specific release candidates."
        exit 459
    }
    capture confirm variable patient_id
    if !_rc quietly isid patient_id
    quietly summarize `weight_variable', meanonly
    local weight_mean = r(mean)
    tempvar table_weight weight_squared cohort_flag
    generate double `table_weight' = `weight_variable'/`weight_mean'
    generate double `weight_squared' = `table_weight'^2
    generate byte `cohort_flag' = 1
    quietly assert !missing(`table_weight', `weight_squared') & `table_weight' > 0 & `weight_squared' > 0
    t1_support `cohort_flag', weight(`table_weight') threshold(`sdc_threshold') miness(`minimum_ess')
    local cohort_unsafe = r(unsafe)
    local cohort_ess = r(ess)

    /* Missing is preserved in all derived indicators. */
    generate byte t1_ageband = .
    replace t1_ageband = floor((age - 50)/10) + 1 if inrange(age, 50, 99)
    assert !missing(t1_ageband) if !missing(age)
    generate byte t1_lowbmi = cal_bmi <= 19 if !missing(cal_bmi)

    generate byte t1_fx_hip = dx_hip_fracture_n > 0 if !missing(dx_hip_fracture_n)
    generate byte t1_fx_vertebral = dx_vertebral_fracture_n > 0 if !missing(dx_vertebral_fracture_n)
    generate byte t1_fx_wrist = dx_wrist_fracture_n > 0 if !missing(dx_wrist_fracture_n)
    generate byte t1_fx_humerus = dx_proximal_humerus_fracture_n > 0 if !missing(dx_proximal_humerus_fracture_n)

    generate byte t1_ra_sle = .
    replace t1_ra_sle = 1 if cal_ra == 1 | b_sle == 1
    replace t1_ra_sle = 0 if cal_ra == 0 & b_sle == 0
    generate byte t1_parental_combined = .
    replace t1_parental_combined = 1 if fh_parental_hip_fracture == 1 | fh_parental_osteoporosis == 1
    replace t1_parental_combined = 0 if fh_parental_hip_fracture == 0 & fh_parental_osteoporosis == 0

    /* Explicit, consistent labels allow Women/Men to be placed side by side. */
    label define t1_age_label 1 "50-59" 2 "60-69" 3 "70-79" 4 "80-89" 5 "90-99", replace
    label values t1_ageband t1_age_label
    label define t1_eth_label 1 "White or not stated" 2 "Indian" 3 "Pakistani" ///
        4 "Bangladeshi" 5 "Other Asian" 6 "Black Caribbean" ///
        7 "Black African" 8 "Chinese" 9 "Other or mixed ethnicity", replace
    label values qf_ethrisk t1_eth_label
    label define t1_smoke_label 0 "Never smoker" 1 "Former smoker" ///
        2 "Current: <10/day" 3 "Current: 10-19/day" 4 "Current: >=20/day", replace
    label values qf_smoke_cat t1_smoke_label
    label define t1_alcohol_label 0 "None" 1 "<1 unit/day" 2 "1-2 units/day" ///
        3 "3-6 units/day" 4 "7-9 units/day" 5 ">9 units/day", replace
    label values qf_alcohol_cat6 t1_alcohol_label
    label define t1_yesno 0 "No" 1 "Yes", replace

    label variable age "Age, years"
    label variable cal_bmi "BMI used for calibration, kg/m2 (uncapped)"
    label variable t1_ageband "Age group, years"
    label variable qf_ethrisk "Ethnicity"
    label variable qf_smoke_cat "Detailed smoking category"
    label variable qf_alcohol_cat6 "Detailed alcohol category"
    label variable cal_prevfx "Previous fracture (calibration definition)"
    label variable t1_fx_hip "Previous hip fracture"
    label variable t1_fx_vertebral "Previous vertebral fracture"
    label variable t1_fx_wrist "Previous wrist fracture"
    label variable t1_fx_humerus "Previous proximal humerus fracture"
    label variable cal_steroid "Oral corticosteroids: >=2 prescriptions in prior 6 months"
    label variable cal_ra "Rheumatoid arthritis"
    label variable cal_smoke "Current smoking, including current intensity unknown"
    label variable cal_alcohol "Alcohol >=3 units/day"
    label variable cal_secondary "FRAX secondary osteoporosis proxy"
    label variable t1_lowbmi "BMI <=19 kg/m2"
    label variable fh_parental_hip_fracture "Parental hip fracture"
    label variable fh_parental_osteoporosis "Parental osteoporosis"
    label variable t1_parental_combined "Parental hip fracture or osteoporosis"
    label variable b_type1 "Type 1 diabetes"
    label variable b_liver "Chronic liver disease"
    label variable b_malabsorption "Malabsorption"
    label variable b_renal "Chronic kidney disease (QFracture definition)"
    label variable b_type2 "Type 2 diabetes"
    label variable b_asthmacopd "Asthma or COPD"
    label variable b_anycancer "Cancer"
    label variable b_cvd "Cardiovascular disease"
    label variable b_epilepsy2 "Epilepsy (QFracture definition)"
    label variable b_parkinsons "Parkinson's disease"
    label variable b_dementia "Dementia"
    label variable b_sle "Systemic lupus erythematosus"
    label variable t1_ra_sle "Rheumatoid arthritis or SLE"
    label variable b_antidepressant "Antidepressant treatment"
    label variable b_carehome "Care-home residency"
    label variable b_falls "History of falls"
    if "`sex'" == "women" {
        label variable b_endocrine "Endocrine disorder (female QFracture predictor)"
        label variable b_hrt_oest "Oestrogen-only HRT"
    }

    assert inrange(qf_ethrisk, 1, 9) & qf_ethrisk == floor(qf_ethrisk) if !missing(qf_ethrisk)
    assert inrange(qf_smoke_cat, 0, 4) & qf_smoke_cat == floor(qf_smoke_cat) if !missing(qf_smoke_cat)
    assert inrange(qf_alcohol_cat6, 0, 5) & qf_alcohol_cat6 == floor(qf_alcohol_cat6) if !missing(qf_alcohol_cat6)

    /* Population N is an actual patient count; the remaining rows are weighted. */
    post `table_post' (0) (0) ("`sex'") ("`cohort_scope'") ///
        ("_population") ("Population") ("Total") ("count") ///
        (`n_cohort') (`n_cohort') (0) (`n_cohort') (`n_cohort') ///
        (.) (.) (.) (`cohort_unsafe') (`cohort_ess')

    /* Same broad order as the supplied Table 1; additional calibration inputs last. */
    local variables age t1_ageband cal_bmi cal_prevfx t1_fx_hip ///
        t1_fx_vertebral t1_fx_wrist t1_fx_humerus qf_smoke_cat ///
        qf_alcohol_cat6 qf_ethrisk b_anycancer b_asthmacopd b_carehome ///
        b_cvd b_dementia b_endocrine b_epilepsy2 b_falls b_liver ///
        b_malabsorption b_parkinsons t1_ra_sle b_renal b_type1 b_type2 ///
        t1_parental_combined b_antidepressant cal_steroid b_hrt_oest ///
        cal_ra cal_secondary cal_smoke cal_alcohol t1_lowbmi ///
        fh_parental_hip_fracture fh_parental_osteoporosis b_sle
    local order = 0
    foreach variable of local variables {
        local order = `order' + 1
        local options "levels(1)"
        if inlist("`variable'", "age", "cal_bmi") local options ""
        if "`variable'" == "t1_ageband" local options "levels(1/5)"
        if "`variable'" == "qf_smoke_cat" local options "levels(0/4)"
        if "`variable'" == "qf_alcohol_cat6" local options "levels(0/5)"
        if "`variable'" == "qf_ethrisk" local options "levels(1/9)"
        /* Female-only predictors are not applicable to men. */
        if "`sex'" == "men" & inlist("`variable'", "b_endocrine", "b_hrt_oest") continue
        if inlist("`options'", "levels(1)", "levels(1 0)") {
            label values `variable' t1_yesno
        }
        t1_weighted_variable `variable', handle(`table_post') sex("`sex'") ///
            weight(`table_weight') scope("`cohort_scope'") order(`order') ///
            threshold(`sdc_threshold') miness(`minimum_ess') `options'
    }
}
postclose `table_post'

/* ------------------------ Disclosure control ------------------------- */
use "`table_results'", clear
isid sex row_order category_order
quietly assert n_category <= n_denominator if !missing(n_category)
quietly assert inrange(weighted_pct, 0, 100.000001) if !missing(weighted_pct)
quietly assert n_observed + n_missing == n_cohort if statistic != "not_applicable"

/* Share suppression across related variables to protect hidden small cells. */
generate str32 sdc_block = variable
quietly replace sdc_block = "age" if inlist(variable, "age", "t1_ageband")
quietly replace sdc_block = "bmi" if inlist(variable, "cal_bmi", "t1_lowbmi")
quietly replace sdc_block = "fracture" if inlist(variable, "cal_prevfx", "t1_fx_hip", "t1_fx_vertebral", "t1_fx_wrist", "t1_fx_humerus")
quietly replace sdc_block = "smoking" if inlist(variable, "qf_smoke_cat", "cal_smoke")
quietly replace sdc_block = "alcohol" if inlist(variable, "qf_alcohol_cat6", "cal_alcohol")
quietly replace sdc_block = "ra_sle" if inlist(variable, "cal_ra", "b_sle", "t1_ra_sle")
quietly replace sdc_block = "family" if inlist(variable, "fh_parental_hip_fracture", "fh_parental_osteoporosis", "t1_parental_combined")
bysort sex sdc_block: egen byte block_unsafe = max(cell_unsafe)
bysort sex: egen byte cohort_unsafe = max(cond(statistic == "count", cell_unsafe, 0))
quietly replace block_unsafe = max(block_unsafe, cohort_unsafe)
sort row_order category_order sex

/* Exact normalized weighted counts stay inside the secure audit. */
generate double weighted_count = (weighted_pct/100) * n_cohort if statistic == "percent"

/* Highly sensitive audit only: no raw fields enter release exports. */
if `write_internal_audit' {
    export delimited using "`internal_dir'/table1_weighted_by_sex_sdc_audit.csv", replace
}

/* Actual n is nonmissing N for continuous rows and category N for categorical rows. */
generate str20 actual_n = ""
generate str20 weighted_n = ""
generate str20 mean_release = ""
generate str20 sd_release = ""
generate str20 weighted_percent = ""

quietly replace actual_n = strtrim(string(round(n_cohort, `round_base'), "%15.0f")) if statistic == "count"
quietly replace actual_n = strtrim(string(round(n_observed, `round_base'), "%15.0f")) if statistic == "mean_sd"
quietly replace actual_n = strtrim(string(round(n_category, `round_base'), "%15.0f")) if statistic == "percent"
quietly replace weighted_n = strtrim(string(round(weighted_count, `round_base'), "%15.0f")) if statistic == "percent"
quietly replace weighted_percent = strtrim(string(100 * ///
    round(weighted_count, `round_base') / round(n_cohort, `round_base'), "%6.1f")) ///
    if statistic == "percent"
quietly replace mean_release = "NA" if statistic == "mean_sd"
quietly replace sd_release = "NA" if statistic == "mean_sd"
quietly replace mean_release = strtrim(string(weighted_mean, "%9.1f")) if statistic == "mean_sd" & !missing(weighted_mean)
quietly replace sd_release = strtrim(string(weighted_sd, "%9.1f")) if statistic == "mean_sd" & !missing(weighted_sd)

/* Apply the existing block-level suppression to every displayed statistic. */
quietly replace actual_n = "[REDACTED]" if block_unsafe == 1 & inlist(statistic, "count", "mean_sd", "percent")
quietly replace weighted_n = "[REDACTED]" if block_unsafe == 1 & statistic == "percent"
quietly replace weighted_percent = "[REDACTED]" if block_unsafe == 1 & statistic == "percent"
quietly replace mean_release = "[REDACTED]" if block_unsafe == 1 & statistic == "mean_sd"
quietly replace sd_release = "[REDACTED]" if block_unsafe == 1 & statistic == "mean_sd"

/* Compact display after all support checks have been done. */
drop if category == "Missing" & inlist(variable, "age", "cal_bmi")
drop if category == "Missing" & n_missing == 0
sort sex row_order category_order

foreach sex in women men {
    preserve
    keep if sex == "`sex'"
    quietly count
    if r(N) > 5000 {
        display as error "Release table exceeds the OpenSAFELY 5000-row limit."
        exit 459
    }
    keep characteristic category actual_n weighted_n mean_release sd_release weighted_percent
    rename mean_release weighted_mean
    rename sd_release weighted_sd
    order characteristic category actual_n weighted_n weighted_mean weighted_sd weighted_percent
    export delimited using "`outdir'/table1_weighted_`sex'.csv", replace
    restore
}

capture program drop t1_support
capture program drop t1_weighted_variable
