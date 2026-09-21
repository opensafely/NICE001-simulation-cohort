/*Sequential conditional models for QFracture variables

Purpose:Fit sex-specific, calibration-weighted conditional models in the
OpenSAFELY DXA-eligible cohort. 

Outputs:
Release package:
output/release/sequential_model_coefficients_release.csv
output/release/sequential_model_manifest.csv
output/release/sequential_model_mapping_rules.csv

Internal review only:
output/internal/sequential_model_diagnostics_internal.csv
output/internal/sequential_model_category_diagnostics_internal.csv
output/internal/sequential_model_cohort_diagnostics_internal.csv
*/

version 16.1
clear all
set more off


/* --------------------------- User settings --------------------------- */
local women_data "output/dxa_eligible_women_calibration_weights.dta"
local men_data   "output/dxa_eligible_men_calibration_weights.dta"
local release_dir  "output/release"
local internal_dir "output/internal"

local weight_variable calibration_weight
local minimum_cell 50

* OpenSAFELY release controls. 
local sdc_threshold 7
local sdc_round_base 5

local ethnicity_levels "1 2 3 4 5 6 7 8 9"
local ethnicity_base 1

capture mkdir "`release_dir'"
capture mkdir "`internal_dir'"

capture erase "`release_dir'/sequential_model_coefficients_release.csv"
capture erase "`release_dir'/sequential_model_manifest.csv"
capture erase "`release_dir'/sequential_model_mapping_rules.csv"


/* ----------------------- Machine-readable outputs -------------------- */
tempname coefficient_post model_post category_post cohort_post mapping_post
tempfile coefficient_results model_results category_results cohort_results mapping_results

postfile `coefficient_post' ///
    str8 sex int sequence str40 model_id str32 outcome ///
    str12 model_type str24 branch str32 equation str100 term ///
    double estimate standard_error n_estimation_rounded df_model ///
    using "`coefficient_results'", replace

postfile `model_post' ///
    str8 sex int sequence str40 model_id str32 outcome ///
    str12 model_type str24 branch str40 expected_levels ///
    double base_category str500 restriction ///
    str2045 predictors str500 logical_rule ///
    double n_estimation n_events n_nonevents n_categories ///
    sum_weight effective_sample_size converged ///
    weighted_observed weighted_predicted brier auc return_code ///
    str28 status ///
    using "`model_results'", replace

postfile `category_post' ///
    str8 sex int sequence str40 model_id str32 outcome ///
    double category n_category weighted_observed weighted_predicted ///
    using "`category_results'", replace

postfile `cohort_post' ///
    str8 sex str32 stage double n n_pct_of_input sum_weight ///
    effective_sample_size ///
    using "`cohort_results'", replace

postfile `mapping_post' ///
    str8 sex str32 target_variable str24 handling ///
    str500 rule ///
    using "`mapping_results'", replace

* Programs need access to the temporary postfile handles.
global SCM_COEFFICIENT_POST `coefficient_post'
global SCM_MODEL_POST       `model_post'
global SCM_CATEGORY_POST    `category_post'
global SCM_MINIMUM_CELL     `minimum_cell'
global SCM_WEIGHT_VARIABLE  `weight_variable'
global SCM_SDC_THRESHOLD    `sdc_threshold'
global SCM_SDC_ROUND_BASE   `sdc_round_base'


/* ---------------------- Coefficient export helper -------------------- */
capture program drop scm_post_coefficients
program define scm_post_coefficients
    version 16.1
    syntax, SEX(string) SEQUENCE(integer) MODELID(string) OUTCOME(string) MODELTYPE(string) [BRANCH(string)]

    if "`branch'" == "" local branch "all"
    tempname coefficient variance omitted
    matrix `coefficient' = e(b)
    matrix `variance' = e(V)

    quietly _ms_omit_info `coefficient'
    matrix `omitted' = r(omit)

    local terms : colnames `coefficient'
    local equations : coleq `coefficient'
    local n_coefficients = colsof(`coefficient')
    local n_estimation = e(N)
    local n_estimation_rounded = round(`n_estimation', $SCM_SDC_ROUND_BASE)
    local df_model = e(df_m)

    * This should be impossible because coefficients are posted only after
    * the final estimation-sample cell checks. Retain as a release safeguard.
    if `n_estimation' <= $SCM_SDC_THRESHOLD {
        display as error "Refusing to post coefficients based on <=7 observations."
        exit 459
    }

    forvalues j = 1/`n_coefficients' {
        if `omitted'[1, `j'] == 0 {
            local term : word `j' of `terms'
            local equation : word `j' of `equations'
            local beta = `coefficient'[1, `j']
            local se = sqrt(`variance'[`j', `j'])

            post $SCM_COEFFICIENT_POST ///
                ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
                ("`modeltype'") ("`branch'") ("`equation'") ///
                ("`term'") (`beta') (`se') ///
                (`n_estimation_rounded') (`df_model')
        }
    }
end


/* -------------------------- Binary models ---------------------------- */
capture program drop scm_fit_logit
program define scm_fit_logit
    version 16.1
    syntax, SEX(string) SEQUENCE(integer) MODELID(string) OUTCOME(name) PREDICTORS(string) [RESTRICTION(string) BRANCH(string) RULE(string)]

    if `"`restriction'"' == "" local restriction "1"
    if "`branch'" == "" local branch "all"
    if `"`rule'"' == "" local rule "Bernoulli draw from predicted probability"

    local weight "$SCM_WEIGHT_VARIABLE"
    local minimum_cell $SCM_MINIMUM_CELL

    capture confirm numeric variable `outcome'
    if _rc {
        post $SCM_MODEL_POST ///
            ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
            ("logit") ("`branch'") ("0 1") (0) ("`restriction'") ///
            (`"`predictors'"') (`"`rule'"') ///
            (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (111) ///
            ("MISSING_OUTCOME")
        exit
    }
    quietly count if (`restriction') & !missing(`outcome') & !inlist(`outcome', 0, 1)
    if r(N) > 0 {
        post $SCM_MODEL_POST ///
            ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
            ("logit") ("`branch'") ("0 1") (0) ("`restriction'") ///
            (`"`predictors'"') (`"`rule'"') ///
            (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (459) ///
            ("NONBINARY_OUTCOME")
        exit
    }

    quietly count if (`restriction') & `outcome' == 0 & !missing(`weight')
    local n0 = r(N)
    quietly count if (`restriction') & `outcome' == 1 & !missing(`weight')
    local n1 = r(N)
    local eligible_n = `n0' + `n1'
    if `n0' < `minimum_cell' | `n1' < `minimum_cell' {
        post $SCM_MODEL_POST ///
            ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
            ("logit") ("`branch'") ("0 1") (0) ("`restriction'") ///
            (`"`predictors'"') (`"`rule'"') ///
            (`eligible_n') (`n1') (`n0') (2) (.) (.) (.) ///
            (.) (.) (.) (.) (2001) ("INSUFFICIENT_CELL")
        exit
    }

    capture noisily logit `outcome' `predictors' if (`restriction') [pweight=`weight'], vce(robust)
    local model_rc = _rc
    if `model_rc' {
        post $SCM_MODEL_POST ///
            ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
            ("logit") ("`branch'") ("0 1") (0) ("`restriction'") ///
            (`"`predictors'"') (`"`rule'"') ///
            (.) (`n1') (`n0') (2) (.) (.) (.) ///
            (.) (.) (.) (.) (`model_rc') ("FIT_FAILED")
        exit
    }

    tempvar model_sample weight_squared predicted squared_error
    generate byte `model_sample' = e(sample)

    quietly count if `model_sample'
    local n_estimation = r(N)
    quietly count if `model_sample' & `outcome' == 1
    local n_events = r(N)
    local n_nonevents = `n_estimation' - `n_events'

    quietly summarize `weight' if `model_sample', meanonly
    local sum_weight = r(sum)
    generate double `weight_squared' = `weight'^2 if `model_sample'
    quietly summarize `weight_squared' if `model_sample', meanonly
    local sum_weight_squared = r(sum)
    local ess = (`sum_weight'^2) / `sum_weight_squared'

    quietly summarize `outcome' if `model_sample' [aw=`weight'], meanonly
    local observed = r(mean)

    predict double `predicted' if `model_sample', pr
    quietly summarize `predicted' if `model_sample' [aw=`weight'], meanonly
    local predicted_mean = r(mean)

    generate double `squared_error' = (`outcome' - `predicted')^2 if `model_sample'
    quietly summarize `squared_error' if `model_sample' [aw=`weight'], meanonly
    local brier = r(mean)

    local auc = .
    capture quietly lroc, nograph
    if !_rc local auc = r(area)

    local converged = e(converged)
    local status "OK"
    if `converged' != 1 local status "NONCONVERGED"
    if `n_events' < `minimum_cell' | `n_nonevents' < `minimum_cell' {
        local status "SMALL_ESTIMATION_CELL"
    }

    post $SCM_MODEL_POST ///
        ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
        ("logit") ("`branch'") ("0 1") (0) ("`restriction'") ///
        (`"`predictors'"') (`"`rule'"') ///
        (`n_estimation') (`n_events') (`n_nonevents') (2) ///
        (`sum_weight') (`ess') (`converged') (`observed') ///
        (`predicted_mean') (`brier') (`auc') (0) ("`status'")

    * Export coefficients only when the converged model also passes the
    * minimum-cell check in the actual estimation sample.
    if "`status'" == "OK" {
        scm_post_coefficients, sex("`sex'") sequence(`sequence') modelid("`modelid'") outcome("`outcome'") modeltype("logit") branch("`branch'")
    }
end


/* ----------------------- Multinomial models -------------------------- */
capture program drop scm_fit_mlogit
program define scm_fit_mlogit
    version 16.1
    syntax, SEX(string) SEQUENCE(integer) MODELID(string) OUTCOME(name) ///
        PREDICTORS(string) BASE(integer) EXPECTEDLEVELS(string) ///
        [RESTRICTION(string) BRANCH(string) RULE(string)]

    if `"`restriction'"' == "" local restriction "1"
    if "`branch'" == "" local branch "all"
    if `"`rule'"' == "" local rule "Categorical draw from softmax probabilities"

    local weight "$SCM_WEIGHT_VARIABLE"
    local minimum_cell $SCM_MINIMUM_CELL
    local n_categories : word count `expectedlevels'

    capture confirm numeric variable `outcome'
    if _rc {
        post $SCM_MODEL_POST ///
            ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
            ("mlogit") ("`branch'") ("`expectedlevels'") (`base') ///
            ("`restriction'") ///
            (`"`predictors'"') (`"`rule'"') ///
            (.) (.) (.) (`n_categories') (.) (.) (.) ///
            (.) (.) (.) (.) (111) ("MISSING_OUTCOME")
        exit
    }

    tempvar expected_category
    generate byte `expected_category' = 0 if (`restriction') & !missing(`outcome')
    local insufficient = 0
    local eligible_n = 0
    foreach category of local expectedlevels {
        replace `expected_category' = 1 if (`restriction') & `outcome' == `category'
        quietly count if (`restriction') & `outcome' == `category' & !missing(`weight')
        local category_n = r(N)
        local eligible_n = `eligible_n' + `category_n'
        if `category_n' < `minimum_cell' local insufficient = 1
    }

    quietly count if (`restriction') & !missing(`outcome') & `expected_category' == 0
    if r(N) > 0 {
        post $SCM_MODEL_POST ///
            ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
            ("mlogit") ("`branch'") ("`expectedlevels'") (`base') ///
            ("`restriction'") ///
            (`"`predictors'"') (`"`rule'"') ///
            (`eligible_n') (.) (.) (`n_categories') (.) (.) (.) ///
            (.) (.) (.) (.) (459) ("UNEXPECTED_CATEGORY")
        exit
    }
    if `insufficient' {
        post $SCM_MODEL_POST ///
            ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
            ("mlogit") ("`branch'") ("`expectedlevels'") (`base') ///
            ("`restriction'") ///
            (`"`predictors'"') (`"`rule'"') ///
            (`eligible_n') (.) (.) (`n_categories') (.) (.) (.) ///
            (.) (.) (.) (.) (2001) ("INSUFFICIENT_CELL")
        exit
    }

    capture noisily mlogit `outcome' `predictors' if (`restriction') [pweight=`weight'], baseoutcome(`base') vce(robust)
    local model_rc = _rc
    if `model_rc' {
        post $SCM_MODEL_POST ///
            ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
            ("mlogit") ("`branch'") ("`expectedlevels'") (`base') ///
            ("`restriction'") ///
            (`"`predictors'"') (`"`rule'"') ///
            (.) (.) (.) (`n_categories') (.) (.) (.) ///
            (.) (.) (.) (.) (`model_rc') ("FIT_FAILED")
        exit
    }

    tempvar model_sample weight_squared multicategory_brier
    generate byte `model_sample' = e(sample)
    quietly count if `model_sample'
    local n_estimation = r(N)

    quietly summarize `weight' if `model_sample', meanonly
    local sum_weight = r(sum)
    generate double `weight_squared' = `weight'^2 if `model_sample'
    quietly summarize `weight_squared' if `model_sample', meanonly
    local sum_weight_squared = r(sum)
    local ess = (`sum_weight'^2) / `sum_weight_squared'

    generate double `multicategory_brier' = 0 if `model_sample'
    local small_estimation_cell = 0

    foreach category of local expectedlevels {
        tempvar observed_category predicted_category category_error
        generate byte `observed_category' = (`outcome' == `category') if `model_sample'
        predict double `predicted_category' if `model_sample', outcome(`category')
        generate double `category_error' = (`observed_category' - `predicted_category')^2 if `model_sample'
        replace `multicategory_brier' = `multicategory_brier' + `category_error' if `model_sample'

        quietly count if `model_sample' & `outcome' == `category'
        local category_n = r(N)
        if `category_n' < `minimum_cell' {
            local small_estimation_cell = 1
        }
        quietly summarize `observed_category' if `model_sample' [aw=`weight'], meanonly
        local category_observed = r(mean)
        quietly summarize `predicted_category' if `model_sample' [aw=`weight'], meanonly
        local category_predicted = r(mean)

        post $SCM_CATEGORY_POST ///
            ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
            (`category') (`category_n') (`category_observed') ///
            (`category_predicted')
    }

    quietly summarize `multicategory_brier' if `model_sample' [aw=`weight'], meanonly
    local brier = r(mean)

    local converged = e(converged)
    local status "OK"
    if `converged' != 1 local status "NONCONVERGED"
    if `small_estimation_cell' == 1 {
        local status "SMALL_ESTIMATION_CELL"
    }

    post $SCM_MODEL_POST ///
        ("`sex'") (`sequence') ("`modelid'") ("`outcome'") ///
        ("mlogit") ("`branch'") ("`expectedlevels'") (`base') ///
        ("`restriction'") ///
        (`"`predictors'"') (`"`rule'"') ///
        (`n_estimation') (.) (.) (`n_categories') ///
        (`sum_weight') (`ess') (`converged') ///
        (.) (.) (`brier') (.) (0) ("`status'")

    * Export coefficients only when every multinomial outcome category in
    * the actual estimation sample passes the minimum-cell safeguard.
    if "`status'" == "OK" {
        scm_post_coefficients, sex("`sex'") sequence(`sequence') modelid("`modelid'") outcome("`outcome'") modeltype("mlogit") branch("`branch'")
    }
end


/* -------------------------- Model sequence --------------------------- */
foreach sex in women men {
    if "`sex'" == "women" local input_data "`women_data'"
    if "`sex'" == "men"   local input_data "`men_data'"
    capture confirm file "`input_data'"
    if _rc {
        display as error "Input file not found: `input_data'"
        exit 601
    }

    use "`input_data'", clear
    capture confirm variable `weight_variable'
    if _rc {
        display as error "Weight variable not found: `weight_variable'"
        exit 111
    }

    local required_variables ///
        age qfracture_calculable cal_bmi cal_prevfx cal_smoke ///
        cal_alcohol cal_steroid cal_ra cal_secondary ///
        fh_parental_hip_fracture qf_ethrisk qf_smoke_cat ///
        qf_alcohol_cat6 fh_parental_osteoporosis ///
        b_type1 b_liver b_malabsorption b_renal b_type2 ///
        b_asthmacopd b_anycancer b_cvd b_epilepsy2 b_parkinsons ///
        b_dementia b_sle b_antidepressant b_falls ///
        dx_hip_fracture_n dx_vertebral_fracture_n ///
        dx_wrist_fracture_n dx_proximal_humerus_fracture_n

    if "`sex'" == "women" {
        local required_variables `required_variables' b_endocrine b_hrt_oest
    }
    if "`sex'" == "men" {
        local required_variables `required_variables' b_carehome
    }

    foreach variable of local required_variables {
        capture confirm variable `variable'
        if _rc {
            display as error "Required variable is missing: `variable'"
            exit 111
        }
    }

    assert inrange(qfracture_calculable, 0, 1) if !missing(qfracture_calculable)
    assert `weight_variable' > 0 if !missing(`weight_variable')

    quietly count
    local n_input = r(N)

    foreach stage in weighted_input model_complete {
        local restriction "`weight_variable' > 0 & !missing(`weight_variable')"
        if "`stage'" == "model_complete" {
            local restriction "!missing(age, cal_bmi, cal_prevfx, cal_smoke, cal_alcohol, cal_steroid, cal_ra, cal_secondary, `weight_variable') & `weight_variable' > 0"
        }

        quietly count if `restriction'
        local stage_n = r(N)
        local stage_pct = 100 * `stage_n' / `n_input'

        quietly summarize `weight_variable' if `restriction', meanonly
        local stage_sum_weight = r(sum)
        tempvar stage_weight_squared
        generate double `stage_weight_squared' = `weight_variable'^2 if `restriction'
        quietly summarize `stage_weight_squared' if `restriction', meanonly
        local stage_sum_weight_squared = r(sum)
        local stage_ess = (`stage_sum_weight'^2) / `stage_sum_weight_squared'
        drop `stage_weight_squared'

        post `cohort_post' ///
            ("`sex'") ("`stage'") (`stage_n') (`stage_pct') ///
            (`stage_sum_weight') (`stage_ess')
    }

    keep if !missing(age, cal_bmi, cal_prevfx, cal_smoke, cal_alcohol, cal_steroid, cal_ra, cal_secondary, `weight_variable') & `weight_variable' > 0
    quietly count
    if r(N) == 0 {
        display as error "No `sex' remain in the modelling cohort."
        exit 2000
    }

    * Portable nonlinear terms.
    generate double age10 = (age - 70) / 10
    generate double age10_sq = age10^2
    generate double bmi5 = (cal_bmi - 27) / 5
    generate double bmi5_sq = bmi5^2

    * Branch outcomes and combined fracture-site indicators.
    generate byte smoke_former = qf_smoke_cat == 1 if inlist(qf_smoke_cat, 0, 1)
    generate byte fx_hip = dx_hip_fracture_n > 0 if !missing(dx_hip_fracture_n)
    generate byte fx_vertebral = dx_vertebral_fracture_n > 0 if !missing(dx_vertebral_fracture_n)
    generate byte fx_wrist = dx_wrist_fracture_n > 0 if !missing(dx_wrist_fracture_n)
    generate byte fx_proximal_humerus = dx_proximal_humerus_fracture_n > 0 if !missing(dx_proximal_humerus_fracture_n)

    * Full-cohort predictor blocks. 
    local continuous "c.age10 c.age10_sq c.bmi5 c.bmi5_sq"
    local fixed_core "`continuous' i.cal_prevfx i.cal_steroid i.cal_ra i.cal_secondary"
    local predictors_common "`fixed_core' i.cal_smoke i.cal_alcohol"
    local predictors_ethnicity "`predictors_common'"
    local predictors_smoking "`fixed_core' i.cal_alcohol"
    local predictors_alcohol "`fixed_core' i.cal_smoke"
    local predictors_behaviour "`predictors_common'"
    local predictors_family "`predictors_common'"

    /* 1. Ethnicity. */
    scm_fit_mlogit, sex("`sex'") sequence(1) ///
        modelid("01_ethnicity") outcome(qf_ethrisk) ///
        predictors(`"`predictors_ethnicity'"') ///
        restriction("!missing(qf_ethrisk)") ///
        base(`ethnicity_base') ///
        expectedlevels(`"`ethnicity_levels'"') ///
        rule("Draw one ethnicity level; use level 1 as the reference")

    /* 2. Smoking detail conditional on fixed FRAX current smoking. */
    scm_fit_logit, sex("`sex'") sequence(2) ///
        modelid("02a_smoking_noncurrent") outcome(smoke_former) ///
        predictors(`"`predictors_smoking'"') ///
        restriction("cal_smoke == 0 & !missing(smoke_former)") ///
        branch("not_current") ///
        rule("If FRAX current smoking=0, draw 0=never or 1=former")

    scm_fit_mlogit, sex("`sex'") sequence(2) ///
        modelid("02b_smoking_current") outcome(qf_smoke_cat) ///
        predictors(`"`predictors_smoking'"') ///
        restriction("cal_smoke == 1 & !missing(qf_smoke_cat)") ///
        branch("current") ///
        base(2) expectedlevels("2 3 4") ///
        rule("If FRAX current smoking=1, draw light, moderate or heavy")

    /* 3. Alcohol detail conditional on fixed FRAX alcohol >=3 units/day. */
    scm_fit_mlogit, sex("`sex'") sequence(3) ///
        modelid("03a_alcohol_below3") outcome(qf_alcohol_cat6) ///
        predictors(`"`predictors_alcohol'"') ///
        restriction("cal_alcohol == 0 & !missing(qf_alcohol_cat6)") ///
        branch("below_3_units") ///
        base(0) expectedlevels("0 1 2") ///
        rule("If FRAX alcohol=0, draw none, trivial or light")

    scm_fit_mlogit, sex("`sex'") sequence(3) ///
        modelid("03b_alcohol_3plus") outcome(qf_alcohol_cat6) ///
        predictors(`"`predictors_alcohol'"') ///
        restriction("cal_alcohol == 1 & !missing(qf_alcohol_cat6)") ///
        branch("3plus_units") ///
        base(3) expectedlevels("3 4 5") ///
        rule("If FRAX alcohol=1, draw moderate, heavy or very heavy")

    /* 4. Additional family history. */
    scm_fit_logit, sex("`sex'") sequence(4) ///
        modelid("04_parental_osteoporosis") ///
        outcome(fh_parental_osteoporosis) ///
        predictors(`"`predictors_behaviour'"') ///
        restriction("fh_parental_hip_fracture == 0") ///
        branch("no_parental_hip") ///
        rule("If fh_parental_hip_fracture == 0, draw parental osteoporosis; set QFracture family history to parental hip OR generated parental osteoporosis")

    /* 5. Components represented within FRAX secondary osteoporosis. */
    scm_fit_logit, sex("`sex'") sequence(5) ///
        modelid("05a_type1_diabetes") outcome(b_type1) ///
        predictors(`"`predictors_family'"') ///
        restriction("cal_secondary == 1") branch("secondary_yes") ///
        rule("Generate only if FRAX secondary osteoporosis=1; otherwise set 0")

    local predictors_secondary_1 "`predictors_family' i.b_type1"
    scm_fit_logit, sex("`sex'") sequence(5) ///
        modelid("05b_liver") outcome(b_liver) ///
        predictors(`"`predictors_secondary_1'"') ///
        restriction("cal_secondary == 1") branch("secondary_yes") ///
        rule("Generate only if FRAX secondary osteoporosis=1; otherwise set 0")

    local predictors_secondary_2 "`predictors_secondary_1' i.b_liver"
    scm_fit_logit, sex("`sex'") sequence(5) ///
        modelid("05c_malabsorption") outcome(b_malabsorption) ///
        predictors(`"`predictors_secondary_2'"') ///
        restriction("cal_secondary == 1") branch("secondary_yes") ///
        rule("Generate only if FRAX secondary osteoporosis=1; otherwise set 0")

    local predictors_secondary_3 "`predictors_secondary_2' i.b_malabsorption"
    scm_fit_logit, sex("`sex'") sequence(5) ///
        modelid("05d_renal") outcome(b_renal) ///
        predictors(`"`predictors_secondary_3'"') ///
        restriction("cal_secondary == 1") branch("secondary_yes") ///
        rule("Generate only if FRAX secondary osteoporosis=1; otherwise set 0")

    local predictors_components "`predictors_secondary_3' i.b_renal"

    /* 6. Female endocrine disease spans both secondary-osteoporosis strata. */
    if "`sex'" == "women" {
        scm_fit_logit, sex("`sex'") sequence(6) ///
            modelid("06a_endocrine_secondary0") outcome(b_endocrine) ///
            predictors(`"`predictors_family'"') ///
            restriction("cal_secondary == 0") ///
            branch("secondary_no") ///
            rule("Generate female QFracture endocrine disease within secondary=0")

        scm_fit_logit, sex("`sex'") sequence(6) ///
            modelid("06b_endocrine_secondary1") outcome(b_endocrine) ///
            predictors(`"`predictors_components'"') ///
            restriction("cal_secondary == 1") ///
            branch("secondary_yes") ///
            rule("Generate female QFracture endocrine disease within secondary=1")

        local predictors_components "`predictors_components' i.b_endocrine"
    }

    /* 7. Other chronic conditions. */
    scm_fit_logit, sex("`sex'") sequence(7) ///
        modelid("07a_type2_diabetes") outcome(b_type2) ///
        predictors(`"`predictors_components'"')

	local predictors_07a "`predictors_components' i.b_type2"
    scm_fit_logit, sex("`sex'") sequence(7) ///
        modelid("07b_asthma_copd") outcome(b_asthmacopd) ///
        predictors(`"`predictors_07a'"')

	local predictors_07b "`predictors_07a' i.b_asthmacopd"
    scm_fit_logit, sex("`sex'") sequence(7) ///
        modelid("07c_any_cancer") outcome(b_anycancer) ///
        predictors(`"`predictors_07b'"')

	local predictors_07c "`predictors_07b' i.b_anycancer"
    scm_fit_logit, sex("`sex'") sequence(7) ///
        modelid("07d_cvd") outcome(b_cvd) ///
        predictors(`"`predictors_07c'"')

	local predictors_07d "`predictors_07c' i.b_cvd"
    scm_fit_logit, sex("`sex'") sequence(7) ///
        modelid("07e_epilepsy") outcome(b_epilepsy2) ///
        predictors(`"`predictors_07d'"')

	local predictors_07e "`predictors_07d' i.b_epilepsy2"
    scm_fit_logit, sex("`sex'") sequence(7) ///
        modelid("07f_parkinsons") outcome(b_parkinsons) ///
        predictors(`"`predictors_07e'"')

	local predictors_07f "`predictors_07e' i.b_parkinsons"
    scm_fit_logit, sex("`sex'") sequence(7) ///
        modelid("07g_dementia") outcome(b_dementia) ///
        predictors(`"`predictors_07f'"')

    /* 8. SLE completes the fixed RA indicator. */
    local predictors_07g "`predictors_07f' i.b_dementia"	
    scm_fit_logit, sex("`sex'") sequence(8) ///
        modelid("08_sle_without_ra") outcome(b_sle) ///
        predictors(`"`predictors_07g'"') ///
        restriction("cal_ra == 0") branch("ra_absent") ///
        rule("If cal_ra == 0, draw b_sle; Set QFracture RA/SLE to fixed RA OR generated SLE")

    /* 9. Medication and residence variables. */
    local predictors_08 "`predictors_07g' i.b_sle"
    scm_fit_logit, sex("`sex'") sequence(9) ///
        modelid("09a_antidepressant") outcome(b_antidepressant) ///
        predictors(`"`predictors_08'"')

	local predictors_09 "`predictors_08' i.b_antidepressant"
    if "`sex'" == "women" {
        scm_fit_logit, sex("`sex'") sequence(9) ///
            modelid("09b_hrt_oestrogen") outcome(b_hrt_oest) ///
            predictors(`"`predictors_09'"')
        local predictors_09 "`predictors_09' i.b_hrt_oest"
    }

    if "`sex'" == "men" {
        scm_fit_logit, sex("`sex'") sequence(9) ///
            modelid("09b_carehome") outcome(b_carehome) ///
            predictors(`"`predictors_09'"')
        local predictors_09 "`predictors_09' i.b_carehome"
    }

    /* 10. Falls after neurological conditions and medicines. */
    scm_fit_logit, sex("`sex'") sequence(10) ///
        modelid("10_falls") outcome(b_falls) ///
        predictors(`"`predictors_09'"')

    /* 11. Non-mutually-exclusive prior-fracture site indicators. */
    local predictors_10 "`predictors_09' i.b_falls"
    scm_fit_logit, sex("`sex'") sequence(11) ///
        modelid("11a_fracture_hip") outcome(fx_hip) ///
        predictors(`"`predictors_10'"') ///
        restriction("cal_prevfx == 1") branch("prior_fracture") ///
        rule("Generate only if previous fracture=1; sites are not mutually exclusive")

    local predictors_fx1 "`predictors_10' i.fx_hip"
    scm_fit_logit, sex("`sex'") sequence(11) ///
        modelid("11b_fracture_vertebral") outcome(fx_vertebral) ///
        predictors(`"`predictors_fx1'"') ///
        restriction("cal_prevfx == 1") branch("prior_fracture") ///
        rule("Generate only if previous fracture=1; sites are not mutually exclusive")

    local predictors_fx2 "`predictors_fx1' i.fx_vertebral"
    scm_fit_logit, sex("`sex'") sequence(11) ///
        modelid("11c_fracture_wrist") outcome(fx_wrist) ///
        predictors(`"`predictors_fx2'"') ///
        restriction("cal_prevfx == 1") branch("prior_fracture") ///
        rule("Generate only if previous fracture=1; sites are not mutually exclusive")

    local predictors_fx3 "`predictors_fx2' i.fx_wrist"
    scm_fit_logit, sex("`sex'") sequence(11) ///
        modelid("11d_fracture_proximal_humerus") ///
        outcome(fx_proximal_humerus) ///
        predictors(`"`predictors_fx3'"') ///
        restriction("cal_prevfx == 1") branch("prior_fracture") ///
        rule("Generate only if previous fracture=1; sites are not mutually exclusive")

    * Direct mappings and post-model construction rules.
    post `mapping_post' ("`sex'") ("age_terms") ("transformation") ///
        ("age10=(age-70)/10; age10_sq=age10^2")
    post `mapping_post' ("`sex'") ("bmi_terms") ("transformation") ///
        ("bmi5=(uncapped_BMI-27)/5; bmi5_sq=bmi5^2; cap BMI to 20-40 only when calculating QFracture")
    post `mapping_post' ("`sex'") ("__logit_formula") ("probability") ///
        ("p=1/(1+exp(-eta)), where eta is the exported intercept plus coefficient times transformed predictor terms")
    post `mapping_post' ("`sex'") ("__mlogit_formula") ("probability") ///
        ("For nonbase k: p_k=exp(eta_k)/(1+sum exp(eta_j)); for the base category: p_base=1/(1+sum exp(eta_j))")
    post `mapping_post' ("`sex'") ("__random_draw") ("simulation") ///
        ("Use a prespecified reproducible random seed; draw Bernoulli or categorical values in sequence order")
    post `mapping_post' ("`sex'") ("b_corticosteroids") ("direct_mapping") ///
        ("Use the fixed FRAX glucocorticoid indicator as QFracture corticosteroid exposure")
    post `mapping_post' ("`sex'") ("qf_smoke_cat") ("conditional_mapping") ///
        ("FRAX current=0 uses model 02a; current=1 uses model 02b")
    post `mapping_post' ("`sex'") ("qf_alcohol_cat6") ("conditional_mapping") ///
        ("FRAX alcohol<3 uses model 03a; alcohol>=3 uses model 03b")
    post `mapping_post' ("`sex'") ("fh_osteoporosis") ("logical_or") ///
        ("Fixed parental hip fracture OR generated parental osteoporosis")
    post `mapping_post' ("`sex'") ("b_ra_sle") ("logical_or") ///
        ("Fixed rheumatoid arthritis OR generated SLE")
    post `mapping_post' ("`sex'") ("secondary_components") ("conditional_zero") ///
        ("Type 1 diabetes, liver, malabsorption and renal are generated when secondary osteoporosis=1 and set to 0 otherwise")
    post `mapping_post' ("`sex'") ("prior_fracture_sites") ("conditional_zero") ///
        ("Generate four site indicators when previous fracture=1 and set all to 0 otherwise")
}

/* --------------------------- Export package -------------------------- */
postclose `coefficient_post'
postclose `model_post'
postclose `category_post'
postclose `cohort_post'
postclose `mapping_post'

preserve
use "`model_results'", clear
sort sex sequence model_id
export delimited using "`internal_dir'/sequential_model_diagnostics_internal.csv", replace
restore

preserve
use "`category_results'", clear
sort sex sequence model_id category
export delimited using "`internal_dir'/sequential_model_category_diagnostics_internal.csv", replace
restore

preserve
use "`cohort_results'", clear
sort sex stage
export delimited using "`internal_dir'/sequential_model_cohort_diagnostics_internal.csv", replace
restore

* Validate the complete package before creating any release files.
preserve
use "`model_results'", clear
quietly count
local manifest_rows = r(N)
quietly count if status != "OK"
local failed_models = r(N)
restore

preserve
use "`coefficient_results'", clear
quietly count
local coefficient_rows = r(N)
capture assert n_estimation_rounded > `sdc_threshold' & mod(n_estimation_rounded, `sdc_round_base') == 0 & !missing(df_model)
local coefficient_sdc_error = _rc
restore

preserve
use "`mapping_results'", clear
quietly count
local mapping_rows = r(N)
restore

if `failed_models' > 0 | `manifest_rows' == 0 | `manifest_rows' > 5000 {
    macro drop SCM_COEFFICIENT_POST SCM_MODEL_POST SCM_CATEGORY_POST SCM_MINIMUM_CELL SCM_WEIGHT_VARIABLE SCM_SDC_THRESHOLD SCM_SDC_ROUND_BASE
    display as error "Release package not created: review internal model diagnostics."
    exit 459
}

if `coefficient_rows' == 0 | `coefficient_rows' > 5000 | `coefficient_sdc_error' != 0 {
    macro drop SCM_COEFFICIENT_POST SCM_MODEL_POST SCM_CATEGORY_POST SCM_MINIMUM_CELL SCM_WEIGHT_VARIABLE SCM_SDC_THRESHOLD SCM_SDC_ROUND_BASE
    display as error "Release package not created: coefficient release checks failed."
    exit 459
}

if `mapping_rows' == 0 | `mapping_rows' > 5000 {
    macro drop SCM_COEFFICIENT_POST SCM_MODEL_POST SCM_CATEGORY_POST SCM_MINIMUM_CELL SCM_WEIGHT_VARIABLE SCM_SDC_THRESHOLD SCM_SDC_ROUND_BASE
    display as error "Release package not created: mapping release checks failed."
    exit 459
}

* File 1 of 3: release-safe model coefficients. Estimation N is rounded to
* the nearest 5; models based on small final outcome cells are not exported.
preserve
use "`coefficient_results'", clear
sort sex sequence model_id equation term
export delimited using "`release_dir'/sequential_model_coefficients_release.csv", replace
restore

* File 2 of 3: model order, predictors, branch rules and reference levels.
preserve
use "`model_results'", clear
sort sex sequence model_id
keep sex sequence model_id outcome model_type branch expected_levels base_category restriction predictors logical_rule status
export delimited using "`release_dir'/sequential_model_manifest.csv", replace
restore

* File 3 of 3: transformations, direct mappings and construction rules.
preserve
use "`mapping_results'", clear
sort sex target_variable
export delimited using "`release_dir'/sequential_model_mapping_rules.csv", replace
restore

macro drop SCM_COEFFICIENT_POST SCM_MODEL_POST SCM_CATEGORY_POST SCM_MINIMUM_CELL SCM_WEIGHT_VARIABLE SCM_SDC_THRESHOLD SCM_SDC_ROUND_BASE
display as result "Sequential conditional modelling completed."
display as result "Review internal diagnostics, then request only the three files in output/release through Airlock."
