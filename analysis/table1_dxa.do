/* Table 1: QFracture predictor characteristics ------------------------- */
version 16.1
clear all
set more off

use "output/qfracture_cohort_with_dxa_flags.dta", clear

*OpenSAFELY disclosure-control programs 

/*Counts of 1–7 are redacted. Other counts are rounded to the nearest 5. Zero is retained.*/
capture program drop sdc_count
program define sdc_count, rclass
    syntax, Count(real)
    if inrange(`count', 1, 7) {
        return local cell "[REDACTED]"
        exit
    }
    local count_rounded = round(`count', 5)
    local cell = strtrim(string(`count_rounded', "%12.0fc"))
    return local cell `"`cell'"'
end

/*For n (%), disclosure control is applied to the raw count first.
The percentage is then calculated from rounded counts.*/
capture program drop sdc_n_pct
program define sdc_n_pct, rclass
    syntax, Count(real) Denominator(real)
    if inrange(`count', 1, 7) | inrange(`denominator', 1, 7) {
        return local cell "[REDACTED]"
        exit
    }

    if `denominator' <= 0 {
        return local cell "NA"
        exit
    }

    local count_rounded = round(`count', 5)
    local denominator_rounded = round(`denominator', 5)
    local percentage = 100 * `count_rounded' / `denominator_rounded'
    local cell = strtrim(string(`count_rounded', "%12.0fc")) + " (" + strtrim(string(`percentage', "%5.1f")) + "%)"
    return local cell `"`cell'"'
end

generate byte prev_hip = inrange(dx_hip_fracture_n, 1, 2)
generate byte prev_vertebral = inrange(dx_vertebral_fracture_n, 1, 2)
generate byte prev_wrist = inrange(dx_wrist_fracture_n, 1, 2)
generate byte prev_humerus = inrange(dx_proximal_humerus_fracture_n, 1, 2)

label variable prev_hip "Previous hip fracture"
label variable prev_vertebral "Previous vertebral fracture"
label variable prev_wrist "Previous wrist fracture"
label variable prev_humerus "Previous proximal humerus fracture"

generate byte table_sex = .
replace table_sex = 1 if lower(strtrim(sex)) == "female"
replace table_sex = 2 if lower(strtrim(sex)) == "male"

label define table_sex_label 1 "Female" 2 "Male"
label values table_sex table_sex_label
label variable table_sex "Sex"

/* QFracture category labels */
label define table_ethnicity_label ///
    1 "White or not stated" ///
    2 "Indian" ///
    3 "Pakistani" ///
    4 "Bangladeshi" ///
    5 "Other Asian" ///
    6 "Black Caribbean" ///
    7 "Black African" ///
    8 "Chinese" ///
    9 "Other or mixed ethnicity"

label values qf_ethrisk table_ethnicity_label
label variable qf_ethrisk "Ethnicity"

label define table_smoking_label ///
    0 "Non-smoker" ///
    1 "Ex-smoker" ///
    2 "Light current smoker (<10/day)" ///
    3 "Moderate current smoker (10-19/day)" ///
    4 "Heavy current smoker (>=20/day)"

label values qf_smoke_cat table_smoking_label
label variable qf_smoke_cat "Smoking status"

label define table_alcohol_label ///
    0 "None" ///
    1 "<1 unit/day" ///
    2 "1-2 units/day" ///
    3 "3-6 units/day" ///
    4 "7-9 units/day" ///
    5 ">9 units/day"

label values qf_alcohol_cat6 table_alcohol_label
label variable qf_alcohol_cat6 "Alcohol consumption"

/* Variable labels for continuous predictors */
label variable age "Age, years"
label variable qf_bmi "BMI, kg/m^2"

/* Variable labels for binary QFracture predictors */
label variable b_anycancer          "Cancer"
label variable b_asthmacopd         "Asthma or COPD"
label variable b_carehome           "Care home"
label variable b_cvd                "Cardiovascular disease"
label variable b_dementia           "Dementia"
label variable b_endocrine          "Endocrine disorder"
label variable b_epilepsy2          "Epilepsy"
label variable b_falls              "History of falls"
label variable b_liver              "Chronic liver disease"
label variable b_malabsorption      "Malabsorption"
label variable b_parkinsons         "Parkinson's disease"
label variable b_ra_sle             "Rheumatoid arthritis or SLE"
label variable b_renal              "Chronic kidney disease"
label variable b_type1              "Type 1 diabetes"
label variable b_type2              "Type 2 diabetes"
label variable fh_osteoporosis      "Parental history of osteoporosis"
label variable b_antidepressant     "Antidepressant treatment"
label variable b_corticosteroids    "Oral corticosteroid treatment"
label variable b_hrt_oest           "Oestrogen-only HRT"

/* Population denominators */
quietly count
local N_overall = r(N)

quietly count if dxa_eligible == 1
local N_high = r(N)

/* Create results dataset */
tempname table1_post
tempfile table1_results

postfile `table1_post' ///
    int row_order ///
    str70 characteristic ///
    str45 category ///
    str30 overall ///
    str30 dxa_eligible ///
    using `table1_results', replace
local order = 1

/* Population size */
sdc_count, count(`N_overall')
local overall_cell `"`r(cell)'"'

sdc_count, count(`N_high')
local high_cell `"`r(cell)'"'

post `table1_post' ///
    (`order') ///
    ("Population, n") ///
    ("") ///
    ("`overall_cell'") ///
    ("`high_cell'")
local order = `order' + 1

/* Continuous-variable helper: mean (SD) */
capture program drop table1_continuous
program define table1_continuous
    syntax varname, Handle(name) Order(integer)
    local row_label : variable label `varlist'

    quietly summarize `varlist'
    local n_overall = r(N)
    local mean_overall = r(mean)
    local sd_overall = r(sd)

    if `n_overall' == 0 {
        local overall_cell = "NA"
    }
    else if inrange(`n_overall', 1, 7) {
        local overall_cell = "[REDACTED]"
    }
    else {
        local overall_cell = strtrim(string(`mean_overall', "%9.1f")) + " (" + strtrim(string(`sd_overall', "%9.1f")) + ")"
    }

    /* DXA-eligible population */
    quietly summarize `varlist' if dxa_eligible == 1
    local n_high = r(N)
    local mean_high = r(mean)
    local sd_high = r(sd)

    if `n_high' == 0 {
        local high_cell = "NA"
    }
    else if inrange(`n_high', 1, 7) {
        local high_cell = "[REDACTED]"
    }
    else {
        local high_cell = strtrim(string(`mean_high', "%9.1f")) + " (" + strtrim(string(`sd_high', "%9.1f")) + ")"
    }

    post `handle' ///
        (`order') ///
        (`"`row_label'"') ///
        ("Mean (SD)") ///
        ("`overall_cell'") ///
        ("`high_cell'")
end

table1_continuous age, handle(`table1_post') order(`order')
local order = `order' + 1
table1_continuous qf_bmi, handle(`table1_post') order(`order')
local order = `order' + 1

/* Categorical-variable helper: n (%) */
capture program drop table1_category
program define table1_category
    syntax varname, Value(integer) Handle(name) Order(integer) NOverall(integer) NHigh(integer)
    local row_label : variable label `varlist'
    local category_label : label (`varlist') `value'

    quietly count if `varlist' == `value'
    local number_overall = r(N)

    sdc_n_pct, count(`number_overall') denominator(`noverall')
    local overall_cell `"`r(cell)'"'

    /* DXA-eligible population */
    quietly count if `varlist' == `value' & dxa_eligible == 1
    local number_high = r(N)

    sdc_n_pct, count(`number_high') denominator(`nhigh')
    local high_cell `"`r(cell)'"'

    post `handle' ///
        (`order') ///
        (`"`row_label'"') ///
        (`"`category_label'"') ///
        ("`overall_cell'") ///
        ("`high_cell'")
end

/* Sex */
forvalues value = 1/2 {
    table1_category table_sex, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        noverall(`N_overall') ///
        nhigh(`N_high')
    local order = `order' + 1
}

/* Ethnicity */
forvalues value = 1/9 {
    table1_category qf_ethrisk, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        noverall(`N_overall') ///
        nhigh(`N_high')
    local order = `order' + 1
}

/* Alcohol */
forvalues value = 0/5 {
    table1_category qf_alcohol_cat6, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        noverall(`N_overall') ///
        nhigh(`N_high')
    local order = `order' + 1
}

/* Smoking */
forvalues value = 0/4 {
    table1_category qf_smoke_cat, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        noverall(`N_overall') ///
        nhigh(`N_high')
    local order = `order' + 1
}

/* Binary QFracture predictors */
local binary_predictors ///
    b_anycancer ///
    b_asthmacopd ///
    b_carehome ///
    b_cvd ///
    b_dementia ///
    b_endocrine ///
    b_epilepsy2 ///
    b_falls ///
    b_liver ///
    b_malabsorption ///
    b_parkinsons ///
    b_ra_sle ///
    b_renal ///
    b_type1 ///
    b_type2 ///
    fh_osteoporosis ///
    b_antidepressant ///
    b_corticosteroids ///
    b_hrt_oest ///
    prev_hip ///
    prev_vertebral ///
    prev_wrist ///
    prev_humerus

foreach variable of local binary_predictors {
    local row_label : variable label `variable'    
    quietly count if `variable' == 1
    local number_overall = r(N)

    sdc_n_pct, count(`number_overall') denominator(`N_overall')
    local overall_cell `"`r(cell)'"'

    /* DXA-eligible population */
    quietly count if `variable' == 1 & dxa_eligible == 1
    local number_high = r(N)

    sdc_n_pct, count(`number_high') denominator(`N_high')
    local high_cell `"`r(cell)'"'

    post `table1_post' ///
        (`order') ///
        (`"`row_label'"') ///
        ("Yes") ///
        ("`overall_cell'") ///
        ("`high_cell'")
    local order = `order' + 1
}

postclose `table1_post'

use `table1_results', clear
sort row_order
drop row_order
export delimited using "output/table1_dxa.csv", replace

