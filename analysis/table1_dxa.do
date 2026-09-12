/* Table 1: QFracture predictor characteristics ------------------------- */
version 16.1
clear all
set more off

use "output/qfracture_cohort_with_dxa_flags.dta", clear
preserve

generate byte table_previous_fracture = .
replace table_previous_fracture = 1 if inrange(dx_hip_fracture_n, 1, 2) 
replace table_previous_fracture = 2 if inrange(dx_vertebral_fracture_n, 1, 2) 
replace table_previous_fracture = 3 if inrange(dx_wrist_fracture_n, 1, 2) 
replace table_previous_fracture = 4 if inrange(dx_proximal_humerus_fracture_n, 1, 2)

label define table_previous_fracture_label ///
	1 "Hip fracture" ///
	2 "Vertebral fracture" ///
	3 "Wrist fracture" ///
	4 "Proximal humerus fracture" 
label values table_previous_fracture table_previous_fracture_label
label variable table_previous_fracture "Previous fracture"

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
local overall_cell = string(`N_overall', "%12.0fc")
local high_cell = string(`N_high', "%12.0fc")

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
    local overall_cell = string(r(mean), "%9.1f") + " (" + string(r(sd), "%9.1f") + ")"
    quietly summarize `varlist' if dxa_eligible == 1
    local high_cell = string(r(mean), "%9.1f") + " (" + string(r(sd), "%9.1f") + ")"

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
    local overall_cell = string(`number_overall', "%12.0fc") + " (" + ///
        string(100 * `number_overall' / `noverall', "%5.1f") + "%)"

    quietly count if `varlist' == `value' & dxa_eligible == 1
    local number_high = r(N)
    if `nhigh' > 0 {
        local high_cell = string(`number_high', "%12.0fc") + " (" + ///
            string(100 * `number_high' / `nhigh', "%5.1f") + "%)"
    }
    else {
        local high_cell = "0 (NA)"
    }

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
    b_hrt_oest 

foreach variable of local binary_predictors {
    quietly count if `variable' == 1
    local number_overall = r(N)

    local overall_cell = string(`number_overall', "%12.0fc") + " (" + ///
        string(100 * `number_overall' / `N_overall', "%5.1f") + "%)"

    quietly count if `variable' == 1 & dxa_eligible == 1
    local number_high = r(N)
    if `N_high' > 0 {
        local high_cell = string(`number_high', "%12.0fc") + " (" + ///
            string(100 * `number_high' / `N_high', "%5.1f") + "%)"
    }
    else {
        local high_cell = "0 (NA)"
    }
    
	local row_label : variable label `variable'
    post `table1_post' ///
        (`order') ///
        (`"`row_label'"') ///
        ("Yes") ///
        ("`overall_cell'") ///
        ("`high_cell'")
    local order = `order' + 1
}

/* Previous fracture */
forvalues value = 1/4 {
    table1_category table_previous_fracture, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        noverall(`N_overall') ///
        nhigh(`N_high')
    local order = `order' + 1
}

postclose `table1_post'

/* Export Table 1 */
use `table1_results', clear
sort row_order
drop row_order

export delimited using "output/table1_dxa.csv", replace
restore
