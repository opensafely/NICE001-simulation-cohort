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

*Derived Table 1 variables
generate byte table_sex = .
replace table_sex = 1 if lower(strtrim(sex)) == "female"
replace table_sex = 2 if lower(strtrim(sex)) == "male"

label define table_sex_label 1 "Female" 2 "Male"
label values table_sex table_sex_label
label variable table_sex "Sex"

generate byte age_band = .
replace age_band = 1 if inrange(age, 50, 59)
replace age_band = 2 if inrange(age, 60, 69)
replace age_band = 3 if inrange(age, 70, 79)
replace age_band = 4 if inrange(age, 80, 89)
replace age_band = 5 if inrange(age, 90, 99)

label define age_band_label ///
    1 "50-59" ///
    2 "60-69" ///
    3 "70-79" ///
    4 "80-89" ///
    5 "90-99"
label values age_band age_band_label
label variable age_band "Age group, years"

generate byte prev_hip = dx_hip_fracture_n > 0 & !missing(dx_hip_fracture_n)
generate byte prev_vertebral = dx_vertebral_fracture_n > 0 & !missing(dx_vertebral_fracture_n)
generate byte prev_wrist = dx_wrist_fracture_n > 0 & !missing(dx_wrist_fracture_n)
generate byte prev_humerus = dx_proximal_humerus_fracture_n > 0 & !missing(dx_proximal_humerus_fracture_n)

generate byte prev_any = ///
    prev_hip == 1 | ///
    prev_vertebral == 1 | ///
    prev_wrist == 1 | ///
    prev_humerus == 1

label define prev_any_label 0 "No" 1 "Yes"
label values prev_any prev_any_label
label variable prev_any "Previous fracture"

label variable prev_hip "Previous hip fracture"
label variable prev_vertebral "Previous vertebral fracture"
label variable prev_wrist "Previous wrist fracture"
label variable prev_humerus "Previous proximal humerus fracture"

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
quietly count if dxa_eligible == 1
local N_dxa = r(N)

quietly count if dxa_eligible == 1 & table_sex == 1
local N_women = r(N)

quietly count if dxa_eligible == 1 & table_sex == 2
local N_men = r(N)

/* Create results dataset */
tempname table1_post
tempfile table1_results

postfile `table1_post' ///
    int row_order ///
    str70 characteristic ///
    str45 category ///
    str30 dxa_eligible_patients ///
    str30 women ///
    str30 men ///
    using `table1_results', replace
local order = 1

/* Population size */
sdc_count, count(`N_dxa')
local dxa_cell `"`r(cell)'"'

sdc_count, count(`N_women')
local women_cell `"`r(cell)'"'

sdc_count, count(`N_men')
local men_cell `"`r(cell)'"'

post `table1_post' ///
    (`order') ///
    ("Population, n") ///
    ("") ///
    ("`dxa_cell'") ///
    ("`women_cell'") ///
    ("`men_cell'")
local order = `order' + 1

/* Continuous-variable helper: mean (SD) */
capture program drop table1_continuous
program define table1_continuous
    syntax varname, Handle(name) Order(integer)
    local row_label : variable label `varlist'

    /* All DXA-eligible patients */
    quietly summarize `varlist' if dxa_eligible == 1
    local n_dxa = r(N)
    local mean_dxa = r(mean)
    local sd_dxa = r(sd)

    if `n_dxa' == 0 {
        local dxa_cell = "NA"
    }
    else if inrange(`n_dxa', 1, 7) {
        local dxa_cell = "[REDACTED]"
    }
    else {
        local dxa_cell = strtrim(string(`mean_dxa', "%9.1f")) + ///
            " (" + strtrim(string(`sd_dxa', "%9.1f")) + ")"
    }

    /* Women */
    quietly summarize `varlist' if dxa_eligible == 1 & table_sex == 1
    local n_women = r(N)
    local mean_women = r(mean)
    local sd_women = r(sd)

    if `n_women' == 0 {
        local women_cell = "NA"
    }
    else if inrange(`n_women', 1, 7) {
        local women_cell = "[REDACTED]"
    }
    else {
        local women_cell = strtrim(string(`mean_women', "%9.1f")) + ///
            " (" + strtrim(string(`sd_women', "%9.1f")) + ")"
    }

    /* Men */
    quietly summarize `varlist' if dxa_eligible == 1 & table_sex == 2
    local n_men = r(N)
    local mean_men = r(mean)
    local sd_men = r(sd)

    if `n_men' == 0 {
        local men_cell = "NA"
    }
    else if inrange(`n_men', 1, 7) {
        local men_cell = "[REDACTED]"
    }
    else {
        local men_cell = strtrim(string(`mean_men', "%9.1f")) + ///
            " (" + strtrim(string(`sd_men', "%9.1f")) + ")"
    }

    post `handle' ///
        (`order') ///
        (`"`row_label'"') ///
        ("Mean (SD)") ///
        ("`dxa_cell'") ///
        ("`women_cell'") ///
        ("`men_cell'")
end

/* Categorical-variable helper: n (%) */
capture program drop table1_category
program define table1_category
    syntax varname, Value(integer) Handle(name) Order(integer) ///
        Ndxa(integer) Nwomen(integer) Nmen(integer) ///
        [Characteristic(string) Category(string)]

    local row_label : variable label `varlist'
    local category_label : label (`varlist') `value'

    if `"`characteristic'"' != "" {
        local row_label `"`characteristic'"'
    }

    if `"`category'"' != "" {
        local category_label `"`category'"'
    }

    /* All DXA-eligible patients */
    quietly count if `varlist' == `value' & dxa_eligible == 1
    local number_dxa = r(N)

    sdc_n_pct, count(`number_dxa') denominator(`ndxa')
    local dxa_cell `"`r(cell)'"'

    /* Women */
    quietly count if `varlist' == `value' & dxa_eligible == 1 & table_sex == 1
    local number_women = r(N)

    sdc_n_pct, count(`number_women') denominator(`nwomen')
    local women_cell `"`r(cell)'"'

    /* Men */
    quietly count if `varlist' == `value' & dxa_eligible == 1 & table_sex == 2
    local number_men = r(N)

    sdc_n_pct, count(`number_men') denominator(`nmen')
    local men_cell `"`r(cell)'"'

    post `handle' ///
        (`order') ///
        (`"`row_label'"') ///
        (`"`category_label'"') ///
        ("`dxa_cell'") ///
        ("`women_cell'") ///
        ("`men_cell'")
end

/* Age: mean (SD), followed by age categories */
table1_continuous age, handle(`table1_post') order(`order')
local order = `order' + 1

forvalues value = 1/5 {
    table1_category age_band, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        ndxa(`N_dxa') ///
        nwomen(`N_women') ///
        nmen(`N_men')
    local order = `order' + 1
}

/* BMI: mean (SD) */
table1_continuous qf_bmi, handle(`table1_post') order(`order')
local order = `order' + 1

/* Previous fracture*/
foreach value in 1 0 {
    table1_category prev_any, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        ndxa(`N_dxa') ///
        nwomen(`N_women') ///
        nmen(`N_men')
    local order = `order' + 1
}

table1_category prev_hip, ///
    value(1) ///
    handle(`table1_post') ///
    order(`order') ///
    ndxa(`N_dxa') ///
    nwomen(`N_women') ///
    nmen(`N_men') ///
    characteristic("Previous fracture site") ///
    category("Hip")
local order = `order' + 1

table1_category prev_vertebral, ///
    value(1) ///
    handle(`table1_post') ///
    order(`order') ///
    ndxa(`N_dxa') ///
    nwomen(`N_women') ///
    nmen(`N_men') ///
    characteristic("Previous fracture site") ///
    category("Vertebral")
local order = `order' + 1

table1_category prev_wrist, ///
    value(1) ///
    handle(`table1_post') ///
    order(`order') ///
    ndxa(`N_dxa') ///
    nwomen(`N_women') ///
    nmen(`N_men') ///
    characteristic("Previous fracture site") ///
    category("Wrist")
local order = `order' + 1

table1_category prev_humerus, ///
    value(1) ///
    handle(`table1_post') ///
    order(`order') ///
    ndxa(`N_dxa') ///
    nwomen(`N_women') ///
    nmen(`N_men') ///
    characteristic("Previous fracture site") ///
    category("Proximal humerus")
local order = `order' + 1

/* Smoking status */
forvalues value = 0/4 {
    table1_category qf_smoke_cat, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        ndxa(`N_dxa') ///
        nwomen(`N_women') ///
        nmen(`N_men')
    local order = `order' + 1
}

/* Alcohol consumption */
forvalues value = 0/5 {
    table1_category qf_alcohol_cat6, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        ndxa(`N_dxa') ///
        nwomen(`N_women') ///
        nmen(`N_men')
    local order = `order' + 1
}

/* Ethnicity */
forvalues value = 1/9 {
    table1_category qf_ethrisk, ///
        value(`value') ///
        handle(`table1_post') ///
        order(`order') ///
        ndxa(`N_dxa') ///
        nwomen(`N_women') ///
        nmen(`N_men')
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
    table1_category `variable', ///
        value(1) ///
        handle(`table1_post') ///
        order(`order') ///
        ndxa(`N_dxa') ///
        nwomen(`N_women') ///
        nmen(`N_men') ///
        category("Yes")
    local order = `order' + 1
}

postclose `table1_post'

use `table1_results', clear
sort row_order
/* Display a characteristic label once*/
generate str70 characteristic_display = characteristic
replace characteristic_display = "" if _n > 1 & characteristic == characteristic[_n - 1]
drop characteristic
rename characteristic_display characteristic
order characteristic category dxa_eligible_patients women men
drop row_order
export delimited using "output/table1_dxa.csv", replace

