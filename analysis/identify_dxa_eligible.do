/*Identify the DXA-eligible population

The patient-level input is the source population already selected in ehrQL:
alive, aged 50-99, and continuously registered through the full baseline year. */

version 16.1
clear all
set more off

import delimited using "output/qfracture_variables.csv", clear varnames(1) stringcols(_all)
local numeric_variables ///
    age bmi_raw alcohol_numeric smoking_numeric ///
    rx_antidepressant_n_6m rx_anticonvulsant_n_6m ///
    rx_hrt_n_6m rx_corticosteroid_n_6m ///
    dx_hip_fracture_n dx_vertebral_fracture_n ///
    dx_wrist_fracture_n dx_proximal_humerus_fracture_n

foreach variable of local numeric_variables {
        destring `variable', replace
}
*

local boolean_variables ///
    rx_osteoporosis_b4 b_carehome carehome_nursing carehome_no_nursing ///
    b_anycancer b_asthmacopd b_cvd b_dementia b_endocrine ///
    dx_epilepsy b_falls b_liver b_malabsorption b_parkinsons ///
    b_ra_sle b_renal b_type1 b_type2 fh_osteoporosis ///
	dx_hip_fracture dx_hip_fracture_hos ///
	dx_vertebral_fracture dx_vertebral_fracture_hos ///
	dx_wrist_fracture dx_wrist_fracture_hos ///
	dx_proximal_humerus_fracture dx_proximal_humerus_fracture_hos 

foreach variable of local boolean_variables {
    tempvar parsed_boolean
    generate byte `parsed_boolean' = .
    replace `parsed_boolean' = 1 if inlist(upper(strtrim(`variable')), "T", "TRUE", "1")
    replace `parsed_boolean' = 0 if inlist(upper(strtrim(`variable')), "F", "FALSE", "0")
    assert !missing(`parsed_boolean') if strtrim(`variable') != ""
    drop `variable'
    rename `parsed_boolean' `variable'
    assert inlist(`variable', 0, 1) | missing(`variable')
}
*

*Construct Route 1 from fracture history
foreach variable in ///
    dx_hip_fracture_n ///
    dx_vertebral_fracture_n ///
    dx_wrist_fracture_n ///
    dx_proximal_humerus_fracture_n { 
    assert inrange(`variable', 0, 2)
}
generate byte route1a = inrange(dx_hip_fracture_n, 1, 2) | inrange(dx_vertebral_fracture_n, 1, 2)
generate byte route1b = (dx_wrist_fracture_n + dx_proximal_humerus_fracture_n >= 2)
generate byte route1 = (rx_osteoporosis_b4 == 0) & (route1a == 1 | route1b == 1)

*Derive medication-based QFracture predictors
generate byte b_antidepressant = rx_antidepressant_n_6m >= 2 if !missing(rx_antidepressant_n_6m)
generate byte b_corticosteroids = rx_corticosteroid_n_6m >= 2 if !missing(rx_corticosteroid_n_6m)
generate byte b_hrt_oest = rx_hrt_n_6m >= 2 if !missing(rx_hrt_n_6m)
generate byte b_epilepsy2 = (dx_epilepsy == 1 | rx_anticonvulsant_n_6m >= 2) if !missing(dx_epilepsy) & !missing(rx_anticonvulsant_n_6m)

*BMI
generate double qf_bmi = bmi_raw
replace qf_bmi = . if qf_bmi <= 0
replace qf_bmi = 20 if qf_bmi < 20 & qf_bmi <.
replace qf_bmi = 40 if qf_bmi > 40 & qf_bmi <.

*Alcohol
generate byte alcohol_no_record = strtrim(alcohol_code) == ""
generate byte qf_alcohol_cat6 = .
rename alcohol_code code 
merge m:1 code using "codelists/uploaded/user-xixiong-alcohol-with-class.dta", keep(1 3) nogen
rename code alcohol_code
rename term alcohol_term
generate byte alcohol_current_unknown = !alcohol_no_record & class == "Current_unknown"
replace qf_alcohol_cat6 = 0 if class == "None"
replace qf_alcohol_cat6 = 1 if class == "<1 unit/day"
replace qf_alcohol_cat6 = 2 if class == "1-2 units/day"
replace qf_alcohol_cat6 = 3 if class == "3-6 units/day"
replace qf_alcohol_cat6 = 4 if class == "7-9 units/day"
replace qf_alcohol_cat6 = 5 if class == ">9 units/day"
drop class

*Smoking 
generate byte smoking_no_record = strtrim(smoking_code) == ""
generate byte qf_smoke_cat = .
rename smoking_code code
merge m:1 code using "codelists/uploaded/user-xixiong-smoking-with-class.dta", keep(1 3) nogen
rename code smoking_code
rename term smoking_term
generate byte smoking_current_unknown = !smoking_no_record & class == 2
replace qf_smoke_cat = 0 if class == 0
replace qf_smoke_cat = 1 if class == 1
replace qf_smoke_cat = 2 if class == 3
replace qf_smoke_cat = 3 if class == 4
replace qf_smoke_cat = 4 if class == 5
drop class

/* Ethnicity mapping from SUS follows the QFracture array order:
      1 White, 2 Indian, 3 Pakistani, 4 Bangladeshi, 5 Other Asian,
      6 Black Caribbean, 7 Black African, 8 Chinese, 9 Other/mixed. */
generate byte qf_ethrisk = .
rename ethnicity_code code
merge m:1 code using "codelists/uploaded/user-xixiong-ethnicity-with-class.dta", keep(1 3) nogen
rename code ethnicity_code

replace qf_ethrisk = 1 if class == "White"
replace qf_ethrisk = 2 if class == "Indian"
replace qf_ethrisk = 3 if class == "Pakistani"
replace qf_ethrisk = 4 if class == "Bangladeshi"
replace qf_ethrisk = 5 if class == "Other Asian"
replace qf_ethrisk = 6 if class == "Black Caribbean"
replace qf_ethrisk = 7 if class == "Black African"
replace qf_ethrisk = 8 if class == "Chinese"
replace qf_ethrisk = 9 if class == "Other"
drop term class

*uses SUS as a fallback
replace qf_ethrisk = 1 if inlist(upper(strtrim(ethnicity_sus_code)), "A", "B", "C") & qf_ethrisk == .
replace qf_ethrisk = 2 if upper(strtrim(ethnicity_sus_code)) == "H" & qf_ethrisk == .
replace qf_ethrisk = 3 if upper(strtrim(ethnicity_sus_code)) == "J" & qf_ethrisk == .
replace qf_ethrisk = 4 if upper(strtrim(ethnicity_sus_code)) == "K" & qf_ethrisk == .
replace qf_ethrisk = 5 if upper(strtrim(ethnicity_sus_code)) == "L" & qf_ethrisk == .
replace qf_ethrisk = 6 if upper(strtrim(ethnicity_sus_code)) == "M" & qf_ethrisk == .
replace qf_ethrisk = 7 if upper(strtrim(ethnicity_sus_code)) == "N" & qf_ethrisk == .
replace qf_ethrisk = 8 if upper(strtrim(ethnicity_sus_code)) == "R" & qf_ethrisk == .
replace qf_ethrisk = 9 if inlist(upper(strtrim(ethnicity_sus_code)), "D", "E", "F", "G", "P", "S") & qf_ethrisk == .
replace qf_ethrisk = 1 if missing(qf_ethrisk)

*carehome
replace b_carehome = 1 if carehome_nursing == 1 | carehome_no_nursing == 1

*Calculate QFracture 2016 10-year MOF risk 
do "analysis/calculate_qfracture_mof.do"

*Construct Route 2 and final DXA eligibility
generate byte route2 = (rx_osteoporosis_b4 == 0 & qfracture_calculable == 1 & qfracture_mof_10y_pct >= 10)
generate byte dxa_eligible = (route1 == 1 | route2 == 1)


/*Display cohort flow and save patient-level outputs --------------- */
* Reasons why QFracture may not be calculable
generate byte qf_age_problem = missing(age) | !inrange(age, 50, 99)
generate byte qf_sex_problem = !inlist(lower(strtrim(sex)), "male", "female")
generate byte qf_bmi_problem = missing(qf_bmi)

quietly count
display "Source population: " r(N)

quietly count if rx_osteoporosis_b4 == 1
display "Excluded for baseline osteoporosis treatment: " r(N)

quietly count if rx_osteoporosis_b4 == 0
display "Treatment-free population: " r(N)

*QFracture calculability
quietly count if rx_osteoporosis_b4 == 0 & qfracture_calculable != 1
display "Excluded because baseline QFracture is not calculable: " r(N)

*Individual reasons: these counts are not mutually exclusive
quietly count if rx_osteoporosis_b4 == 0 & qf_age_problem == 1
display "  Missing or invalid age: " r(N)

quietly count if rx_osteoporosis_b4 == 0 & qf_sex_problem == 1
display "  Sex not male or female (mixed/unknown/missing): " r(N)

quietly count if rx_osteoporosis_b4 == 0 & qf_bmi_problem == 1
display "  Missing or invalid BMI: " r(N)

quietly count if rx_osteoporosis_b4 == 0 & alcohol_no_record == 1
display "  No alcohol record: " r(N)
quietly count if rx_osteoporosis_b4 == 0 & alcohol_current_unknown == 1
display "  Alcohol consumption level unknown: " r(N)

quietly count if rx_osteoporosis_b4 == 0 & smoking_no_record == 1
display "  No smoking record: " r(N)
quietly count if rx_osteoporosis_b4 == 0 & smoking_current_unknown == 1
display "  Smoking level unknown: " r(N)

*Identify patients with at least one listed reason *
generate byte qf_listed_missing_reason = ///
    qf_age_problem == 1 | ///
    qf_sex_problem == 1 | ///
    qf_bmi_problem == 1 | ///
    alcohol_no_record == 1 | ///
    alcohol_current_unknown == 1 | ///
    smoking_no_record == 1 | ///
    smoking_current_unknown == 1

quietly count if rx_osteoporosis_b4 == 0 & ///
    qfracture_calculable != 1 & ///
    qf_listed_missing_reason == 1
display "Non-calculable with at least one listed reason: " r(N)

* Final populations
quietly count if route1 == 1
display "DXA eligible through Route 1: " r(N)

quietly count if route2 == 1
display "DXA eligible through Route 2: " r(N)

quietly count if dxa_eligible == 1
display "Total DXA-eligible population: " r(N)

compress
save "output/qfracture_cohort_with_dxa_flags.dta", replace

preserve
keep if dxa_eligible == 1
save "output/dxa_eligible_population.dta", replace
restore






