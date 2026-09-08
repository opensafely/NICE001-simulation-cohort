"""Patient-level raw measurements needed to construct QFracture 2016 inputs.

This module deliberately extracts raw codes, values, component diagnoses and
prescription counts. 
"""

from ehrql import months
from ehrql.tables.tpp import addresses, ethnicity_from_sus, patients

import codelists
from variables_helpers import (
    ever_recorded,
    fracture_episode_count_capped_at_two,
    hospital_diagnosis_history,
    latest_clinical_event,
    latest_numeric_clinical_event,
    medication_count,
)


RECENT_MEDICATION_LOOKBACK_MONTHS = 6
FRACTURE_EPISODE_WINDOW_DAYS = 180


def _fracture_site_variables(name, snomed_codes, icd10_codes, index_date):
    return {
        f"dx_{name}_fracture": ever_recorded(
            snomed_codes, index_date
        ),
        f"dx_{name}_fracture_hos": hospital_diagnosis_history(
            icd10_codes, index_date
        ),
        f"dx_{name}_fracture_n": (
            fracture_episode_count_capped_at_two(
                snomed_codes,
                icd10_codes,
                index_date,
                episode_window_days=FRACTURE_EPISODE_WINDOW_DAYS,
            )
        ),
    }


def build_qfracture_variables(index_date):
    recent_medication_start = index_date - months(RECENT_MEDICATION_LOOKBACK_MONTHS)

    latest_bmi = latest_numeric_clinical_event(codelists.bmi_codes, index_date)
    latest_alcohol = latest_clinical_event(codelists.alcohol_codes, index_date)
    latest_alcohol_numeric = latest_numeric_clinical_event(codelists.alcohol_codes, index_date)
    latest_smoking = latest_clinical_event(codelists.smoking_codes, index_date)
    latest_smoking_numeric = latest_numeric_clinical_event(codelists.smoking_codes, index_date)
    latest_ethnicity = latest_clinical_event(codelists.ethnicity_codes, index_date)

    address_at_index = addresses.for_patient_on(index_date)

    antidepressant_count = medication_count(codelists.antidepressant_dmd, recent_medication_start, index_date,)
    anticonvulsant_count = medication_count(codelists.anticonvulsant_dmd, recent_medication_start, index_date,)
    hrt_oestrogen_count = medication_count(codelists.hrt_oestrogen_dmd, recent_medication_start, index_date,)
    corticosteroid_count = medication_count(codelists.systemic_corticosteroid_dmd, recent_medication_start, index_date,)

    variables = {
        # Demography and continuous measurements
        "sex": patients.sex,
        "age": patients.age_on(index_date),
        "bmi_raw": latest_bmi.numeric_value,
        "bmi_date": latest_bmi.date,

        # Raw behaviour records. 
        "alcohol_code": latest_alcohol.snomedct_code,
        "alcohol_date": latest_alcohol.date,
        "alcohol_numeric_code": latest_alcohol_numeric.snomedct_code, 
        "alcohol_numeric": latest_alcohol_numeric.numeric_value,
        "alcohol_numeric_date": latest_alcohol_numeric.date,
        "smoking_code": latest_smoking.snomedct_code,
        "smoking_date": latest_smoking.date,
        "smoking_numeric_code": latest_smoking_numeric.snomedct_code,
        "smoking_numeric": latest_smoking_numeric.numeric_value,
        "smoking_numeric_date": latest_smoking_numeric.date,

        # Preserve both sources. The 10-level QFracture mapping is downstream.
        "ethnicity_code": latest_ethnicity.snomedct_code,
        "ethnicity_date": latest_ethnicity.date,
        "ethnicity_sus_code": ethnicity_from_sus.code,

        # Medication counts 
        "rx_antidepressant_n_6m": antidepressant_count,
        "rx_anticonvulsant_n_6m": anticonvulsant_count,
        "rx_hrt_n_6m": hrt_oestrogen_count,
        "rx_corticosteroid_n_6m": corticosteroid_count,

        # Care-home components at index
        "b_carehome": address_at_index.care_home_is_potential_match.when_null_then(False),
        "carehome_nursing": (address_at_index.care_home_requires_nursing.when_null_then(False)),
        "carehome_no_nursing": (address_at_index.care_home_does_not_require_nursing.when_null_then(False)),

        # Diagnoses/history.
        "b_anycancer": ever_recorded(codelists.any_cancer_codes, index_date),
        "b_asthmacopd": ever_recorded(codelists.asthma_or_copd_codes, index_date),
        "b_cvd": ever_recorded(codelists.cvd_codes, index_date),
        "b_dementia": ever_recorded(codelists.dementia_codes, index_date),
        "b_endocrine": ever_recorded(codelists.endocrine_codes, index_date),
        "dx_epilepsy": ever_recorded(codelists.epilepsy_codes, index_date),
        "b_falls": ever_recorded(codelists.falls_codes, index_date),
        "b_liver": ever_recorded(codelists.chronic_liver_disease_codes, index_date),
        "b_malabsorption": ever_recorded(codelists.malabsorption_codes, index_date),
        "b_parkinsons": ever_recorded(codelists.parkinsons_codes, index_date),
        "b_ra_sle": ever_recorded(codelists.ra_or_sle_codes, index_date),
        "b_renal": ever_recorded(codelists.renal_disease_codes, index_date),
        "b_type1": ever_recorded(codelists.type1_diabetes_codes, index_date),
        "b_type2": ever_recorded(codelists.type2_diabetes_codes, index_date),
        "fh_osteoporosis": ever_recorded(codelists.family_history_osteoporosis_codes, index_date),
    }

    variables.update(_fracture_site_variables(
            "hip",
            codelists.hip_fracture_snomed,
            codelists.hip_fracture_icd10,
            index_date,
        )
    )
    variables.update(_fracture_site_variables(
            "vertebral",
            codelists.vertebral_fracture_snomed,
            codelists.vertebral_fracture_icd10,
            index_date,
        )
    )
    variables.update(_fracture_site_variables(
            "wrist",
            codelists.wrist_fracture_snomed,
            codelists.wrist_fracture_icd10,
            index_date,
        )
    )
    variables.update(_fracture_site_variables(
            "proximal_humerus",
            codelists.proximal_humerus_fracture_snomed,
            codelists.proximal_humerus_fracture_icd10,
            index_date,
        )
    )

    return variables
