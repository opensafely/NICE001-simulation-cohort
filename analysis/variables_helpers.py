"""Reusable ehrQL helpers for patient-level QFracture measurements."""

from ehrql import case, days, when, minimum_of
from ehrql.tables.tpp import apcs, clinical_events, medications


def clinical_events_before(codelist, index_date):
    return clinical_events.where(clinical_events.snomedct_code.is_in(codelist)).where(clinical_events.date < index_date)


def latest_clinical_event(codelist, index_date):
    events = clinical_events_before(codelist, index_date)
    return events.sort_by(clinical_events.date).last_for_patient()


def latest_numeric_clinical_event(codelist, index_date):
    events = clinical_events_before(codelist, index_date).where(clinical_events.numeric_value.is_not_null())
    return events.sort_by(clinical_events.date).last_for_patient()


def ever_recorded(codelist, index_date):
    return clinical_events_before(codelist, index_date).exists_for_patient()


def fracture_episode_count_capped_at_two(
    snomed_codelist,
    icd10_codelist,
    index_date,
    episode_window_days=180,
):
    """Return 0, 1 or 2 same-site fracture episodes before index.
    The earliest matching primary-care or hospital record defines the first
    fracture episode. A second episode is identified if any subsequent record
    for the same fracture site occurs more than 180 days after the earliest
    record."""

    primary_care_events = clinical_events_before(snomed_codelist, index_date)
    hospital_events = hospital_admissions_before(icd10_codelist, index_date)

    first_primary_care_date = (primary_care_events.sort_by(clinical_events.date)
        .first_for_patient()
        .date
    )
    first_hospital_date = (hospital_events.sort_by(apcs.admission_date)
        .first_for_patient()
        .admission_date
    )

    first_episode_date = minimum_of(first_primary_care_date, first_hospital_date,)
    second_episode_threshold = first_episode_date + days(episode_window_days)

    later_primary_care_record = primary_care_events.where(
        clinical_events.date > second_episode_threshold
    ).exists_for_patient()

    later_hospital_record = hospital_events.where(
        apcs.admission_date > second_episode_threshold
    ).exists_for_patient()

    return case(
        when(later_primary_care_record | later_hospital_record).then(2),
        when(first_episode_date.is_not_null()).then(1),
        otherwise=0,
    )


def medication_records(codelist, start_date, end_date):
    return medications.where(medications.dmd_code.is_in(codelist)).where(
        (medications.date >= start_date) & (medications.date < end_date)
    )


def medication_count(codelist, start_date, end_date):
    return medication_records(codelist, start_date, end_date).count_for_patient()


def hospital_admissions_before(icd10_codelist, index_date):
    return apcs.where(apcs.admission_date < index_date).where(
        apcs.all_diagnoses.contains_any_of(icd10_codelist)
    )


def hospital_diagnosis_history(icd10_codelist, index_date):
    return hospital_admissions_before(
        icd10_codelist, index_date
    ).exists_for_patient()
