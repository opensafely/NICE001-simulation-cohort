"""ehrQL extraction of the source population and variables required for
QFracture 2016 and DXA-eligibility assessment.

Population:
* alive and aged 50-99 inclusive on 1 January 2025;
* continuously registered at a TPP practice using SystmOne from 
  1 January 2024 through 1 January 2025.

Recorded primary-care osteoporosis-treatment prescriptions during the
12 months before the index date are retained as a baseline flag. Treatment
exclusion, QFracture calculation and final DXA eligibility are applied
downstream in Stata.
"""

from datetime import date
from ehrql import create_dataset
from ehrql.tables.tpp import patients, practice_registrations
import codelists
from qfracture_variables import build_qfracture_variables
from variables_helpers import medication_records


INDEX_DATE = date(2025, 1, 1)
BASELINE_START_DATE = date(2024, 1, 1)
TREATMENT_WASHOUT_START_DATE = date(2024, 1, 1)
MINIMUM_AGE = 50
MAXIMUM_AGE = 99

dataset = create_dataset()
dataset.configure_dummy_data(population_size=1000)

age_at_index = patients.age_on(INDEX_DATE)
registered_for_full_baseline = practice_registrations.spanning_with_systmone(
    BASELINE_START_DATE, INDEX_DATE
).exists_for_patient()

baseline_osteoporosis_treatment_records = medication_records(
    codelists.osteoporosis_treatment_dmd,
    TREATMENT_WASHOUT_START_DATE,
    INDEX_DATE,
)

osteoporosis_treatment_baseline = (
    baseline_osteoporosis_treatment_records.exists_for_patient()
)

dataset.define_population(
    patients.is_alive_on(INDEX_DATE)
    & (age_at_index >= MINIMUM_AGE)
    & (age_at_index <= MAXIMUM_AGE)
    & registered_for_full_baseline
)

dataset.add_column("rx_osteoporosis_b4", osteoporosis_treatment_baseline,)

for variable_name, patient_series in build_qfracture_variables(INDEX_DATE).items():
    dataset.add_column(variable_name, patient_series)
