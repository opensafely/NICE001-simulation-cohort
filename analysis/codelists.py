"""Codelists used by the QFracture 2016 extraction."""

from ehrql import codelist_from_csv

def _snomed(filename): return codelist_from_csv(f"codelists/{filename}", column="code")
def _dmd(filename): return codelist_from_csv(f"codelists/{filename}", column="code")
def _icd10(filename): return codelist_from_csv(f"codelists/{filename}", column="code")


# Measurements and behaviours -------------------------------------------------
bmi_codes = _snomed("nhsd-primary-care-domain-refsets-bmival_cod.csv")
alcohol_codes = codelist_from_csv("codelists/uploaded/user-xixiong-alcohol-without-class.csv", column="code") 
ethnicity_codes = codelist_from_csv("codelists/uploaded/user-xixiong-ethnicity-without-class.csv", column="code") 
smoking_codes = codelist_from_csv("codelists/uploaded/user-xixiong-smoking-without-class.csv", column="code") 

# Prescribing -----------------------------------------------------------------
# Osteoporosis treatment: oral and intravenous bisphosphonates, denosumab, raloxifene, teriparatide and and romosozumab. 
alendronic_acid_dmd = _dmd("opensafely-alendronic-acid.csv")
risedronate_dmd = _dmd("opensafely-risedronate.csv")
ibandronic_acid_dmd = _dmd("opensafely-ibandronic-acid.csv")
zoledronic_acid_dmd = _dmd("opensafely-zoledronic-acid.csv")
teriparatide_dmd = _dmd("opensafely-teriparatide.csv")
denosumab_osteoporosis_dmd = codelist_from_csv("codelists/uploaded/user-xixiong-denosumab-for-osteoporosis.csv", column="code") 
raloxifene_dmd = codelist_from_csv("codelists/uploaded/user-xixiong-raloxifene.csv", column="code") 
romosozumab_dmd = codelist_from_csv("codelists/uploaded/user-xixiong-romosozumab.csv", column="code") 

osteoporosis_treatment_dmd = (alendronic_acid_dmd + risedronate_dmd
                            + ibandronic_acid_dmd + zoledronic_acid_dmd
                            + denosumab_osteoporosis_dmd + raloxifene_dmd
                            + teriparatide_dmd + romosozumab_dmd)

antidepressant_dmd = _dmd("nhs-drug-refsets-antidepdrug_cod.csv")
anticonvulsant_dmd = _dmd("nhs-drug-refsets-epildrug_cod.csv")
systemic_corticosteroid_dmd = codelist_from_csv("codelists/uploaded/user-xixiong-oral-systemic-corticosteroids.csv", column="code") 
hrt_oestrogen_dmd = codelist_from_csv("codelists/uploaded/user-xixiong-oestrogen-only-hrt.csv", column="code") 

# Diagnoses -------------------------------------------------------------------
any_cancer_codes = _snomed("user-ciaranmci-cancer-snomed-ct.csv")
asthma_codes = _snomed("nhsd-primary-care-domain-refsets-ast_cod.csv")
copd_codes = _snomed("nhsd-primary-care-domain-refsets-copd_cod.csv")

mi_codes = _snomed("nhsd-primary-care-domain-refsets-mi_cod.csv")
angina_codes = _snomed("nhsd-primary-care-domain-refsets-angina_cod.csv")
stroke_codes = _snomed("nhsd-primary-care-domain-refsets-strk_cod.csv")
tia_codes = _snomed("nhsd-primary-care-domain-refsets-tia_cod.csv")

dementia_codes = _snomed("nhsd-primary-care-domain-refsets-dem_cod.csv")

hyperparathyroidism_codes = _snomed("opensafely-hyperparathyroidism.csv")
thyrotoxicosis_codes = codelist_from_csv("codelists/uploaded/user-xixiong-thyrotoxicosis.csv", column="code") 
cushing_syndrome_codes = codelist_from_csv("codelists/uploaded/user-xixiong-cushing-syndrome.csv", column="code") 

epilepsy_codes = _snomed("nhsd-primary-care-domain-refsets-epil_cod.csv")
falls_codes = _snomed("nhsd-primary-care-domain-refsets-falls_cod.csv")
chronic_liver_disease_codes = _snomed("nhsd-primary-care-domain-refsets-cldatrisk1_cod.csv")

crohns_codes = _snomed("nhsd-primary-care-domain-refsets-crohns_cod.csv")
ulcerative_colitis_codes = _snomed("nhsd-primary-care-domain-refsets-ulccolitis_cod.csv")
coeliac_codes = _snomed("nhsd-primary-care-domain-refsets-coeliac_cod.csv")
other_malabsorption_codes = _snomed("opensafely-malabsorption-due-to-other-causes.csv")

parkinsons_codes = _snomed("nhsd-primary-care-domain-refsets-pd_cod.csv")

rheumatoid_arthritis_codes = _snomed("nhsd-primary-care-domain-refsets-rarth_cod.csv")
sle_codes = _snomed("nhsd-primary-care-domain-refsets-slupus_cod.csv")

renal_disease_codes = _snomed("nhsd-primary-care-domain-refsets-ckdatrisk1_cod.csv")
type1_diabetes_codes = _snomed("nhsd-primary-care-domain-refsets-dmtype1_cod.csv")
type2_diabetes_codes = _snomed("nhsd-primary-care-domain-refsets-dmtype2_cod.csv")

family_history_osteoporosis_codes = codelist_from_csv("codelists/uploaded/user-xixiong-parental-history-of-osteoporosis-or-hip-fracture.csv", column="code") 

# Fracture sites ---------------------------------------------------------------
hip_fracture_snomed = codelist_from_csv("codelists/uploaded/user-xixiong-hip-fracture.csv", column="code") 
vertebral_fracture_snomed = codelist_from_csv("codelists/uploaded/user-xixiong-vertebral-fracture.csv", column="code") 
wrist_fracture_snomed = codelist_from_csv("codelists/uploaded/user-xixiong-wrist-fracture.csv", column="code") 
proximal_humerus_fracture_snomed = codelist_from_csv("codelists/uploaded/user-xixiong-proximal-humerus-fracture.csv", column="code") 

hip_fracture_icd10 = codelist_from_csv("codelists/uploaded/user-xixiong-hip-fracture-icd10.csv", column="code") 
vertebral_fracture_icd10 = codelist_from_csv("codelists/uploaded/user-xixiong-vertebral-fracture-icd10.csv", column="code") 
wrist_fracture_icd10 = codelist_from_csv("codelists/uploaded/user-xixiong-wrist-fracture-icd10.csv", column="code") 
proximal_humerus_fracture_icd10 = codelist_from_csv("codelists/uploaded/user-xixiong-proximal-humerus-fracture-icd10.csv", column="code") 
osteoporotic_fracture_icd10 = codelist_from_csv("codelists/uploaded/user-xixiong-osteoporotic-fracture-icd10.csv", column="code") 

# Composites ------------------------------------------------------------------
asthma_or_copd_codes = asthma_codes + copd_codes
cvd_codes = mi_codes + angina_codes + stroke_codes + tia_codes
endocrine_codes = (hyperparathyroidism_codes + thyrotoxicosis_codes + cushing_syndrome_codes)
malabsorption_codes = (crohns_codes + ulcerative_colitis_codes + coeliac_codes + other_malabsorption_codes)
ra_or_sle_codes = rheumatoid_arthritis_codes + sle_codes
