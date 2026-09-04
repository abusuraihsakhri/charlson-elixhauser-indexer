# Charlson Comorbidity Index & Elixhauser (van Walraven) Indexer

> **Domain:** Health Services Research, Epidemiology, Risk Adjustment & Biostatistics  
> **Clinical Guidelines & Standards:** Charlson et al. (J Chronic Dis 1987), Quan et al. ICD-10 Coding Algorithm (Med Care 2005, Am J Epidemiol 2011), Elixhauser et al. (Med Care 1998), van Walraven et al. (Med Care 2009), AHRQ Elixhauser Comorbidity Software Refined v2021+

---

## 📖 Clinical & Methodological Overview

The **Charlson & Elixhauser Comorbidity Indexer** extracts and maps secondary ICD-10 diagnosis codes to validated clinical comorbidity frameworks:

1. **Charlson Comorbidity Index (CCI)**: 17 weighted comorbidity categories (e.g., myocardial infarction, metastatic solid tumor, AIDS/HIV, moderate-to-severe liver disease) with age-adjustment and 10-year actuarial survival estimation ($S(10) = 0.983^{\exp(\text{CCI} \times 0.9)} \times 100\%$).
2. **Elixhauser Comorbidity Measures**: 31 distinct chronic condition flags designed for administrative inpatient data.
3. **van Walraven Composite Mortality Weighting**: Empirically derived integer weights (-19 to +89) predicting in-hospital and 30-day mortality.

### Charlson Weight Categories (Quan et al. ICD-10 Mapping)

| Weight | Clinical Conditions |
|:---|:---|
| **1 pt** | Myocardial infarction, Congestive heart failure, Peripheral vascular disease, Cerebrovascular disease, Dementia, Chronic pulmonary disease, Rheumatic disease, Peptic ulcer disease, Mild liver disease, Diabetes without chronic complications |
| **2 pts** | Hemiplegia or paraplegia, Renal disease, Diabetes with chronic complications, Any malignancy (lymphoma, leukemia, solid tumor) |
| **3 pts** | Moderate or severe liver disease |
| **6 pts** | Metastatic solid tumor, AIDS/HIV |

*Age Adjustment:* +1 point per decade above 40 years (50–59: +1, 60–69: +2, 70–79: +3, $\ge 80$: +4).

---

## 💻 CLI Quickstart & Usage

### 1. Score an Individual Patient ICD-10 Profile
```bash
python cli.py eval --patient-id P001 --age 68 --sex M --codes "I21.9, E11.65, I50.9, N18.3, J44.9"
```

### 2. Interactive ICD-10 Questionnaire
```bash
python cli.py interactive
```

### 3. Batch Process Patient Diagnosis Cohort
```bash
python cli.py batch -i sample.csv -o out_results.csv
```

---

## 🧪 Verification & Testing

Execute comprehensive test suite via pytest:
```bash
python -m pytest -p no:zarr
```
