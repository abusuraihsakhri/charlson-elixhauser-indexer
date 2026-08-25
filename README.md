# Charlson & Elixhauser Comorbidity Indexer

A production-grade, zero-dependency Python implementation of the **Charlson Comorbidity Index (CCI)** (Quan et al. 2011 ICD-10 update of Charlson 1987) and the **Elixhauser Comorbidity Index (ECI)** (AHRQ 2021 Refined 31 categories) with **van Walraven (vW 2009) weighted composite scoring** and **10-year estimated survival calculation**.

---

## Clinical Overview & Methodological Background

Comorbidity measurement is essential in epidemiological research, health services outcome studies, risk adjustment, surgical risk prediction, and clinical risk stratification.

### 1. Charlson Comorbidity Index (CCI)
Developed by *Charlson et al. (J Chronic Dis 1987)* and updated for ICD-10 coding by *Quan et al. (Am J Epidemiol 2011)*.
- Evaluates **17 primary comorbid conditions**, each assigned a weight from 1 to 6.
- **Age Adjustment** (Charlson 1994):
  - Age $< 50$: +0
  - Age $50 - 59$: +1
  - Age $60 - 69$: +2
  - Age $70 - 79$: +3
  - Age $\ge 80$: +4
- **10-Year Estimated Survival Probability**:
  $$S(10) = 0.983^{\exp(0.9 \times \text{CCI}_{\text{age-adjusted}})}$$

#### Charlson Condition Weights & Hierarchy:
| Condition | Weight | Hierarchy Exclusion Rules |
| :--- | :---: | :--- |
| Myocardial Infarction | 1 | |
| Congestive Heart Failure | 1 | |
| Peripheral Vascular Disease | 1 | |
| Cerebrovascular Disease | 1 | |
| Dementia | 1 | |
| Chronic Pulmonary Disease | 1 | |
| Rheumatic Disease | 1 | |
| Peptic Ulcer Disease | 1 | |
| Mild Liver Disease | 1 | Suppressed if Moderate/Severe Liver Disease is present |
| Diabetes without complications | 1 | Suppressed if Diabetes with complications is present |
| Diabetes with complications | 2 | |
| Hemiplegia / Paraplegia | 2 | |
| Renal Disease | 2 | |
| Any Malignancy (non-metastatic) | 2 | Suppressed if Metastatic Solid Tumor is present |
| Moderate / Severe Liver Disease | 3 | |
| Metastatic Solid Tumor | 6 | |
| AIDS / HIV | 6 | |

---

### 2. Elixhauser Comorbidity Index (ECI) & van Walraven Score
Developed by *Elixhauser et al. (Med Care 1998)*, refined by AHRQ (2021), and weighted by *van Walraven et al. (Med Care 2009)* to produce a single composite integer index (-19 to +89) strongly predictive of in-hospital mortality.

- Evaluates **31 distinct comorbidity categories**.
- Includes positive weights (e.g. Metastatic Cancer: +12, Liver Disease: +11, Lymphoma: +9, Congestive Heart Failure: +7, Paralysis: +7) and protective/negative weights for specific administrative cohorts (e.g. Obesity: -4, Depression: -3, Drug abuse: -7).

---

## Project Structure

```
charlson-elixhauser-indexer/
├── charlson_elixhauser.py    # Core pure-Python comorbidity engine
├── cli.py                    # Interactive, single, and batch CSV CLI
├── test_charlson_elixhauser.py # Comprehensive unit test suite (27+ tests)
├── benchmark_dataset.json    # Verified clinical reference cases
├── sample_patients.csv       # Sample patient CSV with multi-diagnosis codes
├── Dockerfile                # Production container definition
├── docker-compose.yml        # Compose orchestration file
└── README.md                 # Technical documentation & clinical reference
```

---

## CLI Usage Guide

### 1. Interactive Clinical Consultation
Launch the interactive prompt to enter patient demographics and diagnosis codes:
```bash
python cli.py interactive
```

### 2. Single-Patient Scoring
```bash
# Score patient with MI, CHF, and complicated diabetes at age 68
python cli.py single --id "PT-100" --codes "I21.9; I50.9; E11.22" --age 68 --sex M

# Display detailed breakdown of all 17 Charlson & 31 Elixhauser criteria
python cli.py single --id "PT-100" --codes "I21.9; I50.9; E11.22" --age 68 --detail

# Output JSON for API and EHR integration
python cli.py single --id "PT-100" --codes "I21.9; I50.9; E11.22" --age 68 --json
```

### 3. Batch CSV Processing
Process an entire cohort CSV file:
```bash
python cli.py batch -i sample_patients.csv -o results.csv
```

---

## Programmatic Python API

```python
from charlson_elixhauser import assess_patient

res = assess_patient(
    patient_id="PT-2026",
    raw_codes="I21.9; I50.9; E11.65; N18.3; C34.9",
    age=72,
    sex="M",
)

print(f"Raw Charlson: {res.charlson_score}")
print(f"Age-Adjusted Charlson: {res.charlson_age_adjusted}")
print(f"10-Year Survival: {res.charlson_10yr_survival_pct:.2f}%")
print(f"Elixhauser Count: {res.elix_count}")
print(f"van Walraven Score: {res.elix_van_walraven}")
print(f"Mortality Risk: {res.mortality_risk_tier}")
```

---

## Unit Testing

Run the test suite with standard `unittest`:

```bash
python -m unittest discover -s . -p "test_*.py" -v
```

Test coverage includes:
- All 17 Charlson comorbidity categories and ICD-10 prefix normalization.
- Clinical hierarchy suppression (metastatic vs non-metastatic cancer, severe vs mild liver, complicated vs uncomplicated diabetes/hypertension).
- Charlson age-adjustment step functions and exponential 10-year survival formulas.
- Elixhauser 31 categories, total counts, and van Walraven weighted scores.
- Resilient CSV parsing across varying column headers and delimiters.

---

## References

1. **Charlson ME, et al.** (1987). *A new method of classifying prognostic comorbidity in longitudinal studies: development and validation*. J Chronic Dis; 40(5):373–383.
2. **Quan H, et al.** (2011). *Updating and validating the Charlson comorbidity index and score for risk adjustment in hospital administrative data using ICD-9 and ICD-10*. Am J Epidemiol; 173(6):676–682.
3. **Elixhauser A, et al.** (1998). *Comorbidity measures for use with administrative data*. Med Care; 36(1):8–27.
4. **van Walraven C, et al.** (2009). *A modification of the Elixhauser comorbidity measures for administrative data*. Med Care; 47(6):626–633.

---

## License

MIT License. Developed for clinical informatics and outcomes research.
