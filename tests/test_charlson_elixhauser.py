"""
Unit Test Suite for Charlson & Elixhauser Comorbidity Indexer.
Tests ICD-10 code parsing, Quan 2011 Charlson 17 categories, clinical hierarchy rules,
age-adjusted mortality risk, Elixhauser 31 categories, van Walraven index, and batch workflows.
"""

import csv
import io
import json
import os
import tempfile
import unittest
from unittest.mock import patch

from charlson_elixhauser import (
    assess_patient,
    charlson_10yr_survival,
    charlson_age_adjusted,
    charlson_flags,
    charlson_score,
    classify_diabetes,
    elixhauser_count,
    elixhauser_flags,
    elixhauser_van_walraven,
    normalize_icd,
    process_csv,
)
from cli import main as cli_main


class TestICDNormalization(unittest.TestCase):
    def test_normalize_icd(self):
        self.assertEqual(normalize_icd("i21.9"), "I219")
        self.assertEqual(normalize_icd(" E11.65 "), "E1165")
        self.assertEqual(normalize_icd("c-34.9"), "C349")
        self.assertEqual(normalize_icd(""), "")


class TestCharlsonCategories(unittest.TestCase):
    def test_myocardial_infarction(self):
        flags = charlson_flags(["I21.9"])
        self.assertTrue(flags["Myocardial infarction"])
        self.assertEqual(charlson_score(flags), 1)

    def test_congestive_heart_failure(self):
        flags = charlson_flags(["I50.9"])
        self.assertTrue(flags["Congestive heart failure"])
        self.assertEqual(charlson_score(flags), 1)

    def test_peripheral_vascular_disease(self):
        flags = charlson_flags(["I70.2"])
        self.assertTrue(flags["Peripheral vascular disease"])
        self.assertEqual(charlson_score(flags), 1)

    def test_cerebrovascular_disease(self):
        flags = charlson_flags(["I63.9"])
        self.assertTrue(flags["Cerebrovascular disease"])
        self.assertEqual(charlson_score(flags), 1)

    def test_dementia(self):
        flags = charlson_flags(["G30.9"])
        self.assertTrue(flags["Dementia"])
        self.assertEqual(charlson_score(flags), 1)

    def test_chronic_pulmonary_disease(self):
        flags = charlson_flags(["J44.1"])
        self.assertTrue(flags["Chronic pulmonary disease"])
        self.assertEqual(charlson_score(flags), 1)

    def test_rheumatic_disease(self):
        flags = charlson_flags(["M05.8"])
        self.assertTrue(flags["Rheumatic disease"])
        self.assertEqual(charlson_score(flags), 1)

    def test_peptic_ulcer_disease(self):
        flags = charlson_flags(["K25.0"])
        self.assertTrue(flags["Peptic ulcer disease"])
        self.assertEqual(charlson_score(flags), 1)

    def test_mild_liver_disease(self):
        flags = charlson_flags(["K70.3"])
        self.assertTrue(flags["Mild liver disease"])
        self.assertEqual(charlson_score(flags), 1)

    def test_diabetes_uncomplicated(self):
        for code in ("E11.9", "E11.65"):
            flags = charlson_flags([code])
            self.assertTrue(flags["Diabetes without complications"], code)
            self.assertFalse(flags["Diabetes with complications"], code)
            self.assertEqual(charlson_score(flags), 1, code)

    def test_diabetes_complicated(self):
        flags = charlson_flags(["E11.21"])
        self.assertFalse(flags["Diabetes without complications"])
        self.assertTrue(flags["Diabetes with complications"])
        self.assertEqual(charlson_score(flags), 2)

    def test_diabetes_hierarchy(self):
        # If both uncomplicated (E11.9) and complicated (E11.21) exist, only complicated is scored
        wo, wc = classify_diabetes(["E119", "E1121"])
        self.assertFalse(wo)
        self.assertTrue(wc)
        flags = charlson_flags(["E11.9", "E11.21"])
        self.assertEqual(charlson_score(flags), 2)

    def test_hemiplegia_paraplegia(self):
        flags = charlson_flags(["G81.9"])
        self.assertTrue(flags["Hemiplegia/Paraplegia"])
        self.assertEqual(charlson_score(flags), 2)

    def test_renal_disease(self):
        flags = charlson_flags(["N18.4"])
        self.assertTrue(flags["Renal disease"])
        self.assertEqual(charlson_score(flags), 2)

    def test_solid_tumor_malignancy(self):
        flags = charlson_flags(["C34.9"])
        self.assertTrue(flags["Any malignancy (non-metastatic)"])
        self.assertEqual(charlson_score(flags), 2)

    def test_moderate_severe_liver_and_hierarchy(self):
        # Severe liver (I85.0) scores 3 and cancels mild liver (K70.3)
        flags = charlson_flags(["I85.0", "K70.3"])
        self.assertTrue(flags["Moderate/severe liver disease"])
        self.assertFalse(flags["Mild liver disease"])
        self.assertEqual(charlson_score(flags), 3)

    def test_metastatic_tumor_and_hierarchy(self):
        # Metastatic cancer (C77.0) scores 6 and cancels non-metastatic tumor (C34.9)
        flags = charlson_flags(["C34.9", "C77.0"])
        self.assertTrue(flags["Metastatic solid tumor"])
        self.assertFalse(flags["Any malignancy (non-metastatic)"])
        self.assertEqual(charlson_score(flags), 6)

    def test_aids_hiv(self):
        flags = charlson_flags(["B20"])
        self.assertTrue(flags["AIDS/HIV"])
        self.assertEqual(charlson_score(flags), 6)


class TestCharlsonAgeAdjustmentAndSurvival(unittest.TestCase):
    def test_age_adjustment_tiers(self):
        score = 2
        self.assertEqual(charlson_age_adjusted(score, 45), 2)  # +0 (<50)
        self.assertEqual(charlson_age_adjusted(score, 55), 3)  # +1 (50-59)
        self.assertEqual(charlson_age_adjusted(score, 65), 4)  # +2 (60-69)
        self.assertEqual(charlson_age_adjusted(score, 75), 5)  # +3 (70-79)
        self.assertEqual(charlson_age_adjusted(score, 85), 6)  # +4 (>=80)
        self.assertIsNone(charlson_age_adjusted(score, None))

    def test_charlson_10yr_survival(self):
        # Score 0: 0.983 ^ exp(0) = 0.983 (98.3%)
        surv0 = charlson_10yr_survival(0)
        self.assertAlmostEqual(surv0, 0.983, places=3)

        # Higher score decreases survival probability monotonically
        surv2 = charlson_10yr_survival(2)
        surv5 = charlson_10yr_survival(5)
        self.assertTrue(1.0 > surv0 > surv2 > surv5 > 0.0)


class TestElixhauserAndVanWalraven(unittest.TestCase):
    def test_elixhauser_flags_and_count(self):
        codes = ["I50.9", "E66.0", "C77.0"]  # CHF (+7), Obesity (-4), Metastatic (+12)
        flags = elixhauser_flags(codes)
        self.assertTrue(flags["Congestive heart failure"])
        self.assertTrue(flags["Obesity"])
        self.assertTrue(flags["Metastatic cancer"])
        self.assertFalse(flags["Solid tumor without metastasis"])  # canceled by metastatic

        count = elixhauser_count(flags)
        self.assertEqual(count, 3)

        vw = elixhauser_van_walraven(flags)
        self.assertEqual(vw, 7 + (-4) + 12)  # 15

    def test_elixhauser_hypertension_hierarchy(self):
        codes = ["I10", "I12.9"]  # Uncomplicated + Complicated HTN
        flags = elixhauser_flags(codes)
        self.assertFalse(flags["Hypertension uncomplicated"])
        self.assertTrue(flags["Hypertension complicated"])


class TestPatientAssessmentAndBatch(unittest.TestCase):
    def test_assess_patient_full(self):
        res = assess_patient(
            patient_id="PT-99",
            raw_codes="I21.9; I50.9; E11.22; C34.9",
            age=68,
            sex="M",
        )
        self.assertEqual(res.patient_id, "PT-99")
        self.assertEqual(res.age, 68)
        self.assertEqual(res.n_codes, 4)
        # CCI: MI (1) + CHF (1) + DiabComp (2) + Malignancy (2) = 6
        self.assertEqual(res.charlson_score, 6)
        # Age-adjusted: 6 + 2 (age 68) = 8
        self.assertEqual(res.charlson_age_adjusted, 8)
        self.assertIsNotNone(res.charlson_10yr_survival_pct)
        self.assertEqual(res.mortality_risk_tier, "VERY_HIGH")

    def test_assess_patient_with_garbage_codes(self):
        res = assess_patient("PT-BAD", "INVALID_CODE; 12345; I21.9", age=30)
        self.assertEqual(res.n_codes, 1)
        self.assertEqual(res.charlson_score, 1)
        self.assertTrue(len(res.warnings) > 0)

    def test_implausible_age_is_not_scored(self):
        res = assess_patient("PT-AGE", "I21.9", age=130)
        self.assertIsNone(res.age)
        self.assertIsNone(res.charlson_age_adjusted)
        self.assertIsNone(res.charlson_10yr_survival_pct)
        self.assertTrue(any("age-adjusted scoring was skipped" in w for w in res.warnings))

    def test_batch_csv_processing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "patients.csv")
            out_path = os.path.join(tmpdir, "output.csv")

            with open(in_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["patient_id", "age", "sex", "icd10_codes"])
                writer.writeheader()
                writer.writerow({"patient_id": "P1", "age": "72", "sex": "F", "icd10_codes": "I21.9; I50.9"})
                writer.writerow({"patient_id": "P2", "age": "40", "sex": "M", "icd10_codes": "E11.9"})

            results = process_csv(in_path, out_path)
            self.assertEqual(len(results), 2)
            self.assertEqual(results[0].charlson_score, 2)
            self.assertEqual(results[0].charlson_age_adjusted, 5)  # 2 + 3

            with open(out_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                self.assertEqual(len(rows), 2)
                self.assertEqual(rows[0]["charlson_age_adjusted"], "5")
                self.assertEqual(rows[1]["charlson_score"], "1")

    def test_batch_requires_named_icd_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "patients.csv")
            out_path = os.path.join(tmpdir, "output.csv")
            with open(in_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["patient_id", "age", "notes"])
                writer.writeheader()
                writer.writerow({"patient_id": "P1", "age": "72", "notes": "I21.9"})
            with self.assertRaisesRegex(ValueError, "ICD code column"):
                process_csv(in_path, out_path)


class TestCLI(unittest.TestCase):
    def test_cli_single_json(self):
        buf = io.StringIO()
        with patch("sys.stdout", new=buf):
            ret = cli_main(["single", "--id", "CLI-01", "--codes", "I21.9; I50.9", "--age", "65", "--json"])
            self.assertEqual(ret, 0)
        data = json.loads(buf.getvalue())
        self.assertEqual(data["patient_id"], "CLI-01")
        self.assertEqual(data["charlson_score"], 2)
        self.assertEqual(data["charlson_age_adjusted"], 4)

    def test_cli_single_formatted(self):
        buf = io.StringIO()
        with patch("sys.stdout", new=buf):
            ret = cli_main(["single", "--id", "CLI-02", "--codes", "E11.9; J44.1", "--age", "55", "--detail"])
            self.assertEqual(ret, 0)
        output = buf.getvalue()
        self.assertIn("CHARLSON & ELIXHAUSER COMORBIDITY INDEX REPORT", output)
        self.assertIn("Raw CCI Score:", output)


if __name__ == "__main__":
    unittest.main()
