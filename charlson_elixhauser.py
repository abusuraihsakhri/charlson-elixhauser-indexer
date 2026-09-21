"""
Charlson and Elixhauser Comorbidity Index Engine
=================================================
Calculates Charlson Comorbidity Index (CCI, Quan 2011 ICD-10 update of Charlson 1987)
and Elixhauser Comorbidity Index (AHRQ 2021 refined, 31 categories) with van Walraven
composite score and 10-year estimated survival calculation.

Author: Dr. Abu Suraih Sakhri
License: MIT
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Set, Tuple, Union


# ---------------------------------------------------------------------------
# ICD Normalization & Code Utilities
# ---------------------------------------------------------------------------

def normalize_icd(code: str) -> str:
    """Uppercase, strip whitespace and punctuation dots. E.g., ' i21.9 ' -> 'I219'."""
    return re.sub(r"[\.\s\-]", "", code.strip().upper()) if code else ""


def code_matches(code: str, prefixes: Tuple[str, ...]) -> bool:
    """Checks if normalized code starts with any prefix in the tuple."""
    norm = normalize_icd(code)
    return any(norm.startswith(normalize_icd(p)) for p in prefixes)


# ---------------------------------------------------------------------------
# Charlson 17 Comorbidity Categories (Quan et al. 2011 / Charlson 1987)
# ---------------------------------------------------------------------------

CHARLSON_DEFS: List[Tuple[str, int, Tuple[str, ...]]] = [
    ("Myocardial infarction", 1, ("I21", "I22", "I252")),
    ("Congestive heart failure", 1, ("I099", "I110", "I130", "I132", "I255", "I42", "I43", "I50", "P294")),
    ("Peripheral vascular disease", 1, ("I70", "I71", "I731", "I738", "I739", "I771", "I790", "I792", "K551", "K558", "K559", "Z958", "Z959")),
    ("Cerebrovascular disease", 1, ("G45", "G46", "H340", "I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69")),
    ("Dementia", 1, ("F00", "F01", "F02", "F03", "F051", "G30", "G311")),
    ("Chronic pulmonary disease", 1, ("I270", "I278", "I279", "J40", "J41", "J42", "J43", "J44", "J45", "J46", "J47", "J60", "J61", "J62", "J63", "J64", "J65", "J66", "J67", "J684", "J701", "J703")),
    ("Rheumatic disease", 1, ("M05", "M06", "M315", "M32", "M33", "M34", "M351", "M353", "M360")),
    ("Peptic ulcer disease", 1, ("K25", "K26", "K27", "K28")),
    ("Mild liver disease", 1, ("B18", "K700", "K701", "K702", "K703", "K709", "K713", "K714", "K715", "K717", "K73", "K74", "K760", "K762", "K763", "K764", "K768", "K769", "Z944")),
    ("Diabetes without complications", 1, ()),  # Managed via classify_diabetes
    ("Diabetes with complications", 2, ()),     # Managed via classify_diabetes
    ("Hemiplegia/Paraplegia", 2, ("G041", "G114", "G801", "G802", "G81", "G82", "G830", "G831", "G832", "G833", "G834", "G839")),
    ("Renal disease", 2, ("I120", "I131", "N032", "N033", "N034", "N035", "N036", "N037", "N052", "N053", "N054", "N055", "N056", "N057", "N18", "N19", "N250", "Z490", "Z491", "Z492", "Z940", "Z992")),
    ("Any malignancy (non-metastatic)", 2, (
        "C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14", "C15",
        "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C30", "C31", "C32", "C33", "C34",
        "C37", "C38", "C39", "C40", "C41", "C43", "C45", "C46", "C47", "C48", "C49", "C50", "C51", "C52", "C53", "C54",
        "C55", "C56", "C57", "C58", "C60", "C61", "C62", "C63", "C64", "C65", "C66", "C67", "C68", "C69", "C70", "C71",
        "C72", "C73", "C74", "C75", "C76", "C81", "C82", "C83", "C84", "C85", "C88", "C90", "C91", "C92", "C93", "C94",
        "C95", "C96", "C97"
    )),
    ("Moderate/severe liver disease", 3, ("I850", "I859", "I864", "I982", "K704", "K711", "K721", "K729", "K765", "K766", "K767")),
    ("Metastatic solid tumor", 6, ("C77", "C78", "C79", "C80")),
    ("AIDS/HIV", 6, ("B20", "B21", "B22", "B24")),
]

DIABETES_PREFIXES = ("E10", "E11", "E12", "E13", "E14")
DIAB_WO_SUFFIXES = ("0", "1", "6", "8", "9")
DIAB_WC_SUFFIXES = ("2", "3", "4", "5", "7")


def classify_diabetes(codes: List[str]) -> Tuple[bool, bool]:
    """
    Returns (uncomplicated_diabetes, complicated_diabetes).
    If complicated is true, uncomplicated is suppressed according to hierarchy rules.
    """
    norm_codes = [normalize_icd(c) for c in codes]
    wo = False
    wc = False
    for c in norm_codes:
        if not code_matches(c, DIABETES_PREFIXES):
            continue
        if len(c) == 3:
            wo = True
            continue
        # Quan-style diabetes grouping is determined by the first character
        # after the three-character E1x category. Looking at later digits
        # misclassifies codes such as E11.65 (hyperglycemia) as complicated.
        category_digit = c[3]
        if category_digit in DIAB_WC_SUFFIXES:
            wc = True
        else:
            wo = True

    if wc:
        return False, True
    return wo, False


def charlson_flags(codes: List[str]) -> Dict[str, bool]:
    """Identifies active Charlson categories applying standard clinical hierarchy."""
    norm_codes = [normalize_icd(c) for c in codes if c]
    flags: Dict[str, bool] = {}
    for label, weight, prefixes in CHARLSON_DEFS:
        if label.startswith("Diabetes"):
            continue
        flags[label] = any(code_matches(c, prefixes) for c in norm_codes) if prefixes else False

    wo, wc = classify_diabetes(norm_codes)
    flags["Diabetes without complications"] = wo
    flags["Diabetes with complications"] = wc

    # Hierarchy Rules:
    # 1. Metastatic solid tumor cancels non-metastatic solid tumor
    if flags.get("Metastatic solid tumor"):
        flags["Any malignancy (non-metastatic)"] = False

    # 2. Moderate/severe liver cancels mild liver
    if flags.get("Moderate/severe liver disease"):
        flags["Mild liver disease"] = False

    return flags


def charlson_score(flags: Dict[str, bool]) -> int:
    """Computes raw unadjusted Charlson Comorbidity Index score."""
    total = 0
    weight_map = {label: w for label, w, _ in CHARLSON_DEFS}
    for label, is_present in flags.items():
        if is_present:
            total += weight_map.get(label, 0)
    return total


def charlson_age_adjusted(raw_score: int, age: Optional[float]) -> Optional[int]:
    """
    Computes age-adjusted Charlson Comorbidity Index (CCI):
      - Age < 50: +0
      - Age 50-59: +1
      - Age 60-69: +2
      - Age 70-79: +3
      - Age >= 80: +4
    """
    if age is None:
        return None
    if age < 50:
        add = 0
    elif age < 60:
        add = 1
    elif age < 70:
        add = 2
    elif age < 80:
        add = 3
    else:
        add = 4
    return raw_score + add


def charlson_10yr_survival(age_adjusted_score: Optional[int]) -> Optional[float]:
    """
    Computes estimated 10-year survival probability:
      S(10) = 0.983 ^ exp(0.9 * CCI_age_adjusted)
    """
    if age_adjusted_score is None:
        return None
    try:
        val = 0.983 ** math.exp(0.9 * age_adjusted_score)
        return max(0.0, min(1.0, val))
    except OverflowError:
        return 0.0


# ---------------------------------------------------------------------------
# Elixhauser 31 Categories & van Walraven (vW 2009) Weights
# ---------------------------------------------------------------------------

ELIXHAUSER_DEFS: List[Tuple[str, int, Tuple[str, ...]]] = [
    ("Congestive heart failure", 7, ("I099", "I110", "I130", "I132", "I255", "I420", "I421", "I422", "I423", "I424", "I425", "I426", "I427", "I428", "I429", "I43", "I50", "P294")),
    ("Cardiac arrhythmias", 5, ("I441", "I442", "I443", "I456", "I459", "I47", "I48", "I4900", "I4901", "I491", "I492", "I493", "I494", "I495", "I497", "I498", "R000", "R001", "R008", "T82817", "T82818", "Z450", "Z95810")),
    ("Valvular disease", -1, ("A520", "I05", "I06", "I07", "I08", "I091", "I0981", "I0989", "I34", "I35", "I36", "I37", "I38", "I39", "Q230", "Q231", "Q232", "Q233", "Z952", "Z953", "Z954")),
    ("Pulmonary circulation disorders", 4, ("I26", "I270", "I272", "I278", "I279", "I280", "I281", "I288", "I289")),
    ("Peripheral vascular disorders", 2, ("I70", "I71", "I731", "I738", "I739", "I771", "I790", "I792", "K551", "K558", "K559", "Z958", "Z959")),
    ("Hypertension uncomplicated", 0, ("I10",)),
    ("Hypertension complicated", 0, ("I11", "I12", "I13", "I15")),
    ("Paralysis", 7, ("G041", "G114", "G801", "G802", "G81", "G82", "G830", "G831", "G832", "G833", "G834")),
    ("Other neurological disorders", 6, ("G10", "G11", "G12", "G13", "G20", "G21", "G22", "G25", "G254", "G255", "G312", "G3181", "G3182", "G3183", "G3184", "G3185", "G3189", "G320", "G35", "G36", "G37", "G40", "G41", "G93", "G934", "R4701", "R56")),
    ("Chronic pulmonary disease", 3, ("I278", "I279", "J40", "J41", "J42", "J43", "J44", "J45", "J46", "J47", "J60", "J61", "J62", "J63", "J64", "J65", "J66", "J67", "J684", "J701", "J703")),
    ("Diabetes uncomplicated", 0, ()),
    ("Diabetes complicated", -3, ()),
    ("Hypothyroidism", 0, ("E00", "E01", "E02", "E03", "E890")),
    ("Renal failure", 5, ("I120", "I131", "N18", "N19", "N250", "Z490", "Z491", "Z492", "Z940", "Z992")),
    ("Liver disease", 11, ("B18", "I850", "I859", "I864", "I982", "K700", "K701", "K702", "K703", "K704", "K709", "K713", "K714", "K715", "K717", "K721", "K729", "K73", "K74", "K760", "K762", "K763", "K764", "K765", "K766", "K767", "K768", "K769", "Z944")),
    ("Peptic ulcer disease", 0, ("K25", "K26", "K27", "K28")),
    ("AIDS/HIV", 0, ("B20", "B21", "B22", "B24")),
    ("Lymphoma", 9, ("C81", "C82", "C83", "C84", "C85", "C88", "C900", "C902", "C96")),
    ("Metastatic cancer", 12, ("C77", "C78", "C79", "C80")),
    ("Solid tumor without metastasis", 4, (
        "C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12", "C13", "C14", "C15",
        "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C30", "C31", "C32", "C33", "C34",
        "C37", "C38", "C39", "C40", "C41", "C43", "C45", "C46", "C47", "C48", "C49", "C50", "C51", "C52", "C53", "C54",
        "C55", "C56", "C57", "C58", "C60", "C61", "C62", "C63", "C64", "C65", "C66", "C67", "C68", "C69", "C70", "C71",
        "C72", "C73", "C74", "C75", "C76"
    )),
    ("Rheumatoid arthritis/collagen", 0, ("L940", "L941", "L943", "M05", "M06", "M080", "M120", "M123", "M30", "M31", "M32", "M33", "M34", "M350", "M351", "M353", "M360")),
    ("Coagulopathy", 3, ("D65", "D66", "D67", "D68", "D691", "D693", "D694", "D695", "D696")),
    ("Obesity", -4, ("E66",)),
    ("Weight loss", 6, ("E40", "E41", "E42", "E43", "E44", "E45", "E46", "R634", "R6381", "R6382", "R6383", "R64")),
    ("Fluid and electrolyte disorders", 5, ("E222", "E86", "E87")),
    ("Blood loss anemia", -2, ("D500",)),
    ("Deficiency anemia", -2, ("D508", "D509", "D510", "D511", "D512", "D513", "D518", "D519", "D520", "D521", "D522", "D528", "D529", "D531", "D538", "D539", "D649")),
    ("Alcohol abuse", 0, ("F10", "E52", "G621", "I426", "K292", "K700", "K703", "K709", "T510", "Z502", "Z7141", "Z721")),
    ("Drug abuse", -7, ("F11", "F12", "F13", "F14", "F15", "F16", "F18", "F19", "Z7151", "Z722")),
    ("Psychoses", 0, ("F20", "F22", "F23", "F24", "F25", "F28", "F29", "F302", "F312", "F315")),
    ("Depression", -3, ("F204", "F313", "F314", "F315", "F32", "F33", "F341", "F412", "F432")),
]

ELIX_VW_WEIGHTS: Dict[str, int] = {label: w for label, w, _ in ELIXHAUSER_DEFS}


def elixhauser_flags(codes: List[str]) -> Dict[str, bool]:
    """Identifies active Elixhauser categories applying standard clinical hierarchy."""
    norm_codes = [normalize_icd(c) for c in codes if c]
    flags: Dict[str, bool] = {}
    for label, _, prefixes in ELIXHAUSER_DEFS:
        if label.startswith("Diabetes"):
            continue
        if label == "Depression":
            flags[label] = any(code_matches(c, ("F204", "F313", "F314", "F315", "F32", "F33", "F341", "F412", "F432")) for c in norm_codes)
        else:
            flags[label] = any(code_matches(c, prefixes) for c in norm_codes) if prefixes else False

    # Diabetes hierarchy
    wo, wc = classify_diabetes(norm_codes)
    flags["Diabetes uncomplicated"] = wo and not wc
    flags["Diabetes complicated"] = wc

    # Hypertension hierarchy
    if flags.get("Hypertension complicated"):
        flags["Hypertension uncomplicated"] = False

    # Cancer hierarchy
    if flags.get("Metastatic cancer"):
        flags["Solid tumor without metastasis"] = False

    return flags


def elixhauser_count(flags: Dict[str, bool]) -> int:
    """Calculates total count of active Elixhauser comorbidities."""
    return sum(1 for v in flags.values() if v)


def elixhauser_van_walraven(flags: Dict[str, bool]) -> int:
    """Calculates van Walraven weighted composite score (-19 to +89)."""
    return sum(ELIX_VW_WEIGHTS.get(k, 0) for k, v in flags.items() if v)


# ---------------------------------------------------------------------------
# Patient Assessment & Structured Output
# ---------------------------------------------------------------------------

@dataclass
class ComorbidityResult:
    patient_id: str
    age: Optional[float] = None
    sex: Optional[str] = None
    n_codes: int = 0
    charlson_flags: Dict[str, bool] = field(default_factory=dict)
    charlson_score: int = 0
    charlson_age_adjusted: Optional[int] = None
    charlson_10yr_survival_pct: Optional[float] = None
    charlson_conditions: List[str] = field(default_factory=list)
    elix_flags: Dict[str, bool] = field(default_factory=dict)
    elix_count: int = 0
    elix_van_walraven: int = 0
    elix_conditions: List[str] = field(default_factory=list)
    mortality_risk_tier: str = "LOW"
    warnings: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def assess_patient(
    patient_id: str,
    raw_codes: Union[str, List[str]],
    age: Optional[float] = None,
    sex: Optional[str] = None,
) -> ComorbidityResult:
    """
    Assesses comorbidity burden for an individual patient from ICD diagnostic codes.
    """
    warnings: List[str] = []

    # Age validation
    if age is not None:
        try:
            age = float(age)
            if age < 0 or age > 125:
                warnings.append(f"Age {age} is clinically implausible (expected 0-125); age-adjusted scoring was skipped.")
                age = None
        except (ValueError, TypeError):
            warnings.append(f"Could not parse age: {age!r}")
            age = None

    # Sex validation
    if sex is not None:
        sex = str(sex).strip().upper()
        if sex not in ("M", "F", "MALE", "FEMALE", "", "?", "OTHER", "O", "U"):
            warnings.append(f"Unrecognized sex descriptor: {sex!r}")

    # Parse and normalize codes
    if isinstance(raw_codes, str):
        parts = re.split(r"[;,\|\s]+", raw_codes.strip()) if raw_codes.strip() else []
    elif isinstance(raw_codes, (list, tuple, set)):
        parts = [str(x) for x in raw_codes]
    else:
        parts = []

    codes = [normalize_icd(p) for p in parts if p.strip()]

    # Validate ICD-10/ICD-10-CM structural syntax only; this is not a code-set lookup.
    valid_codes = []
    for c in codes:
        if re.fullmatch(r"[A-Z][0-9][A-Z0-9]{1,5}", c):
            valid_codes.append(c)
        else:
            warnings.append(f"Unrecognized ICD-10 structure: {c!r}")

    c_flags = charlson_flags(valid_codes)
    c_score = charlson_score(c_flags)
    c_age = charlson_age_adjusted(c_score, age)
    surv = charlson_10yr_survival(c_age)
    c_conds = [k for k, v in c_flags.items() if v]

    e_flags = elixhauser_flags(valid_codes)
    e_count = elixhauser_count(e_flags)
    e_vw = elixhauser_van_walraven(e_flags)
    e_conds = [k for k, v in e_flags.items() if v]

    # Backwards-compatible repository heuristic. These thresholds are not a
    # validated mortality prediction model and should not be used for clinical decisions.
    effective_score = c_age if c_age is not None else c_score
    if effective_score >= 6 or e_vw >= 15:
        tier = "VERY_HIGH"
    elif effective_score >= 4 or e_vw >= 6:
        tier = "HIGH"
    elif effective_score >= 2 or e_vw >= 1:
        tier = "MODERATE"
    else:
        tier = "LOW"

    return ComorbidityResult(
        patient_id=str(patient_id),
        age=age,
        sex=sex,
        n_codes=len(valid_codes),
        charlson_flags=c_flags,
        charlson_score=c_score,
        charlson_age_adjusted=c_age,
        charlson_10yr_survival_pct=(round(surv * 100, 2) if surv is not None else None),
        charlson_conditions=c_conds,
        elix_flags=e_flags,
        elix_count=e_count,
        elix_van_walraven=e_vw,
        elix_conditions=e_conds,
        mortality_risk_tier=tier,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Batch CSV Engine
# ---------------------------------------------------------------------------

CSV_OUTPUT_FIELDS = [
    "patient_id", "age", "sex", "n_codes",
    "charlson_score", "charlson_age_adjusted", "charlson_10yr_survival_pct",
    "charlson_conditions",
    "elixhauser_count", "elixhauser_van_walraven", "elixhauser_conditions",
    "mortality_risk_tier",
    "warnings",
]


def process_csv(input_path: str, output_path: str) -> List[ComorbidityResult]:
    """Processes batch CSV dataset of patient diagnostic records."""
    results: List[ComorbidityResult] = []

    with open(input_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("Input CSV file is empty or invalid.")

        lower_map = {h.strip().lower(): h for h in reader.fieldnames}
        pid_col = None
        for cand in ("patient_id", "id", "patient", "mrn", "subject_id"):
            if cand in lower_map:
                pid_col = lower_map[cand]
                break

        icd_col = None
        for cand in ("icd10_codes", "icd_codes", "codes", "icd10", "diagnoses", "icd"):
            if cand in lower_map:
                icd_col = lower_map[cand]
                break

        if icd_col is None:
            raise ValueError(
                "Input CSV must include an ICD code column named one of: "
                "icd10_codes, icd_codes, codes, icd10, diagnoses, icd."
            )

        age_col = lower_map.get("age") or lower_map.get("age_years")
        sex_col = lower_map.get("sex") or lower_map.get("gender")

        for i, row in enumerate(reader, start=1):
            pid = row.get(pid_col, f"PT-{i}") if pid_col else f"PT-{i}"
            raw_codes = row.get(icd_col, "")
            age_val = row.get(age_col) if age_col else None
            sex_val = row.get(sex_col) if sex_col else None
            results.append(assess_patient(patient_id=pid, raw_codes=raw_codes, age=age_val, sex=sex_val))

    with open(output_path, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=CSV_OUTPUT_FIELDS)
        w.writeheader()
        for r in results:
            w.writerow({
                "patient_id": r.patient_id,
                "age": r.age if r.age is not None else "",
                "sex": r.sex or "",
                "n_codes": r.n_codes,
                "charlson_score": r.charlson_score,
                "charlson_age_adjusted": r.charlson_age_adjusted if r.charlson_age_adjusted is not None else "",
                "charlson_10yr_survival_pct": f"{r.charlson_10yr_survival_pct:.2f}" if r.charlson_10yr_survival_pct is not None else "",
                "charlson_conditions": "; ".join(r.charlson_conditions),
                "elixhauser_count": r.elix_count,
                "elixhauser_van_walraven": r.elix_van_walraven,
                "elixhauser_conditions": "; ".join(r.elix_conditions),
                "mortality_risk_tier": r.mortality_risk_tier,
                "warnings": " | ".join(r.warnings),
            })

    return results
