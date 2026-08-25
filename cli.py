#!/usr/bin/env python3
"""
CLI for Charlson & Elixhauser Comorbidity Indexer
Supports single patient evaluation, interactive clinical consultation, batch CSV indexing, and JSON output.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Optional

from charlson_elixhauser import (
    assess_patient,
    process_csv,
    ComorbidityResult,
    ELIX_VW_WEIGHTS,
)


def format_patient_report(r: ComorbidityResult, detail: bool = False) -> str:
    """Formats clinical summary for a patient."""
    lines = []
    lines.append("=" * 74)
    lines.append("  CHARLSON & ELIXHAUSER COMORBIDITY INDEX REPORT")
    lines.append("  Standards: Quan 2011 (ICD-10 CCI) / AHRQ 2021 & van Walraven (2009)")
    lines.append("=" * 74)

    lines.append(f"\n  Patient Identifier:  {r.patient_id}")
    lines.append(f"  Demographics:        Age: {r.age if r.age is not None else 'N/A'} | Sex: {r.sex or 'N/A'}")
    lines.append(f"  Valid ICD Codes:     {r.n_codes}")

    lines.append("\n  [Charlson Comorbidity Index (CCI)]")
    lines.append(f"  * Raw CCI Score:            {r.charlson_score}")
    if r.charlson_age_adjusted is not None:
        lines.append(f"  * Age-Adjusted CCI Score:   {r.charlson_age_adjusted}")
    if r.charlson_10yr_survival_pct is not None:
        lines.append(f"  * Estimated 10-Yr Survival: {r.charlson_10yr_survival_pct:.2f}%")
    cond_str = ", ".join(r.charlson_conditions) if r.charlson_conditions else "None detected"
    lines.append(f"  * CCI Conditions Met:       {cond_str}")

    lines.append("\n  [Elixhauser Comorbidity Index (ECI)]")
    lines.append(f"  * Total Comorbidity Count:  {r.elix_count} / 31")
    lines.append(f"  * van Walraven Index (vW):  {r.elix_van_walraven} (Range: -19 to +89)")
    elix_str = ", ".join(r.elix_conditions) if r.elix_conditions else "None detected"
    lines.append(f"  * ECI Conditions Met:       {elix_str}")

    lines.append(f"\n  [Mortality Risk Tier]:      [{r.mortality_risk_tier}]")

    if r.warnings:
        lines.append("\n  [!] Warnings:")
        for w in r.warnings:
            lines.append(f"      - {w}")

    if detail:
        lines.append("\n  [Detailed Condition Breakdown]:")
        lines.append("  -- Charlson 17 Criteria --")
        for k, v in r.charlson_flags.items():
            lines.append(f"     [{'X' if v else ' '}] {k}")
        lines.append("  -- Elixhauser 31 Criteria & Weights --")
        for k, v in r.elix_flags.items():
            w = ELIX_VW_WEIGHTS.get(k, 0)
            lines.append(f"     [{'X' if v else ' '}] {k:<35} (weight: {w:+d})")

    lines.append("=" * 74)
    return "\n".join(lines)


def interactive_mode():
    """Interactive CLI consultation."""
    print("=" * 74)
    print("  CHARLSON & ELIXHAUSER COMORBIDITY INDEXER - INTERACTIVE MODE")
    print("=" * 74)

    pid = input("Enter Patient ID [PT-001]: ").strip() or "PT-001"
    age_str = input("Enter Patient Age in years (e.g. 68) [optional]: ").strip()
    age = float(age_str) if age_str else None
    sex = input("Enter Patient Sex (M/F) [optional]: ").strip() or None

    print("\nEnter ICD-10 diagnostic codes separated by commas, semicolons, or spaces.")
    print("Example: I21.9, I50.9, E11.65, N18.3, C34.9, J44.9")
    codes_input = input("ICD-10 Codes: ").strip()

    res = assess_patient(patient_id=pid, raw_codes=codes_input, age=age, sex=sex)
    print("\n" + format_patient_report(res, detail=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="charlson-elixhauser-indexer",
        description="Charlson Comorbidity Index & Elixhauser / van Walraven Index Calculator",
    )
    subparsers = parser.add_subparsers(dest="command")

    for cmd_name in ("single", "eval"):
        s = subparsers.add_parser(cmd_name, help="Score an individual patient")
        s.add_argument("--id", "--patient-id", dest="patient_id", default="patient", help="Patient Identifier")
        s.add_argument("--codes", required=True, help='ICD-10 codes separated by ; , or space (e.g., "I21.9; E11.65; I50.9")')
        s.add_argument("--age", type=float, default=None, help="Patient age in years")
        s.add_argument("--sex", default=None, help="Patient sex (M/F)")
        s.add_argument("--detail", action="store_true", help="Display full 17 Charlson & 31 Elixhauser breakdown")
        s.add_argument("--json", action="store_true", help="Output results in JSON format")

    subparsers.add_parser("interactive", help="Run interactive clinical ICD questionnaire")

    b = subparsers.add_parser("batch", help="Process batch CSV dataset")
    b.add_argument("-i", "--input", required=True, help="Input CSV file path")
    b.add_argument("-o", "--output", default="scored_comorbidities.csv", help="Output CSV file path")
    b.add_argument("--detail", action="store_true", help="Print summary of each patient")

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "interactive" or (args.command is None and len(sys.argv) == 1):
        interactive_mode()
        return 0

    if args.command in ("single", "eval"):
        res = assess_patient(patient_id=args.patient_id, raw_codes=args.codes, age=args.age, sex=args.sex)
        if args.json:
            print(json.dumps(res.to_dict(), indent=2))
        else:
            print(format_patient_report(res, detail=args.detail))
        return 0

    if args.command == "batch":
        results = process_csv(args.input, args.output)
        n = len(results)
        avg_c = sum(r.charlson_score for r in results) / n if n else 0
        avg_e = sum(r.elix_count for r in results) / n if n else 0
        print(f"Processed {n} patients -> {args.output}")
        print(f"Mean Charlson Score: {avg_c:.2f} | Mean Elixhauser Count: {avg_e:.2f}")
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
