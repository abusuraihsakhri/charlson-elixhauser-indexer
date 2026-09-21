# Charlson and Elixhauser Indexer

A Python implementation of Charlson Comorbidity Index scoring, a 31-category Elixhauser comorbidity implementation, and van Walraven weighting from ICD-10-style diagnosis codes. The repository provides a command-line interface, batch CSV processing, and a browser interface that runs the same Python module with Pyodide.

## Features

- Single-patient Charlson and Elixhauser scoring from ICD-10 codes
- Age-adjusted Charlson score and legacy 10-year survival estimate
- van Walraven weighted Elixhauser score
- Batch CSV processing with downloadable results
- Browser interface with light and dark themes
- Client-side Python execution; entered patient data is not sent to an application backend
- Automated tests across Python 3.10, 3.11, and 3.12

## Browser use

Open the GitHub Pages site for this repository once Pages is deployed. The browser interface loads Pyodide 314.0.7 from jsDelivr, then loads charlson_elixhauser.py from the same site. Calculations and CSV processing occur in the browser.

The app does not transmit entered patient data to an application server. Loading Pyodide requires network access to jsDelivr.

## Command line

Score one record:

~~~bash
python cli.py eval --patient-id P001 --age 68 --sex M --codes "I21.9; E11.65; I50.9; N18.3; J44.9"
~~~

Interactive entry:

~~~bash
python cli.py interactive
~~~

Batch processing:

~~~bash
python cli.py batch -i sample.csv -o scored.csv
~~~

The batch input must contain an ICD code column named one of: icd10_codes, icd_codes, codes, icd10, diagnoses, or icd. Optional recognized columns include patient_id/id/patient/mrn/subject_id, age/age_years, and sex/gender.

## Methodological notes

The Charlson implementation uses repository-maintained ICD-10 prefix mappings and standard Charlson hierarchy rules. Diabetes grouping is based on the first character after the E10-E14 category, preventing later code digits from changing the complication class.

The Elixhauser implementation in this repository contains 31 categories with van Walraven weights. It should not be described as the current AHRQ Elixhauser Comorbidity Software Refined model, which is a different implementation.

The legacy Charlson 10-year survival formula is retained for compatibility and historical reference. It is not an individualized prognostic model. The mortality_risk_tier output field is also retained for backwards compatibility; its thresholds are repository-defined and are not a validated mortality model.

ICD input checking is structural and prefix-based. This project is not a complete ICD terminology validator and is not intended for direct clinical decision-making.

## Testing

~~~bash
python -m pip install pytest
python -m pytest -v -p no:zarr
python cli.py batch -i sample.csv -o out_smoke.csv
~~~

CI also compiles the Python sources and checks that the static browser assets reference the Python runtime and scoring module.

## Browser compatibility

The browser app requires WebAssembly and a current Chromium-, Firefox-, or Safari-based browser. Pyodide is loaded as a pinned external runtime.

## License

MIT. See LICENSE.
