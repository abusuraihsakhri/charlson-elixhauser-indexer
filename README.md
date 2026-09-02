# Charlson Elixhauser Indexer

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Charlson and Elixhauser Comorbidity Index Engine
=================================================
Calculates Charlson Comorbidity Index (CCI, Quan 2011 ICD-10 update of Charlson 1987)
and Elixhauser Comorbidity Index (AHRQ 2021 refined, 31 categories) with van Walraven
composite score and 10-year estimated survival calculation.

Author: Dr. Abu Suraih Sakhri
License: MIT

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`ComorbidityResult`** — dedicated module for comorbidity result evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  """Calculates total count of active Elixhauser comorbidities."""
  """Calculates van Walraven weighted composite score (-19 to +89)."""
  c_score = charlson_score(c_flags)
  effective_score = c_age if c_age is not None else c_score
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --id <value> --patient-id <value> --codes <value> --age <value>
```

### Parameter Reference
- `--id`: Specifies input measurement or parameter value.
- `--patient-id`: Specifies input measurement or parameter value.
- `--codes`: Specifies input measurement or parameter value.
- `--age`: Specifies input measurement or parameter value.
- `--sex`: Specifies input measurement or parameter value.
- `--detail`: Specifies input measurement or parameter value.
- `--json`: Specifies input measurement or parameter value.
- `--input`: Specifies input measurement or parameter value.
- `--output`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `case_id` | Parameter / observation metric | Required |
| `description` | Parameter / observation metric | Required |
| `inputs` | Parameter / observation metric | Required |
| `expected` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t charlson-elixhauser-indexer .
docker run -p 8000:8000 charlson-elixhauser-indexer
```
