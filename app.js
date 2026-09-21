"use strict";

const PYODIDE_VERSION = "314.0.7";
const PYODIDE_INDEX_URL = "https://cdn.jsdelivr.net/pyodide/v" + PYODIDE_VERSION + "/full/";

const state = {
  pyodide: null,
  ready: false
};

const els = {
  runtimeBadge: document.getElementById("runtime-badge"),
  themeToggle: document.getElementById("theme-toggle"),
  singleTab: document.getElementById("single-tab"),
  batchTab: document.getElementById("batch-tab"),
  singlePane: document.getElementById("single-pane"),
  batchPane: document.getElementById("batch-pane"),
  form: document.getElementById("single-form"),
  patientId: document.getElementById("patient-id"),
  age: document.getElementById("age"),
  sex: document.getElementById("sex"),
  codes: document.getElementById("codes"),
  analyzeButton: document.getElementById("analyze-button"),
  sampleButton: document.getElementById("sample-button"),
  clearButton: document.getElementById("clear-button"),
  csvFile: document.getElementById("csv-file"),
  batchButton: document.getElementById("batch-button"),
  batchSummary: document.getElementById("batch-summary"),
  resultContext: document.getElementById("result-context"),
  resultEmpty: document.getElementById("result-empty"),
  resultContent: document.getElementById("result-content"),
  metricCharlson: document.getElementById("metric-charlson"),
  metricAgeAdjusted: document.getElementById("metric-age-adjusted"),
  metricElixCount: document.getElementById("metric-elix-count"),
  metricVw: document.getElementById("metric-vw"),
  survivalRow: document.getElementById("survival-row"),
  survivalValue: document.getElementById("survival-value"),
  charlsonConditions: document.getElementById("charlson-conditions"),
  elixhauserConditions: document.getElementById("elixhauser-conditions"),
  warningsBox: document.getElementById("warnings-box"),
  warningsList: document.getElementById("warnings-list")
};

function setRuntimeStatus(text, kind) {
  els.runtimeBadge.textContent = text;
  els.runtimeBadge.classList.remove("ready", "error");
  if (kind) {
    els.runtimeBadge.classList.add(kind);
  }
}

function setTheme(theme) {
  document.documentElement.dataset.theme = theme;
  const dark = theme === "dark";
  els.themeToggle.textContent = dark ? "Light" : "Dark";
  els.themeToggle.setAttribute("aria-label", dark ? "Switch to light mode" : "Switch to dark mode");
  localStorage.setItem("ce-theme", theme);
}

function selectTab(mode) {
  const single = mode === "single";
  els.singlePane.hidden = !single;
  els.batchPane.hidden = single;
  els.singleTab.classList.toggle("active", single);
  els.batchTab.classList.toggle("active", !single);
  els.singleTab.setAttribute("aria-selected", String(single));
  els.batchTab.setAttribute("aria-selected", String(!single));
}

function setBusy(busy, label) {
  els.analyzeButton.disabled = busy || !state.ready;
  els.batchButton.disabled = busy || !state.ready || !els.csvFile.files.length;
  if (label) {
    els.resultContext.textContent = label;
  }
}

function addChips(container, values) {
  container.replaceChildren();
  const items = Array.isArray(values) && values.length ? values : ["None detected"];
  items.forEach(function (value) {
    const span = document.createElement("span");
    span.className = "chip" + (value === "None detected" ? " empty" : "");
    span.textContent = value;
    container.appendChild(span);
  });
}

function renderWarnings(warnings) {
  els.warningsList.replaceChildren();
  if (!Array.isArray(warnings) || warnings.length === 0) {
    els.warningsBox.hidden = true;
    return;
  }
  warnings.forEach(function (warning) {
    const item = document.createElement("li");
    item.textContent = warning;
    els.warningsList.appendChild(item);
  });
  els.warningsBox.hidden = false;
}

function renderResult(result) {
  els.resultEmpty.hidden = true;
  els.resultContent.hidden = false;
  els.metricCharlson.textContent = String(result.charlson_score);
  els.metricAgeAdjusted.textContent = result.charlson_age_adjusted == null ? "—" : String(result.charlson_age_adjusted);
  els.metricElixCount.textContent = String(result.elixhauser_count);
  els.metricVw.textContent = String(result.elixhauser_van_walraven);

  if (result.charlson_10yr_survival_pct == null) {
    els.survivalRow.hidden = true;
  } else {
    els.survivalValue.textContent = Number(result.charlson_10yr_survival_pct).toFixed(2) + "%";
    els.survivalRow.hidden = false;
  }

  addChips(els.charlsonConditions, result.charlson_conditions);
  addChips(els.elixhauserConditions, result.elixhauser_conditions);
  renderWarnings(result.warnings);
  els.resultContext.textContent = result.patient_id + " • " + result.n_codes + " valid code" + (result.n_codes === 1 ? "" : "s");
}

async function initPython() {
  try {
    setRuntimeStatus("Loading Python…", "");
    const pyodide = await loadPyodide({ indexURL: PYODIDE_INDEX_URL });

    const response = await fetch("./charlson_elixhauser.py", { cache: "no-store" });
    if (!response.ok) {
      throw new Error("Could not load charlson_elixhauser.py (" + response.status + ").");
    }

    const source = await response.text();
    pyodide.FS.writeFile("/home/pyodide/charlson_elixhauser.py", source);
    await pyodide.runPythonAsync(
      "import sys\n" +
      "if '/home/pyodide' not in sys.path:\n" +
      "    sys.path.insert(0, '/home/pyodide')\n" +
      "import charlson_elixhauser\n"
    );

    state.pyodide = pyodide;
    state.ready = true;
    setRuntimeStatus("Python ready", "ready");
    els.resultContext.textContent = "Ready";
    els.resultEmpty.querySelector("strong").textContent = "Ready to analyze.";
    els.resultEmpty.querySelector("span").textContent = "Enter ICD-10 codes and select Analyze.";
    setBusy(false);
  } catch (error) {
    console.error(error);
    setRuntimeStatus("Runtime failed", "error");
    els.resultContext.textContent = "Runtime error";
    els.resultEmpty.querySelector("strong").textContent = "Python runtime failed to load.";
    els.resultEmpty.querySelector("span").textContent = "Check your network connection and reload the page.";
  }
}

async function analyzeSingle(event) {
  event.preventDefault();
  if (!state.ready) {
    return;
  }

  const codes = els.codes.value.trim();
  if (!codes) {
    els.codes.setCustomValidity("Enter at least one ICD-10 code.");
    els.codes.reportValidity();
    return;
  }
  els.codes.setCustomValidity("");

  const patientId = els.patientId.value.trim() || "patient";
  const age = els.age.value.trim();
  const sex = els.sex.value;

  setBusy(true, "Analyzing…");
  try {
    state.pyodide.globals.set("patient_id_js", patientId);
    state.pyodide.globals.set("codes_js", codes);
    state.pyodide.globals.set("age_js", age);
    state.pyodide.globals.set("sex_js", sex);

    const payload = await state.pyodide.runPythonAsync(
      "import json\n" +
      "from charlson_elixhauser import assess_patient\n" +
      "_result = assess_patient(patient_id_js, codes_js, age_js or None, sex_js or None)\n" +
      "json.dumps(_result.to_dict())"
    );
    renderResult(JSON.parse(payload));
  } catch (error) {
    console.error(error);
    els.resultContext.textContent = "Analysis error";
    els.resultEmpty.hidden = false;
    els.resultContent.hidden = true;
    els.resultEmpty.querySelector("strong").textContent = "Could not calculate the scores.";
    els.resultEmpty.querySelector("span").textContent = "Review the input and try again.";
  } finally {
    ["patient_id_js", "codes_js", "age_js", "sex_js"].forEach(function (name) {
      state.pyodide.globals.delete(name);
    });
    setBusy(false);
  }
}

async function processBatch() {
  if (!state.ready || !els.csvFile.files.length) {
    return;
  }

  const file = els.csvFile.files[0];
  setBusy(true, "Processing CSV…");
  els.batchSummary.textContent = "Processing " + file.name + "…";

  try {
    const csvText = await file.text();
    state.pyodide.FS.writeFile("/tmp/input.csv", csvText);

    const summaryJson = await state.pyodide.runPythonAsync(
      "import json\n" +
      "from charlson_elixhauser import process_csv\n" +
      "_batch = process_csv('/tmp/input.csv', '/tmp/output.csv')\n" +
      "_n = len(_batch)\n" +
      "json.dumps({\n" +
      "  'count': _n,\n" +
      "  'mean_charlson': round(sum(r.charlson_score for r in _batch) / _n, 2) if _n else 0,\n" +
      "  'mean_elixhauser': round(sum(r.elix_count for r in _batch) / _n, 2) if _n else 0\n" +
      "})"
    );

    const summary = JSON.parse(summaryJson);
    const output = state.pyodide.FS.readFile("/tmp/output.csv", { encoding: "utf8" });
    const blob = new Blob([output], { type: "text/csv;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = file.name.replace(/\.csv$/i, "") + "-scored.csv";
    document.body.appendChild(link);
    link.click();
    link.remove();
    URL.revokeObjectURL(url);

    els.batchSummary.textContent =
      "Processed " + summary.count + " rows. Mean CCI " + summary.mean_charlson +
      "; mean Elixhauser count " + summary.mean_elixhauser + ". Download created.";
    els.resultContext.textContent = "Batch complete";
  } catch (error) {
    console.error(error);
    els.batchSummary.textContent = "Batch processing failed: " + (error.message || String(error));
    els.resultContext.textContent = "Batch error";
  } finally {
    try { state.pyodide.FS.unlink("/tmp/input.csv"); } catch (_) {}
    try { state.pyodide.FS.unlink("/tmp/output.csv"); } catch (_) {}
    setBusy(false);
  }
}

els.themeToggle.addEventListener("click", function () {
  setTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark");
});
els.singleTab.addEventListener("click", function () { selectTab("single"); });
els.batchTab.addEventListener("click", function () { selectTab("batch"); });
els.form.addEventListener("submit", analyzeSingle);
els.sampleButton.addEventListener("click", function () {
  els.patientId.value = "P001";
  els.age.value = "68";
  els.sex.value = "M";
  els.codes.value = "I21.9; E11.65; I50.9; N18.3; J44.9";
  els.codes.focus();
});
els.clearButton.addEventListener("click", function () {
  els.patientId.value = "";
  els.age.value = "";
  els.sex.value = "";
  els.codes.value = "";
  els.resultContent.hidden = true;
  els.resultEmpty.hidden = false;
  els.resultContext.textContent = state.ready ? "Ready" : "Waiting for runtime";
  els.codes.focus();
});
els.csvFile.addEventListener("change", function () {
  const hasFile = els.csvFile.files.length > 0;
  els.batchButton.disabled = !state.ready || !hasFile;
  els.batchSummary.textContent = hasFile ? els.csvFile.files[0].name + " selected." : "No file selected.";
});
els.batchButton.addEventListener("click", processBatch);

setTheme(localStorage.getItem("ce-theme") === "dark" ? "dark" : "light");
selectTab("single");
initPython();