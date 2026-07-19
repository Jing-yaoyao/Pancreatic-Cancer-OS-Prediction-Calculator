/* Static OS prediction calculator for M3/M4 Elastic-Net Cox models.
   All computation is client-side; no data leaves the browser. */
(function () {
  "use strict";
  const A = MODEL_ASSETS;
  const TPM = [6, 12, 18];
  let currentModel = "M3";

  // ---- helpers ----
  // interpolate baseline survival S0 at an arbitrary day from the stored grid
  function s0at(model, day) {
    const g = A.models[model].baseline.t_days;
    const s = A.models[model].baseline.S0;
    if (day <= g[0]) return s[0];
    if (day >= g[g.length - 1]) return s[s.length - 1];
    let lo = 0, hi = g.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (g[mid] <= day) lo = mid; else hi = mid; }
    const frac = (day - g[lo]) / (g[hi] - g[lo]);
    return s[lo] + frac * (s[hi] - s[lo]);           // linear interp (curve is fine-grained)
  }
  function rawValue(varName, optIdx) {
    if (A.rawmap[varName]) return A.rawmap[varName][optIdx];  // e.g. T_stage -> 1/2
    return optIdx;                                            // default 0/1
  }
  function linpred(model) {
    const coef = A.models[model].coefficients;
    const used = A.models[model].used_vars;
    let lp = 0;
    used.forEach(function (v) {
      const sel = document.getElementById("f_" + v);
      const x = rawValue(v, parseInt(sel.value, 10));
      lp += coef[v] * x;
    });
    return lp;
  }
  function strataLabel(model, lp) {
    const t = A.models[model].lp_tertiles;
    if (lp < t[0]) return ["Low risk", "low"];
    if (lp < t[1]) return ["Intermediate risk", "med"];
    return ["High risk", "high"];
  }

  // ---- build input panels for the active model ----
  function buildInputs() {
    const used = A.models[currentModel].used_vars;
    // split into two panels: clinical/pathology vs molecular (IHC/RNA)
    const molec = new Set(["IHC_Treg","IHC_epiN","IHC_epiC","RNA_Treg","RNA_epiN","RNA_epiC"]);
    const clin = used.filter(v => !molec.has(v));
    const mol  = used.filter(v => molec.has(v));
    const nice = {
      LVI:"Lymphovascular invasion", PNI:"Perineural invasion", Diff:"Differentiation",
      T_stage:"T stage", N_stage:"N stage", Age:"Age", Sex:"Sex", BMI:"BMI",
      CA199:"CA19-9", CEA:"CEA",
      IHC_Treg:"IHC Treg", IHC_epiN:"IHC epithelial-nuclear", IHC_epiC:"IHC epithelial-cytoplasmic",
      RNA_Treg:"RNAscope Treg", RNA_epiN:"RNAscope epithelial-nuclear", RNA_epiC:"RNAscope epithelial-cytoplasmic"
    };
    function fieldHTML(v) {
      const opts = A.levels[v].map((lab, i) =>
        `<option value="${i}">${lab}</option>`).join("");
      return `<div class="field"><label for="f_${v}">${nice[v] || v}</label>
              <select id="f_${v}">${opts}</select></div>`;
    }
    let html = `<div class="card"><h2>Clinical &amp; histopathology</h2>${clin.map(fieldHTML).join("")}</div>`;
    html += `<div class="card"><h2>Molecular markers${mol.some(v=>v.startsWith("IHC"))?" (IHC + RNAscope)":" (RNAscope)"}</h2>${mol.map(fieldHTML).join("")}</div>`;
    document.getElementById("inputPanels").innerHTML = html;
    used.forEach(v => document.getElementById("f_" + v).addEventListener("change", compute));
  }

  // ---- compute + render ----
  function compute() {
    const model = currentModel;
    const lp = linpred(model);
    const eLp = Math.exp(lp);
    // result cards at 6/12/18 mo
    const tpDays = A.models[model].timepoints_days;
    let cards = "";
    const cls = ["t6","t12","t18"];
    TPM.forEach(function (mo, i) {
      const S = Math.pow(s0at(model, tpDays[i]), eLp);
      const risk = (1 - S) * 100, os = S * 100;
      cards += `<div class="rc ${cls[i]}"><div class="tp">${mo} months</div>
        <div class="risk">${risk.toFixed(1)}%</div>
        <div class="os">death risk<br>OS prob: ${os.toFixed(1)}%</div></div>`;
    });
    document.getElementById("results").innerHTML = cards;
    // strata
    const [slab, scls] = strataLabel(model, lp);
    document.getElementById("strata").innerHTML =
      `Model <b style="background:${model==="M3"?"#7C93C6":"#E070B0"}">${model}</b> &nbsp;
       predicted risk stratum: <b class="${scls}">${slab}</b>
       <span class="muted">(lp = ${lp.toFixed(3)}; tertiles of the training cohort — reference only, not a clinical grade)</span>`;
    drawCurve(model, eLp);
  }

  function drawCurve(model, eLp) {
    const g = A.models[model].baseline.t_days;
    const s0 = A.models[model].baseline.S0;
    const xs = g.map(d => d / 30.4375);                       // days -> months
    const ys = s0.map(v => Math.pow(v, eLp) * 100);
    const col = model === "M3" ? "#7C93C6" : "#E070B0";
    const tpDays = A.models[model].timepoints_days;
    const mx = TPM, my = tpDays.map(d => Math.pow(s0at(model, d), eLp) * 100);
    const line = { x: xs, y: ys, mode: "lines", line: { color: col, width: 3 },
                   name: "Predicted OS", hovertemplate: "%{x:.1f} mo: OS %{y:.1f}%<extra></extra>" };
    const pts = { x: mx, y: my, mode: "markers+text", marker: { color: col, size: 10, line:{color:"#fff",width:1.5} },
                  text: my.map((v,i)=>`${mx[i]}mo`), textposition: "top center",
                  hovertemplate: "%{x} mo: OS %{y:.1f}%<extra></extra>", showlegend:false };
    Plotly.newPlot("curve", [line, pts], {
      margin: { l: 55, r: 18, t: 10, b: 45 },
      xaxis: { title: "Months since surgery", range: [0, Math.max.apply(null, xs) * 1.02], zeroline:false },
      yaxis: { title: "Overall survival probability (%)", range: [0, 100], zeroline:false },
      font: { family: "Liberation Sans, Arimo, Arial, sans-serif", size: 12 },
      paper_bgcolor: "#fff", plot_bgcolor: "#fff", showlegend: false
    }, { displayModeBar: false, responsive: true });
  }

  // ---- model toggle ----
  document.querySelectorAll(".toggle button").forEach(function (btn) {
    btn.addEventListener("click", function () {
      document.querySelectorAll(".toggle button").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      currentModel = btn.getAttribute("data-m");
      const n = A.models[currentModel].used_vars.length;
      document.getElementById("modelHint").textContent = `${n} input variables`;
      buildInputs(); compute();
    });
  });

  // ---- init ----
  document.getElementById("modelHint").textContent =
    `${A.models.M3.used_vars.length} input variables`;
  buildInputs();
  compute();
})();
