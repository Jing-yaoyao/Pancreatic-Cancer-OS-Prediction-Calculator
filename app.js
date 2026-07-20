/* Static OS prediction calculator for M3/M4 Elastic-Net Cox models.
   All computation is client-side; no data leaves the browser. */
(function () {
  "use strict";
  const A = MODEL_ASSETS;
  const TPM = [6, 12, 18];
  const CINDEX = { M3: "0.71", M4: "0.73" };   // aggregated OOB C-index (real data)
  const MCOL = { M3: "#5B5FA2", M4: "#C98AA6" };
  let currentModel = "M3";

  const nice = {
    LVI:"Lymphovascular invasion", PNI:"Perineural invasion", Diff:"Differentiation",
    T_stage:"T stage", N_stage:"N stage", Age:"Age", Sex:"Sex", BMI:"BMI",
    CA199:"CA19-9", CEA:"CEA",
    IHC_Treg:"IHC Treg", IHC_epiN:"IHC epithelial-nuclear", IHC_epiC:"IHC epithelial-cytoplasmic",
    RNA_Treg:"RNAscope Treg", RNA_epiN:"RNAscope epithelial-nuclear", RNA_epiC:"RNAscope epithelial-cytoplasmic"
  };

  // ---- helpers ----
  function s0at(model, day) {
    const g = A.models[model].baseline.t_days;
    const s = A.models[model].baseline.S0;
    if (day <= g[0]) return s[0];
    if (day >= g[g.length - 1]) return s[s.length - 1];
    let lo = 0, hi = g.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (g[mid] <= day) lo = mid; else hi = mid; }
    const frac = (day - g[lo]) / (g[hi] - g[lo]);
    return s[lo] + frac * (s[hi] - s[lo]);
  }
  function rawValue(varName, optIdx) {
    if (A.rawmap[varName]) return A.rawmap[varName][optIdx];
    return optIdx;
  }
  function contributions(model) {
    // per-feature contribution to the linear predictor = beta_i * x_i
    const coef = A.models[model].coefficients;
    const used = A.models[model].used_vars;
    const out = [];
    used.forEach(function (v) {
      const sel = document.getElementById("f_" + v);
      const x = rawValue(v, parseInt(sel.value, 10));
      out.push({ v: v, contrib: coef[v] * x });
    });
    return out;
  }
  function linpredFrom(contribs) {
    return contribs.reduce((a, c) => a + c.contrib, 0);
  }
  function strataLabel(model, lp) {
    const t = A.models[model].lp_tertiles;
    if (lp < t[0]) return ["Low risk", "low"];
    if (lp < t[1]) return ["Intermediate risk", "med"];
    return ["High risk", "high"];
  }

  // ---- build input fields (grouped, no nested cards) ----
  function buildInputs() {
    const used = A.models[currentModel].used_vars;
    const molec = new Set(["IHC_Treg","IHC_epiN","IHC_epiC","RNA_Treg","RNA_epiN","RNA_epiC"]);
    const clin = used.filter(v => !molec.has(v));
    const mol  = used.filter(v => molec.has(v));
    function fieldHTML(v) {
      const opts = A.levels[v].map((lab, i) => `<option value="${i}">${lab}</option>`).join("");
      return `<div class="field"><label for="f_${v}">${nice[v] || v}</label>
              <select id="f_${v}">${opts}</select></div>`;
    }
    const grpStyle = "font-size:12px;font-weight:700;color:#8B91A1;text-transform:uppercase;letter-spacing:.4px;margin:12px 0 2px;";
    let html = `<div style="${grpStyle};margin-top:0">Clinical &amp; histopathology</div>${clin.map(fieldHTML).join("")}`;
    const molTitle = mol.some(v => v.startsWith("IHC")) ? "Molecular (IHC + RNAscope)" : "Molecular (RNAscope)";
    html += `<div style="${grpStyle}">${molTitle}</div>${mol.map(fieldHTML).join("")}`;
    document.getElementById("inputPanels").innerHTML = html;
    used.forEach(v => document.getElementById("f_" + v).addEventListener("change", compute));
  }

  // ---- compute + render ----
  function compute() {
    const model = currentModel;
    const contribs = contributions(model);
    const lp = linpredFrom(contribs);
    const eLp = Math.exp(lp);

    // result cards (OS probability prominent, death risk as sub-line)
    const tpDays = A.models[model].timepoints_days;
    const cls = ["t6","t12","t18"];
    let cards = "";
    TPM.forEach(function (mo, i) {
      const S = Math.pow(s0at(model, tpDays[i]), eLp);
      const os = S * 100, risk = (1 - S) * 100;
      cards += `<div class="rc ${cls[i]}"><div class="tp">${mo}-month OS</div>
        <div class="risk">${os.toFixed(1)}%</div>
        <div class="os">death risk ${risk.toFixed(1)}%</div></div>`;
    });
    document.getElementById("results").innerHTML = cards;

    // strata
    const [slab, scls] = strataLabel(model, lp);
    document.getElementById("strata").innerHTML =
      `Risk stratum: <b class="tag ${scls}">${slab}</b>
       <span class="muted">&nbsp;(lp = ${lp.toFixed(3)}; training-cohort tertiles — reference only)</span>`;

    drawContrib(contribs);
    drawCurve(model, eLp);
  }

  function drawContrib(contribs) {
    // show ALL model features, sorted by absolute contribution for this patient.
    // red = raises risk, violet = protective, faint grey = currently at baseline (0).
    const sorted = contribs.slice()
                           .sort((a, b) => Math.abs(b.contrib) - Math.abs(a.contrib));
    const maxAbs = Math.max.apply(null, sorted.map(c => Math.abs(c.contrib)).concat([1e-9]));
    let html = "";
    sorted.forEach(function (c) {
      const a = Math.abs(c.contrib);
      const zero = a < 1e-9;
      const w = zero ? 0 : Math.max(3, Math.round(a / maxAbs * 100));
      const col = c.contrib > 0 ? "#C98AA6" : (c.contrib < 0 ? "#5B5FA2" : "#D9DDE3");
      const sign = c.contrib > 0 ? "+" : (c.contrib < 0 ? "\u2212" : "");
      const vtxt = zero ? "0" : (sign + a.toFixed(2));
      html += `<div class="crow"><div class="cn" title="${nice[c.v]||c.v}">${nice[c.v]||c.v}</div>
        <div class="cbar"><span style="width:${w}%;background:${col}"></span></div>
        <div class="cv">${vtxt}</div></div>`;
    });
    document.getElementById("contrib").innerHTML = html;
  }

  function drawCurve(model, eLp) {
    const g = A.models[model].baseline.t_days;
    const s0 = A.models[model].baseline.S0;
    const xs = g.map(d => d / 30.4375);
    const ys = s0.map(v => Math.pow(v, eLp) * 100);
    const col = MCOL[model];
    const tpDays = A.models[model].timepoints_days;
    const mx = TPM, my = tpDays.map(d => Math.pow(s0at(model, d), eLp) * 100);
    const line = { x: xs, y: ys, mode: "lines", line: { color: col, width: 3 },
                   name: "Predicted OS", hovertemplate: "%{x:.1f} mo: OS %{y:.1f}%<extra></extra>" };
    const pts = { x: mx, y: my, mode: "markers+text", marker: { color: col, size: 9, line:{color:"#fff",width:1.5} },
                  text: my.map((v,i)=>`  ${mx[i]}mo`), textposition: "middle right",
                  textfont: { size: 10, color: "#8B91A1" },
                  hovertemplate: "%{x} mo: OS %{y:.1f}%<extra></extra>", showlegend:false };
    Plotly.newPlot("curve", [line, pts], {
      margin: { l: 48, r: 14, t: 8, b: 40 },
      xaxis: { title: "Months since surgery", range: [0, Math.max.apply(null, xs) * 1.02], zeroline:false,
               gridcolor:"#EEF0F3" },
      yaxis: { title: "OS probability (%)", range: [0, 100], zeroline:false, gridcolor:"#EEF0F3" },
      font: { family: "Inter, Source Sans 3, Liberation Sans, Arial, sans-serif", size: 11, color:"#2E3340" },
      paper_bgcolor: "#fff", plot_bgcolor: "#fff", showlegend: false
    }, { displayModeBar: false, responsive: true });
  }

  // ---- model toggle ----
  document.querySelectorAll(".mbtn").forEach(function (btn) {
    btn.addEventListener("click", function () {
      document.querySelectorAll(".mbtn").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
      currentModel = btn.getAttribute("data-m");
      document.getElementById("cindexVal").textContent = CINDEX[currentModel];
      document.getElementById("cindexNote").textContent = "Model " + currentModel + ", internal (OOB)";
      buildInputs(); compute();
    });
  });

  // ---- init ----
  buildInputs();
  compute();
})();
