import { useEffect, useMemo, useRef, useState } from "react";
import { Network, Play, Loader2, AlertTriangle, FileText } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow, SectionTitle, StatCard, LimitBanner, ExplainBox } from "../ui/Primitives.jsx";
import { useFastqData } from "../../state/FastqDataContext.jsx";
import { useAssoc } from "../../hooks/useAssoc.js";
import { parseMatrix, parsePhenotypes, matchSamples } from "../../lib/assoc/parseInputs.js";
import { treeClades } from "../../lib/assoc/treeGroups.js";
import { generateExampleGwas } from "../../lib/sampleData/generateExampleGwas.js";
import { buildResultsTsv, buildSignificantTsv } from "../../lib/assoc/exportAssoc.js";
import { downloadBlob } from "../../lib/fastq/exportFastq.js";
import VolcanoPlot from "./VolcanoPlot.jsx";
import ResultsTable, { RankedBar } from "./ResultsTable.jsx";

const CASE_HINTS = /resist|case|positive|^pos$|yes|surviv|pathog|sympt|high/i;

const STAGE_LABEL = {
  tests: "Per-feature association tests",
  correction: "Multiple-testing correction",
};

function saveText(text, fileName) {
  downloadBlob(new Blob([text], { type: "text/plain;charset=utf-8" }), fileName);
}

export default function AssocSection({ explainMode }) {
  const { phyl } = useFastqData();
  const assoc = useAssoc();

  const matrixRef = useRef(null), metaRef = useRef(null);
  const [matrix, setMatrix] = useState(null);   // {genes, samples, variableGenes, format, name}
  const [meta, setMeta] = useState(null);       // {samples, columns, name}
  const [traitColIdx, setTraitColIdx] = useState(-1);
  const [caseLabel, setCaseLabel] = useState("");
  const [strataChoice, setStrataChoice] = useState("none"); // none | meta:<i> | tree
  const [treeCutoff, setTreeCutoff] = useState(0.05);
  const [fdr, setFdr] = useState(0.05);
  const [statPref, setStatPref] = useState(null); // null = auto | "p" | "cmhP"
  const [selectedGene, setSelectedGene] = useState(null);
  const [parseError, setParseError] = useState("");

  const treeAvailable = !!(phyl.status === "done" && phyl.result?.tree);

  /* ------------------------- input handling ------------------------- */

  const handleMatrixFile = (file) => {
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const parsed = parseMatrix(String(reader.result));
        setMatrix({ ...parsed, name: file.name });
        setSelectedGene(null);
        setParseError("");
      } catch (err) { setParseError(`Matrix: ${err.message}`); }
    };
    reader.readAsText(file);
  };

  const handleMetaFile = (file) => {
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const parsed = parsePhenotypes(String(reader.result));
        if (!parsed.columns.some((c) => c.kind !== "text")) throw new Error("No testable column found - need a binary or numeric trait.");
        setMeta({ ...parsed, name: file.name });
        setStrataChoice("none");
        autoPickTrait(parsed);
      } catch (err) { setParseError(`Metadata: ${err.message}`); }
    };
    reader.readAsText(file);
  };

  function autoPickTrait(parsed) {
    let idx = parsed.columns.findIndex((c) => c.kind === "binary");
    if (idx < 0) idx = parsed.columns.findIndex((c) => c.kind === "numeric");
    setTraitColIdx(idx);
    const col = parsed.columns[idx];
    if (col?.kind === "binary") {
      const uniq = [...new Set(col.values.filter(Boolean))];
      const hinted = uniq.find((v) => CASE_HINTS.test(v));
      setCaseLabel(hinted ?? uniq.sort()[uniq.length - 1]);
    }
  }

  function loadExample() {
    const ex = generateExampleGwas(7);
    const m = { genes: ex.genes.map((g) => ({ name: g.name, values: g.values, present: g.present })), samples: ex.samples, variableGenes: ex.genes.length, format: "simulated example panel", name: "Example: clonal outbreak panel", products: Object.fromEntries(ex.genes.map((g) => [g.name, g.product])) };
    const md = { samples: ex.samples, columns: ex.meta, name: "Example metadata" };
    setMatrix(m); setMeta(md); setSelectedGene(null); setStatPref(null); setStrataChoice("meta:1"); setParseError("");
    autoPickTrait(md);
  }

  /* --------------------------- derived data --------------------------- */

  const match = useMemo(() => (
    matrix && meta ? matchSamples(matrix.samples, meta) : null
  ), [matrix, meta]);

  const traitCol = traitColIdx >= 0 ? meta?.columns[traitColIdx] : null;
  const traitType = traitCol?.kind === "numeric" ? "continuous" : traitCol?.kind === "binary" ? "binary" : null;
  const controlLabel = useMemo(() => {
    if (traitCol?.kind !== "binary") return null;
    return [...new Set(traitCol.values.filter(Boolean))].find((v) => v !== caseLabel) ?? null;
  }, [traitCol, caseLabel]);

  const groupCounts = useMemo(() => {
    if (!match || !traitCol || traitType !== "binary") return null;
    const idx = new Map(meta.samples.map((s, i) => [s, i]));
    let nCase = 0, nCtrl = 0;
    for (const s of match.samples) {
      const v = traitCol.values[idx.get(s)];
      if (v === caseLabel) nCase++;
      else if (v != null && v !== "") nCtrl++;
    }
    return { nCase, nCtrl };
  }, [match, traitCol, traitType, caseLabel, meta]);

  const canRun = !!matrix && !!traitCol && !!match && match.keptIdx.length >= 4 &&
    (traitType === "continuous" || (groupCounts && groupCounts.nCase >= 2 && groupCounts.nCtrl >= 2));

  const problems = [];
  if (matrix && meta && match.keptIdx.length < matrix.samples.length) {
    problems.push(`${matrix.samples.length - match.keptIdx.length} matrix sample(s) have no metadata row and will be dropped.`);
  }
  if (matrix && meta && match.unmatchedMeta.length) {
    problems.push(`${match.unmatchedMeta.length} metadata row(s) match no matrix sample.`);
  }
  if (groupCounts && (groupCounts.nCase < 2 || groupCounts.nCtrl < 2)) {
    problems.push(`Need ≥2 samples labelled "${caseLabel}" and ≥2 "${controlLabel}" among matched samples.`);
  }

  /* ------------------------------ run ------------------------------ */

  function run() {
    const kept = match.keptIdx;
    const idxInMeta = new Map(meta.samples.map((s, i) => [s, i]));
    const samples = kept.map((i) => matrix.samples[i]);

    // Trait vector aligned to matched samples.
    let traitValues;
    if (traitType === "binary") {
      traitValues = samples.map((s) => (traitCol.values[idxInMeta.get(s)] === caseLabel ? 1 : 0));
    } else {
      traitValues = samples.map((s) => {
        const v = traitCol.values[idxInMeta.get(s)];
        return Number.isFinite(v) ? v : NaN;
      });
    }
    if (traitValues.some((v) => Number.isNaN(v))) {
      setParseError("Trait column has missing numeric values for matched samples - fill or drop those rows first.");
      return;
    }

    // Flatten presence matrix to matched samples only.
    const flat = new Uint8Array(matrix.genes.length * samples.length);
    for (let g = 0; g < matrix.genes.length; g++) {
      const src = matrix.genes[g].values;
      for (let k = 0; k < kept.length; k++) flat[g * samples.length + k] = src[kept[k]] ? 1 : 0;
    }

    // Strata vector for population-structure correction.
    let strata = null;
    if (strataChoice.startsWith("meta:")) {
      const ci = Number(strataChoice.slice(5));
      const col = meta.columns[ci];
      strata = samples.map((s) => String(col.values[idxInMeta.get(s)] ?? ""));
    } else if (strataChoice === "tree") {
      const groups = treeClades(phyl.result.tree, treeCutoff);
      strata = samples.map((s) => groups.get(s) ?? `unplaced_${s}`);
    }

    setSelectedGene(null);
    setStatPref(null);
    assoc.run({
      samples,
      geneNames: matrix.genes.map((g) => g.name),
      matrix: flat,
      traitValues, traitType,
      strata,
      products: matrix.products ? matrix.genes.map((g) => matrix.products[g.name] ?? null) : null,
      echo: { traitName: traitCol.name, caseLabel, controlLabel },
    });
  }

  useEffect(() => () => assoc.cancel(), []); // eslint-disable-line react-hooks/exhaustive-deps

  const res = assoc.result;
  const rows = useMemo(() => res?.rows ?? [], [res]);
  const hasCMH = !!(res?.meta.hasStrata && res.meta.nStrata >= 2);
  // Displayed statistic: user preference, else auto (stratified when available).
  // `statPref` is cleared whenever a new scan starts, so fresh results re-auto.
  const statView = statPref ?? (hasCMH ? "cmhP" : "p");
  const activeFdr = fdr;

  const selectedRow = rows.find((r) => r.gene === selectedGene);
  const sigCount = useMemo(() => rows.filter((r) => {
    const q = statView === "cmhP" ? r.cmhP : r.q;
    return Number.isFinite(q) && q <= activeFdr;
  }).length, [rows, statView, activeFdr]);

  return (
    <div>
      <SectionTitle icon={Network} color={C.pheno}
        title="Pan-GWAS · Genotype ↔ Phenotype"
        subtitle="Fisher-exact association testing over gene presence/absence, with FDR correction and clade-stratified structure control - runs entirely in your browser" />

      <LimitBanner>
        Bacterial genomes are clonal, so shared ancestry makes unrelated genes travel together: a lineage marker can look
        perfectly associated while carrying no mechanism. Always read the raw Fisher test next to the stratified
        Mantel-Haenszel result (supply clade labels or cut this session's tree), and remember that even a robust
        association is a lead, not causation - confirm by curating context, checking linkage, and experiment.
      </LimitBanner>

      {/* ------------------------------ INPUT ------------------------------ */}
      <Panel style={{ padding: 18, marginTop: 16, marginBottom: 16 }}>
        <Eyebrow color={C.pheno}>Data sources</Eyebrow>
        <div style={{ display: "flex", gap: 10, marginTop: 8, flexWrap: "wrap" }}>
          <SourceCard title="Simulated example panel" color={C.assembly}
            sub="56 strains · 5 clades · causal vs confounded genes with known ground truth"
            onClick={loadExample} />
          <SourceCard title="Presence/absence matrix" color={C.raw}
            sub={matrix ? `${matrix.name} — ${matrix.genes.length.toLocaleString()} genes × ${matrix.samples.length} samples` : "ROARY/Panaroo .Rtab or gene_presence_absence.csv"}
            onClick={() => matrixRef.current.click()} />
          <SourceCard title="Phenotype metadata" color={C.phylo}
            sub={meta ? `${meta.name} — ${meta.samples.length} rows, ${meta.columns.length} column(s)` : "CSV/TSV: sample id + trait (+ optional clade labels)"}
            onClick={() => metaRef.current.click()} />
        </div>
        <input ref={matrixRef} type="file" accept=".rtab,.csv,.tsv,.txt" hidden onChange={(e) => e.target.files[0] && handleMatrixFile(e.target.files[0])} />
        <input ref={metaRef} type="file" accept=".csv,.tsv,.txt" hidden onChange={(e) => e.target.files[0] && handleMetaFile(e.target.files[0])} />

        {parseError && (
          <div style={{ marginTop: 10, fontSize: 12, color: "#f2b3ad", display: "flex", gap: 7 }}>
            <AlertTriangle size={13} style={{ flexShrink: 0, marginTop: 1 }} /> {parseError}
          </div>
        )}

        {matrix && !meta && (
          <div style={{ fontSize: 12, color: C.textDim, marginTop: 10 }}>
            Matrix loaded ({matrix.format}). Now load a phenotype/metadata table to enable testing.
          </div>
        )}

        {/* --------------------- trait & strata configuration --------------------- */}
        {meta && traitCol && (
          <>
            <div style={{ marginTop: 16 }}><Eyebrow color={C.pheno}>Test configuration</Eyebrow></div>
            <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(240px, 1fr))", gap: 18, marginTop: 8 }}>
              <label style={{ fontSize: 12.5, color: C.text }}>
                Trait column
                <select value={traitColIdx} onChange={(e) => { const i = +e.target.value; setTraitColIdx(i); autoPickLabel(meta.columns[i]); }}
                  style={selectStyle}>
                  {meta.columns.map((c, i) => (
                    <option key={c.name} value={i}>{c.name} ({c.kind})</option>
                  ))}
                </select>
              </label>

              {traitType === "binary" && (
                <div style={{ fontSize: 12.5, color: C.text }}>
                  Case group (others become controls)
                  <div style={{ display: "flex", gap: 6, marginTop: 6 }}>
                    {[...new Set(traitCol.values.filter(Boolean))].map((v) => (
                      <Chip key={v} active={v === caseLabel} onClick={() => setCaseLabel(v)}>{v}</Chip>
                    ))}
                  </div>
                </div>
              )}

              <label style={{ fontSize: 12.5, color: C.text }}>
                Structure correction (CMH strata)
                <select value={strataChoice} onChange={(e) => setStrataChoice(e.target.value)} style={selectStyle}>
                  <option value="none">None — pooled test only</option>
                  {meta.columns.map((c, i) => c.kind === "text" && (
                    <option key={c.name} value={`meta:${i}`}>Column: {c.name}</option>
                  ))}
                  <option value="tree" disabled={!treeAvailable}>
                    Session phylogeny{treeAvailable ? "" : " (run Phylogeny first)"}
                  </option>
                </select>
              </label>

              {strataChoice === "tree" && (
                <label style={{ fontSize: 12.5, color: C.text }}>
                  Clade cut-off (subs/site)
                  <input type="number" min={0.0005} max={0.5} step={0.005} value={treeCutoff}
                    onChange={(e) => setTreeCutoff(Math.max(0.0005, +e.target.value || 0.05))}
                    style={selectStyle} />
                </label>
              )}

              <label style={{ fontSize: 12.5, color: C.text }}>
                Significance threshold (FDR)
                <select value={fdr} onChange={(e) => setFdr(+e.target.value)} style={selectStyle}>
                  {[0.01, 0.05, 0.1].map((v) => <option key={v} value={v}>{v}</option>)}
                </select>
              </label>
            </div>

            <div style={{ fontSize: 11.5, color: C.textFaint, marginTop: 10 }}>
              Matched: {match.keptIdx.length} of {matrix?.samples.length ?? 0} matrix samples
              {traitType === "binary" && groupCounts ? <> · cases ({caseLabel}): {groupCounts.nCase} · controls ({controlLabel}): {groupCounts.nCtrl}</> : null}
              {traitType === "continuous" ? <> · quantitative trait — Welch's t-test replaces Fisher exact</> : null}
            </div>
            {problems.map((p, i) => (
              <div key={i} style={{ fontSize: 11.5, color: C.qc, marginTop: 5 }}>⚠ {p}</div>
            ))}
          </>
        )}

        <div style={{ display: "flex", alignItems: "center", gap: 12, marginTop: 16, flexWrap: "wrap" }}>
          <button onClick={run} disabled={!canRun || assoc.status === "running"}
            style={{
              all: "unset", cursor: canRun && assoc.status !== "running" ? "pointer" : "default",
              display: "inline-flex", alignItems: "center", gap: 7, padding: "9px 16px", borderRadius: 8,
              border: `1px solid ${C.pheno}66`, color: C.pheno, fontSize: 13,
              opacity: canRun && assoc.status !== "running" ? 1 : 0.55,
            }}>
            {assoc.status === "running" ? <Loader2 size={14} className="spin" /> : <Play size={14} />}
            {assoc.status === "running" ? "Testing…" : matrix?.genes.length > 3000 ? "Run genome-wide scan" : "Run association scan"}
          </button>
          {assoc.status === "running" && (
            <button onClick={assoc.cancel} style={{ all: "unset", cursor: "pointer", fontSize: 12, color: C.bad, border: `1px solid ${C.bad}66`, borderRadius: 6, padding: "7px 12px" }}>
              Cancel
            </button>
          )}
          {!canRun && matrix && traitCol && (
            <span style={{ fontSize: 11.5, color: C.textFaint }}>resolve the warnings above to enable testing</span>
          )}
        </div>

        <style>{`.spin { animation: spin 1s linear infinite; } @keyframes spin { from{transform:rotate(0)} to{transform:rotate(360deg)} }`}</style>

        {assoc.status === "running" && (
          <div style={{ marginTop: 14 }}>
            <div style={{ fontSize: 12, color: C.textDim, marginBottom: 6 }}>
              {STAGE_LABEL[assoc.stage] || "Working…"}{assoc.detail ? ` — ${assoc.detail}` : ""}
            </div>
            <div style={{ height: 6, background: "#05070a", borderRadius: 3, overflow: "hidden" }}>
              <div style={{ width: `${assoc.pct}%`, height: "100%", background: C.pheno, transition: "width .2s" }} />
            </div>
          </div>
        )}
        {assoc.error && (
          <div style={{ marginTop: 12, display: "flex", gap: 8, fontSize: 12.5, color: "#f2b3ad" }}>
            <AlertTriangle size={14} /> {assoc.error}
          </div>
        )}

        <ExplainBox explainMode={explainMode} color={C.pheno}>
          Each gene becomes a 2×2 table (present/absent × case/control). Fisher's exact test asks how surprising that
          split is if gene carriage and phenotype were independent; the volcano plots effect size against that surprise.
          The stratified test then repeats the comparison inside each clade separately - a confounder that only looked
          good because resistant strains are relatives loses its shine here, while a real determinant carried across
          lineages keeps it.
        </ExplainBox>
      </Panel>

      {/* ------------------------------ RESULTS ------------------------------ */}
      {assoc.status === "done" && res && (
        <>
          <StatsRow res={res} sigCount={sigCount} statView={statView} fdr={activeFdr} />

          <Panel style={{ padding: 16, marginBottom: 16 }}>
            <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", gap: 8, flexWrap: "wrap", marginBottom: 8 }}>
              <Eyebrow color={C.pheno}>Effect vs evidence</Eyebrow>
              <div style={{ display: "flex", gap: 6, alignItems: "center" }}>
                <Chip active={statView === "p"} onClick={() => setStatPref("p")}>
                  {res.meta.traitType === "binary" ? "Fisher (pooled)" : "Welch t-test"}
                </Chip>
                <Chip active={statView === "cmhP"} onClick={() => hasCMH && setStatPref("cmhP")} dim={!hasCMH}>
                  Stratified CMH
                </Chip>
              </div>
            </div>
            {statView === "cmhP" && (
              <div style={{ fontSize: 11.5, color: C.textFaint, marginBottom: 6 }}>
                Y axis uses the Cochran-Mantel-Haenszel p across {res.meta.nStrata} strata (x keeps the pooled odds
                ratio so the two views stay comparable). Non-informative strata - clades where every strain matches on
                the gene or the phenotype - contribute nothing; that is the point.
              </div>
            )}
            {res.meta.traitType === "binary"
              ? <VolcanoPlot rows={rows} stat={statView} fdr={activeFdr} selected={selectedGene} onSelect={setSelectedGene} />
              : <RankedBar rows={rows} stat={statView} fdr={activeFdr} />}
            {selectedRow && <SelectedStrip row={selectedRow} traitType={res.meta.traitType} />}
          </Panel>

          <Panel style={{ padding: 16, marginBottom: 16 }}>
            <Eyebrow color={C.pheno}>All tested features</Eyebrow>
            <div style={{ marginTop: 8 }}>
              <ResultsTable rows={rows} traitType={res.meta.traitType} stat={statView} fdr={activeFdr}
                selected={selectedGene} onSelect={setSelectedGene} />
            </div>
          </Panel>

          <div style={{ display: "flex", gap: 8, flexWrap: "wrap", margin: "4px 0 14px" }}>
            <DownloadBtn label="Full results TSV" onClick={() => saveText(buildResultsTsv(rows, res.meta.traitType), "pan_gwas_results.tsv")} />
            <DownloadBtn label={`Significant hits (q ≤ ${activeFdr})`} onClick={() => saveText(buildSignificantTsv(rows, res.meta.traitType, activeFdr), "pan_gwas_significant.tsv")} />
          </div>

          <ExplainBox explainMode={explainMode} color={C.pheno}>
            Compare the two views before believing any hit: genes that survive both the pooled and the stratified test
            are your strongest candidates; genes significant only when pooled deserve suspicion of being inherited
            lineage markers. The exported TSV carries every intermediate number so reviewers can audit the calls.
          </ExplainBox>
        </>
      )}

      {assoc.status !== "done" && !canRun && (
        <Panel style={{ padding: 16, fontSize: 12.5, color: C.textFaint }}>
          No scan yet. Load a presence/absence matrix plus phenotype metadata (or start from the simulated panel),
          choose a trait, then press Run. ROARY/Panaroo outputs drop straight in.
        </Panel>
      )}
    </div>
  );

  function autoPickLabel(col) {
    if (col.kind === "binary") {
      const uniq = [...new Set(col.values.filter(Boolean))];
      const hinted = uniq.find((v) => CASE_HINTS.test(v));
      setCaseLabel(hinted ?? uniq.sort()[uniq.length - 1]);
    }
  }
}

/* ------------------------------- pieces ------------------------------- */

const selectStyle = {
  display: "block", width: "100%", marginTop: 6, background: "#05070a", border: `1px solid ${C.border}`,
  borderRadius: 6, padding: "6px 9px", color: C.text, fontFamily: FONT_DISPLAY, fontSize: 12.5,
};

function SourceCard({ title, sub, onClick }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", flex: "1 1 220px", padding: "12px 14px", borderRadius: 8,
      background: "#05070a", border: `1px solid ${C.border}`,
    }}>
      <div style={{ fontSize: 12.5, color: C.text }}>{title}</div>
      <div style={{ fontSize: 11, color: C.textFaint, marginTop: 3 }}>{sub}</div>
    </button>
  );
}

function Chip({ children, active, onClick, dim }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: dim ? "default" : "pointer", fontSize: 11, padding: "4px 10px", borderRadius: 6,
      background: active ? `${C.pheno}22` : "transparent",
      border: `1px solid ${active ? C.pheno : C.border}`,
      color: dim ? C.textFaint : active ? C.pheno : C.textDim,
      opacity: dim ? 0.55 : 1,
    }}>{children}</button>
  );
}

function StatsRow({ res, sigCount, statView, fdr }) {
  const m = res.meta;
  const top = res.rows.find((r) => {
    const q = statView === "cmhP" ? r.cmhP : r.q;
    return Number.isFinite(q);
  });
  return (
    <div style={{ display: "flex", gap: 10, flexWrap: "wrap", marginBottom: 14 }}>
      <StatCard label="Features tested" value={m.nGenes.toLocaleString()} color={C.pheno} />
      <StatCard label="Samples" value={m.nSamples} color={C.raw} />
      {m.traitType === "binary"
        ? <StatCard label={`Cases / controls`} value={`${m.nCase} / ${m.nCtrl}`} color={C.qc} />
        : <StatCard label="Design" value="quant." unit="trait" color={C.qc} />}
      <StatCard label={statView === "cmhP" ? `Significant (CMH p ≤ ${fdr})` : `Significant (q ≤ ${fdr})`}
        value={sigCount} color={sigCount ? C.bad : C.good} />
      <PairCard text={top?.gene ?? "-"} sub={top?.product ?? ""} />
      {m.hasStrata && <StatCard label="CMH strata" value={m.nStrata} color={C.phylo} />}
      <StatCard label="Runtime" value={(m.ms / 1000).toFixed(2)} unit="s" color={C.textDim} />
    </div>
  );
}

function PairCard({ text, sub }) {
  return (
    <div style={{ background: C.bgPanel2, border: `1px solid ${C.border}`, borderRadius: 8, padding: "10px 14px", minWidth: 150 }}>
      <div style={{ fontSize: 10.5, color: C.textFaint, textTransform: "uppercase", letterSpacing: "0.08em" }}>Top-ranked feature</div>
      <div style={{ fontFamily: FONT_DISPLAY, fontSize: 13, color: C.pheno, marginTop: 3 }}>{text}</div>
      {sub && <div style={{ fontSize: 10.5, color: C.textDim, marginTop: 1 }}>{sub}</div>}
    </div>
  );
}

function SelectedStrip({ row, traitType }) {
  return (
    <div style={{
      marginTop: 10, padding: "10px 12px", borderRadius: 8, background: `${C.pheno}10`,
      border: `1px solid ${C.pheno}44`, fontFamily: FONT_DISPLAY, fontSize: 12, color: C.text,
      display: "flex", gap: 18, flexWrap: "wrap",
    }}>
      <strong style={{ color: C.pheno }}>{row.gene}</strong>
      {row.product && <span style={{ color: C.textDim }}>{row.product}</span>}
      {traitType === "binary" ? (
        <>
          <span>present {row.a}/{row.a + row.b} cases · {row.c}/{row.c + row.d} controls</span>
          <span>OR {Number.isFinite(row.or) ? row.or.toFixed(2) : "∞"} [{Number.isFinite(row.lo) ? row.lo.toFixed(2) : "-"}, {Number.isFinite(row.hi) ? row.hi.toFixed(2) : "-"}]</span>
        </>
      ) : (
        <span>Δ mean {Number.isFinite(row.meanDiff) ? row.meanDiff.toFixed(3) : "-"} (t = {Number.isFinite(row.t) ? row.t.toFixed(2) : "-"})</span>
      )}
      <span>p {fmt(row.p)}</span>
      <span>q {fmt(row.q)}</span>
      <span>CMH p {Number.isFinite(row.cmhP) ? fmt(row.cmhP) : "—"}</span>
    </div>
  );
}

const fmt = (v) => (Number.isFinite(v) ? v.toExponential(2) : "-");

function DownloadBtn({ label, onClick }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 6,
      fontSize: 12, color: C.good, border: `1px solid ${C.good}66`, borderRadius: 6, padding: "6px 12px",
    }}>
      <FileText size={12} /> {label}
    </button>
  );
}
