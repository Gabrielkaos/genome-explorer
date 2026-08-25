import { useMemo, useRef, useState } from "react";
import { GitBranch, Play, Loader2, AlertTriangle } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow, SectionTitle, StatCard, LimitBanner, ExplainBox } from "../ui/Primitives.jsx";
import { useFastqData } from "../../state/FastqDataContext.jsx";
import { generateSampleTreeData } from "../../lib/sampleData/generateSampleTreeData.js";
import { MODELS, GAP_MODES } from "../../lib/phylo/distances.js";
import TreeView from "./TreeView.jsx";
import DistanceMatrix from "./DistanceMatrix.jsx";
import NewickPanel from "./NewickPanel.jsx";
import {
  SourceCard, Segmented, Slider, NumberField, ParamNote, DownloadBtn,
} from "./controls.jsx";
import { save, buildNexusTrees, buildMatrixTsv } from "../../lib/phylo/exportPhylo.js";

const STAGE_LABEL = {
  encode: "Encoding alignment",
  distance: "Estimating pairwise distances",
  tree: "Building tree topology",
  bootstrap: "Bootstrap replicates",
  export: "Serializing Newick",
};

function parseFastaText(text) {
  const records = [];
  let cur = null;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (cur) records.push(cur);
      const header = line.slice(1).trim();
      cur = { id: header.split(/\s+/)[0] || `seq_${records.length + 1}`, seq: "" };
    } else if (cur) cur.seq += line.trim();
  }
  if (cur) records.push(cur);
  return records;
}

function parsePasted(text) {
  if (text.includes(">")) return parseFastaText(text);
  return text.split(/\s+/)
    .map((t) => t.replace(/[^A-Za-z*?.~-]/g, ""))
    .filter((t) => t.length > 0)
    .map((seq, i) => ({ id: `pasted_${i + 1}`, seq }));
}

/** Uppercase; strip whitespace/digits; keep gaps & ambiguity codes; unique-ify ids. */
function sanitize(records) {
  const seen = new Map();
  return records
    .map((r) => ({ ...r, seq: r.seq.toUpperCase().replace(/[^A-Z*?.~-]/g, "") }))
    .filter((r) => r.seq.length > 0)
    .map((r) => {
      const k = seen.get(r.id) ?? 0;
      seen.set(r.id, k + 1);
      return { ...r, id: k ? `${r.id}_${k + 1}` : r.id };
    });
}

export default function PhyloSection({ explainMode }) {
  const { msa, phyl } = useFastqData();

  const [source, setSource] = useState("example"); // example | msa | upload | paste
  const [uploadRecords, setUploadRecords] = useState(null);
  const [pasteText, setPasteText] = useState("");
  const [method, setMethod] = useState("nj");
  const [model, setModel] = useState("jc");
  const [gapMode, setGapMode] = useState("pairwise");
  const [rooting, setRooting] = useState("midpoint");
  const [bootstrap, setBootstrap] = useState(100);
  const [seed, setSeed] = useState(42);
  // Viewer display state.
  const [layout, setLayout] = useState("rect");
  const [mode, setMode] = useState("phylogram");
  const [showSupports, setShowSupports] = useState(true);
  const [imported, setImported] = useState({ result: null, tree: null });
  const fileRef = useRef(null);

  const exampleRecords = useMemo(() => sanitize(generateSampleTreeData(11)), []);
  const msaRecords = useMemo(() => (
    msa.status === "done" && msa.result
      ? msa.result.ids.map((id, i) => ({ id, seq: msa.result.rows[i] }))
      : null
  ), [msa.status, msa.result]);
  const pasteRecords = useMemo(() => parsePasted(pasteText), [pasteText]);

  const sources = {
    msa: msaRecords,
    upload: uploadRecords,
    paste: pasteRecords.length ? pasteRecords : null,
    example: exampleRecords,
  };
  const activeRecords = sources[source];

  const lenInfo = useMemo(() => {
    if (!activeRecords?.length) return null;
    const lens = activeRecords.map((r) => r.seq.length);
    const min = Math.min(...lens), max = Math.max(...lens);
    return { min, max, aligned: min === max };
  }, [activeRecords]);

  const proteinish = useMemo(() => {
    if (!activeRecords?.length) return false;
    const sample = activeRecords.slice(0, 5).map((r) => r.seq.slice(0, 400)).join("");
    let nt = 0;
    for (const ch of sample) if ("ACGTUN-*?.".includes(ch)) nt++;
    return sample.length > 50 && nt / sample.length < 0.85;
  }, [activeRecords]);

  const canRun = !!activeRecords && activeRecords.length >= 3 && lenInfo?.aligned && phyl.status !== "running";

  function handleUpload(file) {
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const recs = sanitize(parseFastaText(String(reader.result)));
        if (!recs.length) throw new Error("No sequences found in file.");
        if (recs.length < 3) throw new Error(`Only ${recs.length} sequence(s) found - tree inference needs at least three.`);
        setUploadRecords(recs);
        setSource("upload");
      } catch (err) {
        alert(`Could not use this FASTA: ${err.message}`);
      }
    };
    reader.readAsText(file);
  }

  function handleRun() {
    const clean = sanitize(activeRecords);
    if (clean.length < 3) return;
    phyl.run(clean.map((r) => ({ id: r.id, seq: r.seq })), {
      method, model, gapMode,
      rooting: method === "upgma" ? "midpoint" : rooting,
      bootstrap, seed,
    });
    setImported({ result: null, tree: null });
  }

  const res = phyl.result;
  const displayedTree = imported.result === res ? imported.tree : res?.tree ?? null;

  return (
    <div>
      <SectionTitle icon={GitBranch} color={C.phylo}
        title="Phylogenetic Inference"
        subtitle="Neighbor-Joining & UPGMA over model-corrected distances, with nonparametric bootstraps - runs entirely in your browser" />

      <LimitBanner>
        These are genuine distance-matrix methods: JC69/K2P-corrected pairwise distances feed Neighbor-Joining or
        UPGMA, and support values come from real column-resampling bootstrap replicates. This is not maximum
        likelihood - there is no substitution-model selection (ModelFinder), no ML branch-length optimization, and
        bootstrap percentages are not posterior probabilities. UPGMA additionally assumes a molecular clock - prefer
        NJ unless you specifically need a clock tree. Very short internodes may stay unresolved under any method.
      </LimitBanner>

      {/* ------------------------------ INPUT ------------------------------ */}
      <Panel style={{ padding: 18, marginTop: 16, marginBottom: 16 }}>
        <Eyebrow color={C.phylo}>Aligned sequences</Eyebrow>
        <div style={{ fontSize: 11.5, color: C.textFaint, marginTop: 2, marginBottom: 10 }}>
          Distances require an ALREADY-ALIGNED matrix (equal-length rows). Align first in the Alignment section,
          then infer the tree on its output here.
        </div>
        <div style={{ display: "flex", gap: 10, flexWrap: "wrap" }}>
          <SourceCard title="This session's Alignment" color={C.raw}
            sub={msaRecords ? `${msaRecords.length} seqs × ${msa.result.length} cols` : "run the Alignment section first"}
            active={source === "msa"} enabled={!!msaRecords} onBrowse={() => setSource("msa")} />
          <SourceCard title="Uploaded aligned FASTA" color={C.phylo}
            sub={uploadRecords ? `${uploadRecords.length} sequence(s)` : ".fasta / .aln / plain text"}
            active={source === "upload"} enabled onBrowse={() => fileRef.current.click()} />
          <input ref={fileRef} type="file" accept=".fasta,.fa,.fna,.fas,.aln,.txt,.nex" hidden
            onChange={(e) => e.target.files[0] && handleUpload(e.target.files[0])} />
          <SourceCard title="Pasted sequences" color={C.qc}
            sub={pasteRecords.length ? `${pasteRecords.length} sequence(s) parsed` : "FASTA or one sequence per line"}
            active={source === "paste"} enabled onBrowse={() => setSource("paste")} />
          <SourceCard title="Synthetic example" color={C.assembly}
            sub={`${exampleRecords.length} strains × 1200 bp · true topology known`}
            active={source === "example"} enabled onBrowse={() => setSource("example")} />
        </div>

        {source === "paste" && (
          <textarea value={pasteText} onChange={(e) => setPasteText(e.target.value)} rows={5} spellCheck={false}
            placeholder={"Paste an ALIGNED multi-FASTA…\n>strainA\nATGGCGAAG…\n>strainB\nATGGCAAAA…"}
            style={{
              width: "100%", marginTop: 10, background: "#05070a", border: `1px solid ${C.border}`,
              borderRadius: 2, padding: "10px 12px", color: C.text, fontFamily: FONT_DISPLAY, fontSize: 12.5,
            }} />
        )}

        {activeRecords && lenInfo && (
          <div style={{ fontSize: 11.5, marginTop: 8, color: C.textFaint }}>
            Ready: {activeRecords.length} sequence(s) ×{" "}
            {lenInfo.aligned ? `${lenInfo.min} bp` : `UNEQUAL LENGTHS (${lenInfo.min}-${lenInfo.max} bp)`}
            {!lenInfo.aligned && <span style={{ color: C.bad }}> — align these first; distance methods need equal-length rows.</span>}
            {activeRecords.length < 3 && <span style={{ color: C.bad }}> — need at least 3 sequences.</span>}
          </div>
        )}
        {proteinish && (
          <div style={{ fontSize: 11.5, color: C.qc, marginTop: 6 }}>
            Heads-up: this input looks like PROTEIN sequences. The models here are nucleotide models - results will be meaningless.
          </div>
        )}
        {!activeRecords && source === "upload" && (
          <div style={{ fontSize: 11.5, color: C.textFaint, marginTop: 8 }}>
            Upload an aligned multi-FASTA, or load the synthetic example to explore the tool.
          </div>
        )}

        {/* --------------------------- PARAMETERS --------------------------- */}
        <div style={{ marginTop: 16 }}><Eyebrow color={C.phylo}>Inference parameters</Eyebrow></div>
        <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(230px, 1fr))", gap: 18, marginTop: 8 }}>
          <Segmented label="Tree-building method" value={method} onChange={setMethod}
            options={[
              { v: "nj", label: "Neighbor-Joining", hint: "clock-free; usually the right default" },
              { v: "upgma", label: "UPGMA", hint: "molecular-clock assumption; rooted" },
            ]} />
          <Segmented label="Substitution model" value={model} onChange={setModel}
            options={[
              { v: "jc", label: MODELS.jc.short, hint: "corrects multiple hits, equal rates" },
              { v: "k2p", label: MODELS.k2p.short, hint: "separate transition/transversion rates" },
              { v: "p", label: MODELS.p.short, hint: "uncorrected; fine for very close taxa" },
            ]} />
          <Segmented label="Gap handling" value={gapMode} onChange={setGapMode}
            options={[
              { v: "pairwise", label: GAP_MODES.pairwise.label, hint: "skip gapped columns per pair" },
              { v: "complete", label: GAP_MODES.complete.label, hint: "drop any column gapped anywhere" },
            ]} />
          {method === "nj" ? (
            <Segmented label="Rooting" value={rooting} onChange={setRooting}
              options={[
                { v: "midpoint", label: "Midpoint root", hint: "root halfway along the longest path" },
                { v: "none", label: "As inferred", hint: "leave unrooted (trifurcation)" },
              ]} />
          ) : (
            <ParamNote>UPGMA yields a rooted ultrametric tree by construction - every tip is equidistant from the root.</ParamNote>
          )}
          <Slider label="Bootstrap replicates" value={bootstrap} min={0} max={500} step={25}
            fmt={(v) => (v === 0 ? "off" : `${v} reps`)} onChange={setBootstrap} color={C.phylo} />
          <NumberField label="Bootstrap RNG seed" value={seed} onChange={setSeed} hint="same seed = same resampling" />
        </div>

        <div style={{ display: "flex", alignItems: "center", gap: 12, marginTop: 16, flexWrap: "wrap" }}>
          <button onClick={handleRun} disabled={!canRun}
            style={{
              all: "unset", cursor: canRun ? "pointer" : "default", display: "inline-flex",
              alignItems: "center", gap: 7, padding: "9px 16px", borderRadius: 2,
              border: `1px solid ${C.phylo}66`, color: C.phylo, fontSize: 13,
              fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
              opacity: canRun ? 1 : 0.55,
              textShadow: canRun ? `0 0 8px ${C.phylo}44` : "none",
            }}>
            {phyl.status === "running" ? <Loader2 size={14} className="spin" /> : <Play size={14} />}
            {phyl.status === "running" ? "[ inferring… ]" : "[ infer_tree ]"}
          </button>
          {phyl.status === "running" && (
            <button onClick={phyl.cancel}
              style={{
                all: "unset", cursor: "pointer", fontSize: 12, color: C.bad,
                border: `1px solid ${C.bad}66`, borderRadius: 2, padding: "7px 12px",
                fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
                background: "#05070a",
              }}>
              [ cancel ]
            </button>
          )}
          {canRun && bootstrap > 100 && (
            <span style={{ fontSize: 11.5, color: C.qc }}>this rebuilds the full pipeline {bootstrap} times</span>
          )}
        </div>

        <style>{`.spin { animation: spin 1s linear infinite; } @keyframes spin { from{transform:rotate(0)} to{transform:rotate(360deg)} }`}</style>

        {phyl.status === "running" && (
          <div style={{ marginTop: 14 }}>
            <div style={{ fontSize: 12, color: C.textDim, marginBottom: 6 }}>
              {STAGE_LABEL[phyl.stage] || "Working…"}{phyl.detail ? ` — ${phyl.detail}` : ""}
            </div>
            <div style={{ height: 6, background: "#05070a", borderRadius: 1, overflow: "hidden" }}>
              <div style={{ width: `${phyl.pct}%`, height: "100%", background: C.phylo, transition: "width .2s" }} />
            </div>
          </div>
        )}
        {phyl.error && (
          <div style={{ marginTop: 12, display: "flex", gap: 8, fontSize: 12.5, color: "#f2b3ad" }}>
            <AlertTriangle size={14} /> {phyl.error}
          </div>
        )}

        <ExplainBox explainMode={explainMode} color={C.phylo}>
          A phylogeny is a hypothesis about descent: tips that share a recent branch share a recent common ancestor.
          The bootstrap answers a different question - "if evolution replayed this alignment, would I recover the
          same clade again?" - which is why a beautiful tree with 30% supports deserves distrust, and why strain
          relationships should be read together with branch lengths, not just topology.
        </ExplainBox>
      </Panel>

      {/* ------------------------------ RESULTS ------------------------------ */}
      {phyl.status === "done" && res && displayedTree && (
        <>
          <StatsRow res={res} />

          <Panel style={{ padding: 16, marginBottom: 16 }}>
            <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", flexWrap: "wrap", gap: 8, marginBottom: 10 }}>
              <Eyebrow color={C.phylo}>{imported.result === res ? "Imported tree viewer" : "Inferred tree"}</Eyebrow>
              <div style={{ display: "flex", gap: 6, alignItems: "center", flexWrap: "wrap" }}>
                <ViewBtn active={layout === "rect"} onClick={() => setLayout("rect")}>Rectangular</ViewBtn>
                <ViewBtn active={layout === "radial"} onClick={() => setLayout("radial")}>Radial</ViewBtn>
                <div style={{ width: 1, height: 16, background: C.border }} />
                <ViewBtn active={mode === "phylogram"} onClick={() => setMode("phylogram")}>Phylogram</ViewBtn>
                <ViewBtn active={mode === "cladogram"} onClick={() => setMode("cladogram")}>Cladogram</ViewBtn>
                <div style={{ width: 1, height: 16, background: C.border }} />
                <ViewBtn active={showSupports} onClick={() => setShowSupports((s) => !s)}>Support values</ViewBtn>
                {imported.result !== res && <RootChip res={res} />}
              </div>
            </div>
            <TreeView
              key={`${res.ms}-${res.paramsEcho.method}`}
              tree={displayedTree}
              layout={layout}
              mode={mode}
              showSupports={showSupports && (res.bootstrapReps > 0 || imported.result === res)}
            />
            <LimitBanner>
              {res.bootstrapReps > 0
                ? `Support dots come from ${res.bootstrapReps} bootstrap replicates (seed ${res.params.seed}): they measure how often each split reappears when alignment columns are resampled - splits below ~70% are unreliable. Branch lengths are in expected substitutions/site under ${MODELS[res.paramsEcho.model].label}.`
                : "Bootstrapping was disabled for this run - no support values are shown. Turn replicates above 0 to quantify confidence per split."}
            </LimitBanner>
          </Panel>

          <div style={{ display: "grid", gridTemplateColumns: "minmax(0, 1.15fr) minmax(0, 1fr)", gap: 16, alignItems: "start" }}>
            <DistanceMatrix ids={res.ids} matrix={res.matrix}
              transitions={res.transitions} transversions={res.transversions} />

            <div style={{ display: "grid", gap: 16 }}>
              <NewickPanel
                newick={res.newick}
                importActive={imported.result === res}
                onImport={(tree) => setImported({ result: res, tree })}
                onClearImport={() => setImported({ result: null, tree: null })}
              />

              <Panel style={{ padding: 16 }}>
                <Eyebrow color={C.phylo}>Exports</Eyebrow>
                <div style={{ display: "flex", gap: 8, flexWrap: "wrap", marginTop: 8 }}>
                  <DownloadBtn label="Newick (.tre)" onClick={() => save(res.newick, "tree.tre")} />
                  <DownloadBtn label="NEXUS" onClick={() => save(buildNexusTrees(res.newick), "tree.nex")} />
                  <DownloadBtn label="Distance matrix TSV" onClick={() => save(buildMatrixTsv(res.ids, res.matrix), "distances.tsv")} />
                </div>
                <div style={{ fontSize: 10.5, color: C.textFaint, marginTop: 8 }}>
                  Newick drops straight into IQ-TREE's tree-viewer steps, iTOL annotation, or MEGA.
                </div>
              </Panel>
            </div>
          </div>

          <ExplainBox explainMode={explainMode} color={C.phylo}>
            Read the distance matrix against the tree: NJ builds exactly the topology whose implied path lengths best
            match those numbers, so pairs sharing a shallow split should also sit close in the matrix. When they
            disagree - e.g. one taxon looks equidistant to everything - suspect saturation (too many multiple hits)
            or a different ortholog sneaked into the alignment.
          </ExplainBox>
        </>
      )}

      {phyl.status !== "done" && !canRun && (
        <Panel style={{ padding: 16, fontSize: 12.5, color: C.textFaint }}>
          No inference yet. Load at least three ALIGNED sequences (session alignment, upload, paste, or the synthetic
          example), pick a method, and press Infer tree.
        </Panel>
      )}
    </div>
  );
}

/* ------------------------------- pieces ------------------------------- */

function StatsRow({ res }) {
  const s = res.summary;
  const pair = (p) => (p ? `${res.ids[p[0]]} ↔ ${res.ids[p[1]]}` : "-");
  return (
    <div style={{ display: "flex", gap: 10, flexWrap: "wrap", marginBottom: 14 }}>
      <StatCard label="Taxa" value={res.numTaxa} color={C.phylo} />
      <StatCard label="Columns used" value={res.columnsUsed.toLocaleString()} unit={`/ ${res.columns}`} color={C.phylo} />
      <StatCard label="Mean distance" value={s.meanDistance.toFixed(4)} unit="subs/site" color={C.raw} />
      <PairCard text={pair(s.furthestPair)} d={s.furthestPair ? res.matrix[s.furthestPair[0] * res.numTaxa + s.furthestPair[1]] : 0} />
      {Number.isFinite(s.tiTvRatio) && <StatCard label="Ti/Tv ratio" value={s.tiTvRatio.toFixed(2)} color={C.qc} />}
      {s.saturatedPairs > 0 && <StatCard label="Saturated pairs" value={s.saturatedPairs} color={C.bad} />}
      <StatCard label="Bootstrap reps" value={res.bootstrapReps || "off"} color={C.good} />
      <StatCard label="Runtime" value={(res.ms / 1000).toFixed(2)} unit="s" color={C.textDim} />
    </div>
  );
}

function PairCard({ text, d }) {
  return (
    <div style={{ background: C.bgPanel, border: `1px solid ${C.border}`, borderRadius: 2, padding: "10px 14px" }}>
      <div style={{ fontSize: 10.5, color: C.textFaint, textTransform: "uppercase", letterSpacing: "0.08em", fontFamily: FONT_DISPLAY }}>Most divergent pair</div>
      <div style={{ fontFamily: FONT_DISPLAY, fontSize: 13, color: C.bad, marginTop: 3, textShadow: `0 0 8px ${C.bad}44` }}>
        {text} <span style={{ fontSize: 10, color: C.textDim }}>d = {d.toFixed(3)}</span>
      </div>
    </div>
  );
}

function ViewBtn({ active, onClick, children }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", fontSize: 11, padding: "3px 9px", borderRadius: 2,
      background: active ? C.bgPanel2 : "transparent",
      border: `1px solid ${active ? C.phylo : C.border}`,
      color: active ? C.phylo : C.textDim,
      fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
      textShadow: active ? `0 0 6px ${C.phylo}44` : "none",
    }}>[ {children} ]</button>
  );
}

function RootChip({ res }) {
  const p = res.paramsEcho;
  const label = p.method === "upgma"
    ? "rooted (UPGMA)"
    : p.rooting === "midpoint" ? "midpoint-rooted" : "unrooted";
  const color = p.method === "upgma" ? C.assembly : p.rooting === "midpoint" ? C.raw : C.qc;
  return (
    <span style={{
      fontSize: 10.5, color, border: `1px solid ${color}44`, borderRadius: 2, padding: "2px 8px",
      fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
      textShadow: `0 0 6px ${color}44`,
    }}>[ {label} ]</span>
  );
}
