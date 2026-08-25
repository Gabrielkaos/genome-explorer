import { useMemo, useRef, useState } from "react";
import { Table2, Play, Loader2, AlertTriangle, FileText } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow, SectionTitle, StatCard, LimitBanner, ExplainBox } from "../ui/Primitives.jsx";
import { useFastqData } from "../../state/FastqDataContext.jsx";
import { downloadBlob } from "../../lib/fastq/exportFastq.js";
import { generateSampleAlignment } from "../../lib/sampleData/generateSampleAlignment.js";
import AlignmentViewer from "./AlignmentViewer.jsx";
import ColumnInspector from "./ColumnInspector.jsx";
import IdentityMatrix from "./IdentityMatrix.jsx";
import VariantTable from "./VariantTable.jsx";

const STAGE_LABEL = {
  distances: "Estimating pairwise k-mer distances",
  tree: "Building UPGMA guide tree",
  align: "Progressive profile-profile alignment",
  stats: "Computing alignment statistics",
  exports: "Building FASTA / Clustal / NEXUS / PHYLIP exports",
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
    .map((t) => t.replace(/[^A-Za-z*-]/g, ""))
    .filter((t) => t.length > 0)
    .map((seq, i) => ({ id: `pasted_${i + 1}`, seq }));
}

/** Uppercase, strip whitespace/digits; unique-ify ids. */
function sanitize(records) {
  const seen = new Map();
  return records
    .map((r) => ({ ...r, seq: r.seq.toUpperCase().replace(/[^A-Z-]/g, "") }))
    .filter((r) => r.seq.length > 0)
    .map((r) => {
      const k = seen.get(r.id) ?? 0;
      seen.set(r.id, k + 1);
      return { ...r, id: k ? `${r.id}_${k + 1}` : r.id };
    });
}

export default function MsaSection({ explainMode }) {
  const { asm, msa } = useFastqData();

  const [source, setSource] = useState("upload"); // upload | paste | contigs | example
  const [uploadRecords, setUploadRecords] = useState(null);
  const [pasteText, setPasteText] = useState("");
  const [match, setMatch] = useState(2);
  const [mismatch, setMismatch] = useState(1);
  const [gapOpen, setGapOpen] = useState(8);
  const [gapExtend, setGapExtend] = useState(1);
  const [kmerSize, setKmerSize] = useState(6);
  const [scrollTarget, setScrollTarget] = useState(null);
  const [colorMode, setColorMode] = useState("nucleotide");
  const fileRef = useRef(null);
  const viewerRef = useRef(null);

  const pasteRecords = useMemo(() => parsePasted(pasteText), [pasteText]);
  const contigRecords = useMemo(
    () => (asm.status === "done" && asm.contigs?.length
      ? asm.contigs.map((c) => ({ id: c.id, seq: c.seq }))
      : null),
    [asm.status, asm.contigs]
  );
  const exampleRecords = useMemo(() => sanitize(generateSampleAlignment(7)), []);

  const sources = {
    upload: uploadRecords,
    paste: pasteRecords.length ? pasteRecords : null,
    contigs: contigRecords,
    example: exampleRecords,
  };
  const activeRecords = sources[source];
  const canRun = !!activeRecords && activeRecords.length >= 2 && msa.status !== "running";
  const totalBases = activeRecords?.reduce((a, r) => a + r.seq.length, 0) ?? 0;

  function handleUpload(file) {
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const recs = sanitize(parseFastaText(String(reader.result)));
        if (!recs.length) throw new Error("No sequences found in file.");
        if (recs.length < 2) throw new Error("Only one sequence found — a multiple alignment needs at least two homologous sequences.");
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
    if (clean.length < 2) return;
    msa.run(clean.map((r) => ({ id: r.id, seq: r.seq })), {
      match,
      mismatch,
      gapOpen,
      gapExtend,
      kmerSize,
    });
    setScrollTarget(null);
  }

  const res = msa.result;
  const sites = useMemo(() => (res ? res.stats.variableSites : []), [res]);

  // Selection is scoped to the result it was made on: a new alignment run
  // automatically clears it without needing a state-reset effect.
  const [selection, setSelection] = useState({ result: null, col: null });
  const selectedCol = selection.result === res ? selection.col : null;

  function selectColumn(col) {
    setSelection({ result: res, col });
    setScrollTarget({ col, nonce: Date.now() });
  }
  function stepSite(dir) {
    const idx = sites.findIndex((s) => s.pos === selectedCol);
    const next = sites[Math.max(0, Math.min(sites.length - 1, (idx < 0 ? 0 : idx) + dir))];
    if (next) selectColumn(next.pos);
  }
  function clearSelection() {
    setSelection({ result: null, col: null });
  }

  return (
    <div>
      <SectionTitle icon={Table2} color={C.phylo} title="Multiple Sequence Alignment"
        subtitle="Progressive alignment with a k-mer guide tree and profile-profile DP — runs entirely in your browser" />

      <LimitBanner>
        This is a genuine ClustalW-style progressive aligner (approximate k-mer distances → UPGMA guide tree →
        profile-profile global alignment with affine gap penalties). It is single-pass: no iterative refinement
        rounds (à la MAFFT FFT-NS-2 or ClustalΩ), no consistency scoring (T-Coffee), and N/ambiguity codes are
        treated as gaps. Expect good results on close-to-moderately diverged sequences (genes, viral genomes,
        plasmids); heavily diverged or rearranged sequences need a server-grade tool.
      </LimitBanner>

      {/* ------------------------- INPUT ------------------------- */}
      <Panel style={{ padding: 18, marginTop: 16, marginBottom: 16 }}>
        <Eyebrow color={C.phylo}>Homologous sequences to align</Eyebrow>
        <div style={{ display: "flex", gap: 10, marginTop: 8, flexWrap: "wrap" }}>
          <SourceCard
            title="Uploaded multi-FASTA" sub={uploadRecords
              ? `${uploadRecords.length} sequence(s)`
              : "the same gene/locus from several strains"} active={source === "upload"}
            enabled={!!uploadRecords} onBrowse={() => fileRef.current.click()} color={C.phylo} />
          <input ref={fileRef} type="file" accept=".fasta,.fa,.fna,.fas,.txt" hidden
            onChange={(e) => e.target.files[0] && handleUpload(e.target.files[0])} />
          <SourceCard
            title="Pasted sequences" sub={pasteRecords.length
              ? `${pasteRecords.length} sequence(s) parsed` : "raw sequence lines or FASTA text"}
            active={source === "paste"} enabled onBrowse={() => setSource("paste")} color={C.raw} />
          <SourceCard
            title="This session's Assembly contigs" sub={contigRecords
              ? `${contigRecords.length} contig(s)` : "run the Assembly section first"}
            active={source === "contigs"} enabled={!!contigRecords}
            onBrowse={() => setSource("contigs")} color={C.assembly} />
          <SourceCard
            title="Synthetic example" sub={`${exampleRecords.length} simulated strains × ~900 bp`}
            active={source === "example"} enabled onBrowse={() => setSource("example")} color={C.qc} />
        </div>

        {source === "paste" && (
          <textarea value={pasteText} onChange={(e) => setPasteText(e.target.value)} rows={5} spellCheck={false}
            placeholder={"Paste FASTA or one sequence per line…\n>strainA\nATGGCGAAG…\n>strainB\nATGGCAAAA…"}
            style={{
              width: "100%", marginTop: 10, background: "#05070a", border: `1px solid ${C.border}`,
              borderRadius: 2, padding: "10px 12px", color: C.text, fontFamily: FONT_DISPLAY, fontSize: 12.5,
            }} />
        )}

        {activeRecords && (
          <div style={{ fontSize: 11.5, color: C.textFaint, marginTop: 8 }}>
            Ready: {activeRecords.length} sequence(s) · {(totalBases / 1000).toFixed(1)} kb total · lengths{" "}
            {Math.min(...activeRecords.map((r) => r.seq.length))}-{Math.max(...activeRecords.map((r) => r.seq.length))} bp
            {activeRecords.length < 2 && " · need at least 2 sequences"}
          </div>
        )}
        {!activeRecords && source === "upload" && (
          <div style={{ fontSize: 11.5, color: C.textFaint, marginTop: 8 }}>
            Upload a multi-FASTA (e.g. one gene exported from several annotations), or try the synthetic example.
          </div>
        )}

        <div style={{ marginTop: 14 }}>
          <Eyebrow color={C.phylo}>Scoring parameters</Eyebrow>
        </div>
        <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(170px, 1fr))", gap: 18, marginTop: 8 }}>
          <Slider label="Match score" value={match} min={1} max={5} step={0.5} fmt={(v) => `+${v}`} onChange={setMatch} color={C.good} />
          <Slider label="Mismatch penalty" value={mismatch} min={0.5} max={5} step={0.5} fmt={(v) => `−${v}`} onChange={setMismatch} color={C.bad} />
          <Slider label="Gap open" value={gapOpen} min={2} max={20} step={1} fmt={(v) => `${v}`} onChange={setGapOpen} color={C.pheno} />
          <Slider label="Gap extend" value={gapExtend} min={0.5} max={4} step={0.5} fmt={(v) => `${v}`} onChange={setGapExtend} color={C.pheno} />
          <Slider label="Guide-tree k-mer size" value={kmerSize} min={4} max={8} step={1} fmt={(v) => `${v}-mers`} onChange={setKmerSize} color={C.phylo} />
        </div>

        <button onClick={handleRun} disabled={!canRun}
          style={{
            all: "unset", cursor: canRun ? "pointer" : "default", display: "inline-flex",
            alignItems: "center", gap: 7, marginTop: 16, padding: "9px 16px", borderRadius: 2,
            border: `1px solid ${C.phylo}66`, color: C.phylo, fontSize: 13,
            fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
            opacity: canRun ? 1 : 0.55,
            textShadow: canRun ? `0 0 8px ${C.phylo}44` : "none",
          }}>
          {msa.status === "running" ? <Loader2 size={14} className="spin" /> : <Play size={14} />}
          {msa.status === "running" ? "[ aligning… ]" : "[ align ]"}
        </button>
        {totalBases > 1_500_000 && canRun && (
          <span style={{ marginLeft: 12, fontSize: 11.5, color: C.qc }}>
            large dataset — may take tens of seconds
          </span>
        )}

        <style>{`.spin { animation: spin 1s linear infinite; } @keyframes spin { from{transform:rotate(0)} to{transform:rotate(360deg)} }`}</style>

        {msa.status === "running" && (
          <div style={{ marginTop: 14 }}>
            <div style={{ fontSize: 12, color: C.textDim, marginBottom: 6 }}>
              {STAGE_LABEL[msa.stage] || "Working…"}{msa.detail ? ` — ${msa.detail}` : ""}
            </div>
            <div style={{ height: 6, background: "#05070a", borderRadius: 1, overflow: "hidden" }}>
              <div style={{ width: `${msa.pct}%`, height: "100%", background: C.phylo, transition: "width .2s" }} />
            </div>
          </div>
        )}
        {msa.error && (
          <div style={{ marginTop: 12, display: "flex", gap: 8, fontSize: 12.5, color: "#f2b3ad" }}>
            <AlertTriangle size={14} /> {msa.error}
          </div>
        )}

        <ExplainBox explainMode={explainMode} color={C.phylo}>
          The guide tree decides WHO merges first (closest pairs make small, reliable profiles early); the
          profile DP then decides WHERE the gaps go. Once a merge is made it is never revisited — which is why
          alignment order matters, and why a bad early merge propagates through everything downstream.
        </ExplainBox>
      </Panel>

      {/* ------------------------- RESULTS ------------------------- */}
      {msa.status === "done" && res && (
        <>
          <StatsRow stats={res.stats} />

          <div style={{ display: "flex", gap: 8, flexWrap: "wrap", margin: "14px 0" }}>
            <span style={{ fontSize: 11, color: C.textFaint, alignSelf: "center", marginRight: 4 }}>
              Exports — feed the aligned matrix straight into IQ-TREE, RAxML, MEGA or any MSA-aware tool:
            </span>
            <DownloadBtn label="FASTA (aligned)" onClick={() => save(res.exports.fasta, "alignment.fasta")} />
            <DownloadBtn label="Clustal (.aln)" onClick={() => save(res.exports.clustal, "alignment.aln")} />
            <DownloadBtn label="NEXUS" onClick={() => save(res.exports.nexus, "alignment.nex")} />
            <DownloadBtn label="PHYLIP" onClick={() => save(res.exports.phylip, "alignment.phylip")} />
            <DownloadBtn label="Variants TSV" onClick={() => save(res.exports.variantsTsv, "variants.tsv")} />
            <DownloadBtn label="Consensus FASTA" onClick={() => save(res.exports.consensusFasta, "consensus.fasta")} />
          </div>

          <Panel style={{ padding: 16 }}>
            <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", flexWrap: "wrap", gap: 8, marginBottom: 10 }}>
              <Eyebrow color={C.phylo}>Aligned matrix</Eyebrow>
              <div style={{ display: "flex", gap: 6, alignItems: "center" }}>
                <LegendChip letter="A" color="#68c98f" /><LegendChip letter="C" color="#5ec8d8" />
                <LegendChip letter="G" color="#e8c15a" /><LegendChip letter="T" color="#ef7fa3" />
                <LegendChip letter="-" color="#39424e" />
                <div style={{ width: 1, height: 18, background: C.border, margin: "0 6px" }} />
                <ModeBtn active={colorMode === "nucleotide"} onClick={() => setColorMode("nucleotide")}>Nucleotide</ModeBtn>
                <ModeBtn active={colorMode === "identity"} onClick={() => setColorMode("identity")}>Identity</ModeBtn>
                <ZoomBtn onClick={() => viewerRef.current?.zoomBy(-1)}>−</ZoomBtn>
                <ZoomBtn onClick={() => viewerRef.current?.zoomBy(1)}>+</ZoomBtn>
              </div>
            </div>
            <AlignmentViewer
              ref={viewerRef}
              ids={res.ids} rows={res.rows} length={res.length}
              consensusCodes={res.stats.consensusCodes}
              selectedCol={selectedCol} onSelectCol={selectColumn}
              scrollTarget={scrollTarget} colorMode={colorMode} />
            {selectedCol != null && (
              <div style={{ marginTop: 14, paddingTop: 12, borderTop: `1px solid ${C.border}` }}>
                <ColumnInspector col={selectedCol} ids={res.ids} rows={res.rows} sites={sites}
                  onClose={clearSelection} onStep={stepSite} />
              </div>
            )}
          </Panel>

          <div className="grid-2-col-asym1" style={{ gap: 16, marginTop: 16 }}>
            <VariantTable sites={sites} totalSites={res.stats.totalVariableSites}
              selectedCol={selectedCol} onSelect={selectColumn} />
            <IdentityMatrix ids={res.ids} matrix={res.identityMatrix} />
          </div>

          <LimitBanner>
            The identity matrix and variable-site list are computed from THIS alignment only — change the input
            set or scoring and every number changes. Parsimony-informative sites (≥2 alleles each seen twice or
            more) are what tree-building actually uses; singleton changes carry no shared-history signal.
          </LimitBanner>
          <ExplainBox explainMode={explainMode} color={C.phylo}>
            The consensus row shows the majority base per column (ties/N-majority shown as N in exports). The
            identity shading makes recombinant-looking tracts pop visually: long runs where one strain switches
            from clustering with its relatives to clustering elsewhere hint at horizontal transfer.
          </ExplainBox>
        </>
      )}

      {msa.status !== "done" && !canRun && (
        <Panel style={{ padding: 16, fontSize: 12.5, color: C.textFaint }}>
          No input yet. Upload a multi-FASTA of homologous sequences, paste them directly, or load the synthetic
          example to see the aligner work end-to-end.
        </Panel>
      )}
    </div>
  );
}

function StatsRow({ stats }) {
  const t = stats.timings;
  return (
    <div style={{ display: "flex", gap: 10, flexWrap: "wrap", marginBottom: 4 }}>
      <StatCard label="Sequences" value={stats.numSeqs.toLocaleString()} color={C.phylo} />
      <StatCard label="Alignment length" value={stats.length.toLocaleString()} unit="col" color={C.phylo} />
      <StatCard label="Mean pair identity" value={(stats.meanPairIdentity * 100).toFixed(1)} unit="%" color={C.good} />
      <StatCard label="Conserved cols" value={stats.conservedColumns.toLocaleString()} color={C.good} />
      <StatCard label="Variable cols" value={stats.variableColumns.toLocaleString()} color={C.qc} />
      <StatCard label="Informative" value={stats.informativeColumns.toLocaleString()} color={C.raw} />
      <StatCard label="Singletons" value={stats.singletonColumns.toLocaleString()} color={C.textDim} />
      <StatCard label="Gap fraction" value={(stats.gapFraction * 100).toFixed(1)} unit="%" color={C.pheno} />
      <StatCard label="Runtime" value={(t.totalMs / 1000).toFixed(2)} unit="s" color={C.textDim} />
    </div>
  );
}

function SourceCard({ title, sub, active, enabled, onBrowse, color }) {
  return (
    <button onClick={() => { if (enabled) onBrowse(); }} style={{
      all: "unset", cursor: enabled ? "pointer" : "default", flex: "1 1 200px", padding: "12px 14px", borderRadius: 2,
      background: active ? C.bgPanel2 : "#05070a",
      border: `1px solid ${active ? color : C.border}`,
      opacity: enabled ? 1 : 0.45,
    }}>
      <div style={{ fontSize: 12.5, color: active ? color : C.textDim, fontFamily: FONT_DISPLAY, textShadow: active ? `0 0 8px ${color}44` : "none" }}>{title}</div>
      <div style={{ fontSize: 11, color: C.textFaint, marginTop: 3, fontFamily: FONT_DISPLAY }}>{sub}</div>
    </button>
  );
}

function Slider({ label, value, min, max, step, fmt, onChange, color }) {
  return (
    <div>
      <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 4 }}>
        <span style={{ fontSize: 12.5, color: C.text, fontFamily: FONT_DISPLAY }}>{label}</span>
        <span style={{ fontFamily: FONT_DISPLAY, color, fontSize: 13, textShadow: `0 0 6px ${color}44` }}>{fmt(value)}</span>
      </div>
      <input type="range" min={min} max={max} step={step} value={value}
        onChange={(e) => onChange(+e.target.value)} style={{ width: "100%", accentColor: color }} />
    </div>
  );
}

function DownloadBtn({ label, onClick }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 6,
      fontSize: 12, color: C.good, border: `1px solid ${C.good}66`, borderRadius: 2, padding: "6px 12px",
      fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
      background: "#05070a",
    }}>
      <FileText size={12} /> [ {label} ]
    </button>
  );
}

function LegendChip({ letter, color }) {
  return (
    <span style={{
      fontFamily: FONT_DISPLAY, fontSize: 10.5, color, border: `1px solid ${color}44`,
      borderRadius: 2, padding: "1px 5px", textShadow: `0 0 6px ${color}44`,
    }}>{letter}</span>
  );
}

function ModeBtn({ active, onClick, children }) {
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

function ZoomBtn({ onClick, children }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", fontSize: 13, lineHeight: 1, padding: "3px 9px", borderRadius: 2,
      border: `1px solid ${C.border}`, color: C.textDim, fontFamily: FONT_DISPLAY,
      background: "#05070a",
    }}>{children}</button>
  );
}

function save(text, fileName) {
  downloadBlob(new Blob([text], { type: "text/plain" }), fileName);
}
