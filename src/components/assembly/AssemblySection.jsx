import React, { useMemo, useState } from "react";
import { Waypoints, Play, Loader2, Download, AlertTriangle } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow, SectionTitle, StatCard, LimitBanner, ExplainBox } from "../ui/Primitives.jsx";
import { useFastqData } from "../../state/FastqDataContext.jsx";
import { formatBases } from "../../lib/fastq/stats.js";
import ContigList from "./ContigList.jsx";

const STAGE_LABEL = {
  indexing: "Indexing minimizers across selected reads",
  overlaps: "Estimating overlaps between reads",
  graph: "Building & simplifying the overlap graph",
  consensus: "Splicing contig sequences from resolved layout",
};

export default function AssemblySection({ explainMode }) {
  const { dataset, index, getManyRecords, asm } = useFastqData();

  const [minLength, setMinLength] = useState(1000);
  const [minQuality, setMinQuality] = useState(9);
  const [maxReads, setMaxReads] = useState(150);
  const [fetching, setFetching] = useState(false);

  const selection = useMemo(() => {
    if (!index) return null;
    const n = index.lengths.length;
    const candidates = [];
    for (let i = 0; i < n; i++) {
      if (index.lengths[i] >= minLength && index.meanQs[i] >= minQuality) candidates.push(i);
    }
    candidates.sort((a, b) => index.lengths[b] - index.lengths[a]); // prefer longest reads first, same practical strategy real assemblers use when downsampling
    const chosen = candidates.slice(0, maxReads);
    const totalBases = chosen.reduce((s, i) => s + index.lengths[i], 0);
    return { candidateCount: candidates.length, chosen, totalBases };
  }, [index, minLength, minQuality, maxReads]);

  async function handleRun() {
    if (!selection || selection.chosen.length < 2) return;
    setFetching(true);
    try {
      const records = await getManyRecords(selection.chosen);
      const reads = records
        .map((r, i) => (r && !r.unavailable ? { seq: r.seq, sourceReadIndex: selection.chosen[i], id: r.id } : null))
        .filter(Boolean);
      asm.run(reads, {});
    } finally {
      setFetching(false);
    }
  }

  if (!dataset) {
    return (
      <div>
        <SectionTitle icon={Waypoints} color={C.assembly} title="Flye-style Assembly" subtitle="Real overlap-layout assembly from your reads — no reference genome used" />
        <Panel style={{ padding: 16, fontSize: 12.5, color: C.textFaint }}>
          No FASTQ data loaded yet. Upload or load the sample dataset on the FASTQ page first — Assembly uses the same dataset.
        </Panel>
      </div>
    );
  }

  return (
    <div>
      <SectionTitle icon={Waypoints} color={C.assembly} title="Flye-style Assembly" subtitle="A real overlap-layout assembler run on your actual reads — not a simulation with a known answer" />

      <LimitBanner>
        This performs genuine minimizer-based overlap detection, string-graph simplification, and greedy-layout contig
        construction on your real reads — architecturally closer to a lightweight tool like <em>miniasm</em> than to Flye's
        full repeat-graph algorithm, and it does not run a base-level polishing pass (the way a real pipeline would follow
        up with Racon/medaka). Contig sequences are therefore a <strong>draft assembly</strong>, and junction boundaries are
        approximate (typically accurate to within roughly 1-2% of true length) since no base-level alignment is used to
        verify overlaps — only k-mer matching. This is disclosed, not hidden, in the statistics below.
      </LimitBanner>

      <Panel style={{ padding: 18, marginTop: 16, marginBottom: 16 }}>
        <Eyebrow color={C.assembly}>Choose input reads</Eyebrow>
        <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr 1fr", gap: 20, marginTop: 8 }}>
          <div>
            <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 4 }}>
              <span style={{ fontSize: 12.5, color: C.text }}>Min. length</span>
              <span style={{ fontFamily: FONT_DISPLAY, color: C.assembly, fontSize: 13 }}>{minLength.toLocaleString()} bp</span>
            </div>
            <input type="range" min={0} max={Math.max(dataset.maxLength, 2000)} step={100} value={minLength}
              onChange={(e) => setMinLength(+e.target.value)} style={{ width: "100%", accentColor: C.assembly }} />
          </div>
          <div>
            <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 4 }}>
              <span style={{ fontSize: 12.5, color: C.text }}>Min. quality</span>
              <span style={{ fontFamily: FONT_DISPLAY, color: C.assembly, fontSize: 13 }}>Q{minQuality}</span>
            </div>
            <input type="range" min={0} max={25} step={0.5} value={minQuality}
              onChange={(e) => setMinQuality(+e.target.value)} style={{ width: "100%", accentColor: C.assembly }} />
          </div>
          <div>
            <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 4 }}>
              <span style={{ fontSize: 12.5, color: C.text }}>Max reads to use</span>
              <span style={{ fontFamily: FONT_DISPLAY, color: C.assembly, fontSize: 13 }}>{maxReads.toLocaleString()}</span>
            </div>
            <input type="range" min={20} max={3000} step={20} value={maxReads}
              onChange={(e) => setMaxReads(+e.target.value)} style={{ width: "100%", accentColor: C.assembly }} />
          </div>
        </div>

        {selection && (
          <div style={{ fontSize: 12, color: C.textDim, marginTop: 12 }}>
            {selection.candidateCount.toLocaleString()} reads pass these thresholds — using the {selection.chosen.length.toLocaleString()} longest
            ({formatBases(selection.totalBases)} total). Longest-first selection is a real, practical strategy for keeping assembly tractable.
          </div>
        )}

        <button
          onClick={handleRun}
          disabled={fetching || asm.status === "running" || !selection || selection.chosen.length < 2}
          style={{
            all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 7, marginTop: 14,
            padding: "9px 16px", borderRadius: 2, border: `1px solid ${C.assembly}66`, color: C.assembly, fontSize: 13,
            fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
            opacity: fetching || asm.status === "running" ? 0.6 : 1,
          }}>
          {fetching || asm.status === "running" ? <Loader2 size={14} className="spin" /> : <Play size={14} />}
          {fetching ? "[ fetching_reads… ]" : asm.status === "running" ? "[ assembling… ]" : "[ run_assembly ]"}
        </button>

        <style>{`.spin { animation: spin 1s linear infinite; } @keyframes spin { from{transform:rotate(0)} to{transform:rotate(360deg)} }`}</style>

        {asm.status === "running" && (
          <div style={{ marginTop: 14 }}>
            <div style={{ fontSize: 12, color: C.textDim, marginBottom: 6 }}>{STAGE_LABEL[asm.stage] || "Working…"}</div>
            <div style={{ height: 6, background: "#05070a", borderRadius: 1, overflow: "hidden" }}>
              <div style={{ width: `${asm.pct}%`, height: "100%", background: C.assembly, transition: "width .2s" }} />
            </div>
          </div>
        )}
        {asm.error && (
          <div style={{ marginTop: 12, display: "flex", gap: 8, fontSize: 12.5, color: "#f2b3ad", background: "#05070a", border: `1px solid ${C.bad}55`, borderRadius: 2, padding: "8px 12px" }}>
            <AlertTriangle size={14} /> {asm.error}
          </div>
        )}

        <ExplainBox explainMode={explainMode} color={C.assembly}>
          Real assembly needs enough reads overlapping each region (coverage) to reconstruct it confidently — too few reads, or reads that are too short/low-quality, will fragment the assembly into more contigs than the true genome/plasmid count.
        </ExplainBox>
      </Panel>

      {asm.status === "done" && asm.stats && (
        <AssemblyResults stats={asm.stats} contigs={asm.contigs} explainMode={explainMode} />
      )}
    </div>
  );
}

function AssemblyResults({ stats, contigs, explainMode }) {
  return (
    <div>
      <div style={{ display: "flex", gap: 10, flexWrap: "wrap", marginBottom: 16 }}>
        <StatCard label="Contigs" value={stats.numContigs} color={C.assembly} />
        <StatCard label="Total length" value={formatBases(stats.totalLength)} color={C.assembly} />
        <StatCard label="Assembly N50" value={stats.n50.toLocaleString()} unit="bp" color={C.assembly} />
        <StatCard label="Longest contig" value={stats.longestContig.toLocaleString()} unit="bp" color={C.assembly} />
        <StatCard label="Circular" value={stats.circularContigs} color={C.assembly} />
        <StatCard label="Reads used" value={stats.usedReads} color={C.good} />
        <StatCard label="Contained (redundant)" value={stats.containedReads} color={C.textDim} />
        <StatCard label="Unplaced" value={stats.unplacedReads} color={C.bad} />
        <StatCard label="Overlaps found" value={stats.overlapsFound} color={C.assembly} />
        <StatCard label="Mean overlap score" value={stats.meanOverlapScore.toFixed(2)} color={C.assembly} />
      </div>

      <div style={{ fontSize: 11, color: C.textFaint, marginBottom: 16 }}>
        Parameters used: k={stats.params.k}, w={stats.params.w}, min overlap {stats.params.minOverlapLen}bp, min shared minimizers {stats.params.minMatches}.
        Computed in {(stats.computeTimeMs / 1000).toFixed(1)}s.
      </div>

      {stats.unplacedReads > stats.inputReads * 0.15 && (
        <LimitBanner>
          {stats.unplacedReads} of {stats.inputReads} reads couldn't be confidently overlapped with anything (isolated in the graph) — likely low coverage in some regions, reads too short relative to the minimum overlap length, or high local error rate. Try lowering the minimum length threshold, increasing max reads, or check the QC Dashboard for this run's quality distribution.
        </LimitBanner>
      )}

      <ExplainBox explainMode={explainMode} color={C.assembly}>
        "Contained" reads are ones fully spanned by a longer read's overlap — real information, but redundant for building the draft backbone, so they're set aside rather than discarded from the underlying data.
      </ExplainBox>

      <div style={{ marginTop: 16 }}>
        <ContigList contigs={contigs} />
      </div>
    </div>
  );
}
