import React from "react";
import { C } from "../../theme.js";
import { StatCard, Eyebrow } from "../ui/Primitives.jsx";
import { formatBases, formatBytes } from "../../lib/fastq/stats.js";

export default function FastqSummaryStats({ dataset, passCount, passPct }) {
  if (!dataset) return null;
  return (
    <div>
      <div style={{ display: "flex", gap: 10, flexWrap: "wrap" }}>
        <StatCard label="Reads" value={dataset.totalReads.toLocaleString()} color={C.raw} />
        <StatCard label="Total bases" value={formatBases(dataset.totalBases)} color={C.raw} />
        <StatCard label="Mean length" value={Math.round(dataset.meanLength).toLocaleString()} unit="bp" color={C.raw} />
        <StatCard label="Median length" value={Math.round(dataset.medianLength).toLocaleString()} unit="bp" color={C.raw} />
        <StatCard label="N50" value={dataset.n50.toLocaleString()} unit="bp" color={C.raw} />
        <StatCard label="N90" value={dataset.n90.toLocaleString()} unit="bp" color={C.raw} />
        <StatCard label="Mean quality" value={dataset.datasetMeanQ.toFixed(1)} unit="Q" color={C.raw} />
        <StatCard label="GC content" value={dataset.gcPercent.toFixed(1)} unit="%" color={C.raw} />
        <StatCard label="Passing filter" value={`${passPct.toFixed(1)}`} unit="%" color={C.good} />
      </div>
      <div style={{ fontSize: 11, color: C.textFaint, marginTop: 8, display: "flex", gap: 16, flexWrap: "wrap" }}>
        <span>File: {dataset.fileName} ({formatBytes(dataset.fileSize)}{dataset.compressed ? ", gzip" : ""})</span>
        <span>Length range: {dataset.minLength.toLocaleString()}–{dataset.maxLength.toLocaleString()} bp</span>
        <span>Parsed in {(dataset.parseTimeMs / 1000).toFixed(1)}s</span>
        <span>Base composition: A {dataset.baseCounts.A.toLocaleString()} · C {dataset.baseCounts.C.toLocaleString()} · G {dataset.baseCounts.G.toLocaleString()} · T {dataset.baseCounts.T.toLocaleString()}{dataset.baseCounts.N ? ` · N ${dataset.baseCounts.N.toLocaleString()}` : ""}</span>
      </div>
      {dataset.warnings?.length > 0 && (
        <div style={{ marginTop: 10, fontSize: 11.5, color: "#f2b3ad" }}>
          <Eyebrow color={C.bad}>Parse warnings ({dataset.warnings.length})</Eyebrow>
          <ul style={{ margin: 0, paddingLeft: 18 }}>
            {dataset.warnings.slice(0, 6).map((w, i) => <li key={i}>{w}</li>)}
          </ul>
        </div>
      )}
      {dataset.parseMode === "buffered" && (
        <div style={{ marginTop: 10, fontSize: 11.5, color: C.qc }}>
          Gzip input: full per-read sequence/quality text is only kept in memory for the first {dataset.previewCap.toLocaleString()} reads. All statistics above (N50, GC%, histograms, mean quality) are computed from the entire file regardless.
        </div>
      )}
    </div>
  );
}
