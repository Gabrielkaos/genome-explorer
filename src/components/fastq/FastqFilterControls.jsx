import React, { useState } from "react";
import { Download, Loader2 } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow, ExplainBox } from "../ui/Primitives.jsx";
import { buildFilteredFastqBlob, downloadBlob } from "../../lib/fastq/exportFastq.js";

export default function FastqFilterControls({
  dataset, index, file, previewRecords,
  qualityThreshold, setQualityThreshold,
  minLength, setMinLength,
  maxLength, setMaxLength,
  passCount, passPct, explainMode,
}) {
  const [exporting, setExporting] = useState(false);
  const canExportFull = dataset?.parseMode === "offset";

  async function handleExport() {
    setExporting(true);
    try {
      const passingIndices = [];
      for (let i = 0; i < index.lengths.length; i++) {
        const len = index.lengths[i], q = index.meanQs[i];
        if (q >= qualityThreshold && len >= minLength && (!maxLength || len <= maxLength)) passingIndices.push(i);
      }
      const result = buildFilteredFastqBlob({ file, dataset, index, previewRecords, passingIndices });
      const skipped = result?.skipped;
      downloadBlob(result, `filtered_${dataset.fileName.replace(/\.gz$/i, "")}`);
      if (skipped) {
        // eslint-disable-next-line no-alert
        alert(`${skipped} passing read(s) beyond the in-memory preview window for this gzip file could not be included. Export from an uncompressed .fastq for a complete filtered file.`);
      }
    } finally {
      setExporting(false);
    }
  }

  return (
    <Panel style={{ padding: 18 }}>
      <Eyebrow color={C.qc}>Filter thresholds</Eyebrow>

      <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 20, marginTop: 8 }}>
        <div>
          <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 4 }}>
            <span style={{ fontSize: 12.5, color: C.text }}>Min. mean quality</span>
            <span style={{ fontFamily: FONT_DISPLAY, color: C.qc, fontSize: 14 }}>Q{qualityThreshold}</span>
          </div>
          <input type="range" min={0} max={30} step={0.5} value={qualityThreshold}
            onChange={(e) => setQualityThreshold(+e.target.value)} style={{ width: "100%", accentColor: C.qc }} />
        </div>
        <div>
          <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 4 }}>
            <span style={{ fontSize: 12.5, color: C.text }}>Min. length</span>
            <span style={{ fontFamily: FONT_DISPLAY, color: C.qc, fontSize: 14 }}>{minLength.toLocaleString()} bp</span>
          </div>
          <input type="range" min={0} max={Math.max(dataset?.maxLength || 5000, 1000)} step={50} value={minLength}
            onChange={(e) => setMinLength(+e.target.value)} style={{ width: "100%", accentColor: C.qc }} />
        </div>
      </div>

      <div style={{ display: "flex", gap: 14, alignItems: "center", marginTop: 14, flexWrap: "wrap" }}>
        <label style={{ fontSize: 12, color: C.textDim, display: "flex", alignItems: "center", gap: 6 }}>
          Max length (optional)
          <input type="number" placeholder="none" value={maxLength || ""}
            onChange={(e) => setMaxLength(e.target.value ? +e.target.value : null)}
            style={{ width: 90, background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 6, color: C.text, padding: "4px 8px", fontFamily: FONT_DISPLAY, fontSize: 12 }} />
        </label>

        <div style={{ marginLeft: "auto", fontSize: 13, color: C.text }}>
          <span style={{ color: C.good, fontFamily: FONT_DISPLAY }}>{passCount.toLocaleString()}</span> / {dataset?.totalReads.toLocaleString()} reads pass
          <span style={{ color: C.textFaint }}> ({passPct.toFixed(1)}%)</span>
        </div>

        <button onClick={handleExport} disabled={exporting || !dataset}
          style={{
            all: "unset", cursor: exporting ? "default" : "pointer", display: "flex", alignItems: "center", gap: 6,
            fontSize: 12.5, color: C.good, border: `1px solid ${C.good}66`, borderRadius: 6, padding: "7px 13px",
            opacity: exporting ? 0.6 : 1,
          }}>
          {exporting ? <Loader2 size={13} className="spin" /> : <Download size={13} />}
          Export filtered FASTQ
        </button>
      </div>

      {!canExportFull && (
        <div style={{ fontSize: 11, color: C.textFaint, marginTop: 8 }}>
          This file was gzip-compressed, so export is limited to reads within the first {dataset?.previewCap?.toLocaleString()} (in-memory preview). Upload an uncompressed .fastq for a complete filtered export of arbitrarily large files.
        </div>
      )}

      <ExplainBox explainMode={explainMode} color={C.qc}>
        This filters whole reads out of the dataset - it does not touch individual bases within a read that passes. That's a deliberate, separate step from QC (which only measures).
      </ExplainBox>
      <style>{`.spin { animation: spin 1s linear infinite; } @keyframes spin { from{transform:rotate(0)} to{transform:rotate(360deg)} }`}</style>
    </Panel>
  );
}
