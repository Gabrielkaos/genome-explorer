import React, { useCallback, useRef, useState } from "react";
import { Upload, FileUp, X, Sparkles, Loader2 } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel } from "../ui/Primitives.jsx";
import { formatBytes } from "../../lib/fastq/stats.js";
import { generateSampleFastqFile } from "../../lib/sampleData/generateSampleFastq.js";

export default function FastqUploader({ status, progress, onFile, onCancel, error }) {
  const [dragOver, setDragOver] = useState(false);
  const inputRef = useRef(null);

  const handleFiles = useCallback((files) => {
    const file = files?.[0];
    if (!file) return;
    if (!/\.(fastq|fq|fastq\.gz|fq\.gz|txt)$/i.test(file.name)) {
      // Not a hard block - some pipelines use non-standard extensions -
      // but warn visually via the uploader's own state if needed later.
    }
    onFile(file);
  }, [onFile]);

  const isParsing = status === "parsing";
  const pct = progress.totalBytes ? Math.min(100, (progress.bytesProcessed / progress.totalBytes) * 100) : 0;

  return (
    <Panel style={{ padding: 20 }}>
      <div
        onDragOver={(e) => { e.preventDefault(); setDragOver(true); }}
        onDragLeave={() => setDragOver(false)}
        onDrop={(e) => { e.preventDefault(); setDragOver(false); handleFiles(e.dataTransfer.files); }}
        onClick={() => !isParsing && inputRef.current?.click()}
        style={{
          border: `1px solid ${dragOver ? C.raw : C.border}`, borderRadius: 2, padding: "26px 20px",
          textAlign: "center", cursor: isParsing ? "default" : "pointer",
          background: dragOver ? `${C.raw}08` : "transparent", transition: "border-color .15s, background .15s",
        }}
      >
        <input ref={inputRef} type="file" accept=".fastq,.fq,.fastq.gz,.fq.gz,.txt,.gz" hidden
          onChange={(e) => handleFiles(e.target.files)} />

        {!isParsing ? (
          <>
            <Upload size={22} color={C.raw} style={{ marginBottom: 8 }} />
            <div style={{ fontSize: 14, color: C.text }}>Drop a .fastq or .fastq.gz file here, or click to browse</div>
            <div style={{ fontSize: 12, color: C.textFaint, marginTop: 4 }}>
              Parsed entirely in your browser via a background worker — files never leave your machine. Phred+33 encoding assumed.
            </div>
          </>
        ) : (
          <>
            <Loader2 size={22} color={C.raw} className="spin" style={{ marginBottom: 8 }} />
            <div style={{ fontSize: 14, color: C.text }}>Parsing… {progress.readsProcessed.toLocaleString()} reads processed</div>
            <div style={{ width: "100%", maxWidth: 420, height: 6, background: "#05070a", borderRadius: 1, margin: "12px auto 6px", overflow: "hidden" }}>
              <div style={{ width: `${pct}%`, height: "100%", background: C.raw, transition: "width .2s" }} />
            </div>
            <div style={{ fontSize: 11.5, color: C.textFaint, fontFamily: FONT_DISPLAY }}>
              {formatBytes(progress.bytesProcessed)} / {formatBytes(progress.totalBytes)}
            </div>
            <button onClick={(e) => { e.stopPropagation(); onCancel(); }}
              style={{ all: "unset", cursor: "pointer", marginTop: 12, fontSize: 12, color: C.bad, display: "inline-flex", alignItems: "center", gap: 5, fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em" }}>
              <X size={13} /> [ cancel ]
            </button>
          </>
        )}
      </div>

      {error && (
        <div style={{ marginTop: 12, padding: "10px 12px", borderRadius: 2, fontSize: 12.5, background: "#05070a", border: `1px solid ${C.bad}55`, color: "#f2b3ad" }}>
          {error}
        </div>
      )}

      <div style={{ display: "flex", alignItems: "center", gap: 10, marginTop: 14 }}>
        <div style={{ flex: 1, height: 1, background: C.border }} />
        <span style={{ fontSize: 11, color: C.textFaint }}>or</span>
        <div style={{ flex: 1, height: 1, background: C.border }} />
      </div>
      <button
        disabled={isParsing}
        onClick={() => onFile(generateSampleFastqFile())}
        style={{
          all: "unset", cursor: isParsing ? "default" : "pointer", display: "flex", alignItems: "center", justifyContent: "center", gap: 7,
          marginTop: 14, padding: "9px 14px", borderRadius: 2, border: `1px solid ${C.raw}55`, color: C.raw, fontSize: 12.5,
          fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
          opacity: isParsing ? 0.5 : 1,
        }}>
        <Sparkles size={13} /> [ load_sample ]
      </button>

      <style>{`
        .spin { animation: spin 1s linear infinite; }
        @keyframes spin { from { transform: rotate(0deg); } to { transform: rotate(360deg); } }
      `}</style>
    </Panel>
  );
}
