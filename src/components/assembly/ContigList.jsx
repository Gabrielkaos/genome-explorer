import React, { useState } from "react";
import { ChevronDown, Download, RotateCcw } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow } from "../ui/Primitives.jsx";

function downloadFasta(contigs) {
  const text = contigs.map((c) => `>${c.id} length=${c.length}${c.circular ? " circular=true" : ""}\n${(c.seq.match(/.{1,80}/g) || [c.seq]).join("\n")}`).join("\n");
  const blob = new Blob([text], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = "assembly.fasta"; document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

function ReadLayoutDiagram({ contig }) {
  const W = 680, H = 16 + contig.members.length * 9 + 20;
  const scaleX = (W - 20) / contig.length;
  return (
    <svg viewBox={`0 0 ${W} ${H}`} style={{ width: "100%", height: Math.min(H, 260) }}>
      <line x1={10} y1={H - 14} x2={10 + contig.length * scaleX} y2={H - 14} stroke={C.border} strokeWidth={2} />
      {contig.members.map((m, i) => {
        const y = 10 + i * 9;
        const x = 10 + m.contigStart * scaleX;
        const w = Math.max(1.5, (m.contigEnd - m.contigStart) * scaleX);
        return <rect key={i} x={x} y={y} width={w} height={4} rx={1} fill={m.strand === 1 ? C.assembly : C.pheno} opacity={0.8} />;
      })}
      <text x={10} y={H - 4} fontSize="9" fill={C.textFaint} fontFamily={FONT_DISPLAY}>{contig.length.toLocaleString()} bp · {contig.readCount} reads</text>
    </svg>
  );
}

export default function ContigList({ contigs }) {
  const [expanded, setExpanded] = useState(contigs?.[0]?.id ?? null);
  if (!contigs || contigs.length === 0) {
    return <Panel style={{ padding: 16, fontSize: 12.5, color: C.textFaint }}>No contigs were formed from this read set.</Panel>;
  }

  const longest = contigs[0].length; // buildAllContigs already sorts by length desc
  const primary = contigs.filter((c) => c.length >= longest * 0.2);
  const fragments = contigs.filter((c) => c.length < longest * 0.2);

  return (
    <div>
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 10 }}>
        <Eyebrow color={C.assembly}>assembly.fasta ({contigs.length} contig{contigs.length !== 1 ? "s" : ""})</Eyebrow>
        <button onClick={() => downloadFasta(contigs)}
          style={{
            all: "unset", cursor: "pointer", display: "flex", alignItems: "center", gap: 6,
            fontSize: 12, color: C.good, border: `1px solid ${C.good}66`, borderRadius: 2, padding: "6px 12px",
            fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
          }}>
          <Download size={12} /> [ download_fasta ]
        </button>
      </div>

      {fragments.length > 0 && (
        <div style={{ fontSize: 11.5, color: C.textFaint, marginBottom: 10, lineHeight: 1.5 }}>
          Showing {primary.length} primary contig{primary.length !== 1 ? "s" : ""} (≥20% of the longest) plus {fragments.length} smaller fragment{fragments.length !== 1 ? "s" : ""}.
          Under high redundant coverage, small fragments often come from reads that genuinely overlap an already-assembled region but weren't consolidated into it by this simplified graph-matching step — a known limitation, not necessarily separate genomic content. All contigs are included in the FASTA download either way.
        </div>
      )}

      <ContigTable contigs={primary} expanded={expanded} setExpanded={setExpanded} />
      {fragments.length > 0 && (
        <details style={{ marginTop: 10 }}>
          <summary style={{ cursor: "pointer", fontSize: 12, color: C.textDim, fontFamily: FONT_DISPLAY }}>Show {fragments.length} smaller fragment{fragments.length !== 1 ? "s" : ""}</summary>
          <div style={{ marginTop: 8 }}>
            <ContigTable contigs={fragments} expanded={expanded} setExpanded={setExpanded} />
          </div>
        </details>
      )}
    </div>
  );
}

function ContigTable({ contigs, expanded, setExpanded }) {
  return (
    <Panel style={{ padding: 4 }}>
      {contigs.map((c) => {
        const isOpen = expanded === c.id;
        return (
          <div key={c.id} style={{ borderBottom: `1px solid ${C.border}` }}>
            <div onClick={() => setExpanded(isOpen ? null : c.id)}
              style={{ display: "flex", alignItems: "center", gap: 12, padding: "10px 14px", cursor: "pointer" }}>
              <span style={{ fontFamily: FONT_DISPLAY, fontSize: 12.5, color: C.text, width: 90 }}>{c.id}</span>
              <span style={{ fontSize: 11.5, color: C.textDim, width: 100, fontFamily: FONT_DISPLAY }}>{c.length.toLocaleString()} bp</span>
              <span style={{ fontSize: 11.5, color: C.textDim, width: 90, fontFamily: FONT_DISPLAY }}>{c.readCount} reads</span>
              {c.circular && <span style={{ fontSize: 10, color: C.good, border: `1px solid ${C.good}55`, borderRadius: 2, padding: "1px 6px", display: "flex", alignItems: "center", gap: 3, fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.04em" }}><RotateCcw size={9} /> circular</span>}
              <ChevronDown size={14} color={C.textFaint} style={{ marginLeft: "auto", transform: isOpen ? "rotate(180deg)" : "none", transition: "transform .15s" }} />
            </div>
            {isOpen && (
              <div style={{ padding: "4px 14px 16px" }}>
                <div style={{ fontSize: 10.5, color: C.textFaint, marginBottom: 6 }}>Read layout (bar color = strand: teal-purple forward, orange reverse)</div>
                <ReadLayoutDiagram contig={c} />
                <div style={{ fontFamily: FONT_DISPLAY, fontSize: 10.5, lineHeight: 1.7, background: "#05070a", border: `1px solid ${C.border}`, padding: 10, borderRadius: 2, marginTop: 10, maxHeight: 160, overflowY: "auto", wordBreak: "break-all" }}>
                  <div style={{ color: C.assembly, textShadow: `0 0 6px ${C.assembly}33` }}>&gt;{c.id} length={c.length}{c.circular ? " circular=true" : ""}</div>
                  <div style={{ color: C.textDim }}>{(c.seq.match(/.{1,80}/g) || [c.seq]).slice(0, 6).join("\n")}{c.seq.length > 480 ? "\n…" : ""}</div>
                </div>
              </div>
            )}
          </div>
        );
      })}
    </Panel>
  );
}
