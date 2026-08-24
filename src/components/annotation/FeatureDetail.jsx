import React from "react";
import { Copy } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow } from "../ui/Primitives.jsx";

function copy(text) {
  navigator.clipboard?.writeText(text).catch(() => {});
}

function Row({ label, value, mono = false }) {
  return (
    <div style={{ display: "flex", justifyContent: "space-between", gap: 10, fontSize: 12.5, padding: "3px 0" }}>
      <span style={{ color: C.textFaint }}>{label}</span>
      <span style={{ color: C.text, fontFamily: mono ? FONT_DISPLAY : undefined, textAlign: "right" }}>{value}</span>
    </div>
  );
}

export default function FeatureDetail({ gene }) {
  if (!gene) {
    return (
      <Panel style={{ padding: 18, fontSize: 12.5, color: C.textFaint }}>
        Select a feature on the map or in the table to inspect it.
      </Panel>
    );
  }

  const { start, end } = { start: gene.start + 1, end: gene.end };

  // paint the protein sequence with signal/TM spans highlighted
  const spans = [];
  if (gene.signal) spans.push({ a: 1, b: gene.signal.cleavageAfter, color: `${C.good}55`, label: "signal peptide" });
  for (const tm of gene.tmSegments) spans.push({ a: tm.start, b: tm.end, color: `${C.pheno}55`, label: "TM helix" });

  const aaCells = Array.from(gene.protSeq).map((aa, i) => {
    const span = spans.find((s) => i + 1 >= s.a && i + 1 <= s.b);
    return (
      <span key={i}
        title={span ? span.label : `residue ${i + 1}`}
        style={{
          background: span?.color,
          color: aa === "X" ? C.textFaint : C.text,
          borderRadius: 2,
        }}>{aa}</span>
    );
  });

  return (
    <Panel style={{ padding: 18 }}>
      <Eyebrow color={C.annotation}>CDS · {gene.strand === "+" ? "forward" : "reverse"} strand</Eyebrow>
      <div style={{ fontFamily: FONT_DISPLAY, fontSize: 16, color: C.text }}>{gene.locusTag}</div>
      <div style={{ fontSize: 12.5, color: C.textDim, margin: "4px 0 12px" }}>{gene.product}{gene.partial ? ` — truncated at contig edge (${gene.partial})` : ""}</div>

      <Row label="Location" value={`${gene.contigId}:${start.toLocaleString()}–${end.toLocaleString()} (${(end - start + 1).toLocaleString()} bp)`} mono />
      <Row label="Start codon / RBS" value={`${gene.startCodon ?? "—"} · RBS ${gene.rbsScore.toFixed(2)}${gene.rbsSpacer !== null && gene.rbsSpacer !== undefined ? ` (spacer ${gene.rbsSpacer} nt)` : ""}`} />
      <Row label="Coding score" value={`${gene.score.toFixed(2)} total · coding term ${gene.codingTerm.toFixed(2)} · mean hexamer LLR ${gene.meanHexLR.toFixed(2)}`} mono />
      <Row label="Protein size" value={`${gene.lengthAa.toLocaleString()} aa`} mono />
      <Row label="Molecular weight" value={`${(gene.mw / 1000).toFixed(1)} kDa`} mono />
      <Row label="Theoretical pI" value={gene.pi ?? "—"} mono />
      <Row label="GRAVY hydropathy" value={gene.gravy.toFixed(3)} mono />
      <Row label="Aromaticity" value={`${(gene.aromaticity * 100).toFixed(1)}% F+W+Y`} mono />
      <Row label="Gene GC%" value={`${gene.gcContent.toFixed(1)}%`} mono />
      <Row label="Localization hint" value={gene.localization} />

      {gene.tmSegments.length > 0 && (
        <div style={{ marginTop: 8, fontSize: 11.5, color: C.textDim }}>
          Transmembrane helices (heuristic): {gene.tmSegments.map((t) => `residues ${t.start}–${t.end}`).join(", ")}
        </div>
      )}
      {gene.signal && (
        <div style={{ marginTop: 6, fontSize: 11.5, color: C.good }}>
          Signal peptide predicted (score {gene.signal.score}) — cleavage after residue {gene.signal.cleavageAfter}. {gene.signal.note}
        </div>
      )}

      <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", marginTop: 14, marginBottom: 4 }}>
        <Eyebrow color={C.annotation}>Protein ({gene.lengthAa} aa)</Eyebrow>
        <button onClick={() => copy(`>${gene.locusTag}\n${gene.protSeq}`)} style={copyBtn}><Copy size={11} /> copy</button>
      </div>
      <div style={{ background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 6, padding: 8, maxHeight: 120, overflowY: "auto", wordBreak: "break-all", lineHeight: 1.8, fontSize: 11.5, fontFamily: FONT_DISPLAY }}>
        {aaCells}
      </div>
      <div style={{ fontSize: 10, color: C.textFaint, marginTop: 3 }}>highlighted = predicted signal peptide (green) / TM helices (orange)</div>

      <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", marginTop: 12, marginBottom: 4 }}>
        <Eyebrow color={C.annotation}>Nucleotide ({gene.lengthNt.toLocaleString()} nt)</Eyebrow>
        <button onClick={() => copy(`>${gene.locusTag}\n${gene.dnaSeq}`)} style={copyBtn}><Copy size={11} /> copy</button>
      </div>
      <div style={{ background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 6, padding: 8, maxHeight: 90, overflowY: "auto", wordBreak: "break-all", lineHeight: 1.7, fontSize: 10.5, fontFamily: FONT_DISPLAY, color: C.textDim }}>
        {gene.dnaSeq}
      </div>
    </Panel>
  );
}

const copyBtn = {
  all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 5,
  fontSize: 11, color: C.textDim, border: `1px solid ${C.border}`, borderRadius: 6, padding: "3px 9px",
};
