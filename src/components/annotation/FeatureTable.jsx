import React, { useMemo, useState } from "react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow } from "../ui/Primitives.jsx";

const RENDER_CAP = 400;

/**
 * Sortable/filterable table of predicted CDS features. Filtering is
 * client-side over the full result set; rendering is capped so a large
 * annotation doesn't stall the DOM.
 */
export default function FeatureTable({ genes, selectedGeneId, onSelectGene }) {
  const [search, setSearch] = useState("");
  const [sortKey, setSortKey] = useState("locus");
  const [sortDir, setSortDir] = useState(1);
  const [onlyTm, setOnlyTm] = useState(false);
  const [onlySignal, setOnlySignal] = useState(false);
  const [onlyPartial, setOnlyPartial] = useState(false);

  const filtered = useMemo(() => {
    const q = search.trim().toLowerCase();
    let out = genes.filter((g) => {
      if (q && !(g.locusTag.toLowerCase().includes(q) || (g.protSeq || "").toLowerCase().includes(q.slice(0, 12)))) return false;
      if (onlyTm && !g.tmSegments?.length) return false;
      if (onlySignal && !g.signal) return false;
      if (onlyPartial && !g.partial) return false;
      return true;
    });
    const keyFn = {
      locus: (g) => g.locusTag,
      start: (g) => g.start,
      len: (g) => g.lengthAa,
      score: (g) => g.score,
      rbs: (g) => g.rbsScore,
      gc: (g) => g.gcContent,
    }[sortKey];
    if (keyFn) out.sort((a, b) => {
      const va = keyFn(a), vb = keyFn(b);
      return (typeof va === "string" ? va.localeCompare(vb) : va - vb) * sortDir;
    });
    return out;
  }, [genes, search, sortKey, sortDir, onlyTm, onlySignal, onlyPartial]);

  const shown = filtered.slice(0, RENDER_CAP);
  const chip = (active, onClick, children) => (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", fontSize: 11, padding: "3px 9px", borderRadius: 6,
      background: active ? `${C.annotation}22` : "transparent",
      border: `1px solid ${active ? C.annotation : C.border}`,
      color: active ? C.annotation : C.textDim,
    }}>{children}</button>
  );
  const th = (key, label, width) => (
    <th onClick={() => { setSortKey(key); setSortDir(sortKey === key ? -sortDir : 1); }}
      style={{ textAlign: "left", padding: "5px 8px", color: C.textFaint, fontSize: 10.5, cursor: "pointer", userSelect: "none", width }}>
      {label}{sortKey === key ? (sortDir === 1 ? " ▲" : " ▼") : ""}
    </th>
  );

  return (
    <Panel style={{ padding: 14 }}>
      <div style={{ display: "flex", gap: 10, alignItems: "center", flexWrap: "wrap", marginBottom: 10 }}>
        <Eyebrow color={C.annotation}>Feature table</Eyebrow>
        <input value={search} onChange={(e) => setSearch(e.target.value)} placeholder="filter by locus tag or protein sequence…"
          style={{ flex: 1, minWidth: 160, background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 6, padding: "5px 10px", color: C.text, fontSize: 12 }} />
        {chip(onlyTm, () => setOnlyTm(!onlyTm), "TM helices")}
        {chip(onlySignal, () => setOnlySignal(!onlySignal), "signal peptide")}
        {chip(onlyPartial, () => setOnlyPartial(!onlyPartial), "partial")}
      </div>
      <div style={{ maxHeight: 380, overflowY: "auto" }}>
        <table style={{ borderCollapse: "collapse", width: "100%", fontSize: 11.5 }}>
          <thead>
            <tr>
              {th("locus", "Locus tag", 110)}
              {th("start", "Start")}
              {th("len", "Length (aa)", 80)}
              {th("score", "Score", 60)}
              {th("rbs", "RBS", 55)}
              {th("gc", "GC%", 55)}
              <th style={{ textAlign: "left", padding: "5px 8px", color: C.textFaint, fontSize: 10.5 }}>Flags</th>
            </tr>
          </thead>
          <tbody>
            {shown.map((g) => (
              <tr key={g.locusTag} onClick={() => onSelectGene(g.locusTag)}
                style={{ cursor: "pointer", background: g.locusTag === selectedGeneId ? `${C.annotation}16` : "transparent" }}>
                <td style={{ padding: "4px 8px", fontFamily: FONT_DISPLAY, color: C.text, borderTop: `1px solid ${C.border}` }}>{g.locusTag}</td>
                <td style={{ padding: "4px 8px", fontFamily: FONT_DISPLAY, color: C.textDim, borderTop: `1px solid ${C.border}` }}>
                  {(g.start + 1).toLocaleString()}{g.strand === "-" ? " ←" : " →"}
                </td>
                <td style={{ padding: "4px 8px", color: C.textDim, borderTop: `1px solid ${C.border}` }}>{g.lengthAa.toLocaleString()}</td>
                <td style={{ padding: "4px 8px", fontFamily: FONT_DISPLAY, color: C.annotation, borderTop: `1px solid ${C.border}` }}>{g.score.toFixed(2)}</td>
                <td style={{ padding: "4px 8px", fontFamily: FONT_DISPLAY, borderTop: `1px solid ${C.border}`, color: g.rbsScore >= 0.6 ? C.good : C.textFaint }}>{g.rbsScore.toFixed(2)}</td>
                <td style={{ padding: "4px 8px", fontFamily: FONT_DISPLAY, color: C.qc, borderTop: `1px solid ${C.border}` }}>{g.gcContent.toFixed(1)}</td>
                <td style={{ padding: "4px 8px", borderTop: `1px solid ${C.border}`, whiteSpace: "nowrap" }}>
                  {g.tmSegments?.length > 0 && <Flag color={C.pheno}>TM×{g.tmSegments.length}</Flag>}
                  {g.signal && <Flag color={C.good}>signal</Flag>}
                  {g.partial && <Flag color={C.bad}>partial {g.partial}</Flag>}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
        {filtered.length > RENDER_CAP && (
          <div style={{ fontSize: 11, color: C.textFaint, padding: "8px" }}>
            Showing first {RENDER_CAP.toLocaleString()} of {filtered.length.toLocaleString()} features — refine the filters to see more.
          </div>
        )}
        {filtered.length === 0 && <div style={{ fontSize: 12, color: C.textFaint, padding: 12 }}>No features match.</div>}
      </div>
    </Panel>
  );
}

function Flag({ color, children }) {
  return (
    <span style={{ fontSize: 9.5, fontFamily: FONT_DISPLAY, color, border: `1px solid ${color}55`, borderRadius: 4, padding: "1px 5px", marginRight: 4 }}>{children}</span>
  );
}
