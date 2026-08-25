import { useEffect, useMemo, useState } from "react";
import { Search as SearchIcon, X } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Eyebrow } from "../ui/Primitives.jsx";
import { findMotif, searchFeatures } from "../../lib/explorer/search.js";

/**
 * Two search modes over the loaded contig:
 *  - feature text: substring across locus tag / product / type
 *  - DNA motif: IUPAC-aware (RYSWKMBDHVN), both strands
 * Clicking a hit recenters the view; motif hits also get temporary markers.
 */
export default function SearchPanel({ contig, features, onJump, onMarkHits }) {
  const [mode, setMode] = useState("feature"); // feature | motif
  const [query, setQuery] = useState("");

  const featureHits = useMemo(() => (
    mode === "feature" ? searchFeatures(features, query).slice(0, 60) : []
  ), [features, query, mode]);

  const motifHits = useMemo(() => {
    if (mode !== "motif") return [];
    const p = query.trim().toUpperCase();
    if (!p || p.length < 3 || p.length > 40) return [];
    return findMotif(contig.seq, p).slice(0, 120);
  }, [contig.seq, query, mode]);

  const hits = mode === "feature" ? featureHits : motifHits;

  // Persistently highlight current motif hits on the map.
  useEffect(() => {
    onMarkHits(mode === "motif" ? motifHits : []);
  }, [mode, motifHits, onMarkHits]);

  function jump(pos0, pos1) {
    const span = Math.max(200, (pos1 - pos0) * 4);
    onJump({ start: Math.max(0, (pos0 + pos1) / 2 - span / 2), end: Math.min(contig.seq.length, (pos0 + pos1) / 2 + span / 2) });
  }

  return (
    <div>
      <Eyebrow color={C.raw}>Search</Eyebrow>
      <div style={{ display: "flex", gap: 6, marginTop: 6, marginBottom: 8 }}>
        <Chip active={mode === "feature"} onClick={() => setMode("feature")}>Features</Chip>
        <Chip active={mode === "motif"} onClick={() => setMode("motif")}>DNA motif</Chip>
      </div>
      <div style={{ position: "relative" }}>
        <SearchIcon size={13} style={{ position: "absolute", left: 9, top: 9, color: C.textFaint }} />
        <input value={query} onChange={(e) => setQuery(e.target.value)}
          spellCheck={false}
          placeholder={mode === "feature" ? "locus tag or product…" : "IUPAC motif e.g. GCRWTTAT…"}
          style={{
            width: "100%", background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 2,
            padding: "7px 26px 7px 28px", color: C.text, fontFamily: FONT_DISPLAY, fontSize: 12.5,
          }} />
        {query && (
          <button onClick={() => setQuery("")} style={{ all: "unset", cursor: "pointer", position: "absolute", right: 8, top: 9, color: C.textFaint }}>
            <X size={12} />
          </button>
        )}
      </div>

      <div style={{ fontSize: 10.5, color: C.textFaint, margin: "7px 0" }}>
        {mode === "motif"
          ? `${motifHits.length} hit(s)${query.trim().length > 0 && query.trim().length < 3 ? " — pattern too short" : ""} · both strands`
          : `${featureHits.length} match(es)`}
      </div>

      <div style={{ maxHeight: 260, overflowY: "auto" }}>
        {hits.map((h, i) => (
          <button key={`${h.locusTag ?? h.start0}_${i}`}
            onClick={() => mode === "feature"
              ? jump(h.start, h.end)
              : jump(h.start0, h.end)}
            style={{
              all: "unset", cursor: "pointer", display: "block", width: "100%", boxSizing: "border-box",
              padding: "5px 8px", borderRadius: 2, borderBottom: `1px solid ${C.border}`, fontSize: 11.5,
              background: "#05070a",
            }}>
            {mode === "feature" ? (
              <>
                <span style={{ fontFamily: FONT_DISPLAY, color: C.text }}>{h.locusTag}</span>
                <span style={{ color: C.textFaint }}> {h.type} </span>
                <span style={{ color: C.textDim, display: "block", whiteSpace: "nowrap", overflow: "hidden", textOverflow: "ellipsis" }}>{h.product}</span>
              </>
            ) : (
              <span style={{ fontFamily: FONT_DISPLAY, color: C.text }}>
                {(h.start0 + 1).toLocaleString()}–{h.end.toLocaleString()}{" "}
                <span style={{ color: h.strand === "+" ? C.good : C.bad }}>{h.strand === "+" ? "→" : "←"}</span>{" "}
                <span style={{ color: C.qc, textShadow: `0 0 6px ${C.qc}44` }}>{h.match}</span>
              </span>
            )}
          </button>
        ))}
      </div>
    </div>
  );
}

function Chip({ children, active, onClick }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", fontSize: 11, padding: "3px 10px", borderRadius: 2,
      background: active ? `${C.raw}22` : "#05070a",
      border: `1px solid ${active ? C.raw : C.border}`,
      color: active ? C.raw : C.textDim,
      fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
      textShadow: active ? `0 0 6px ${C.raw}44` : "none",
    }}>{children}</button>
  );
}
