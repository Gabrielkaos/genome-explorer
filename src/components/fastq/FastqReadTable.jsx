import { useEffect, useMemo, useRef, useState } from "react";
import { Search, ChevronDown, ArrowUpDown } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow } from "../ui/Primitives.jsx";

const ROW_HEIGHT = 32;
const VIEWPORT_HEIGHT = 420;
const BUFFER_ROWS = 8;

function QualityStrip({ qual, width = 200 }) {
  return (
    <div style={{ display: "flex", height: 12, borderRadius: 2, overflow: "hidden", width, flexShrink: 0 }}>
      {Array.from(qual).map((c, i) => {
        const q = c.charCodeAt(0) - 33;
        const t = Math.max(0, Math.min(1, q / 30));
        return <div key={i} style={{ flex: 1, background: `hsl(${t * 130}, 70%, ${38 + t * 18}%)` }} />;
      })}
    </div>
  );
}

export default function FastqReadTable({ index, qualityThreshold, minLength, maxLength, getRecord }) {
  const [search, setSearch] = useState("");
  const [sortKey, setSortKey] = useState("index");
  const [sortDir, setSortDir] = useState("asc");
  const [expanded, setExpanded] = useState(null);
  const [expandedRecord, setExpandedRecord] = useState(null);
  const [loadingRecord, setLoadingRecord] = useState(false);
  const [scrollTop, setScrollTop] = useState(0);
  const scrollRef = useRef(null);

  const order = useMemo(() => {
    if (!index) return [];
    const n = index.lengths.length;
    let arr = [];
    const q = search.trim().toLowerCase();
    if (q) {
      for (let i = 0; i < n; i++) if (index.ids[i].toLowerCase().includes(q)) arr.push(i);
    } else {
      arr = Array.from({ length: n }, (_, i) => i);
    }
    if (sortKey !== "index") {
      const key = sortKey === "length" ? index.lengths : index.meanQs;
      arr.sort((a, b) => (key[a] - key[b]) * (sortDir === "asc" ? 1 : -1));
    } else if (sortDir === "desc") {
      arr.reverse();
    }
    return arr;
  }, [index, search, sortKey, sortDir]);

  useEffect(() => {
    if (expanded === null) return;
    let cancelled = false;
    // eslint-disable-next-line react-hooks/set-state-in-effect
    setLoadingRecord(true);
    getRecord(expanded).then((rec) => { if (!cancelled) { setExpandedRecord(rec); setLoadingRecord(false); } });
    return () => { cancelled = true; };
  }, [expanded, getRecord]);

  if (!index) return null;

  const total = order.length;
  const startIdx = Math.max(0, Math.floor(scrollTop / ROW_HEIGHT) - BUFFER_ROWS);
  const endIdx = Math.min(total, Math.ceil((scrollTop + VIEWPORT_HEIGHT) / ROW_HEIGHT) + BUFFER_ROWS);
  const visible = order.slice(startIdx, endIdx);

  function toggleSort(key) {
    if (sortKey === key) setSortDir((d) => (d === "asc" ? "desc" : "asc"));
    else { setSortKey(key); setSortDir("desc"); }
    setScrollTop(0);
    if (scrollRef.current) scrollRef.current.scrollTop = 0;
  }

  function handleSearch(e) {
    setSearch(e.target.value);
    setScrollTop(0);
    if (scrollRef.current) scrollRef.current.scrollTop = 0;
  }

  const expandedOrderIdx = expanded !== null ? order.indexOf(expanded) : -1;
  const EXPANDED_HEIGHT = 140;

  return (
    <Panel style={{ padding: 16 }}>
      <div style={{ display: "flex", alignItems: "center", gap: 10, marginBottom: 10, flexWrap: "wrap" }}>
        <Eyebrow color={C.raw}>Reads ({total.toLocaleString()}{search ? ` matching "${search}"` : ""})</Eyebrow>
        <div style={{ marginLeft: "auto", display: "flex", alignItems: "center", gap: 8 }}>
          <Search size={13} color={C.textFaint} />
          <input value={search} onChange={handleSearch} placeholder="Search read ID…"
            style={{ background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 2, padding: "5px 10px", fontSize: 12, color: C.text, width: 180, fontFamily: FONT_DISPLAY }} />
        </div>
      </div>

      <div style={{ display: "flex", gap: 12, padding: "0 4px 6px", fontSize: 10.5, color: C.textFaint, borderBottom: `1px solid ${C.border}`, fontFamily: FONT_DISPLAY, letterSpacing: "0.06em" }}>
        <span style={{ width: 90 }}>READ ID</span>
        <button onClick={() => toggleSort("length")} style={sortBtnStyle(sortKey === "length")}>LENGTH <ArrowUpDown size={10} /></button>
        <button onClick={() => toggleSort("quality")} style={{ ...sortBtnStyle(sortKey === "quality"), width: 60 }}>QUALITY <ArrowUpDown size={10} /></button>
        <span style={{ flex: 1 }}>QUALITY PROFILE (first 60 bp)</span>
        <span style={{ width: 50 }}>STATUS</span>
      </div>

      <div ref={scrollRef} onScroll={(e) => setScrollTop(e.currentTarget.scrollTop)}
        style={{ height: VIEWPORT_HEIGHT, overflowY: "auto", position: "relative" }}>
        <div style={{ height: total * ROW_HEIGHT + (expanded !== null ? EXPANDED_HEIGHT : 0), position: "relative" }}>
          {visible.map((i, vi) => {
            const rowOrderIdx = startIdx + vi;
            const rowTop = rowOrderIdx * ROW_HEIGHT + (expandedOrderIdx !== -1 && rowOrderIdx > expandedOrderIdx ? EXPANDED_HEIGHT : 0);
            const len = index.lengths[i], q = index.meanQs[i], id = index.ids[i];
            const pass = q >= qualityThreshold && len >= minLength && (!maxLength || len <= maxLength);
            const isOpen = expanded === i;
            return (
              <div key={i} style={{ position: "absolute", top: rowTop, left: 0, right: 0 }}>
                <div onClick={() => {
                  const willOpen = !isOpen;
                  setExpanded(willOpen ? i : null);
                  if (!willOpen) setExpandedRecord(null);
                }}
                  style={{ display: "flex", alignItems: "center", gap: 12, height: ROW_HEIGHT, padding: "0 4px", cursor: "pointer", borderBottom: `1px solid ${C.border}`, fontFamily: FONT_DISPLAY }}>
                  <span style={{ width: 8, height: 8, borderRadius: 2, background: pass ? C.good : C.bad, flexShrink: 0 }} />
                  <span style={{ fontFamily: FONT_DISPLAY, fontSize: 11.5, color: C.text, width: 82, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>{id}</span>
                  <span style={{ fontSize: 11, color: C.textDim, width: 70, fontFamily: FONT_DISPLAY }}>{len.toLocaleString()} bp</span>
                  <span style={{ fontSize: 11, color: C.textDim, width: 46, fontFamily: FONT_DISPLAY }}>Q{q.toFixed(1)}</span>
                  <span style={{ fontSize: 10.5, color: pass ? C.good : C.bad, marginLeft: "auto", fontFamily: FONT_DISPLAY }}>{pass ? "PASS" : "FAIL"}</span>
                  <ChevronDown size={13} color={C.textFaint} style={{ transform: isOpen ? "rotate(180deg)" : "none", transition: "transform .15s" }} />
                </div>
                {isOpen && (
                  <div style={{ padding: "10px 12px", background: "#05070a", fontSize: 11.5, fontFamily: FONT_DISPLAY, color: C.textDim }}>
                    {loadingRecord && "Loading record…"}
                    {!loadingRecord && expandedRecord?.unavailable && (
                      <span style={{ color: C.qc, fontFamily: "inherit" }}>
                        Full sequence not available: this read is beyond the in-memory preview window for this gzip file.
                      </span>
                    )}
                    {!loadingRecord && expandedRecord && !expandedRecord.unavailable && (
                      <>
                        <div style={{ color: C.raw }}>@{expandedRecord.id}{expandedRecord.desc ? ` ${expandedRecord.desc}` : ""}</div>
                        <div style={{ color: C.text, wordBreak: "break-all", margin: "4px 0" }}>{expandedRecord.seq.slice(0, 240)}{expandedRecord.seq.length > 240 ? "…" : ""}</div>
                        <QualityStrip qual={expandedRecord.qual.slice(0, 60)} width={260} />
                      </>
                    )}
                  </div>
                )}
              </div>
            );
          })}
        </div>
      </div>
    </Panel>
  );
}

function sortBtnStyle(active) {
  return {
    all: "unset", cursor: "pointer", display: "flex", alignItems: "center", gap: 3, width: 70,
    fontFamily: FONT_DISPLAY, letterSpacing: "0.06em",
    color: active ? C.qc : "inherit",
  };
}
