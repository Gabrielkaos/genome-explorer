import { useMemo } from "react";
import { ChevronLeft, ChevronRight, X } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Eyebrow } from "../ui/Primitives.jsx";
import { BASE_COLORS, BASE_LETTERS } from "./colors.js";

/**
 * Detail panel for one selected alignment column: allele composition,
 * which strains carry the minority alleles, and prev/next stepping through
 * variable sites.
 */
export default function ColumnInspector({ col, ids, rows, sites, onClose, onStep }) {
  const info = useMemo(() => {
    const counts = [0, 0, 0, 0, 0];
    for (const row of rows) {
      const ch = row[col];
      counts[ch === "-" ? 4 : "ACGT".indexOf(ch)]++;
    }
    let majCode = 0, majCount = -1;
    for (let a = 0; a < 5; a++) if (counts[a] > majCount) { majCount = counts[a]; majCode = a; }
    const differing = [];
    for (let r = 0; r < ids.length; r++) {
      const ch = rows[r][col];
      const code = ch === "-" ? 4 : "ACGT".indexOf(ch);
      if (code !== majCode) differing.push({ id: ids[r], char: ch });
    }
    return { counts, majCode, majCount, differing };
  }, [rows, ids, col]);

  // ungapped position within the first sequence (if that base is real)
  let ungapped = null;
  const row0 = rows[0];
  if (row0 && row0[col] !== "-") {
    ungapped = 0;
    for (let c = 0; c <= col; c++) if (row0[c] !== "-") ungapped++;
  }

  const siteIdx = sites.findIndex((s) => s.pos === col);
  const total = rows.length;

  return (
    <div>
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 6 }}>
        <Eyebrow color={C.phylo}>
          Column {col + 1}{siteIdx >= 0 ? ` · variant ${siteIdx + 1}/${sites.length}` : " · conserved"}
          {ungapped != null && ` · pos ${ungapped} in ${ids[0].slice(0, 18)}`}
        </Eyebrow>
        <div style={{ display: "flex", gap: 4 }}>
          <IconMini disabled={siteIdx <= 0} onClick={() => onStep(-1)}><ChevronLeft size={13} /></IconMini>
          <IconMini disabled={siteIdx < 0 || siteIdx >= sites.length - 1} onClick={() => onStep(1)}><ChevronRight size={13} /></IconMini>
          <IconMini onClick={onClose}><X size={13} /></IconMini>
        </div>
      </div>

      <div style={{ display: "flex", gap: 10, flexWrap: "wrap" }}>
        {[0, 1, 2, 3, 4].map((a) => (
          <div key={a} style={{ minWidth: 92 }}>
            <div style={{ fontSize: 11, color: a === 4 ? C.textFaint : BASE_COLORS[a], fontFamily: FONT_DISPLAY }}>
              {BASE_LETTERS[a]} · {info.counts[a]}
            </div>
            <div style={{ height: 5, background: "#05070a", borderRadius: 3, marginTop: 3, overflow: "hidden" }}>
              <div style={{
                width: `${(info.counts[a] / Math.max(1, total)) * 100}%`, height: "100%",
                background: a === 4 ? "#39424e" : BASE_COLORS[a],
              }} />
            </div>
          </div>
        ))}
      </div>

      {info.differing.length > 0 && (
        <div style={{ marginTop: 12, fontSize: 12, color: C.textDim }}>
          Minority alleles in{" "}
          <span style={{ color: C.text }}>{info.differing.length}</span> of {total} sequences:{" "}
          {info.differing.slice(0, 14).map((d, i) => (
            <span key={d.id} style={{ marginLeft: i ? 7 : 2, fontFamily: FONT_DISPLAY, color: d.char === "-" ? C.textFaint : C.text }}>
              {d.id.slice(0, 22)}:<span style={{ color: d.char === "-" ? C.textFaint : "#ffd9d2" }}>{d.char}</span>
            </span>
          ))}
          {info.differing.length > 14 && <span style={{ color: C.textFaint }}> …+{info.differing.length - 14} more</span>}
        </div>
      )}
    </div>
  );
}

function IconMini({ children, onClick, disabled }) {
  return (
    <button onClick={onClick} disabled={disabled} style={{
      all: "unset", cursor: disabled ? "default" : "pointer", padding: 5, borderRadius: 6,
      border: `1px solid ${C.border}`, color: disabled ? C.borderStrong : C.textDim, display: "flex",
    }}>{children}</button>
  );
}
