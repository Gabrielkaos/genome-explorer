import { C, FONT_DISPLAY } from "../../theme.js";
import { Eyebrow, Panel } from "../ui/Primitives.jsx";

/** Short label for tight matrix edges; hover any cell/label for the full id. */
const shortCol = (id) => (id.length > 9 ? `${id.slice(0, 8)}…` : id);
const shortRow = (id) => (id.length > 13 ? `${id.slice(0, 12)}…` : id);

/**
 * Lower-triangle heatmap of model-corrected pairwise distances
 * (substitutions/site). Hover a cell for the raw transition/transversion
 * counts behind that estimate.
 */
export default function DistanceMatrix({ ids, matrix, transitions, transversions }) {
  const n = ids.length;
  let min = Infinity, max = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const v = matrix[i * n + j];
      if (v < min) min = v;
      if (v > max) max = v;
    }
  }
  const span = Math.max(1e-9, max - min);

  return (
    <Panel>
      <Eyebrow color={C.phylo}>Pairwise distances (substitutions/site)</Eyebrow>
      <div style={{ overflow: "auto", maxHeight: 340, maxWidth: "100%", marginTop: 8 }}>
        <table style={{ borderCollapse: "collapse", fontFamily: FONT_DISPLAY, fontSize: 10.5 }}>
          <thead>
            <tr>
              <th style={corner} />
              {ids.map((id) => <th key={id} style={headCell} title={id}>{shortCol(id)}</th>)}
            </tr>
          </thead>
          <tbody>
            {ids.map((rowId, i) => (
              <tr key={rowId}>
                <td style={rowHead} title={rowId}>{shortRow(rowId)}</td>
                {ids.map((colId, j) => {
                  if (j > i) return <td key={colId} style={{ ...cell, background: "rgba(255,255,255,0.02)" }} />;
                  if (i === j) return <td key={colId} style={{ ...cell, background: "rgba(255,255,255,0.07)", color: C.textFaint }}>—</td>;
                  const v = matrix[i * n + j];
                  const t = (v - min) / span;
                  const ti = transitions?.[i * n + j];
                  const tv = transversions?.[i * n + j];
                  return (
                    <td key={colId}
                      title={`${rowId} vs ${colId}: d = ${v.toFixed(4)}${Number.isFinite(ti) ? ` (${ti} ti, ${tv} tv)` : ""}`}
                      style={{
                        ...cell,
                        background: heat(t),
                        color: t > 0.55 ? "#0a0f0c" : "#cfe3d8",
                        fontWeight: 500,
                      }}>
                      {v >= 10 ? "—" : v.toFixed(3)}
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <div style={{ display: "flex", alignItems: "center", gap: 8, flexWrap: "wrap", marginTop: 8, fontSize: 10.5, color: C.textFaint }}>
        <span>close</span>
        {[0, 0.25, 0.5, 0.75, 1].map((t) => (
          <span key={t} style={{ width: 26, height: 10, borderRadius: 2, background: heat(t), display: "inline-block" }} />
        ))}
        <span>divergent</span>
        <span style={{ marginLeft: "auto" }}>range {min.toFixed(3)} – {max.toFixed(3)}</span>
      </div>
    </Panel>
  );
}

function heat(t) {
  // green (close) -> amber -> red (divergent)
  const r = Math.round(104 + (230 - 104) * t);
  const g = Math.round(201 - (201 - 90) * t);
  const b = Math.round(143 - (143 - 70) * t);
  return `rgba(${r},${g},${b},${0.35 + t * 0.45})`;
}

const clip = { overflow: "hidden", whiteSpace: "nowrap", textOverflow: "ellipsis" };
const corner = { padding: "4px 6px", position: "sticky", top: 0, left: 0, zIndex: 2, background: "#10151d" };
const headCell = {
  padding: "4px 6px", color: C.textFaint, textAlign: "right", fontWeight: 400,
  position: "sticky", top: 0, zIndex: 1, background: "#10151d",
  minWidth: 46, maxWidth: 58, ...clip,
};
const rowHead = {
  padding: "4px 6px", color: C.textDim, textAlign: "left",
  position: "sticky", left: 0, zIndex: 1, background: "#10151d",
  minWidth: 92, maxWidth: 92, ...clip,
};
const cell = { padding: "3px 5px", textAlign: "right", minWidth: 46 };
