import { C, FONT_DISPLAY } from "../../theme.js";
import { Eyebrow, Panel } from "../ui/Primitives.jsx";

/** Short label for tight matrix edges; hover any cell/label for the full id. */
const shortCol = (id) => (id.length > 9 ? `${id.slice(0, 8)}…` : id);
const shortRow = (id) => (id.length > 13 ? `${id.slice(0, 12)}…` : id);

/** Lower-triangle pairwise identity matrix with heat coloring. */
export default function IdentityMatrix({ ids, matrix }) {
  const n = ids.length;
  return (
    <Panel>
      <Eyebrow color={C.phylo}>Pairwise identity (%)</Eyebrow>
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
                  if (j > i) return <td key={colId} style={{ ...cell, background: "transparent" }} />;
                  const id = matrix[i * n + j];
                  const isDiag = i === j;
                  return (
                    <td key={colId} title={`${rowId} vs ${colId}: ${(id * 100).toFixed(2)}%`}
                      style={{
                        ...cell,
                        background: isDiag ? "#0e1320" : heat(id),
                        color: isDiag ? C.textFaint : id > 0.55 ? "#0a0f0c" : C.text,
                        fontWeight: !isDiag && id < 0.9 ? 600 : 400,
                      }}>
                      {isDiag ? "—" : (id * 100).toFixed(1)}
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </Panel>
  );
}

function heat(id) {
  if (id >= 0.5) {
    const t = (id - 0.5) / 0.5;
    return `rgba(104,201,143,${0.15 + t * 0.75})`;
  }
  const t = Math.max(0, id) / 0.5;
  return `rgba(230,104,95,${0.35 - t * 0.25})`;
}

const clip = { overflow: "hidden", whiteSpace: "nowrap", textOverflow: "ellipsis" };
const corner = { padding: "4px 6px", position: "sticky", top: 0, left: 0, zIndex: 2, background: C.bgPanel };
const headCell = {
  padding: "4px 6px", color: C.textFaint, textAlign: "right", fontWeight: 400,
  position: "sticky", top: 0, zIndex: 1, background: C.bgPanel,
  minWidth: 46, maxWidth: 58, ...clip,
};
const rowHead = {
  padding: "4px 6px", color: C.textDim, textAlign: "left",
  position: "sticky", left: 0, zIndex: 1, background: C.bgPanel,
  minWidth: 92, maxWidth: 92, ...clip,
};
const cell = { padding: "3px 5px", textAlign: "right", minWidth: 46 };
