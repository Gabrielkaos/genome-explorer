import { C, FONT_DISPLAY } from "../../theme.js";
import { Eyebrow, Panel } from "../ui/Primitives.jsx";

const MAX_ROWS = 300;

/** Variable-site table; clicking a row jumps the viewer to that column. */
export default function VariantTable({ sites, totalSites, selectedCol, onSelect }) {
  const shown = sites.slice(0, MAX_ROWS);
  return (
    <Panel>
      <Eyebrow color={C.phylo}>
        Variable sites {totalSites > shown.length ? `(showing ${shown.length} of ${totalSites})` : `(${sites.length})`}
      </Eyebrow>
      <div style={{ overflow: "auto", maxHeight: 340, marginTop: 8 }}>
        <table style={{ borderCollapse: "collapse", fontFamily: FONT_DISPLAY, fontSize: 11.5, width: "100%" }}>
          <thead>
            <tr>
              {["Col", "Consensus", "Alleles", "Type", "Informative"].map((h) => (
                <th key={h} style={head}>{h}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {shown.map((s) => (
              <tr key={s.pos} onClick={() => onSelect(s.pos)}
                style={{
                  cursor: "pointer",
                  background: s.pos === selectedCol ? "rgba(255,255,255,0.08)" : "transparent",
                }}>
                <td style={cellTd}>{s.pos + 1}</td>
                <td style={{ ...cellTd, color: BASE_COLOR_MAP[s.consensus] || C.textDim }}>{s.consensus}</td>
                <td style={{ ...cellTd, color: C.text }}>{s.alleles.join("  ")}</td>
                <td style={{ ...cellTd, color: s.indel ? C.qc : C.textDim }}>{s.indel ? "indel" : "SNP"}</td>
                <td style={{ ...cellTd, color: s.informative ? C.good : C.textFaint }}>
                  {s.informative ? "yes" : "singleton"}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
        {totalSites > shown.length && (
          <div style={{ fontSize: 10.5, color: C.textFaint, padding: "8px 4px" }}>
            …{totalSites - shown.length} more — export the variants TSV for the full list.
          </div>
        )}
      </div>
    </Panel>
  );
}

const BASE_COLOR_MAP = { A: "#68c98f", C: "#5ec8d8", G: "#e8c15a", T: "#ef7fa3" };

const head = { textAlign: "left", padding: "5px 10px", color: C.textFaint, position: "sticky", top: 0, background: "#10151d", letterSpacing: "0.04em" };
const cellTd = { padding: "4px 10px", borderBottom: `1px solid ${C.border}` };
