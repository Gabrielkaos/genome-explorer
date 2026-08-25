import { useMemo, useState } from "react";
import { BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine, Cell } from "recharts";
import { ChevronDown, ChevronUp } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";

const fmtSci = (v) => (Number.isFinite(v) ? v.toExponential(2) : "—");
const fmtNum = (v, d = 2) => (Number.isFinite(v) ? v.toFixed(d) : "—");

/**
 * Sortable, searchable results table. For binary traits every statistic of
 * the 2x2 test is exposed; for continuous traits the Welch t-test columns.
 * `stat` picks which p-value column drives sorting default + significance
 * coloring (raw Fisher/Welch vs clade-stratified CMH).
 */
export default function ResultsTable({ rows, traitType, stat = "p", fdr, selected, onSelect }) {
  const [sortKey, setSortKey] = useState(stat === "cmhP" ? "cmhP" : "p");
  const [sortDir, setSortDir] = useState(1);
  const [query, setQuery] = useState("");
  const [limit, setLimit] = useState(150);

  const view = useMemo(() => {
    const q = query.trim().toLowerCase();
    let r = q ? rows.filter((row) => row.gene.toLowerCase().includes(q) || (row.product ?? "").toLowerCase().includes(q)) : rows.slice();
    const keyOf = (row) => {
      if (sortKey === "or") return Number.isFinite(row.or) ? Math.abs(Math.log2(row.or)) : -Infinity;
      return Number.isFinite(row[sortKey]) ? row[sortKey] : Infinity;
    };
    r.sort((a, b) => (keyOf(a) - keyOf(b)) * sortDir);
    return r;
  }, [rows, query, sortKey, sortDir]);

  const shown = view.slice(0, limit);

  function header(label, key, title) {
    const active = sortKey === key;
    return (
      <th onClick={() => { if (active) setSortDir((d) => -d); else { setSortKey(key); setSortDir(1); } }}
        title={title}
        style={{
          padding: "6px 9px", cursor: "pointer", userSelect: "none", textAlign: "right",
          color: active ? C.pheno : C.textFaint, borderBottom: `1px solid ${active ? C.pheno : C.border}`,
          whiteSpace: "nowrap",
        }}>
        {label} {active && (sortDir > 0 ? <ChevronUp size={10} style={{ verticalAlign: "-1px" }} /> : <ChevronDown size={10} style={{ verticalAlign: "-1px" }} />)}
      </th>
    );
  }

  const sigOf = (r) => {
    const p = stat === "cmhP" ? r.cmhP : r.p;
    const q = stat === "cmhP" ? null : r.q;
    return Number.isFinite(p) && (q != null ? q <= fdr : p <= 0.01);
  };

  return (
    <div>
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", gap: 10, marginBottom: 8, flexWrap: "wrap" }}>
        <input value={query} onChange={(e) => setQuery(e.target.value)} placeholder="Filter by gene or product…"
          spellCheck={false}
          style={{
            flex: "0 1 280px", background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 2,
            padding: "6px 10px", color: C.text, fontFamily: FONT_DISPLAY, fontSize: 12,
          }} />
        <span style={{ fontSize: 11, color: C.textFaint }}>
          {view.length.toLocaleString()} feature(s){query && ` match "${query}"`}{view.length > limit && ` — showing ${limit}`}
        </span>
      </div>

      <div style={{ overflowX: "auto", maxHeight: 460, overflowY: "auto", border: `1px solid ${C.border}`, borderRadius: 2 }}>
        <table style={{ borderCollapse: "collapse", fontFamily: FONT_DISPLAY, fontSize: 11.5, width: "100%" }}>
          <thead style={{ position: "sticky", top: 0, background: "#0c1018", zIndex: 1 }}>
            <tr>
              <th onClick={() => { if (sortKey === "gene" && sortDir === -1) setSortDir(1); else { setSortKey("gene"); setSortDir(-1); } }}
                style={{ padding: "6px 9px", cursor: "pointer", textAlign: "left", color: sortKey === "gene" ? C.pheno : C.textFaint, borderBottom: `1px solid ${C.border}` }}>
                gene
              </th>
              {traitType === "binary" ? (
                <>
                  {header("case+", "a", "present among cases")}
                  {header("ctrl+", "b", "present among controls")}
                  {header("OR", "or", "odds ratio (Haldane-corrected)")}
                  {header("95% CI", "lo", "Woolf log confidence interval")}
                  {header("p Fisher", "p")}
                  {header("q FDR", "q", "Benjamini-Hochberg adjusted")}
                  {header("CMH p", "cmhP", "clade-stratified Mantel-Haenszel")}
                </>
              ) : (
                <>
                  {header("mean present", "meanPresent")}
                  {header("mean absent", "meanAbsent")}
                  {header("Δ mean", "meanDiff")}
                  {header("t", "t")}
                  {header("p Welch", "p")}
                  {header("q FDR", "q")}
                  {header("CMH p", "cmhP")}
                </>
              )}
            </tr>
          </thead>
          <tbody>
            {shown.map((r) => {
              const isSel = selected === r.gene;
              const sig = sigOf(r);
              return (
                <tr key={r.gene} onClick={() => onSelect?.(isSel ? null : r.gene)}
                  style={{ cursor: "pointer", borderBottom: `1px solid ${C.border}`,
                    background: isSel ? `${C.pheno}22` : "transparent" }}>
                  <td style={{ padding: "5px 9px", color: C.text }} title={r.product ?? ""}>
                    {r.gene}
                    {sig && <span style={{ color: C.pheno, marginLeft: 5 }}>●</span>}
                  </td>
                  {traitType === "binary" ? (
                    <>
                      <NumCell>{r.a}</NumCell>
                      <NumCell>{r.b}</NumCell>
                      <NumCell>{fmtNum(r.or)}</NumCell>
                      <NumCell>[{fmtNum(r.lo)}, {fmtNum(r.hi)}]</NumCell>
                      <NumCell>{fmtSci(r.p)}</NumCell>
                      <NumCell>{fmtSci(r.q)}</NumCell>
                      <NumCell>{Number.isFinite(r.cmhP) ? fmtSci(r.cmhP) : "—"}</NumCell>
                    </>
                  ) : (
                    <>
                      <NumCell>{fmtNum(r.meanPresent)}</NumCell>
                      <NumCell>{fmtNum(r.meanAbsent)}</NumCell>
                      <NumCell>{fmtNum(r.meanDiff)}</NumCell>
                      <NumCell>{fmtNum(r.t)}</NumCell>
                      <NumCell>{fmtSci(r.p)}</NumCell>
                      <NumCell>{fmtSci(r.q)}</NumCell>
                      <NumCell>{Number.isFinite(r.cmhP) ? fmtSci(r.cmhP) : "—"}</NumCell>
                    </>
                  )}
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      {view.length > limit && (
        <button onClick={() => setLimit((l) => l + 500)}
          style={{
            all: "unset", cursor: "pointer", display: "block", margin: "10px auto 0",
            fontSize: 12, color: C.raw, border: `1px solid ${C.border}`, borderRadius: 2, padding: "6px 14px",
            fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
            textShadow: `0 0 6px ${C.raw}44`,
          }}>
          [ show_more ({(view.length - limit).toLocaleString()}_hidden) ]
        </button>
      )}
    </div>
  );
}

function NumCell({ children }) {
  return <td style={{ padding: "5px 9px", textAlign: "right", color: C.textDim, whiteSpace: "nowrap" }}>{children}</td>;
}

/**
 * Manhattan-style ranked bar chart of -log10 p for continuous traits
 * (no odds ratio exists there, so a volcano is meaningless).
 */
export function RankedBar({ rows, stat = "p", fdr }) {
  const points = useMemo(() => (
    rows
      .map((r) => {
        const p = stat === "cmhP" ? r.cmhP : r.p;
        if (!Number.isFinite(p)) return null;
        return { gene: r.gene, y: -Math.log10(Math.max(p, 1e-300)), q: r.q };
      })
      .filter(Boolean)
      .slice(0, 400)
  ), [rows, stat]);
  if (!points.length) return null;
  const yLine = -Math.log10(Math.max(fdr, 1e-300));
  return (
    <div style={{ width: "100%", height: 260 }}>
      <ResponsiveContainer>
        <BarChart data={points} margin={{ top: 12, right: 18, bottom: 40, left: 4 }}>
          <XAxis dataKey="gene" tick={false} label={{ value: "features ranked by evidence (first 400)", position: "insideBottom", offset: -16, fill: C.textDim, fontSize: 11.5 }} />
          <YAxis tick={{ fill: C.textDim, fontSize: 11 }} width={46}
            label={{ value: "-log10 p", angle: -90, position: "insideLeft", fill: C.textDim, fontSize: 11.5 }} />
          <Tooltip contentStyle={{ background: "#0c1018", border: `1px solid ${C.border}`, borderRadius: 2, fontSize: 12, color: C.text, fontFamily: FONT_DISPLAY }}
            formatter={(v) => [`-log10 p = ${v.toFixed(2)}`]} />
          <ReferenceLine y={yLine} stroke={C.bad} strokeDasharray="5 4" />
          <Bar dataKey="y">
            {points.map((pt, i) => (
              <Cell key={i} fill={pt.y >= yLine ? C.pheno : C.textFaint} fillOpacity={0.85} />
            ))}
          </Bar>
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
}
