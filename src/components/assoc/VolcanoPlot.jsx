import { tooltipStyle } from "../../theme.js";
import { useMemo } from "react";
import {
  ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer,
  ReferenceLine, Cell,
} from "recharts";
import { C, FONT_DISPLAY } from "../../theme.js";


const X_CLAMP = 9; // log2 OR clamp so perfectly-segregating features stay on-plot

function xOf(row) {
  const or = row.or;
  if (!Number.isFinite(or)) return null;
  const v = Math.log2(Math.max(or, 1e-6));
  return Math.max(-X_CLAMP, Math.min(X_CLAMP, v));
}

/**
 * Volcano of association strength (log2 odds ratio) vs evidence (-log10 p).
 * `stat` selects which p-value column drives the Y axis: the raw Fisher test
 * or the clade-stratified CMH test.
 */
export default function VolcanoPlot({ rows, stat = "p", fdr, selected, onSelect }) {
  const points = useMemo(() => (
    rows
      .map((r) => {
        const p = stat === "cmhP" ? r.cmhP : r.p;
        const x = xOf(r);
        if (x === null || !Number.isFinite(p)) return null;
        return { ...r, x, y: -Math.log10(Math.max(p, 1e-300)), _p: p };
      })
      .filter(Boolean)
  ), [rows, stat]);

  const yMax = useMemo(() => {
    let m = 4;
    for (const pt of points) if (pt.y > m) m = pt.y;
    return Math.ceil(m) + 0.5;
  }, [points]);

  const yLine = -Math.log10(Math.max(fdr, 1e-300));

  const significant = (pt) => Number.isFinite(pt.q ?? NaN) ? pt.q <= fdr : false;

  return (
    <div style={{ width: "100%", height: 340 }}>
      <ResponsiveContainer>
        <ScatterChart margin={{ top: 12, right: 18, bottom: 22, left: 4 }}>
          <CartesianGrid stroke={C.border} strokeDasharray="3 3" />
          <XAxis
            type="number" dataKey="x" domain={[-X_CLAMP, X_CLAMP]}
            ticks={[-8, -6, -4, -2, 0, 2, 4, 6, 8]}
            tick={{ fill: C.textDim, fontSize: 11 }}
            label={{ value: "log2 odds ratio", position: "insideBottom", offset: -12, fill: C.textDim, fontSize: 11.5 }}
          />
          <YAxis
            type="number" dataKey="y" domain={[0, yMax]} allowDataOverflow
            tick={{ fill: C.textDim, fontSize: 11 }} width={52}
            label={{ value: "-log10 p", angle: -90, position: "insideLeft", fill: C.textDim, fontSize: 11.5 }}
          />
          <ReferenceLine y={yLine} stroke={C.bad} strokeDasharray="5 4"
            label={{ value: `FDR ${fdr}`, fill: C.bad, fontSize: 10, position: "right" }} />
          <ReferenceLine x={0} stroke={C.borderStrong} />
          <Tooltip
            content={<VolcanoTip fdr={fdr} />}
            contentStyle={tooltipStyle}
            cursor={{ strokeDasharray: "3 3", stroke: C.borderStrong }}
          />
          <Scatter data={points.filter((pt) => !significant(pt))} onClick={(pt) => onSelect?.(pt.gene)}>
            {points.filter((pt) => !significant(pt)).map((pt, i) => (
              <Cell key={i} fill={selected === pt.gene ? "#ffffff" : `${C.textFaint}`} fillOpacity={selected === pt.gene ? 1 : 0.55} />
            ))}
          </Scatter>
          <Scatter data={points.filter((pt) => significant(pt))} onClick={(pt) => onSelect?.(pt.gene)}>
            {points.filter((pt) => significant(pt)).map((pt, i) => (
              <Cell key={i} fill={pt.x >= 0 ? C.pheno : C.raw} fillOpacity={selected === pt.gene ? 1 : 0.85} />
            ))}
          </Scatter>
        </ScatterChart>
      </ResponsiveContainer>
    </div>
  );
}

function VolcanoTip({ active, payload, fdr }) {
  if (!active || !payload?.length) return null;
  const r = payload[0].payload;
  const fmtP = (v) => (Number.isFinite(v) ? v.toExponential(2) : "-");
  const qTxt = Number.isFinite(r.q) ? r.q.toExponential(2) : "-";
  return (
    <div style={{ ...tooltipStyle, padding: "8px 11px", maxWidth: 280 }}>
      <div style={{ fontFamily: FONT_DISPLAY, color: C.text, fontSize: 12.5 }}>{r.gene}</div>
      {r.product && <div style={{ color: C.textDim, fontSize: 11, marginTop: 2 }}>{r.product}</div>}
      <div style={{ color: C.textDim, fontSize: 11, marginTop: 5, lineHeight: 1.6 }}>
        {r.a}/{r.a + r.b} case vs {r.c}/{r.c + r.d} control present<br />
        OR {Number.isFinite(r.or) ? r.or.toFixed(2) : "∞"} [{Number.isFinite(r.lo) ? r.lo.toFixed(2) : "-"}, {Number.isFinite(r.hi) ? r.hi.toFixed(2) : "-"}]<br />
        p = {fmtP(r._p)} · q = {qTxt}
        {r.strataUsed > 0 && <> · CMH p = {fmtP(r.cmhP)}</>}
        {Number.isFinite(r.q) && r.q <= fdr && <span style={{ color: C.pheno, textShadow: `0 0 6px ${C.pheno}44` }}> · significant</span>}
      </div>
    </div>
  );
}
