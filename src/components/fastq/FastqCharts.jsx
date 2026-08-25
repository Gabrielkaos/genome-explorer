import { useMemo } from "react";
import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer,
  ScatterChart, Scatter, Cell, ZAxis,
} from "recharts";
import { C } from "../../theme.js";
import { Panel, Eyebrow} from "../ui/Primitives.jsx";
import { tooltipStyle } from "../../theme.js";

function ChartPanel({ title, color, children }) {
  return (
    <Panel style={{ padding: 14 }}>
      <Eyebrow color={color}>{title}</Eyebrow>
      <div style={{ width: "100%", height: 190 }}><ResponsiveContainer>{children}</ResponsiveContainer></div>
    </Panel>
  );
}

export default function FastqCharts({ dataset, index, qualityThreshold, minLength, maxLength, neutral = false }) {
  // Scatter is sampled (not the full dataset) to keep the SVG light even
  // for million-read runs; histograms below use the exact full-dataset bins.
  const scatterSample = useMemo(() => {
    if (!index) return [];
    const n = index.lengths.length;
    const step = Math.max(1, Math.floor(n / 3000));
    const pts = [];
    for (let i = 0; i < n; i += step) {
      const len = index.lengths[i], q = index.meanQs[i];
      const pass = neutral ? true : (q >= qualityThreshold && len >= minLength && (maxLength ? len <= maxLength : true));
      pts.push({ x: len, y: q, pass });
    }
    return pts;
  }, [index, qualityThreshold, minLength, maxLength, neutral]);

  if (!dataset) return null;

  return (
    <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(300px, 1fr))", gap: 14 }}>
      <ChartPanel title="Read length distribution" color={C.qc}>
        <BarChart data={dataset.lengthHistogram}>
          <CartesianGrid stroke={C.border} vertical={false} />
          <XAxis dataKey="label" tick={{ fill: C.textFaint, fontSize: 10 }} />
          <YAxis tick={{ fill: C.textFaint, fontSize: 10 }} />
          <Tooltip contentStyle={tooltipStyle} formatter={(v, n, p) => [v, `${p.payload.rangeStart}-${p.payload.rangeEnd} bp`]} />
          <Bar dataKey="count" fill={C.qc} radius={[2, 2, 0, 0]} />
        </BarChart>
      </ChartPanel>

      <ChartPanel title="Quality distribution" color={C.qc}>
        <BarChart data={dataset.qualityHistogram}>
          <CartesianGrid stroke={C.border} vertical={false} />
          <XAxis dataKey="label" tick={{ fill: C.textFaint, fontSize: 10 }} label={{ value: "mean Q", fill: C.textFaint, fontSize: 10, position: "insideBottom", dy: 10 }} />
          <YAxis tick={{ fill: C.textFaint, fontSize: 10 }} />
          <Tooltip contentStyle={tooltipStyle} />
          <Bar dataKey="count" fill={C.raw} radius={[2, 2, 0, 0]} />
        </BarChart>
      </ChartPanel>

      <ChartPanel title={`Quality vs. length (${scatterSample.length.toLocaleString()}-read sample)`} color={C.qc}>
        <ScatterChart>
          <CartesianGrid stroke={C.border} />
          <XAxis dataKey="x" name="length" tick={{ fill: C.textFaint, fontSize: 10 }} />
          <YAxis dataKey="y" name="quality" tick={{ fill: C.textFaint, fontSize: 10 }} />
          <ZAxis range={[16, 16]} />
          <Tooltip contentStyle={tooltipStyle} cursor={{ strokeDasharray: "3 3" }} />
          <Scatter data={scatterSample}>
            {scatterSample.map((p, i) => <Cell key={i} fill={neutral ? C.qc : (p.pass ? C.good : C.bad)} fillOpacity={0.65} />)}
          </Scatter>
        </ScatterChart>
      </ChartPanel>

      <Panel style={{ padding: 14 }}>
        <Eyebrow color={C.qc}>Base composition</Eyebrow>
        <div style={{ display: "flex", height: 190, alignItems: "flex-end", gap: 10, padding: "0 10px" }}>
          {["A", "C", "G", "T", "N"].map((b) => {
            const count = dataset.baseCounts[b] || 0;
            const total = dataset.totalBases || 1;
            const pct = (count / total) * 100;
            const colors = { A: C.raw, C: C.qc, G: C.assembly, T: C.annotation, N: C.textFaint };
            return (
              <div key={b} style={{ flex: 1, display: "flex", flexDirection: "column", alignItems: "center", justifyContent: "flex-end", height: "100%" }}>
                <div style={{ fontSize: 10.5, color: C.textDim, marginBottom: 4 }}>{pct.toFixed(1)}%</div>
                <div style={{ width: "60%", height: `${Math.max(2, pct * 1.6)}px`, background: colors[b], borderRadius: "2px 2px 0 0" }} />
                <div style={{ fontSize: 12, color: C.text, marginTop: 6 }}>{b}</div>
              </div>
            );
          })}
        </div>
      </Panel>
    </div>
  );
}
