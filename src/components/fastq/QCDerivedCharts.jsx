import { useMemo } from "react";
import {
  AreaChart, Area, LineChart, Line, BarChart, Bar, XAxis, YAxis, CartesianGrid,
  Tooltip, ResponsiveContainer, ReferenceDot, ReferenceLine,
} from "recharts";
import { C } from "../../theme.js";
import { Panel, Eyebrow} from "../ui/Primitives.jsx";
import { tooltipStyle } from "../../theme.js";
import { computeYieldByLength, computeGCHistogram, computePositionalQuality } from "../../lib/fastq/qcDerived.js";
import { formatBases } from "../../lib/fastq/stats.js";

function ChartPanel({ title, color, children, note }) {
  return (
    <Panel style={{ padding: 14 }}>
      <Eyebrow color={color}>{title}</Eyebrow>
      <div style={{ width: "100%", height: 200 }}><ResponsiveContainer>{children}</ResponsiveContainer></div>
      {note && <div style={{ fontSize: 10.5, color: C.textFaint, marginTop: 6 }}>{note}</div>}
    </Panel>
  );
}

export default function QCDerivedCharts({ dataset, index, previewRecords }) {
  const yieldData = useMemo(() => computeYieldByLength(index.lengths, dataset.n50), [index, dataset.n50]);
  const gcHist = useMemo(() => computeGCHistogram(index.gcPerRead), [index]);
  const posQuality = useMemo(() => computePositionalQuality(previewRecords), [previewRecords]);

  return (
    <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(300px, 1fr))", gap: 14 }}>
      <ChartPanel title="Yield by length (N50 marker)" color={C.qc}
        note={`Reading right-to-left: total bases contributed by reads at or above each length. The marked point is N50 (${dataset.n50.toLocaleString()} bp) — half of all sequenced bases come from reads at least this long.`}>
        <AreaChart data={yieldData.points} margin={{ left: 4, right: 12 }}>
          <CartesianGrid stroke={C.border} vertical={false} />
          <XAxis dataKey="length" type="number" reversed tick={{ fill: C.textFaint, fontSize: 10 }}
            label={{ value: "min. read length (bp)", fill: C.textFaint, fontSize: 10, position: "insideBottom", dy: 10 }} />
          <YAxis tickFormatter={(v) => formatBases(v)} tick={{ fill: C.textFaint, fontSize: 10 }} width={54} />
          <Tooltip contentStyle={tooltipStyle} formatter={(v) => formatBases(v)} labelFormatter={(l) => `≥ ${l.toLocaleString()} bp`} />
          <Area type="monotone" dataKey="cumulativeBases" stroke={C.qc} fill={`${C.qc}22`} strokeWidth={2} />
          <ReferenceLine x={dataset.n50} stroke={C.pheno} strokeDasharray="4 3" />
          <ReferenceDot x={dataset.n50} y={yieldData.n50CumBases} r={5} fill={C.pheno} stroke="#fff" strokeWidth={1.5} />
        </AreaChart>
      </ChartPanel>

      <ChartPanel title="Per-read GC content" color={C.qc}
        note="Distribution of GC% within individual reads — a second peak far from your organism's expected GC% can flag contamination or a mixed sample.">
        <BarChart data={gcHist}>
          <CartesianGrid stroke={C.border} vertical={false} />
          <XAxis dataKey="label" tick={{ fill: C.textFaint, fontSize: 10 }} label={{ value: "% GC", fill: C.textFaint, fontSize: 10, position: "insideBottom", dy: 10 }} />
          <YAxis tick={{ fill: C.textFaint, fontSize: 10 }} />
          <Tooltip contentStyle={tooltipStyle} formatter={(v, n, p) => [v, `${p.payload.rangeStart}-${p.payload.rangeEnd}%`]} />
          <Bar dataKey="count" fill={C.qc} radius={[2, 2, 0, 0]} />
        </BarChart>
      </ChartPanel>

      <Panel style={{ padding: 14, gridColumn: "1 / -1" }}>
        <Eyebrow color={C.qc}>Mean quality by position in read</Eyebrow>
        <div style={{ width: "100%", height: 200 }}>
          <ResponsiveContainer>
            <LineChart data={posQuality} margin={{ left: 4, right: 12 }}>
              <CartesianGrid stroke={C.border} vertical={false} />
              <XAxis dataKey="position" tick={{ fill: C.textFaint, fontSize: 10 }} label={{ value: "position in read (bp)", fill: C.textFaint, fontSize: 10, position: "insideBottom", dy: 10 }} />
              <YAxis tick={{ fill: C.textFaint, fontSize: 10 }} domain={["dataMin - 2", "dataMax + 2"]} />
              <Tooltip contentStyle={tooltipStyle} formatter={(v, n, p) => [`Q${v.toFixed(1)} (n=${p.payload.readsAtPosition})`, "mean quality"]} />
              <Line type="monotone" dataKey="meanQ" stroke={C.raw} dot={false} strokeWidth={2} />
            </LineChart>
          </ResponsiveContainer>
        </div>
        <div style={{ fontSize: 10.5, color: C.textFaint, marginTop: 6 }}>
          Computed from a sample of up to {previewRecords?.length.toLocaleString() ?? 0} reads with full quality text in memory (see FASTQ page for why very large gzip files only keep a bounded sample) — not the entire dataset. A declining tail is normal and typical for nanopore reads.
        </div>
      </Panel>
    </div>
  );
}
