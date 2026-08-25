import React from "react";
import { CheckCircle2, AlertTriangle, XCircle } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow } from "../ui/Primitives.jsx";

/**
 * A handful of common, well-known long-read QC heuristics (not a validated
 * clinical/production spec - thresholds are conventional rules of thumb
 * used when eyeballing a Nanopore run, and are labeled as such).
 */
function evaluateChecks(dataset) {
  const checks = [];

  checks.push({
    label: "N50",
    value: `${dataset.n50.toLocaleString()} bp`,
    level: dataset.n50 >= 3000 ? "good" : dataset.n50 >= 1000 ? "warn" : "bad",
    note: dataset.n50 >= 3000
      ? "Healthy for most long-read assembly workflows."
      : dataset.n50 >= 1000
      ? "Usable, but longer reads would help assembly contiguity."
      : "Quite short for long-read assembly - check size selection / library prep.",
  });

  checks.push({
    label: "Mean quality",
    value: `Q${dataset.datasetMeanQ.toFixed(1)}`,
    level: dataset.datasetMeanQ >= 15 ? "good" : dataset.datasetMeanQ >= 10 ? "warn" : "bad",
    note: dataset.datasetMeanQ >= 15
      ? "Consistent with a modern high-accuracy basecalling model."
      : dataset.datasetMeanQ >= 10
      ? "Typical of a fast/older basecalling model - consider re-basecalling with a high-accuracy model if available."
      : "Low even for raw nanopore reads - check flow cell health and basecaller settings.",
  });

  checks.push({
    label: "Total yield",
    value: dataset.totalBases >= 1e9 ? `${(dataset.totalBases / 1e9).toFixed(2)} Gb` : `${(dataset.totalBases / 1e6).toFixed(1)} Mb`,
    level: dataset.totalBases >= 3e8 ? "good" : dataset.totalBases >= 5e7 ? "warn" : "bad",
    note: "Whether this is \"enough\" depends entirely on your target genome size and desired coverage - use the yield-by-length plot below to estimate coverage at a given length cutoff.",
  });

  const gcDev = Math.abs(dataset.gcPercent - 50);
  checks.push({
    label: "GC content",
    value: `${dataset.gcPercent.toFixed(1)}%`,
    level: gcDev < 15 ? "good" : gcDev < 25 ? "warn" : "bad",
    note: "Compare against your organism's known genomic GC% - a large mismatch can indicate contamination or a mislabeled sample, not just \"good\" or \"bad\" in isolation.",
  });

  if (dataset.baseCounts.N > 0) {
    const nPct = (dataset.baseCounts.N / dataset.totalBases) * 100;
    checks.push({
      label: "Ambiguous bases (N)",
      value: `${nPct.toFixed(2)}%`,
      level: nPct < 0.5 ? "good" : nPct < 2 ? "warn" : "bad",
      note: "Higher N content than usual for basecalled Nanopore data - worth double-checking the basecaller output.",
    });
  }

  return checks;
}

const LEVEL_STYLE = {
  good: { color: C.good, Icon: CheckCircle2 },
  warn: { color: C.qc, Icon: AlertTriangle },
  bad: { color: C.bad, Icon: XCircle },
};

export default function QCVerdictPanel({ dataset }) {
  const checks = evaluateChecks(dataset);
  return (
    <Panel style={{ padding: 18 }}>
      <Eyebrow color={C.qc}>Run quality — heuristic guidance</Eyebrow>
      <div style={{ fontSize: 11.5, color: C.textFaint, marginBottom: 12 }}>
        Conventional rules of thumb for eyeballing a long-read run, not a validated pass/fail spec. What counts as "good enough" always depends on your specific experiment.
      </div>
      <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(220px, 1fr))", gap: 10 }}>
        {checks.map((c) => {
          const { color, Icon } = LEVEL_STYLE[c.level];
          return (
            <div key={c.label} style={{ display: "flex", gap: 10, padding: "10px 12px", borderRadius: 2, border: `1px solid ${color}44`, background: "#05070a" }}>
              <Icon size={16} color={color} style={{ flexShrink: 0, marginTop: 1 }} />
              <div>
                <div style={{ fontSize: 12.5, color: C.text, display: "flex", gap: 6, alignItems: "baseline", fontFamily: FONT_DISPLAY }}>
                  <span>{c.label}</span><span style={{ color, fontWeight: 600 }}>{c.value}</span>
                </div>
                <div style={{ fontSize: 11, color: C.textDim, marginTop: 2, lineHeight: 1.4 }}>{c.note}</div>
              </div>
            </div>
          );
        })}
      </div>
    </Panel>
  );
}
