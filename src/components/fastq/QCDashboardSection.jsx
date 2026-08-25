import { Microscope, Link2 } from "lucide-react";
import { C } from "../../theme.js";
import { SectionTitle, ExplainBox, Panel } from "../ui/Primitives.jsx";
import { useFastqData } from "../../state/FastqDataContext.jsx";
import FastqUploader from "./FastqUploader.jsx";
import FastqSummaryStats from "./FastqSummaryStats.jsx";
import FastqCharts from "./FastqCharts.jsx";
import QCDerivedCharts from "./QCDerivedCharts.jsx";
import QCVerdictPanel from "./QCVerdictPanel.jsx";

export default function QCDashboardSection({ explainMode }) {
  const { status, progress, dataset, index, previewRecords, error, parseFile, cancel } = useFastqData();

  return (
    <div>
      <SectionTitle icon={Microscope} color={C.qc} title="QC Dashboard"
        subtitle="Read-only diagnostics — measurement only. Nothing here filters or modifies your reads." />

      <ExplainBox explainMode={explainMode} color={C.qc}>
        This mirrors what a tool like NanoPlot does: it describes the dataset you already have. If you want to actually remove low-quality or short reads, use the filter controls on the FASTQ page — that's a deliberately separate step.
      </ExplainBox>

      {dataset && (
        <div style={{ display: "flex", alignItems: "center", gap: 6, fontSize: 11.5, color: C.textFaint, margin: "14px 0" }}>
          <Link2 size={12} /> Showing the same dataset loaded on the FASTQ page ({dataset.fileName}). Load a different file below to replace it everywhere.
        </div>
      )}

      <div style={{ marginBottom: 20 }}>
        <FastqUploader status={status} progress={progress} onFile={parseFile} onCancel={cancel} error={error} />
      </div>

      {status === "done" && dataset && (
        <>
          <div style={{ marginBottom: 20 }}>
            <FastqSummaryStats dataset={dataset} passCount={dataset.totalReads} passPct={100} />
          </div>

          <div style={{ marginBottom: 20 }}>
            <QCVerdictPanel dataset={dataset} />
          </div>

          <div style={{ marginBottom: 20 }}>
            <FastqCharts dataset={dataset} index={index} neutral />
          </div>

          <QCDerivedCharts dataset={dataset} index={index} previewRecords={previewRecords} />
        </>
      )}

      {status !== "done" && !dataset && (
        <Panel style={{ padding: 16, fontSize: 12.5, color: C.textFaint }}>
          No data loaded yet. Upload a file above (or on the FASTQ page — they share the same dataset).
        </Panel>
      )}
    </div>
  );
}
