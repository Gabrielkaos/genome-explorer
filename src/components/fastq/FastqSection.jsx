import React, { useState } from "react";
import { ScrollText } from "lucide-react";
import { C } from "../../theme.js";
import { SectionTitle } from "../ui/Primitives.jsx";
import { useFastqData } from "../../state/FastqDataContext.jsx";
import FastqUploader from "./FastqUploader.jsx";
import FastqSummaryStats from "./FastqSummaryStats.jsx";
import FastqFilterControls from "./FastqFilterControls.jsx";
import FastqCharts from "./FastqCharts.jsx";
import FastqReadTable from "./FastqReadTable.jsx";

export default function FastqSection({ explainMode }) {
  const { status, progress, dataset, index, previewRecords, error, parseFile, cancel, getRecord, file } = useFastqData();

  const [qualityThreshold, setQualityThreshold] = useState(9);
  const [minLength, setMinLength] = useState(0);
  const [maxLength, setMaxLength] = useState(null);

  function handleFile(f) {
    parseFile(f);
  }

  let passCount = 0;
  if (index) {
    for (let i = 0; i < index.lengths.length; i++) {
      if (index.meanQs[i] >= qualityThreshold && index.lengths[i] >= minLength && (!maxLength || index.lengths[i] <= maxLength)) passCount++;
    }
  }
  const passPct = dataset?.totalReads ? (passCount / dataset.totalReads) * 100 : 0;

  return (
    <div>
      <SectionTitle icon={ScrollText} color={C.raw} title="FASTQ" subtitle="Upload and inspect real long-read sequencing data — parsed entirely client-side" />

      <div style={{ marginBottom: 18 }}>
        <FastqUploader status={status} progress={progress} onFile={handleFile} onCancel={cancel} error={error} />
      </div>

      {status === "done" && dataset && (
        <>
          <div style={{ marginBottom: 20 }}>
            <FastqSummaryStats dataset={dataset} passCount={passCount} passPct={passPct} />
          </div>

          <div style={{ marginBottom: 20 }}>
            <FastqFilterControls
              dataset={dataset} index={index} file={file} previewRecords={previewRecords}
              qualityThreshold={qualityThreshold} setQualityThreshold={setQualityThreshold}
              minLength={minLength} setMinLength={setMinLength}
              maxLength={maxLength} setMaxLength={setMaxLength}
              passCount={passCount} passPct={passPct} explainMode={explainMode}
            />
          </div>

          <div style={{ marginBottom: 20 }}>
            <FastqCharts dataset={dataset} index={index} qualityThreshold={qualityThreshold} minLength={minLength} maxLength={maxLength} />
          </div>

          <FastqReadTable
            dataset={dataset} index={index}
            qualityThreshold={qualityThreshold} minLength={minLength} maxLength={maxLength}
            getRecord={getRecord}
          />
        </>
      )}

      {status === "idle" && !dataset && (
        <div style={{ fontSize: 12.5, color: C.textFaint, marginTop: 4 }}>
          No file loaded yet. Upload a real .fastq / .fastq.gz run, or load the synthetic sample to try the tool first.
        </div>
      )}
    </div>
  );
}
