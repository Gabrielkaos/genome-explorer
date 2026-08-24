import { useCallback, useRef, useState } from "react";

/**
 * Manages the lifecycle of parsing one FASTQ file: spins up the worker,
 * tracks progress, stores the resulting lightweight per-read index (NOT
 * full sequence text for every read - see fastqParser.worker.js for why),
 * and exposes a lazy `getRecord(index)` for on-demand full-text lookup.
 */
export function useFastqParser() {
  const [status, setStatus] = useState("idle"); // idle | parsing | done | error
  const [progress, setProgress] = useState({ bytesProcessed: 0, totalBytes: 0, readsProcessed: 0 });
  const [dataset, setDataset] = useState(null); // aggregate stats object
  const [index, setIndex] = useState(null);     // { lengths, meanQs, ids, offsets }
  const [previewRecords, setPreviewRecords] = useState(null);
  const [error, setError] = useState(null);

  const workerRef = useRef(null);
  const fileRef = useRef(null);
  const recordCacheRef = useRef(new Map());

  const reset = useCallback(() => {
    setStatus("idle");
    setProgress({ bytesProcessed: 0, totalBytes: 0, readsProcessed: 0 });
    setDataset(null);
    setIndex(null);
    setPreviewRecords(null);
    setError(null);
    recordCacheRef.current = new Map();
  }, []);

  const cancel = useCallback(() => {
    if (workerRef.current) {
      workerRef.current.terminate();
      workerRef.current = null;
    }
    setStatus("idle");
  }, []);

  const parseFile = useCallback((file) => {
    cancel();
    reset();
    fileRef.current = file;
    setStatus("parsing");
    setProgress({ bytesProcessed: 0, totalBytes: file.size, readsProcessed: 0 });

    const worker = new Worker(new URL("../workers/fastqParser.worker.js", import.meta.url), { type: "module" });
    workerRef.current = worker;

    worker.onmessage = (e) => {
      const msg = e.data;
      if (msg.type === "progress") {
        setProgress({ bytesProcessed: msg.bytesProcessed, totalBytes: msg.totalBytes, readsProcessed: msg.readsProcessed });
      } else if (msg.type === "done") {
        setDataset(msg.dataset);
        setIndex({ lengths: msg.lengths, meanQs: msg.meanQs, gcPerRead: msg.gcPerRead, ids: msg.ids, offsets: msg.offsets });
        setPreviewRecords(msg.previewRecords);
        setStatus("done");
        worker.terminate();
        workerRef.current = null;
      } else if (msg.type === "error") {
        setError(msg.message);
        setStatus("error");
        worker.terminate();
        workerRef.current = null;
      }
    };
    worker.onerror = (err) => {
      setError(err.message || "Worker error while parsing file.");
      setStatus("error");
    };
    worker.postMessage({ type: "parse", file });
  }, [cancel, reset]);

  /**
   * Fetch the full sequence/quality text for read `i` (0-indexed, in file
   * order). Uses byte-offset slicing for plain files (scales to any file
   * size) or the bounded in-memory preview for gzip files.
   */
  const getRecord = useCallback(async (i) => {
    if (recordCacheRef.current.has(i)) return recordCacheRef.current.get(i);
    if (!dataset) return null;

    if (dataset.parseMode === "offset" && index?.offsets && fileRef.current) {
      const start = index.offsets[i * 2];
      const len = index.offsets[i * 2 + 1];
      const blob = fileRef.current.slice(start, start + len);
      const text = await blob.text();
      const lines = text.split(/\r?\n/).filter((l) => l.length > 0 || l === "");
      const header = lines[0]?.replace(/^@/, "") || "";
      const spaceIdx = header.indexOf(" ");
      const id = spaceIdx === -1 ? header : header.slice(0, spaceIdx);
      const desc = spaceIdx === -1 ? "" : header.slice(spaceIdx + 1);
      const record = { id, desc, seq: lines[1] || "", qual: lines[3] || "" };
      recordCacheRef.current.set(i, record);
      return record;
    }

    if (previewRecords && i < previewRecords.length) {
      const record = previewRecords[i];
      recordCacheRef.current.set(i, record);
      return record;
    }

    return { unavailable: true };
  }, [dataset, index, previewRecords]);

  return { status, progress, dataset, index, previewRecords, error, parseFile, cancel, reset, getRecord };
}
