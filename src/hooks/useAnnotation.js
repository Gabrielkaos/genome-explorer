import { useCallback, useRef, useState } from "react";

export function useAnnotation() {
  const [status, setStatus] = useState("idle"); // idle | running | done | error
  const [stage, setStage] = useState(null);
  const [pct, setPct] = useState(0);
  const [genes, setGenes] = useState(null);
  const [stats, setStats] = useState(null);
  const [gcTracks, setGcTracks] = useState(null);
  const [modelInfo, setModelInfo] = useState(null);
  const [exports, setExports] = useState(null);
  const [error, setError] = useState(null);
  const workerRef = useRef(null);

  const cancel = useCallback(() => {
    if (workerRef.current) { workerRef.current.terminate(); workerRef.current = null; }
    setStatus("idle");
  }, []);

  const run = useCallback((contigs, params) => {
    cancel();
    setStatus("running"); setStage("training"); setPct(0);
    setGenes(null); setStats(null); setGcTracks(null); setModelInfo(null); setExports(null); setError(null);

    const worker = new Worker(new URL("../workers/annotation.worker.js", import.meta.url), { type: "module" });
    workerRef.current = worker;
    worker.onmessage = (e) => {
      const msg = e.data;
      if (msg.type === "progress") { setStage(msg.stage); setPct(msg.pct); }
      else if (msg.type === "done") {
        setGenes(msg.genes); setStats(msg.stats); setGcTracks(msg.gcTracks);
        setModelInfo(msg.modelInfo); setExports(msg.exports);
        setStatus("done");
        worker.terminate(); workerRef.current = null;
      }
      else if (msg.type === "error") { setError(msg.message); setStatus("error"); worker.terminate(); workerRef.current = null; }
    };
    worker.onerror = (err) => { setError(err.message || "Worker error during annotation."); setStatus("error"); };
    worker.postMessage({ type: "annotate", contigs, params });
  }, [cancel]);

  return { status, stage, pct, genes, stats, gcTracks, modelInfo, exports, error, run, cancel };
}
