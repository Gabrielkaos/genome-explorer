import { useCallback, useRef, useState } from "react";

/**
 * Runs the pan-GWAS worker and holds its result. Local to the Association
 * section (nothing downstream consumes it), mirroring usePhylo's protocol.
 */
export function useAssoc() {
  const [status, setStatus] = useState("idle"); // idle | running | done | error
  const [stage, setStage] = useState(null);
  const [pct, setPct] = useState(0);
  const [detail, setDetail] = useState("");
  const [result, setResult] = useState(null);
  const [error, setError] = useState(null);
  const workerRef = useRef(null);

  const cancel = useCallback(() => {
    if (workerRef.current) { workerRef.current.terminate(); workerRef.current = null; }
    setStatus("idle");
    setPct(0);
    setStage(null);
    setDetail("");
  }, []);

  const run = useCallback((data) => {
    cancel();
    setStatus("running"); setStage("tests"); setPct(0); setDetail("");
    setResult(null); setError(null);

    const worker = new Worker(new URL("../workers/assoc.worker.js", import.meta.url), { type: "module" });
    workerRef.current = worker;
    worker.onmessage = (e) => {
      const msg = e.data;
      if (msg.type === "progress") { setStage(msg.stage); setPct(msg.pct); setDetail(msg.detail || ""); }
      else if (msg.type === "done") {
        setResult({ ...msg, params: data });
        setStatus("done");
        setPct(100); setStage("done"); setDetail("");
        worker.terminate(); workerRef.current = null;
      }
      else if (msg.type === "error") { setError(msg.message); setStatus("error"); worker.terminate(); workerRef.current = null; }
    };
    worker.onerror = (err) => { setError(err.message || "Worker error during association testing."); setStatus("error"); };
    worker.postMessage({ type: "assoc", data });
  }, [cancel]);

  return { status, stage, pct, detail, result, error, run, cancel };
}
