import { useCallback, useRef, useState } from "react";

export function useAssembly() {
  const [status, setStatus] = useState("idle"); // idle | running | done | error
  const [stage, setStage] = useState(null);
  const [pct, setPct] = useState(0);
  const [contigs, setContigs] = useState(null);
  const [stats, setStats] = useState(null);
  const [error, setError] = useState(null);
  const workerRef = useRef(null);

  const cancel = useCallback(() => {
    if (workerRef.current) { workerRef.current.terminate(); workerRef.current = null; }
    setStatus("idle");
  }, []);

  const run = useCallback((reads, params) => {
    cancel();
    setStatus("running"); setStage("indexing"); setPct(0); setContigs(null); setStats(null); setError(null);

    const worker = new Worker(new URL("../workers/assembly.worker.js", import.meta.url), { type: "module" });
    workerRef.current = worker;
    worker.onmessage = (e) => {
      const msg = e.data;
      if (msg.type === "progress") { setStage(msg.stage); setPct(msg.pct); }
      else if (msg.type === "done") { setContigs(msg.contigs); setStats(msg.stats); setStatus("done"); worker.terminate(); workerRef.current = null; }
      else if (msg.type === "error") { setError(msg.message); setStatus("error"); worker.terminate(); workerRef.current = null; }
    };
    worker.onerror = (err) => { setError(err.message || "Worker error during assembly."); setStatus("error"); };
    worker.postMessage({ type: "assemble", reads, params });
  }, [cancel]);

  return { status, stage, pct, contigs, stats, error, run, cancel };
}
