import { runPhylo } from "../lib/engine/phyloEngine.js";

self.onmessage = async (e) => {
  if (e.data?.type !== "phylo") return;
  try {
    const records = e.data.records;
    const params = e.data.params || {};
    const report = (m) => self.postMessage({ type: "progress", ...m });
    const result = await runPhylo(records, params, report);
    report({ stage: "done", pct: 100 });
    self.postMessage({ type: "done", ...result });
  } catch (err) {
    self.postMessage({ type: "error", message: String(err?.message || err) });
  }
};
