import { runMsa } from "../lib/engine/msaEngine.js";
import { decodeChars } from "../lib/msa/alphabet.js";
import {
  buildAlignedFasta, buildClustal, buildNexus, buildPhylip,
  buildVariantsTsv, buildConsensusFasta,
} from "../lib/msa/exportMsa.js";

self.onmessage = async (e) => {
  if (e.data?.type !== "msa") return;
  try {
    const records = e.data.records;
    const params = e.data.params || {};
    const report = (m) => self.postMessage({ type: "progress", ...m });

    const { length, rows, stats } = await runMsa(records, params, report);

    report({ stage: "exports", pct: 96, detail: "Building FASTA / Clustal / NEXUS / PHYLIP exports" });
    const ids = records.map((r) => r.id);
    const consensusStr = decodeChars(stats.consensusCodes).replace(/-/g, "N");

    const { identityMatrix, variableSites, ...reportableStats } = stats;
    const exports = {
      fasta: buildAlignedFasta(ids, rows),
      clustal: buildClustal(ids, rows),
      nexus: buildNexus(ids, rows),
      phylip: buildPhylip(ids, rows),
      variantsTsv: buildVariantsTsv(variableSites),
      consensusFasta: buildConsensusFasta(ids, consensusStr),
    };

    report({ stage: "done", pct: 100 });
    self.postMessage({
      type: "done",
      ids,
      rows,
      length,
      stats: {
        ...reportableStats,
        variableSites: variableSites.slice(0, 5000),
        totalVariableSites: stats.totalVariableSites ?? variableSites.length,
        kUsed: stats.kUsed,
        timings: stats.timings,
      },
      identityMatrix: Array.from(identityMatrix),
      exports,
    });
  } catch (err) {
    self.postMessage({ type: "error", message: String(err?.message || err) });
  }
};
