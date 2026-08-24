import { annotateGenome } from "../lib/annotation/genecaller.js";
import { buildGff3, buildGenbank, buildFaa, buildFfn, buildTsv } from "../lib/annotation/exportAnnotation.js";

self.onmessage = async (e) => {
  if (e.data?.type !== "annotate") return;
  try {
    await runAnnotation(e.data.contigs, e.data.params || {});
  } catch (err) {
    self.postMessage({ type: "error", message: String(err?.message || err) });
  }
};

async function runAnnotation(contigs, params) {
  const t0 = performance.now();

  // Yield to the event loop between heavy stages so progress messages flush.
  const tick = () => new Promise((r) => setTimeout(r, 0));

  const report = (m) => self.postMessage({ type: "progress", ...m });

  report({ stage: "training", pct: 5 });
  await tick();
  const result = annotateGenome(contigs, params, report);

  report({ stage: "exports", pct: 92 });
  await tick();
  const genesByContig = new Map();
  for (const g of result.genes) {
    if (!genesByContig.has(g.contigId)) genesByContig.set(g.contigId, []);
    genesByContig.get(g.contigId).push(g);
  }
  const exports = {
    gff3: buildGff3(contigs, result.genes),
    genbank: buildGenbank(contigs, genesByContig, { locusPrefix: params.locusPrefix ?? "GE" }),
    faa: buildFaa(result.genes),
    ffn: buildFfn(result.genes),
    tsv: buildTsv(result.genes),
  };

  result.stats.computeTimeMs = performance.now() - t0;
  report({ stage: "done", pct: 100 });
  self.postMessage({ type: "done", ...result, exports });
}
