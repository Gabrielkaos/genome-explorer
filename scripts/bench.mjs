/*
 * Benchmark: WASM kernels vs pure-JS engines on realistic data sizes.
 *
 *   npm run bench
 */
import { gxConfig } from "../src/lib/wasm/bridge.js";
import { runAssembly } from "../src/lib/engine/assemblyEngine.js";
import { runMsa } from "../src/lib/engine/msaEngine.js";
import { runPhylo } from "../src/lib/engine/phyloEngine.js";

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const BASES = "ACGT";
const revcomp = (s) => [...s].reverse().map((c) => ({ A: "T", C: "G", G: "C", T: "A" }[c] ?? "N")).join("");

function makeReads(nReads, readLen, genomeLen, errRate = 0.02) {
  const rng = mulberry32(11);
  let genome = "";
  for (let i = 0; i < genomeLen; i++) genome += BASES[Math.floor(rng() * 4)];
  const reads = [];
  for (let r = 0; r < nReads; r++) {
    const start = Math.floor(rng() * (genomeLen - readLen));
    const err = [...genome.slice(start, start + readLen)];
    for (let i = 0; i < err.length; i++) if (rng() < errRate) err[i] = BASES[Math.floor(rng() * 4)];
    const seq = err.join("");
    reads.push({ id: `r${r}`, seq: rng() < 0.5 ? revcomp(seq) : seq });
  }
  return reads;
}

function makeAlignment(nSeq, len) {
  const rng = mulberry32(22);
  let ancestor = "";
  for (let i = 0; i < len; i++) ancestor += BASES[Math.floor(rng() * 4)];
  const records = [];
  for (let s = 0; s < nSeq; s++) {
    const rate = 0.02 + 0.08 * (s / nSeq);
    let seq = "";
    for (let i = 0; i < len; i++) seq += rng() < rate ? BASES[Math.floor(rng() * 4)] : ancestor[i];
    records.push({ id: `s${s}`, seq });
  }
  return records;
}

function makePhylo(nTaxa, len) {
  const rng = mulberry32(33);
  let ancestor = "";
  for (let i = 0; i < len; i++) ancestor += BASES[Math.floor(rng() * 4)];
  const records = [];
  for (let t = 0; t < nTaxa; t++) {
    const rate = 0.03 + 0.12 * (t / nTaxa);
    let seq = "";
    for (let i = 0; i < len; i++) seq += rng() < rate ? BASES[Math.floor(rng() * 4)] : ancestor[i];
    records.push({ id: `t${t}`, seq });
  }
  return records;
}

async function timeBoth(label, fnJs, fnWasm) {
  gxConfig.disableWasm = true;
  const tJs0 = performance.now();
  await fnJs();
  const jsMs = performance.now() - tJs0;

  gxConfig.disableWasm = false;
  const tW0 = performance.now();
  await fnWasm();
  const wMs = performance.now() - tW0;

  console.log(`${label.padEnd(34)} JS ${jsMs.toFixed(0).padStart(7)} ms   WASM ${wMs.toFixed(0).padStart(7)} ms   ${(jsMs / wMs).toFixed(1)}x`);
}

async function main() {
  console.log("benchmark (node; browser workers add transfer overhead but same compute)\n");

  // Assembly: ~1.6 Mbp of reads (200 x 3kb over a 300kb genome)
  {
    const reads = makeReads(200, 3000, 300_000);
    await timeBoth(
      `assembly ${reads.length} reads x 3kb`,
      () => runAssembly(reads, {}, () => {}),
      () => runAssembly(reads, {}, () => {})
    );
  }

  // MSA: 40 sequences x 1200 bp
  {
    const recs = makeAlignment(40, 1200);
    await timeBoth(
      `msa ${recs.length} seqs x 1200bp`,
      () => runMsa(recs, {}, () => {}),
      () => runMsa(recs, {}, () => {})
    );
  }

  // Phylo: 60 taxa x 2000 bp with 100 bootstrap replicates
  {
    const recs = makePhylo(60, 2000);
    const params = { bootstrap: 100, seed: 42 };
    await timeBoth(
      `phylo ${recs.length} taxa x 2000bp +100 reps`,
      () => runPhylo(recs, params, () => {}),
      () => runPhylo(recs, params, () => {})
    );
  }
}

main();
