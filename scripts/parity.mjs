/*
 * Parity harness: runs every engine twice (WASM kernel vs pure-JS fallback)
 * over deterministic generated datasets and asserts the results match.
 *
 *   npm run parity
 */
import { gxConfig } from "../src/lib/wasm/bridge.js";
import { runAssembly } from "../src/lib/engine/assemblyEngine.js";
import { runMsa } from "../src/lib/engine/msaEngine.js";
import { runPhylo } from "../src/lib/engine/phyloEngine.js";

/* deterministic RNG (same generator the phylo worker uses) */
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const BASES = "ACGT";
let failures = 0;

function check(name, ok, detail = "") {
  if (!ok) {
    failures++;
    console.error(`  FAIL ${name}${detail ? ` - ${detail}` : ""}`);
  } else {
    console.log(`  ok   ${name}`);
  }
}

function closeEnough(a, b, tol = 1e-9) {
  if (a === b) return true;
  if (Number.isNaN(a) && Number.isNaN(b)) return true;
  if (!Number.isFinite(a) || !Number.isFinite(b)) return false;
  return Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b));
}

function revcomp(s) {
  const map = { A: "T", C: "G", G: "C", T: "A" };
  return [...s].reverse().map((c) => map[c] ?? "N").join("");
}

/* ------------------------------------------------------------ datasets */

function makeReads() {
  const rng = mulberry32(101);
  const genomeLen = 12000;
  let genome = "";
  for (let i = 0; i < genomeLen; i++) genome += BASES[Math.floor(rng() * 4)];

  const reads = [];
  const readLen = 1400;
  const step = 220;
  for (let start = 0; start + readLen < genomeLen; start += step) {
    let seq = genome.slice(start, start + readLen);
    // sprinkle sequencing errors
    const err = [...seq];
    for (let i = 0; i < err.length; i++) {
      if (rng() < 0.015) err[i] = BASES[Math.floor(rng() * 4)];
    }
    seq = err.join("");
    const rev = rng() < 0.5;
    reads.push({ id: `read_${reads.length}`, seq: rev ? revcomp(seq) : seq });
  }
  return reads;
}

function makeAlignment(nSeq = 10, len = 260) {
  const rng = mulberry32(202);
  let ancestor = "";
  for (let i = 0; i < len; i++) ancestor += BASES[Math.floor(rng() * 4)];
  const records = [];
  for (let s = 0; s < nSeq; s++) {
    const rate = 0.02 + 0.05 * (s / nSeq);
    let seq = "";
    for (let i = 0; i < len; i++) {
      if (rng() < rate) seq += BASES[Math.floor(rng() * 4)];
      else seq += ancestor[i];
    }
    records.push({ id: `taxon_${s}`, seq });
  }
  // add a shared gap column block to exercise gap handling
  for (const r of records) r.seq = r.seq.slice(0, 40) + "--" + r.seq.slice(42);
  return records;
}

function makePhyloRecords(nTaxa = 7, len = 320) {
  const rng = mulberry32(303);
  let ancestor = "";
  for (let i = 0; i < len; i++) ancestor += BASES[Math.floor(rng() * 4)];
  const records = [];
  for (let t = 0; t < nTaxa; t++) {
    const rate = 0.03 + 0.09 * (t / nTaxa);
    let seq = "";
    for (let i = 0; i < len; i++) {
      if (rng() < rate) seq += BASES[Math.floor(rng() * 4)];
      else seq += ancestor[i];
    }
    records.push({ id: `t${t}`, seq });
  }
  return records;
}

/* ---------------------------------------------------------------- tests */

async function testAssembly() {
  console.log("assembly");
  const reads = makeReads();
  const params = {};

  gxConfig.disableWasm = true;
  const jsRes = await runAssembly(reads, params, () => {});
  gxConfig.disableWasm = false;
  const wmRes = await runAssembly(reads, params, () => {});

  check("contig count", jsRes.contigs.length === wmRes.contigs.length,
        `js=${jsRes.contigs.length} wasm=${wmRes.contigs.length}`);
  check("total length", jsRes.stats.totalLength === wmRes.stats.totalLength);
  const nCmp = Math.min(jsRes.contigs.length, wmRes.contigs.length);
  let seqsEq = true, membersEq = true, circEq = true;
  for (let i = 0; i < nCmp; i++) {
    const a = jsRes.contigs[i], b = wmRes.contigs[i];
    seqsEq &&= a.seq === b.seq;
    circEq &&= a.circular === b.circular;
    membersEq &&= a.members.length === b.members.length &&
      a.members.every((m, k) =>
        m.readIdx === b.members[k].readIdx && m.strand === b.members[k].strand &&
        m.contigStart === b.members[k].contigStart && m.contigEnd === b.members[k].contigEnd);
  }
  check("sequences byte-identical", seqsEq);
  check("circularity flags", circEq);
  check("member layouts", membersEq);

  const statKeys = ["inputReads", "containedReads", "usedReads", "unplacedReads",
                    "overlapsFound", "numContigs", "totalLength", "longestContig",
                    "circularContigs", "meanOverlapScore", "n50"];
  let statsEq = true, badKey = "";
  for (const key of statKeys) {
    if (!closeEnough(jsRes.stats[key], wmRes.stats[key])) { statsEq = false; badKey = key; break; }
  }
  check("stats", statsEq, badKey);
}

async function testMsa() {
  console.log("msa");
  const records = makeAlignment();
  const paramSets = [
    {},
    { kmerSize: 3 },
    { gapOpen: 12, gapExtend: 2, match: 3, mismatch: 2, gapGapBonus: 1 },
  ];

  for (let pi = 0; pi < paramSets.length; pi++) {
    const params = paramSets[pi];
    gxConfig.disableWasm = true;
    const jsRes = await runMsa(records, params, () => {});
    gxConfig.disableWasm = false;
    const wmRes = await runMsa(records, params, () => {});

    const tag = `[params ${pi}]`;
    check(`${tag} length`, jsRes.length === wmRes.length, `js=${jsRes.length} wasm=${wmRes.length}`);
    check(`${tag} rows identical`, jsRes.rows.join("#") === wmRes.rows.join("#"));
    check(`${tag} kUsed`, jsRes.stats.kUsed === wmRes.stats.kUsed);
    const ctrKeys = ["conservedColumns", "variableColumns", "informativeColumns",
                     "singletonColumns", "gapColumns", "gapFraction",
                     "meanPairIdentity", "minPairIdentity", "maxPairIdentity"];
    let statsEq = true, badKey = "";
    for (const key of ctrKeys) {
      if (!closeEnough(jsRes.stats[key], wmRes.stats[key])) { statsEq = false; badKey = key; break; }
    }
    check(`${tag} counters`, statsEq, badKey);
    check(`${tag} variable site count`,
          jsRes.stats.totalVariableSites === wmRes.stats.totalVariableSites);
    check(`${tag} consensus codes`,
          String.fromCharCode(...jsRes.stats.consensusCodes) === String.fromCharCode(...wmRes.stats.consensusCodes));
    let sitesEq = JSON.stringify(jsRes.stats.variableSites.slice(0, 50)) ===
                 JSON.stringify(wmRes.stats.variableSites.slice(0, 50));
    check(`${tag} variable sites (first 50)`, sitesEq);
  }
}

async function testPhylo() {
  console.log("phylo");
  const records = makePhyloRecords();
  const cases = [
    { method: "nj", model: "jc", bootstrap: 20, seed: 7 },
    { method: "nj", model: "k2p", gapMode: "complete", bootstrap: 10, seed: 99 },
    { method: "upgma", model: "p", bootstrap: 8, seed: 5 },
    { method: "nj", model: "jc", bootstrap: 0, rooting: "none" },
  ];

  for (let ci = 0; ci < cases.length; ci++) {
    const params = cases[ci];
    gxConfig.disableWasm = true;
    const jsRes = await runPhylo(records, params, () => {});
    gxConfig.disableWasm = false;
    const wmRes = await runPhylo(records, params, () => {});

    const tag = `[case ${ci}]`;
    check(`${tag} newick identical`, jsRes.newick === wmRes.newick,
          `\n    js  : ${jsRes.newick}\n    wasm: ${wmRes.newick}`);
    check(`${tag} columnsUsed`, jsRes.columnsUsed === wmRes.columnsUsed);

    const n = jsRes.numTaxa;
    let matEq = true;
    for (let i = 0; i < n * n && matEq; i++)
      matEq = closeEnough(jsRes.matrix[i], wmRes.matrix[i], 1e-10);
    check(`${tag} distance matrix`, matEq);

    let sumEq = closeEnough(jsRes.summary.meanDistance, wmRes.summary.meanDistance, 1e-10) &&
                closeEnough(jsRes.summary.tiTvRatio ?? NaN, wmRes.summary.tiTvRatio ?? NaN, 1e-10) &&
                jsRes.summary.saturatedPairs === wmRes.summary.saturatedPairs &&
                jsRes.summary.noDataPairs === wmRes.summary.noDataPairs &&
                jsRes.summary.clampedLimbs === wmRes.summary.clampedLimbs;
    check(`${tag} summary`, sumEq);
  }
}

/* ----------------------------------------------------------------- main */

async function main() {
  console.log("WASM/JS parity checks");
  try {
    await testAssembly();
    await testMsa();
    await testPhylo();
  } catch (err) {
    failures++;
    console.error("UNEXPECTED ERROR:", err);
  }
  if (failures) {
    console.error(`\n${failures} check(s) FAILED`);
    process.exit(1);
  }
  console.log("\nall parity checks passed");
  process.exit(0);
}

main();
