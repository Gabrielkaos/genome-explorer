/*
 * MSA engine: WASM-first with automatic fallback to the original JS
 * progressive-alignment pipeline. Both paths produce identical worker
 * messages and result objects.
 */
import { gxWasm, mallocIntoHeap, takeBlob, BlobReader } from "../wasm/bridge.js";
import { encodeSeq, decodeChars, CODE_CHARS } from "../msa/alphabet.js";
import { kmerDistanceMatrix } from "../msa/kmer.js";
import { upgma } from "../msa/guidetree.js";
import { singleProfile, mergeProfiles } from "../msa/profile.js";
import { computeAlnStats } from "../msa/alnstats.js";

/* ------------------------------------------------------------------ JS path */

export function progressiveAlignJs(records, params = {}, report = () => {}) {
  const t0 = performance.now();
  const n = records.length;
  if (n < 2) throw new Error("Need at least two sequences to align.");

  report({ stage: "distances", pct: 4, detail: `Encoding ${n} sequences` });
  const encoded = records.map((r) => encodeSeq(r.seq));

  const tK = performance.now();
  report({ stage: "distances", pct: 8, detail: "Estimating pairwise k-mer distances" });
  const { dist, k } = kmerDistanceMatrix(encoded, params.kmerSize ?? 6);
  const distMs = performance.now() - tK;

  report({ stage: "tree", pct: 12, detail: "Building UPGMA guide tree" });
  const tree = upgma(dist, n);

  let mergesDone = 0;
  const totalMerges = n - 1;
  const tA = performance.now();

  function build(node) {
    if (node.leaf) return singleProfile(encoded[node.idx], node.idx);
    const left = build(node.left);
    const right = build(node.right);
    mergesDone++;
    if (mergesDone % Math.max(1, Math.ceil(totalMerges / 40)) === 0 || mergesDone === totalMerges) {
      report({
        stage: "align",
        pct: 14 + Math.round((mergesDone / totalMerges) * 72),
        detail: `Progressive merge ${mergesDone}/${totalMerges}`,
      });
    }
    return mergeProfiles(left, right, params);
  }

  const root = build(tree);
  const alignMs = performance.now() - tA;

  const chars = new Array(n);
  for (const mem of root.members) chars[mem.origIdx] = mem.chars;

  return {
    length: root.len,
    chars,
    cols: root.cols,
    dist,
    kUsed: k,
    timings: {
      distancesMs: distMs,
      alignmentMs: alignMs,
      totalMs: performance.now() - t0,
    },
  };
}

function runMsaJs(records, params, report) {
  const result = progressiveAlignJs(records, params, report);
  report({ stage: "stats", pct: 90, detail: "Computing alignment statistics" });
  const stats = computeAlnStats(result);
  stats.kUsed = result.kUsed;
  stats.timings = result.timings;
  stats.totalVariableSites = stats.variableSites.length;
  return { length: result.length, rows: result.chars.map(decodeChars), stats };
}

/* ---------------------------------------------------------------- WASM path */

const DEFAULTS = { match: 2, mismatch: 1, gapOpen: 8, gapExtend: 1, gapGapBonus: 0.5 };

/*
 * Blob layout (little-endian), written by gx_msa_run:
 *   "GXMR" u32 version u32 n u32 L
 *   f64 distancesMs f64 alignmentMs
 *   i32 kUsed (-1 = null) u32 reserved
 *   f64 dist[n*n]
 *   u32 cols[5*L]
 *   u8 rows[n*L] (input order), pad to 4
 *   f64 identity[n*n]
 *   u8 consensus[L] u8 colClass[L] pad to 8
 *   f64 perGapPct[n]
 *   f64 gapFraction meanPairIdentity minPairIdentity maxPairIdentity
 *   u32 conservedColumns variableColumns informativeColumns singletonColumns
 *       gapColumns variableSitesTotal
 *
 * colClass bits: 1 conserved, 2 variable, 4 informative, 8 singleton, 16 all-gap.
 */
async function runMsaWasm(records, params, report) {
  const mod = await gxWasm();
  const n = records.length;

  const encoded = records.map((r) => encodeSeq(r.seq));
  const flatLen = encoded.reduce((a, e) => a + e.length, 0);
  const flat = new Uint8Array(flatLen);
  const lens = new Uint32Array(n);
  let p = 0;
  encoded.forEach((e, i) => {
    flat.set(e, p);
    p += e.length;
    lens[i] = e.length;
  });

  const {
    kmerSize = 6,
    match = DEFAULTS.match,
    mismatch = DEFAULTS.mismatch,
    gapOpen = DEFAULTS.gapOpen,
    gapExtend = DEFAULTS.gapExtend,
    gapGapBonus = DEFAULTS.gapGapBonus,
  } = params;

  const t0 = performance.now();
  report({ stage: "distances", pct: 4, detail: `Encoding ${n} sequences` });
  report({ stage: "distances", pct: 8, detail: "Estimating pairwise k-mer distances" });
  report({ stage: "tree", pct: 12, detail: "Building UPGMA guide tree" });

  const fp = mallocIntoHeap(mod, flat);
  const lp = mallocIntoHeap(mod, lens);
  try {
    const rc = mod._gx_msa_run(fp, lp, n, kmerSize, match, mismatch, gapOpen, gapExtend, gapGapBonus);
    if (rc === -1) throw new Error("Cannot align an empty profile.");
    if (rc === -2) throw new Error(
      "Profile DP too large for the WASM engine. Try fewer or shorter sequences."
    );
    if (rc !== 0) throw new Error(`WASM alignment failed (code ${rc}).`);
  } finally {
    // blob stays valid until the next _gx_*_run; copy before freeing inputs
  }
  const outBlob = takeBlob(mod, mod._gx_msa_result_ptr(), mod._gx_msa_result_len());
  mod._gx_msa_free();
  mod._free(fp);
  mod._free(lp);

  report({ stage: "stats", pct: 90, detail: "Computing alignment statistics" });

  const r = new BlobReader(outBlob);
  r.magic("GXMR");
  r.u32(); // version
  r.u32(); // n
  const L = r.u32();
  const distancesMs = r.f64();
  const alignmentMs = r.f64();
  const kUsedRaw = r.i32();
  r.u32(); // reserved

  r.skip(8 * n * n); // distance matrix (not consumed by the worker)
  const colsOff = r.pos;
  r.skip(20 * L);
  const rowsOff = r.pos;
  r.skip(n * L);
  r.align(8); /* identity matrix is f64-aligned by the writer */
  const identOff = r.pos;
  r.skip(8 * n * n);
  const consensusCodes = r.bytes(L).slice();
  const colClass = r.bytes(L).slice();
  r.align(8);
  const perGapOff = r.pos;
  r.skip(8 * n);
  const gapFraction = r.f64();
  const meanPairIdentity = r.f64();
  const minPairIdentity = r.f64();
  const maxPairIdentity = r.f64();
  const conservedColumns = r.u32();
  const variableColumns = r.u32();
  const informativeColumns = r.u32();
  const singletonColumns = r.u32();
  const gapColumns = r.u32();
  const variableSitesTotal = r.u32();

  const heapU32 = new Uint32Array(outBlob.buffer, outBlob.byteOffset, Math.ceil(outBlob.length / 4));

  // per-sequence aligned rows as strings, input order preserved by the kernel
  const rows = new Array(n);
  for (let i = 0; i < n; i++) {
    const seg = outBlob.subarray(rowsOff + i * L, rowsOff + (i + 1) * L);
    rows[i] = decodeChars(seg);
  }

  // rebuild the (capped) variable-site list exactly like alnstats.js did
  const MAX_SITES = 5000;
  const variableSites = [];
  for (let c = 0; c < L && variableSites.length < MAX_SITES; c++) {
    if (!(colClass[c] & 2)) continue;
    const counts = [heapU32[(colsOff >> 2) + c * 5], heapU32[(colsOff >> 2) + c * 5 + 1],
                    heapU32[(colsOff >> 2) + c * 5 + 2], heapU32[(colsOff >> 2) + c * 5 + 3]];
    const gapCount = heapU32[(colsOff >> 2) + c * 5 + 4];
    const present = [];
    for (let a = 0; a < 4; a++) if (counts[a]) present.push(a);
    variableSites.push({
      pos: c,
      consensus: CODE_CHARS[consensusCodes[c]],
      alleles: present.map((a) => `${CODE_CHARS[a]}:${counts[a]}`),
      informative: !!(colClass[c] & 4),
      indel: gapCount > 0,
    });
  }

  const identityMatrix = new Float64Array(n * n);
  identityMatrix.set(new Float64Array(outBlob.buffer, outBlob.byteOffset + identOff, n * n));
  const perSeqGapPct = Array.from(new Float64Array(outBlob.buffer, outBlob.byteOffset + perGapOff, n));

  const stats = {
    length: L,
    numSeqs: n,
    conservedColumns,
    variableColumns,
    informativeColumns,
    singletonColumns,
    gapColumns,
    gapFraction,
    meanPairIdentity,
    minPairIdentity,
    maxPairIdentity,
    identityMatrix,
    perSeqGapPct,
    consensusCodes,
    variableSites,
    totalVariableSites: variableSitesTotal,
    kUsed: kUsedRaw < 0 ? null : kUsedRaw,
    timings: {
      distancesMs,
      alignmentMs,
      totalMs: performance.now() - t0,
    },
  };

  return { length: L, rows, stats };
}

/* ------------------------------------------------------------------ public */

/**
 * Returns { length, rows: string[], stats } where stats carries every field
 * msa.worker.js posts (identityMatrix included; variableSites capped at 5000
 * on the WASM path, with totalVariableSites always the true count).
 */
export async function runMsa(records, params, report) {
  report = report || (() => {});
  try {
    return await runMsaWasm(records, params, report);
  } catch (err) {
    console.warn("[msa] WASM engine unavailable, using JS fallback:", err?.message || err);
    return runMsaJs(records, params, report);
  }
}

/** Distance matrix access for parity tests (JS path returns Float64Array too). */
export function jsProgressiveAlign(records, params, report) {
  return progressiveAlignJs(records, params, report);
}
