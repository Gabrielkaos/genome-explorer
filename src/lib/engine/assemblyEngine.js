/*
 * Assembly engine: WASM-first with automatic fallback to the original JS
 * pipeline. Both paths emit identical progress messages and identical
 * result objects, so the worker and UI cannot tell which one ran.
 */
import { gxWasm, mallocIntoHeap, takeBlob, withProgress, BlobReader } from "../wasm/bridge.js";
import { buildMinimizerIndex, findCandidatePairs, estimateOverlap } from "../assembly/overlap.js";
import { buildGraph, extractPaths } from "../assembly/graph.js";
import { buildAllContigs } from "../assembly/consensus.js";
import { computeNStat } from "../fastq/stats.js";

const DEFAULTS = { k: 10, w: 5, maxOccurrence: 40, diagTolerance: 60, minMatches: 4, minOverlapLen: 200 };

/* ------------------------------------------------------------------ JS path */

async function runAssemblyJs(reads, params, report) {
  report({ stage: "indexing", pct: 5 });
  const index = buildMinimizerIndex(reads, params);

  report({ stage: "indexing", pct: 35 });
  const pairs = findCandidatePairs(index);

  report({ stage: "overlaps", pct: 45 });
  const overlaps = [];
  let processed = 0;
  const total = pairs.size;
  for (const [key, matches] of pairs) {
    const [i, j] = key.split(",").map(Number);
    const ov = estimateOverlap(matches, reads[i].seq.length, reads[j].seq.length, params);
    if (ov) overlaps.push({ a: i, b: j, ...ov });
    processed++;
    if (processed % 200 === 0) {
      report({ stage: "overlaps", pct: 45 + Math.round((processed / Math.max(1, total)) * 30) });
    }
  }

  report({ stage: "graph", pct: 80 });
  const graph = buildGraph(reads, overlaps);
  const paths = extractPaths(reads, graph);

  report({ stage: "consensus", pct: 90 });
  return { contigs: buildAllContigs(reads, paths), overlaps, graph };
}

function summarize(reads, params, contigs, overlaps, graph, t0) {
  const n = reads.length;
  const lengths = contigs.map((c) => c.length);
  const totalLength = lengths.reduce((a, b) => a + b, 0);
  const usedReads = contigs.reduce((a, c) => a + c.readCount, 0);
  return {
    inputReads: n,
    containedReads: graph.contained.size,
    usedReads,
    unplacedReads: n - graph.contained.size - usedReads,
    overlapsFound: overlaps.length,
    numContigs: contigs.length,
    totalLength,
    n50: computeNStat(lengths, 0.5),
    longestContig: lengths.length ? Math.max(...lengths) : 0,
    circularContigs: contigs.filter((c) => c.circular).length,
    meanOverlapScore: overlaps.length ? overlaps.reduce((a, o) => a + o.score, 0) / overlaps.length : 0,
    computeTimeMs: performance.now() - t0,
    params: {
      k: params.k ?? DEFAULTS.k,
      w: params.w ?? DEFAULTS.w,
      minOverlapLen: params.minOverlapLen ?? DEFAULTS.minOverlapLen,
      minMatches: params.minMatches ?? DEFAULTS.minMatches,
    },
  };
}

/* ---------------------------------------------------------------- WASM path */

/*
 * Blob layout (little-endian), written by gx_assembly_run:
 *   "GXAR" u32 version=1 u32 nContigs
 *   per contig: u32 seqLen, u32 circular, u32 nMembers, u32 reserved,
 *               seqLen bytes, pad to 4, then nMembers x
 *               { u32 readIdx, i32 strand, u32 contigStart, u32 contigEnd,
 *                 u32 trimmedFromStart }
 *   tail: f64 n50, f64 meanOverlapScore,
 *         u32 inputReads, containedReads, usedReads, unplacedReads,
 *             overlapsFound, numContigs(again), totalLength, longestContig,
 *         u32 circularContigs
 */
const textDecoder = new TextDecoder("ascii");

function decodeAssemblyBlob(buf, reads, params, t0) {
  const r = new BlobReader(buf);
  r.magic("GXAR");
  r.u32(); // version
  const nContigs = r.u32();

  const contigs = [];
  for (let i = 0; i < nContigs; i++) {
    const seqLen = r.u32();
    const circular = r.u32();
    const nMembers = r.u32();
    r.u32(); // reserved
    const seq = textDecoder.decode(r.bytes(seqLen));
    r.align(4);
    const members = [];
    for (let mm = 0; mm < nMembers; mm++) {
      members.push({
        readIdx: r.u32(),
        strand: r.i32(),
        contigStart: r.u32(),
        contigEnd: r.u32(),
        // trimmedFromStart: consumed to keep the stream aligned
      });
      r.u32();
    }
    contigs.push({
      seq,
      length: seqLen,
      circular: circular === 1,
      readCount: members.length,
      members,
    });
  }

  const n50 = r.f64();
  const meanOverlapScore = r.f64();
  const inputReads = r.u32();
  const containedReads = r.u32();
  const usedReads = r.u32();
  const unplacedReads = r.u32();
  const overlapsFound = r.u32();
  r.u32(); // numContigs duplicate
  const totalLength = r.u32();
  const longestContig = r.u32();
  const circularContigs = r.u32();

  const stats = {
    inputReads,
    containedReads,
    usedReads,
    unplacedReads,
    overlapsFound,
    numContigs: nContigs,
    totalLength,
    n50,
    longestContig,
    circularContigs,
    meanOverlapScore,
    computeTimeMs: performance.now() - t0,
    params: {
      k: params.k ?? DEFAULTS.k,
      w: params.w ?? DEFAULTS.w,
      minOverlapLen: params.minOverlapLen ?? DEFAULTS.minOverlapLen,
      minMatches: params.minMatches ?? DEFAULTS.minMatches,
    },
  };
  return { contigs, stats };
}

async function runAssemblyWasm(reads, params, report) {
  const mod = await gxWasm();

  const encoder = new TextEncoder();
  let flatLen = 0;
  for (const rd of reads) flatLen += rd.seq.length;
  const flat = new Uint8Array(flatLen);
  const offs = new Uint32Array(reads.length + 1);
  let p = 0;
  reads.forEach((rd, i) => {
    flat.set(encoder.encode(rd.seq), p);
    p += rd.seq.length;
    offs[i + 1] = p;
  });

  const t0 = performance.now();
  const { k = DEFAULTS.k, w = DEFAULTS.w, maxOccurrence = DEFAULTS.maxOccurrence,
          diagTolerance = DEFAULTS.diagTolerance, minMatches = DEFAULTS.minMatches,
          minOverlapLen = DEFAULTS.minOverlapLen } = params;

  const out = await withProgress(() => {
    const fp = mallocIntoHeap(mod, flat);
    const op = mallocIntoHeap(mod, offs);
    try {
      mod._gx_assembly_run(fp, op, reads.length, k, w, maxOccurrence, diagTolerance, minMatches, minOverlapLen);
      return takeBlob(mod, mod._gx_assembly_result_ptr(), mod._gx_assembly_result_len());
    } finally {
      mod._gx_assembly_free();
      mod._free(fp);
      mod._free(op);
    }
  }, (code, pct) => {
    if (code <= 2 || code === 4 || code === 5) {
      report({ stage: ["indexing", "indexing", "overlaps", "overlaps", "graph", "consensus"][code], pct });
    } else if (code === 3) {
      report({ stage: "overlaps", pct });
    }
  });

  return decodeAssemblyBlob(out, reads, params, t0);
}

/* ------------------------------------------------------------------ public */

export async function runAssembly(reads, params, report) {
  report = report || (() => {});
  try {
    const res = await runAssemblyWasm(reads, params, report);
    return res;
  } catch (err) {
    console.warn("[assembly] WASM engine unavailable, using JS fallback:", err?.message || err);
    const t0 = performance.now();
    const { contigs, overlaps, graph } = await runAssemblyJs(reads, params, report);
    return { contigs, stats: summarize(reads, params, contigs, overlaps, graph, t0) };
  }
}
