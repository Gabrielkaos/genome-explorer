/*
 * Phylo engine: WASM-first with automatic fallback to the original JS
 * pipeline (distance models -> NJ/UPGMA -> bootstrap -> Newick). Both paths
 * emit identical progress messages and identical result objects.
 *
 * The WASM bootstrap runs in chunks (gx_phylo_bootstrap_chunk) so the worker
 * can yield between batches and keep the progress bar moving - replicating
 * the tick()-based yielding the JS pipeline relied on.
 */
import { gxWasm, mallocIntoHeap, takeBlob, BlobReader } from "../wasm/bridge.js";
import { computeDistanceMatrix, encodeRows, MODELS } from "../phylo/distances.js";
import { neighborJoining } from "../phylo/nj.js";
import { upgmaTree } from "../phylo/upgma.js";
import { annotateSupports, collectSplits, indexAndAssignMembers, midpointRoot } from "../phylo/splits.js";
import { toNewick } from "../phylo/newick.js";

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const MODEL_IDS = { p: 0, jc: 1, k2p: 2 };

/* ------------------------------------------------------------------ shared */

function resolveOpts(params) {
  const method = params.method === "upgma" ? "upgma" : "nj";
  const opts = {
    model: MODELS[params.model] ? params.model : "jc",
    gapMode: params.gapMode === "complete" ? "complete" : "pairwise",
  };
  const rooting = method === "nj" && params.rooting === "none" ? "none" : "midpoint";
  return { method, opts, rooting };
}

function validateRecords(records) {
  const ids = records.map((r) => r.id);
  const n = ids.length;
  if (n < 3) throw new Error(`Phylogenetic inference needs at least 3 sequences (got ${n}).`);
  const L = records[0].seq.length;
  for (const r of records) {
    if (r.seq.length !== L) {
      throw new Error(
        `Sequences are not aligned: "${r.id}" is ${r.seq.length} bp but the first sequence is ${L} bp. ` +
        `Build a multiple alignment first (Alignment section), then infer the tree on its output.`
      );
    }
  }
  return { ids, n, L };
}

function clearSupports(node) {
  node.support = undefined;
  node.children.forEach(clearSupports);
}

/* ------------------------------------------------------------------ JS path */

async function runPhyloJs(records, params, report, opts, rooting, ids, n) {
  const tick = () => new Promise((r) => setTimeout(r, 0));

  report({ stage: "encode", pct: 4, detail: "Encoding alignment" });
  await tick();
  const codes = encodeRows(records.map((r) => r.seq));

  report({
    stage: "distance", pct: 9,
    detail: `${((n * (n - 1)) / 2)} pairwise distances · ${MODELS[opts.model].short} · ${opts.gapMode} deletion`,
  });
  await tick();
  const distRes = computeDistanceMatrix(codes, opts);

  report({ stage: "tree", pct: 20, detail: opts.method === "upgma" ? "UPGMA clustering" : "Neighbor-Joining search" });
  await tick();
  const main = opts.method === "upgma"
    ? upgmaTree(distRes.matrix, n, ids)
    : neighborJoining(distRes.matrix, n, ids);

  let tree = main.root;
  if (opts.method === "nj" && rooting === "midpoint") tree = midpointRoot(tree);
  indexAndAssignMembers(tree);

  const reps = Math.max(0, Math.min(1000, Math.round(params.bootstrap ?? 100)));
  const tally = new Map();
  let clampedLimbs = main.clampedLimbs;
  if (reps > 0) {
    const rng = mulberry32(Number.isFinite(+params.seed) ? +params.seed : 42);
    for (let r = 0; r < reps; r++) {
      const sample = new Uint32Array(distRes.columnsUsed);
      for (let k = 0; k < sample.length; k++) sample[k] = Math.floor(rng() * distRes.columnsUsed);
      const rep = opts.method === "upgma"
        ? upgmaTree(computeDistanceMatrix(codes, opts, sample).matrix, n, ids)
        : neighborJoining(computeDistanceMatrix(codes, opts, sample).matrix, n, ids);
      clampedLimbs += rep.clampedLimbs;
      for (const key of collectSplits(rep.root)) tally.set(key, (tally.get(key) || 0) + 1);
      report({
        stage: "bootstrap",
        pct: 24 + Math.round(((r + 1) / reps) * 64),
        detail: `replicate ${r + 1}/${reps}`,
      });
      if (r % 2 === 0) await tick();
    }
  }
  if (reps > 0) annotateSupports(tree, tally, reps);
  else clearSupports(tree);

  return {
    columnsUsed: distRes.columnsUsed,
    matrix: distRes.matrix,
    transitions: distRes.transitions,
    transversions: distRes.transversions,
    comparable: distRes.comparable,
    tree,
    clampedLimbs,
    reps,
    summary: {
      meanDistance: distRes.meanDistance,
      closestPair: distRes.closestPair,
      furthestPair: distRes.furthestPair,
      tiTvRatio: distRes.tiTvRatio,
      saturatedPairs: distRes.saturatedPairs,
      noDataPairs: distRes.noDataPairs,
    },
  };
}

/* ---------------------------------------------------------------- WASM path */

/*
 * Blob layout (little-endian), written by gx_phylo_finish():
 *   "GXPR" u32 version u32 n u32 L u32 columnsUsed i32 numNodes(preorder)
 *   i32 parent[numNodes] i32 nameIdx[numNodes] u8 childCount[numNodes] pad8
 *   f64 branchLen[numNodes*3]   (per-node child slots; parent-relative)
 *   f64 meanDistance f64 tiTvRatio f64 clampedLimbsMain
 *   i32 closestPair[2] i32 furthestPair[2]
 *   u32 saturatedPairs u32 noDataPairs u32 clampedLimbsReps
 *   f64 d[n*n] f64 tiM[n*n] f64 tvM[n*n] f64 cmpM[n*n]
 *   u32 nsplits u32 maskBytes u32 counts[nsplits] u8 masks[nsplits*maskBytes]
 */
function decodePhyloBlob(buf, ids) {
  const r = new BlobReader(buf);
  r.magic("GXPR");
  r.u32(); // version
  const n = r.u32();
  const L = r.u32();
  const columnsUsed = r.u32();
  const nn = r.i32();

  for (let i = 0; i < nn; i++) r.i32(); // parent[] (implicit in preorder; kept for debuggers)
  const nameIdx = new Int32Array(nn);
  for (let i = 0; i < nn; i++) nameIdx[i] = r.i32();
  const childCount = new Uint8Array(nn);
  for (let i = 0; i < nn; i++) childCount[i] = r.readU8();
  r.align(8);
  const blens = new Float64Array(nn * 3);
  for (let i = 0; i < nn * 3; i++) blens[i] = r.f64();

  const meanDistance = r.f64();
  const tiTvRatioRaw = r.f64();
  const clampedLimbsMain = r.f64();
  const closestPair = [r.i32(), r.i32()];
  const furthestPair = [r.i32(), r.i32()];
  const saturatedPairs = r.u32();
  const noDataPairs = r.u32();
  const clampedLimbsReps = r.u32();

  r.align(8); /* writer pads to keep matrices f64-aligned */
  const f64 = (count) => {
    const off = r.pos - buf.byteOffset;
    r.skip(8 * count);
    return new Float64Array(buf.buffer, off, count);
  };
  const d = f64(n * n);
  const tiM = f64(n * n);
  const tvM = f64(n * n);
  const cmpM = f64(n * n);

  const nsplits = r.u32();
  const maskBytes = r.u32();
  const counts = new Uint32Array(nsplits);
  for (let i = 0; i < nsplits; i++) counts[i] = r.u32();
  const masks = r.bytes(nsplits * maskBytes).slice();

  // rebuild the preorder-serialized tree, preserving child order exactly
  let cursor = 0;
  function readNode(parentBl, slot) {
    const pos = cursor++;
    const isLeaf = nameIdx[pos] >= 0;
    const node = {
      isLeaf,
      name: isLeaf ? ids[nameIdx[pos]] : null,
      branchLen: parentBl ? parentBl[slot] : 0,
      children: [],
    };
    const cc = childCount[pos];
    const myBl = [blens[pos * 3], blens[pos * 3 + 1], blens[pos * 3 + 2]];
    for (let c = 0; c < cc; c++) node.children.push(readNode(myBl, c));
    return node;
  }

  return {
    n, L, columnsUsed, meanDistance, tiTvRatioRaw, clampedLimbsMain,
    closestPair, furthestPair, saturatedPairs, noDataPairs, clampedLimbsReps,
    d, tiM, tvM, cmpM, nsplits, maskBytes, counts, masks, readNode,
  };
}

async function runPhyloWasm(records, params, report, opts, rooting, ids, n, L, tick) {
  const mod = await gxWasm();

  const encoded = encodeRows(records.map((r) => r.seq));
  const flat = new Uint8Array(n * L);
  encoded.forEach((row, i) => flat.set(row, i * L));

  report({ stage: "encode", pct: 4, detail: "Encoding alignment" });
  report({
    stage: "distance", pct: 9,
    detail: `${((n * (n - 1)) / 2)} pairwise distances · ${MODELS[opts.model].short} · ${opts.gapMode} deletion`,
  });

  const reps = Math.max(0, Math.min(1000, Math.round(params.bootstrap ?? 100)));
  const seed = Number.isFinite(+params.seed) ? +params.seed : 42;

  report({ stage: "tree", pct: 20, detail: opts.method === "upgma" ? "UPGMA clustering" : "Neighbor-Joining search" });

  const fp = mallocIntoHeap(mod, flat);
  let outBlob;
  try {
    mod._gx_phylo_init(fp, n, L, MODEL_IDS[opts.model], opts.gapMode === "complete" ? 1 : 0,
                       opts.method === "upgma" ? 1 : 0, reps, seed);

    if (reps > 0) {
      while (mod._gx_phylo_reps_done() < reps) {
        const did = mod._gx_phylo_bootstrap_chunk(Math.min(4, reps - mod._gx_phylo_reps_done()));
        const done = mod._gx_phylo_reps_done();
        report({
          stage: "bootstrap",
          pct: 24 + Math.round((done / reps) * 64),
          detail: `replicate ${done}/${reps}`,
        });
        await tick();
        if (did === 0) break;
      }
    }

    mod._gx_phylo_finish();
    outBlob = takeBlob(mod, mod._gx_phylo_result_ptr(), mod._gx_phylo_result_len());
  } finally {
    mod._gx_phylo_abort();
    mod._free(fp);
  }

  const dec = decodePhyloBlob(outBlob, ids);

  const tree = dec.readNode(null, 0);
  if (opts.method === "nj") tree.unrooted = true;

  // split tally keyed by canonical comma-list strings (smaller side), exactly
  // what splits.js canonKeyOf produces for the same taxon set
  const tally = new Map();
  for (let s = 0; s < dec.nsplits; s++) {
    const members = [];
    for (let t = 0; t < n; t++) {
      if ((dec.masks[s * dec.maskBytes + (t >> 3)] >> (t & 7)) & 1) members.push(t);
    }
    tally.set(members.join(","), dec.counts[s]);
  }

  let finalTree = tree;
  if (opts.method === "nj" && rooting === "midpoint") finalTree = midpointRoot(tree);
  indexAndAssignMembers(finalTree);

  if (reps > 0) annotateSupports(finalTree, tally, reps);
  else clearSupports(finalTree);

  return {
    columnsUsed: dec.columnsUsed,
    matrix: dec.d,
    transitions: dec.tiM,
    transversions: dec.tvM,
    comparable: dec.cmpM,
    tree: finalTree,
    clampedLimbs: dec.clampedLimbsMain + dec.clampedLimbsReps,
    reps,
    summary: {
      meanDistance: dec.meanDistance,
      closestPair: dec.closestPair[0] >= 0 ? dec.closestPair : null,
      furthestPair: dec.furthestPair[0] >= 0 ? dec.furthestPair : null,
      tiTvRatio: Number.isFinite(dec.tiTvRatioRaw) ? dec.tiTvRatioRaw : null,
      saturatedPairs: dec.saturatedPairs,
      noDataPairs: dec.noDataPairs,
    },
  };
}

/* ------------------------------------------------------------------ public */

/**
 * Returns the exact fields phylo.worker.js posts in its "done" message
 * (minus type/ms, which stay in the worker).
 */
export async function runPhylo(records, params, report) {
  report = report || (() => {});
  const tick = () => new Promise((r) => setTimeout(r, 0));
  const { ids, n, L } = validateRecords(records);
  const { method, opts, rooting } = resolveOpts(params);
  const t0 = performance.now();

  let core;
  let engineUsed = "wasm";
  try {
    core = await runPhyloWasm(records, params, report, opts, rooting, ids, n, L, tick);
  } catch (err) {
    console.warn("[phylo] WASM engine unavailable, using JS fallback:", err?.message || err);
    core = await runPhyloJs(records, params, report, opts, rooting, ids, n);
    engineUsed = "js";
  }

  const newick = toNewick(core.tree, { includeSupports: core.reps > 0 });

  return {
    ids,
    numTaxa: n,
    columns: L,
    columnsUsed: core.columnsUsed,
    matrix: Array.from(core.matrix),
    transitions: Array.from(core.transitions),
    transversions: Array.from(core.transversions),
    comparable: Array.from(core.comparable),
    tree: core.tree,
    unrooted: method === "nj",
    newick,
    bootstrapReps: core.reps,
    summary: { ...core.summary, clampedLimbs: core.clampedLimbs },
    paramsEcho: { method, model: opts.model, gapMode: opts.gapMode, rooting },
    ms: performance.now() - t0,
    engine: engineUsed,
  };
}
