/**
 * Progressive multiple sequence alignment - the ClustalW architecture:
 *
 *   1. all-against-all k-mer distance estimates (fast, approximate)
 *   2. UPGMA guide tree over those distances
 *   3. depth-first merge: sequences pair up at the leaves, then profiles
 *      align against profiles via global DP with affine gaps until the
 *      root profile IS the multiple alignment
 *
 * Deterministic throughout: tie-breaking in clustering and traceback is
 * index-stable, so the same input always yields the same alignment.
 */
import { encodeSeq } from "./alphabet.js";
import { kmerDistanceMatrix } from "./kmer.js";
import { upgma } from "./guidetree.js";
import { singleProfile, mergeProfiles } from "./profile.js";

export function progressiveAlign(records, params = {}, report = () => {}) {
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

  // Depth-first progressive merge along the guide tree.
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

  // Members carry origIdx, so reassemble rows in input order.
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
