import { computeMinimizers } from "./minimizers.js";

/**
 * Index every read's minimizers into hash -> [{readIdx,pos,strand}, ...].
 * Minimizers shared by more than `maxOccurrence` distinct reads are dropped
 * (repeat-region filtering - the same reason minimap2/miniasm do this:
 * without it, one repetitive element creates a combinatorial explosion of
 * spurious "overlaps" between otherwise-unrelated reads). Minimizers unique
 * to a single read are dropped too since they can never support an overlap.
 */
export function buildMinimizerIndex(reads, { k = 10, w = 5, maxOccurrence = 40 } = {}) {
  const index = new Map();
  reads.forEach((r, readIdx) => {
    const mins = computeMinimizers(r.seq, { k, w });
    for (const m of mins) {
      let bucket = index.get(m.hash);
      if (!bucket) { bucket = []; index.set(m.hash, bucket); }
      bucket.push({ readIdx, pos: m.pos, strand: m.strand });
    }
  });
  for (const [hash, bucket] of index) {
    const distinct = new Set(bucket.map((b) => b.readIdx));
    if (distinct.size > maxOccurrence || distinct.size < 2) index.delete(hash);
  }
  return index;
}

/** Tally shared-minimizer matches between every pair of reads that share at least one. */
export function findCandidatePairs(index) {
  const pairMatches = new Map(); // "i,j" (i<j) -> [{posA,posB,strandA,strandB}]
  for (const bucket of index.values()) {
    for (let a = 0; a < bucket.length; a++) {
      for (let b = a + 1; b < bucket.length; b++) {
        const X = bucket[a], Y = bucket[b];
        if (X.readIdx === Y.readIdx) continue;
        const first = X.readIdx < Y.readIdx ? X : Y;
        const second = X.readIdx < Y.readIdx ? Y : X;
        const key = `${first.readIdx},${second.readIdx}`;
        let arr = pairMatches.get(key);
        if (!arr) { arr = []; pairMatches.set(key, arr); }
        arr.push({ posA: first.pos, posB: second.pos, strandA: first.strand, strandB: second.strand });
      }
    }
  }
  return pairMatches;
}

/**
 * Given the shared-minimizer matches between read A and read B, estimate
 * whether they genuinely overlap: which orientation, where, how long, and
 * with what confidence.
 *
 * Method: a k-mer match where both reads picked the SAME strand as
 * canonical means A and B are in the same orientation at that point; a
 * match where they picked OPPOSITE strands means B needs to be reverse
 * complemented to align with A. We convert every match into a unified
 * "diagonal" coordinate (position-in-A minus position-in-B, after
 * transforming B's coordinate into whichever orientation the match
 * implies), then take the most heavily supported diagonal as the overlap.
 * This is a simplified stand-in for full collinear chaining (as used in
 * minimap2) - fast, but approximate; it does not verify the overlap with a
 * base-level alignment, so `score` is a density-based confidence estimate,
 * not a true percent identity.
 */
export function estimateOverlap(matches, lenA, lenB, { k = 10, w = 5, diagTolerance = 60, minMatches = 4, minOverlapLen = 200 } = {}) {
  if (matches.length < minMatches) return null;

  // Convert every match into a unified diagonal coordinate per orientation
  // hypothesis, then chain SEQUENTIALLY (sorted by posA), allowing the
  // diagonal to drift by up to `diagTolerance` between *consecutive*
  // matches rather than requiring the whole overlap to sit in one fixed
  // bin. This matters a lot for real noisy long reads: each indel error
  // shifts the diagonal by about 1bp, and those small shifts accumulate
  // over a long overlap, so a rigid global bin would fragment one true
  // overlap into several under-sized pieces.
  function chainBestChain(list) {
    if (list.length < minMatches) return null;
    const sorted = [...list].sort((a, b) => a.posA - b.posA);
    let best = null, chain = [sorted[0]];
    for (let i = 1; i < sorted.length; i++) {
      const prevDiag = chain[chain.length - 1].posA - chain[chain.length - 1].posBt;
      const curDiag = sorted[i].posA - sorted[i].posBt;
      if (Math.abs(curDiag - prevDiag) <= diagTolerance) {
        chain.push(sorted[i]);
      } else {
        if (!best || chain.length > best.length) best = chain;
        chain = [sorted[i]];
      }
    }
    if (!best || chain.length > best.length) best = chain;
    return best.length >= minMatches ? best : null;
  }

  const same = [], rev = [];
  for (const m of matches) {
    const sameOrientation = m.strandA === m.strandB;
    const posBt = sameOrientation ? m.posB : (lenB - m.posB - k);
    (sameOrientation ? same : rev).push({ posA: m.posA, posBt });
  }

  let best = null;
  for (const [orientation, list] of [["same", same], ["rev", rev]]) {
    const chain = chainBestChain(list);
    if (!chain) continue;
    const aVals = chain.map((x) => x.posA), bVals = chain.map((x) => x.posBt);
    const aMin = Math.min(...aVals), aMax = Math.max(...aVals) + k;
    const bMin = Math.min(...bVals), bMax = Math.max(...bVals) + k;
    const rawSpan = Math.max(aMax - aMin, bMax - bMin);
    const overlapLen = rawSpan + 2 * w;
    if (overlapLen < minOverlapLen) continue;
    const expectedMinimizers = rawSpan / w;
    const score = Math.min(1, chain.length / Math.max(1, expectedMinimizers));
    if (!best || chain.length > best.matchCount) {
      best = { orientation, matchCount: chain.length, overlapLen, score, aMin, aMax, bMin, bMax };
    }
  }
  if (!best) return null;

  // Classify the overlap geometry: does A extend past B, B past A, or is
  // one essentially contained within the other? Slack accounts for the
  // same edge-sampling / error-dropout effect as above - a fixed, moderate
  // margin (not scaled to read length, which would misclassify large
  // partial overlaps on short reads as full containment).
  const slack = Math.max(2 * w + 15, 60);
  const aStartsAt0 = best.aMin <= slack, aEndsAtLen = best.aMax >= lenA - slack;
  const bStartsAt0 = best.bMin <= slack, bEndsAtLen = best.bMax >= lenB - slack;

  let type;
  if (aStartsAt0 && aEndsAtLen) type = "A_contained_in_B";
  else if (bStartsAt0 && bEndsAtLen) type = "B_contained_in_A";
  else if (aEndsAtLen && bStartsAt0) type = "A_to_B"; // A's suffix overlaps B's prefix (in B's effective orientation)
  else if (bEndsAtLen && aStartsAt0) type = "B_to_A"; // B's suffix overlaps A's prefix
  else return null; // internal/spurious match, not a valid dovetail or containment overlap

  return {
    orientation: best.orientation, // "same" | "rev" (B reverse-complemented relative to A)
    type,
    matchCount: best.matchCount,
    overlapLen: best.overlapLen,
    score: best.score,
    aRange: [best.aMin, best.aMax],
    bRange: [best.bMin, best.bMax], // in transformed ("B as if same-orientation as A") coordinates
  };
}

/** Full pipeline: reads -> accepted overlap edge list. */
export function detectOverlaps(reads, params = {}, onProgress) {
  const index = buildMinimizerIndex(reads, params);
  const pairs = findCandidatePairs(index);
  const overlaps = [];
  let processed = 0;
  const total = pairs.size;
  for (const [key, matches] of pairs) {
    const [i, j] = key.split(",").map(Number);
    const ov = estimateOverlap(matches, reads[i].seq.length, reads[j].seq.length, params);
    if (ov) overlaps.push({ a: i, b: j, ...ov });
    processed++;
    if (onProgress && processed % 500 === 0) onProgress(processed, total);
  }
  return overlaps;
}
