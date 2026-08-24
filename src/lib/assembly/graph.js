/**
 * Turns a list of pairwise overlaps into contigs.
 *
 * Model: every read has two ends, 'L' (5', its own native start) and 'R'
 * (3', its own native end). A dovetail overlap connects exactly one end of
 * read A to exactly one end of read B, optionally with a "flip" (read B
 * must be used in reverse-complement relative to A for this connection to
 * make sense). Since each end can only truly be glued to one other end in
 * a correct layout, we resolve conflicts with mutual-best-match: an edge
 * is confirmed only if it is the highest-scoring candidate for BOTH ends
 * it touches. This is a simplification of full string-graph transitive
 * reduction (as used in real assemblers like miniasm) - simpler to reason
 * about and implement correctly, at the cost of being more conservative
 * (it can leave some real overlaps unused rather than risk a wrong join).
 *
 * Once conflicts are resolved, every end has degree <= 1, so the whole
 * graph is trivially a disjoint union of simple paths and simple cycles -
 * paths become linear contigs, cycles become circular contigs (e.g. a
 * plasmid).
 */

function edgeForOverlap(ov) {
  // Returns { endA, endB, flip } - which end of A connects to which end
  // of B, and whether B needs to be used in reverse-complement.
  const flip = ov.orientation === "rev";
  if (ov.type === "A_to_B") return { endA: "R", endB: flip ? "R" : "L", flip };
  if (ov.type === "B_to_A") return { endA: "L", endB: flip ? "L" : "R", flip };
  return null; // containment overlaps don't create dovetail edges
}

export function buildGraph(reads, overlaps) {
  const contained = new Set();
  for (const ov of overlaps) {
    if (ov.type === "A_contained_in_B") contained.add(ov.a);
    else if (ov.type === "B_contained_in_A") contained.add(ov.b);
  }

  const dovetail = overlaps.filter(
    (ov) => (ov.type === "A_to_B" || ov.type === "B_to_A") && !contained.has(ov.a) && !contained.has(ov.b)
  );

  // candidates per (readIdx:end) node
  const candidates = new Map(); // key -> [{otherKey, flip, score, ov}]
  function addCandidate(key, otherKey, flip, score, ov) {
    let arr = candidates.get(key);
    if (!arr) { arr = []; candidates.set(key, arr); }
    arr.push({ otherKey, flip, score, ov });
  }
  const edgeList = []; // {keyA, keyB, flip, score, ov}
  for (const ov of dovetail) {
    const e = edgeForOverlap(ov);
    if (!e) continue;
    const keyA = `${ov.a}:${e.endA}`, keyB = `${ov.b}:${e.endB}`;
    addCandidate(keyA, keyB, e.flip, ov.matchCount, ov);
    addCandidate(keyB, keyA, e.flip, ov.matchCount, ov);
    edgeList.push({ keyA, keyB, flip: e.flip, score: ov.matchCount, ov });
  }

  // Greedy maximum-weight matching: each (read,end) node can carry at most
  // one confirmed edge (that's what keeps the graph a simple union of paths
  // and cycles - see extractPaths). Strict "reciprocal best choice only"
  // was too conservative under real coverage depth: with many reads having
  // similar-strength overlaps with several neighbors, a read's single best
  // partner often doesn't return the favor, so most real, valid overlaps
  // never got confirmed. Sorting all candidate edges by strength and
  // greedily taking each one whose *both* ends are still free is the
  // standard greedy approximation to maximum-weight matching, and uses
  // far more of the real signal already computed.
  edgeList.sort((a, b) => b.score - a.score);
  const usedEnds = new Set();
  const confirmed = new Map(); // key -> {otherKey, flip, ov}
  for (const e of edgeList) {
    if (usedEnds.has(e.keyA) || usedEnds.has(e.keyB)) continue;
    usedEnds.add(e.keyA); usedEnds.add(e.keyB);
    confirmed.set(e.keyA, { otherKey: e.keyB, flip: e.flip, ov: e.ov });
    confirmed.set(e.keyB, { otherKey: e.keyA, flip: e.flip, ov: e.ov });
  }

  return { contained, confirmed, dovetail };
}

function otherEnd(end) { return end === "L" ? "R" : "L"; }

/** Walk the simplified end-graph into an ordered list of contigs (paths and cycles). */
export function extractPaths(reads, graph) {
  const { contained, confirmed } = graph;
  const visitedEnds = new Set();
  const paths = [];

  function degree(readIdx) {
    let d = 0;
    if (confirmed.has(`${readIdx}:L`)) d++;
    if (confirmed.has(`${readIdx}:R`)) d++;
    return d;
  }

  const allReadIdx = reads.map((_, i) => i).filter((i) => !contained.has(i));

  // Pass 1: linear chains, starting from reads with degree 0 or 1 (a true endpoint).
  for (const startIdx of allReadIdx) {
    if (visitedEnds.has(`${startIdx}:L`) || visitedEnds.has(`${startIdx}:R`)) continue;
    if (degree(startIdx) === 2) continue; // internal to some chain - visit from its real endpoint
    walkFrom(startIdx);
  }
  // Pass 2: whatever's left must be pure cycles (every node degree 2, no free ends).
  for (const startIdx of allReadIdx) {
    if (visitedEnds.has(`${startIdx}:L`) || visitedEnds.has(`${startIdx}:R`)) continue;
    walkFrom(startIdx, true);
  }

  function walkFrom(startIdx, isCycle = false) {
    // Free end to start from: prefer the end with no confirmed edge (true
    // chain start); if both are free (isolated read) or both connected
    // (cycle), start at 'L' arbitrarily.
    // Exit via whichever end actually has a confirmed connection (that's
    // the direction the chain continues in); if neither does, this is an
    // isolated single-read contig and the choice is arbitrary.
    const startEnd = confirmed.has(`${startIdx}:R`) ? "R" : confirmed.has(`${startIdx}:L`) ? "L" : "L";

    const steps = [{ readIdx: startIdx, strand: 1, entryEnd: otherEnd(startEnd) }];
    visitedEnds.add(`${startIdx}:L`); visitedEnds.add(`${startIdx}:R`);

    let curIdx = startIdx, curExitEnd = startEnd, curStrand = 1;
    let circular = false;
    let guard = 0;
    while (guard++ < reads.length + 5) {
      const edge = confirmed.get(`${curIdx}:${curExitEnd}`);
      if (!edge) break; // reached a free end - path complete
      const [nextIdxStr, nextEntryEndNative] = edge.otherKey.split(":");
      const nextIdx = Number(nextIdxStr);
      if (nextIdx === startIdx && steps.length > 1) { circular = true; break; } // closed the loop
      const nextStrand = edge.flip ? -curStrand : curStrand;
      // If we're using the next read in reverse strand, its "effective"
      // entry end (in the orientation we'll actually traverse it) is the
      // opposite of its native end label.
      const effectiveEntryEnd = nextStrand === 1 ? nextEntryEndNative : otherEnd(nextEntryEndNative);
      steps.push({ readIdx: nextIdx, strand: nextStrand, entryEnd: effectiveEntryEnd, overlapWithPrev: edge.ov });
      visitedEnds.add(`${nextIdx}:L`); visitedEnds.add(`${nextIdx}:R`);
      curIdx = nextIdx;
      curExitEnd = otherEnd(nextEntryEndNative); // native end labels are strand-independent
      curStrand = nextStrand;
    }
    paths.push({ steps, circular: circular || isCycle, closingOverlap: circular ? confirmed.get(`${curIdx}:${curExitEnd}`)?.ov : null });
  }

  return paths;
}
