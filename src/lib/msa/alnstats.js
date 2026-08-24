/**
 * Alignment-level statistics computed from the merged root profile's column
 * counts (authoritative composition per column) plus the per-sequence code
 * rows (needed for pairwise identity). Column classification follows the
 * standard phylogenetics vocabulary:
 *
 *   conserved   - every sequence carries the same non-gap character
 *   variable    - at least two distinct non-gap characters present
 *   informative - parsimony-informative: >=2 distinct non-gap characters,
 *                 each seen in >=2 sequences (shared derived states)
 *   singleton   - variable but NOT informative (one sequence's private change)
 */
import { CODE_CHARS } from "./alphabet.js";

export function computeAlnStats({ cols, length, chars }) {
  const nSeq = chars.length;

  let conserved = 0, variable = 0, informative = 0, singleton = 0;
  let gapCols = 0, totalGapCells = 0;
  const consensus = new Uint8Array(length);
  const sites = [];

  for (let c = 0; c < length; c++) {
    const counts = [cols[c * 5], cols[c * 5 + 1], cols[c * 5 + 2], cols[c * 5 + 3]];
    const gapCount = cols[c * 5 + 4];
    totalGapCells += gapCount;

    let bestCode = 4, bestCount = -1, distinct = 0, resTotal = 0;
    const present = [];
    for (let a = 0; a < 4; a++) {
      if (!counts[a]) continue;
      resTotal += counts[a];
      distinct++;
      present.push({ code: a, count: counts[a] });
      if (counts[a] > bestCount) { bestCount = counts[a]; bestCode = a; }
    }
    if (resTotal === 0) { consensus[c] = 4; gapCols++; continue; }

    consensus[c] = bestCount / nSeq >= 0.5 ? bestCode : 4;

    if (distinct === 1 && resTotal === nSeq) conserved++;
    else {
      variable++;
      let infStates = 0;
      for (const p of present) if (p.count >= 2) infStates++;
      if (infStates >= 2) informative++;
      else singleton++;
      sites.push({
        pos: c,
        consensus: CODE_CHARS[bestCode],
        alleles: present.map((p) => `${CODE_CHARS[p.code]}:${p.count}`),
        informative: infStates >= 2,
        indel: gapCount > 0,
      });
    }
  }

  // Pairwise identity over positions where both sequences have a real base.
  const identity = new Float64Array(nSeq * nSeq);
  let sumId = 0, pairs = 0, minId = 1, maxId = 0;
  for (let i = 0; i < nSeq; i++) identity[i * nSeq + i] = 1;
  for (let i = 0; i < nSeq; i++) {
    for (let j = i + 1; j < nSeq; j++) {
      const A = chars[i], B = chars[j];
      let match = 0, compared = 0;
      const L = Math.min(A.length, B.length);
      for (let c = 0; c < L; c++) {
        if (A[c] > 3 || B[c] > 3) continue;
        compared++;
        if (A[c] === B[c]) match++;
      }
      const id = compared ? match / compared : 0;
      identity[i * nSeq + j] = identity[j * nSeq + i] = id;
      sumId += id; pairs++;
      if (id < minId) minId = id;
      if (id > maxId) maxId = id;
    }
  }

  const perSeqGapPct = chars.map((row) => {
    let g = 0;
    for (let c = 0; c < row.length; c++) if (row[c] === 4) g++;
    return row.length ? g / row.length : 0;
  });

  return {
    length,
    numSeqs: nSeq,
    conservedColumns: conserved,
    variableColumns: variable,
    informativeColumns: informative,
    singletonColumns: singleton,
    gapColumns: gapCols,
    gapFraction: length ? totalGapCells / (length * nSeq) : 0,
    meanPairIdentity: pairs ? sumId / pairs : 1,
    minPairIdentity: pairs ? minId : 1,
    maxPairIdentity: pairs ? maxId : 1,
    identityMatrix: identity,
    perSeqGapPct,
    consensusCodes: consensus,
    variableSites: sites,
  };
}
