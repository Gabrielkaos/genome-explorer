import { reverseComplement, canonicalKmerHash } from "./sequence.js";

/**
 * Compute a minimizer sketch for one sequence: for every window of `w`
 * consecutive k-mers, keep the one with the smallest canonical hash. This
 * is the standard trick (Roberts et al. 2004) that lets long-read overlap
 * tools like minimap2/miniasm avoid comparing every k-mer of every read
 * against every other read - only a sparse, deterministic sample of each
 * read's k-mers needs to be indexed and compared.
 *
 * Returns an array of { hash, pos, strand } - `pos` is the k-mer's start
 * position in the read's forward-strand coordinates, `strand` is +1 if the
 * k-mer's own sequence was the canonical (smaller-hash) orientation, -1 if
 * its reverse complement was. Consecutive duplicate minimizers are
 * collapsed (only the first occurrence within a run is kept).
 */
export function computeMinimizers(seq, { k = 10, w = 5 } = {}) {
  const n = seq.length;
  if (n < k) return [];

  // Precompute canonical hash + strand for every k-mer position.
  // (Straightforward string-slicing implementation, not bit-packed rolling
  // hashes - simpler to verify correct, and fast enough at the read
  // lengths / read counts this in-browser tool targets.)
  const numKmers = n - k + 1;
  const hashes = new Array(numKmers);
  const strands = new Int8Array(numKmers);
  for (let i = 0; i < numKmers; i++) {
    const km = seq.substring(i, i + k);
    if (/[^ACGT]/i.test(km)) { hashes[i] = null; continue; } // skip ambiguous bases
    const rc = reverseComplement(km);
    const { hash, strand } = canonicalKmerHash(km, rc);
    hashes[i] = hash;
    strands[i] = strand;
  }

  const minimizers = [];
  let lastMinPos = -1;
  for (let i = 0; i + w <= numKmers; i++) {
    let minIdx = -1, minHash = Infinity;
    for (let j = i; j < i + w; j++) {
      if (hashes[j] === null) continue;
      if (hashes[j] < minHash) { minHash = hashes[j]; minIdx = j; }
    }
    if (minIdx !== -1 && minIdx !== lastMinPos) {
      minimizers.push({ hash: minHash, pos: minIdx, strand: strands[minIdx] });
      lastMinPos = minIdx;
    }
  }
  return minimizers;
}
