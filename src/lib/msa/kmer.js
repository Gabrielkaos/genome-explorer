/**
 * K-mer based pairwise distance estimation, used only to order the guide
 * tree for progressive alignment (the same role the 10-mer/11-mer tuple
 * score plays in ClustalW). It is deliberately approximate: only the tree
 * TOPOLOGY matters downstream, never the absolute distances.
 *
 * Sequences are reduced to dense word-count vectors (2-bit packed k-mers,
 * so k<=8 keeps the vector at <=65536 Int32 slots) and each pair is scored
 * with a Jaccard-style fraction of shared words, then Jukes-Cantor
 * corrected so saturation does not compress divergent pairs together.
 * Pairs with p >= 0.75 (JC's singularity) are capped at distance 3.
 */

export function kmerCounts(encoded, k) {
  const size = 1 << (2 * k);
  const mask = size - 1;
  const vec = new Int32Array(size);
  let code = 0;
  let run = 0;
  let total = 0;
  for (let i = 0; i < encoded.length; i++) {
    const b = encoded[i];
    if (b > 3) { run = 0; continue; }
    code = ((code << 2) | b) & mask;
    if (++run >= k) { vec[code]++; total++; }
  }
  return { vec, total };
}

/** Jukes-Cantor corrected distance from two k-mer count vectors.
 *  Quantized to 1e-9 so the JS and WASM kernels agree despite libm
 *  last-bit differences in log() - the guide tree only needs ordering. */
export function kmerDistance(a, b) {
  let inter = 0;
  let sumA = 0;
  let sumB = 0;
  for (let i = 0; i < a.vec.length; i++) {
    const x = a.vec[i], y = b.vec[i];
    inter += x < y ? x : y;
    sumA += x;
    sumB += y;
  }
  const union = sumA + sumB - inter;
  if (union === 0) return 1;
  const p = 1 - inter / union;
  if (p <= 0) return 0;
  if (p >= 0.75) return 3;
  const d = -0.75 * Math.log(1 - (4 * p) / 3);
  return Math.round(d * 1e9) / 1e9;
}

/**
 * Full symmetric distance matrix in O(n^2). `k` shrinks automatically to
 * fit the shortest sequence; sequences too short for even a 3-mer fall
 * back to prefix identity.
 */
export function kmerDistanceMatrix(encodedSeqs, requestedK = 6) {
  const n = encodedSeqs.length;
  const minLen = Math.min(...encodedSeqs.map((s) => s.length));
  const k = Math.max(3, Math.min(requestedK, minLen - 1));

  const dist = new Float64Array(n * n);
  if (k < 3) {
    // degenerate short-input fallback: crude prefix identity
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const L = Math.min(encodedSeqs[i].length, encodedSeqs[j].length);
        let diff = 0;
        for (let p = 0; p < L; p++) if (encodedSeqs[i][p] !== encodedSeqs[j][p]) diff++;
        const d = L ? diff / L : 1;
        dist[i * n + j] = dist[j * n + i] = d;
      }
    }
    return { dist, k: null };
  }

  const vectors = encodedSeqs.map((s) => kmerCounts(s, k));
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const d = kmerDistance(vectors[i], vectors[j]);
      dist[i * n + j] = dist[j * n + i] = d;
    }
  }
  return { dist, k };
}
