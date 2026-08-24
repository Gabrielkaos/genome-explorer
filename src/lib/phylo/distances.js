/**
 * Pairwise evolutionary distance estimation from an ALIGNED matrix
 * (equal-length rows). Three standard models:
 *
 *   p    - uncorrected proportion of differing sites
 *   jc   - Jukes-Cantor (1969) single-parameter correction:
 *          d = -3/4 · ln(1 - 4p/3), which corrects for multiple hits
 *          and back-mutations assuming equal base frequencies and rates.
 *   k2p  - Kimura (1980) two-parameter correction that counts
 *          transitions (A<->G, C<->T) separately from transversions:
 *          d = -1/2·ln(1-2P-Q) - 1/4·ln(1-2Q)
 *
 * Gaps / ambiguous codes are excluded pairwise (only columns where BOTH
 * sequences carry a resolved ACGT are compared) or completely (columns
 * with any unresolved character in ANY sequence are dropped first).
 *
 * An optional `colIdx` (e.g. a bootstrap resample of column positions)
 * restricts the computation to those columns without materializing
 * resampled copies of the alignment.
 */

export const MODELS = {
  p: { id: "p", label: "Uncorrected p-distance", short: "p-dist" },
  jc: { id: "jc", label: "Jukes-Cantor (JC69)", short: "JC69" },
  k2p: { id: "k2p", label: "Kimura 2-parameter (K2P)", short: "K2P" },
};

export const GAP_MODES = {
  pairwise: { id: "pairwise", label: "Pairwise deletion" },
  complete: { id: "complete", label: "Complete deletion" },
};

/** Pairs whose corrected distance exceeds this cap are flagged saturated. */
export const SATURATION_CAP = 1.5; // substitutions/site

const CODE = { A: 0, C: 1, G: 2, T: 3 };
const TRANSITION = [
  //  A      C      G      T     (from \ to)
  [false, false, true, false],
  [false, false, false, true],
  [true, false, false, false],
  [false, true, false, false],
];

/** Normalize one aligned row into base codes: 0-3 = ACGT, 4 = gap/ambiguity. */
export function encodeRow(seq) {
  const out = new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) {
    let ch = seq[i];
    if (ch === "U") ch = "T"; // tolerate RNA input
    const c = CODE[ch];
    out[i] = c === undefined ? 4 : c;
  }
  return out;
}

export function encodeRows(rows) {
  return rows.map(encodeRow);
}

/**
 * rows: array of equal-length strings (aligned).
 * Returns flat n*n matrices (row-major) plus summary flags.
 */
export function computeDistanceMatrix(rowsOrCodes, { model = "jc", gapMode = "pairwise" } = {}, colIdx = null) {
  const codes = Array.isArray(rowsOrCodes) && typeof rowsOrCodes[0] === "string"
    ? encodeRows(rowsOrCodes)
    : rowsOrCodes;
  const n = codes.length;
  if (!n) throw new Error("No sequences supplied.");
  const L = codes[0].length;

  const d = new Float64Array(n * n);
  const tiM = new Float64Array(n * n);
  const tvM = new Float64Array(n * n);
  const cmpM = new Float64Array(n * n);

  // Effective column universe.
  let cols;
  if (colIdx) {
    cols = colIdx;
  } else {
    cols = new Uint32Array(L);
    for (let c = 0; c < L; c++) cols[c] = c;
  }

  // Complete deletion: drop every column where any sequence is unresolved.
  if (gapMode === "complete") {
    const keep = [];
    for (let k = 0; k < cols.length; k++) {
      const c = cols[k];
      let ok = true;
      for (let r = 0; r < n; r++) {
        if (codes[r][c] > 3) { ok = false; break; }
      }
      if (ok) keep.push(c);
    }
    cols = Uint32Array.from(keep);
  }

  let saturatedPairs = 0, noDataPairs = 0;
  const saturatedList = [];

  for (let i = 0; i < n; i++) {
    const ri = codes[i];
    d[i * n + i] = 0;
    for (let j = i + 1; j < n; j++) {
      const rj = codes[j];
      let diffs = 0, ti = 0, tv = 0, comp = 0;
      for (let k = 0; k < cols.length; k++) {
        const a = ri[cols[k]], b = rj[cols[k]];
        if (a > 3 || b > 3) continue;
        comp++;
        if (a !== b) {
          diffs++;
          if (TRANSITION[a][b]) ti++;
          else tv++;
        }
      }
      const raw = comp ? diffs / comp : 0;
      let dist = raw;
      if (model === "jc") {
        if (raw > 0 && raw < 0.75) dist = -0.75 * Math.log(1 - (4 * raw) / 3);
        else if (raw >= 0.75) dist = Infinity;
      } else if (model === "k2p") {
        const P = comp ? ti / comp : 0;
        const Q = comp ? tv / comp : 0;
        const s = 1 - 2 * P - Q;
        const v = 1 - 2 * Q;
        if ((s <= 0 || v <= 0) && (P > 0 || Q > 0)) dist = Infinity;
        else dist = -0.5 * Math.log(Math.max(s, Number.MIN_VALUE)) - 0.25 * Math.log(Math.max(v, Number.MIN_VALUE));
      }
      /* quantize corrected distances to 1e-9 so JS and WASM kernels agree
       * despite libm last-bit differences in log() */
      if (Number.isFinite(dist) && dist !== raw) dist = Math.round(dist * 1e9) / 1e9;
      if (!Number.isFinite(dist)) {
        saturatedPairs++;
        saturatedList.push([i, j]);
        dist = SATURATION_CAP;
      }
      if (comp === 0) noDataPairs++;
      d[i * n + j] = d[j * n + i] = dist;
      tiM[i * n + j] = tiM[j * n + i] = ti;
      tvM[i * n + j] = tvM[j * n + i] = tv;
      cmpM[i * n + j] = cmpM[j * n + i] = comp;
    }
  }

  let tiSum = 0, tvSum = 0, sum = 0, count = 0;
  let minD = Infinity, maxD = 0, minPair = null, maxPair = null;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const v = d[i * n + j];
      tiSum += tiM[i * n + j];
      tvSum += tvM[i * n + j];
      sum += v; count++;
      if (v < minD) { minD = v; minPair = [i, j]; }
      if (v > maxD) { maxD = v; maxPair = [i, j]; }
    }
  }

  return {
    n,
    columnsUsed: cols.length,
    matrix: d,
    transitions: tiM,
    transversions: tvM,
    comparable: cmpM,
    saturatedPairs,
    saturatedList,
    noDataPairs,
    meanDistance: count ? sum / count : 0,
    closestPair: minPair,
    furthestPair: maxPair,
    tiTvRatio: tvSum > 0 ? tiSum / tvSum : null,
  };
}
