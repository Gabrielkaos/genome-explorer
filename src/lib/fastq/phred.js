/**
 * Phred quality score utilities.
 *
 * Assumption: Phred+33 encoding (Sanger / Illumina 1.8+ / all modern Nanopore
 * basecallers including Guppy & Dorado). Legacy Phred+64 is not supported.
 *
 * IMPORTANT correctness note: quality scores are on a *logarithmic* scale, so
 * naively averaging Q values (sum/n) is mathematically wrong and biases the
 * result. The correct approach (used by NanoPlot, FastQC, etc.) is:
 *   1. Convert each Q to a linear error probability: p = 10^(-Q/10)
 *   2. Average the probabilities
 *   3. Convert the average probability back to Q: Q = -10 * log10(meanP)
 * This module implements that method throughout.
 */

export const PHRED_OFFSET = 33;

/** Single quality character -> Phred Q score. */
export function qualCharToScore(charCode) {
  return charCode - PHRED_OFFSET;
}

/** Phred Q score -> linear error probability. */
export function scoreToErrorProb(q) {
  return Math.pow(10, -q / 10);
}

/** Linear error probability -> Phred Q score. */
export function errorProbToScore(p) {
  return -10 * Math.log10(Math.max(p, 1e-12));
}

/**
 * Compute the statistically-correct mean quality for a run of raw quality
 * bytes (Uint8Array of ASCII char codes, Phred+33 encoded).
 * Returns { meanQ, sumErrorProb, count }.
 */
export function meanQualityFromBytes(qualBytes, start = 0, end = qualBytes.length) {
  let sumP = 0;
  const n = end - start;
  for (let i = start; i < end; i++) {
    const q = qualBytes[i] - PHRED_OFFSET;
    sumP += scoreToErrorProb(q);
  }
  const meanP = n > 0 ? sumP / n : 0;
  return { meanQ: errorProbToScore(meanP), sumErrorProb: sumP, count: n };
}

/** Same as above but for a JS string of quality characters. */
export function meanQualityFromString(qualStr) {
  let sumP = 0;
  for (let i = 0; i < qualStr.length; i++) {
    const q = qualStr.charCodeAt(i) - PHRED_OFFSET;
    sumP += scoreToErrorProb(q);
  }
  const n = qualStr.length;
  const meanP = n > 0 ? sumP / n : 0;
  return { meanQ: errorProbToScore(meanP), sumErrorProb: sumP, count: n };
}

/**
 * Dataset-level mean quality: base-weighted average across every base in
 * every read (NOT the average of per-read means, which would over-weight
 * short reads). Feed it running totals accumulated during parsing.
 */
export function datasetMeanQuality(totalErrorProbSum, totalBases) {
  if (totalBases === 0) return 0;
  return errorProbToScore(totalErrorProbSum / totalBases);
}
