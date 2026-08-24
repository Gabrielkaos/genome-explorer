import { createHistogramAccumulator } from "./stats.js";
import { scoreToErrorProb, errorProbToScore } from "./phred.js";

/**
 * "Yield by length" curve (the classic NanoPlot chart used to visually
 * locate N50): for each length threshold L, how many total bases come from
 * reads at least L long. Downsampled to a fixed number of plot points so
 * this stays fast even for million-read runs.
 */
export function computeYieldByLength(lengths, n50, maxPoints = 400) {
  const sorted = Array.from(lengths).sort((a, b) => b - a); // longest first
  const n = sorted.length;
  if (n === 0) return { points: [], n50CumBases: 0 };

  const points = [];
  let running = 0;
  const step = Math.max(1, Math.floor(n / maxPoints));
  let n50CumBases = 0;
  let n50Recorded = false;

  for (let i = 0; i < n; i++) {
    running += sorted[i];
    if (!n50Recorded && sorted[i] <= n50) { n50CumBases = running; n50Recorded = true; }
    if (i % step === 0 || i === n - 1) {
      points.push({ length: sorted[i], cumulativeBases: running });
    }
  }
  if (!n50Recorded) n50CumBases = running;
  return { points, n50CumBases, totalBases: running };
}

/** Per-read GC% histogram, 0-100 in bins of 5. */
export function computeGCHistogram(gcPerRead) {
  const edges = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 101];
  const hist = createHistogramAccumulator(edges);
  for (let i = 0; i < gcPerRead.length; i++) hist.add(gcPerRead[i]);
  return hist.toArray();
}

/**
 * Mean quality by position-in-read, averaged (correctly, via error
 * probability) across a sample of reads. Only previewRecords have full
 * per-base quality text available in memory (see worker docs), so this is
 * explicitly a sample-based statistic, not a whole-dataset one - the UI
 * must say so.
 */
export function computePositionalQuality(previewRecords, maxPosition = 3000, minReadsPerPosition = 20) {
  if (!previewRecords || previewRecords.length === 0) return [];
  const sumErrByPos = new Float64Array(maxPosition);
  const countByPos = new Int32Array(maxPosition);

  for (const rec of previewRecords) {
    const qual = rec.qual;
    const len = Math.min(qual.length, maxPosition);
    for (let p = 0; p < len; p++) {
      const q = qual.charCodeAt(p) - 33;
      sumErrByPos[p] += scoreToErrorProb(q);
      countByPos[p] += 1;
    }
  }

  const out = [];
  // stride the output so charts don't render thousands of points for long reads
  const strideCandidatePositions = [];
  for (let p = 0; p < maxPosition; p++) {
    if (countByPos[p] >= minReadsPerPosition) strideCandidatePositions.push(p);
  }
  const stride = Math.max(1, Math.floor(strideCandidatePositions.length / 300));
  for (let i = 0; i < strideCandidatePositions.length; i += stride) {
    const p = strideCandidatePositions[i];
    const meanQ = errorProbToScore(sumErrByPos[p] / countByPos[p]);
    out.push({ position: p, meanQ, readsAtPosition: countByPos[p] });
  }
  return out;
}
