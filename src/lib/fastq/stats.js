/**
 * Pure statistics helpers used by the FASTQ tool. No parsing or DOM logic
 * lives here so these are easy to unit-test in isolation.
 */

/**
 * Generic "N-statistic" (N50, N90, ...): sort lengths descending, walk the
 * list accumulating total length, and report the length of the entry at
 * which the running total first reaches `fraction` of the overall total.
 * `lengths` may be a plain array or typed array; it is sorted descending
 * in-place on a copy (original is not mutated).
 */
export function computeNStat(lengths, fraction = 0.5) {
  if (!lengths || lengths.length === 0) return 0;
  const sorted = Array.from(lengths).sort((a, b) => b - a);
  const total = sorted.reduce((a, b) => a + b, 0);
  const target = total * fraction;
  let running = 0;
  for (const len of sorted) {
    running += len;
    if (running >= target) return len;
  }
  return sorted[sorted.length - 1] || 0;
}

export function median(values) {
  if (!values || values.length === 0) return 0;
  const sorted = Array.from(values).sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 0 ? (sorted[mid - 1] + sorted[mid]) / 2 : sorted[mid];
}

/**
 * Build log-ish, read-length-appropriate histogram bin edges. Nanopore
 * length distributions are long-tailed, so evenly-spaced bins over the full
 * range tend to waste most bins on the tail. We use fixed, human-legible
 * breakpoints instead, extending automatically if data exceeds them.
 */
export function buildLengthBinEdges(maxLength) {
  const base = [0, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 10000, 15000, 20000, 30000, 50000];
  const edges = base.filter((e) => e <= Math.max(maxLength, 200));
  const last = edges[edges.length - 1];
  if (maxLength > last) edges.push(Math.ceil(maxLength / 1000) * 1000 + 1);
  else edges.push(last + 1);
  return edges;
}

/** Fixed quality bin edges, Q0-Q40+ in steps of 2. */
export function buildQualityBinEdges(maxQ = 40) {
  const edges = [];
  for (let q = 0; q <= Math.max(maxQ, 30); q += 2) edges.push(q);
  edges.push(edges[edges.length - 1] + 2);
  return edges;
}

/**
 * Streaming-friendly histogram accumulator: given fixed bin edges, returns
 * an object with an `add(value)` method and a `toArray()` reader. Designed
 * so the parser worker can call `.add()` once per read without holding
 * every raw value in memory.
 */
export function createHistogramAccumulator(edges) {
  const counts = new Array(edges.length - 1).fill(0);
  function add(value) {
    // binary search for the bin (edges are sorted ascending)
    let lo = 0, hi = counts.length - 1;
    if (value < edges[0]) { counts[0]++; return; }
    if (value >= edges[edges.length - 1]) { counts[counts.length - 1]++; return; }
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (edges[mid] <= value) lo = mid; else hi = mid - 1;
    }
    counts[lo]++;
  }
  function toArray() {
    return counts.map((count, i) => ({
      label: edges[i] >= 1000 ? `${(edges[i] / 1000).toFixed(edges[i] % 1000 === 0 ? 0 : 1)}k` : `${edges[i]}`,
      rangeStart: edges[i],
      rangeEnd: edges[i + 1],
      count,
    }));
  }
  return { add, toArray, counts };
}

/** GC content (%) from accumulated base counts, excluding ambiguous bases. */
export function gcPercent(baseCounts) {
  const { A = 0, C = 0, G = 0, T = 0 } = baseCounts;
  const acgt = A + C + G + T;
  if (acgt === 0) return 0;
  return ((G + C) / acgt) * 100;
}

export function formatBases(n) {
  if (n >= 1e9) return `${(n / 1e9).toFixed(2)} Gb`;
  if (n >= 1e6) return `${(n / 1e6).toFixed(2)} Mb`;
  if (n >= 1e3) return `${(n / 1e3).toFixed(1)} kb`;
  return `${n} bp`;
}

export function formatBytes(n) {
  if (n >= 1e9) return `${(n / 1e9).toFixed(2)} GB`;
  if (n >= 1e6) return `${(n / 1e6).toFixed(2)} MB`;
  if (n >= 1e3) return `${(n / 1e3).toFixed(1)} KB`;
  return `${n} B`;
}
