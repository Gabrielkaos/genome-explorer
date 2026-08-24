/**
 * Statistical core for the pan-GWAS / genotype-phenotype association tool.
 * Everything here is standard, deterministic biostatistics suitable for
 * gene presence/absence (binary) or quantitative traits:
 *
 *  - Fisher's exact test (two-sided) on the 2x2 contingency table
 *  - Cochran-Armitage-free alternative for continuous traits: Welch's t-test
 *  - Odds ratio with Haldane-Anscombe correction + Woolf log 95% CI
 *  - Multiple-testing correction: Benjamini-Hochberg FDR and Bonferroni
 *  - Population-structure correction: Cochran-Mantel-Haenszel stratified test
 *    (each clade/lineage is a stratum; tests association WITHIN strata),
 *    plus the Mantel-Haenszel common odds-ratio estimate.
 */

/* ------------------------------- special fns ------------------------------- */

const LG_C = [
  76.18009172947146, -86.5053203294168, 24.01409824083091,
  -1.23173957245, 0.0012086509738662, -0.000005395239384953,
];

export function logGamma(x) {
  if (x <= 0) return NaN;
  let y = x;
  const tmp = x + 5.5 - (x + 0.5) * Math.log(x + 5.5);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) ser += LG_C[j] / ++y;
  return -tmp + Math.log((SQRT_2PI * ser) / x);
}

const SQRT_2PI = Math.sqrt(2 * Math.PI);

const lchoose = (n, k) =>
  lgammaCache(n + 1) - lgammaCache(k + 1) - lgammaCache(n - k + 1);

const lgammaCacheMap = new Map();
function lgammaCache(x) {
  let v = lgammaCacheMap.get(x);
  if (v === undefined) { v = logGamma(x); lgammaCacheMap.set(x, v); }
  return v;
}

/** Regularized incomplete beta function I_x(a,b) (continued fractions). */
export function betai(a, b, x) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const lbeta = logGamma(a + b) - logGamma(a) - logGamma(b);
  const front = Math.exp(lbeta + a * Math.log(x) + b * Math.log(1 - x));
  // Use symmetry when x is close to 1 for better convergence.
  if (x < (a + 1) / (a + b + 2)) return (front * betacf(a, b, x)) / a;
  return 1 - (front * betacf(b, a, 1 - x)) / b;
}

function betacf(a, b, x) {
  const MAXIT = 200, EPS = 3e-14, FPMIN = 1e-300;
  const qab = a + b, qap = a + 1, qam = a - 1;
  let c = 1, d = 1 - (qab * x) / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d; h *= d * c;
    aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return h;
}

/** Upper-tail chi-square survival for df=1 via erfc. */
export function chiSquareSFdf1(x) {
  if (!Number.isFinite(x)) return NaN;
  if (x <= 0) return 1;
  return erfc(Math.sqrt(x / 2));
}

export function erfc(x) {
  // Numerical Recipes erfc with fractional-error < 1.2e-7 everywhere.
  const z = Math.abs(x), t = 2 / (2 + z);
  const ty = 4 * t - 2;
  const cof = [
    -1.30265371978171, 0.64196979235649, 1.9476473204185836e-2,
    -9.561514786808631e-3, -9.46595344482036e-4, 3.66839497852761e-4,
    4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6,
    1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8, 6.529054439e-9,
    5.059343495e-9, -9.91364156e-10, -2.27365122e-10, 9.6467911e-11,
    2.394038e-12, -6.886027e-12, 8.94487e-13, 3.13092e-13, -1.12708e-13,
    3.81e-16, 7.106e-15,
  ];
  let d = 0, dd = 0;
  for (let j = cof.length - 1; j > 0; j--) {
    const tmp = d;
    d = ty * d - dd + cof[j];
    dd = tmp;
  }
  const ans = t * Math.exp(-z * z + 0.5 * (cof[0] + ty * d) - dd);
  return x >= 0 ? ans : 2 - ans;
}

/* ----------------------------- Fisher exact ------------------------------ */

/**
 * Two-sided Fisher exact test for the 2x2 table:
 *   rows    = feature present / absent
 *   columns = case / control
 *   [[a, b],
 *    [c, d]]
 * Sums all hypergeometric probabilities <= observed under fixed margins.
 */
export function fisherExactTwoSided(a, b, c, d) {
  const r1 = a + b, r2 = c + d, c1 = a + c, c2 = b + d, n = a + b + c + d;
  if (r1 === 0 || r2 === 0 || c1 === 0 || c2 === 0) return 1; // degenerate margins

  const kLo = Math.max(0, r1 - r2), kHi = Math.min(r1, c1);
  if (kHi < kLo) return 1;

  const pObs = hyperLogP(a, r1, c1, c2, n);
  let p = 0;
  for (let k = kLo; k <= kHi; k++) {
    const pk = hyperLogP(k, r1, c1, c2, n);
    // Additive epsilon: these are LOG probabilities (negative), so a
    // multiplicative slack would move the cut the wrong way.
    if (pk <= pObs + 1e-9) p += Math.exp(pk);
  }
  return Math.min(1, Math.max(p, Number.MIN_VALUE));
}

/** log P(X=k): X = #cases among the r1 feature-present samples. */
function hyperLogP(k, r1, c1, c2, n) {
  return lchoose(c1, k) + lchoose(c2, r1 - k) - lchoose(n, r1);
}

/**
 * Odds ratio with Haldane-Anscombe 0.5 correction (always applied so zero
 * cells stay finite) plus Woolf-log 95% CI.
 */
export function oddsRatioCI(a, b, c, d) {
  const A = a + 0.5, B = b + 0.5, C = c + 0.5, D = d + 0.5;
  const logOR = Math.log((A * D) / (B * C));
  const se = Math.sqrt(1 / A + 1 / B + 1 / C + 1 / D);
  const z = 1.959963984540054;
  return {
    or: Math.exp(logOR),
    lo: Math.exp(logOR - z * se),
    hi: Math.exp(logOR + z * se),
  };
}

/* ------------------------- Welch t (continuous) --------------------------- */

/**
 * Welch's unequal-variance t-test, two-sided, between trait values of
 * samples with vs without the feature. Returns t, df, p, mean difference.
 */
export function welchTTest(xs, ys) {
  const n1 = xs.length, n2 = ys.length;
  if (n1 < 2 || n2 < 2) return null;
  const m1 = xs.reduce((s, v) => s + v, 0) / n1;
  const m2 = ys.reduce((s, v) => s + v, 0) / n2;
  const v1 = xs.reduce((s, v) => s + (v - m1) * (v - m1), 0) / (n1 - 1);
  const v2 = ys.reduce((s, v) => s + (v - m2) * (v - m2), 0) / (n2 - 1);
  const seD = Math.sqrt(v1 / n1 + v2 / n2);
  if (!(seD > 0)) return null;
  const t = (m1 - m2) / seD;
  const num = (v1 / n1 + v2 / n2) ** 2;
  const den = (v1 * v1) / (n1 * n1 * (n1 - 1)) + (v2 * v2) / (n2 * n2 * (n2 - 1));
  const df = num / den;
  const p = betai(df / 2, 0.5, df / (df + t * t));
  return { t, df, meanDiff: m1 - m2, p };
}

/* ------------------- Cochran-Mantel-Haenszel (stratified) ----------------- */

/**
 * CMH test over K strata (e.g. clades/lineages). Each stratum contributes a
 * 2x2 table; the test asks whether, holding stratum constant, feature
 * presence associates with phenotype. Returns chi-square(1) p-value and the
 * Mantel-Haenszel common odds ratio estimate.
 *
 * tables: array of [a,b,c,d] per stratum (same layout as Fisher above).
 */
export function cmhTest(tables, continuityCorrection = true) {
  let sumA = 0, sumE = 0, sumV = 0, numOR = 0, denOR = 0, usable = 0;
  for (const [a, b, c, d] of tables) {
    const nT = a + b + c + d;
    if (nT === 0) continue;
    const col1 = a + c, row1 = a + b;
    if (col1 === 0 || col1 === nT || row1 === 0 || row1 === nT) continue; // non-informative stratum
    sumA += a;
    sumE += (row1 * col1) / nT;
    sumV += (row1 * col1 * (nT - row1) * (nT - col1)) / (nT * nT * (nT - 1));
    numOR += (a * d) / nT;
    denOR += (b * c) / nT;
    usable++;
  }
  if (usable === 0 || sumV <= 0) return { p: NaN, orMH: NaN, strataUsed: 0, chi2: NaN };
  const diff = Math.abs(sumA - sumE);
  const cc = continuityCorrection && diff >= 0.5 ? 0.5 : 0;
  const chi2 = (diff - cc) ** 2 / sumV;
  return {
    p: chiSquareSFdf1(chi2),
    orMH: denOR > 0 ? numOR / denOR : Infinity,
    strataUsed: usable,
    chi2,
  };
}

/* ----------------------- multiple-testing correction ---------------------- */

/**
 * Benjamini-Hochberg FDR adjustment. Input array of raw p-values (NaN
 * allowed; NaN outputs stay NaN). Returns new array of adjusted q-values.
 */
export function benjaminiHochberg(pvals) {
  const n = pvals.length;
  const out = new Array(n).fill(NaN);
  const order = [];
  for (let i = 0; i < n; i++) if (Number.isFinite(pvals[i])) order.push(i);
  order.sort((i, j) => pvals[i] - pvals[j]);
  let prev = 1;
  for (let rank = order.length; rank >= 1; rank--) {
    const idx = order[rank - 1];
    const q = Math.min(prev, (pvals[idx] * order.length) / rank);
    prev = q;
    out[idx] = Math.min(1, q);
  }
  return out;
}

export function bonferroni(pvals) {
  const n = pvals.filter(Number.isFinite).length || 1;
  return pvals.map((p) => (Number.isFinite(p) ? Math.min(1, p * n) : NaN));
}
