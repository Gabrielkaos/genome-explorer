import { fisherExactTwoSided, oddsRatioCI, welchTTest, cmhTest, benjaminiHochberg, bonferroni } from "../lib/assoc/stats.js";

/**
 * Pan-GWAS worker: one association test per feature row.
 * Binary traits -> Fisher exact (+ odds ratio/CI); continuous -> Welch t.
 * If a per-sample stratum vector is supplied (clade/lineage labels), also
 * runs the Cochran-Mantel-Haenszel stratified test as the structure-aware
 * statistic. Multiple-testing correction over all tested features.
 */
self.onmessage = async (e) => {
  if (e.data?.type !== "assoc") return;
  try {
    await runAssoc(e.data.data);
  } catch (err) {
    self.postMessage({ type: "error", message: String(err?.message || err) });
  }
};

async function runAssoc(data) {
  const t0 = performance.now();
  const tick = () => new Promise((r) => setTimeout(r, 0));
  const report = (m) => self.postMessage({ type: "progress", ...m });

  const { samples, geneNames, matrix, traitValues, traitType, strata, products } = data;
  const nS = samples.length;
  const nG = geneNames.length;

  report({ stage: "tests", pct: 3 });
  await tick();

  // Precompute strata layout.
  const stratumIds = [];
  const stratumOf = [];
  if (strata && strata.length === nS) {
    const map = new Map();
    for (let s = 0; s < nS; s++) {
      const key = strata[s];
      if (!map.has(key)) { map.set(key, stratumIds.length); stratumIds.push(key); }
      stratumOf[s] = map.get(key);
    }
  }

  const rows = new Array(nG);
  let nCase = 0, nCtrl = 0;

  for (let g = 0; g < nG; g++) {
    const off = g * nS;
    const row = {
      gene: geneNames[g],
      product: products?.[g] ?? null,
      a: 0, b: 0, c: 0, d: 0,
      or: NaN, lo: NaN, hi: NaN,
      p: NaN,
      meanDiff: null,
      cmhP: NaN, orMH: NaN, strataUsed: 0,
    };

    if (traitType === "binary") {
      for (let s = 0; s < nS; s++) {
        const present = matrix[off + s] === 1;
        const isCase = traitValues[s] === 1;
        if (present && isCase) row.a++;
        else if (present) row.b++;
        else if (isCase) row.c++;
        else row.d++;
      }
      nCase = Math.max(nCase, row.a + row.c);
      nCtrl = Math.max(nCtrl, row.b + row.d);
      const f = fisherExactTwoSided(row.a, row.b, row.c, row.d);
      row.p = f;
      ({ or: row.or, lo: row.lo, hi: row.hi } = oddsRatioCI(row.a, row.b, row.c, row.d));

      if (stratumIds.length) {
        const tables = stratumIds.map(() => [0, 0, 0, 0]);
        for (let s = 0; s < nS; s++) {
          const st = tables[stratumOf[s]];
          const present = matrix[off + s] === 1;
          const isCase = traitValues[s] === 1;
          if (present && isCase) st[0]++;
          else if (present) st[1]++;
          else if (isCase) st[2]++;
          else st[3]++;
        }
        const cmh = cmhTest(tables);
        row.cmhP = cmh.p; row.orMH = cmh.orMH; row.strataUsed = cmh.strataUsed;
      }
    } else {
      const xs = [], ys = [];
      for (let s = 0; s < nS; s++) (matrix[off + s] === 1 ? xs : ys).push(traitValues[s]);
      const wt = welchTTest(xs, ys);
      if (wt) {
        row.p = wt.p;
        row.meanDiff = wt.meanDiff;
        row.t = wt.t;
        row.df = wt.df;
        row.meanPresent = xs.length ? wtMean(xs) : NaN;
        row.meanAbsent = ys.length ? wtMean(ys) : NaN;
        // Effect direction shown via signed mean difference; CI not applicable.
        row.or = NaN; row.lo = NaN; row.hi = NaN;
      }
    }

    rows[g] = row;

    if (g % 400 === 0) {
      report({ stage: "tests", pct: 4 + Math.round((g / nG) * 82), detail: `${g.toLocaleString()} / ${nG.toLocaleString()} features` });
      await tick();
    }
  }

  report({ stage: "correction", pct: 92, detail: "Benjamini-Hochberg FDR" });
  await tick();
  const raw = rows.map((r) => r.p);
  const qs = benjaminiHochberg(raw);
  const qb = bonferroni(raw);
  for (let i = 0; i < nG; i++) { rows[i].q = qs[i]; rows[i].qBonf = qb[i]; }

  rows.sort((x, y) => (Number.isFinite(x.p) && Number.isFinite(y.p) ? x.p - y.p : Number.isFinite(x.p) ? -1 : 1));

  report({ stage: "done", pct: 100 });
  self.postMessage({
    type: "done",
    rows,
    meta: {
      nSamples: nS,
      nGenes: nG,
      nCase,
      nCtrl,
      traitType,
      hasStrata: !!stratumIds.length,
      nStrata: stratumIds.length,
      ms: performance.now() - t0,
    },
  });
}

function wtMean(arr) {
  return arr.reduce((s, v) => s + v, 0) / arr.length;
}
