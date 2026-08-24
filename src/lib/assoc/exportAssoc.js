/**
 * Export builders for pan-GWAS results: a full per-feature TSV (every raw and
 * corrected statistic) and a filtered significant-hit list suitable for
 * downstream annotation or plotting elsewhere.
 */

const COLS_BINARY = [
  "gene", "product", "present_case", "present_control", "absent_case", "absent_control",
  "odds_ratio", "or_ci_low", "or_ci_high",
  "p_fisher", "q_bh_fdr", "q_bonferroni",
  "cmh_p_stratified", "cmh_or_common", "strata_used",
];

const COLS_CONTINUOUS = [
  "gene", "product", "mean_present", "mean_absent", "mean_difference",
  "t_statistic", "df",
  "p_welch", "q_bh_fdr", "q_bonferroni",
  "cmh_p_stratified", "cmh_or_common", "strata_used",
];

const num = (v, digits = 6) => (Number.isFinite(v) ? String(+(v.toPrecision(digits))) : "");

export function buildResultsTsv(rows, traitType) {
  const cols = traitType === "binary" ? COLS_BINARY : COLS_CONTINUOUS;
  const out = [cols.join("\t")];
  for (const r of rows) {
    const cells = traitType === "binary"
      ? [r.gene, r.product ?? "", r.a, r.b, r.c, r.d,
         num(r.or, 4), num(r.lo, 4), num(r.hi, 4),
         r.p.toExponential(4), r.q?.toExponential(4) ?? "", r.qBonf?.toExponential(4) ?? "",
         Number.isFinite(r.cmhP) ? r.cmhP.toExponential(4) : "", num(r.orMH, 4), r.strataUsed]
      : [r.gene, r.product ?? "", num(r.meanPresent, 5), num(r.meanAbsent, 5), num(r.meanDiff, 5),
         num(r.t, 5), num(r.df, 4),
         r.p?.toExponential(4) ?? "", r.q?.toExponential(4) ?? "", r.qBonf?.toExponential(4) ?? "",
         Number.isFinite(r.cmhP) ? r.cmhP.toExponential(4) : "", num(r.orMH, 4), r.strataUsed];
    out.push(cells.join("\t"));
  }
  return out.join("\n") + "\n";
}

export function buildSignificantTsv(rows, traitType, fdr) {
  const sig = rows.filter((r) => Number.isFinite(r.q) && r.q <= fdr);
  return buildResultsTsv(sig, traitType);
}
