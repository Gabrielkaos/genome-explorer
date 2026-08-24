/**
 * Protein-level analysis for predicted CDS features.
 * All computations here are standard, deterministic bioinformatics:
 * genetic-code translation, average-residue molecular weight, pKa-based
 * iterative pI, Kyte-Doolittle hydropathy, and published heuristics for
 * transmembrane helices (sliding-window hydropathy) and signal peptides
 * (positive-inside rule + hydrophobic core + cleavage-site pattern).
 */

// Standard genetic code (NCBI transl_table=11 differs from table 1 only in
// start-codon usage, not in elongation; translation of internal codons is
// identical). Codons containing non-ACGT translate to "X".
const CODON_TABLE = (() => {
  const bases = ["T", "C", "A", "G"];
  const aas =
    "FFLLSSSSYY**CC*W" +
    "LLLLPPPPHHQQRRRR" +
    "IIIMTTTTNNKKSSRR" +
    "VVVVAAAADDEEGGGG";
  const table = {};
  let i = 0;
  for (const b1 of bases) {
    for (const b2 of bases) {
      for (const b3 of bases) {
        table[b1 + b2 + b3] = aas[i++];
      }
    }
  }
  return table;
})();

export function translate(seq) {
  let prot = "";
  for (let i = 0; i + 3 <= seq.length; i += 3) {
    const codon = seq.slice(i, i + 3);
    prot += CODON_TABLE[codon] ?? "X";
  }
  return prot.replace(/\*$/, ""); // drop terminal stop only
}

/** Average residue masses (Da), residue-only (i.e. post water-loss). */
const RESIDUE_MASS = {
  A: 71.0788, R: 156.1875, N: 114.1038, D: 115.0886, C: 103.1388,
  E: 129.1155, Q: 128.1307, G: 57.0519, H: 137.1411, I: 113.1594,
  L: 113.1594, K: 128.1741, M: 131.1926, F: 147.1766, P: 97.1167,
  S: 87.0782, T: 101.1051, W: 186.2132, Y: 163.1760, V: 99.1326,
};

export function molecularWeight(prot) {
  let m = 18.01524; // free termini contribute one H2O equivalent net
  for (const aa of prot) m += RESIDUE_MASS[aa] ?? 110.0;
  return m;
}

// pKa values (Expasy/EMBOSS convention).
const PKA = { nTerm: 9.69, cTerm: 2.34, C: 8.18, D: 3.65, E: 4.25, H: 6.0, K: 10.53, R: 12.48, Y: 10.07 };

function netCharge(prot, pH, counts) {
  const pos =
    Math.pow(10, PKA.nTerm - pH) / (1 + Math.pow(10, PKA.nTerm - pH)) +
    counts.K * (Math.pow(10, PKA.K - pH) / (1 + Math.pow(10, PKA.K - pH))) +
    counts.R * (Math.pow(10, PKA.R - pH) / (1 + Math.pow(10, PKA.R - pH))) +
    counts.H * (Math.pow(10, PKA.H - pH) / (1 + Math.pow(10, PKA.H - pH)));
  const neg =
    counts.D * (Math.pow(10, pH - PKA.D) / (1 + Math.pow(10, pH - PKA.D))) +
    counts.E * (Math.pow(10, pH - PKA.E) / (1 + Math.pow(10, pH - PKA.E))) +
    counts.C * (Math.pow(10, pH - PKA.C) / (1 + Math.pow(10, pH - PKA.C))) +
    counts.Y * (Math.pow(10, pH - PKA.Y) / (1 + Math.pow(10, pH - PKA.Y))) +
    Math.pow(10, pH - PKA.cTerm) / (1 + Math.pow(10, pH - PKA.cTerm));
  return pos - neg;
}

export function theoreticalPI(prot) {
  if (!prot.length) return null;
  const counts = { K: 0, R: 0, H: 0, D: 0, E: 0, C: 0, Y: 0 };
  for (const aa of prot) if (counts[aa] !== undefined) counts[aa]++;
  let lo = 0, hi = 14;
  for (let it = 0; it < 40; it++) {
    const mid = (lo + hi) / 2;
    if (netCharge(prot, mid, counts) > 0) lo = mid;
    else hi = mid;
  }
  return (lo + hi) / 2;
}

// Kyte-Doolittle hydropathy index per residue.
const KD = {
  A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5,
  Q: -3.5, E: -3.5, G: -0.4, H: -3.2, I: 4.5,
  L: 3.8, K: -3.9, M: 1.9, F: 2.8, P: -1.6,
  S: -0.8, T: -0.7, W: -0.9, Y: -1.3, V: 4.2,
};

export function gravy(prot) {
  if (!prot.length) return 0;
  let sum = 0;
  for (const aa of prot) sum += KD[aa] ?? 0;
  return sum / prot.length;
}

export function aromaticity(prot) {
  if (!prot.length) return 0;
  let n = 0;
  for (const aa of prot) if (aa === "F" || aa === "W" || aa === "Y") n++;
  return n / prot.length;
}

/**
 * Transmembrane-helix segments via sliding-window mean hydropathy.
 * Window 19 residues; a window is "transmembrane-like" when mean KD >= 1.6
 * (typical TM helices sit ~1.6-2.5 on this scale). Positive windows are merged
 * across gaps of up to 5 aa, then kept if the resulting segment is >= 17 aa.
 * This is deliberately simpler than a full HMM (TMHMM/Phobos); results are a
 * strong heuristic, not a membrane-topology prediction.
 */
export function findTmSegments(prot, { window = 19, threshold = 1.6, minLen = 17, maxGap = 5 } = {}) {
  if (prot.length < window) return [];
  const means = new Array(prot.length - window + 1);
  let running = 0;
  for (let i = 0; i < prot.length; i++) {
    running += KD[prot[i]] ?? 0;
    if (i >= window) running -= KD[prot[i - window]] ?? 0;
    if (i >= window - 1) means[i - window + 1] = running / window;
  }
  const raw = [];
  let start = null;
  for (let i = 0; i < means.length; i++) {
    if (means[i] >= threshold && start === null) start = i;
    else if (means[i] < threshold && start !== null) { raw.push([start, i + window - 1]); start = null; }
  }
  if (start !== null) raw.push([start, means.length + window - 1]);

  // merge close positive runs
  const merged = [];
  for (const seg of raw) {
    const last = merged[merged.length - 1];
    if (last && seg[0] - last[1] <= maxGap) last[1] = seg[1];
    else merged.push([...seg]);
  }
  return merged.filter(([a, b]) => b - a + 1 >= minLen).map(([a, b]) => ({ start: a + 1, end: b + 1 }));
}

const HYDRO_STARTS = new Set(["A", "G", "S", "T", "C", "V"]); // small residues allowed at cleavage -3/-1

/**
 * Signal-peptide heuristic combining the three classical signals:
 *  1) an N-region with at least 2 Lys/Arg within the first 8 aa (positive-inside rule)
 *  2) a hydrophobic core in roughly positions 6-24 (mean KD >= 1.4 over any 15-aa span)
 *  3) a plausible cleavage site: small residue at position -3 and -1 for some cut point between aa 12 and 32
 * Returns null or { score, cleavageAfter, note }.
 */
export function predictSignalPeptide(prot) {
  if (prot.length < 35) return null;

  let basicN = 0;
  for (let i = 0; i < Math.min(8, prot.length); i++) {
    if (prot[i] === "K" || prot[i] === "R") basicN++;
  }

  let bestCore = -Infinity;
  for (let s = 2; s <= Math.min(10, prot.length - 15); s++) {
    let sum = 0;
    for (let i = s; i < s + 15; i++) sum += KD[prot[i]] ?? 0;
    bestCore = Math.max(bestCore, sum / 15);
  }

  let cleavage = -1;
  for (let cut = 12; cut <= Math.min(32, prot.length - 2); cut++) {
    const m3 = prot[cut - 3], m1 = prot[cut - 1];
    if (HYDRO_STARTS.has(m3) && HYDRO_STARTS.has(m1)) { cleavage = cut; break; }
  }

  if (basicN < 2 || bestCore < 1.4 || cleavage < 0) return null;

  const score =
    Math.min(1, basicN / 3) * 0.3 +
    Math.min(1, Math.max(0, bestCore - 1.4) / 1.1) * 0.45 +
    (cleavage <= 26 ? 0.25 : 0.12);

  return {
    score: +score.toFixed(2),
    cleavageAfter: cleavage,
    note: `n-region +${basicN} K/R · core mean KD ${bestCore.toFixed(2)} · predicted cleavage after residue ${cleavage}`,
  };
}

/**
 * Subcellular-localization hint derived purely from the two predictions above.
 * Explicitly a coarse heuristic: bacteria can also carry lipoproteins,
 * beta-barrel outer-membrane proteins and intracellular-inclusion proteins
 * that these two signals do not capture.
 */
export function localizationHint(tmSegments, signal) {
  if (signal) return "secretory pathway (signal peptide detected)";
  if (tmSegments.length > 0) return `membrane (${tmSegments.length} TM helix${tmSegments.length > 1 ? "es" : ""})`;
  return "likely cytosolic";
}
