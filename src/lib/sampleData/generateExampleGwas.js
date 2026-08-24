/**
 * Synthetic but structurally realistic example dataset for the association
 * tool: N strains across 5 phylogenetic clades, an antibiotic-resistance
 * phenotype, and three classes of genes:
 *
 *   - "true" causal determinants carried by MULTIPLE lineages (so their signal
 *     survives population-structure correction),
 *   - lineage marker genes present in entire clades (pure confounders whose
 *     naive signal VANISHES once you stratify by clade),
 *   - background/HGT genes at low frequency (nulls).
 *
 * The seeded RNG makes every regeneration identical for a given seed.
 */

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const CLADES = [
  { id: "ST11", size: 13, resistP: 0.92 },
  { id: "ST131", size: 12, resistP: 0.85 },
  { id: "ST10", size: 12, resistP: 0.45 },
  { id: "ST73", size: 10, resistP: 0.12 },
  { id: "ST95", size: 9, resistP: 0.08 },
];

const MARKER_PRODUCTS = [
  ["outer membrane protein OmpA", "core"],
  ["DNA polymerase III subunit alpha", "core"],
  ["chaperonin GroEL", "core"],
  ["ribosomal protein S15", "core"],
  ["ATP synthase subunit B", "core"],
];

export function generateExampleGwas(seed = 7) {
  const rng = mulberry32(seed);

  // ---- samples, clades, phenotypes --------------------------------------
  const samples = [];
  const sampleClade = [];
  const phenotypes = [];
  CLADES.forEach((clade, ci) => {
    for (let i = 1; i <= clade.size; i++) {
      samples.push(`${clade.id}_${String(i).padStart(2, "0")}`);
      sampleClade.push(ci);
      phenotypes.push(rng() < clade.resistP ? "Resistant" : "Susceptible");
    }
  });
  const n = samples.length;

  // ---- genes --------------------------------------------------------------
  const genes = [];

  const push = (name, values, product, klass) =>
    genes.push({ name, values: Uint8Array.from(values), present: values.reduce((a, b) => a + b, 0), product, klass });

  // TRUE determinants: carried across several lineages; phenotype is driven
  // by blaCTX-M-gene with ~6% misclassification noise.
  addCrossLineageGene("blaCTX-M-15", 0.75, 0.55, 0.30);
  addCrossLineageGene("aac(6')-Ib-cr", 0.45, 0.35, 0.18);

  // LINEAGE MARKERS (confounders): fixed per clade with rare dropout.
  let mi = 0;
  for (let ci = 0; ci < CLADES.length; ci++) {
    const [product] = MARKER_PRODUCTS[mi++ % MARKER_PRODUCTS.length];
    const name = `${CLADES[ci].id}-marker_${ci + 1}`;
    const values = sampleClade.map((c) => (c === ci && rng() > 0.04 ? 1 : 0));
    push(name, values, product, "marker");
  }
  // One extra marker shared by the two resistant-heavy clades (classic
  // clonal-expansion confounder).
  push(
    "ICEAbsence-region_1",
    sampleClade.map((c) => (c <= 1 && rng() > 0.03 ? 1 : 0)),
    "integrative conjugative element, hypothetical cargo",
    "marker",
  );

  // BACKGROUND / null genes at low frequency, unlinked to anything.
  for (let g = 1; g <= 40; g++) {
    const values = Array.from({ length: n }, () => (rng() < 0.16 ? 1 : 0));
    push(`hypothetical_bg${String(g).padStart(2, "0")}`, values, "hypothetical protein", "background");
  }

  function addCrossLineageGene(name, p1, p2, pOther) {
    const probs = [p1, p2, pOther * 1.4, pOther * 0.5, pOther * 0.4];
    const carrier = sampleClade.map((c) => rng() < probs[c]);
    const values = carrier.map((carr) => (carr ? 1 : 0));
    // Phenotype follows this class of genes with modest noise.
    push(name, values,
      name.startsWith("bla") ? "class A extended-spectrum beta-lactamase" : "aminoglycoside acetyltransferase",
      "causal");
  }

  // Re-derive phenotype from causal carriage (dominant, noisy) so the ground
  // truth is genuinely encoded in the data rather than asserted.
  const causalIdx = genes.filter((g) => g.klass === "causal");
  for (let i = 0; i < n; i++) {
    const dose = causalIdx.reduce((s, g) => s + g.values[i], 0);
    const pRes = Math.min(0.97, dose === 0 ? 0.06 : dose === 1 ? 0.72 : 0.94);
    phenotypes[i] = rng() < pRes ? "Resistant" : "Susceptible";
  }

  return {
    samples,
    genes,
    meta: [
      { name: "phenotype", kind: "binary", values: phenotypes },
      { name: "clade", kind: "text", values: sampleClade.map((c) => CLADES[c].id) },
    ],
    nCausal: causalIdx.length,
    note: "Simulated outbreak panel: two cross-lineage determinants drive resistance; clade markers are inherited confounders.",
  };
}
