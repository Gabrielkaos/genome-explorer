/**
 * Deterministic synthetic ALIGNED locus (equal-length rows, no indels) with
 * genuine hierarchical structure for trying the phylogenetics tool without
 * an upload: an outgroup on a long branch, two 4-taxon clades each with
 * sister-pair substructure, and private substitutions per strain. The true
 * topology is known by construction - useful for judging how many bootstrap
 * replicates it takes to recover which internodes.
 */

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const BASES = ["A", "C", "G", "T"];

function mutate(b, rng) {
  const opts = BASES.filter((x) => x !== b);
  return opts[Math.floor(rng() * opts.length)];
}

function pickPositions(L, count, rng, used) {
  const out = [];
  while (out.length < count) {
    const p = Math.floor(rng() * L);
    if (!used.has(p)) { used.add(p); out.push(p); }
  }
  return out;
}

export function generateSampleTreeData(seed = 11) {
  const rng = mulberry32(seed);
  const L = 1200;
  const ancestor = Array.from({ length: L }, () => BASES[Math.floor(rng() * 4)]);

  // Which mutation layers apply to each strain.
  const plan = [
    { id: "outgroup_X", clade: null, pairGroup: null, privates: 42 },
    { id: "clade1_A1", clade: "clade1", pairGroup: "A12", privates: 3 },
    { id: "clade1_A2", clade: "clade1", pairGroup: "A12", privates: 2 },
    { id: "clade1_A3", clade: "clade1", pairGroup: "A34", privates: 2 },
    { id: "clade1_A4", clade: "clade1", pairGroup: "A34", privates: 4 },
    { id: "clade2_B1", clade: "clade2", pairGroup: "B12", privates: 2 },
    { id: "clade2_B2", clade: "clade2", pairGroup: "B12", privates: 3 },
    { id: "clade2_B3", clade: "clade2", pairGroup: "B34", privates: 2 },
    { id: "clade2_B4", clade: "clade2", pairGroup: "B34", privates: 3 },
  ];

  const used = new Set();
  const layers = {
    clade1: pickPositions(L, 18, rng, used),
    clade2: pickPositions(L, 15, rng, used),
    A12: pickPositions(L, 7, rng, used),
    A34: pickPositions(L, 6, rng, used),
    B12: pickPositions(L, 8, rng, used),
    B34: pickPositions(L, 5, rng, used),
  };

  const records = plan.map((spec) => {
    const seq = ancestor.slice();
    for (const layerKey of [spec.clade, spec.pairGroup]) {
      if (!layerKey) continue;
      for (const p of layers[layerKey]) seq[p] = mutate(seq[p], rng);
    }
    for (const p of pickPositions(L, spec.privates, rng, used)) seq[p] = mutate(seq[p], rng);
    return { id: spec.id, seq: seq.join("") };
  });
  return records;
}
