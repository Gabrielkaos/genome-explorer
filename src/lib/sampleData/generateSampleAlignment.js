/**
 * Deterministic synthetic multi-strain locus for trying the aligner without
 * an upload. Clearly labeled as simulated everywhere it appears in the UI.
 * Two clades (A-D vs E-H) share clade-specific SNPs, each strain carries
 * private substitutions plus occasional short indels - enough structure for
 * the guide tree to recover something tree-like from real signal.
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

export function generateSampleAlignment(seed = 7) {
  const rng = mulberry32(seed);
  const L = 900;

  const ancestor = Array.from({ length: L }, () => BASES[Math.floor(rng() * 4)]);

  const cladeSnps = { early: [], late: [] };
  for (let i = 0; i < 14; i++) cladeSnps.early.push(Math.floor(rng() * L));
  for (let i = 0; i < 9; i++) cladeSnps.late.push(Math.floor(rng() * L));

  const records = [];
  for (let s = 0; s < 8; s++) {
    const id = `strain_${String.fromCharCode(65 + s)}${s < 4 ? "_clade1" : "_clade2"}`;
    const seq = ancestor.slice();
    const shared = s < 4 ? cladeSnps.early : cladeSnps.late;
    for (const p of shared) seq[p] = mutate(seq[p], rng);

    const nPrivate = 3 + Math.floor(rng() * 5);
    for (let k = 0; k < nPrivate; k++) {
      const p = Math.floor(rng() * L);
      seq[p] = mutate(seq[p], rng);
    }
    if (rng() < 0.65) {
      const p = 50 + Math.floor(rng() * (L - 100));
      const len = 1 + Math.floor(rng() * 5);
      if (rng() < 0.5) seq.splice(p, len); // deletion
      else seq.splice(p, 0, ...Array.from({ length: len }, () => BASES[Math.floor(rng() * 4)])); // insertion
    }
    records.push({ id, seq: seq.join("") });
  }
  return records;
}

function mutate(b, rng) {
  const opts = BASES.filter((x) => x !== b);
  return opts[Math.floor(rng() * opts.length)];
}
