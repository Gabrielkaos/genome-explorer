/**
 * Deterministic example genome for the browser: one circular ~64 kb replicon
 * carrying CDS features on both strands, two RNA features, and a handful of
 * planted IUPAC-searchable motifs so the search panel has something to find.
 */

function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const PRODUCTS = [
  "DNA-directed RNA polymerase subunit beta", "ABC transporter substrate-binding protein",
  "hypothetical protein", "hypothetical protein", "ribosomal protein L2",
  "site-specific recombinase", "phage tail family protein", "transposase IS30 family",
  "two-component system response regulator", "sigma-70 factor domain protein",
  "methyl-accepting chemotaxis protein", "short-chain dehydrogenase",
  "TonB-dependent receptor", "cold shock protein", "integrase core domain",
  "paratox regulatory protein", "DUF411 domain-containing protein", "peptidase M23",
];

export function generateExampleGenome(seed = 11) {
  const rng = mulberry32(seed);
  const LEN = 64000;
  const bases = "ACGT";
  let seq = "";
  for (let i = 0; i < LEN; i++) seq += bases[Math.floor(rng() * 4)];

  // Planted motifs: a palindromic repeat present at several sites.
  const motifSites = [4200, 9800, 15550, 27300, 40100, 52900];
  const MOTIF = "GAATTCGCGAATTC"; // self-reverse-complement-ish 14mer
  for (const s of motifSites) seq = seq.slice(0, s) + MOTIF + seq.slice(s + MOTIF.length);

  // Promoter-like -10 boxes scattered around.
  for (let i = 0; i < 40; i++) {
    const s = Math.floor(rng() * (LEN - 20));
    seq = seq.slice(0, s) + "TATAAT" + seq.slice(s + 6);
  }

  const features = [];
  let pos = 240;
  let gi = 0;
  while (pos < LEN - 400) {
    const kind = rng();
    if (kind < 0.055 && gi % 7 === 3) {
      const len = 1542;
      features.push(mkFeature(++gi, pos, pos + len, "+", "rRNA", "16S ribosomal RNA"));
      pos += len + 90;
      continue;
    }
    if (kind < 0.085 && gi % 5 === 2) {
      const len = 72;
      features.push(mkFeature(++gi, pos, pos + len, "-", "tRNA", `tRNA-${["Lys", "Gly", "Met"][gi % 3]}(${["TTT", "CCC", "CAT"][gi % 3]})`));
      pos += len + 55;
      continue;
    }
    const strand = rng() < 0.5 ? "+" : "-";
    const codons = Math.max(120, Math.round(160 + rng() * 460)); // ~40-200 aa
    const lenNt = codons * 3;
    const product = PRODUCTS[gi % PRODUCTS.length];
    features.push(mkFeature(++gi, pos, pos + lenNt, strand, "CDS", product));
    pos += lenNt + Math.floor(30 + rng() * 220);
  }

  function mkFeature(n, start, end, strand, type, product) {
    return {
      contigId: "pGEK_64kb",
      type,
      start,
      end: Math.min(end, LEN),
      strand,
      phase: ".",
      locusTag: `GEK_${String(n).padStart(4, "0")}`,
      product,
    };
  }

  return {
    id: "pGEK_64kb",
    desc: "Genome Explorer example plasmid (simulated, circular)",
    seq,
    circular: true,
    features,
  };
}
