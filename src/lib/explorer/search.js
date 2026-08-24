/**
 * Sequence-level utilities for the browser: IUPAC-aware motif search on both
 * strands, feature text search, and reverse-complement handling.
 */

const IUPAC = {
  A: "A", C: "C", G: "G", T: "T",
  R: "AG", Y: "CT", S: "GC", W: "AT", K: "GT", M: "AC",
  B: "CGT", D: "AGT", H: "ACT", V: "ACG", N: "ACGT",
};

export const COMPLEMENT = { A: "T", C: "G", G: "C", T: "A", N: "N" };

export function revcomp(seq) {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) out += COMPLEMENT[seq[i]] ?? "N";
  return out;
}

/** Regex source for an IUPAC motif pattern (already uppercased). */
function motifRegexSource(pattern) {
  let src = "";
  for (const ch of pattern) src += IUPAC[ch] ? `[${IUPAC[ch]}]` : "[N]";
  return src;
}

/**
 * Find all occurrences of `pattern` in both strands of `seq`.
 * Returns [{start0, end, strand, match}] sorted by position.
 */
export function findMotif(seq, pattern) {
  const p = pattern.toUpperCase().replace(/[^A-Z]/g, "");
  if (!p) return [];
  const re = new RegExp(motifRegexSource(p), "g");
  const hits = [];
  let m;
  while ((m = re.exec(seq)) !== null) {
    hits.push({ start0: m.index, end: m.index + p.length, strand: "+", match: m[0] });
    if (m.index === re.lastIndex) re.lastIndex++; // safety against zero-length
  }
  const rc = revcomp(seq);
  while ((m = re.exec(rc)) !== null) {
    // Map back to forward coordinates.
    const end1 = seq.length - m.index;          // exclusive on forward strand
    const start0 = end1 - p.length;
    hits.push({ start0: Math.max(0, start0), end: Math.min(seq.length, end1), strand: "-", match: revcomp(m[0]) });
    if (m.index === re.lastIndex) re.lastIndex++;
  }
  hits.sort((a, b) => a.start0 - b.start0 || (a.strand < b.strand ? -1 : 1));
  return hits;
}

/** Case-insensitive substring search across locus tag / product / type. */
export function searchFeatures(features, query) {
  const q = query.trim().toLowerCase();
  if (!q) return [];
  return features.filter((f) =>
    f.locusTag?.toLowerCase().includes(q) ||
    f.product?.toLowerCase().includes(q) ||
    f.type.toLowerCase() === q
  );
}
