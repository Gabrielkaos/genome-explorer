/**
 * Low-level sequence utilities shared by the assembly engine.
 */

const COMPLEMENT = { A: "T", C: "G", G: "C", T: "A", a: "t", c: "g", g: "c", t: "a", N: "N", n: "n" };

export function reverseComplement(seq) {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) out += COMPLEMENT[seq[i]] || "N";
  return out;
}

/** FNV-1a string hash, 32-bit. Fast, well-distributed, deterministic. */
export function hashKmer(km) {
  let h = 2166136261;
  for (let i = 0; i < km.length; i++) {
    h ^= km.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  return h >>> 0;
}

/**
 * Canonical hash for a k-mer: a k-mer and its reverse complement represent
 * the same underlying DNA duplex, so we always hash whichever orientation
 * sorts first (by hash value) and remember which strand that was. This is
 * what lets us find overlaps between reads regardless of which strand each
 * one happened to be sequenced from.
 */
export function canonicalKmerHash(km, revCompKm) {
  const hFwd = hashKmer(km);
  const hRev = hashKmer(revCompKm);
  return hFwd <= hRev ? { hash: hFwd, strand: 1 } : { hash: hRev, strand: -1 };
}
