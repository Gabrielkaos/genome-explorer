/**
 * Shared sequence alphabet for the MSA pipeline. Sequences are encoded into
 * compact Uint8Arrays once at the boundary (0=A, 1=C, 2=G, 3=T, 4=gap) and
 * every downstream stage works on codes. Anything that is not ACGT or an
 * explicit gap character - N, lowercase letters, whitespace, IUPAC ambiguity
 * codes - collapses to gap, which is the honest way to treat "base unknown"
 * in a sum-of-pairs scorer.
 */
export const CODE_CHARS = ["A", "C", "G", "T", "-"];

const LOOKUP = new Int8Array(128).fill(4);
LOOKUP[65] = 0; // A
LOOKUP[67] = 1; // C
LOOKUP[71] = 2; // G
LOOKUP[84] = 3; // T
LOOKUP[45] = 4; // '-'

export function encodeSeq(seq) {
  const out = new Uint8Array(seq.length);
  for (let i = 0; i < seq.length; i++) {
    const c = seq.charCodeAt(i) & 127;
    out[i] = c === 45 ? 4 : LOOKUP[c];
  }
  return out;
}

export function decodeChars(codes) {
  let s = "";
  for (let i = 0; i < codes.length; i++) s += CODE_CHARS[codes[i]];
  return s;
}
