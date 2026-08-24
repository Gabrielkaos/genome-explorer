/**
 * Profiles and profile-profile global alignment - the engine of progressive
 * MSA. A profile is an aligned set of sequences summarized as per-column
 * residue counts; aligning two profiles is a global DP (Gotoh affine gaps)
 * whose cell scores are the AVERAGE sum-of-pairs score between the two
 * columns involved:
 *
 *   diag(x,y) = sum_a f_x(a) * t_y(a) + gapGapBonus * f_x(-) * f_y(-)
 *               t_y(a) = match*f_y(a) - mismatch*(residues_y - f_y(a))
 *
 * This is ClustalW's profile-score formulation: only residue-residue pairs
 * are charged on the diagonal, gap-gap pairs get a small BONUS (aligning
 * gaps with gaps is evidence of a shared indel event), and residue-gap
 * pairs contribute nothing there because those gaps were already priced
 * when their run was opened. Insertion runs cost (gapOpen + gapExtend*k)
 * scaled by the fraction of REAL residues in the opposing column, so new
 * gaps slide toward pre-existing gap columns instead of splitting them.
 */
import { encodeSeq } from "./alphabet.js";

export function singleProfile(seq, origIdx) {
  const enc = typeof seq === "string" ? encodeSeq(seq) : seq;
  const len = enc.length;
  const cols = new Uint32Array(5 * len);
  for (let i = 0; i < len; i++) cols[i * 5 + enc[i]]++;
  return { len, cols, members: [{ origIdx, chars: enc }] };
}

function fractions(profile) {
  const f = new Float32Array(5 * profile.len);
  const n = profile.members.length;
  for (let c = 0; c < profile.len * 5; c++) f[c] = profile.cols[c] / n;
  return f;
}

/** Global profile-profile alignment. Returns a NEW merged profile. */
export function mergeProfiles(P, Q, params) {
  const m = P.len, n = Q.len;
  if (!m || !n) throw new Error("Cannot align an empty profile.");

  if (m * n > 160e6) {
    throw new Error(
      `Profile DP too large (${m} x ${n} cells). Try fewer or shorter sequences - ` +
      `the in-browser aligner targets gene-to-small-genome scale inputs.`
    );
  }

  const { match = 2, mismatch = 1, gapOpen = 8, gapExtend = 1, gapGapBonus = 0.5 } = params;
  const coef = match + mismatch;

  const fx = fractions(P);
  const fy = fractions(Q);

  // Per-column score vectors: tx[i*5+a] = score of one P-residue `a`
  // against everything in Q-column j (and vice versa).
  const tx = new Float32Array(5 * m);
  const ty = new Float32Array(5 * n);
  for (const [f, t, len] of [[fx, tx, m], [fy, ty, n]]) {
    for (let c = 0; c < len; c++) {
      let resTot = 0;
      for (let a = 0; a < 4; a++) resTot += f[c * 5 + a];
      for (let a = 0; a < 4; a++) t[c * 5 + a] = f[c * 5 + a] * coef - resTot * mismatch;
    }
  }

  // Real-residue fraction per column: the weight charged when the OTHER
  // profile inserts a gap run opposite this column.
  const wx = new Float32Array(m), wy = new Float32Array(n);
  for (let c = 0; c < m; c++) wx[c] = 1 - fx[c * 5 + 4];
  for (let c = 0; c < n; c++) wy[c] = 1 - fy[c * 5 + 4];

  const NEG = -Infinity;
  const trace = new Uint8Array(m * n); // bits0-1 M-pred | 2-3 Ix | 4-5 Iy (3 = extended)
  let prevM = new Float64Array(n + 1), prevX = new Float64Array(n + 1), prevY = new Float64Array(n + 1);
  let curM = new Float64Array(n + 1), curX = new Float64Array(n + 1), curY = new Float64Array(n + 1);

  prevM[0] = 0;
  prevM.fill(NEG, 1);
  prevX.fill(NEG);
  prevY.fill(NEG);

  for (let i = 1; i <= m; i++) {
    curM[0] = NEG;
    curY[0] = NEG;
    curX[0] = 0; // free leading-gap chain down the left edge

    const ti = (i - 1) * 5;
    for (let j = 1; j <= n; j++) {
      const tj = (j - 1) * 5;
      const cell = (i - 1) * n + (j - 1);

      let s = 0;
      for (let a = 0; a < 4; a++) s += fx[ti + a] * ty[tj + a];
      s += gapGapBonus * fx[ti + 4] * fy[tj + 4];

      let best = prevM[j - 1], arg = 0;
      if (prevX[j - 1] > best) { best = prevX[j - 1]; arg = 1; }
      if (prevY[j - 1] > best) { best = prevY[j - 1]; arg = 2; }
      curM[j] = best === NEG ? NEG : best + s;
      trace[cell] |= arg;

      const oSrcX = prevY[j] >= prevM[j] ? 2 : 0;
      const oValX = prevY[j] >= prevM[j] ? prevY[j] : prevM[j];
      const openX = oValX === NEG ? NEG : oValX - gapOpen * wy[j - 1];
      const extX = prevX[j] === NEG ? NEG : prevX[j] - gapExtend * wy[j - 1];
      if (extX !== NEG && extX >= openX) { curX[j] = extX; trace[cell] |= (3 << 2); }
      else if (openX !== NEG) { curX[j] = openX; trace[cell] |= (oSrcX << 2); }
      else curX[j] = NEG;

      const oSrcY = curX[j - 1] >= curM[j - 1] ? 1 : 0;
      const oValY = curX[j - 1] >= curM[j - 1] ? curX[j - 1] : curM[j - 1];
      const openY = oValY === NEG ? NEG : oValY - gapOpen * wx[i - 1];
      const extY = curY[j - 1] === NEG ? NEG : curY[j - 1] - gapExtend * wx[i - 1];
      if (extY !== NEG && extY >= openY) { curY[j] = extY; trace[cell] |= (3 << 4); }
      else if (openY !== NEG) { curY[j] = openY; trace[cell] |= (oSrcY << 4); }
      else curY[j] = NEG;
    }

    [prevM, curM] = [curM, prevM];
    [prevX, curX] = [curX, prevX];
    [prevY, curY] = [curY, prevY];
  }

  let state = 0;
  let bestEnd = prevM[n];
  if (prevX[n] > bestEnd) { bestEnd = prevX[n]; state = 1; }
  if (prevY[n] > bestEnd) state = 2;

  const ops = [];
  let i = m, j = n;
  while (i > 0 && j > 0) {
    const cell = (i - 1) * n + (j - 1);
    if (state === 0) {
      ops.push(0);
      state = trace[cell] & 3;
      i--; j--;
    } else if (state === 1) {
      ops.push(1);
      const info = (trace[cell] >> 2) & 3;
      state = info === 3 ? 1 : info;
      i--;
    } else {
      ops.push(2);
      const info = (trace[cell] >> 4) & 3;
      state = info === 3 ? 2 : info;
      j--;
    }
  }
  while (i > 0) { ops.push(1); i--; }
  while (j > 0) { ops.push(2); j--; }
  ops.reverse();

  return buildMerged(P, Q, ops);
}

function buildMerged(P, Q, ops) {
  const L = ops.length;
  const cols = new Uint32Array(5 * L);
  let pI = 0, qI = 0;
  for (let k = 0; k < L; k++) {
    if (ops[k] === 0) {
      for (let a = 0; a < 5; a++) cols[k * 5 + a] = P.cols[pI * 5 + a] + Q.cols[qI * 5 + a];
      pI++; qI++;
    } else if (ops[k] === 1) {
      cols.set(P.cols.subarray(pI * 5, pI * 5 + 5), k * 5);
      pI++;
    } else {
      cols.set(Q.cols.subarray(qI * 5, qI * 5 + 5), k * 5);
      qI++;
    }
  }

  const mapP = new Int32Array(P.len).fill(-1);
  const mapQ = new Int32Array(Q.len).fill(-1);
  pI = 0; qI = 0;
  for (let k = 0; k < L; k++) {
    if (ops[k] === 0) { mapP[pI++] = k; mapQ[qI++] = k; }
    else if (ops[k] === 1) mapP[pI++] = k;
    else mapQ[qI++] = k;
  }

  const members = [];
  for (const mem of P.members) members.push({ origIdx: mem.origIdx, chars: remap(mem.chars, mapP, L) });
  for (const mem of Q.members) members.push({ origIdx: mem.origIdx, chars: remap(mem.chars, mapQ, L) });
  return { len: L, cols, members };
}

function remap(oldChars, mapArr, L) {
  const out = new Uint8Array(L).fill(4);
  for (let c = 0; c < oldChars.length; c++) {
    const nc = mapArr[c];
    if (nc >= 0) out[nc] = oldChars[c];
  }
  return out;
}
