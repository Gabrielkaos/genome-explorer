/**
 * Neighbor-Joining (Saitou & Nei 1987), the standard distance-matrix
 * method. Greedily joins the pair of taxa that minimizes the Q criterion
 *
 *   Q(i,j) = (m-2)·d(i,j) - r_i - r_j ,  r_x = Σ_k d(x,k)
 *
 * which is guaranteed to recover the correct topology for additive
 * distance matrices. Branch lengths are computed so total path lengths
 * between leaves match the input distances; the result is an UNROOTED
 * binary tree whose root carries a trifurcation of the last three nodes.
 * Negative limb solutions (possible with non-additive data) are clamped
 * to zero and counted.
 */

export function neighborJoining(dist, n, ids) {
  let D = Float64Array.from(dist);
  const clampCount = { n: 0 };
  let clusters = [];
  for (let i = 0; i < n; i++) {
    clusters.push({ isLeaf: true, name: ids[i], branchLen: 0, children: [] });
  }

  let m = n;
  while (m > 3) {
    // r_i = row sums
    const r = new Float64Array(m);
    for (let i = 0; i < m; i++) {
      let s = 0;
      for (let j = 0; j < m; j++) if (j !== i) s += D[i * m + j];
      r[i] = s;
    }
    let bi = -1, bj = -1, bq = Infinity;
    for (let i = 0; i < m; i++) {
      for (let j = i + 1; j < m; j++) {
        const q = (m - 2) * D[i * m + j] - r[i] - r[j];
        if (q < bq) { bq = q; bi = i; bj = j; }
      }
    }
    const dij = D[bi * m + bj];
    const li = Math.max(0, 0.5 * dij + (r[bi] - r[bj]) / (2 * (m - 2)));
    const lj = Math.max(0, dij - li);
    if (li !== 0.5 * dij + (r[bi] - r[bj]) / (2 * (m - 2))) clampCount.n++;
    if (lj !== dij - li) clampCount.n++;

    const nodeA = clusters[bi], nodeB = clusters[bj];
    const u = {
      isLeaf: false,
      name: null,
      branchLen: 0,
      children: [
        Object.assign({}, nodeA, { branchLen: li }),
        Object.assign({}, nodeB, { branchLen: lj }),
      ],
    };

    // Build the reduced matrix: replace bi with u, drop bj.
    const next = new Float64Array((m - 1) * (m - 1));
    const mapTo = []; // old index -> new index (-1 for dropped bj)
    let t = 0;
    for (let k = 0; k < m; k++) {
      if (k === bj) { mapTo.push(-1); continue; }
      mapTo.push(k === bi ? t++ : t);
      if (k !== bi) t++;
    }
    for (let a = 0; a < m; a++) {
      if (a === bj) continue;
      for (let b = 0; b < m; b++) {
        if (b === bj || b === a) continue;
        next[mapTo[a] * (m - 1) + mapTo[b]] = D[a * m + b];
      }
    }
    // u's distances to everyone else.
    for (let a = 0; a < m; a++) {
      if (a === bi || a === bj) continue;
      const dau = (D[a * m + bi] + D[a * m + bj] - dij) / 2;
      next[mapTo[bi] * (m - 1) + mapTo[a]] = dau;
      next[mapTo[a] * (m - 1) + mapTo[bi]] = dau;
    }

    const nextClusters = [];
    for (let k = 0; k < m; k++) {
      if (k === bj) continue;
      nextClusters.push(k === bi ? u : clusters[k]);
    }
    clusters = nextClusters;
    D = next;
    m--;
  }

  // Final three nodes form the trifurcating (unrooted) root:
  // each limb = half the sum of its two incident pairwise distances minus
  // the opposite edge.
  const [x, y, z] = [0, 1, 2];
  const dxy = D[x * 3 + y], dxz = D[x * 3 + z], dyz = D[y * 3 + z];
  const lx = Math.max(0, (dxy + dxz - dyz) / 2);
  const ly = Math.max(0, (dxy + dyz - dxz) / 2);
  const lz = Math.max(0, (dxz + dyz - dxy) / 2);
  const root = {
    isLeaf: false,
    name: null,
    branchLen: 0,
    unrooted: true,
    children: [
      Object.assign({}, clusters[x], { branchLen: lx }),
      Object.assign({}, clusters[y], { branchLen: ly }),
      Object.assign({}, clusters[z], { branchLen: lz }),
    ],
  };
  return { root, clampedLimbs: clampCount.n };
}
