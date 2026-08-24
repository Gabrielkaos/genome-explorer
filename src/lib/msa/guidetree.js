/**
 * UPGMA (average-linkage) clustering over the k-mer distance matrix,
 * producing the binary guide tree that fixes the order in which sequences
 * and profiles are merged during progressive alignment. This mirrors
 * ClustalW's "guide tree from pairwise distances" stage, with cluster-pair
 * distances updated via the standard Lance-Williams average-linkage rule.
 */

export function upgma(dist, n) {
  const clusters = [];
  for (let i = 0; i < n; i++) {
    clusters.push({ id: i, size: 1, members: [i], alive: true, node: { leaf: true, idx: i } });
  }
  let nextId = n;
  let aliveCount = n;

  while (aliveCount > 1) {
    let bi = -1, bj = -1, bd = Infinity;
    for (let i = 0; i < clusters.length; i++) {
      const A = clusters[i];
      if (!A.alive) continue;
      for (let j = i + 1; j < clusters.length; j++) {
        const B = clusters[j];
        if (!B.alive) continue;
        // average linkage between the two clusters
        const d = pairAverage(dist, n, A.members, B.members);
        if (d < bd) { bd = d; bi = i; bj = j; }
      }
    }
    const A = clusters[bi], B = clusters[bj];
    A.alive = false; B.alive = false;
    clusters.push({
      id: nextId++,
      size: A.size + B.size,
      members: A.members.concat(B.members),
      node: { leaf: false, left: A.node, right: B.node },
      alive: true,
    });
    aliveCount--;
  }

  const root = clusters.find((c) => c.alive);
  return root ? root.node : null;

  function pairAverage(d, size, memA, memB) {
    let sum = 0;
    for (let a = 0; a < memA.length; a++) {
      const row = memA[a] * size;
      for (let b = 0; b < memB.length; b++) sum += d[row + memB[b]];
    }
    return sum / (memA.length * memB.length);
  }
}
