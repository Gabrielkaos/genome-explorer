/**
 * UPGMA (Unweighted Pair Group Method with Arithmetic mean, Sokal & Michener
 * 1958) over a pairwise distance matrix. Assumes a MOLECULAR CLOCK: the two
 * closest clusters merge, the new cluster's height is half its average
 * inter-cluster distance, and each child's branch length is the difference
 * between the new height and the child's own height. The result is a rooted,
 * ultrametric tree - which is also its main caveat: real data often violates
 * the clock assumption, and UPGMA will then infer wrong topologies that NJ
 * would get right.
 */

export function upgmaTree(dist, n, ids) {
  const clusters = [];
  for (let i = 0; i < n; i++) {
    clusters.push({
      id: i, members: [i], alive: true, height: 0,
      node: { isLeaf: true, name: ids[i], branchLen: 0, children: [] },
    });
  }
  let nextId = n;
  let aliveCount = n;
  let clamped = 0;

  while (aliveCount > 1) {
    let bi = -1, bj = -1, bd = Infinity;
    for (let i = 0; i < clusters.length; i++) {
      const A = clusters[i];
      if (!A.alive) continue;
      for (let j = i + 1; j < clusters.length; j++) {
        const B = clusters[j];
        if (!B.alive) continue;
        // arithmetic-mean distance over all member pairs (UPGMA linkage)
        let sum = 0;
        for (const a of A.members) {
          const row = a * n;
          for (const b of B.members) sum += dist[row + b];
        }
        const dAvg = sum / (A.members.length * B.members.length);
        if (dAvg < bd) { bd = dAvg; bi = i; bj = j; }
      }
    }
    const A = clusters[bi], B = clusters[bj];
    A.alive = false; B.alive = false;

    const height = bd / 2;
    const branchA = Math.max(0, height - A.height);
    const branchB = Math.max(0, height - B.height);
    if (height - A.height < 0 || height - B.height < 0) clamped++;

    clusters.push({
      id: nextId++,
      members: A.members.concat(B.members),
      alive: true,
      height,
      node: {
        isLeaf: false, name: null, branchLen: 0,
        children: [
          Object.assign({}, A.node, { branchLen: branchA }),
          Object.assign({}, B.node, { branchLen: branchB }),
        ],
      },
    });
    aliveCount--;
  }

  return { root: clusters.find((c) => c.alive).node, clampedLimbs: clamped };
}
