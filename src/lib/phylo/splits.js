/**
 * Split (bipartition) handling - the machinery behind nonparametric
 * bootstrap support values.
 *
 * Every internal edge of a tree divides the taxa into two sets. We
 * canonicalize each edge to its SMALLER side, which makes splits
 * comparable across trees regardless of where each tree happens to be
 * rooted - essential for NJ trees, whose root is arbitrary
 * (a trifurcation). Bootstrap replicates rebuild the whole pipeline
 * (resample columns -> distances -> tree) and every replicate votes for
 * the splits it contains; a split recovered in 87% of replicates gets
 * support 87. Trivial splits (one taxon vs everything else) carry no
 * shared-history signal and are excluded, per convention.
 */

/** Canonical representation of an edge: the sorted smaller-side member list.
 *  When both sides have exactly n/2 members the tie is broken
 *  lexicographically so complementary edges always map to the SAME key. */
export function canonKeyOf(membersSorted, total) {
  const ownKey = membersSorted.join(",");
  if (membersSorted.length * 2 < total) return { key: ownKey, size: membersSorted.length };
  if (membersSorted.length * 2 === total) {
    const mask = new Uint8Array(total);
    for (const m of membersSorted) mask[m] = 1;
    const comp = [];
    for (let i = 0; i < total; i++) if (!mask[i]) comp.push(i);
    const compKey = comp.join(",");
    return compKey < ownKey
      ? { key: compKey, size: comp.length }
      : { key: ownKey, size: membersSorted.length };
  }
  const mask = new Uint8Array(total);
  for (const m of membersSorted) mask[m] = 1;
  const comp = [];
  for (let i = 0; i < total; i++) if (!mask[i]) comp.push(i);
  return { key: comp.join(","), size: comp.length };
}

/**
 * Index leaves 0..k-1 in DFS order, then fill node._members (sorted leaf
 * indices under each node) for the whole tree.
 */
export function indexAndAssignMembers(tree) {
  let cursor = 0;
  const indexOf = new Map();
  (function assignIdx(node) {
    if (node.isLeaf) { indexOf.set(node, cursor++); return; }
    node.children.forEach(assignIdx);
  })(tree);

  (function assignMembers(node) {
    if (node.isLeaf) { node._members = [indexOf.get(node)]; return; }
    node.children.forEach(assignMembers);
    const merged = [].concat(...node.children.map((c) => c._members));
    merged.sort((a, b) => a - b);
    node._members = merged;
  })(tree);

  return { total: cursor };
}

/** All non-trivial split keys present in the tree (one per internal edge). */
export function collectSplits(tree) {
  const { total } = indexAndAssignMembers(tree);
  const keys = [];
  (function visitEdges(node) {
    if (node.isLeaf) return;
    for (const child of node.children) {
      const m = child._members;
      if (m.length > 1 && m.length < total) keys.push(canonKeyOf(m, total).key);
      visitEdges(child);
    }
  })(tree);
  return keys;
}

/**
 * Annotate every internal node except the root with `support`
 * (0-100) read off the replicate tally. Rooting-invariant because keys are
 * canonicalized to the smaller side.
 */
export function annotateSupports(tree, tally, reps) {
  const { total } = indexAndAssignMembers(tree);
  (function visit(node, isRoot) {
    if (node.isLeaf) return;
    node.children.forEach((ch) => visit(ch, false));
    if (!isRoot) {
      const hits = tally.get(canonKeyOf(node._members, total).key);
      node.support = hits != null ? Math.round((hits / reps) * 100) : 0;
    }
  })(tree, true);
  return tree;
}

/**
 * Midpoint rooting: place the root halfway along the longest path between
 * any two leaves (the tree diameter) - the standard way to root an
 * unrooted distance tree when no outgroup is available.
 */
export function midpointRoot(root) {
  const adj = new Map();
  function link(a, b, len) {
    if (!adj.has(a)) adj.set(a, []);
    if (!adj.has(b)) adj.set(b, []);
    adj.get(a).push({ node: b, len });
    adj.get(b).push({ node: a, len });
  }
  (function build(node) {
    for (const ch of node.children) link(node, ch, Math.max(0, ch.branchLen || 0));
    node.children.forEach(build);
  })(root);

  const anyNode = [...adj.keys()][0];
  if (!anyNode) return root;

  function farthestFrom(start) {
    const dist = new Map([[start, 0]]);
    const parent = new Map([[start, null]]);
    const stack = [start];
    let best = start;
    while (stack.length) {
      const cur = stack.pop();
      if (dist.get(cur) > dist.get(best)) best = cur;
      for (const { node, len } of adj.get(cur)) {
        if (dist.has(node)) continue;
        dist.set(node, dist.get(cur) + len);
        parent.set(node, cur);
        stack.push(node);
      }
    }
    return { best, dist, parent };
  }

  const sweepB = farthestFrom(farthestFrom(anyNode).best);
  const endB = sweepB.best;
  const diameter = sweepB.dist.get(endB);
  if (!(diameter > 0)) return root;

  // Path endA .. endB (via parents recorded from endB's sweep).
  const path = [];
  for (let cur = endB; cur != null; cur = sweepB.parent.get(cur)) path.push(cur);
  path.reverse(); // now endA ... endB

  const half = diameter / 2;
  let acc = 0, i = 0;
  while (i < path.length - 1) {
    const w = adj.get(path[i]).find((e) => e.node === path[i + 1])?.len ?? 0;
    if (acc + w >= half) break;
    acc += w;
    i++;
  }

  const cloneComponent = (start, exclude, limb) => {
    function rec(cur, cameFrom) {
      const children = [];
      for (const { node, len } of adj.get(cur)) {
        if (node === cameFrom || (cur === start && node === exclude)) continue;
        const child = rec(node, cur);
        child.branchLen = len;
        children.push(child);
      }
      return {
        isLeaf: children.length === 0,
        name: cur.name ?? null,
        branchLen: 0,
        support: undefined,
        children,
      };
    }
    const out = rec(start, null);
    out.branchLen = limb;
    return out;
  };

  // The cut may land inside an edge or exactly on a node; either way we
  // split that edge into two limbs (one possibly zero-length).
  const nodeA = path[i], nodeB = path[i + 1];
  const wAB = adj.get(nodeA).find((e) => e.node === nodeB)?.len ?? 0;
  const leftLen = half - acc;
  const rightLen = Math.max(0, wAB - leftLen);

  return {
    isLeaf: false,
    name: null,
    branchLen: 0,
    rootedByMidpoint: true,
    children: [cloneComponent(nodeA, nodeB, leftLen), cloneComponent(nodeB, nodeA, rightLen)],
  };
}
