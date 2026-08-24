/**
 * Derive per-sample lineage groups from a session-inferred phylogeny by
 * cutting the tree at a cumulative branch-length threshold. Every subtree
 * whose root sits at or beyond the cut becomes one stratum (clade) - the
 * standard "condition on lineages" trick that makes a stratified
 * Mantel-Haenszel test a usable population-structure control in-browser.
 */

function collectLeaves(node, out) {
  if (node.isLeaf) { out.push(node.name); return; }
  for (const c of node.children) collectLeaves(c, out);
}

export function treeClades(tree, cutoff) {
  const groups = [];
  (function walk(node, dist) {
    const d = dist + Math.max(0, node.branchLen || 0);
    if (node.isLeaf || d >= cutoff) {
      const leaves = [];
      collectLeaves(node, leaves);
      groups.push(leaves);
      return;
    }
    for (const c of node.children) walk(c, d);
  })(tree, 0);
  // A single group covering everything carries no stratification.
  const nameToGroup = new Map();
  groups.forEach((members, gi) => members.forEach((m) => nameToGroup.set(m, `clade_${gi + 1} (${members.length})`)));
  return nameToGroup;
}
