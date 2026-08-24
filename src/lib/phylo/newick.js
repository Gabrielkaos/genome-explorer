/**
 * Newick string serialization & parsing - the lingua franca of tree
 * exchange between this tool and IQ-TREE / RAxML / MEGA / BEAST /
 * FigTree / iTOL. Internal numeric labels are interpreted as bootstrap
 * support values (the near-universal convention), quoted names are
 * supported, branch lengths are optional.
 */

/** Serialize a tree object ({isLeaf,name,branchLen,children,support}) to Newick. */
export function toNewick(root, { includeSupports = true, digits = 6 } = {}) {
  function fmtLen(v) {
    if (v == null || !Number.isFinite(v)) return "";
    const s = Number(v.toFixed(digits)).toString();
    return `:${s}`;
  }
  function rec(node) {
    if (node.isLeaf) return `${maybeQuote(node.name ?? "")}${fmtLen(node.branchLen)}`;
    const kids = node.children.map(rec).join(",");
    const label = includeSupports && Number.isFinite(node.support) ? String(Math.round(node.support)) : "";
    return `(${kids})${label}${fmtLen(node.branchLen)}`;
  }
  return `${rec(root)};`;
}

function maybeQuote(name) {
  if (/^[A-Za-z0-9_\-.|/[\]]+$/.test(name)) return name;
  return `"${String(name).replace(/"/g, '""')}"`;
}

/**
 * Parse a Newick string into the same tree shape used everywhere here.
 * Numeric internal labels <= 100 become `support`. Throws Error with
 * position info on malformed input.
 */
export function parseNewick(text) {
  const src = text.trim().replace(/;\s*$/, "");
  let pos = 0;
  const n = src.length;

  function fail(msg) {
    throw new Error(`Invalid Newick at character ${pos + 1}: ${msg}`);
  }
  function skipWs() {
    while (pos < n && /\s/.test(src[pos])) pos++;
  }
  function peek() {
    skipWs();
    return pos < n ? src[pos] : null;
  }
  function readName() {
    skipWs();
    if (pos >= n) fail("unexpected end of input");
    const ch = src[pos];
    if (ch === '"' || ch === "'") {
      const quote = ch;
      pos++;
      let out = "";
      while (pos < n) {
        if (src[pos] === quote) {
          if (quote === '"' && src[pos + 1] === '"') { out += '"'; pos += 2; continue; }
          pos++;
          return out;
        }
        out += src[pos++];
      }
      fail("unterminated quoted name");
    }
    let out = "";
    while (pos < n && !"(),:;".includes(src[pos]) && !/\s/.test(src[pos])) out += src[pos++];
    return out;
  }
  function readLength() {
    skipWs();
    let num = "";
    while (pos < n && /[0-9eE+\-.]/.test(src[pos])) num += src[pos++];
    const v = parseFloat(num);
    if (!Number.isFinite(v)) fail(`bad branch length "${num}"`);
    return v;
  }

  function parseNode() {
    const c = peek();
    if (c === "(") {
      pos++; // consume (
      const children = [];
      for (;;) {
        children.push(parseNode());
        const d = peek();
        if (d === ",") { pos++; continue; }
        if (d === ")") { pos++; break; }
        fail("expected ',' or ')' inside internal node");
      }
      let label = null;
      if (peek() && !"(),:;".includes(peek())) label = readName();
      let len = 0;
      if (peek() === ":") { pos++; len = readLength(); }
      const support = label != null && /^\d+(\.\d+)?$/.test(label) ? Math.round(parseFloat(label)) : undefined;
      return { isLeaf: false, name: null, branchLen: len, support, children };
    }
    const name = readName();
    let len = 0;
    if (peek() === ":") { pos++; len = readLength(); }
    if (!name) fail("empty leaf name");
    return { isLeaf: true, name, branchLen: len, children: [] };
  }

  const root = parseNode();
  skipWs();
  if (pos < n) fail(`unexpected trailing text "${src.slice(pos, pos + 12)}"`);

  // Normalize: a degree-2 root (some programs emit "(A,B);") is kept as-is;
  // anything with >=3 children at the root is treated as unrooted.
  if (!root.isLeaf) root.unrooted = root.children.length > 2;
  return root;
}
