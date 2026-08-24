/**
 * Downstream-format builders for an aligned set of sequences. All take
 * plain id[] + string[] rows (equal length, gaps as '-') and emit text in
 * formats that IQ-TREE / RAxML / MrBayes / MEGA accept directly.
 */
export function buildAlignedFasta(ids, rows, width = 80) {
  const out = [];
  for (let i = 0; i < ids.length; i++) {
    out.push(`>${ids[i]}`);
    const s = rows[i];
    for (let p = 0; p < s.length; p += width) out.push(s.slice(p, p + width));
  }
  return out.join("\n") + "\n";
}

/** Purine/pyrimidine transition pair, used for Clustal ':' marks. */
function isTransition(a, b) {
  const s = [a, b].sort().join("");
  return s === "AG" || s === "CT";
}

export function buildClustal(ids, rows, blockWidth = 60) {
  const nameWidth = Math.max(...ids.map((id) => Math.min(id.length, 24))) + 4;
  const len = rows[0]?.length ?? 0;
  const out = [
    "CLUSTAL W style alignment produced by Genome Explorer (in-browser progressive MSA)",
    "",
  ];

  for (let start = 0; start < len; start += blockWidth) {
    const end = Math.min(start + blockWidth, len);
    for (let i = 0; i < ids.length; i++) {
      out.push(ids[i].slice(0, 24).padEnd(nameWidth) + rows[i].slice(start, end));
    }
    let marks = "";
    for (let c = start; c < end; c++) {
      const chars = new Set();
      for (const row of rows) if (row[c] !== "-") chars.add(row[c]);
      if (chars.size === 0) marks += " ";
      else if (chars.size === 1) marks += "*";
      else if (chars.size === 2 && isTransition(...chars)) marks += ":";
      else marks += " ";
    }
    out.push(" ".repeat(nameWidth) + marks);
    out.push("");
  }
  return out.join("\n");
}

export function buildNexus(ids, rows) {
  const ntax = ids.length;
  const nchar = rows[0]?.length ?? 0;
  const lines = [
    "#NEXUS",
    "BEGIN DATA;",
    `    DIMENSIONS NTAX=${ntax} NCHAR=${nchar};`,
    "    FORMAT DATATYPE=DNA GAP=- MISSING=? ;",
    "    MATRIX",
  ];
  ids.forEach((id, i) => lines.push(`    ${safeName(id).padEnd(Math.max(...ids.map(safeName).map((s) => s.length)))} ${rows[i]}`));
  lines.push("    ;", "END;");
  return lines.join("\n") + "\n";
}

export function buildPhylip(ids, rows) {
  // PHYLIP wants <=10-char unique names on the same line as the sequence.
  const used = new Set();
  const names = ids.map((id) => {
    let base = safeName(id).slice(0, 10).padEnd(10);
    let name = base, k = 1;
    while (used.has(name)) { base = base.slice(0, 8); name = base + String(k).padStart(2, "0"); k++; }
    used.add(name);
    return name;
  });
  const nchar = rows[0]?.length ?? 0;
  const lines = [` ${ids.length} ${nchar}`];
  for (let i = 0; i < names.length; i++) lines.push(names[i] + rows[i]);
  return lines.join("\n") + "\n";
}

export function buildVariantsTsv(sites) {
  const head = ["alignment_column\tconsensus\talleles\tparsimony_informative\tinvolves_gap"];
  const body = sites.map((s) =>
    [s.pos + 1, s.consensus, s.alleles.join(" "), s.informative ? "yes" : "no", s.indel ? "yes" : "no"].join("\t")
  );
  return head.concat(body).join("\n") + "\n";
}

export function buildConsensusFasta(ids, consensus, width = 80) {
  const name = ids[0] ? `${ids[0]}_consensus` : "consensus";
  return buildAlignedFasta([name], [consensus], width);
}

function safeName(id) {
  return id.replace(/[^A-Za-z0-9_.-]/g, "_") || "seq";
}
