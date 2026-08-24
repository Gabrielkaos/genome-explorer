/**
 * Export builders for the genome browser: region/feature FASTA, feature TSV,
 * and a GFF3 dump of whatever is loaded.
 */

export function buildRegionFasta(seq, header, start0, end, strand = "+") {
  const s = Math.max(0, Math.floor(start0));
  const e = Math.min(seq.length, Math.ceil(end));
  let sub = seq.slice(s, e);
  if (strand === "-") {
    const COMP = { A: "T", C: "G", G: "C", T: "A", N: "N" };
    let rc = "";
    for (let i = sub.length - 1; i >= 0; i--) rc += COMP[sub[i]] ?? "N";
    sub = rc;
  }
  const lines = sub.match(/.{1,70}/g) ?? [""];
  return `>${header}\n${lines.join("\n")}\n`;
}

const escGff = (s) => String(s ?? "").replace(/([&;,=\n])/g, (m) => ({ "&": "&amp;", ",": "%2C", ";": "%3B", "=": "%3D", "\n": " " }[m]));

export function buildFeaturesGff3(features, contigLengths) {
  const lines = ["##gff-version 3"];
  for (const [id, len] of contigLengths) lines.push(`##sequence-region ${id} 1 ${len}`);
  for (const f of features) {
    lines.push([
      f.contigId, "genome-explorer", f.type, f.start + 1, f.end, ".", f.strand, ".",
      `ID=${escGff(f.locusTag)};locus_tag=${escGff(f.locusTag)}${f.product ? `;product=${escGff(f.product)}` : ""}`,
    ].join("\t"));
  }
  return lines.join("\n") + "\n";
}

export function buildFeaturesTsv(features) {
  const rows = [["locus_tag", "contig", "start", "end", "strand", "type", "length_nt", "product"].join("\t")];
  for (const f of features) {
    rows.push([f.locusTag, f.contigId, f.start + 1, f.end, f.strand, f.type, f.end - f.start, f.product ?? ""].join("\t"));
  }
  return rows.join("\n") + "\n";
}
