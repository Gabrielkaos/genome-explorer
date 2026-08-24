/**
 * Standards-compliant export builders for annotation results.
 * GFF3, GenBank flat file, protein FASTA (.faa), nucleotide FASTA (.ffn),
 * and a plain TSV summary. Gene coordinates are stored internally as
 * 0-based half-open [start, end); all exported formats use 1-based closed.
 */

const fmt1 = (g) => ({ start: g.start + 1, end: g.end });
const esc = (s) => String(s).replace(/([&;,=\n])/g, (m) => ({ "&": "&amp;", ",": "%2C", ";": "%3B", "=": "%3D", "\n": " " }[m]));

export function buildGff3(contigs, genes) {
  const lines = ["##gff-version 3"];
  for (const c of contigs) {
    if (c.seq?.length) lines.push(`##sequence-region ${c.id} 1 ${c.seq.length}`);
  }
  for (const g of genes) {
    const { start, end } = fmt1(g);
    // GFF phase: bases to skip from the feature start to reach whole codons.
    // Complete genes start on a start codon -> phase 0; partial 5' genes may
    // begin mid-codon when the ORF was trimmed at a contig edge.
    const phase = g.partial ? (3 - (g.lengthNt % 3)) % 3 : 0;
    const attrs = [
      `ID=${esc(g.locusTag)}`,
      `locus_tag=${esc(g.locusTag)}`,
      `product=${esc(g.product)}`,
      `start_codon=${g.startCodon ?? "none"}`,
      `rbs_score=${g.rbsScore}`,
      `score=${g.score}`,
    ];
    if (g.partial) attrs.push(`partial=${g.partial}`);
    if (g.tmSegments?.length) attrs.push(`tm_helices=${g.tmSegments.length}`);
    if (g.signal) attrs.push(`signal_peptide=predicted`);
    lines.push([g.contigId, "genome-explorer", "CDS", start, end, ".", g.strand, phase, attrs.join(";")].join("\t"));
  }
  return lines.join("\n") + "\n";
}

function wrapTranslation(prot, firstIndent, contIndent) {
  const out = [];
  let first = true;
  for (let i = 0; i < prot.length; i += 60) {
    out.push((first ? firstIndent : contIndent) + prot.slice(i, i + 60));
    first = false;
  }
  return out;
}

const MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"];

export function buildGenbank(contigs, genesByContig, { locusPrefix = "GE" } = {}) {
  const d = new Date();
  const date = `${String(d.getDate()).padStart(2, "0")}-${MONTHS[d.getMonth()]}-${d.getFullYear()}`;
  const blocks = [];

  for (const c of contigs) {
    const genes = genesByContig.get(c.id) ?? [];
    const topology = c.circular ? "circular" : "linear";
    const name = c.id.slice(0, 16).padEnd(16);
    const lenStr = String(c.seq.length).padStart(11);

    const b = [];
    b.push(`LOCUS       ${name}${lenStr} bp    DNA     ${topology} BCT ${date}`);
    b.push(`DEFINITION  ${c.id}`);
    b.push(`ACCESSION   ${c.id}`);
    b.push(`VERSION     ${c.id}`);
    b.push(`KEYWORDS    .`);
    b.push(`SOURCE      .`);
    b.push(`  ORGANISM  .`);
    b.push(`COMMENT     Annotated in-browser by Genome Explorer (de novo bacterial`);
    b.push(`            gene calling; transl_table=11). Products are hypothetical:`);
    b.push(`            no similarity search against reference databases was run.`);
    b.push(`FEATURES             Location/Qualifiers`);
    b.push(`     source          1..${c.seq.length}`);
    b.push(`                     /organism="${locusPrefix} sample"`);
    b.push(`                     /mol_type="genomic DNA"`);
    for (const g of genes) {
      const { start, end } = fmt1(g);
      const loc = g.strand === "-" ? `complement(${start}..${end})` : `${start}..${end}`;
      b.push(`     CDS             ${loc}`);
      b.push(`                     /locus_tag="${g.locusTag}"`);
      b.push(`                     /product="${g.product}"`);
      b.push(`                     /codon_start=1`);
      b.push(`                     /transl_table=11`);
      if (g.startCodon) b.push(`                     /start_codon="${g.startCodon}"`);
      if (g.partial) b.push(`                     /partial="true"`);
      const trans = wrapTranslation(g.protSeq, `                     /translation="`, `                      `);
      trans[trans.length - 1] += '"';
      b.push(...trans);
    }
    b.push(`ORIGIN      `);
    let idx = 1;
    const seqLower = c.seq.toLowerCase();
    for (; idx <= seqLower.length; idx += 60) {
      const chunk = seqLower.slice(idx - 1, idx + 59);
      const groups = chunk.match(/.{1,10}/g).join(" ");
      b.push(`${String(idx).padStart(9)} ${groups}`);
    }
    b.push(`//`);
    blocks.push(b.join("\n"));
  }
  return blocks.join("\n") + "\n";
}

export function buildFaa(genes) {
  const recs = genes.map((g) => {
    const { start, end } = fmt1(g);
    const head = [
      g.locusTag,
      g.product,
      `contig=${g.contigId}`,
      `coords=${start}..${end}`,
      `strand=${g.strand}`,
      g.startCodon ? `start=${g.startCodon}` : null,
      g.partial ? `partial=${g.partial}` : null,
      `rbs=${g.rbsScore}`,
      g.tmSegments?.length ? `tm_helices=${g.tmSegments.length}` : null,
      g.signal ? "signal_peptide=yes" : null,
    ].filter(Boolean).join(" ");
    return `>${head}\n${(g.protSeq.match(/.{1,60}/g) || []).join("\n")}`;
  });
  return recs.join("\n") + "\n";
}

export function buildFfn(genes) {
  const recs = genes.map((g) => {
    const { start, end } = fmt1(g);
    const head = [
      g.locusTag,
      g.product,
      `contig=${g.contigId}`,
      `coords=${start}..${end}`,
      `strand=${g.strand}`,
    ].join(" ");
    return `>${head}\n${(g.dnaSeq.match(/.{1,70}/g) || []).join("\n")}`;
  });
  return recs.join("\n") + "\n";
}

export function buildTsv(genes) {
  const cols = ["locus_tag", "contig", "start", "end", "strand", "length_nt", "length_aa", "start_codon", "rbs_score", "gc_percent", "mw_kda", "pi", "gravy", "tm_helices", "signal_peptide", "localization", "partial", "product"];
  const rows = [cols.join("\t")];
  for (const g of genes) {
    rows.push([
      g.locusTag,
      g.contigId,
      g.start + 1,
      g.end,
      g.strand,
      g.lengthNt,
      g.lengthAa,
      g.startCodon ?? "",
      g.rbsScore,
      g.gcContent,
      (g.mw / 1000).toFixed(1),
      g.pi ?? "",
      g.gravy,
      g.tmSegments?.length ?? 0,
      g.signal ? "yes" : "no",
      g.localization,
      g.partial ?? "",
      g.product,
    ].join("\t"));
  }
  return rows.join("\n") + "\n";
}
