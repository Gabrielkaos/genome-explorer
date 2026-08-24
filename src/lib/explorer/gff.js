/**
 * Annotation parsers for the genome browser: GFF3 (and tolerant GFF2) plus
 * GenBank flat files. Coordinates are normalized to 0-based half-open
 * [start, end) internally, matching the rest of the app.
 *
 * GFF3:  seqid source type start end score strand phase attributes
 * GenBank: FEATURES table CDS/gene/rRNA/tRNA/etc. entries + ORIGIN sequence,
 * parsed per LOCUS record so multi-locus .gbff files load as multi-contig.
 */

const GFF_TYPES = new Set(["gene", "CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "mobile_genetic_element", "region"]);

export function parseGff3(text, fallbackSeqId = "sequence") {
  const bySeq = new Map();
  let count = 0;

  for (const rawLine of text.split(/\r?\n/)) {
    if (!rawLine || rawLine.startsWith("#")) continue;
    const cols = rawLine.split("\t");
    if (cols.length < 9) continue;
    const [seqidRaw, , typeRaw, s, e, , strandRaw, phaseRaw, attrRaw] = cols;
    const type = typeRaw.trim();
    // Keep the feature classes a bacterial genome browser cares about;
    // everything else is noise for this view.
    if (!GFF_TYPES.has(type) && !/^(pseudogene|exon|five_prime_UTR|three_prime_UTR)$/.test(type)) continue;

    const start1 = parseInt(s, 10);
    const end1 = parseInt(e, 10);
    if (!Number.isFinite(start1) || !Number.isFinite(end1) || end1 < start1) continue;

    const attrs = {};
    for (const pair of attrRaw.split(";")) {
      const i = pair.indexOf("=");
      if (i < 0) continue;
      const k = pair.slice(0, i).trim();
      const v = pair.slice(i + 1).trim().replace(/%3D/gi, "=").replace(/%2C/gi, ",").replace(/%3B/gi, ";");
      if (!(k in attrs)) attrs[k] = v;
    }

    const seqid = seqidRaw.trim() || fallbackSeqId;
    if (!bySeq.has(seqid)) bySeq.set(seqid, []);
    bySeq.get(seqid).push({
      contigId: seqid,
      type: type === "gene" && attrs.product ? "CDS" : type,
      start: start1 - 1,
      end: end1,
      strand: strandRaw.trim() === "-" ? "-" : "+",
      phase: phaseRaw.trim(),
      locusTag: attrs.locus_tag || attrs.ID || attrs.Name || `feat_${++count}`,
      product: attrs.product || attrs.Note || attrs.gene || "",
    });
    count++;
  }
  return { featuresByContig: bySeq, numFeatures: count };
}

/* --------------------------------- GenBank --------------------------------- */

/**
 * Parse one or more LOCUS records. Returns { records:[{id, desc, seq}],
 * featuresByContig: Map }. Only the feature kinds above are kept; qualifiers
 * map onto the same display fields as the GFF path.
 */
export function parseGenbank(text) {
  const lines = text.split(/\r?\n/);
  const records = [];
  const featuresByContig = new Map();

  // Locate record boundaries first (multi-locus files).
  const starts = [];
  lines.forEach((l, i) => { if (/^LOCUS\s/.test(l)) starts.push(i); });
  if (!starts.length) throw new Error("No LOCUS line found - is this a GenBank file?");
  starts.push(lines.length);

  for (let r = 0; r < starts.length - 1; r++) {
    const recLines = lines.slice(starts[r], starts[r + 1]);
    const locus = recLines[0].replace(/^LOCUS\s+/, "").trim().split(/\s+/)[0];
    const definition = recLines.find((l) => /^DEFINITION/.test(l))?.replace(/^DEFINITION\s+/, "").trim() ?? "";
    const seqParts = [];
    const feats = [];

    // Feature table spans from "FEATURES" to "ORIGIN" (or end).
    let i = recLines.findIndex((l) => /^FEATURES/.test(l));
    const originIdx = recLines.findIndex((l) => /^ORIGIN/.test(l));
    const featEnd = originIdx >= 0 ? originIdx : recLines.length;

    while (i >= 0 && i < featEnd) {
      const l = recLines[i];
      const m = /^ {5}(\S+)\s+(.+?)\s*$/.exec(l); // top-level feature key at col 6
      if (m) {
        const type = m[1];
        if (GFF_TYPES.has(type)) {
          // Location may wrap onto numeric continuation lines.
          let loc = m[2];
          let j = i + 1;
          while (j < featEnd && /^\s+\d/.test(recLines[j]) && !/\//.test(recLines[j])) {
            loc += recLines[j].trim();
            j++;
          }
          const strand = /complement/.test(loc) ? "-" : "+";
          // Strip complement()/join()/<>/ etc. down to numeric spans and take
          // the outermost bounds (good enough for a browser view).
          const nums = [...loc.matchAll(/(\d+)/g)].map((mm) => parseInt(mm[1], 10)).filter(Number.isFinite);
          if (nums.length >= 2) {
            const start = Math.min(...nums);
            const end = Math.max(...nums);
            // Collect qualifiers until next top-level feature.
            const quals = {};
            for (let q = i + 1; q < featEnd; q++) {
              const ql = recLines[q];
              if (/^ {5}\S/.test(ql)) break;
              const qm = /^\s*\/([^=]+)(?:=(.*))?$/.exec(ql);
              if (qm) {
                let key = qm[1].trim();
                let val = (qm[2] ?? "").replace(/^"|"$/g, "");
                // Continuation lines for this qualifier are indented further
                // and not themselves qualifiers.
                let k = q + 1;
                while (k < featEnd && !/^\s*\//.test(recLines[k]) && !/^ {5}\S/.test(recLines[k])) {
                  val += recLines[k].trim();
                  k++;
                }
                if (!(key in quals)) quals[key] = val;
              }
            }
            feats.push({
              contigId: locus,
              type: type === "gene" ? "CDS" : type,
              start: start - 1,
              end,
              strand,
              phase: ".",
              locusTag: quals.locus_tag || `feat_${feats.length + 1}`,
              product: quals.product || quals.gene || "",
            });
          }
        }
      }
      i++;
    }

    if (originIdx >= 0) {
      for (let s = originIdx + 1; s < recLines.length; s++) {
        if (/^\/\//.test(recLines[s])) break;
        seqParts.push(recLines[s].replace(/[^A-Za-z]/g, "").toUpperCase());
      }
    }

    const seq = seqParts.join("");
    records.push({ id: locus, desc: definition, seq });
    if (feats.length) featuresByContig.set(locus, feats);
  }

  return { records, featuresByContig };
}
