/**
 * Annotation parsers for the genome browser: GFF3 (and tolerant GFF2) plus
 * GenBank flat files. Coordinates are normalized to 0-based half-open
 * [start, end) internally, matching the rest of the app.
 *
 * GFF3:  seqid source type start end score strand phase attributes
 * GenBank: FEATURES table CDS/gene/rRNA/tRNA/etc. entries + ORIGIN sequence,
 * parsed per LOCUS record so multi-locus .gbff files load as multi-contig.
 */

const GFF_TYPES = new Set(["gene", "CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "mobile_genetic_element", "pseudogene", "region", "repeat_region", "exon", "mRNA", "transcript"]);

function safeDecode(val) {
  if (!val) return "";
  try {
    return decodeURIComponent(val.replace(/\+/g, "%20"));
  } catch {
    return val.replace(/%3D/gi, "=").replace(/%2C/gi, ",").replace(/%3B/gi, ";").replace(/%20/gi, " ");
  }
}

function parseAttributes(attrRaw) {
  const attrs = {};
  if (!attrRaw || attrRaw === ".") return attrs;
  const parts = attrRaw.split(/;\s*/);
  for (const part of parts) {
    if (!part.trim()) continue;
    const eqIdx = part.indexOf("=");
    if (eqIdx >= 0) {
      const k = part.slice(0, eqIdx).trim();
      const v = safeDecode(part.slice(eqIdx + 1).trim().replace(/^["']|["']$/g, ""));
      if (k && !(k in attrs)) attrs[k] = v;
    } else {
      const m = part.trim().match(/^([A-Za-z0-9_.-]+)\s+["']?(.*?)["']?$/);
      if (m) {
        const k = m[1].trim();
        const v = safeDecode(m[2].trim());
        if (k && !(k in attrs)) attrs[k] = v;
      }
    }
  }
  return attrs;
}

function parseFastaPart(text) {
  const records = [];
  let cur = null;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (cur) records.push(cur);
      const header = line.slice(1).trim();
      cur = { id: header.split(/\s+/)[0] || `seq_${records.length + 1}`, desc: header, seq: "" };
    } else if (cur) cur.seq += line.trim().toUpperCase();
  }
  if (cur) records.push(cur);
  return records.filter((r) => r.seq.length > 0);
}

const NORM_TYPES = {
  cds: "CDS",
  gene: "gene",
  pseudogene: "pseudogene",
  trna: "tRNA",
  rrna: "rRNA",
  tmrna: "tmRNA",
  ncrna: "ncRNA",
  misc_rna: "misc_RNA",
  snrna: "ncRNA",
  snorna: "ncRNA",
  mirna: "ncRNA",
  mrna: "mRNA",
  transcript: "mRNA",
  exon: "exon",
  mobile_genetic_element: "mobile_genetic_element",
  insertion_sequence: "mobile_genetic_element",
  transposon: "mobile_genetic_element",
  prophage: "mobile_genetic_element",
  repeat_region: "repeat_region",
  crispr: "repeat_region",
  region: "region",
  source: "region",
};

export function parseGff3(text, fallbackSeqId = "sequence") {
  const bySeq = new Map();
  const sequenceRegions = new Map();
  let fastaRecords = [];
  let count = 0;

  // 1. Separate annotations and ##FASTA if present
  let annotText = text;
  const fastaMatch = /(?:^|\r?\n)##FASTA\s*(?:\r?\n|$)/i.exec(text);
  if (fastaMatch) {
    annotText = text.slice(0, fastaMatch.index);
    const fastaPart = text.slice(fastaMatch.index + fastaMatch[0].length);
    fastaRecords = parseFastaPart(fastaPart).map((r) => ({
      ...r,
      seq: r.seq.toUpperCase().replace(/[^ACGTN]/gi, "N"),
      circular: false,
    }));
  }

  // Raw features collection for hierarchy resolution
  const rawFeatures = [];

  for (const rawLine of annotText.split(/\r?\n/)) {
    const line = rawLine.trim();
    if (!line) continue;

    if (line.startsWith("#")) {
      // Parse ##sequence-region directives
      const seqRegMatch = /^##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)/i.exec(line);
      if (seqRegMatch) {
        const id = seqRegMatch[1];
        const s = parseInt(seqRegMatch[2], 10);
        const e = parseInt(seqRegMatch[3], 10);
        if (Number.isFinite(s) && Number.isFinite(e) && e >= s) {
          sequenceRegions.set(id, { start: s, end: e, length: e - s + 1 });
        }
      }
      continue;
    }

    const cols = line.split("\t");
    if (cols.length < 9) continue;
    const [seqidRaw, , typeRaw, s, e, , strandRaw, phaseRaw, attrRaw] = cols;
    const rawType = typeRaw.trim();
    const normType = NORM_TYPES[rawType.toLowerCase()] || rawType;

    const start1 = parseInt(s, 10);
    const end1 = parseInt(e, 10);
    if (!Number.isFinite(start1) || !Number.isFinite(end1) || end1 < start1) continue;

    const attrs = parseAttributes(attrRaw);
    const seqid = seqidRaw.trim() || fallbackSeqId;

    if (normType === "region") {
      const isCirc = /true|yes|circular/i.test(attrs.Is_circular || attrs.is_circular || attrs.circular || "");
      const prev = sequenceRegions.get(seqid) || { start: start1, end: end1, length: end1 - start1 + 1 };
      sequenceRegions.set(seqid, {
        ...prev,
        start: Math.min(prev.start, start1),
        end: Math.max(prev.end, end1),
        length: Math.max(prev.length, end1 - start1 + 1),
        circular: isCirc || prev.circular || false,
      });
      continue; // Skip full-contig region features from cluttering the feature tracks
    }

    rawFeatures.push({
      seqid,
      rawType,
      type: normType,
      start: start1 - 1,
      end: end1,
      strand: strandRaw.trim() === "-" ? "-" : "+",
      phase: phaseRaw.trim(),
      attrs,
      id: attrs.ID || attrs.transcript_id || attrs.gene_id || null,
      parent: attrs.Parent || (attrs.transcript_id && normType !== "mRNA" && normType !== "transcript" ? attrs.transcript_id : attrs.gene_id && normType !== "gene" ? attrs.gene_id : null),
      locusTag: attrs.locus_tag || attrs.Locus_tag || attrs.locusTag || attrs.gene_id || attrs.protein_id || null,
      name: attrs.Name || attrs.gene || attrs.Gene || attrs.gene_name || null,
      product: attrs.product || attrs.Product || attrs.description || attrs.Description || attrs.Note || attrs.note || attrs.function || attrs.gene_desc || "",
      isPseudo: /^(true|yes|1|unknown)$/i.test(attrs.pseudo || attrs.pseudogene || "") || normType === "pseudogene",
    });
  }

  // 2. Hierarchy resolution: link CDS/tRNA/rRNA/etc. to parent gene
  const genesById = new Map();
  const genesByLocus = new Map();
  const genesBySpan = new Map();

  for (const f of rawFeatures) {
    if (f.type === "gene") {
      if (f.id) genesById.set(f.id, f);
      if (f.locusTag) genesByLocus.set(f.locusTag, f);
      genesBySpan.set(`${f.seqid}:${f.start}-${f.end}:${f.strand}`, f);
    }
  }

  const parentsWithChildren = new Set();
  for (const f of rawFeatures) {
    if (f.type !== "gene" && f.type !== "region") {
      if (f.parent && genesById.has(f.parent)) {
        parentsWithChildren.add(genesById.get(f.parent));
      } else if (f.locusTag && genesByLocus.has(f.locusTag)) {
        parentsWithChildren.add(genesByLocus.get(f.locusTag));
      } else {
        const spanKey = `${f.seqid}:${f.start}-${f.end}:${f.strand}`;
        if (genesBySpan.has(spanKey)) {
          parentsWithChildren.add(genesBySpan.get(spanKey));
        }
      }
    }
  }

  // 3. Build unified features list
  for (const f of rawFeatures) {
    if (f.type === "gene" && parentsWithChildren.has(f)) {
      // Skip redundant parent gene container when specific child CDS/RNA features exist
      continue;
    }

    let parentGene = null;
    if (f.parent && genesById.has(f.parent)) parentGene = genesById.get(f.parent);
    else if (f.locusTag && genesByLocus.has(f.locusTag)) parentGene = genesByLocus.get(f.locusTag);
    else {
      const spanKey = `${f.seqid}:${f.start}-${f.end}:${f.strand}`;
      if (genesBySpan.has(spanKey)) parentGene = genesBySpan.get(spanKey);
    }

    const locusTag = f.locusTag || (parentGene?.locusTag) || f.name || (parentGene?.name) || f.id || `feat_${++count}`;
    const product = f.product || (parentGene?.product) || (f.name ? `gene ${f.name}` : "") || "";
    const displayType = f.isPseudo ? "pseudogene" : f.type === "gene" ? (product ? "CDS" : "gene") : f.type;

    if (!bySeq.has(f.seqid)) bySeq.set(f.seqid, []);
    bySeq.get(f.seqid).push({
      contigId: f.seqid,
      type: displayType,
      start: f.start,
      end: f.end,
      strand: f.strand,
      phase: f.phase,
      locusTag,
      product,
    });
    count++;
  }

  // Apply circularity to fastaRecords if defined in sequenceRegions
  for (const r of fastaRecords) {
    const reg = sequenceRegions.get(r.id);
    if (reg?.circular) r.circular = true;
  }

  return {
    featuresByContig: bySeq,
    numFeatures: count,
    records: fastaRecords,
    sequenceRegions,
  };
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
