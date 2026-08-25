/**
 * Parsers for pan-GWAS inputs.
 *
 * Presence/absence matrix — accepts the standard outputs of ROARY/Panaroo/
 * Panacota pipelines plus any generic 0/1 table:
 *   - ".Rtab" style: TAB-separated, header "Gene<TAB>sampleA<TAB>sampleB...",
 *     one row per gene with 0/1 cells.
 *   - "gene_presence_absence.csv" style: COMMA-separated, leading annotation
 *     columns ("Gene", optionally "Non-unique Gene name", "Annotation"),
 *     then one column per sample whose cells hold locus tags when present
 *     and are empty when absent.
 *   - Any other delimited 0/1 / yes-no matrix with a header row.
 *
 * Phenotype metadata — delimited table with a header row; one row per
 * sample; columns may be binary categories, numbers, or labels (usable as
 * CMH strata such as clade/lineage/source).
 */

function sniffDelimiter(line) {
  const counts = { "\t": 0, ",": 0, ";": 0 };
  for (const ch of line) if (counts[ch] !== undefined) counts[ch]++;
  return Object.entries(counts).sort((a, b) => b[1] - a[1])[0][0] || "\t";
}

function splitDelim(line, delim) {
  const out = [];
  let cur = "", inQ = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i];
    if (ch === '"') { inQ = !inQ; continue; }
    if (ch === delim && !inQ) { out.push(cur.trim()); cur = ""; continue; }
    cur += ch;
  }
  out.push(cur.trim());
  return out;
}

const isBinaryCell = (v) => v === "" || /^[01]$/.test(v) || /^(true|false|yes|no|present|absent)$/i.test(v);
const binaryToOne = (v) => {
  if (v === "") return 0;
  if (v === "1" || /^true|yes|present$/i.test(v)) return 1;
  return 0;
};

/**
 * Parse a presence/absence matrix. Returns
 * { genes: [{name, values:Uint8Array}], samples:[string], format:string }
 */
export function parseMatrix(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0);
  if (lines.length < 2) throw new Error("Matrix needs a header row plus at least one data row.");
  const delim = sniffDelimiter(lines[0]);
  const header = splitDelim(lines[0], delim).map((h) => h.replace(/^"|"$/g, ""));

  // Locate the feature-name column (ROARY files call it exactly "Gene").
  let nameCol = header.findIndex((h) => h.toLowerCase() === "gene");
  if (nameCol < 0) nameCol = 0;
  const sampleCols = [];
  for (let i = 0; i < header.length; i++) {
    if (i === nameCol) continue;
    // Skip ROARY's optional annotation columns.
    if (i !== nameCol && /^(non-unique gene name|annotation|avg sequences? per gene|min|max)/i.test(header[i])) continue;
    sampleCols.push(i);
  }
  if (sampleCols.length < 2) throw new Error("Could not find at least two sample columns in the matrix header.");

  // Decide cell semantics from the first data rows: strict 0/1 (Rtab-ish)
  // vs locus-tag-or-empty (gene_presence_absence.csv-ish). Genuine Rtab
  // files never contain anything but 0/1, so a single foreign string
  // switches us to csv mode (presence = non-empty cell).
  let tagCells = 0, rowsPeeked = 0;
  outer: for (let li = 1; li < lines.length && rowsPeeked < 20; li++, rowsPeeked++) {
    const cells = splitDelim(lines[li], delim);
    for (const ci of sampleCols) {
      const v = cells[ci] ?? "";
      if (!isBinaryCell(v) && v.length > 0) { tagCells++; break outer; }
    }
  }
  const tagMode = tagCells > 0;

  const samples = sampleCols.map((i) => header[i]);
  const genes = [];
  const seenNames = new Set();

  for (let li = 1; li < lines.length; li++) {
    const cells = splitDelim(lines[li], delim);
    const name = (cells[nameCol] ?? "").replace(/^"|"$/g, "");
    if (!name || name.toLowerCase() === "gene") continue;
    const values = new Uint8Array(samples.length);
    let present = 0;
    for (let s = 0; s < sampleCols.length; s++) {
      const raw = (cells[sampleCols[s]] ?? "");
      const v = tagMode ? (raw ? 1 : 0) : binaryToOne(raw);
      if (v) { values[s] = 1; present++; }
    }
    // Constant features (present in every or no sample) are retained; their
    // Fisher margins are degenerate so they simply return p = 1 downstream,
    // and the summary exposes how many of them there are.
    if (!seenNames.has(name)) {
      seenNames.add(name);
      genes.push({ name, values, present });
    }
  }

  if (!genes.length) throw new Error("No gene rows could be parsed from the matrix.");

  return {
    genes,
    samples,
    variableGenes: genes.filter((g) => g.present > 0 && g.present < samples.length).length,
    format: tagMode ? "gene_presence_absence.csv (locus-tag cells)" : "presence/absence 0-1 matrix",
  };
}

/**
 * Parse a phenotype/metadata table. Returns
 * { samples:[string], columns:[{name, kind:'binary'|'numeric'|'text', values:[]}] }
 */
export function parsePhenotypes(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0);
  if (lines.length < 2) throw new Error("Metadata needs a header row plus at least one data row.");
  const delim = sniffDelimiter(lines[0]);
  const header = splitDelim(lines[0], delim).map((h) => h.replace(/^"|"$/g, ""));
  if (header.length < 2) throw new Error("Metadata table needs an ID column and at least one trait column.");

  const rows = lines.slice(1).map((l) => splitDelim(l, delim));
  const nCols = header.length;
  const cols = [];

  for (let c = 0; c < nCols; c++) {
    const rawVals = rows.map((r) => (r[c] ?? "").replace(/^"|"$/g, "")).filter((v) => v !== "");
    const numericVals = rawVals.map((v) => Number(v)).filter(Number.isFinite);
    const uniq = [...new Set(rawVals)];
    let kind, values;
    if (uniq.length === 2) {
      kind = "binary";
      values = rows.map((r) => (r[c] ?? "").replace(/^"|"$/g, ""));
    } else if (rawVals.length >= 2 && numericVals.length === rawVals.length) {
      kind = "numeric";
      values = rows.map((r) => {
        const v = Number((r[c] ?? "").replace(/^"|"$/g, ""));
        return Number.isFinite(v) ? v : null;
      });
    } else {
      kind = uniq.length < 2 ? "constant" : "text";
      values = rows.map((r) => (r[c] ?? "").replace(/^"|"$/g, ""));
    }
    cols.push({ name: header[c], kind, values });
  }

  const samples = rows.map((r, i) => ((r[0] ?? "").trim() || `row_${i + 1}`));

  // Prefer a clearly-named ID column for display purposes only; we keep the
  // first column convention but upgrade the sample ids if a better key exists.
  const idIdx = header.findIndex((h) => /^(id|sample|isolate|strain|name)$/i.test(h));
  if (idIdx > 0) {
    for (let i = 0; i < samples.length; i++) {
      const better = (rows[i][idIdx] ?? "").trim();
      if (better) samples[i] = better;
    }
  }

  return { samples, columns: cols.filter((c) => c.kind !== "constant") };
}

/**
 * Intersect matrix samples with metadata samples. Returns aligned arrays +
 * diagnostics about dropped/unmatched ids.
 */
export function matchSamples(matrixSamples, meta) {
  const metaIndex = new Map(meta.samples.map((s, i) => [s, i]));
  const kept = [], missingMeta = [];
  matrixSamples.forEach((s, i) => {
    if (metaIndex.has(s)) kept.push(i);
    else missingMeta.push(s);
  });
  const unmatchedMeta = meta.samples.filter((s) => !matrixSamples.includes(s));
  const samples = kept.map((i) => matrixSamples[i]);
  return { keptIdx: kept, samples, missingMeta, unmatchedMeta };
}
