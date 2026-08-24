/**
 * Builds a filtered FASTQ file for download. For offset-indexed (plain)
 * files this is essentially free: Blob parts can be File.slice() ranges,
 * which the browser handles lazily without materializing the whole file
 * in JS memory, even for files far larger than available RAM would allow
 * as a string. Gzip input is exported uncompressed (re-gzipping large
 * output client-side is possible but out of scope for now).
 */
export function buildFilteredFastqBlob({ file, dataset, index, previewRecords, passingIndices }) {
  if (dataset.parseMode === "offset" && index?.offsets) {
    const parts = [];
    for (const i of passingIndices) {
      const start = index.offsets[i * 2];
      const len = index.offsets[i * 2 + 1];
      parts.push(file.slice(start, start + len));
    }
    return new Blob(parts, { type: "text/plain" });
  }

  // buffered / gzip mode: only records within the in-memory preview window
  // have full text available.
  const lines = [];
  let skipped = 0;
  for (const i of passingIndices) {
    const rec = previewRecords?.[i];
    if (!rec) { skipped++; continue; }
    const header = rec.desc ? `${rec.id} ${rec.desc}` : rec.id;
    lines.push(`@${header}`, rec.seq, "+", rec.qual);
  }
  return { blob: new Blob([lines.join("\n") + "\n"], { type: "text/plain" }), skipped };
}

export function downloadBlob(blobOrResult, fileName) {
  const blob = blobOrResult instanceof Blob ? blobOrResult : blobOrResult.blob;
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = fileName;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}
