/**
 * FASTQ parsing worker.
 *
 * Runs entirely off the main thread so large files (hundreds of MB - GB)
 * don't freeze the UI. Handles both plain .fastq and gzip-compressed
 * .fastq.gz input.
 *
 * Random-access strategy:
 *  - Plain files: decompressed byte offsets == original file byte offsets,
 *    so we record {start,len} per record and the main thread can later
 *    File.slice(start, start+len) to re-fetch any single record's exact
 *    text on demand. This scales to arbitrarily large files.
 *  - Gzip files: decompressed offsets do NOT map back to compressed file
 *    offsets, so true random access isn't cheaply available in-browser.
 *    We instead keep full sequence/quality text for a bounded "preview"
 *    window (the first PREVIEW_CAP records) and compute exact aggregate
 *    stats (length/quality/GC/N50/etc.) for the *entire* file regardless.
 *    This limitation is reported back to the UI explicitly.
 */
import { Inflate } from "pako";
import { createLineScanner, bytesToString } from "../lib/fastq/lineScanner.js";
import { meanQualityFromBytes } from "../lib/fastq/phred.js";
import {
  computeNStat, median, buildLengthBinEdges, buildQualityBinEdges,
  createHistogramAccumulator, gcPercent,
} from "../lib/fastq/stats.js";

const PREVIEW_CAP = 5000;
const PROGRESS_EVERY_BYTES = 4_000_000;

self.onmessage = async (e) => {
  if (e.data?.type !== "parse") return;
  const { file } = e.data;
  try {
    await parseFile(file);
  } catch (err) {
    self.postMessage({ type: "error", message: String(err?.message || err) });
  }
};

async function parseFile(file) {
  const t0 = performance.now();
  const isGzip = /\.gz$/i.test(file.name);
  const parseMode = isGzip ? "buffered" : "offset";

  const warnings = [];
  let totalReads = 0;
  let totalBases = 0;
  let totalErrorProbSum = 0; // for correct base-weighted dataset mean quality
  let minLength = Infinity, maxLength = 0;
  const lengths = []; // Int32 later
  const meanQs = []; // Float32 later
  const gcPerRead = []; // Float32 later — per-read GC%, cheap to compute since we already loop bases
  const ids = [];
  const offsets = parseMode === "offset" ? [] : null; // flat [start0,len0,start1,len1,...]
  const previewRecords = []; // {seq, qual} for first PREVIEW_CAP reads
  const baseCounts = { A: 0, C: 0, G: 0, T: 0, N: 0 };

  // --- FASTQ 4-line record state machine, driven by the byte-level line scanner ---
  let lineIdx = 0; // 0=header 1=seq 2=plus 3=qual
  let curId = "", curDesc = "";
  let curSeqBytes = null;
  let recordStartOffset = 0;
  let malformedCount = 0;

  function handleLine(bytes, startOffset, endOffset) {
    if (lineIdx === 0) {
      if (bytes.length === 0) return; // tolerate stray blank lines between records
      if (bytes[0] !== 0x40 /* '@' */) {
        malformedCount++;
        if (warnings.length < 20) warnings.push(`Line at byte ${startOffset}: expected '@' header, found something else - skipping to resync.`);
        return; // stay at lineIdx 0, try next line as a header
      }
      recordStartOffset = startOffset;
      const headerText = bytesToString(bytes.subarray(1));
      const spaceIdx = headerText.indexOf(" ");
      curId = spaceIdx === -1 ? headerText : headerText.slice(0, spaceIdx);
      curDesc = spaceIdx === -1 ? "" : headerText.slice(spaceIdx + 1);
      lineIdx = 1;
    } else if (lineIdx === 1) {
      curSeqBytes = bytes.slice(); // copy, since underlying buffer gets reused
      lineIdx = 2;
    } else if (lineIdx === 2) {
      // '+' separator line, content ignored
      lineIdx = 3;
    } else if (lineIdx === 3) {
      const qualBytes = bytes;
      if (qualBytes.length !== curSeqBytes.length) {
        malformedCount++;
        if (warnings.length < 20) {
          warnings.push(`Read "${curId}": sequence length (${curSeqBytes.length}) != quality length (${qualBytes.length}) - record skipped.`);
        }
        lineIdx = 0;
        return;
      }
      const len = curSeqBytes.length;
      const { meanQ, sumErrorProb } = meanQualityFromBytes(qualBytes);

      // base composition (both dataset-wide and per-read, since we're looping anyway)
      let ra = 0, rc = 0, rg = 0, rt = 0;
      for (let i = 0; i < len; i++) {
        const b = curSeqBytes[i];
        if (b === 65 || b === 97) { baseCounts.A++; ra++; }       // A/a
        else if (b === 67 || b === 99) { baseCounts.C++; rc++; }  // C/c
        else if (b === 71 || b === 103) { baseCounts.G++; rg++; } // G/g
        else if (b === 84 || b === 116) { baseCounts.T++; rt++; } // T/t
        else baseCounts.N++;
      }
      const racgt = ra + rc + rg + rt;
      gcPerRead.push(racgt > 0 ? ((rg + rc) / racgt) * 100 : 0);

      totalReads++;
      totalBases += len;
      totalErrorProbSum += sumErrorProb;
      if (len < minLength) minLength = len;
      if (len > maxLength) maxLength = len;
      lengths.push(len);
      meanQs.push(meanQ);
      ids.push(curId);

      if (parseMode === "offset") {
        const recLen = endOffset - recordStartOffset;
        offsets.push(recordStartOffset, recLen);
      }
      if (previewRecords.length < PREVIEW_CAP) {
        previewRecords.push({ seq: bytesToString(curSeqBytes), qual: bytesToString(qualBytes), id: curId, desc: curDesc });
      }
      lineIdx = 0;
    }
  }

  const scanner = createLineScanner(handleLine);

  const totalBytes = file.size;
  let bytesForProgress = 0; // counts compressed bytes read (progress proxy)
  let lastProgressPost = 0;

  function postProgress(force = false) {
    if (force || bytesForProgress - lastProgressPost > PROGRESS_EVERY_BYTES) {
      lastProgressPost = bytesForProgress;
      self.postMessage({ type: "progress", bytesProcessed: bytesForProgress, totalBytes, readsProcessed: totalReads });
    }
  }

  if (isGzip) {
    await new Promise((resolve, reject) => {
      const inflator = new Inflate();
      inflator.onData = (chunk) => scanner.push(chunk);
      inflator.onEnd = (status) => {
        if (status !== 0) reject(new Error(`gzip decompression failed (status ${status})`));
        else resolve();
      };
      const reader = file.stream().getReader();
      function pump() {
        reader.read().then(({ done, value }) => {
          if (done) { inflator.push(new Uint8Array(0), true); return; }
          bytesForProgress += value.byteLength;
          postProgress();
          inflator.push(value, false);
          pump();
        }).catch(reject);
      }
      pump();
    });
  } else {
    const reader = file.stream().getReader();
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      bytesForProgress += value.byteLength;
      postProgress();
      scanner.push(value);
    }
  }
  scanner.flush();
  postProgress(true);

  if (totalReads === 0) {
    self.postMessage({ type: "error", message: "No valid FASTQ records were found in this file. Check that it's a standard 4-line-per-record FASTQ file (and that a .gz file is really gzip-compressed)." });
    return;
  }
  if (malformedCount > 0) {
    warnings.unshift(`${malformedCount} malformed record(s) were skipped during parsing.`);
  }

  // ---- aggregate stats ----
  const lengthEdges = buildLengthBinEdges(maxLength);
  const qualEdges = buildQualityBinEdges(Math.max(...meanQs, 30));
  const lenHist = createHistogramAccumulator(lengthEdges);
  const qHist = createHistogramAccumulator(qualEdges);
  for (let i = 0; i < lengths.length; i++) { lenHist.add(lengths[i]); qHist.add(meanQs[i]); }

  const n50 = computeNStat(lengths, 0.5);
  const n90 = computeNStat(lengths, 0.9);
  const medianLength = median(lengths);
  const meanLength = totalBases / totalReads;
  const datasetMeanQErrProb = totalErrorProbSum / totalBases;
  const datasetMeanQ = -10 * Math.log10(Math.max(datasetMeanQErrProb, 1e-12));

  const lengthsTyped = Int32Array.from(lengths);
  const meanQsTyped = Float32Array.from(meanQs);
  const gcPerReadTyped = Float32Array.from(gcPerRead);
  const offsetsTyped = offsets ? Uint32Array.from(offsets) : null;

  const dataset = {
    fileName: file.name,
    fileSize: file.size,
    compressed: isGzip,
    parseMode,
    previewCap: PREVIEW_CAP,
    totalReads,
    totalBases,
    minLength, maxLength, meanLength, medianLength, n50, n90,
    datasetMeanQ,
    baseCounts,
    gcPercent: gcPercent(baseCounts),
    lengthHistogram: lenHist.toArray(),
    qualityHistogram: qHist.toArray(),
    warnings,
    parseTimeMs: performance.now() - t0,
  };

  const transferables = [lengthsTyped.buffer, meanQsTyped.buffer, gcPerReadTyped.buffer];
  if (offsetsTyped) transferables.push(offsetsTyped.buffer);

  self.postMessage({
    type: "done",
    dataset,
    lengths: lengthsTyped,
    meanQs: meanQsTyped,
    gcPerRead: gcPerReadTyped,
    ids,
    offsets: offsetsTyped,
    previewRecords,
  }, transferables);
}
