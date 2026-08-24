/**
 * Produces a synthetic Nanopore-like FASTQ *File* object (not a pre-parsed
 * dataset). The point is that "load sample data" exercises exactly the
 * same worker/parser code path as a real upload, instead of maintaining a
 * second, divergent data shape just for the demo.
 */
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const BASES = ["A", "C", "G", "T"];
function randSeq(len, rng, gc = 0.507) {
  let s = "";
  for (let i = 0; i < len; i++) {
    const r = rng();
    s += r < gc / 2 ? "G" : r < gc ? "C" : r < gc + (1 - gc) / 2 ? "A" : "T";
  }
  return s;
}
function gaussian(rng, mean, sd) {
  const u = Math.max(rng(), 1e-9), v = Math.max(rng(), 1e-9);
  const z = Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  return mean + z * sd;
}
function qChar(q) { return String.fromCharCode(33 + Math.max(2, Math.min(60, Math.round(q)))); }
function mutate(rng, b) { const o = BASES.filter((x) => x !== b); return o[Math.floor(rng() * o.length)]; }

export function generateSampleFastqFile({ seed = 42, readCount = 500, genomeLength = 20000, fileName = "sample_nanopore_reads.fastq" } = {}) {
  // Safety cap: this is a demo-data generator, not the real-file pipeline (which
  // streams and scales fine - see fastqParser.worker.js). Building output as one
  // giant retained JS string array does NOT scale past tens of thousands of reads
  // (GC pressure grows superlinearly with total live string data), so we both cap
  // and stream the construction in bounded batches to stay safe at any size.
  if (readCount > 20000) {
    console.warn(`generateSampleFastqFile: readCount ${readCount} capped to 20000 to avoid excessive memory use.`);
    readCount = 20000;
  }

  const refRng = mulberry32(seed);
  const reference = randSeq(genomeLength, refRng, 0.507);
  const rng = mulberry32(seed + 1000);

  const BATCH_SIZE = 2000;
  const blobParts = [];
  let batch = [];

  for (let i = 0; i < readCount; i++) {
    let len = Math.max(180, Math.min(genomeLength - 5, Math.round(gaussian(rng, 1450, 620))));
    const start = Math.floor(rng() * Math.max(1, genomeLength - len));
    const trueSeq = reference.slice(start, start + len);
    const errRate = 0.06 + rng() * 0.05;
    let seq = "", qual = "";
    for (let p = 0; p < trueSeq.length; p++) {
      const edgeFactor = Math.min(p, trueSeq.length - p) < trueSeq.length * 0.06 ? 1.6 : 1;
      let base = trueSeq[p];
      let q = gaussian(rng, 20, 5) / edgeFactor;
      if (rng() < errRate) {
        const roll = rng();
        if (roll < 0.7) base = mutate(rng, base);
        else if (roll < 0.85) { qual += qChar(Math.max(2, q - 10)); continue; }
        else { seq += mutate(rng, base); qual += qChar(Math.max(2, q - 8)); }
        q -= 9;
      }
      seq += base;
      qual += qChar(Math.max(2, q));
    }
    qual = qual.slice(0, seq.length).padEnd(seq.length, "#");
    const runId = "a1b2c3d4e5f6a7b8c9d0e1f2a3b4c5d6";
    const header = `read_${String(i + 1).padStart(4, "0")} runid=${runId} ch=${1 + (i % 512)} start_time=2026-08-1${(i % 9)}T0${(i % 9)}:00:00Z`;
    batch.push(`@${header}\n${seq}\n+\n${qual}\n`);

    if (batch.length >= BATCH_SIZE) {
      blobParts.push(new Blob(batch, { type: "text/plain" }));
      batch = []; // drop the batch's string references so GC can reclaim them now
    }
  }
  if (batch.length) blobParts.push(new Blob(batch, { type: "text/plain" }));

  const blob = new Blob(blobParts, { type: "text/plain" });
  return new File([blob], fileName, { type: "text/plain" });
}
