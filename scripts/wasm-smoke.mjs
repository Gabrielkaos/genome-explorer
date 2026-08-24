// Smoke test: instantiate the WASM module under Node and run tiny jobs.
import init from "../src/wasm/genome_core.js";

const M = await init();
console.log("wasm ready");

const enc = new TextEncoder();

function packReads(seqs) {
  const total = seqs.reduce((a, s) => a + s.length, 0);
  const flat = new Uint8Array(total);
  const offs = new Uint32Array(seqs.length + 1);
  let p = 0;
  seqs.forEach((s, i) => {
    flat.set(enc.encode(s), p);
    p += s.length;
    offs[i + 1] = p;
  });
  return { flat, offs };
}

function pushHeap(u8) {
  const ptr = M._malloc(u8.byteLength);
  if (!ptr) throw new Error("OOM");
  const bytes = u8 instanceof Uint8Array ? u8 : new Uint8Array(u8.buffer, u8.byteOffset, u8.byteLength);
  M.HEAPU8.set(bytes, ptr);
  return ptr;
}

// ---- assembly smoke ----
{
  // two overlapping reads
  const g = "ACGT".repeat(200);
  const reads = [g, g.slice(50) + "TTTT"];
  const { flat, offs } = packReads(reads);
  globalThis.__gxWasmProgress = (code, pct) => console.log("progress", code, pct);
  const fp = pushHeap(flat), op = pushHeap(offs);
  const rc = M._gx_assembly_run(fp, op, reads.length, 10, 5, 40, 60, 4, 200);
  const ptr = M._gx_assembly_result_ptr(), len = M._gx_assembly_result_len();
  console.log("assembly rc", rc, "blob bytes", len);
  const blob = M.HEAPU8.slice(ptr, ptr + len);
  console.log("magic", String.fromCharCode(...blob.slice(0, 4)));
  M._gx_assembly_free(); M._free(fp); M._free(op);
}

// ---- msa smoke ----
{
  const seqs = ["ACGTACGTAC", "ACGTTCGTAC", "ACGAACGT--"];
  const LOOKUP = new Int8Array(128).fill(4);
  LOOKUP[65] = 0; LOOKUP[67] = 1; LOOKUP[71] = 2; LOOKUP[84] = 3; LOOKUP[45] = 4;
  const enc2 = seqs.map((s) => Uint8Array.from(s, (c) => LOOKUP[c.charCodeAt(0) & 127]));
  const flat = new Uint8Array(enc2.reduce((a, e) => a + e.length, 0));
  const lens = new Uint32Array(enc2.map((e) => e.length));
  let p = 0;
  enc2.forEach((e) => { flat.set(e, p); p += e.length; });
  const fp = pushHeap(flat), lp = pushHeap(lens);
  const rc = M._gx_msa_run(fp, lp, seqs.length, 3, 2, 1, 8, 1, 0.5);
  const ptr = M._gx_msa_result_ptr(), len = M._gx_msa_result_len();
  const blob = M.HEAPU8.slice(ptr, ptr + len);
  console.log("msa rc", rc, "bytes", len, "magic", String.fromCharCode(...blob.slice(0, 4)));
  M._gx_msa_free(); M._free(fp); M._free(lp);
}

// ---- phylo smoke ----
{
  const n = 5, L = 100;
  const rows = [];
  for (let i = 0; i < n; i++) {
    let s = "";
    for (let c = 0; c < L; c++) s += "ACGT"[(i * 7 + c * 3) % 4];
    rows.push(s);
  }
  const flat = new Uint8Array(n * L);
  rows.forEach((r, i) => flat.set(enc.encode(r), i * L));
  const fp = pushHeap(flat);
  M._gx_phylo_init(fp, n, L, 1 /*jc*/, 0, 0 /*nj*/, 5, 42);
  while (M._gx_phylo_reps_done() < 5) M._gx_phylo_bootstrap_chunk(2);
  M._gx_phylo_finish();
  const ptr = M._gx_phylo_result_ptr(), len = M._gx_phylo_result_len();
  const blob = M.HEAPU8.slice(ptr, ptr + len);
  console.log("phylo ok bytes", len, "magic", String.fromCharCode(...blob.slice(0, 4)));
  M._gx_phylo_abort(); M._free(fp);
}
console.log("SMOKE OK");
