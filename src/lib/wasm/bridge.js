/*
 * Loader + binary plumbing for the WASM kernels (src/wasm/genome_core.js).
 *
 * The bundle is a single self-contained ES module (WASM embedded as base64),
 * built by scripts/build-wasm.sh from src/wasm/src/*.c. It exposes three
 * engines - assembly, MSA and phylo - each following the same protocol:
 *
 *   ..._run / ..._init+bootstrap_chunk+finish   compute into an internal blob
 *   ..._result_ptr() / ..._result_len()          locate the blob on the heap
 *   ..._free() / ..._abort()                     release everything
 *
 * Blob layouts are documented next to each decoder in src/lib/engine/*.js.
 * Every engine falls back to the original pure-JS implementation whenever
 * this module cannot be instantiated, so the app keeps working everywhere.
 */

let modPromise = null;

/** Set true to force the pure-JS fallback everywhere (tests / debugging). */
export const gxConfig = { disableWasm: false };

/** Resolves the Emscripten module (singleton). Throws if unavailable. */
export function gxWasm() {
  if (gxConfig.disableWasm) throw new Error("WASM disabled by configuration");
  if (!modPromise) {
    modPromise = import("../../wasm/genome_core.js")
      .then((m) => m.default())
      .catch((err) => {
        modPromise = null;
        throw err;
      });
  }
  return modPromise;
}

/** True when the WASM kernels are usable in this environment. */
export async function hasGxWasm() {
  try {
    await gxWasm();
    return true;
  } catch {
    return false;
  }
}

export function mallocIntoHeap(mod, u8) {
  const bytes = u8 instanceof Uint8Array ? u8 : new Uint8Array(u8.buffer, u8.byteOffset, u8.byteLength);
  const ptr = mod._malloc(bytes.length);
  if (!ptr) throw new Error("Out of WebAssembly memory while packing input data.");
  mod.HEAPU8.set(bytes, ptr);
  return ptr;
}

/** Copies a result blob out of the heap before freeing it. */
export function takeBlob(mod, ptr, len) {
  return mod.HEAPU8.slice(ptr, ptr + len);
}

/**
 * Runs `fn` with globalThis.__gxWasmProgress bound to `onProgress(code, pct)`.
 * The C kernels call gx_progress() at the same checkpoints the JS engines do.
 */
export function withProgress(fn, onProgress) {
  const prev = globalThis.__gxWasmProgress;
  globalThis.__gxWasmProgress = onProgress;
  return Promise.resolve()
    .then(fn)
    .finally(() => {
      globalThis.__gxWasmProgress = prev ?? null;
    });
}

export class BlobReader {
  constructor(buf) {
    this.u8 = buf;
    this.view = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
    this.pos = 0;
  }
  magic(expected) {
    const got = String.fromCharCode(this.u8[0], this.u8[1], this.u8[2], this.u8[3]);
    if (got !== expected) throw new Error(`Corrupt WASM result blob (expected ${expected}, got ${got})`);
    this.pos = 4;
  }
  u32() { const v = this.view.getUint32(this.pos, true); this.pos += 4; return v; }
  i32() { const v = this.view.getInt32(this.pos, true); this.pos += 4; return v; }
  f64() { const v = this.view.getFloat64(this.pos, true); this.pos += 8; return v; }
  readU8() { const v = this.u8[this.pos]; this.pos += 1; return v; }
  skip(n) { this.pos += n; }
  align(a) { const r = this.pos % a; if (r) this.pos += a - r; }
  bytes(n) { const v = this.u8.subarray(this.pos, this.pos + n); this.pos += n; return v; }
}
