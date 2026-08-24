# Genome Explorer

In-browser genomics workbench: FASTQ QC, genome assembly, multiple sequence
alignment, phylogenetic inference, annotation, and association mapping.

## WASM compute kernels (Assembly / Alignment / Phylo)

The three compute-heavy pipelines run as **C compiled to WebAssembly**:

| Section   | C kernel                | What it does                                              |
|-----------|-------------------------|-----------------------------------------------------------|
| Assembly  | `src/wasm/src/assembly.c` | minimizer index → overlap chains → string graph → contigs |
| Alignment | `src/wasm/src/msa.c`      | k-mer distances → UPGMA guide tree → progressive profile DP |
| Phylo     | `src/wasm/src/phylo.c`    | p/JC69/K2P distances → NJ/UPGMA → chunked bootstrap        |

- The kernels are bit-compatible ports of the original JS engines
  (`src/lib/{assembly,msa,phylo}`), including tie-breaking order. A parity
  harness (`npm run parity`) asserts both engines produce identical output.
- Every worker uses WASM when available and **falls back to the pure-JS
  engine automatically** (old browsers, WASM failure, or
  `gxConfig.disableWasm = true` in `src/lib/wasm/bridge.js`).
- The C engines raise the profile-DP capacity from 160M to 900M cells and use
  explicit memory (growable to ~2 GB), so large / real datasets no longer hit
  JS GC limits.

Measured speedups (node, single-threaded, see `npm run bench`):
assembly ~2–3x, phylo ~1.5–2x, MSA ~1.3x (the affine-gap DP is branchy and
V8 already JITs it well; the win there is capacity + stability).

### Rebuilding the WASM bundle

The compiled bundle is committed at `src/wasm/genome_core.js` (a
self-contained ES module with the WASM embedded), so `npm run dev` / `build`
work without any native toolchain. To change the C sources:

```bash
# one-time: install Emscripten (no sudo required)
git clone --depth 1 https://github.com/emscripten-core/emsdk.git ~/.local/emsdk
~/.local/emsdk/emsdk install latest && ~/.local/emsdk/emsdk activate latest

# rebuild after editing src/wasm/src/*.c
npm run build:wasm

# optional debug build with assertions + SAFE_HEAP
bash scripts/build-debug.sh
```

### Verification

```bash
npm run parity   # WASM vs JS engines must produce identical results
npm run bench    # side-by-side timings
```

## Development

```bash
npm install
npm run dev      # start dev server
npm run build    # production build
npm run lint
```
