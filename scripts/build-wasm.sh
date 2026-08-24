#!/usr/bin/env bash
# Compiles the C kernels (assembly / MSA / phylo) to a single ES-module WASM
# bundle consumed by src/lib/wasm/bridge.js. Requires Emscripten (emsdk).
#
#   source ~/.local/emsdk/emsdk_env.sh && npm run build:wasm
#
set -euo pipefail
cd "$(dirname "$0")/.."

if ! command -v emcc >/dev/null 2>&1; then
  if [ -f "$HOME/.local/emsdk/emsdk_env.sh" ]; then
    # shellcheck disable=SC1091
    source "$HOME/.local/emsdk/emsdk_env.sh"
  fi
fi
command -v emcc >/dev/null 2>&1 || { echo "error: emcc not found (install emsdk)" >&2; exit 1; }

mkdir -p src/wasm

emcc src/wasm/src/gx_glue.c src/wasm/src/assembly.c src/wasm/src/msa.c src/wasm/src/phylo.c \
  -o src/wasm/genome_core.js \
  -O3 \
  --no-entry \
  -sMODULARIZE=1 \
  -sEXPORT_ES6=1 \
  -sSINGLE_FILE=1 \
  -sENVIRONMENT=web,worker,node \
  -sALLOW_MEMORY_GROWTH=1 \
  -sINITIAL_MEMORY=134217728 \
  -sSTACK_SIZE=4194304 \
  -sABORTING_MALLOC=0 \
  -sEXPORTED_RUNTIME_METHODS=HEAPU8,HEAPU32,HEAPF64 \
  -sEXPORTED_FUNCTIONS=_malloc,_free,\
_gx_assembly_run,_gx_assembly_result_ptr,_gx_assembly_result_len,_gx_assembly_free,\
_gx_msa_run,_gx_msa_result_ptr,_gx_msa_result_len,_gx_msa_free,\
_gx_phylo_init,_gx_phylo_bootstrap_chunk,_gx_phylo_reps_done,_gx_phylo_finish,_gx_phylo_result_ptr,_gx_phylo_result_len,_gx_phylo_abort

echo "built src/wasm/genome_core.js ($(du -h src/wasm/genome_core.js | cut -f1))"
