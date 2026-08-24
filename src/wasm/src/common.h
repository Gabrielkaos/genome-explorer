/*
 * Shared plumbing for the Genome Explorer WASM kernels (assembly, MSA, phylo).
 * All kernels are deterministic ports of the reference JS implementations in
 * src/lib/{assembly,msa,phylo} - same algorithms, same tie-breaking, so both
 * engines produce identical output on identical input.
 */
#ifndef GX_COMMON_H
#define GX_COMMON_H

#define _POSIX_C_SOURCE 199309L

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
/* Implemented via EM_JS in gx_glue.c (must live in a single TU). */
void gx_progress(int code, int pct);
double gx_now(void);
#else
/* Non-Emscripten parse context (clangd/native tests): silent no-op / POSIX clock. */
#include <time.h>
static void gx_progress(int code, int pct) { (void)code; (void)pct; }
static double gx_now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec * 1000.0 + ts.tv_nsec / 1e6;
}
#endif

static void *xmalloc(size_t n) {
  void *p = malloc(n ? n : 1);
  if (!p) abort();
  return p;
}

static void *xcalloc(size_t n, size_t sz) {
  void *p = calloc(n ? n : 1, sz ? sz : 1);
  if (!p) abort();
  return p;
}

static void *xrealloc(void *p, size_t n) {
  void *q = realloc(p, n ? n : 1);
  if (!q) abort();
  return q;
}

/* Growable byte buffer. */
typedef struct { uint8_t *data; size_t len, cap; } Buf;

static void buf_init(Buf *b) { b->data = NULL; b->len = b->cap = 0; }
static void buf_reserve(Buf *b, size_t extra) {
  if (b->len + extra > b->cap) {
    size_t cap = b->cap ? b->cap * 2 : 256;
    while (cap < b->len + extra) cap *= 2;
    b->data = (uint8_t *)xrealloc(b->data, cap);
    b->cap = cap;
  }
}
static void buf_push(Buf *b, uint8_t c) { buf_reserve(b, 1); b->data[b->len++] = c; }
static void buf_append(Buf *b, const void *src, size_t n) {
  buf_reserve(b, n); memcpy(b->data + b->len, src, n); b->len += n;
}
static void buf_free(Buf *b) { free(b->data); buf_init(b); }

/* FNV-1a 32-bit over bytes - must match hashKmer() in sequence.js exactly. */
static inline uint32_t fnv1a(const uint8_t *s, size_t n) {
  uint32_t h = 2166136261u;
  for (size_t i = 0; i < n; i++) { h ^= s[i]; h *= 16777619u; }
  return h;
}

static inline char comp_base(char c) {
  switch (c) {
    case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
    case 'a': return 't'; case 'c': return 'g'; case 'g': return 'c'; case 't': return 'a';
    default: return 'N';
  }
}

/* Stable sort of n indices by (key asc, index asc). keys are doubles. */
static void stable_sort_idx(const double *keys, uint32_t *idx, uint32_t n, uint32_t *tmp) {
  if (n < 2) return;
  for (uint32_t width = 1; width < n; width *= 2) {
    /* snapshot current permutation, then merge one level into idx */
    memcpy(tmp, idx, n * sizeof(uint32_t));
    for (uint32_t lo = 0; lo < n; lo += 2 * width) {
      uint32_t mid = lo + width < n ? lo + width : n;
      uint32_t hi = lo + 2 * width < n ? lo + 2 * width : n;
      uint32_t i = lo, j = mid, k = lo;
      while (i < mid && j < hi) {
        /* strict < keeps earlier index first on ties => stable */
        if (keys[tmp[j]] < keys[tmp[i]]) idx[k++] = tmp[j++];
        else idx[k++] = tmp[i++];
      }
      while (i < mid) idx[k++] = tmp[i++];
      while (j < hi) idx[k++] = tmp[j++];
    }
  }
}

#endif
