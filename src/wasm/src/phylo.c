/*
 * Phylogenetic inference kernel (WASM port of src/lib/phylo/distances|nj|
 * upgma|splits). Distance models (p / JC69 / K2P), Neighbor-Joining and
 * UPGMA mirror the JS reference implementations bit-for-bit, including
 * saturation caps, limb-clamp accounting and Q-criterion tie-breaking.
 *
 * The nonparametric bootstrap runs through a chunked API (init -> chunks ->
 * finish) so the calling worker can yield to the event loop between
 * replicates and keep the progress bar alive. Split voting tallies the
 * canonical (smaller-side) taxon bitmask of every non-trivial internal edge,
 * matching splits.js canonKeyOf() semantics exactly - including the
 * lexicographic comma-list tie-break when a split divides taxa 50/50.
 */
#include "common.h"
#include <stdio.h>
#include <float.h>

typedef struct {
  int32_t model;    /* 0=p, 1=jc, 2=k2p */
  int32_t gapMode;  /* 0=pairwise, 1=complete */
  int32_t method;   /* 0=nj, 1=upgma */
} PhyloOpts;

#define GX_SATURATION_CAP 1.5

/* --------------------------------------------------------------- RNG */

/* mulberry32 - bit-exact port of the PRNG in phylo.worker.js */
static uint32_t g_rngState;

static void rng_seed(uint32_t a) { g_rngState = a; }

static double rng_next(void) {
  uint32_t a = g_rngState;
  a += 0x6D2B79F5u;
  uint32_t t = (uint32_t)((a ^ (a >> 15)) * (1u | a));
  t = ((uint32_t)((t ^ (t >> 7)) * (61u | t)) + t) ^ t;
  g_rngState = a;
  return (double)(t ^ (t >> 14)) / 4294967296.0;
}

/* ----------------------------------------------------- distance matrices */

typedef struct {
  double *d, *tiM, *tvM, *cmpM;
  uint32_t columnsUsed;
  uint32_t saturatedPairs, noDataPairs;
  double meanDistance;
  double tiTvRatio; /* NaN means null */
  int32_t closestPair[2], furthestPair[2];
} DistResult;

static int is_transition(int a, int b) {
  return (a == 0 && b == 2) || (a == 2 && b == 0) || (a == 1 && b == 3) || (a == 3 && b == 1);
}

/*
 * rowsFlat: n x L encoded codes. colIdx (may be NULL): column subset.
 * Mirrors computeDistanceMatrix() including its order of operations: the
 * column universe becomes colIdx when provided, THEN complete-deletion
 * filtering is applied over that universe.
 */
static void dist_matrix(const uint8_t *rowsFlat, uint32_t n, uint32_t L,
                        const PhyloOpts *o, const uint32_t *colIdxIn, uint32_t nColIn,
                        DistResult *R) {
  R->d = xcalloc((size_t)n * n, sizeof(double));
  R->tiM = xcalloc((size_t)n * n, sizeof(double));
  R->tvM = xcalloc((size_t)n * n, sizeof(double));
  R->cmpM = xcalloc((size_t)n * n, sizeof(double));

  uint32_t *cols;
  if (colIdxIn) {
    cols = xmalloc(nColIn * sizeof(uint32_t));
    memcpy(cols, colIdxIn, nColIn * sizeof(uint32_t));
  } else {
    cols = xmalloc(L * sizeof(uint32_t));
    for (uint32_t c = 0; c < L; c++) cols[c] = c;
  }
  uint32_t nCols = colIdxIn ? nColIn : L;

  if (o->gapMode == 1) {
    uint32_t *keep = xmalloc((size_t)nCols * sizeof(uint32_t));
    uint32_t nk = 0;
    for (uint32_t kk = 0; kk < nCols; kk++) {
      uint32_t c = cols[kk];
      int ok = 1;
      for (uint32_t r = 0; r < n; r++)
        if (rowsFlat[(size_t)r * L + c] > 3) { ok = 0; break; }
      if (ok) keep[nk++] = c;
    }
    free(cols);
    cols = keep;
    nCols = nk;
  }
  R->columnsUsed = nCols;

  R->saturatedPairs = 0;
  R->noDataPairs = 0;
  R->meanDistance = 0;
  R->tiTvRatio = NAN;
  R->closestPair[0] = R->closestPair[1] = -1;
  R->furthestPair[0] = R->furthestPair[1] = -1;

  for (uint32_t i = 0; i < n; i++) {
    const uint8_t *ri = rowsFlat + (size_t)i * L;
    R->d[(size_t)i * n + i] = 0;
    for (uint32_t j = i + 1; j < n; j++) {
      const uint8_t *rj = rowsFlat + (size_t)j * L;
      uint64_t diffs = 0, ti = 0, tv = 0, comp = 0;
      for (uint32_t kk = 0; kk < nCols; kk++) {
        uint8_t a = ri[cols[kk]], b = rj[cols[kk]];
        if (a > 3 || b > 3) continue;
        comp++;
        if (a != b) {
          diffs++;
          if (is_transition(a, b)) ti++;
          else tv++;
        }
      }
      double raw = comp ? (double)diffs / (double)comp : 0;
      double dist = raw;
      if (o->model == 1) {
        if (raw > 0 && raw < 0.75) dist = -0.75 * log(1 - (4 * raw) / 3);
        else if (raw >= 0.75) dist = INFINITY;
      } else if (o->model == 2) {
        double Pp = comp ? (double)ti / (double)comp : 0;
        double Qq = comp ? (double)tv / (double)comp : 0;
        double s = 1 - 2 * Pp - Qq;
        double v = 1 - 2 * Qq;
        if ((s <= 0 || v <= 0) && (Pp > 0 || Qq > 0)) dist = INFINITY;
        else dist = -0.5 * log(fmax(s, DBL_MIN)) - 0.25 * log(fmax(v, DBL_MIN));
      }
      if (!isfinite(dist)) {
        R->saturatedPairs++;
        dist = GX_SATURATION_CAP;
      } else if (dist != raw) {
        /* quantize like the JS kernel so both engines agree despite libm
         * last-bit differences in log(); raw p-distances are already exact */
        dist = floor(dist * 1e9 + 0.5) / 1e9;
      }
      if (comp == 0) R->noDataPairs++;
      R->d[(size_t)i * n + j] = R->d[(size_t)j * n + i] = dist;
      R->tiM[(size_t)i * n + j] = R->tiM[(size_t)j * n + i] = (double)ti;
      R->tvM[(size_t)i * n + j] = R->tvM[(size_t)j * n + i] = (double)tv;
      R->cmpM[(size_t)i * n + j] = R->cmpM[(size_t)j * n + i] = (double)comp;
    }
  }
  free(cols);

  double tiSum = 0, tvSum = 0, sum = 0;
  uint64_t count = 0;
  double minD = INFINITY, maxD = 0;
  R->closestPair[0] = R->closestPair[1] = -1;
  R->furthestPair[0] = R->furthestPair[1] = -1;
  for (uint32_t i = 0; i < n; i++) {
    for (uint32_t j = i + 1; j < n; j++) {
      double v = R->d[(size_t)i * n + j];
      tiSum += R->tiM[(size_t)i * n + j];
      tvSum += R->tvM[(size_t)i * n + j];
      sum += v; count++;
      if (v < minD) { minD = v; R->closestPair[0] = (int32_t)i; R->closestPair[1] = (int32_t)j; }
      if (v > maxD) { maxD = v; R->furthestPair[0] = (int32_t)i; R->furthestPair[1] = (int32_t)j; }
    }
  }
  R->meanDistance = count ? sum / (double)count : 0;
  R->tiTvRatio = tvSum > 0 ? tiSum / tvSum : NAN;
}

static void dist_result_free(DistResult *R) {
  free(R->d); free(R->tiM); free(R->tvM); free(R->cmpM);
  memset(R, 0, sizeof(*R));
}

/* ------------------------------------------------------------------ trees */

/* Small-degree tree: up to 3 children (the NJ root trifurcates). */
typedef struct {
  uint8_t leaf;
  int32_t nameId;   /* leaf: taxon index; internal: unused */
  int32_t child[3];
  int32_t nChild;
} PNode;

typedef struct {
  PNode *nodes;
  int32_t count, cap;
  int32_t root;
  uint64_t clampedLimbs;
} PTree;

/* per-node child branch lengths (child index -> length) */
typedef struct { double bl[3]; } BLTable;

static void ptree_init(PTree *T, uint32_t nTaxa) {
  T->cap = (int32_t)(2 * nTaxa);
  T->nodes = xcalloc((size_t)T->cap, sizeof(PNode));
  T->count = (int32_t)nTaxa; /* leaves preallocated as slots 0..n-1 */
  T->root = -1;
  T->clampedLimbs = 0;
  for (int32_t i = 0; i < T->count; i++) {
    T->nodes[i].leaf = 1;
    T->nodes[i].nameId = i;
    T->nodes[i].nChild = 0;
    T->nodes[i].child[0] = T->nodes[i].child[1] = T->nodes[i].child[2] = -1;
  }
}

static int32_t ptree_new_internal(PTree *T, int32_t a, int32_t b, int32_t c) {
  if (T->count >= T->cap) {
    T->cap *= 2;
    T->nodes = xrealloc(T->nodes, (size_t)T->cap * sizeof(PNode));
  }
  PNode *nd = &T->nodes[T->count];
  nd->leaf = 0;
  nd->nameId = -1;
  nd->nChild = c >= 0 ? 3 : 2;
  nd->child[0] = a; nd->child[1] = b; nd->child[2] = c;
  return T->count++;
}

/* --- Neighbor-Joining (Saitou & Nei), topology + branch lengths together --- */

static void nj_build(const double *distIn, uint32_t n, PTree *T, BLTable **blOut) {
  ptree_init(T, n);
  BLTable *bl = xcalloc((size_t)T->cap, sizeof(BLTable));
  double *D = xmalloc((size_t)n * n * sizeof(double));
  memcpy(D, distIn, (size_t)n * n * sizeof(double));
  double *r = xmalloc(n * sizeof(double));

  /* live positions map matrix indices to node slots; starts as identity */
  int32_t *live = xmalloc(n * sizeof(int32_t));
  for (uint32_t i = 0; i < n; i++) live[i] = (int32_t)i;

  int32_t m = (int32_t)n;
  while (m > 3) {
    for (int32_t i = 0; i < m; i++) {
      double s = 0;
      for (int32_t j = 0; j < m; j++) if (j != i) s += D[(size_t)i * m + j];
      r[i] = s;
    }
    int32_t bi = -1, bj = -1;
    double bq = INFINITY;
    for (int32_t i = 0; i < m; i++) {
      for (int32_t j = i + 1; j < m; j++) {
        double q = (double)(m - 2) * D[(size_t)i * m + j] - r[i] - r[j];
        if (q < bq) { bq = q; bi = i; bj = j; }
      }
    }
    double dij = D[(size_t)bi * m + bj];
    double liIdeal = 0.5 * dij + (r[bi] - r[bj]) / (2.0 * (double)(m - 2));
    double li = fmax(0, liIdeal);
    double lj = fmax(0, dij - li);
    if (li != liIdeal) T->clampedLimbs++;
    if (lj != dij - li) T->clampedLimbs++;

    int32_t uSlot = ptree_new_internal(T, live[bi], live[bj], -1);
    bl[uSlot].bl[0] = li;
    bl[uSlot].bl[1] = lj;

    /* reduced matrix: replace position bi with u, drop bj (mirrors nj.js) */
    double *next = xcalloc((size_t)(m - 1) * (m - 1), sizeof(double));
    int32_t *mapTo = xmalloc(m * sizeof(int32_t));
    int32_t t2 = 0;
    for (int32_t kk = 0; kk < m; kk++) {
      if (kk == bj) { mapTo[kk] = -1; continue; }
      mapTo[kk] = (kk == bi) ? t2++ : t2;
      if (kk != bi) t2++;
    }
    for (int32_t a = 0; a < m; a++) {
      if (a == bj) continue;
      for (int32_t b = 0; b < m; b++) {
        if (b == bj || b == a) continue;
        next[(size_t)mapTo[a] * (m - 1) + mapTo[b]] = D[(size_t)a * m + b];
      }
    }
    for (int32_t a = 0; a < m; a++) {
      if (a == bi || a == bj) continue;
      double dau = (D[(size_t)a * m + bi] + D[(size_t)a * m + bj] - dij) / 2;
      next[(size_t)mapTo[bi] * (m - 1) + mapTo[a]] = dau;
      next[(size_t)mapTo[a] * (m - 1) + mapTo[bi]] = dau;
    }

    int32_t w = 0;
    for (int32_t kk = 0; kk < m; kk++) {
      if (kk == bj) continue;
      live[w++] = (kk == bi) ? uSlot : live[kk];
    }

    free(D);
    D = next;
    free(mapTo);
    m--;
  }
  free(r);

  /* final three nodes form the trifurcating (unrooted) root; limbs are
   * clamped without incrementing clampedLimbs, like nj.js */
  double dxy = D[0 * 3 + 1], dxz = D[0 * 3 + 2], dyz = D[1 * 3 + 2];
  double lx = fmax(0, (dxy + dxz - dyz) / 2);
  double ly = fmax(0, (dxy + dyz - dxz) / 2);
  double lz = fmax(0, (dxz + dyz - dxy) / 2);
  T->root = ptree_new_internal(T, live[0], live[1], live[2]);
  bl[T->root].bl[0] = lx;
  bl[T->root].bl[1] = ly;
  bl[T->root].bl[2] = lz;

  free(D);
  free(live);
  *blOut = bl;
}

/* --- UPGMA (rooted, molecular clock) --- */

typedef struct {
  uint32_t size;
  int32_t *members;
  uint8_t alive;
  int32_t nodeSlot;
  double height;
} UCluster;

static void upgma_build(const double *dist, uint32_t n, PTree *T, BLTable **blOut) {
  ptree_init(T, n);
  BLTable *bl = xcalloc((size_t)T->cap, sizeof(BLTable));

  UCluster *cl = xcalloc((size_t)(2 * n - 1 ? 2 * n - 1 : 1), sizeof(UCluster));
    /* leaves + every merged cluster's full member list:
   n + sum(k, k=2..n) = n + (n*n+n-2)/2 entries */
  int32_t *memPool = xmalloc(((size_t)n * (n + 3) / 2 + 16) * sizeof(int32_t));
  int32_t poolTop = (int32_t)n;
  for (uint32_t i = 0; i < n; i++) {
    cl[i].size = 1; cl[i].members = memPool + i; cl[i].members[0] = (int32_t)i;
    cl[i].alive = 1; cl[i].nodeSlot = (int32_t)i; cl[i].height = 0;
  }
  int32_t nextSlot = (int32_t)n;
  int32_t aliveCount = (int32_t)n;

  while (aliveCount > 1) {
    int32_t bi = -1, bj = -1;
    double bd = INFINITY;
    for (int32_t i = 0; i < nextSlot; i++) {
      if (!cl[i].alive) continue;
      for (int32_t j = i + 1; j < nextSlot; j++) {
        if (!cl[j].alive) continue;
        double sum = 0;
        for (uint32_t a = 0; a < cl[i].size; a++) {
          const double *row = dist + (size_t)cl[i].members[a] * n;
          for (uint32_t b = 0; b < cl[j].size; b++) sum += row[cl[j].members[b]];
        }
        double dAvg = sum / ((double)cl[i].size * (double)cl[j].size);
        if (dAvg < bd) { bd = dAvg; bi = i; bj = j; }
      }
    }
    double height = bd / 2;
    double branchA = fmax(0, height - cl[bi].height);
    double branchB = fmax(0, height - cl[bj].height);
    if (height - cl[bi].height < 0 || height - cl[bj].height < 0) T->clampedLimbs++;

    cl[bi].alive = 0; cl[bj].alive = 0;
    int32_t slot = nextSlot++;
    T->nodes[slot].leaf = 0;
    T->nodes[slot].nChild = 2;
    T->nodes[slot].child[0] = cl[bi].nodeSlot;
    T->nodes[slot].child[1] = cl[bj].nodeSlot;
    T->nodes[slot].child[2] = -1;
    bl[slot].bl[0] = branchA;
    bl[slot].bl[1] = branchB;

    uint32_t nsz = cl[bi].size + cl[bj].size;
    cl[slot].size = nsz;
    cl[slot].members = memPool + poolTop;
    memcpy(cl[slot].members, cl[bi].members, cl[bi].size * sizeof(int32_t));
    memcpy(cl[slot].members + cl[bi].size, cl[bj].members, cl[bj].size * sizeof(int32_t));
    poolTop += (int32_t)nsz;
    cl[slot].alive = 1;
    cl[slot].nodeSlot = slot;
    cl[slot].height = height;
    aliveCount--;
  }

  for (int32_t i = 0; i < nextSlot; i++)
    if (cl[i].alive) { T->root = cl[i].nodeSlot; break; }
  free(memPool); free(cl);
  *blOut = bl;
}

/* ------------------------------------------------------- split collection */

typedef struct {
  uint8_t **masks;   /* representative canonical mask per distinct split */
  uint32_t *counts;
  size_t maskBytes;
  size_t count;
  size_t cap;
  int32_t *slot;     /* open addressing on FNV(mask bytes) */
  uint8_t *used;
  size_t smask;
} SplitTally;

static uint32_t fnv_mask(const uint8_t *m, size_t nBytes) {
  uint32_t h = 2166136261u;
  for (size_t i = 0; i < nBytes; i++) { h ^= m[i]; h *= 16777619u; }
  return h;
}

static void tally_init(SplitTally *S, size_t expect, size_t maskBytes) {
  size_t sz = 16;
  while (sz < expect * 2) sz *= 2;
  S->smask = sz - 1;
  S->slot = xmalloc(sz * sizeof(int32_t));
  S->used = xcalloc(sz, 1);
  S->cap = expect + 16;
  S->masks = xmalloc(S->cap * sizeof(uint8_t *));
  S->counts = xmalloc(S->cap * sizeof(uint32_t));
  S->maskBytes = maskBytes;
  S->count = 0;
}
static void tally_rehash(SplitTally *S) {
  free(S->slot); free(S->used);
  size_t sz = (S->smask + 1) * 2;
  S->smask = sz - 1;
  S->slot = xmalloc(sz * sizeof(int32_t));
  S->used = xcalloc(sz, 1);
  for (size_t i = 0; i < S->count; i++) {
    size_t h = fnv_mask(S->masks[i], S->maskBytes) & S->smask;
    while (S->used[h]) h = (h + 1) & S->smask;
    S->used[h] = 1; S->slot[h] = (int32_t)i;
  }
}
static void tally_add(SplitTally *S, const uint8_t *canonMask) {
  if ((S->count + 1) * 2 > S->smask + 1) tally_rehash(S);
  size_t h = fnv_mask(canonMask, S->maskBytes) & S->smask;
  while (S->used[h]) {
    int32_t id = S->slot[h];
    if (memcmp(S->masks[id], canonMask, S->maskBytes) == 0) { S->counts[id]++; return; }
    h = (h + 1) & S->smask;
  }
  if ((size_t)S->count >= S->cap) {
    S->cap *= 2;
    S->masks = xrealloc(S->masks, S->cap * sizeof(uint8_t *));
    S->counts = xrealloc(S->counts, S->cap * sizeof(uint32_t));
  }
  uint8_t *copy = xmalloc(S->maskBytes);
  memcpy(copy, canonMask, S->maskBytes);
  S->masks[S->count] = copy;
  S->counts[S->count] = 1;
  S->used[h] = 1;
  S->slot[h] = (int32_t)S->count;
  S->count++;
}
static void tally_free(SplitTally *S) {
  for (size_t i = 0; i < S->count; i++) free(S->masks[i]);
  free(S->masks); free(S->counts); free(S->slot); free(S->used);
  memset(S, 0, sizeof(*S));
}

/* ascending comma-list of set bits, e.g. "0,3,7" */
static void mask_to_csv(const uint8_t *mask, uint32_t total, char *out) {
  size_t op = 0;
  int first = 1;
  for (uint32_t i = 0; i < total; i++) {
    if ((mask[i >> 3] >> (i & 7)) & 1) {
      if (first) { op += (size_t)sprintf(out + op, "%u", i); first = 0; }
      else op += (size_t)sprintf(out + op, ",%u", i);
    }
  }
  out[op] = 0;
}

/*
 * Canonicalize one edge given its own-side membership mask, replicating
 * canonKeyOf(): smaller side wins; at exactly n/2 the side whose ascending
 * comma-list compares lexicographically SMALLER wins.
 */
static void tally_edge(SplitTally *S, const uint8_t *ownMask, uint32_t ownSize,
                       uint32_t total) {
  if (ownSize <= 1 || ownSize >= total) return;

  if ((uint64_t)ownSize * 2 < (uint64_t)total) {
    tally_add(S, ownMask);
    return;
  }

  uint8_t *comp = xcalloc(S->maskBytes, 1);
  for (uint32_t i = 0; i < total; i++)
    if (!((ownMask[i >> 3] >> (i & 7)) & 1)) comp[i >> 3] |= (uint8_t)(1u << (i & 7));

  if ((uint64_t)ownSize * 2 == (uint64_t)total) {
    size_t bufCap = (size_t)total * 13 + 16;
    char *ownS = xmalloc(bufCap), *compS = xmalloc(bufCap);
    mask_to_csv(ownMask, total, ownS);
    mask_to_csv(comp, total, compS);
    if (strcmp(compS, ownS) < 0) {
      tally_add(S, comp);
    } else {
      tally_add(S, ownMask);
    }
    free(ownS); free(compS);
  } else {
    tally_add(S, comp);
  }
  free(comp);
}

/* Vote every non-trivial internal-edge split of the tree. Iterative DFS. */
static void collect_splits(const PTree *T, SplitTally *S, uint32_t total) {
  int32_t cap = T->cap;

  /* leaf ids in DFS pre-order (children visited in stored order) */
  int32_t *stk = xmalloc((size_t)cap * sizeof(int32_t));
  int32_t *leafOf = xmalloc((size_t)cap * sizeof(int32_t));
  for (int32_t i = 0; i < cap; i++) leafOf[i] = -1;
  int32_t sp = 0, cursor = 0;
  stk[sp++] = T->root;
  while (sp > 0) {
    int32_t nid = stk[--sp];
    const PNode *nd = &T->nodes[nid];
    if (nd->leaf) { leafOf[nid] = cursor++; continue; }
    for (int32_t c = nd->nChild - 1; c >= 0; c--) stk[sp++] = nd->child[c];
  }

  typedef struct { int32_t node; uint8_t expanded; } Frame;
  Frame *frames = xmalloc((size_t)cap * sizeof(Frame));
  uint8_t **masksByNode = xcalloc((size_t)cap, sizeof(uint8_t *));
  uint32_t *sizeByNode = xcalloc((size_t)cap, sizeof(uint32_t));

  sp = 0;
  frames[sp].node = T->root; frames[sp].expanded = 0; sp++;
  while (sp > 0) {
    Frame *fr = &frames[sp - 1];
    const PNode *nd = &T->nodes[fr->node];
    if (nd->leaf) {
      uint8_t *mk = xcalloc(S->maskBytes, 1);
      mk[leafOf[fr->node] >> 3] |= (uint8_t)(1u << (leafOf[fr->node] & 7));
      masksByNode[fr->node] = mk;
      sizeByNode[fr->node] = 1;
      sp--;
      continue;
    }
    if (!fr->expanded) {
      fr->expanded = 1;
      for (int32_t c = 0; c < nd->nChild; c++) {
        frames[sp].node = nd->child[c];
        frames[sp].expanded = 0;
        sp++;
      }
      continue;
    }
    uint8_t *mk = xcalloc(S->maskBytes, 1);
    uint32_t szSum = 0;
    for (int32_t c = 0; c < nd->nChild; c++) {
      int32_t cid = nd->child[c];
      const uint8_t *cm = masksByNode[cid];
      tally_edge(S, cm, sizeByNode[cid], total);
      for (size_t b = 0; b < S->maskBytes; b++) mk[b] |= cm[b];
      szSum += sizeByNode[cid];
    }
    masksByNode[fr->node] = mk;
    sizeByNode[fr->node] = szSum;
    for (int32_t c = 0; c < nd->nChild; c++) {
      free(masksByNode[T->nodes[fr->node].child[c]]);
      masksByNode[T->nodes[fr->node].child[c]] = NULL;
    }
    sp--;
  }

  free(frames); free(stk); free(leafOf); free(masksByNode); free(sizeByNode);
}

/* ------------------------------------------------------------- driver */

static Buf g_result;

/* Persistent state across init / bootstrap chunks / finish */
static struct {
  int active;
  uint32_t n, L;
  PhyloOpts opts;
  DistResult mainDist;
  PTree mainTree;
  BLTable *mainBL;
  uint8_t *rowsFlat;
  uint32_t seed;
  int32_t repsTarget, repsDone;
  uint64_t clampedReps;
  SplitTally tally;
} ST;

int32_t gx_phylo_init(const uint8_t *rowsFlat, uint32_t n, uint32_t L,
                      int32_t model, int32_t gapMode, int32_t method,
                      int32_t reps, uint32_t seed) {
  memset(&ST, 0, sizeof(ST));
  buf_init(&g_result);
  ST.rowsFlat = xmalloc((size_t)n * L);
  memcpy(ST.rowsFlat, rowsFlat, (size_t)n * L);
  ST.n = n; ST.L = L;
  ST.opts.model = model; ST.opts.gapMode = gapMode; ST.opts.method = method;
  ST.seed = seed;
  ST.repsTarget = reps;
  ST.repsDone = 0;
  ST.clampedReps = 0;

  dist_matrix(ST.rowsFlat, n, L, &ST.opts, NULL, 0, &ST.mainDist);

  if (method == 1) upgma_build(ST.mainDist.d, n, &ST.mainTree, &ST.mainBL);
  else nj_build(ST.mainDist.d, n, &ST.mainTree, &ST.mainBL);

  size_t maskBytes = ((size_t)n + 7) / 8;
  tally_init(&ST.tally, (size_t)n * 4 + 16, maskBytes);
  rng_seed(seed); /* JS seeds once before the replicate loop */
  ST.active = 1;
  return 0;
}

/* Runs up to `howMany` bootstrap replicates; returns replicates completed.
 * RNG state persists across chunks, matching the JS generator semantics. */
int32_t gx_phylo_bootstrap_chunk(int32_t howMany) {
  if (!ST.active || howMany <= 0) return 0;
  uint32_t n = ST.n, L = ST.L;
  uint32_t cu = ST.mainDist.columnsUsed;
  int32_t doneHere = 0;

  uint32_t *sample = xmalloc((cu ? cu : 1) * sizeof(uint32_t));
  while (ST.repsDone < ST.repsTarget && doneHere < howMany) {
    for (uint32_t kk = 0; kk < cu; kk++)
      sample[kk] = (uint32_t)(rng_next() * (double)cu);

    DistResult rd;
    dist_matrix(ST.rowsFlat, n, L, &ST.opts, sample, cu, &rd);
    PTree rt;
    BLTable *rbl = NULL;
    if (ST.opts.method == 1) upgma_build(rd.d, n, &rt, &rbl);
    else nj_build(rd.d, n, &rt, &rbl);
    ST.clampedReps += rt.clampedLimbs;

    collect_splits(&rt, &ST.tally, n);

    free(rt.nodes); free(rbl);
    dist_result_free(&rd);

    ST.repsDone++;
    doneHere++;
  }
  free(sample);
  return doneHere;
}

int32_t gx_phylo_reps_done(void) { return ST.repsDone; }

static void serialize_and_cleanup(void) {
  Buf *o = &g_result;
  uint32_t n = ST.n, L = ST.L;
  uint32_t version = 1;
  buf_append(o, "GXPR", 4);
  buf_append(o, &version, 4);
  buf_append(o, &n, 4);
  buf_append(o, &L, 4);
  uint32_t cu = ST.mainDist.columnsUsed;
  buf_append(o, &cu, 4);

  /* main tree, serialized in DFS preorder with children in stored order */
  PTree *T = &ST.mainTree;
  int32_t nn = T->count;
  buf_append(o, &nn, 4);
  int32_t *order = xmalloc((size_t)nn * sizeof(int32_t));
  int32_t *pos = xmalloc((size_t)nn * sizeof(int32_t));
  int32_t *parent = xmalloc((size_t)nn * sizeof(int32_t));
  int32_t *nameIdx = xmalloc((size_t)nn * sizeof(int32_t));
  uint8_t *childCount = xmalloc((size_t)nn);
  int32_t *stk = xmalloc((size_t)nn * sizeof(int32_t));
  int32_t sp = 0, oi = 0;
  stk[sp++] = T->root;
  while (sp > 0) {
    int32_t nid = stk[--sp];
    order[oi] = nid; pos[nid] = oi; oi++;
    const PNode *nd = &T->nodes[nid];
    for (int32_t c = nd->nChild - 1; c >= 0; c--) stk[sp++] = nd->child[c];
  }
  for (int32_t idx = 0; idx < nn; idx++) parent[idx] = -1;
  for (int32_t idx = 0; idx < nn; idx++) {
    const PNode *nd = &T->nodes[order[idx]];
    for (int32_t c = 0; c < nd->nChild; c++) parent[pos[nd->child[c]]] = idx;
    nameIdx[idx] = nd->leaf ? nd->nameId : -1;
    childCount[idx] = (uint8_t)nd->nChild;
  }
  buf_append(o, parent, (size_t)nn * sizeof(int32_t));
  buf_append(o, nameIdx, (size_t)nn * sizeof(int32_t));
  buf_append(o, childCount, nn);
  while (o->len % 8) buf_push(o, 0);
  double *blens = xmalloc((size_t)nn * 3 * sizeof(double));
  for (int32_t idx = 0; idx < nn; idx++) {
    int32_t nid = order[idx];
    const PNode *nd = &T->nodes[nid];
    for (int32_t c = 0; c < 3; c++)
      blens[idx * 3 + c] = c < nd->nChild ? ST.mainBL[nid].bl[c] : 0;
  }
  buf_append(o, blens, (size_t)nn * 3 * sizeof(double));
  free(blens);
  free(stk); free(order); free(pos); free(parent); free(nameIdx); free(childCount);

  /* summary numbers */
  double fs[3] = { ST.mainDist.meanDistance, ST.mainDist.tiTvRatio, (double)ST.mainTree.clampedLimbs };
  buf_append(o, fs, sizeof(fs));
  int32_t pairs[4] = {
    ST.mainDist.closestPair[0], ST.mainDist.closestPair[1],
    ST.mainDist.furthestPair[0], ST.mainDist.furthestPair[1]
  };
  buf_append(o, pairs, sizeof(pairs));
  uint32_t us[3] = { ST.mainDist.saturatedPairs, ST.mainDist.noDataPairs,
                     (uint32_t)ST.clampedReps };
  buf_append(o, us, sizeof(us));
  while (o->len % 8) buf_push(o, 0); /* f64-align the matrices */

  /* matrices */
  buf_append(o, ST.mainDist.d, (size_t)n * n * sizeof(double));
  buf_append(o, ST.mainDist.tiM, (size_t)n * n * sizeof(double));
  buf_append(o, ST.mainDist.tvM, (size_t)n * n * sizeof(double));
  buf_append(o, ST.mainDist.cmpM, (size_t)n * n * sizeof(double));

  /* bootstrap split tally */
  uint32_t nsplits = (uint32_t)ST.tally.count;
  uint32_t maskBytes = (uint32_t)ST.tally.maskBytes;
  buf_append(o, &nsplits, 4);
  buf_append(o, &maskBytes, 4);
  for (size_t i = 0; i < ST.tally.count; i++) buf_append(o, &ST.tally.counts[i], 4);
  for (size_t i = 0; i < ST.tally.count; i++) buf_append(o, ST.tally.masks[i], ST.tally.maskBytes);

  /* cleanup */
  dist_result_free(&ST.mainDist);
  free(ST.mainTree.nodes);
  free(ST.mainBL);
  free(ST.rowsFlat);
  tally_free(&ST.tally);
  ST.active = 0;
}

int32_t gx_phylo_finish(void) {
  if (ST.active) serialize_and_cleanup();
  return 0;
}

const uint8_t *gx_phylo_result_ptr(void) { return g_result.data; }
uint32_t gx_phylo_result_len(void) { return (uint32_t)g_result.len; }

void gx_phylo_abort(void) {
  if (ST.active) {
    dist_result_free(&ST.mainDist);
    free(ST.mainTree.nodes);
    free(ST.mainBL);
    free(ST.rowsFlat);
    tally_free(&ST.tally);
    ST.active = 0;
  }
  buf_free(&g_result);
}
