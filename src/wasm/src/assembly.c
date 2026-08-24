/*
 * Whole-genome-overlap assembler kernel (WASM port of src/lib/assembly).
 *
 * Pipeline, identical to the JS reference implementation:
 *   minimizer sketch -> shared-minimizer index -> candidate pair tallying ->
 *   diagonal-chain overlap acceptance -> greedy end matching -> path walk ->
 *   greedy-layout consensus.
 *
 * Enumeration order, tie-breaking and float arithmetic mirror the JS code so
 * both engines emit byte-identical contigs for the same reads and params.
 */
#include "common.h"

typedef struct {
  int32_t k, w, maxOccurrence, diagTolerance, minMatches, minOverlapLen;
} AsmParams;

/* ------------------------------------------------------------- minimizers */

typedef struct { uint32_t hash, pos; int8_t strand; } Min;

static Min *compute_minimizers(const uint8_t *seq, size_t len, const AsmParams *p, size_t *outN) {
  *outN = 0;
  int32_t k = p->k;
  size_t nk;
  if (len < (size_t)k || k > 64) return NULL;
  nk = len - (size_t)k + 1;

  uint32_t *hashes = xmalloc(nk * sizeof(uint32_t));
  int8_t *strands = xmalloc(nk);
  uint8_t *valid = xcalloc(nk, 1);
  uint8_t km[64], rc[64];

  for (size_t i = 0; i < nk; i++) {
    memcpy(km, seq + i, (size_t)k);
    int ok = 1;
    for (int32_t c = 0; c < k; c++) {
      uint8_t ch = km[c] | 0x20;
      if (!(ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')) { ok = 0; break; }
    }
    if (!ok) continue;
    valid[i] = 1;
    for (int32_t c = 0; c < k; c++) rc[c] = (uint8_t)comp_base((char)km[k - 1 - c]);
    uint32_t hf = fnv1a(km, (size_t)k), hr = fnv1a(rc, (size_t)k);
    if (hf <= hr) { hashes[i] = hf; strands[i] = 1; }
    else          { hashes[i] = hr; strands[i] = -1; }
  }

  Min *out = xmalloc(nk * sizeof(Min));
  size_t cnt = 0;
  int64_t lastMinPos = -1;
  for (size_t i = 0; i + (size_t)p->w <= nk; i++) {
    int64_t minIdx = -1;
    uint32_t minHash = UINT32_MAX;
    for (size_t j = i; j < i + (size_t)p->w; j++) {
      if (!valid[j]) continue;
      if (hashes[j] < minHash) { minHash = hashes[j]; minIdx = (int64_t)j; }
    }
    if (minIdx != -1 && minIdx != lastMinPos) {
      out[cnt].hash = minHash;
      out[cnt].pos = (uint32_t)minIdx;
      out[cnt].strand = strands[minIdx];
      cnt++;
      lastMinPos = minIdx;
    }
  }
  free(hashes); free(strands); free(valid);
  *outN = cnt;
  return out;
}

/* -------------------------------------------------- hash map (u32 -> bucket) */

typedef struct { uint32_t *rd, *pos; int8_t *st; size_t n, cap; } Bucket;

typedef struct {
  uint32_t *keys; int32_t *slot; uint8_t *used;
  size_t mask, count;
  Bucket *buckets;
  uint32_t *dirKeys; /* keys in first-seen order */
} U32Map;

static void map_init(U32Map *m, size_t expect) {
  size_t sz = 16;
  while (sz < expect * 2) sz *= 2;
  m->keys = xcalloc(sz, sizeof(uint32_t));
  m->slot = xmalloc(sz * sizeof(int32_t));
  m->used = xcalloc(sz, 1);
  m->mask = sz - 1;
  m->count = 0;
  m->buckets = xcalloc(sz / 2 + 1, sizeof(Bucket));
  m->dirKeys = xcalloc(sz / 2 + 1, sizeof(uint32_t));
}
static void map_grow(U32Map *m) {
  U32Map n;
  size_t oldCap = m->mask + 1;
  map_init(&n, oldCap); /* sz' >= 2*oldCap; entry arrays sized sz'/2+1 >= oldCap+1 */
  n.count = m->count;
  memcpy(n.dirKeys, m->dirKeys, m->count * sizeof(uint32_t));
  for (size_t i = 0; i < m->count; i++) {
    uint32_t key = m->dirKeys[i];
    size_t h = ((size_t)(key * 2654435769u)) & n.mask;
    while (n.used[h]) h = (h + 1) & n.mask;
    n.used[h] = 1; n.keys[h] = key; n.slot[h] = (int32_t)i;
  }
  for (size_t i = 0; i < m->count; i++) n.buckets[i] = m->buckets[i];
  free(m->keys); free(m->slot); free(m->used); free(m->buckets); free(m->dirKeys);
  *m = n;
}
static Bucket *map_get_or_add(U32Map *m, uint32_t key) {
  if ((m->count + 1) * 2 > m->mask + 1) map_grow(m);
  size_t h = ((size_t)(key * 2654435769u)) & m->mask;
  while (m->used[h]) {
    if (m->keys[h] == key) return &m->buckets[m->slot[h]];
    h = (h + 1) & m->mask;
  }
  m->used[h] = 1; m->keys[h] = key; m->slot[h] = (int32_t)m->count;
  m->dirKeys[m->count] = key;
  memset(&m->buckets[m->count], 0, sizeof(Bucket));
  return &m->buckets[m->count++];
}
static void map_free(U32Map *m) {
  for (size_t i = 0; i < m->count; i++) {
    free(m->buckets[i].rd); free(m->buckets[i].pos); free(m->buckets[i].st);
  }
  free(m->buckets); free(m->dirKeys); free(m->keys); free(m->slot); free(m->used);
  memset(m, 0, sizeof(*m));
}

/* ------------------------------------------------------------ pair matches */

typedef struct { int32_t posA, posB; int8_t strandA, strandB; } Match;

static int cmp_u32_asc(const void *a, const void *b) {
  uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
  return x < y ? -1 : (x > y ? 1 : 0);
}

typedef struct {
  uint64_t *hkeys; int32_t *slot; uint8_t *used;
  size_t mask, count;
  Match **lists; size_t *lens, *caps;
  uint64_t *pairIds; /* (i<<32)|j by discovery order */
} PairMap;

static void pmap_init(PairMap *m, size_t expect) {
  size_t sz = 16;
  while (sz < expect * 2) sz *= 2;
  m->hkeys = xcalloc(sz, sizeof(uint64_t));
  m->slot = xmalloc(sz * sizeof(int32_t));
  m->used = xcalloc(sz, 1);
  m->mask = sz - 1; m->count = 0;
  m->lists = xcalloc(sz / 2 + 1, sizeof(Match *));
  m->lens = xcalloc(sz / 2 + 1, sizeof(size_t));
  m->caps = xcalloc(sz / 2 + 1, sizeof(size_t));
  m->pairIds = xcalloc(sz / 2 + 1, sizeof(uint64_t));
}
static void pmap_grow(PairMap *m) {
  PairMap n;
  size_t oldCap = m->mask + 1;
  pmap_init(&n, oldCap); /* sz' >= 2*oldCap; per-pair arrays sized sz'/2+1 */
  n.count = m->count;
  for (size_t i = 0; i < m->count; i++) {
    uint64_t key = m->pairIds[i];
    size_t h = (size_t)((key * 0x9E3779B97F4A7C15ull) >> 33) & n.mask;
    while (n.used[h]) h = (h + 1) & n.mask;
    n.used[h] = 1; n.hkeys[h] = key; n.slot[h] = (int32_t)i;
    n.pairIds[i] = key;
  }
  for (size_t i = 0; i < m->count; i++) {
    n.lists[i] = m->lists[i]; n.lens[i] = m->lens[i]; n.caps[i] = m->caps[i];
  }
  free(m->hkeys); free(m->slot); free(m->used);
  free(m->lists); free(m->lens); free(m->caps); free(m->pairIds);
  *m = n;
}
static int32_t pmap_get_or_add(PairMap *m, uint64_t key) {
  if ((m->count + 1) * 2 > m->mask + 1) pmap_grow(m);
  size_t h = (size_t)((key * 0x9E3779B97F4A7C15ull) >> 33) & m->mask;
  while (m->used[h]) {
    if (m->hkeys[h] == key) return m->slot[h];
    h = (h + 1) & m->mask;
  }
  m->used[h] = 1; m->hkeys[h] = key; m->slot[h] = (int32_t)m->count;
  m->pairIds[m->count] = key;
  m->lists[m->count] = NULL; m->lens[m->count] = 0; m->caps[m->count] = 0;
  return (int32_t)m->count++;
}
static void pmap_push(PairMap *m, int32_t id, Match mt) {
  if (m->lens[id] == m->caps[id]) {
    m->caps[id] = m->caps[id] ? m->caps[id] * 2 : 4;
    m->lists[id] = xrealloc(m->lists[id], m->caps[id] * sizeof(Match));
  }
  m->lists[id][m->lens[id]++] = mt;
}
static void pmap_free(PairMap *m) {
  for (size_t i = 0; i < m->count; i++) free(m->lists[i]);
  free(m->lists); free(m->lens); free(m->caps); free(m->pairIds);
  free(m->hkeys); free(m->slot); free(m->used);
  memset(m, 0, sizeof(*m));
}

/* --------------------------------------------------------------- overlaps */

enum { OT_SAME = 0, OT_REV = 1 };
enum { TT_A_IN_B, TT_B_IN_A, TT_A_TO_B, TT_B_TO_A };

typedef struct {
  uint32_t a, b;
  int8_t orientation; /* OT_SAME | OT_REV */
  int8_t type;
  uint32_t matchCount, overlapLen;
  double score;
  int32_t aMin, aMax, bMin, bMax;
} Overlap;

typedef struct { double posA, posBt; uint32_t ord; } ChainPt;

static int cmp_chainpt(const void *pa, const void *pb) {
  const ChainPt *a = (const ChainPt *)pa, *b = (const ChainPt *)pb;
  if (a->posA < b->posA) return -1;
  if (a->posA > b->posA) return 1;
  if (a->ord < b->ord) return -1;
  if (a->ord > b->ord) return 1;
  return 0;
}

/* chainBestChain port; fills best run bounds, returns its length (0 = reject). */
static size_t chain_best(ChainPt *arr, size_t n, const AsmParams *p,
                         size_t *bestStart, size_t *bestLenOut) {
  if (n < (size_t)p->minMatches) return 0;
  qsort(arr, n, sizeof(ChainPt), cmp_chainpt); /* ord tiebreak => stable */
  size_t bestLen = 0, bestS = 0, runStart = 0;
  for (size_t i = 1; i <= n; i++) {
    int brk = 0;
    if (i == n) brk = 1;
    else {
      double prevDiag = arr[i - 1].posA - arr[i - 1].posBt;
      double curDiag = arr[i].posA - arr[i].posBt;
      if (fabs(curDiag - prevDiag) > (double)p->diagTolerance) brk = 1;
    }
    if (brk) {
      size_t runLen = i - runStart;
      if (runLen > bestLen) { bestLen = runLen; bestS = runStart; }
      runStart = i;
    }
  }
  if (bestLen < (size_t)p->minMatches) return 0;
  *bestStart = bestS; *bestLenOut = bestLen;
  return bestLen;
}

static int estimate_overlap(const Match *ms, size_t nm, int32_t lenA, int32_t lenB,
                            const AsmParams *p, Overlap *ov) {
  if ((int64_t)nm < p->minMatches) return 0;

  ChainPt *same = xmalloc(nm * sizeof(ChainPt));
  ChainPt *rev = xmalloc(nm * sizeof(ChainPt));
  size_t ns = 0, nr = 0;
  for (size_t t = 0; t < nm; t++) {
    int sameOri = ms[t].strandA == ms[t].strandB;
    double posBt = sameOri ? (double)ms[t].posB : (double)(lenB - ms[t].posB - p->k);
    if (sameOri) { same[ns].posA = (double)ms[t].posA; same[ns].ord = (uint32_t)t; same[ns++].posBt = posBt; }
    else         { rev[nr].posA = (double)ms[t].posA; rev[nr].ord = (uint32_t)t; rev[nr++].posBt = posBt; }
  }

  int haveBest = 0;
  Overlap best;
  ChainPt *lists[2] = { same, rev };
  size_t nlens[2] = { ns, nr };

  for (int o = 0; o < 2; o++) {
    size_t st = 0, ln = 0;
    if (!chain_best(lists[o], nlens[o], p, &st, &ln)) continue;
    const ChainPt *ch = lists[o] + st;
    int32_t aMin = INT32_MAX, aMax = INT32_MIN, bMin = INT32_MAX, bMax = INT32_MIN;
    for (size_t q = 0; q < ln; q++) {
      int32_t av = (int32_t)ch[q].posA, bv = (int32_t)ch[q].posBt;
      if (av < aMin) aMin = av;
      if (av > aMax) aMax = av;
      if (bv < bMin) bMin = bv;
      if (bv > bMax) bMax = bv;
    }
    aMax += p->k; bMax += p->k;
    int32_t spanA = aMax - aMin, spanB = bMax - bMin;
    int32_t rawSpan = spanA > spanB ? spanA : spanB;
    int32_t overlapLen = rawSpan + 2 * p->w;
    if (overlapLen < p->minOverlapLen) continue;
    double expected = (double)rawSpan / (double)p->w;
    if (expected < 1) expected = 1;
    double score = (double)ln / expected;
    if (score > 1) score = 1;
    if (!haveBest || ln > (size_t)best.matchCount) {
      haveBest = 1;
      best.orientation = (int8_t)o;
      best.matchCount = (uint32_t)ln;
      best.overlapLen = (uint32_t)overlapLen;
      best.score = score;
      best.aMin = aMin; best.aMax = aMax; best.bMin = bMin; best.bMax = bMax;
    }
  }
  free(same); free(rev);
  if (!haveBest) return 0;

  int32_t slack = (2 * p->w + 15 > 60) ? 2 * p->w + 15 : 60;
  int aStartsAt0 = best.aMin <= slack, aEndsAtLen = best.aMax >= lenA - slack;
  int bStartsAt0 = best.bMin <= slack, bEndsAtLen = best.bMax >= lenB - slack;

  int type;
  if (aStartsAt0 && aEndsAtLen) type = TT_A_IN_B;
  else if (bStartsAt0 && bEndsAtLen) type = TT_B_IN_A;
  else if (aEndsAtLen && bStartsAt0) type = TT_A_TO_B;
  else if (bEndsAtLen && aStartsAt0) type = TT_B_TO_A;
  else return 0;

  best.type = (int8_t)type;
  *ov = best;
  return 1;
}

/* ------------------------------------------------------- graph & layouts */

typedef struct {
  int32_t endA, endB; /* node ids: readIdx*2 + (0=L, 1=R) */
  int flip;
  uint32_t ovIdx;
} Dovetail;

typedef struct {
  uint32_t readIdx;
  int32_t strand;
  int32_t entryEnd;
  uint32_t ovIdxPlus1; /* overlapWithPrev; 0 = none */
} Step;

typedef struct {
  Step *steps; uint32_t nSteps;
  int circular;
  uint32_t closingOvIdxPlus1;
} Path;

typedef struct {
  uint8_t *contained;
  Dovetail *dovetails; size_t nDovetail;
  int32_t *confTo;   /* per node: matched node or -1 */
  int32_t *confFlip;
  uint32_t *confOv;  /* ovIdx+1 of confirmed edge */
  Path *paths; size_t nPaths;
} Layout;

static void layout_build(Layout *L, uint32_t n, const Overlap *ovs, size_t nOv) {
  L->contained = xcalloc(n, 1);
  for (size_t t = 0; t < nOv; t++) {
    if (ovs[t].type == TT_A_IN_B) L->contained[ovs[t].a] = 1;
    else if (ovs[t].type == TT_B_IN_A) L->contained[ovs[t].b] = 1;
  }

  Dovetail *dt = xmalloc((nOv ? nOv : 1) * sizeof(Dovetail));
  size_t nd = 0;
  for (size_t t = 0; t < nOv; t++) {
    const Overlap *o = &ovs[t];
    if (!(o->type == TT_A_TO_B || o->type == TT_B_TO_A)) continue;
    if (L->contained[o->a] || L->contained[o->b]) continue;
    int flip = o->orientation == OT_REV;
    int endA, endB;
    if (o->type == TT_A_TO_B) { endA = 1; endB = flip ? 1 : 0; }
    else                      { endA = 0; endB = flip ? 0 : 1; }
    dt[nd].endA = (int32_t)((size_t)o->a * 2 + (size_t)endA);
    dt[nd].endB = (int32_t)((size_t)o->b * 2 + (size_t)endB);
    dt[nd].flip = flip;
    dt[nd].ovIdx = (uint32_t)t;
    nd++;
  }
  L->dovetails = dt; L->nDovetail = nd;

  /* Greedy maximum-weight matching: all edges sorted by strength desc
   * (stable), take each whose both ends are still free. */
  uint32_t *order = xmalloc((nd ? nd : 1) * sizeof(uint32_t));
  uint32_t *tmpo = xmalloc((nd ? nd : 1) * sizeof(uint32_t));
  double *keys = xmalloc((nd ? nd : 1) * sizeof(double));
  for (size_t t = 0; t < nd; t++) { order[t] = (uint32_t)t; keys[t] = -(double)ovs[dt[t].ovIdx].matchCount; }
  stable_sort_idx(keys, order, (uint32_t)nd, tmpo);
  free(keys); free(tmpo);

  size_t nn = (size_t)n * 2;
  L->confTo = xmalloc(nn * sizeof(int32_t));
  L->confFlip = xcalloc(nn, sizeof(int32_t));
  L->confOv = xcalloc(nn, sizeof(uint32_t));
  uint8_t *usedEnds = xcalloc(nn, 1);
  for (size_t i = 0; i < nn; i++) L->confTo[i] = -1;
  for (size_t t = 0; t < nd; t++) {
    const Dovetail *e = &dt[order[t]];
    if (usedEnds[e->endA] || usedEnds[e->endB]) continue;
    usedEnds[e->endA] = usedEnds[e->endB] = 1;
    L->confTo[e->endA] = e->endB; L->confFlip[e->endA] = e->flip; L->confOv[e->endA] = e->ovIdx + 1;
    L->confTo[e->endB] = e->endA; L->confFlip[e->endB] = e->flip; L->confOv[e->endB] = e->ovIdx + 1;
  }
  free(order); free(usedEnds);

  /* ---- extractPaths ---- */
  uint8_t *visitedEnds = xcalloc(nn, 1);
  Path *paths = xcalloc(n ? n : 1, sizeof(Path));
  size_t np = 0;

  for (int pass = 0; pass < 2; pass++) {
    for (uint32_t r = 0; r < n; r++) {
      if (L->contained[r]) continue;
      if (visitedEnds[r * 2] || visitedEnds[r * 2 + 1]) continue;
      int deg = (L->confTo[r * 2] >= 0) + (L->confTo[r * 2 + 1] >= 0);
      if (pass == 0 && deg == 2) continue;

      int32_t startEnd = L->confTo[r * 2 + 1] >= 0 ? 1 : 0;
      /* up to n+1 steps (first step + one per walk iteration) */
      Step *steps = xcalloc((size_t)n + 2, sizeof(Step));
      uint32_t ns = 0;
      steps[ns].readIdx = r; steps[ns].strand = 1;
      steps[ns].entryEnd = startEnd == 1 ? 0 : 1;
      steps[ns].ovIdxPlus1 = 0;
      ns++;
      visitedEnds[r * 2] = visitedEnds[r * 2 + 1] = 1;

      int32_t curIdx = (int32_t)r, curExitEnd = startEnd, curStrand = 1;
      int circular = 0;
      uint32_t closingOv = 0;
      int32_t guard = 0;
      while (guard++ < (int32_t)n + 5) {
        size_t exitNode = (size_t)curIdx * 2 + (size_t)curExitEnd;
        if (L->confTo[exitNode] < 0) break;
        int32_t nextIdx = L->confTo[exitNode] / 2;
        int32_t nextEntryEndNative = L->confTo[exitNode] % 2;
        int32_t flip = L->confFlip[exitNode];
        uint32_t ovP1 = L->confOv[exitNode];
        if (nextIdx == (int32_t)r && ns > 1) {
          circular = 1;
          closingOv = ovP1;
          break;
        }
        int32_t nextStrand = flip ? -curStrand : curStrand;
        steps[ns].readIdx = (uint32_t)nextIdx;
        steps[ns].strand = nextStrand;
        steps[ns].entryEnd = nextStrand == 1 ? nextEntryEndNative : (nextEntryEndNative == 0 ? 1 : 0);
        steps[ns].ovIdxPlus1 = ovP1;
        ns++;
        visitedEnds[nextIdx * 2] = visitedEnds[nextIdx * 2 + 1] = 1;
        curIdx = nextIdx;
        curExitEnd = nextEntryEndNative == 0 ? 1 : 0;
        curStrand = nextStrand;
      }
      paths[np].steps = xrealloc(steps, (ns ? ns : 1) * sizeof(Step));
      paths[np].nSteps = ns;
      paths[np].circular = circular || pass == 1;
      paths[np].closingOvIdxPlus1 = circular ? closingOv : 0;
      np++;
    }
  }
  free(visitedEnds);
  L->paths = xrealloc(paths, (np ? np : 1) * sizeof(Path));
  L->nPaths = np;
}

static void layout_free(Layout *L) {
  free(L->contained);
  free(L->dovetails);
  free(L->confTo); free(L->confFlip); free(L->confOv);
  for (size_t i = 0; i < L->nPaths; i++) free(L->paths[i].steps);
  free(L->paths);
  memset(L, 0, sizeof(*L));
}

/* -------------------------------------------------------------- consensus */

typedef struct {
  uint32_t readIdx;
  int32_t strand;
  uint32_t contigStart, contigEnd, trimmedFromStart;
} MemberRec;

typedef struct {
  uint8_t *seq; size_t len;
  MemberRec *members; uint32_t nMembers;
  int circular;
} ContigBuild;

/* boundary range of `readIdx` in the coordinate system of a step using it */
static void overlap_range_in_step(const Overlap *ov, uint32_t readIdx, int32_t len, int32_t strand,
                                  int32_t outStart[1], int32_t outEnd[1]) {
  int32_t native[2];
  if (readIdx == ov->a) {
    native[0] = ov->aMin; native[1] = ov->aMax;
  } else if (ov->orientation == OT_REV) {
    native[0] = len - ov->bMax; native[1] = len - ov->bMin;
  } else {
    native[0] = ov->bMin; native[1] = ov->bMax;
  }
  if (strand == 1) { *outStart = native[0]; *outEnd = native[1]; }
  else { *outStart = len - native[1]; *outEnd = len - native[0]; }
}

static void append_read_full(Buf *sb, const uint8_t *seqs, const uint32_t *offs,
                             uint32_t ri, int32_t strand) {
  const uint8_t *s = seqs + offs[ri];
  size_t l = offs[ri + 1] - offs[ri];
  buf_reserve(sb, l);
  if (strand == 1) {
    memcpy(sb->data + sb->len, s, l);
    sb->len += l;
  } else {
    for (size_t q = 0; q < l; q++) sb->data[sb->len++] = (uint8_t)comp_base((char)s[l - 1 - q]);
  }
}

static void append_read_suffix(Buf *sb, const uint8_t *seqs, const uint32_t *offs,
                               uint32_t ri, int32_t strand, int32_t trimStart) {
  const uint8_t *s = seqs + offs[ri];
  size_t l = offs[ri + 1] - offs[ri];
  buf_reserve(sb, l);
  if (strand == 1) {
    memcpy(sb->data + sb->len, s + trimStart, l - (size_t)trimStart);
    sb->len += l - (size_t)trimStart;
  } else {
    /* revcomp(read).slice(trimStart) */
    for (size_t q = (size_t)trimStart; q < l; q++)
      sb->data[sb->len++] = (uint8_t)comp_base((char)s[l - 1 - q]);
  }
}

static void build_contig(const uint8_t *seqs, const uint32_t *offs, const Overlap *ovs,
                         const Path *path, ContigBuild *out) {
  const Step *steps = path->steps;
  uint32_t nSteps = path->nSteps;

  Buf sb; buf_init(&sb);
  MemberRec *mem = xcalloc(nSteps ? nSteps : 1, sizeof(MemberRec));
  uint32_t nMem = 0;

  {
    uint32_t ri = steps[0].readIdx;
    size_t l = offs[ri + 1] - offs[ri];
    append_read_full(&sb, seqs, offs, ri, steps[0].strand);
    mem[nMem].readIdx = ri; mem[nMem].strand = steps[0].strand;
    mem[nMem].contigStart = 0; mem[nMem].contigEnd = (uint32_t)l;
    mem[nMem].trimmedFromStart = 0;
    nMem++;
  }

  for (uint32_t i = 1; i < nSteps; i++) {
    const Step *st = &steps[i];
    int32_t len = (int32_t)(offs[st->readIdx + 1] - offs[st->readIdx]);
    const Overlap *ov = &ovs[st->ovIdxPlus1 - 1];
    int32_t rs, re;
    overlap_range_in_step(ov, st->readIdx, len, st->strand, &rs, &re);
    int32_t trimStart = re < 0 ? 0 : re;
    if (trimStart > len) trimStart = len;

    size_t before = sb.len;
    append_read_suffix(&sb, seqs, offs, st->readIdx, st->strand, trimStart);
    mem[nMem].readIdx = st->readIdx; mem[nMem].strand = st->strand;
    mem[nMem].contigStart = (uint32_t)before;
    mem[nMem].contigEnd = (uint32_t)sb.len;
    mem[nMem].trimmedFromStart = (uint32_t)trimStart;
    nMem++;
  }

  int circular = path->circular;
  if (path->circular && path->closingOvIdxPlus1) {
    const Step *lastStep = &steps[nSteps - 1];
    int32_t lastLen = (int32_t)(offs[lastStep->readIdx + 1] - offs[lastStep->readIdx]);
    const Overlap *cov = &ovs[path->closingOvIdxPlus1 - 1];
    int32_t ovS, ovE;
    overlap_range_in_step(cov, lastStep->readIdx, lastLen, lastStep->strand, &ovS, &ovE);
    MemberRec *lm = &mem[nMem - 1];
    int32_t trimTo = (int32_t)lm->contigStart + (ovS - (int32_t)lm->trimmedFromStart);
    if (trimTo > (int32_t)lm->contigStart && trimTo < (int32_t)sb.len) {
      sb.len = (size_t)trimTo;
      lm->contigEnd = (uint32_t)trimTo;
    } else {
      circular = 0;
    }
  }

  out->seq = sb.data ? sb.data : xmalloc(1);
  out->len = sb.len;
  out->members = mem;
  out->nMembers = nMem;
  out->circular = circular;
}

/* ------------------------------------------------------------------ driver */

static Buf g_result;

int gx_assembly_run(const uint8_t *seqsFlat, const uint32_t *offsets, uint32_t nReads,
                    int32_t k, int32_t w, int32_t maxOccurrence, int32_t diagTolerance,
                    int32_t minMatches, int32_t minOverlapLen) {
  AsmParams P;
  P.k = k; P.w = w; P.maxOccurrence = maxOccurrence;
  P.diagTolerance = diagTolerance; P.minMatches = minMatches; P.minOverlapLen = minOverlapLen;
  buf_init(&g_result);

  /* ---- stage 1: minimizer sketch + shared index ---- */
  gx_progress(0, 5);
  U32Map index;
  {
    size_t estBuckets = 1024;
    Min **perRead = xcalloc(nReads ? nReads : 1, sizeof(Min *));
    size_t *perReadN = xcalloc(nReads ? nReads : 1, sizeof(size_t));
    for (uint32_t r = 0; r < nReads; r++) {
      perRead[r] = compute_minimizers(seqsFlat + offsets[r], offsets[r + 1] - offsets[r], &P, &perReadN[r]);
      estBuckets += perReadN[r] / 8 + 16;
    }
    map_init(&index, estBuckets);
    for (uint32_t r = 0; r < nReads; r++) {
      const Min *ms = perRead[r];
      for (size_t t = 0; t < perReadN[r]; t++) {
        Bucket *bk = map_get_or_add(&index, ms[t].hash);
        if (bk->n == bk->cap) {
          bk->cap = bk->cap ? bk->cap * 2 : 4;
          bk->rd = xrealloc(bk->rd, bk->cap * sizeof(uint32_t));
          bk->pos = xrealloc(bk->pos, bk->cap * sizeof(uint32_t));
          bk->st = xrealloc(bk->st, bk->cap);
        }
        bk->rd[bk->n] = r; bk->pos[bk->n] = ms[t].pos; bk->st[bk->n] = ms[t].strand;
        bk->n++;
      }
      free(perRead[r]);
    }
    free(perRead); free(perReadN);
  }

  /* repeat filtering: drop buckets seen in >maxOccurrence or <2 distinct reads */
  {
    uint32_t *tmpRd = NULL;
    size_t tmpCap = 0;
    uint8_t *keep = xcalloc(index.count ? index.count : 1, 1);
    for (size_t bi = 0; bi < index.count; bi++) {
      Bucket *bk = &index.buckets[bi];
      if ((size_t)bk->n > tmpCap) {
        tmpCap = bk->n;
        free(tmpRd);
        tmpRd = xmalloc(tmpCap * sizeof(uint32_t));
      }
      memcpy(tmpRd, bk->rd, (size_t)bk->n * sizeof(uint32_t));
      /* count distinct read ids */
      qsort(tmpRd, bk->n, sizeof(uint32_t), cmp_u32_asc);
      size_t distinct = 0;
      for (size_t q = 0; q < bk->n; q++)
        if (q == 0 || tmpRd[q] != tmpRd[q - 1]) distinct++;
      keep[bi] = (distinct <= (size_t)P.maxOccurrence && distinct >= 2) ? 1 : 0;
    }
    free(tmpRd);
    /* compact, preserving discovery order; then rebuild the hash table */
    size_t w2 = 0;
    for (size_t bi = 0; bi < index.count; bi++) {
      if (!keep[bi]) {
        free(index.buckets[bi].rd); free(index.buckets[bi].pos); free(index.buckets[bi].st);
        continue;
      }
      index.dirKeys[w2] = index.dirKeys[bi];
      index.buckets[w2] = index.buckets[bi];
      w2++;
    }
    index.count = w2;
    memset(index.used, 0, index.mask + 1);
    for (size_t i = 0; i < w2; i++) {
      uint32_t key = index.dirKeys[i];
      size_t h = ((size_t)(key * 2654435769u)) & index.mask;
      while (index.used[h]) h = (h + 1) & index.mask;
      index.used[h] = 1; index.keys[h] = key; index.slot[h] = (int32_t)i;
    }
    free(keep);
  }

  gx_progress(1, 35);

  /* ---- stage 2: candidate pairs ---- */
  gx_progress(2, 45);
  PairMap pm;
  pmap_init(&pm, index.count * 4 + 64);
  for (size_t bi = 0; bi < index.count; bi++) {
    const Bucket *bk = &index.buckets[bi];
    for (size_t a = 0; a < bk->n; a++) {
      for (size_t b = a + 1; b < bk->n; b++) {
        uint32_t xr = bk->rd[a], yr = bk->rd[b];
        if (xr == yr) continue;
        uint32_t firstR = xr < yr ? xr : yr, secondR = xr < yr ? yr : xr;
        int fa = xr < yr; /* entry a belongs to `firstR` when xr<yr */
        uint32_t pa = fa ? bk->pos[a] : bk->pos[b];
        int8_t sa = (int8_t)(fa ? bk->st[a] : bk->st[b]);
        uint32_t pb = fa ? bk->pos[b] : bk->pos[a];
        int8_t sb2 = (int8_t)(fa ? bk->st[b] : bk->st[a]);
        uint64_t key = ((uint64_t)firstR << 32) | secondR;
        Match mt;
        mt.posA = (int32_t)pa; mt.posB = (int32_t)pb;
        mt.strandA = sa; mt.strandB = sb2;
        pmap_push(&pm, pmap_get_or_add(&pm, key), mt);
      }
    }
  }

  gx_progress(3, 45);

  /* ---- stage 3: overlap estimation per candidate pair ---- */
  Overlap *ovs = xmalloc((pm.count ? pm.count : 1) * sizeof(Overlap));
  size_t nOv = 0;
  for (size_t pi = 0; pi < pm.count; pi++) {
    uint32_t i = (uint32_t)(pm.pairIds[pi] >> 32);
    uint32_t j = (uint32_t)(pm.pairIds[pi] & 0xffffffffu);
    int32_t lenI = (int32_t)(offsets[i + 1] - offsets[i]);
    int32_t lenJ = (int32_t)(offsets[j + 1] - offsets[j]);
    Overlap ov;
    if (estimate_overlap(pm.lists[pi], pm.lens[pi], lenI, lenJ, &P, &ov)) {
      ov.a = i; ov.b = j;
      ovs[nOv++] = ov;
    }
    /* same progress cadence as the JS worker (every 200 pairs) */
    if (((pi + 1) % 200) == 0) {
      double frac = (double)(pi + 1) / (double)(pm.count ? pm.count : 1);
      gx_progress(3, 45 + (int)(frac * 30.0 + 0.5));
    }
  }
  pmap_free(&pm);

  Layout L;
  layout_build(&L, nReads, ovs, nOv);
  gx_progress(4, 80);
  size_t nn2 = (size_t)nReads * 2;

  /* ---- stage 4: consensus ---- */
  gx_progress(5, 90);
  ContigBuild *contigs = xcalloc(L.nPaths ? L.nPaths : 1, sizeof(ContigBuild));
  for (size_t pi = 0; pi < L.nPaths; pi++) {
    build_contig(seqsFlat, offsets, ovs, &L.paths[pi], &contigs[pi]);
  }

  /* sort contigs by length desc (stable) */
  if (L.nPaths) {
    uint32_t *ord = xmalloc(L.nPaths * sizeof(uint32_t));
    uint32_t *tmp = xmalloc(L.nPaths * sizeof(uint32_t));
    double *keys = xmalloc(L.nPaths * sizeof(double));
    for (size_t t = 0; t < L.nPaths; t++) { ord[t] = (uint32_t)t; keys[t] = -(double)contigs[t].len; }
    stable_sort_idx(keys, ord, (uint32_t)L.nPaths, tmp);
    free(keys); free(tmp);
    ContigBuild *sorted = xcalloc(L.nPaths, sizeof(ContigBuild));
    for (size_t t = 0; t < L.nPaths; t++) sorted[t] = contigs[ord[t]];
    free(ord); free(contigs);
    contigs = sorted;
  }

  /* ---- serialize result blob ---- */
  {
    Buf *o = &g_result;
    uint32_t nContigs = (uint32_t)L.nPaths;
    uint32_t version = 1;
    buf_append(o, "GXAR", 4);
    buf_append(o, &version, 4);
    buf_append(o, &nContigs, 4);
    for (uint32_t t = 0; t < nContigs; t++) {
      ContigBuild *c = &contigs[t];
      uint32_t slen = (uint32_t)c->len, circ = (uint32_t)c->circular, resv = 0;
      buf_append(o, &slen, 4);
      buf_append(o, &circ, 4);
      buf_append(o, &c->nMembers, 4);
      buf_append(o, &resv, 4);
      buf_append(o, c->seq, c->len);
      while (o->len % 4) buf_push(o, 0);
      for (uint32_t mm = 0; mm < c->nMembers; mm++) {
        buf_append(o, &c->members[mm].readIdx, 4);
        buf_append(o, &c->members[mm].strand, 4);
        buf_append(o, &c->members[mm].contigStart, 4);
        buf_append(o, &c->members[mm].contigEnd, 4);
        buf_append(o, &c->members[mm].trimmedFromStart, 4);
      }
    }
    /* stats tail */
    uint64_t totalLength = 0, longest = 0, circCount = 0, used = 0;
    for (uint32_t t = 0; t < nContigs; t++) {
      totalLength += contigs[t].len;
      if (contigs[t].len > longest) longest = contigs[t].len;
      circCount += contigs[t].circular ? 1 : 0;
      used += contigs[t].nMembers;
    }
    double n50 = 0;
    if (nContigs) {
      uint32_t *lens = xmalloc(nContigs * sizeof(uint32_t));
      for (uint32_t t = 0; t < nContigs; t++) lens[t] = (uint32_t)contigs[t].len;
      /* sort descending (insertion; contig counts are modest) */
      for (uint32_t a = 1; a < nContigs; a++) {
        uint32_t v = lens[a];
        int64_t b = (int64_t)a - 1;
        while (b >= 0 && lens[b] < v) { lens[b + 1] = lens[b]; b--; }
        lens[b + 1] = v;
      }
      uint64_t tot = 0;
      for (uint32_t t = 0; t < nContigs; t++) tot += lens[t];
      double target = (double)tot * 0.5;
      uint64_t running = 0;
      n50 = (double)lens[nContigs - 1];
      for (uint32_t t = 0; t < nContigs; t++) {
        running += lens[t];
        if ((double)running >= target) { n50 = (double)lens[t]; break; }
      }
      free(lens);
    }
    double meanScore = 0;
    for (size_t t = 0; t < nOv; t++) meanScore += ovs[t].score;
    if (nOv) meanScore /= (double)nOv;
    uint32_t containedCount = 0;
    for (uint32_t r = 0; r < nReads; r++) containedCount += L.contained[r];

    double f64s[2] = { n50, meanScore };
    buf_append(o, f64s, sizeof(f64s));
    uint32_t u32s[8] = {
      nReads, containedCount, (uint32_t)used,
      nReads - containedCount - (uint32_t)used,
      (uint32_t)nOv, nContigs, (uint32_t)totalLength, (uint32_t)longest
    };
    buf_append(o, u32s, sizeof(u32s));
    uint32_t circV = (uint32_t)circCount;
    buf_append(o, &circV, 4);

  }

  /* cleanup */
  for (size_t t = 0; t < L.nPaths; t++) {
    free(contigs[t].seq);
    free(contigs[t].members);
  }
  free(contigs);
  layout_free(&L);
  free(ovs);
  return 0;
}

const uint8_t *gx_assembly_result_ptr(void) { return g_result.data; }
uint32_t gx_assembly_result_len(void) { return (uint32_t)g_result.len; }
void gx_assembly_free(void) { buf_free(&g_result); }
