/*
 * Progressive multiple-sequence-alignment kernel (WASM port of src/lib/msa).
 *
 * Stages, identical to the JS reference implementation:
 *   k-mer distance estimation -> UPGMA guide tree -> depth-first progressive
 *   merge with profile-profile global DP (Gotoh affine gaps, ClustalW-style
 *   column scoring) -> alignment statistics.
 *
 * The DP traceback, tie-breaking order and Float32 rounding points mirror
 * profile.js exactly; sequences are compared bit-for-bit against the JS
 * engine by scripts/parity.mjs.
 */
#include "common.h"
#include <time.h>

#define GX_MSA_MAX_CELLS 900000000.0 /* trace bytes cap; JS reference used 160e6 */

/* ------------------------------------------------------------ k-mer stage */

static int32_t msa_k_for(const uint32_t *lens, uint32_t n, int32_t reqK) {
  uint32_t minLen = lens[0];
  for (uint32_t i = 1; i < n; i++) if (lens[i] < minLen) minLen = lens[i];
  int32_t k = (int32_t)minLen - 1;
  if (reqK < k) k = reqK;
  if (k < 3) k = 3;
  return k;
}

static void kmer_counts(const uint8_t *seq, uint32_t len, int32_t k, int32_t *vec, int32_t *total) {
  int32_t mask = (1 << (2 * k)) - 1;
  memset(vec, 0, ((size_t)1 << (2 * k)) * sizeof(int32_t));
  int32_t code = 0, run = 0, tot = 0;
  for (uint32_t i = 0; i < len; i++) {
    uint8_t b = seq[i];
    if (b > 3) { run = 0; continue; }
    code = ((code << 2) | b) & mask;
    if (++run >= k) { vec[code]++; tot++; }
  }
  *total = tot;
}

static double kmer_distance(const int32_t *av, const int32_t *bv, int32_t size) {
  int64_t inter = 0, sumA = 0, sumB = 0;
  for (int32_t i = 0; i < size; i++) {
    int64_t x = av[i], y = bv[i];
    inter += x < y ? x : y;
    sumA += x; sumB += y;
  }
  double dUnion = (double)(sumA + sumB - inter);
  if (dUnion == 0) return 1;
  double p = 1 - (double)inter / dUnion;
  if (p <= 0) return 0;
  if (p >= 0.75) return 3;
  double d = -0.75 * log(1 - (4 * p) / 3);
  /* quantize like the JS kernel so guide trees agree despite libm ULPs */
  return floor(d * 1e9 + 0.5) / 1e9;
}

/* --------------------------------------------------------- UPGMA guide tree */

typedef struct {
  uint8_t leaf;
  int32_t idx;        /* leaf sequence index */
  int32_t left, right;/* internal children (cluster slots) */
} GTNode;

typedef struct {
  uint32_t size;
  int32_t *members;
  uint8_t alive;
  int32_t node;
} GCluster;

/* Builds guide tree; returns root node id. Node ids into *nodesOut (count n*2-1). */
static int32_t msa_upgma(const double *dist, uint32_t n, GTNode **nodesOut) {
  GTNode *nodes = xcalloc(n * 2 - 1 ? n * 2 - 1 : 1, sizeof(GTNode));
  GCluster *cl = xcalloc(n * 2 - 1 ? n * 2 - 1 : 1, sizeof(GCluster));
    /* leaves + every merged cluster's full member list:
   n + sum(k, k=2..n) = n + (n*n+n-2)/2 entries */
  int32_t *memPool = xmalloc(((size_t)n * (n + 3) / 2 + 16) * sizeof(int32_t));

  for (uint32_t i = 0; i < n; i++) {
    nodes[i].leaf = 1; nodes[i].idx = (int32_t)i;
    cl[i].size = 1; cl[i].members = memPool + i; cl[i].members[0] = (int32_t)i; cl[i].alive = 1; cl[i].node = (int32_t)i;
  }
  int32_t nextNode = (int32_t)n;
  int32_t poolTop = (int32_t)n;
  int32_t aliveCount = (int32_t)n;

  while (aliveCount > 1) {
    int32_t bi = -1, bj = -1;
    double bd = INFINITY;
    for (int32_t i = 0; i < nextNode; i++) {
      if (!cl[i].alive) continue;
      for (int32_t j = i + 1; j < nextNode; j++) {
        if (!cl[j].alive) continue;
        /* average linkage over member pairs */
        double sum = 0;
        for (uint32_t a = 0; a < cl[i].size; a++) {
          const double *row = dist + (size_t)cl[i].members[a] * n;
          for (uint32_t b = 0; b < cl[j].size; b++) sum += row[cl[j].members[b]];
        }
        double dAvg = sum / ((double)cl[i].size * (double)cl[j].size);
        if (dAvg < bd) { bd = dAvg; bi = i; bj = j; }
      }
    }
    cl[bi].alive = 0; cl[bj].alive = 0;
    uint32_t nsz = cl[bi].size + cl[bj].size;
    int32_t slot = nextNode++;
    nodes[slot].leaf = 0; nodes[slot].left = cl[bi].node; nodes[slot].right = cl[bj].node;
    cl[slot].size = nsz;
    cl[slot].members = memPool + poolTop;
    memcpy(cl[slot].members, cl[bi].members, cl[bi].size * sizeof(int32_t));
    memcpy(cl[slot].members + cl[bi].size, cl[bj].members, cl[bj].size * sizeof(int32_t));
    poolTop += (int32_t)nsz;
    cl[slot].alive = 1;
    cl[slot].node = slot;
    aliveCount--;
  }

  int32_t root = -1;
  for (int32_t i = 0; i < nextNode; i++) if (cl[i].alive) { root = cl[i].node; break; }
  free(memPool); free(cl);
  *nodesOut = nodes;
  return root;
}

/* ------------------------------------------------------------- profiles */

typedef struct {
  uint32_t len;
  uint32_t *cols;      /* 5 * len residue counts */
  uint32_t nMem;
  int32_t *origIdx;
  uint8_t *chars;      /* nMem * len codes */
} Profile;

static void profile_free(Profile *p) {
  free(p->cols); free(p->origIdx); free(p->chars);
  memset(p, 0, sizeof(*p));
}

static Profile *profile_single(const uint8_t *seq, uint32_t len, int32_t origIdx) {
  Profile *p = xcalloc(1, sizeof(Profile));
  p->len = len;
  p->cols = xcalloc((size_t)5 * (len ? len : 1), sizeof(uint32_t));
  p->nMem = 1;
  p->origIdx = xmalloc(sizeof(int32_t));
  p->origIdx[0] = origIdx;
  p->chars = xmalloc(len ? len : 1);
  memcpy(p->chars, seq, len);
  for (uint32_t i = 0; i < len; i++) p->cols[(size_t)i * 5 + seq[i]]++;
  return p;
}

/*
 * Global profile-profile alignment (Gotoh, ClustalW column scores).
 * ops[k]: 0 = diag, 1 = gap in Q (consume P column), 2 = gap in P.
 */
static int merge_profiles(const Profile *P, const Profile *Q, double match, double mismatch,
                          double gapOpen, double gapExtend, double gapGapBonus,
                          Profile **outOps) {
  uint32_t m = P->len, n = Q->len;
  if (!m || !n) return -1;
  if ((double)m * (double)n > GX_MSA_MAX_CELLS) return -2;

  double coef = match + mismatch;

  float *fx = xmalloc((size_t)5 * m * sizeof(float));
  float *fy = xmalloc((size_t)5 * n * sizeof(float));
  for (uint32_t c = 0; c < m * 5; c++) fx[c] = (float)((double)P->cols[c] / (double)P->nMem);
  for (uint32_t c = 0; c < n * 5; c++) fy[c] = (float)((double)Q->cols[c] / (double)Q->nMem);

  float *tx = xmalloc((size_t)5 * m * sizeof(float));
  float *ty = xmalloc((size_t)5 * n * sizeof(float));
  const float *fs[2] = { fx, fy };
  float *ts[2] = { tx, ty };
  uint32_t lens[2] = { m, n };
  for (int which = 0; which < 2; which++) {
    const float *f = fs[which];
    float *t = ts[which];
    for (uint32_t c = 0; c < lens[which]; c++) {
      double resTot = 0;
      for (int a = 0; a < 4; a++) resTot += f[c * 5 + a];
      for (int a = 0; a < 4; a++)
        t[c * 5 + a] = (float)((double)f[c * 5 + a] * coef - resTot * mismatch);
    }
  }

  float *wx = xmalloc(m * sizeof(float)), *wy = xmalloc(n * sizeof(float));
  for (uint32_t c = 0; c < m; c++) wx[c] = (float)(1.0 - (double)fx[c * 5 + 4]);
  for (uint32_t c = 0; c < n; c++) wy[c] = (float)(1.0 - (double)fy[c * 5 + 4]);

  const double NEG = -INFINITY;
  uint8_t *trace = xcalloc((size_t)m * n, 1);
  double *prevM = xmalloc((n + 1) * sizeof(double));
  double *prevX = xmalloc((n + 1) * sizeof(double));
  double *prevY = xmalloc((n + 1) * sizeof(double));
  double *curM = xmalloc((n + 1) * sizeof(double));
  double *curX = xmalloc((n + 1) * sizeof(double));
  double *curY = xmalloc((n + 1) * sizeof(double));

  prevM[0] = 0;
  for (uint32_t j = 1; j <= n; j++) prevM[j] = NEG;
  for (uint32_t j = 0; j <= n; j++) { prevX[j] = NEG; prevY[j] = NEG; }

  for (uint32_t i = 1; i <= m; i++) {
    curM[0] = NEG;
    curY[0] = NEG;
    curX[0] = 0; /* free leading-gap chain down the left edge */

    uint32_t ti = (i - 1) * 5;
    for (uint32_t j = 1; j <= n; j++) {
      uint32_t tj = (j - 1) * 5;
      size_t cell = ((size_t)(i - 1)) * n + (j - 1);

      double s = 0;
      for (int a = 0; a < 4; a++) s += (double)fx[ti + a] * (double)ty[tj + a];
      s += gapGapBonus * (double)fx[ti + 4] * (double)fy[tj + 4];

      double best = prevM[j - 1];
      int arg = 0;
      if (prevX[j - 1] > best) { best = prevX[j - 1]; arg = 1; }
      if (prevY[j - 1] > best) { best = prevY[j - 1]; arg = 2; }
      curM[j] = (best == NEG) ? NEG : best + s;
      trace[cell] |= (uint8_t)arg;

      int oSrcX = prevY[j] >= prevM[j] ? 2 : 0;
      double oValX = prevY[j] >= prevM[j] ? prevY[j] : prevM[j];
      double openX = (oValX == NEG) ? NEG : oValX - gapOpen * (double)wy[j - 1];
      double extX = (prevX[j] == NEG) ? NEG : prevX[j] - gapExtend * (double)wy[j - 1];
      if (extX != NEG && extX >= openX) { curX[j] = extX; trace[cell] |= (3 << 2); }
      else if (openX != NEG) { curX[j] = openX; trace[cell] |= (oSrcX << 2); }
      else curX[j] = NEG;

      int oSrcY = curX[j - 1] >= curM[j - 1] ? 1 : 0;
      double oValY = curX[j - 1] >= curM[j - 1] ? curX[j - 1] : curM[j - 1];
      double openY = (oValY == NEG) ? NEG : oValY - gapOpen * (double)wx[i - 1];
      double extY = (curY[j - 1] == NEG) ? NEG : curY[j - 1] - gapExtend * (double)wx[i - 1];
      if (extY != NEG && extY >= openY) { curY[j] = extY; trace[cell] |= (3 << 4); }
      else if (openY != NEG) { curY[j] = openY; trace[cell] |= (oSrcY << 4); }
      else curY[j] = NEG;
    }

    double *sw;
    sw = prevM; prevM = curM; curM = sw;
    sw = prevX; prevX = curX; curX = sw;
    sw = prevY; prevY = curY; curY = sw;
  }

  /* end-state selection */
  int state = 0;
  {
    double bestEnd = prevM[n];
    if (prevX[n] > bestEnd) { bestEnd = prevX[n]; state = 1; }
    if (prevY[n] > bestEnd) state = 2;
  }

  /* traceback (reverse-order op collection, then reverse) */
  uint8_t *ops = xmalloc((size_t)m + n);
  size_t nOps = 0;
  {
    uint32_t i = m, j = n;
    while (i > 0 && j > 0) {
      size_t cell = ((size_t)(i - 1)) * n + (j - 1);
      if (state == 0) {
        ops[nOps++] = 0;
        state = trace[cell] & 3;
        i--; j--;
      } else if (state == 1) {
        ops[nOps++] = 1;
        int info = (trace[cell] >> 2) & 3;
        state = info == 3 ? 1 : info;
        i--;
      } else {
        ops[nOps++] = 2;
        int info = (trace[cell] >> 4) & 3;
        state = info == 3 ? 2 : info;
        j--;
      }
    }
    while (i > 0) { ops[nOps++] = 1; i--; }
    while (j > 0) { ops[nOps++] = 2; j--; }
    for (size_t a = 0; a < nOps / 2; a++) {
      uint8_t t2 = ops[a]; ops[a] = ops[nOps - 1 - a]; ops[nOps - 1 - a] = t2;
    }
  }

  /* buildMerged */
  Profile *R = xcalloc(1, sizeof(Profile));
  uint32_t Lo = (uint32_t)nOps;
  R->len = Lo;
  R->cols = xcalloc((size_t)5 * Lo, sizeof(uint32_t));
  R->nMem = P->nMem + Q->nMem;
  R->origIdx = xmalloc((size_t)R->nMem * sizeof(int32_t));
  R->chars = xmalloc((size_t)R->nMem * Lo);

  int32_t *mapP = xmalloc(m * sizeof(int32_t));
  int32_t *mapQ = xmalloc(n * sizeof(int32_t));
  for (uint32_t c = 0; c < m; c++) mapP[c] = -1;
  for (uint32_t c = 0; c < n; c++) mapQ[c] = -1;

  {
    uint32_t pI = 0, qI = 0;
    for (size_t kk = 0; kk < nOps; kk++) {
      if (ops[kk] == 0) {
        for (int a = 0; a < 5; a++) R->cols[kk * 5 + a] = P->cols[pI * 5 + a] + Q->cols[qI * 5 + a];
        mapP[pI] = (int32_t)kk; mapQ[qI] = (int32_t)kk;
        pI++; qI++;
      } else if (ops[kk] == 1) {
        memcpy(R->cols + kk * 5, P->cols + pI * 5, 5 * sizeof(uint32_t));
        mapP[pI] = (int32_t)kk;
        pI++;
      } else {
        memcpy(R->cols + kk * 5, Q->cols + qI * 5, 5 * sizeof(uint32_t));
        mapQ[qI] = (int32_t)kk;
        qI++;
      }
    }
  }

  /* remap member rows through the column maps (unmapped -> gap code 4) */
  {
    uint32_t wRow = 0;
    for (uint32_t mi = 0; mi < P->nMem; mi++, wRow++) {
      R->origIdx[wRow] = P->origIdx[mi];
      uint8_t *dst = R->chars + (size_t)wRow * Lo;
      memset(dst, 4, Lo);
      const uint8_t *src = P->chars + (size_t)mi * P->len;
      for (uint32_t c = 0; c < P->len; c++) {
        int32_t nc = mapP[c];
        if (nc >= 0) dst[nc] = src[c];
      }
    }
    for (uint32_t mi = 0; mi < Q->nMem; mi++, wRow++) {
      R->origIdx[wRow] = Q->origIdx[mi];
      uint8_t *dst = R->chars + (size_t)wRow * Lo;
      memset(dst, 4, Lo);
      const uint8_t *src = Q->chars + (size_t)mi * Q->len;
      for (uint32_t c = 0; c < Q->len; c++) {
        int32_t nc = mapQ[c];
        if (nc >= 0) dst[nc] = src[c];
      }
    }
  }

  free(mapP); free(mapQ); free(ops);
  free(trace); free(fx); free(fy); free(tx); free(ty); free(wx); free(wy);
  free(prevM); free(prevX); free(prevY); free(curM); free(curX); free(curY);
  *outOps = R;
  return 0;
}

/* ----------------------------------------------------------- statistics */

typedef struct {
  uint32_t conservedColumns, variableColumns, informativeColumns, singletonColumns, gapColumns;
  uint32_t variableSitesTotal;
  double gapFraction, meanPairIdentity, minPairIdentity, maxPairIdentity;
} AlnStats;

/* colClass bits: 1 conserved, 2 variable, 4 informative, 8 singleton, 16 gapCol */
static void aln_stats(const Profile *root, const uint8_t *rowsByOrig, AlnStats *st,
                      uint8_t *colClass, uint8_t *consensusCodes, double *identity,
                      double *perSeqGapPct) {
  uint32_t length = root->len;
  uint32_t nSeq = root->nMem;

  uint32_t conserved = 0, variable = 0, informative = 0, singleton = 0, gapCols = 0;
  uint64_t totalGapCells = 0;
  uint32_t varTotal = 0;

  for (uint32_t c = 0; c < length; c++) {
    const uint32_t *cnt = root->cols + (size_t)c * 5;
    uint32_t gapCount = cnt[4];
    totalGapCells += gapCount;

    uint32_t bestCount = 0;
    int32_t bestCode = 4;
    uint32_t distinct = 0, resTotal = 0;
    for (int a = 0; a < 4; a++) {
      if (!cnt[a]) continue;
      resTotal += cnt[a];
      distinct++;
      if (cnt[a] > bestCount) { bestCount = cnt[a]; bestCode = a; }
    }
    if (resTotal == 0) {
      consensusCodes[c] = 4;
      gapCols++;
      colClass[c] = 16;
      continue;
    }
    consensusCodes[c] = (uint8_t)((double)bestCount / (double)nSeq >= 0.5 ? bestCode : 4);

    if (distinct == 1 && resTotal == nSeq) {
      conserved++;
      colClass[c] = 1;
    } else {
      variable++;
      varTotal++;
      uint32_t infStates = 0;
      for (int a = 0; a < 4; a++) if (cnt[a] >= 2) infStates++;
      uint8_t cls = 2;
      if (infStates >= 2) { informative++; cls |= 4; }
      else { singleton++; cls |= 8; }
      colClass[c] = cls;
    }
  }

  for (uint32_t i = 0; i < nSeq; i++) {
    uint32_t g = 0;
    const uint8_t *row = rowsByOrig + (size_t)i * length;
    for (uint32_t c = 0; c < length; c++) if (row[c] == 4) g++;
    perSeqGapPct[i] = length ? (double)g / (double)length : 0;
  }

  for (uint32_t i = 0; i < nSeq; i++) identity[(size_t)i * nSeq + i] = 1;
  double sumId = 0;
  uint32_t pairs = 0;
  double minId = 1, maxId = 0;
  for (uint32_t i = 0; i < nSeq; i++) {
    const uint8_t *A = rowsByOrig + (size_t)i * length;
    for (uint32_t j = i + 1; j < nSeq; j++) {
      const uint8_t *B = rowsByOrig + (size_t)j * length;
      uint64_t matchN = 0, compared = 0;
      for (uint32_t c = 0; c < length; c++) {
        if (A[c] > 3 || B[c] > 3) continue;
        compared++;
        if (A[c] == B[c]) matchN++;
      }
      double id = compared ? (double)matchN / (double)compared : 0;
      identity[(size_t)i * nSeq + j] = identity[(size_t)j * nSeq + i] = id;
      sumId += id; pairs++;
      if (id < minId) minId = id;
      if (id > maxId) maxId = id;
    }
  }

  st->conservedColumns = conserved;
  st->variableColumns = variable;
  st->informativeColumns = informative;
  st->singletonColumns = singleton;
  st->gapColumns = gapCols;
  st->variableSitesTotal = varTotal;
  st->gapFraction = length ? (double)totalGapCells / ((double)length * (double)nSeq) : 0;
  st->meanPairIdentity = pairs ? sumId / (double)pairs : 1;
  st->minPairIdentity = pairs ? minId : 1;
  st->maxPairIdentity = pairs ? maxId : 1;
}

/* ---------------------------------------------------------------- driver */

static Buf g_result;

/*
 * encFlat: concatenated encoded sequences (codes 0..4), lens: per-seq lengths.
 * Returns 0 on success, -1 empty profile, -2 DP too large.
 */
int gx_msa_run(const uint8_t *encFlat, const uint32_t *lens, uint32_t n,
               int32_t kmerSize, double match, double mismatch, double gapOpen,
               double gapExtend, double gapGapBonus) {
  buf_init(&g_result);
  if (n < 2) return -3;

  double distancesMs, alignmentMs;

  /* offsets */
  uint32_t *offs = xmalloc((n + 1) * sizeof(uint32_t));
  offs[0] = 0;
  for (uint32_t i = 0; i < n; i++) offs[i + 1] = offs[i] + lens[i];

  /* ---- pairwise k-mer distances ---- */
  double tK = gx_now();
  int32_t k = msa_k_for(lens, n, kmerSize);
  double *dist = xmalloc((size_t)n * n * sizeof(double));

  if (k < 3) { /* unreachable given msa_k_for, kept for exact parity */
    for (uint32_t i = 0; i < n; i++) {
      for (uint32_t j = i + 1; j < n; j++) {
        uint32_t Lo = lens[i] < lens[j] ? lens[i] : lens[j];
        uint64_t diff = 0;
        for (uint32_t p = 0; p < Lo; p++)
          if (encFlat[offs[i] + p] != encFlat[offs[j] + p]) diff++;
        double d = Lo ? (double)diff / (double)Lo : 1;
        dist[(size_t)i * n + j] = dist[(size_t)j * n + i] = d;
      }
    }
    k = -1;
  } else {
    int32_t size = 1 << (2 * k);
    int32_t *vecs = xmalloc((size_t)n * size * sizeof(int32_t));
    int32_t *tots = xcalloc(n, sizeof(int32_t));
    for (uint32_t i = 0; i < n; i++)
      kmer_counts(encFlat + offs[i], lens[i], k, vecs + (size_t)i * size, &tots[i]);
    for (uint32_t i = 0; i < n; i++) {
      dist[(size_t)i * n + i] = 0;
      for (uint32_t j = i + 1; j < n; j++) {
        double d = kmer_distance(vecs + (size_t)i * size, vecs + (size_t)j * size, size);
        dist[(size_t)i * n + j] = dist[(size_t)j * n + i] = d;
      }
    }
    free(vecs); free(tots);
  }
  distancesMs = gx_now() - tK;

  /* ---- guide tree ---- */
  GTNode *nodes = NULL;
  int32_t root = msa_upgma(dist, n, &nodes);

  /* ---- progressive merge (iterative post-order DFS) ---- */
  double tA = gx_now();

  /* stack frames: node id + phase; profiles stored per cluster slot */
  Profile **prof = xcalloc((size_t)n * 2 - 1 ? (size_t)n * 2 - 1 : 1, sizeof(Profile *));
  typedef struct { int32_t node; uint8_t expanded; } Frame;
  Frame *stack = xmalloc(((size_t)n * 2 + 8) * sizeof(Frame));
  int32_t sp = 0;
  stack[sp].node = root; stack[sp].expanded = 0; sp++;

  while (sp > 0) {
    Frame *fr = &stack[sp - 1];
    const GTNode *nd = &nodes[fr->node];
    if (nd->leaf) {
      prof[fr->node] = profile_single(encFlat + offs[nd->idx], lens[nd->idx], nd->idx);
      sp--;
      continue;
    }
    if (!fr->expanded) {
      fr->expanded = 1;
      stack[sp].node = nd->left; stack[sp].expanded = 0; sp++;
      stack[sp].node = nd->right; stack[sp].expanded = 0; sp++;
      continue;
    }
    Profile *LP = prof[nd->left], *RP = prof[nd->right];
    Profile *merged = NULL;
    int rc = merge_profiles(LP, RP, match, mismatch, gapOpen, gapExtend, gapGapBonus, &merged);
    profile_free(LP); free(LP);
    profile_free(RP); free(RP);
    prof[nd->left] = prof[nd->right] = NULL;
    if (rc != 0) {
      free(stack); free(prof); free(nodes); free(dist); free(offs);
      while (sp > 0) sp--;
      return rc;
    }
    prof[fr->node] = merged;
    sp--;
  }
  free(stack);

  Profile *rootProf = prof[root];
  alignmentMs = gx_now() - tA;

  /* rows reassembled in input order */
  size_t rowsSz = (size_t)n * rootProf->len;
  uint8_t *rowsByOrig = xcalloc(rowsSz ? rowsSz : 1, 1);
  for (uint32_t mi = 0; mi < rootProf->nMem; mi++) {
    memcpy(rowsByOrig + (size_t)rootProf->origIdx[mi] * rootProf->len,
           rootProf->chars + (size_t)mi * rootProf->len,
           rootProf->len);
  }

  /* ---- stats ---- */
  AlnStats st;
  double *identity = xmalloc((size_t)n * n * sizeof(double));
  double *perGap = xmalloc(n * sizeof(double));
  uint8_t *colClass = xmalloc(rootProf->len);
  uint8_t *consensus = xmalloc(rootProf->len);
  aln_stats(rootProf, rowsByOrig, &st, colClass, consensus, identity, perGap);

  /* ---- serialize blob ---- */
  {
    Buf *o = &g_result;
    uint32_t version = 1, Lout = rootProf->len, kOut = k;
    buf_append(o, "GXMR", 4);
    buf_append(o, &version, 4);
    buf_append(o, &n, 4);
    buf_append(o, &Lout, 4);
    double times[2] = { distancesMs, alignmentMs };
    buf_append(o, times, sizeof(times));
    buf_append(o, &kOut, 4);
    buf_append(o, &(uint32_t){ 0 }, 4); /* reserved */
    buf_append(o, dist, (size_t)n * n * sizeof(double));
    buf_append(o, rootProf->cols, (size_t)5 * Lout * sizeof(uint32_t));
    buf_append(o, rowsByOrig, (size_t)n * Lout);
    while (o->len % 8) buf_push(o, 0); /* f64-align the identity matrix */
    buf_append(o, identity, (size_t)n * n * sizeof(double));
    buf_append(o, consensus, Lout);
    buf_append(o, colClass, Lout);
    while (o->len % 8) buf_push(o, 0);
    buf_append(o, perGap, n * sizeof(double));
    double fs[4] = { st.gapFraction, st.meanPairIdentity, st.minPairIdentity, st.maxPairIdentity };
    buf_append(o, fs, sizeof(fs));
    uint32_t us[6] = { st.conservedColumns, st.variableColumns, st.informativeColumns,
                       st.singletonColumns, st.gapColumns, st.variableSitesTotal };
    buf_append(o, us, sizeof(us));

  }

  profile_free(rootProf); free(rootProf);
  free(rowsByOrig); free(identity); free(perGap); free(colClass); free(consensus);
  free(prof); free(nodes); free(dist); free(offs);
  return 0;
}

const uint8_t *gx_msa_result_ptr(void) { return g_result.data; }
uint32_t gx_msa_result_len(void) { return (uint32_t)g_result.len; }
void gx_msa_free(void) { buf_free(&g_result); }
