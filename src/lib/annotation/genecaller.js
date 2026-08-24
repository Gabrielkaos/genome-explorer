/**
 * De novo bacterial CDS prediction, browser-side.
 *
 * Architecture follows the same broad strategy as Prodigal's anonymous mode,
 * reduced to what can run comfortably client-side:
 *   1. Unsupervised training: the longest stop-to-stop ORFs on both strands
 *      are treated as presumptive coding sequence and used to build a
 *      hexamer (dicodon) log-likelihood-ratio model vs. genomic background.
 *   2. Six-frame scanning for stop-to-stop ORFs with bacterial start codons
 *      (ATG / GTG / TTG), including partial ORFs truncated by contig edges.
 *   3. Start-codon selection per ORF by combining coding potential
 *      (hexamer log-ratio over the candidate's span), Shine-Dalgarno RBS
 *      strength in the upstream region, a start-type prior learned from the
 *      training set, and a mild length prior.
 *   4. Overlap resolution across frames/strands with strand-aware tolerance
 *      (bacterial operons legitimately allow short same-strand overlaps).
 *
 * It deliberately does NOT call RNA features (rRNA/tRNA need covariance-model
 * search like Infernal/tRNAscan-SE/Barrnap) and assigns no product names from
 * reference databases (needs BLAST/DIAMOND vs UniProt etc.) — everything it
 * calls "hypothetical" honestly is hypothetical until similarity evidence
 * exists. Those limits are surfaced in the UI rather than papered over.
 */

import { reverseComplement } from "../assembly/sequence.js";
import { translate, molecularWeight, theoreticalPI, gravy, aromaticity, findTmSegments, predictSignalPeptide, localizationHint } from "./protein.js";

const START_CODONS = ["ATG", "GTG", "TTG"];
const STOPS = new Set(["TAA", "TAG", "TGA"]);
const SD_CONSENSUS = "TAAGGAGGT";

/* ------------------------------------------------------------------ */
/* Six-frame ORF scanning                                              */
/* ------------------------------------------------------------------ */

/**
 * Finds all stop-to-stop ORF intervals >= minLen on one oriented sequence.
 * Returns intervals in oriented coordinates (end exclusive):
 *   partial5 — no in-frame stop upstream, i.e. the gene's start region runs
 *              off the contig edge and no complete start context exists
 *   partial3 — no downstream stop before the contig edge (stop codon missing)
 */
export function findOrfIntervals(seq, minLen) {
  const out = [];
  const n = seq.length;
  for (let frame = 0; frame < 3; frame++) {
    let orfStart = frame;
    let seenStop = false;
    for (let pos = frame; pos + 3 <= n; pos += 3) {
      if (!STOPS.has(seq.slice(pos, pos + 3))) continue;
      const end = pos + 3;
      if (end - orfStart >= minLen) out.push({ start: orfStart, end, partial5: !seenStop, partial3: false });
      orfStart = end;
      seenStop = true;
    }
    const lastFullCodonEnd = n - ((n - orfStart) % 3);
    if (lastFullCodonEnd - orfStart >= minLen) {
      out.push({ start: orfStart, end: lastFullCodonEnd, partial5: !seenStop, partial3: true });
    }
  }
  return out;
}

/** Candidate start codons inside an ORF (ORF-relative offsets, step 3). */
function findCandidateStarts(orfSeq, allowedStarts) {
  const wanted = new Set(allowedStarts);
  const cands = [];
  for (let pos = 0; pos + 3 <= orfSeq.length; pos += 3) {
    const c = orfSeq.slice(pos, pos + 3);
    if (wanted.has(c)) cands.push({ rel: pos, codon: c });
  }
  return cands;
}

/* ------------------------------------------------------------------ */
/* Coding model                                                        */
/* ------------------------------------------------------------------ */

function hexIndex(h) {
  let v = 0;
  for (let i = 0; i < 6; i++) {
    const c = h.charCodeAt(i);
    const b2 = c === 65 ? 0 : c === 67 ? 1 : c === 71 ? 2 : c === 84 ? 3 : -1;
    if (b2 < 0) return -1;
    v = v * 4 + b2;
  }
  return v;
}

/**
 * Builds the hexamer log-ratio table from training ORFs (longest ones),
 * against whole-genome hexamer frequencies as background. Returns
 * Float64Array(4096) of log2 odds (0 for unseen/ambiguous entries means
 * neutral, not informative).
 */
export function trainCodingModel(seqs, { maxTrainingGenes = 600 } = {}) {
  const bgCounts = new Float64Array(4096);
  for (const s of seqs) {
    for (let i = 0; i + 6 <= s.length; i++) {
      const idx = hexIndex(s.slice(i, i + 6));
      if (idx >= 0) bgCounts[idx]++;
    }
  }

  // Training set: longest ORFs across all contigs, both strands. On small
  // inputs the 300nt floor yields too few examples, so relax adaptively.
  const collectOrfs = (minLen) => {
    const arr = [];
    for (const s of seqs) {
      for (const strandSeq of [s, reverseComplement(s)]) {
        for (const orf of findOrfIntervals(strandSeq, minLen)) {
          arr.push({ seq: strandSeq.slice(orf.start, orf.end), len: orf.end - orf.start });
        }
      }
    }
    return arr;
  };
  let trainOrfs = collectOrfs(300);
  let minUsed = 300;
  if (trainOrfs.length < 40) { trainOrfs = collectOrfs(180); minUsed = 180; }
  if (trainOrfs.length < 20) { trainOrfs = collectOrfs(120); minUsed = 120; }
  trainOrfs.sort((a, b) => b.len - a.len);
  const chosen = trainOrfs.slice(0, maxTrainingGenes);

  const codingCounts = new Float64Array(4096);
  let codingHexamers = 0;
  for (const t of chosen) {
    for (let i = 0; i + 6 <= t.seq.length; i += 3) {
      const idx = hexIndex(t.seq.slice(i, i + 6));
      if (idx >= 0) { codingCounts[idx]++; codingHexamers++; }
    }
  }

  const bgTotal = bgCounts.reduce((a, b) => a + b, 0);
  const pseudo = 1.0;
  const CLAMP = 3.0; // cap per-hexamer evidence so an undertrained table can't dominate
  const logOdds = new Float64Array(4096);
  for (let i = 0; i < 4096; i++) {
    if (bgCounts[i] === 0 && codingCounts[i] === 0) continue;
    const pC = (codingCounts[i] + pseudo) / (codingHexamers + 4096 * pseudo);
    const pB = (bgCounts[i] + pseudo) / (bgTotal + 4096 * pseudo);
    logOdds[i] = Math.max(-CLAMP, Math.min(CLAMP, Math.log2(pC / pB)));
  }

  return {
    logOdds,
    trained: chosen.length >= 30,
    trainingGenes: chosen.length,
    minTrainingLen: minUsed,
    startUsage: learnStartUsage(chosen),
  };
}

/** Presume the true start of each training ORF is its strongest-RBS candidate; tally usage. */
function learnStartUsage(trainingOrfs) {
  const counts = { ATG: 0, GTG: 0, TTG: 0 };
  let usable = 0;
  for (const t of trainingOrfs) {
    let best = null;
    for (let pos = 0; pos + 3 <= t.seq.length - 3; pos += 3) {
      const codon = t.seq.slice(pos, pos + 3);
      if (!START_CODONS.includes(codon)) continue;
      const rbs = scoreRBS(t.seq, pos).score;
      const cand = { codon, rbs };
      if (!best || cand.rbs > best.rbs || (cand.rbs === best.rbs && codon === "ATG")) best = cand;
    }
    if (best) { counts[best.codon]++; usable++; }
  }
  if (usable < 20) return { ATG: 0.42, GTG: 0.34, TTG: 0.24 }; // literature-typical bacterial mix fallback
  const smoothed = {};
  for (const c of START_CODONS) smoothed[c] = 0.8 * (counts[c] / usable) + 0.2 * ({ ATG: 0.42, GTG: 0.34, TTG: 0.24 })[c];
  return smoothed;
}

/* ------------------------------------------------------------------ */
/* Ribosome-binding-site scoring                                       */
/* ------------------------------------------------------------------ */

/**
 * Best Shine-Dalgarno-like match upstream of a candidate start (oriented
 * coordinates). The consensus TAAGGAGGT is aligned at every offset whose
 * full 9-mer fits in [startPos-22, startPos); spacer = distance from motif
 * end to the start codon, with a mild preference for the biological ~4-12 nt
 * range. Score normalized to [0,1]; ~0.6+ indicates a convincing RBS.
 */
export function scoreRBS(seq, startPos) {
  let bestMatches = 0, bestOffset = -1;
  const firstOff = Math.max(0, startPos - 22);
  for (let off = firstOff; off + SD_CONSENSUS.length <= startPos; off++) {
    let m = 0;
    for (let j = 0; j < SD_CONSENSUS.length; j++) {
      if (seq[off + j] === SD_CONSENSUS[j]) m++;
    }
    if (m > bestMatches) { bestMatches = m; bestOffset = off; }
  }

  let score = 0, spacer = null;
  if (bestOffset >= 0) {
    spacer = startPos - (bestOffset + SD_CONSENSUS.length);
    score = Math.max(0, Math.min(1, (bestMatches - 4) / 5));
    if (spacer < 3) score *= 0.55;
    else if (spacer > 14) score *= 0.6;
    else if (spacer > 12) score *= 0.85;
  }
  const strong = bestMatches >= 7 && spacer !== null && spacer >= 4 && spacer <= 13;
  return { score: +score.toFixed(3), matches: bestMatches, spacer, strong };
}

/* ------------------------------------------------------------------ */
/* Start selection                                                     */
/* ------------------------------------------------------------------ */

const sigmoid = (x) => 1 / (1 + Math.exp(-x));

/**
 * Picks the best start codon for one ORF. All offsets here are RELATIVE to
 * the ORF slice (orfSeq), which keeps prefix-sum indexing and RBS windowing
 * consistent; callers map `rel` back to genomic coordinates.
 *
 * Weights adapt to model confidence: when the hexamer table trained on few
 * ORFs (small contigs), coding statistics are noisy, so Shine-Dalgarno
 * evidence and the start-type prior carry more weight — mirroring how a
 * careful human annotator leans on conserved context for draft genomes.
 */
function chooseStartForOrf(orfSeq, logOdds, startUsage, useRbs, trained, allowedStarts = START_CODONS) {
  let cum = 0;
  const cumAt = new Float64Array(Math.floor(orfSeq.length / 3) + 1);
  let ci = 1;
  for (let pos = 0; pos + 6 <= orfSeq.length; pos += 3) {
    const idx = hexIndex(orfSeq.slice(pos, pos + 6));
    cum += idx >= 0 ? logOdds[idx] : 0;
    cumAt[ci++] = cum;
  }
  const totalCells = ci - 1;

  // Coding statistics are informative but noisy, especially when trained on
  // few ORFs; a strong Shine-Dalgarno motif at the correct spacer distance is
  // rarer and more decisive evidence for a true start (this mirrors how
  // Prodigal's start model treats SD strength as a first-class feature).
  const wCoding = trained ? 0.9 : 0.6;
  const wRbs = useRbs ? (trained ? 1.0 : 1.2) : 0;
  const strongBonus = useRbs ? (trained ? 0.3 : 0.4) : 0;

  let best = null;
  for (const cand of findCandidateStarts(orfSeq, allowedStarts)) {
    const rel = cand.rel;
    const cellsFrom = Math.floor(rel / 3);
    const denom = Math.max(1, totalCells - cellsFrom);
    const meanLR = totalCells > cellsFrom ? (cum - cumAt[cellsFrom]) / denom : 0;
    const codingTerm = sigmoid(meanLR * 1.3); // squashes log-ratio scale into 0..1
    const rbs = useRbs ? scoreRBS(orfSeq, rel) : { score: 0.35, matches: 0, spacer: null, strong: false }; // neutral prior when disabled
    const prior = startUsage[cand.codon] ?? 0.33;
    const lenAhead = orfSeq.length - rel;
    const lenBonus = 0.3 * Math.min(1, lenAhead / 1200);

    const total =
      wCoding * codingTerm +
      wRbs * rbs.score +
      (rbs.strong ? strongBonus : 0) +
      0.5 * prior +
      lenBonus;

    if (!best || total > best.total) {
      best = { rel, codon: cand.codon, total, codingTerm, rbs, meanLR };
    }
  }
  return best;
}

/* ------------------------------------------------------------------ */
/* Overlap resolution                                                  */
/* ------------------------------------------------------------------ */

/**
 * Greedy resolution: accept highest-scoring genes first; reject a gene whose
 * overlap with any accepted gene exceeds the strand-aware allowance.
 */
export function resolveOverlaps(genes, { oppositeStrandMax = 25, sameStrandMax = 60 } = {}) {
  const sorted = [...genes].sort((a, b) => b.score - a.score);
  const accepted = [];
  const rejected = [];
  for (const g of sorted) {
    let ok = true;
    for (const acc of accepted) {
      if (acc.contigId !== g.contigId) continue;
      const ov = Math.min(acc.end, g.end) - Math.max(acc.start, g.start);
      if (ov <= 0) continue;
      const limit = acc.strand === g.strand ? sameStrandMax : oppositeStrandMax;
      if (ov > limit) { ok = false; break; }
    }
    (ok ? accepted : rejected).push(g);
  }
  return { accepted, rejectedCount: rejected.length };
}

/* ------------------------------------------------------------------ */
/* Genome-wide orchestration                                           */
/* ------------------------------------------------------------------ */

export function annotateGenome(contigs, params, reportProgress = () => {}) {
  const t0 = performance.now();
  const { minGeneLen = 90 } = params;

  const clean = contigs
    .map((c) => ({ ...c, seq: (c.seq || "").toUpperCase().replace(/[^ACGTN]/g, "N") }))
    .filter((c) => c.seq.length >= minGeneLen + 60);
  const skippedTiny = contigs.length - clean.length;

  reportProgress({ stage: "training", pct: 8 });
  const model = trainCodingModel(clean.map((c) => c.seq));
  const allowed = (params.allowedStarts?.length ? params.allowedStarts : [...START_CODONS]).filter((c) => START_CODONS.includes(c));

  reportProgress({ stage: "scanning", pct: 30 });
  const rawGenes = [];
  clean.forEach((contig, ci) => {
    const L = contig.seq.length;
    // forward strand
    for (const orf of findOrfIntervals(contig.seq, minGeneLen)) {
      collectGene(rawGenes, contig, ci, contig.seq, orf, "+", L, model, params, minGeneLen, allowed);
    }
    // reverse strand (work in revcomp coords, map back)
    const rc = reverseComplement(contig.seq);
    for (const orf of findOrfIntervals(rc, minGeneLen)) {
      collectGene(rawGenes, contig, ci, rc, orf, "-", L, model, params, minGeneLen, allowed);
    }
    reportProgress({ stage: "scanning", pct: 30 + Math.round(((ci + 1) / clean.length) * 35) });
  });

  reportProgress({ stage: "resolving", pct: 72 });
  const { accepted, rejectedCount } = resolveOverlaps(rawGenes);

  // deterministic locus tags ordered by contig then coordinate
  accepted.sort((a, b) => (a.contigIdx - b.contigIdx) || (a.start - b.start));
  accepted.forEach((g, i) => {
    g.locusTag = `${params.locusPrefix ?? "GE"}_${String(i + 1).padStart(5, "0")}`;
  });

  reportProgress({ stage: "proteins", pct: 82 });
  for (const g of accepted) finalizeProtein(g);

  const stats = buildStats(accepted, clean, rejectedCount, skippedTiny, model, performance.now() - t0);
  const gcTracks = clean.map((c) => gcTrack(c));

  reportProgress({ stage: "done", pct: 100 });
  return { genes: accepted, stats, gcTracks, modelInfo: { trained: model.trained, trainingGenes: model.trainingGenes, startUsage: model.startUsage } };
}

function collectGene(out, contig, contigIdx, orientedSeq, orf, strand, contigLen, model, params, minGeneLen, allowedStarts) {
  const orfSeq = orientedSeq.slice(orf.start, orf.end);
  const best = chooseStartForOrf(orfSeq, model.logOdds, model.startUsage, params.useRbs ?? true, model.trained, allowedStarts);
  if (!best || best.rel < 3) return; // no usable start codon inside the ORF

  const geneLen = orfSeq.length - best.rel;
  if (geneLen < minGeneLen) return;

  // Short-ORF evidence gate: random sequence produces many short open frames,
  // so ORFs under 240 nt need independent support (a convincing RBS, or an
  // ATG start with strong coding statistics). Long ORFs pass unconditionally —
  // chance stop-to-stop spans of that size are rare enough to be worth
  // reporting. Edge-truncated genes are exempt from the RBS requirement since
  // their upstream context may simply be missing.
  if (geneLen < 240 && (params.useRbs ?? true) && !orf.partial3) {
    const hasEvidence = best.rbs.score >= 0.45 || (best.codon === "ATG" && best.codingTerm >= 0.78);
    if (!hasEvidence) return;
  }

  // map oriented -> forward coordinates (best.rel is ORF-relative)
  const absStart = orf.start + best.rel;
  let fStart, fEnd;
  if (strand === "+") { fStart = absStart; fEnd = orf.end; }
  else { fStart = contigLen - orf.end; fEnd = contigLen - absStart; }

  const partial5 = !!(orf.partial5 && best.rel === 0);
  const partial3 = !!orf.partial3;

  out.push({
    contigId: contig.id,
    contigIdx,
    strand,
    start: fStart, // 0-based half-open, converted to 1-based closed at export/render time
    end: fEnd,
    score: +best.total.toFixed(3),
    codingTerm: +best.codingTerm.toFixed(3),
    meanHexLR: +best.meanLR.toFixed(3),
    rbsScore: best.rbs.score,
    rbsMatches: best.rbs.matches,
    rbsSpacer: best.rbs.spacer,
    startCodon: best.codon,
    partial: partial5 && partial3 ? "both" : partial5 ? "5'" : partial3 ? "3'" : null,
    _orientedSeq: orfSeq.slice(best.rel),
  });
}

function finalizeProtein(gene) {
  const dna = gene._orientedSeq;
  const prot = translate(dna);
  gene.dnaSeq = dna;
  gene.protSeq = prot;
  gene.lengthNt = dna.length;
  gene.lengthAa = prot.length;
  gene.mw = +molecularWeight(prot).toFixed(1);
  const pi = theoreticalPI(prot);
  gene.pi = pi !== null ? +pi.toFixed(2) : null;
  gene.gravy = +gravy(prot).toFixed(3);
  gene.aromaticity = +aromaticity(prot).toFixed(3);
  gene.tmSegments = findTmSegments(prot);
  gene.signal = predictSignalPeptide(prot);
  gene.localization = localizationHint(gene.tmSegments, gene.signal);
  let gc = 0;
  for (const b of dna) if (b === "G" || b === "C") gc++;
  gene.gcContent = dna.length ? +(gc / dna.length * 100).toFixed(1) : 0;
  gene.product = gene.partial ? "hypothetical protein (partial)" : "hypothetical protein";
  delete gene._orientedSeq;
}

function buildStats(genes, contigs, rejectedCount, skippedTiny, model, ms) {
  const totalBases = contigs.reduce((a, c) => a + c.seq.length, 0);
  const codingBases = genes.reduce((a, g) => a + (g.end - g.start), 0);
  const starts = { ATG: 0, GTG: 0, TTG: 0 };
  let partials = 0, withSignal = 0, withTm = 0;
  for (const g of genes) {
    starts[g.startCodon]++;
    if (g.partial) partials++;
    if (g.signal) withSignal++;
    if (g.tmSegments.length) withTm++;
  }
  const lens = genes.map((g) => g.lengthAa).sort((a, b) => a - b);
  const q = (p) => (lens.length ? lens[Math.min(lens.length - 1, Math.floor(p * lens.length))] : 0);
  return {
    numContigs: contigs.length,
    skippedTinyContigs: skippedTiny,
    totalBases,
    numGenes: genes.length,
    codingDensity: totalBases ? codingBases / totalBases : 0,
    meanGeneLenAa: genes.length ? Math.round(lens.reduce((a, b) => a + b, 0) / genes.length) : 0,
    medianGeneLenAa: q(0.5),
    shortestAa: lens[0] ?? 0,
    longestAa: lens[lens.length - 1] ?? 0,
    starts,
    partialGenes: partials,
    signalPeptides: withSignal,
    tmProteins: withTm,
    resolvedOverlaps: rejectedCount,
    modelTrained: model.trained,
    trainingGenes: model.trainingGenes,
    computeTimeMs: ms,
  };
}

function gcTrack(contig) {
  const win = Math.max(200, Math.floor(contig.seq.length / 400));
  const points = [];
  let gc = 0;
  const seq = contig.seq;
  for (let i = 0; i < seq.length; i++) {
    if (seq[i] === "G" || seq[i] === "C") gc++;
    if (i >= win && (seq[i - win] === "G" || seq[i - win] === "C")) gc--;
    if (i >= win - 1 && (i - win + 1) % Math.floor(win / 2) === 0) {
      points.push({ pos: i - win + 1 + Math.floor(win / 2), gc: +(gc / win * 100).toFixed(2) });
    }
  }
  return { contigId: contig.id, window: win, points };
}
