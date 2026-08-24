import { reverseComplement } from "./sequence.js";

function mirrorRange([s, e], len) { return [len - e, len - s]; }

/** The overlap boundary for `readIdx`, expressed in the coordinate system of `seqForStep` (native if strand=+1, reverse-complemented if strand=-1). */
function overlapRangeInStepCoords(ov, readIdx, len, strand) {
  let nativeRange;
  if (readIdx === ov.a) {
    nativeRange = ov.aRange; // A's position is always native in our overlap model
  } else {
    // B's range was computed in "transformed" coordinates during detection
    // (native B coords if same-orientation, reverseComplement(B) coords if
    // reverse-orientation) - convert back to B's own native coordinates first.
    nativeRange = ov.orientation === "rev" ? mirrorRange(ov.bRange, len) : ov.bRange;
  }
  return strand === 1 ? nativeRange : mirrorRange(nativeRange, len);
}

function seqForStep(reads, step) {
  const s = reads[step.readIdx].seq;
  return step.strand === 1 ? s : reverseComplement(s);
}

/**
 * Build one contig's draft sequence from its resolved read path. This is a
 * "greedy layout" consensus: each read after the first contributes only
 * its non-overlapping suffix, trimmed at the estimated (minimizer-based,
 * approximate) overlap boundary. There is no base-level multiple-sequence
 * consensus / error correction here - the output is a draft assembly, the
 * same as what a tool like miniasm produces before a polishing pass
 * (Racon/medaka) is run on it. This is disclosed to the user, not hidden.
 */
export function buildContigFromPath(reads, path) {
  const { steps } = path;
  let seq = seqForStep(reads, steps[0]);
  const members = [{ readIdx: steps[0].readIdx, strand: steps[0].strand, contigStart: 0, contigEnd: seq.length }];

  for (let i = 1; i < steps.length; i++) {
    const step = steps[i];
    const len = reads[step.readIdx].seq.length;
    const [, ovEnd] = overlapRangeInStepCoords(step.overlapWithPrev, step.readIdx, len, step.strand);
    const stepSeq = seqForStep(reads, step);
    const trimStart = Math.max(0, Math.min(ovEnd, stepSeq.length));
    const addition = stepSeq.slice(trimStart);
    const start = seq.length;
    seq += addition;
    members.push({ readIdx: step.readIdx, strand: step.strand, contigStart: start, contigEnd: seq.length, trimmedFromStart: trimStart });
  }

  let circular = path.circular;
  if (path.circular && path.closingOverlap) {
    const lastStep = steps[steps.length - 1];
    const lastLen = reads[lastStep.readIdx].seq.length;
    // Position (within the last read's own sequence) where it starts
    // overlapping back with the contig's start.
    const [ovStart] = overlapRangeInStepCoords(path.closingOverlap, lastStep.readIdx, lastLen, lastStep.strand);
    const lastMember = members[members.length - 1];
    // Convert that read-local position into a contig position: contigStart
    // corresponds to the read's local position `trimmedFromStart` (where
    // its non-overlapping contribution began), not local position 0.
    const trimTo = lastMember.contigStart + (ovStart - (lastMember.trimmedFromStart ?? 0));
    if (trimTo > lastMember.contigStart && trimTo < seq.length) {
      seq = seq.slice(0, trimTo);
      lastMember.contigEnd = trimTo;
    } else {
      circular = false; // closure trim didn't make sense - report as linear rather than guess
    }
  }

  return { seq, length: seq.length, members, circular, readCount: steps.length };
}

export function buildAllContigs(reads, paths) {
  return paths
    .map((p) => buildContigFromPath(reads, p))
    .sort((a, b) => b.length - a.length);
}
