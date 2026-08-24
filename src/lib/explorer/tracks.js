/**
 * Per-contig precomputed tracks for the genome browser: sliding-window GC
 * content and GC skew ((G-C)/(G+C)), plus cumulative-skew minima as a rough
 * origin-of-replication hint (linear chromosomes/plasmids excluded).
 */

export function gcContent(seq, from = 0, to = seq.length) {
  let gc = 0, n = 0;
  for (let i = from; i < to && i < seq.length; i++) {
    const b = seq[i];
    if (b === "G" || b === "C" || b === "g" || b === "c") gc++;
    if (b !== "-" && b) n++;
  }
  return n ? gc / n : 0;
}

/**
 * Windowed tracks. window ~ max(32, floor(len/1200)) keeps arrays small
 * enough to redraw cheaply while staying smooth when zoomed in.
 */
export function buildTracks(seq) {
  const len = seq.length;
  const win = Math.max(32, Math.floor(len / 1400));
  const step = Math.max(1, Math.floor(win / 2));
  const starts = [];
  const gc = [];
  const skew = [];
  for (let s = 0; s + win <= len; s += step) {
    let g = 0, c = 0, n = 0;
    const e = Math.min(s + win, len);
    for (let i = s; i < e; i++) {
      const b = seq[i];
      if (b === "G") g++;
      else if (b === "C") c++;
      if (b === "A" || b === "T" || b === "G" || b === "C") n++;
    }
    if (n < win * 0.5) continue; // skip windows dominated by N
    starts.push(s);
    gc.push(n ? (g + c) / n : 0);
    skew.push(g + c ? (g - c) / (g + c) : 0);
  }
  return { window: win, step, starts, gc, skew };
}

/** Mean GC over an arbitrary region using the raw sequence. */
export function regionGC(seq, start0, end) {
  return gcContent(seq, Math.max(0, start0), Math.min(seq.length, end));
}

/**
 * Rough ori/ter estimate from cumulative GC-skew extrema (DoriC-style trick).
 * Only meaningful for circular replicons; returns null otherwise.
 */
export function skewExtrema(tracks) {
  if (!tracks.skew.length) return null;
  let cum = 0;
  let minV = Infinity, maxV = -Infinity;
  let minI = 0, maxI = 0;
  for (let i = 0; i < tracks.skew.length; i++) {
    cum += tracks.skew[i];
    if (cum < minV) { minV = cum; minI = i; }
    if (cum > maxV) { maxV = cum; maxI = i; }
  }
  return { originPos: tracks.starts[minI], terminusPos: tracks.starts[maxI] };
}
