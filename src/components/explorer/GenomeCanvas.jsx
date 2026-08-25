import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { regionGC } from "../../lib/explorer/tracks.js";

const TYPE_COLORS = {
  CDS: "#ef7fa3",
  tRNA: "#e8c15a",
  rRNA: "#e8a95a",
  tmRNA: "#e8c15a",
  ncRNA: "#c9b3f5",
  misc_RNA: "#c9b3f5",
  mobile_genetic_element: "#a08cf0",
  pseudogene: "#93a0ae",
};
const BASE_COLORS = { A: "#68c98f", T: "#e6685f", G: "#e8c15a", C: "#5ec8d8" };

const RULER_H = 26;
const FWD_TOP = 34, BAND_H = 52, AXIS_Y = FWD_TOP + BAND_H + 12; // axis center line
const REV_TOP = AXIS_Y + 6;
const GC_TOP = REV_TOP + BAND_H + 14;
const GC_H = 44;
const SKEW_TOP = GC_TOP + GC_H + 12;
const SKEW_H = 36;
const SEQ_TOP = SKEW_TOP + SKEW_H + 14;

/**
 * Canvas-based genome browser for one contig:
 * ruler / forward & reverse feature bands / GC content / GC skew /
 * base-level sequence strip when zoomed in. Wheel zooms around the cursor,
 * drag pans (or rubber-band selects), click selects features.
 */
export default function GenomeCanvas({
  contig, features, tracks, oriPos,
  view, onViewChange,
  selection, onSelectionChange,
  mode = "pan",
  selectedFeatureId, onSelectFeature,
  motifHits = [],
}) {
  const wrapRef = useRef(null);
  const canvasRef = useRef(null);
  const miniRef = useRef(null);
  const dragRef = useRef(null);
  const [size, setSize] = useState({ w: 900, h: SEQ_TOP + 34 });
  const [hover, setHover] = useState({ x: null });

  const seq = contig.seq;
  const len = seq.length;

  /* ------------------------------ sizing ------------------------------ */
  useEffect(() => {
    const el = wrapRef.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => {
      const w = Math.max(320, Math.floor(entries[0].contentRect.width));
      setSize((s) => ({ ...s, w }));
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  /* --------------------------- interactions --------------------------- */
  const gOf = useCallback((x) => {
    const span = view.end - view.start;
    return view.start + (x / size.w) * span;
  }, [view, size.w]);

  const clampView = useCallback((start, end) => {
    let span = end - start;
    span = Math.max(12, Math.min(len, span));
    let s = Math.max(0, Math.min(start, len - span));
    return { start: s, end: s + span };
  }, [len]);

  function handleWheel(e) {
    e.preventDefault();
    const rect = canvasRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const anchor = gOf(x);
    const factor = e.deltaY > 0 ? 1.18 : 1 / 1.18;
    const span = (view.end - view.start) * factor;
    const frac = (anchor - view.start) / (view.end - view.start);
    const nv = clampView(anchor - frac * span, anchor + (1 - frac) * span);
    onViewChange(nv);
  }

  function handleMouseDown(e) {
    const rect = canvasRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    dragRef.current = { x0: x, moved: false };
  }

  function handleMouseMove(e) {
    const rect = canvasRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    setHover({ x });
    const d = dragRef.current;
    if (!d) return;
    const dx = x - d.x0;
    if (!d.moved && Math.abs(dx) < 3) return;
    d.moved = true;
    if (mode === "pan") {
      const span = view.end - view.start;
      if (d.baseStart === undefined) d.baseStart = view.start;
      const targetStart = d.baseStart - ((x - d.x0) / size.w) * span;
      onViewChange(clampView(targetStart, targetStart + span));
    } else {
      const a = Math.min(gOf(d.x0), gOf(x));
      const b = Math.max(gOf(d.x0), gOf(x));
      onSelectionChange({ start0: Math.floor(a), end: Math.ceil(b) });
    }
  }

  function handleMouseUp(e) {
    const d = dragRef.current;
    dragRef.current = null;
    if (!d || d.moved) return;
    const rect = canvasRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    // Click: hit-test features first.
    const hit = hitTest(x, y);
    onSelectFeature(hit?.locusTag ?? null);
  }

  function handleDoubleClick(e) {
    const rect = canvasRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    const hit = hitTest(x, y);
    if (hit) {
      const pad = Math.max(40, (hit.end - hit.start) * 0.15);
      onViewChange(clampView(hit.start - pad, hit.end + pad));
    }
  }

  function handleMouseLeave() {
    setHover({ x: null });
    dragRef.current = null;
  }

  const hitTest = useCallback((x, y) => {
    const g = gOf(x);
    // Features are sorted by start; scan candidates whose span could cover g.
    let lo = 0, hi = features.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (features[mid].end <= g) lo = mid + 1; else hi = mid; }
    for (let i = lo; i >= 0 && i < features.length && features[i].start <= g; i--) {
      const f = features[i];
      if (g >= f.start && g <= f.end) {
        const inFwd = f.strand !== "-" && y >= FWD_TOP - 4 && y <= FWD_TOP + BAND_H;
        const inRev = f.strand === "-" && y >= REV_TOP - 4 && y <= REV_TOP + BAND_H;
        if (inFwd || inRev) return f;
      }
      if (i === 0) break;
    }
    return null;
  }, [features, gOf]);

  /* ------------------------------ drawing ------------------------------ */
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !len) return;
    const dpr = window.devicePixelRatio || 1;
    const W = size.w, H = size.h;
    canvas.width = W * dpr; canvas.height = H * dpr;
    canvas.style.width = `${W}px`; canvas.style.height = `${H}px`;
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);

    const span = view.end - view.start;
    const xOf = (g) => ((g - view.start) / span) * W;

    // ---- background bands
    ctx.fillStyle = "rgba(255,255,255,0.02)";
    ctx.fillRect(0, FWD_TOP - 4, W, BAND_H + 4);
    ctx.fillRect(0, REV_TOP - 4, W, BAND_H + 4);
    ctx.fillStyle = "rgba(94,200,216,0.03)";
    ctx.fillRect(0, GC_TOP - 4, W, GC_H + 4);
    ctx.fillStyle = "rgba(160,140,240,0.03)";
    ctx.fillRect(0, SKEW_TOP - 4, W, SKEW_H + 4);

    // ---- ruler
    ctx.font = `10.5px ${FONT_DISPLAY}`;
    ctx.textBaseline = "middle";
    const rawStep = span / (W / 90);
    const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
    const steps = [1, 2, 5, 10].map((m) => m * mag).filter((st) => st >= rawStep / 2);
    const step = steps[0] ?? mag * 10;
    const fmtTick = (v) => (v >= 10000 ? `${(v / 1000).toFixed(v % 1000 ? 1 : 0)} kb` : String(Math.round(v)));
    for (let t = Math.ceil(view.start / step) * step; t <= view.end; t += step) {
      const x = xOf(t);
      ctx.strokeStyle = "rgba(255,255,255,0.16)";
      ctx.beginPath(); ctx.moveTo(x, RULER_H - 6); ctx.lineTo(x, H); ctx.stroke();
      ctx.fillStyle = C.textDim;
      ctx.textAlign = "center";
      ctx.fillText(fmtTick(t), x, RULER_H / 2 - 2);
    }
    // Axis line
    ctx.strokeStyle = "rgba(255,255,255,0.22)";
    ctx.beginPath(); ctx.moveTo(0, AXIS_Y); ctx.lineTo(W, AXIS_Y); ctx.stroke();

    // ---- features
    ctx.font = `10px ${FONT_DISPLAY}`;
    for (const f of features) {
      if (f.end < view.start || f.start > view.end) continue;
      const x1 = xOf(f.start), x2 = xOf(f.end);
      const wPx = Math.max(2.5, x2 - x1);
      const fwd = f.strand !== "-";
      const color = TYPE_COLORS[f.type] ?? "#93a0ae";
      const sel = selectedFeatureId === f.locusTag;
      const h = Math.min(BAND_H - 8, 26);
      const yTop = fwd ? AXIS_Y - 8 - h : AXIS_Y + 8;

      ctx.fillStyle = color;
      ctx.globalAlpha = sel ? 1 : 0.82;
      roundBar(ctx, x1, yTop, wPx, h, fwd, color);
      ctx.globalAlpha = 1;
      if (sel) {
        ctx.strokeStyle = "#ffffff";
        ctx.lineWidth = 1.6;
        roundBarStroke(ctx, x1, yTop, wPx, h, fwd);
        ctx.lineWidth = 1;
      }
      if (wPx > 52) {
        ctx.fillStyle = "#05070a";
        ctx.textAlign = "left";
        const label = f.locusTag.length > Math.floor(wPx / 6.2) ? f.locusTag.slice(0, Math.max(3, Math.floor(wPx / 6.2) - 1)) + "…" : f.locusTag;
        ctx.fillText(label, x1 + 4, yTop + h / 2 + 1);
      }
    }

    // ---- origin marker (cumulative-skew minimum)
    if (Number.isFinite(oriPos)) {
      const x = xOf(oriPos);
      if (x >= 0 && x <= W) {
        ctx.strokeStyle = "#68c98f";
        ctx.setLineDash([4, 4]);
        ctx.beginPath(); ctx.moveTo(x, RULER_H); ctx.lineTo(x, SEQ_TOP - 4); ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = "#68c98f";
        ctx.textAlign = x > W - 70 ? "right" : "left";
        ctx.fillText("ori?", x + (x > W - 70 ? -4 : 4), RULER_H + 8);
      }
    }

    // ---- GC content track
    ctx.strokeStyle = "#5ec8d8";
    ctx.lineWidth = 1.3;
    ctx.beginPath();
    const meanGC = meanGCOf(seq);
    let started = false;
    for (let i = 0; i < tracks.starts.length; i++) {
      const gMid = tracks.starts[i] + tracks.window / 2;
      if (gMid < view.start || gMid > view.end) continue;
      const x = xOf(gMid);
      const y = GC_TOP + GC_H - (tracks.gc[i]) * GC_H;
      if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.strokeStyle = "rgba(255,255,255,0.18)";
    ctx.setLineDash([3, 4]);
    ctx.beginPath(); ctx.moveTo(0, GC_TOP + GC_H - meanGC * GC_H); ctx.lineTo(W, GC_TOP + GC_H - meanGC * GC_H); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = C.textFaint;
    ctx.textAlign = "left";
    ctx.fillText(`GC% (win ${tracks.window.toLocaleString()} bp · mean ${(meanGC * 100).toFixed(1)}%)`, 6, GC_TOP - 8);

    // ---- GC skew track
    ctx.strokeStyle = "#a08cf0";
    ctx.beginPath();
    started = false;
    for (let i = 0; i < tracks.skew.length; i++) {
      const gMid = tracks.starts[i] + tracks.window / 2;
      if (gMid < view.start || gMid > view.end) continue;
      const x = xOf(gMid);
      const y = SKEW_TOP + SKEW_H / 2 - tracks.skew[i] * (SKEW_H / 2 - 3) * 1.6;
      if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
    }
    ctx.stroke();
    ctx.strokeStyle = "rgba(255,255,255,0.18)";
    ctx.beginPath(); ctx.moveTo(0, SKEW_TOP + SKEW_H / 2); ctx.lineTo(W, SKEW_TOP + SKEW_H / 2); ctx.stroke();
    ctx.fillStyle = C.textFaint;
    ctx.fillText("GC skew", 6, SKEW_TOP - 6);

    // ---- sequence strip (only when bases stay legible)
    const seqSpanLimit = Math.max(30, Math.floor(W / 5));
    if (span <= seqSpanLimit) {
      ctx.font = `12px ${FONT_DISPLAY}`;
      ctx.textAlign = "center";
      const bps = Array.from({ length: span }, (_, k) => Math.floor(view.start) + k);
      for (const g of bps) {
        if (g < 0 || g >= len) continue;
        const x = xOf(g + 0.5);
        const b = seq[g];
        ctx.fillStyle = BASE_COLORS[b] ?? C.textDim;
        ctx.fillText(b, x, SEQ_TOP + 10);
        if (span <= 40) ctx.fillRect(x - 3.5, SEQ_TOP + 17, 7, 2.5);
      }
      ctx.textAlign = "left";
      ctx.fillStyle = C.textFaint;
      ctx.font = `10.5px ${FONT_DISPLAY}`;
      ctx.fillText("sequence (forward strand)", 6, SEQ_TOP - 6);
    }

    // ---- selection overlay
    if (selection && selection.end > selection.start0) {
      const x1 = xOf(selection.start0), x2 = xOf(selection.end);
      if (x2 > 0 && x1 < W) {
        ctx.fillStyle = "rgba(232,193,90,0.13)";
        ctx.fillRect(Math.max(0, x1), RULER_H, Math.min(W, x2) - Math.max(0, x1), SEQ_TOP + 26 - RULER_H);
        ctx.strokeStyle = "rgba(232,193,90,0.75)";
        for (const bx of [x1, x2]) {
          if (bx < 0 || bx > W) continue;
          ctx.beginPath(); ctx.moveTo(bx, RULER_H); ctx.lineTo(bx, SEQ_TOP + 26); ctx.stroke();
        }
      }
    }

    // ---- motif-hit markers
    ctx.fillStyle = "#e6685f";
    for (const h of motifHits) {
      const x = xOf(h.start0 + (h.end - h.start0) / 2);
      if (x < 0 || x > W) continue;
      ctx.beginPath();
      ctx.moveTo(x, RULER_H - 4);
      ctx.lineTo(x - 3.5, RULER_H + 2);
      ctx.lineTo(x + 3.5, RULER_H + 2);
      ctx.closePath();
      ctx.fill();
    }

    if (hover.x != null && !dragRef.current?.moved) {
      ctx.strokeStyle = "rgba(231,237,243,0.35)";
      ctx.setLineDash([2, 3]);
      ctx.beginPath(); ctx.moveTo(hover.x, RULER_H - 6); ctx.lineTo(hover.x, SEQ_TOP + 26); ctx.stroke();
      ctx.setLineDash([]);
      const g = Math.round(gOf(hover.x));
      if (g >= 0 && g < len) {
        const label = `${(g + 1).toLocaleString()}`;
        ctx.font = `10.5px ${FONT_DISPLAY}`;
        const tw = ctx.measureText(label).width + 8;
        const bx = Math.min(W - tw, Math.max(0, hover.x - tw / 2));
        ctx.fillStyle = "#0e131a";
        ctx.fillRect(bx, 0, tw, 14);
        ctx.strokeStyle = C.borderStrong; ctx.strokeRect(bx, 0, tw, 14);
        ctx.fillStyle = C.text;
        ctx.textAlign = "center";
        ctx.fillText(label, bx + tw / 2, 7.5);
      }
    }

    // ---- minimap
    const mini = miniRef.current;
    if (mini) {
      mini.width = W * dpr; mini.height = 42 * dpr;
      mini.style.width = `${W}px`; mini.style.height = "42px";
      const mctx = mini.getContext("2d");
      mctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      mctx.clearRect(0, 0, W, 42);
      mctx.fillStyle = "rgba(255,255,255,0.03)";
      mctx.fillRect(0, 0, W, 42);
      for (const f of features) {
        const x1 = (f.start / len) * W, x2 = (f.end / len) * W;
        mctx.fillStyle = TYPE_COLORS[f.type] ?? "#93a0ae";
        mctx.globalAlpha = 0.65;
        mctx.fillRect(x1, f.strand === "-" ? 24 : 10, Math.max(1, x2 - x1), 8);
        mctx.globalAlpha = 1;
      }
      const vx1 = (view.start / len) * W, vx2 = (view.end / len) * W;
      mctx.strokeStyle = "#ffffff";
      mctx.lineWidth = 1.4;
      mctx.strokeRect(vx1, 1.5, Math.max(3, vx2 - vx1), 39);
      mctx.lineWidth = 1;
    }
  }, [seq, len, features, tracks, view, size, selectedFeatureId, selection, hover.x, oriPos, motifHits, gOf]);

  const stats = useMemo(() => {
    if (!selection || selection.end <= selection.start0) return null;
    const gc = regionGC(seq, selection.start0, selection.end);
    const nf = features.filter((f) => f.start < selection.end && f.end > selection.start0).length;
    return { lenBp: selection.end - selection.start0, gc, nf };
  }, [selection, seq, features]);

  return (
    <div>
      <div style={{ fontSize: 11.5, color: C.textFaint, marginBottom: 4, display: "flex", justifyContent: "space-between", flexWrap: "wrap", gap: 6 }}>
        <span>
          {contig.id} · {len.toLocaleString()} bp{contig.circular ? " · circular" : ""}
        </span>
        <span style={{ fontFamily: FONT_DISPLAY }}>
          viewing {(view.start + 1).toLocaleString()}–{view.end.toLocaleString()} ({Math.round(view.end - view.start).toLocaleString()} bp)
        </span>
      </div>
      <div ref={wrapRef} style={{ border: `1px solid ${C.border}`, borderRadius: 2, overflow: "hidden", background: "#05070a" }}>
        <canvas
          ref={canvasRef}
          style={{ display: "block", cursor: mode === "pan" ? "grab" : "crosshair" }}
          onWheel={handleWheel}
          onMouseDown={handleMouseDown}
          onMouseMove={handleMouseMove}
          onMouseUp={handleMouseUp}
          onMouseLeave={handleMouseLeave}
          onDoubleClick={handleDoubleClick}
        />
        <canvas
          ref={miniRef}
          style={{ display: "block", borderTop: `1px solid ${C.border}`, cursor: "pointer" }}
          onMouseDown={(e) => {
            const rect = e.currentTarget.getBoundingClientRect();
            const frac = (e.clientX - rect.left) / rect.width;
            const span = view.end - view.start;
            onViewChange(clampView(frac * len - span / 2, frac * len + span / 2));
          }}
        />
      </div>
      {stats && (
        <div style={{
          marginTop: 8, padding: "8px 12px", borderRadius: 2, background: "#0c1018",
          border: "1px solid rgba(232,193,90,0.35)", fontFamily: FONT_DISPLAY, fontSize: 12, color: C.text,
          display: "flex", gap: 18, flexWrap: "wrap",
        }}>
          <strong style={{ color: C.qc, textShadow: `0 0 6px ${C.qc}44` }}>Selected region</strong>
          <span>{stats.lenBp.toLocaleString()} bp</span>
          <span>GC {(stats.gc * 100).toFixed(1)}%</span>
          <span>{stats.nf} feature(s)</span>
          <span style={{ color: C.textDim }}>{(view.start + 1).toLocaleString()} … context</span>
        </div>
      )}
    </div>
  );
}

function roundBar(ctx, x, y, w, h, fwd) {
  const r = Math.min(5, h / 2, Math.abs(w) / 2);
  ctx.beginPath();
  if (fwd) {
    ctx.moveTo(x, y);
    ctx.lineTo(x + Math.max(0, w - r), y);
    ctx.lineTo(x + w, y + h / 2);
    ctx.lineTo(x + Math.max(0, w - r), y + h);
    ctx.lineTo(x, y + h);
  } else {
    ctx.moveTo(x + w, y);
    ctx.lineTo(x + r, y);
    ctx.lineTo(x, y + h / 2);
    ctx.lineTo(x + r, y + h);
    ctx.lineTo(x + w, y + h);
  }
  ctx.closePath();
  ctx.fill();
}

function roundBarStroke(ctx, x, y, w, h, fwd) {
  const r = Math.min(5, h / 2, Math.abs(w) / 2);
  ctx.beginPath();
  if (fwd) {
    ctx.moveTo(x, y); ctx.lineTo(x + w - r, y); ctx.lineTo(x + w, y + h / 2);
    ctx.lineTo(x + w - r, y + h); ctx.lineTo(x, y + h);
  } else {
    ctx.moveTo(x + w, y); ctx.lineTo(x + r, y); ctx.lineTo(x, y + h / 2);
    ctx.lineTo(x + r, y + h); ctx.lineTo(x + w, y + h);
  }
  ctx.closePath();
  ctx.stroke();
}

let cachedMeanGC = { key: "", v: 0.5 };
function meanGCOf(seq) {
  if (cachedMeanGC.key !== seq) {
    cachedMeanGC = { key: seq, v: regionGC(seq, 0, seq.length) };
  }
  return cachedMeanGC.v;
}
