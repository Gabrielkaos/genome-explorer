import { forwardRef, useEffect, useImperativeHandle, useMemo, useRef, useState } from "react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { BASE_COLORS } from "./colors.js";

const LABEL_W = 148;
const RULER_H = 16;
const CONS_H = 28;
const CONSENSUS_H = 18;
const CANVAS_H = 430;
const ZOOMS = [1.5, 2.5, 4, 6, 8, 11, 15];

/**
 * Canvas MSA viewer: colored bases (or identity shading), a per-column
 * conservation track, ruler and consensus row. Everything is virtualized -
 * only the visible column/row window is drawn - so alignments of many
 * thousands of columns stay smooth. Wheel scrolls vertically, shift+wheel
 * or trackpad-X pans horizontally, drag pans both axes, click selects the
 * nearest column.
 */
const AlignmentViewer = forwardRef(function AlignmentViewer({
  ids, rows, length, consensusCodes,
  selectedCol, onSelectCol, scrollTarget, colorMode = "nucleotide",
}, ref) {
  const wrapRef = useRef(null);
  const canvasRef = useRef(null);
  const [width, setWidth] = useState(900);
  const [zoomIdx, setZoomIdx] = useState(() => {
    const fit = (900 - LABEL_W) / Math.max(length, 1);
    let idx = 0;
    for (let i = ZOOMS.length - 1; i >= 0; i--) if (ZOOMS[i] <= fit) { idx = i; break; }
    return idx;
  });
  const view = useRef({ x: 0, y: 0 });
  const drag = useRef(null);
  const [, forceTick] = useState(0);

  useImperativeHandle(ref, () => ({
    zoomBy(dir) {
      setZoomIdx((z) => Math.max(0, Math.min(ZOOMS.length - 1, z + dir)));
    },
    reset() {
      setZoomIdx(0);
      view.current.x = 0;
      view.current.y = 0;
      forceTick((t) => t + 1);
    },
  }), []);

  const pxBase = ZOOMS[zoomIdx];
  const rowH = Math.max(9, Math.min(24, pxBase * 1.7));
  const fontSize = Math.max(7, Math.min(rowH - 2, pxBase - 0.5));

  // Per-column majority frequency (conservation track + identity mode).
  const colMajority = useMemo(() => {
    const out = new Float32Array(length);
    const n = rows.length;
    const counts = new Int32Array(5);
    for (let c = 0; c < length; c++) {
      counts.fill(0);
      for (let r = 0; r < n; r++) {
        const ch = rows[r][c];
        counts[ch === "-" ? 4 : "ACGT".indexOf(ch)]++;
      }
      let max = 0;
      for (let a = 0; a < 4; a++) if (counts[a] > max) max = counts[a];
      out[c] = n ? max / n : 0;
    }
    return out;
  }, [rows, length]);

  useEffect(() => {
    const el = wrapRef.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => setWidth(Math.max(320, entries[0].contentRect.width)));
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  useEffect(() => { view.current.x = 0; view.current.y = 0; }, [rows]);

  useEffect(() => {
    if (!scrollTarget) return;
    const span = (width - LABEL_W) / pxBase;
    view.current.x = Math.max(0, Math.min(length - span, scrollTarget.col - span / 2));
    forceTick((t) => t + 1);
  }, [scrollTarget, width, pxBase, length]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const maxX = Math.max(0, length - (width - LABEL_W) / pxBase);
    const onWheel = (e) => {
      e.preventDefault();
      if (e.shiftKey || Math.abs(e.deltaX) > Math.abs(e.deltaY)) {
        view.current.x = Math.max(0, Math.min(maxX, view.current.x + (e.deltaX !== 0 ? e.deltaX : e.deltaY) / pxBase));
      } else {
        view.current.y += e.deltaY;
      }
      forceTick((t) => t + 1);
    };
    canvas.addEventListener("wheel", onWheel, { passive: false });
    return () => canvas.removeEventListener("wheel", onWheel);
  }, [pxBase, width, length]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const dpr = window.devicePixelRatio || 1;
    const W = width, H = CANVAS_H;
    canvas.width = W * dpr; canvas.height = H * dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    ctx.fillStyle = "#0b0f15";
    ctx.fillRect(0, 0, W, H);

    const gridW = W - LABEL_W;
    const c0 = Math.max(0, Math.floor(view.current.x));
    const c1 = Math.min(length, Math.ceil(view.current.x + gridW / pxBase) + 1);
    const bodyTop = RULER_H + CONS_H;
    const availH = H - bodyTop - CONSENSUS_H;
    const r0 = Math.max(0, Math.floor(view.current.y));
    const r1 = Math.min(rows.length, r0 + Math.ceil(availH / rowH));
    const xOf = (c) => LABEL_W + (c - view.current.x) * pxBase;

    // conservation track
    for (let c = c0; c < c1; c++) {
      const frac = colMajority[c];
      if (frac <= 0) continue;
      const h = frac * (CONS_H - 6);
      ctx.fillStyle = `rgba(104,201,143,${0.25 + frac * 0.75})`;
      ctx.fillRect(xOf(c), CONS_H - h, Math.max(1, pxBase - 0.5), h);
    }

    // ruler
    ctx.font = `9px ${FONT_DISPLAY}`;
    ctx.fillStyle = C.textFaint;
    let step = 2;
    for (const s of [2, 5, 10, 20, 50, 100, 200, 500, 1000, 2500]) { step = s; if (s * pxBase >= 46) break; }
    for (let c = Math.ceil(c0 / step) * step; c < c1; c += step) {
      const x = xOf(c);
      ctx.fillRect(x, CONS_H + 3, 1, 4);
      ctx.fillText(String(c + 1), x + 2, CONS_H + 13);
    }

    // cells
    const textMode = pxBase >= 5;
    ctx.textAlign = "center";
    for (let r = r0; r < r1; r++) {
      const y = bodyTop + (r - r0) * rowH;
      const row = rows[r];
      for (let c = c0; c < c1; c++) {
        const ch = row[c];
        if (!ch) continue;
        const isGap = ch === "-";
        const code = isGap ? 4 : "ACGT".indexOf(ch);
        const x = xOf(c);
        if (colorMode === "identity") {
          const matchesConsensus = consensusCodes ? consensusCodes[c] === code : false;
          const conserved = colMajority[c] >= 0.999;
          if (isGap) ctx.fillStyle = "#141a22";
          else if (matchesConsensus && conserved) ctx.fillStyle = "rgba(104,201,143,0.55)";
          else if (matchesConsensus) ctx.fillStyle = "rgba(104,201,143,0.18)";
          else ctx.fillStyle = "rgba(230,104,95,0.45)";
          ctx.fillRect(x, y, Math.max(1, pxBase - 0.5), rowH - 1);
          if (textMode && !isGap && !conserved) {
            ctx.fillStyle = "#ffd9d2";
            ctx.font = `${fontSize}px ${FONT_DISPLAY}`;
            ctx.fillText(ch, x + pxBase / 2, y + rowH - 4);
          }
        } else if (isGap) {
          ctx.fillStyle = "#161c25";
          ctx.fillRect(x, y, Math.max(1, pxBase - 0.5), rowH - 1);
        } else if (textMode) {
          ctx.fillStyle = BASE_COLORS[code];
          ctx.font = `${fontSize}px ${FONT_DISPLAY}`;
          ctx.fillText(ch, x + pxBase / 2, y + rowH - 4);
        } else {
          ctx.fillStyle = BASE_COLORS[code];
          ctx.fillRect(x, y, Math.max(1, pxBase - 0.5), rowH - 1);
        }
      }
    }

    // selection overlay
    if (selectedCol != null && selectedCol >= c0 && selectedCol < c1) {
      const x = xOf(selectedCol);
      ctx.fillStyle = "rgba(255,255,255,0.10)";
      ctx.fillRect(x, CONS_H, pxBase, availH + CONSENSUS_H);
      ctx.strokeStyle = "rgba(255,255,255,0.7)";
      ctx.lineWidth = 1;
      ctx.strokeRect(x - 0.5, CONS_H + 0.5, pxBase + 1, availH);
    }

    // labels
    ctx.textAlign = "left";
    ctx.font = `10px ${FONT_DISPLAY}`;
    ctx.fillStyle = "rgba(11,15,21,0.85)";
    ctx.fillRect(0, 0, LABEL_W - 8, H);
    for (let r = r0; r < r1; r++) {
      const y = bodyTop + (r - r0) * rowH;
      ctx.fillStyle = r % 2 ? C.textDim : C.text;
      ctx.fillText(ids[r].slice(0, Math.floor((LABEL_W - 14) / 6)), 6, y + rowH - 4);
      ctx.fillStyle = "rgba(255,255,255,0.05)";
      ctx.fillRect(LABEL_W - 6, y, 2, rowH - 1);
    }
    ctx.fillStyle = "rgba(255,255,255,0.08)";
    ctx.fillRect(LABEL_W - 4, 0, 1, H);
    ctx.fillStyle = "rgba(255,255,255,0.06)";
    ctx.fillRect(LABEL_W, CONS_H, gridW, 1);
    ctx.fillStyle = "rgba(255,255,255,0.06)";
    ctx.fillRect(LABEL_W, bodyTop + availH, gridW, 1);

    // consensus
    if (consensusCodes && consensusCodes.length === length) {
      ctx.font = `${Math.min(fontSize, 12)}px ${FONT_DISPLAY}`;
      ctx.textAlign = "center";
      for (let c = c0; c < c1; c++) {
        const code = consensusCodes[c];
        if (code > 3) continue;
        ctx.fillStyle = BASE_COLORS[code];
        ctx.fillText("ACGT"[code], xOf(c) + pxBase / 2, H - 5);
      }
    }
  });

  const clampX = (x) => Math.max(0, Math.min(Math.max(0, length - (width - LABEL_W) / pxBase), x));

  function onPointerDown(e) {
    drag.current = { sx: e.clientX, sy: e.clientY, ox: view.current.x, oy: view.current.y, moved: false };
    e.currentTarget.setPointerCapture(e.pointerId);
  }
  function onPointerMove(e) {
    if (!drag.current) return;
    const dx = e.clientX - drag.current.sx, dy = e.clientY - drag.current.sy;
    if (Math.abs(dx) + Math.abs(dy) > 4) drag.current.moved = true;
    if (drag.current.moved) {
      view.current.x = clampX(drag.current.ox - dx / pxBase);
      view.current.y = Math.max(0, drag.current.oy - dy);
      forceTick((t) => t + 1);
    }
  }
  function onPointerUp(e) {
    const wasClick = drag.current && !drag.current.moved;
    drag.current = null;
    if (!wasClick || !onSelectCol) return;
    const rect = e.currentTarget.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const col = Math.floor(view.current.x + (mx - LABEL_W) / pxBase);
    if (col >= 0 && col < length) onSelectCol(col);
  }

  return (
    <div ref={wrapRef}>
      <canvas
        ref={canvasRef}
        style={{ width: "100%", height: CANVAS_H, borderRadius: 2, border: `1px solid ${C.border}`, cursor: drag.current?.moved ? "grabbing" : "crosshair", display: "block" }}
        onPointerDown={onPointerDown}
        onPointerMove={onPointerMove}
        onPointerUp={onPointerUp}
      />
      <div style={{ fontSize: 10.5, color: C.textFaint, marginTop: 6 }}>
        Click a column to inspect it · drag to pan · wheel scrolls strains · shift+wheel scrolls bases
      </div>
    </div>
  );
});

AlignmentViewer.displayName = "AlignmentViewer";
export default AlignmentViewer;
