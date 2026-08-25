import React, { useMemo, useState } from "react";
import { ChevronLeft, ChevronRight, ZoomIn, ZoomOut, RotateCcw } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow } from "../ui/Primitives.jsx";

const MIN_WINDOW = 800;

function niceStep(windowLen, targetTicks = 9) {
  const rough = windowLen / targetTicks;
  const mag = Math.pow(10, Math.floor(Math.log10(rough)));
  for (const m of [1, 2, 5, 10]) {
    if (m * mag >= rough) return m * mag;
  }
  return 10 * mag;
}

/**
 * Linear genome map for one contig: forward-strand CDS arrows above the axis,
 * reverse-strand below, with a zoomable viewport, an overview strip, position
 * ruler, and a windowed GC-content track underneath.
 */
export default function GenomeMap({ contig, genes, gcPoints, selectedGeneId, onSelectGene }) {
  const [windowLen, setWindowLen] = useState(contig.seq.length);
  const [viewStart, setViewStart] = useState(0);

  const clampedWindow = Math.min(windowLen, contig.seq.length);
  const viewEnd = Math.min(viewStart + clampedWindow, contig.seq.length);

  const visible = useMemo(() => genes.filter((g) => g.end > viewStart && g.start < viewEnd), [genes, viewStart, viewEnd]);

  // greedy lane packing so overlapping arrows don't collide when zoomed out
  const lanes = useMemo(() => {
    const fwd = [], rev = [];
    const sorted = [...visible].sort((a, b) => a.start - b.start);
    const pack = (list, arr) => {
      for (const g of list) {
        let placed = false;
        for (const lane of arr) {
          if (g.start >= lane.lastEnd) { lane.items.push(g); lane.lastEnd = g.end; placed = true; break; }
        }
        if (!placed) arr.push({ items: [g], lastEnd: g.end });
      }
    };
    pack(sorted.filter((g) => g.strand === "+"), fwd);
    pack(sorted.filter((g) => g.strand === "-"), rev);
    return { fwd, rev };
  }, [visible]);

  const W = 900, H = 150 + lanes.fwd.length * 16 + lanes.rev.length * 16;
  const padL = 12, plotW = W - padL * 2;
  const midY = 24 + lanes.fwd.length * 16;
  const scaleX = plotW / clampedWindow;
  const x = (pos) => padL + (pos - viewStart) * scaleX;

  const step = niceStep(clampedWindow);
  const ticks = [];
  for (let t = Math.ceil(viewStart / step) * step; t <= viewEnd; t += step) ticks.push(t);

  const gcInView = (gcPoints ?? []).filter((p) => p.pos >= viewStart - 2000 && p.pos <= viewEnd + 2000);
  const gcBaseY = H - 26, gcAmp = 20;
  const gcPath = gcInView.map((p, i) => `${i === 0 ? "M" : "L"} ${x(p.pos).toFixed(1)} ${(gcBaseY - ((p.gc - 35) / 30) * gcAmp).toFixed(1)}`).join(" ");

  function zoom(f) {
    const c = viewStart + clampedWindow / 2;
    const w = Math.max(MIN_WINDOW, Math.min(contig.seq.length, Math.round(clampedWindow * f)));
    setWindowLen(w);
    setViewStart(Math.max(0, Math.min(c - w / 2, contig.seq.length - w)));
  }
  function shift(dir) {
    const d = Math.round(clampedWindow * 0.4) * dir;
    setViewStart(Math.max(0, Math.min(d + viewStart, contig.seq.length - clampedWindow)));
  }

  const labelGenes = clampedWindow < 24000;

  function renderArrow(g, laneIdx, strandDir) {
    const s = Math.max(g.start, viewStart - 500), e = Math.min(g.end, viewEnd + 500);
    const x0 = x(s), x1 = x(e);
    const h = 11, tip = Math.min(7, Math.max(3, (x1 - x0) * 0.25));
    const y = strandDir === 1 ? midY - 14 - laneIdx * 16 : midY + 3 + laneIdx * 16;
    const active = g.id === selectedGeneId || g.locusTag === selectedGeneId;
    const d = strandDir === 1
      ? `M ${x0} ${y} H ${x1 - tip} L ${x1} ${y + h / 2} L ${x1 - tip} ${y + h} H ${x0} Z`
      : `M ${x1} ${y} H ${x0 + tip} L ${x0} ${y + h / 2} L ${x0 + tip} ${y + h} H ${x1} Z`;
    return (
      <g key={g.locusTag} onClick={() => onSelectGene(g.locusTag)} style={{ cursor: "pointer" }}>
        <path d={d} fill={active ? "#ffffff" : "#ef7fa3"} fillOpacity={active ? 1 : 0.78}
          stroke={active ? C.annotation : "none"} strokeWidth={active ? 2 : 0} />
        {labelGenes && x1 - x0 > 64 && (
          <text x={(x0 + x1) / 2} y={y + h - 2.5} textAnchor="middle" fontSize={7.5} fontFamily={FONT_DISPLAY}
            fill={active ? "#000" : "#fff"} style={{ pointerEvents: "none" }}>{g.locusTag}</text>
        )}
      </g>
    );
  }

  return (
    <Panel style={{ padding: 16, marginBottom: 16 }}>
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 8 }}>
        <Eyebrow color={C.annotation}>{contig.id} · {contig.seq.length.toLocaleString()} bp · {genes.length} CDS</Eyebrow>
        <div style={{ display: "flex", gap: 6 }}>
          <IconBtn onClick={() => shift(-1)} title="Pan left"><ChevronLeft size={13} /></IconBtn>
          <IconBtn onClick={() => shift(1)} title="Pan right"><ChevronRight size={13} /></IconBtn>
          <IconBtn onClick={() => zoom(0.5)} title="Zoom in"><ZoomIn size={13} /></IconBtn>
          <IconBtn onClick={() => zoom(2)} title="Zoom out"><ZoomOut size={13} /></IconBtn>
          <IconBtn onClick={() => { setWindowLen(contig.seq.length); setViewStart(0); }} title="Reset view"><RotateCcw size={13} /></IconBtn>
        </div>
      </div>

      {/* overview strip */}
      <div style={{ height: 14, background: "#05070a", borderRadius: 2, border: `1px solid ${C.border}`, position: "relative", marginBottom: 6, cursor: "pointer" }}
        onClick={(e) => {
          const r = e.currentTarget.getBoundingClientRect();
          const frac = (e.clientX - r.left) / r.width;
          const w = clampedWindow;
          setViewStart(Math.max(0, Math.min(frac * contig.seq.length - w / 2, contig.seq.length - w)));
        }}>
        <div style={{
          position: "absolute", top: 0, bottom: 0,
          left: `${(viewStart / contig.seq.length) * 100}%`,
          width: `${(clampedWindow / contig.seq.length) * 100}%`,
          background: `${C.annotation}33`, border: `1px solid ${C.annotation}`,
          borderRadius: 1,
        }} />
      </div>

      <svg viewBox={`0 0 ${W} ${H}`} style={{ width: "100%" }}>
        <line x1={padL} y1={midY - 5.5} x2={padL + plotW} y2={midY - 5.5} stroke={C.borderStrong} strokeWidth={1.5} />
        {ticks.map((t) => (
          <g key={t}>
            <line x1={x(t)} y1={midY - 9} x2={x(t)} y2={midY - 2} stroke={C.textFaint} strokeWidth={1} />
            <text x={x(t)} y={midY - 13} fontSize={8} fill={C.textFaint} textAnchor="middle" fontFamily={FONT_DISPLAY}>
              {t.toLocaleString()}
            </text>
          </g>
        ))}
        {lanes.fwd.map((lane, i) => lane.items.map((g) => renderArrow(g, i, 1)))}
        {lanes.rev.map((lane, i) => lane.items.map((g) => renderArrow(g, i, -1)))}

        {/* GC track */}
        <line x1={padL} y1={gcBaseY} x2={padL + plotW} y2={gcBaseY} stroke={C.border} strokeWidth={1} strokeDasharray="3 3" />
        <path d={gcPath} fill="none" stroke={C.qc} strokeWidth={1.2} opacity={0.85} />
        <text x={padL} y={H - 6} fontSize={8.5} fill={C.textFaint} fontFamily={FONT_DISPLAY}>GC% (windowed, dashed = 35%)</text>
        <text x={padL + plotW} y={H - 6} fontSize={8.5} fill={C.textFaint} textAnchor="end" fontFamily={FONT_DISPLAY}>
          viewing {viewStart.toLocaleString()}–{viewEnd.toLocaleString()} bp
        </text>
      </svg>
      <div style={{ fontSize: 11, color: C.textFaint, marginTop: 4 }}>
        Arrows point in transcription direction (top rows = forward strand, bottom = reverse). Click a feature to inspect it; click the strip above to jump.
      </div>
    </Panel>
  );
}

function IconBtn({ children, onClick, title }) {
  return (
    <button onClick={onClick} title={title}
      style={{
        all: "unset", cursor: "pointer", padding: 5, borderRadius: 2, border: `1px solid ${C.border}`,
        color: C.textDim, display: "flex", background: "#05070a",
      }}>
      {children}
    </button>
  );
}
