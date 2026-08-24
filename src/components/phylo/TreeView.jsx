import { useEffect, useMemo, useRef, useState } from "react";
import { ZoomIn, ZoomOut, RotateCcw } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";

/**
 * SVG phylogenetic tree viewer.
 *   layout: "rect" | "radial"
 *   mode:   "phylogram" (branch lengths drawn to scale) | "cladogram"
 * Supports bootstrap dots + labels, hover inspection, click-to-highlight
 * (leaf -> lineage to root, internal dot -> whole clade), wheel zoom and
 * drag panning. Pure rendering: the tree comes precomputed from the worker.
 */

export default function TreeView({ tree, layout = "rect", mode = "phylogram", showSupports = true }) {
  const svgRef = useRef(null);
  const drag = useRef(null);

  // Interaction state is scoped to a "session" (tree+layout+mode combo):
  // switching either resets zoom, pan and selection via render-phase reset.
  const session = useMemo(() => ({ tree, layout, mode }), [tree, layout, mode]);
  const [st, setSt] = useState(() => ({
    session,
    view: { k: 1, x: 0, y: 0 },
    selLeaf: null,
    selClade: null,
    hover: null,
    dragging: false,
  }));
  if (st.session !== session) {
    setSt({
      session,
      view: { k: 1, x: 0, y: 0 },
      selLeaf: null,
      selClade: null,
      hover: null,
      dragging: false,
    });
  }
  const { view, selLeaf, selClade, hover, dragging } = st;
  const setView = (updater) => setSt((s) => ({ ...s, view: typeof updater === "function" ? updater(s.view) : updater }));
  const setHover = (v) => setSt((s) => ({ ...s, hover: v }));
  const setSelLeaf = (v) => setSt((s) => ({ ...s, selLeaf: typeof v === "function" ? v(s.selLeaf) : v }));
  const setSelClade = (v) => setSt((s) => ({ ...s, selClade: typeof v === "function" ? v(s.selClade) : v }));

  const geo = useMemo(() => buildGeometry(tree, layout, mode), [tree, layout, mode]);

  // Non-passive wheel handler for cursor-centered zoom.
  useEffect(() => {
    const el = svgRef.current;
    if (!el) return;
    const onWheel = (e) => {
      e.preventDefault();
      const rect = el.getBoundingClientRect();
      const cx = ((e.clientX - rect.left) / rect.width) * geo.W;
      const cy = ((e.clientY - rect.top) / rect.height) * geo.H;
      setView((v) => {
        const factor = e.deltaY < 0 ? 1.15 : 1 / 1.15;
        const k = Math.max(0.4, Math.min(12, v.k * factor));
        const s = k / v.k;
        return { k, x: cx - (cx - v.x) * s, y: cy - (cy - v.y) * s };
      });
    };
    el.addEventListener("wheel", onWheel, { passive: false });
    return () => el.removeEventListener("wheel", onWheel);
  }, [geo.W, geo.H]);

  if (!geo.ok) {
    return (
      <div style={{ padding: 16, fontSize: 12.5, color: C.textFaint }}>
        No tree to display yet.
      </div>
    );
  }

  const { W, H, padL, scale, edges, nodes, leaves, totalDepth } = geo;

  function edgeVisual(P, Ch) {
    const inLineage = selLeaf != null && idxContains(Ch._lo, Ch._hi, selLeaf);
    if (inLineage) return { stroke: C.pheno, w: 2.4, opacity: 1 };
    if (selClade) {
      if (Ch === selClade || rangeWithin(selClade, Ch)) return { stroke: C.pheno, w: 2.2, opacity: 1 };
      if (rangeWithin(Ch, selClade)) return { stroke: C.phylo, w: 1.3, opacity: 1 };
      return { stroke: C.phylo, w: 1.3, opacity: 0.22 };
    }
    const isHover = hover?.child === Ch || hover?.parent === P;
    return { stroke: C.phylo, w: isHover ? 2.4 : 1.4, opacity: 1 };
  }

  function renderEdges() {
    const els = [];
    for (const e of edges) {
      const { P, Ch } = e;
      const sty = edgeVisual(P, Ch);
      els.push(
        <g key={`e${e.id}`}>
          <path d={e.d} fill="none" stroke={sty.stroke} strokeWidth={sty.w} opacity={sty.opacity} strokeLinecap="round" />
          <path d={e.d} fill="none" stroke="transparent" strokeWidth={11}
            onMouseEnter={() => setHover(e)} onMouseLeave={() => setHover(null)} style={{ cursor: "pointer" }} />
        </g>
      );
    }
    return els;
  }

  const infoText = (() => {
    const e = hover;
    if (!e) {
      if (selClade) {
        const names = descNames(selClade);
        return `Selected clade (${names.length} taxon${names.length === 1 ? "" : "a"}): ${names.slice(0, 10).join(", ")}${names.length > 10 ? ` +${names.length - 10} more` : ""} — click the dot again or empty space to clear.`;
      }
      if (selLeaf != null) return `Highlighted lineage of ${leaves[selLeaf].name} — click again to clear.`;
      return "Hover a branch to inspect it · click a tip to trace its lineage · click an internal dot to light up that whole clade.";
    }
    const sup = Number.isFinite(e.child.support) ? ` · bootstrap ${e.child.support}%` : "";
    if (!e.parent) {
      const names = e.child._memberNames ?? [];
      return `Internal node${sup} · ${names.length} descendants: ${names.slice(0, 8).join(", ")}${names.length > 8 ? ` +${names.length - 8} more` : ""}`;
    }
    const desc = e.child.isLeaf ? "" : ` · clade of ${e.child._members.length}: ${(e.child._memberNames ?? []).slice(0, 6).join(", ")}${e.child._members.length > 6 ? ` +${e.child._members.length - 6} more` : ""}`;
    return `Branch → ${e.child.name ?? "internal node"} · length ${(e.child.branchLen ?? 0).toFixed(4)} substitutions/site${sup}${desc}`;
  })();

  return (
    <div>
      <div style={{ display: "flex", justifyContent: "flex-end", gap: 6, marginBottom: 6 }}>
        <IconBtn onClick={() => setView((v) => ({ ...v, k: Math.min(12, v.k * 1.25) }))}><ZoomIn size={14} /></IconBtn>
        <IconBtn onClick={() => setView((v) => ({ ...v, k: Math.max(0.4, v.k / 1.25) }))}><ZoomOut size={14} /></IconBtn>
        <IconBtn onClick={() => { setView({ k: 1, x: 0, y: 0 }); }}><RotateCcw size={14} /></IconBtn>
      </div>
      <svg ref={svgRef} viewBox={`0 0 ${W} ${H}`} style={{ width: "100%", maxHeight: H + 40, display: "block",
        borderRadius: 8, border: `1px solid ${C.border}`, background: "#0b0f15",
        cursor: dragging ? "grabbing" : "grab" }}
        onPointerDown={(e) => {
          drag.current = { sx: e.clientX, sy: e.clientY, vx: view.x, vy: view.y };
          setSt((s) => ({ ...s, dragging: true }));
          e.currentTarget.setPointerCapture(e.pointerId);
        }}
        onPointerMove={(e) => {
          if (!drag.current) return;
          const rect = svgRef.current.getBoundingClientRect();
          const dx = ((e.clientX - drag.current.sx) / rect.width) * W;
          const dy = ((e.clientY - drag.current.sy) / rect.height) * H;
          setView((v) => ({ ...v, x: drag.current.vx + dx, y: drag.current.vy + dy }));
        }}
        onPointerUp={() => { drag.current = null; setSt((s) => ({ ...s, dragging: false })); }}
        onClick={(e) => {
          if (e.target.tagName === "svg") { setSelLeaf(null); setSelClade(null); }
        }}
      >
        <g transform={`translate(${view.x},${view.y}) scale(${view.k})`}>
          {renderEdges()}

          {/* internal support dots */}
          {nodes.map((nd) => nd !== tree && Number.isFinite(nd.support) && (
            <g key={`n${nd._id}`}>
              <circle cx={nd._px} cy={nd._py} r={4.2} fill={supportColor(nd.support)} stroke="#0b0f15" strokeWidth={1}
                style={{ cursor: "pointer" }}
                onClick={() => { setSelClade((s) => (s === nd ? null : nd)); setSelLeaf(null); }}
                onMouseEnter={() => setHover({ parent: null, child: nd })}
                onMouseLeave={() => setHover(null)} />
              {showSupports && view.k > 0.7 && (
                <text x={nd._px - 6} y={nd._py - 7} fontSize={9.5} fontFamily={FONT_DISPLAY}
                  fill={supportColor(nd.support)} textAnchor="middle">{Math.round(nd.support)}</text>
              )}
            </g>
          ))}

          {/* root marker */}
          {!tree.unrooted && (
            <rect x={tree._px - 3.5} y={tree._py - 3.5} width={7} height={7}
              fill="none" stroke={C.textDim} strokeWidth={1.2} />
          )}

          {/* leaves */}
          {leaves.map((lf) => {
            const dimmed = selClade && !(lf._lo >= selClade._lo && lf._hi <= selClade._hi);
            const active = selLeaf === lf._idx;
            return (
              <g key={`l${lf._id}`} opacity={dimmed ? 0.45 : 1}>
                {layout === "rect" ? (
                  <>
                    <circle cx={lf._px} cy={lf._py} r={active ? 5 : 3.6} fill={active ? C.pheno : C.phylo} stroke={active ? "#fff" : "none"} strokeWidth={1.2}
                      style={{ cursor: "pointer" }} onClick={() => setSelLeaf(active ? null : lf._idx)} />
                    <text x={lf._px + 8} y={lf._py + 4} fontSize={12} fontFamily={FONT_DISPLAY}
                      fill={active ? C.pheno : C.text} style={{ cursor: "pointer", userSelect: "none" }}
                      onClick={() => setSelLeaf(active ? null : lf._idx)}>{trimName(lf.name)}</text>
                  </>
                ) : (
                  <>
                    <circle cx={lf._px} cy={lf._py} r={active ? 5 : 3.6} fill={active ? C.pheno : C.phylo} stroke={active ? "#fff" : "none"} strokeWidth={1.2}
                      style={{ cursor: "pointer" }} onClick={() => setSelLeaf(active ? null : lf._idx)} />
                    <text x={lf._lx} y={lf._ly} fontSize={11.5} fontFamily={FONT_DISPLAY}
                      textAnchor={lf._anchor} dominantBaseline="middle"
                      fill={active ? C.pheno : C.text} style={{ cursor: "pointer", userSelect: "none" }}
                      onClick={() => setSelLeaf(active ? null : lf._idx)}>{trimName(lf.name)}</text>
                  </>
                )}
              </g>
            );
          })}
        </g>

        {/* scale bar (fixed, outside zoom) */}
        {mode === "phylogram" && geo.scaleStep != null && (layout === "rect" ? (
          <g fontSize={10} fontFamily={FONT_DISPLAY} fill={C.textDim}>
            <line x1={padL} y1={H - 18} x2={padL + geo.scaleStep * scale} y2={H - 18} stroke={C.textDim} strokeWidth={1.5} />
            <line x1={padL} y1={H - 21} x2={padL} y2={H - 15} stroke={C.textDim} strokeWidth={1.5} />
            <line x1={padL + geo.scaleStep * scale} y1={H - 21} x2={padL + geo.scaleStep * scale} y2={H - 15} stroke={C.textDim} strokeWidth={1.5} />
            <text x={padL + (geo.scaleStep * scale) / 2} y={H - 24} textAnchor="middle">{geo.scaleStep} subs/site</text>
          </g>
        ) : (
          <text x={14} y={H - 12} fontSize={10} fontFamily={FONT_DISPLAY} fill={C.textDim}>
            ring spacing ≈ {geo.scaleStep} substitutions/site
          </text>
        ))}
      </svg>

      <div style={{ marginTop: 8, fontSize: 11.5, color: C.textDim, minHeight: 32, lineHeight: 1.5 }}>{infoText}</div>
      <div style={{ fontSize: 11, color: C.textFaint }}>
        Drag to pan · scroll to zoom · {showSupports ? "dots are colored by bootstrap support (<70 red, 70-89 amber, ≥90 green)" : "support labels hidden"}
        {totalDepth > 0 && mode === "phylogram" ? ` · deepest tip ≈ ${totalDepth.toFixed(3)} substitutions/site` : ""}
      </div>
    </div>
  );

  function descNames(node) {
    return node._memberNames ?? [];
  }
}

/* ------------------------------ geometry ------------------------------ */

function buildGeometry(tree, layout, mode) {
  if (!tree || typeof tree !== "object") return { ok: false };
  const leaves = [];
  let counter = 0;
  let nodeId = 0;

  function prep(node) {
    node._id = nodeId++;
    node._len = mode === "cladogram" ? 1 : Math.max(0, node.branchLen || 0);
    if (!node.children || !node.children.length) node.isLeaf = true;
    if (node.isLeaf) {
      node._idx = counter;
      node._lo = node._hi = counter;
      node._members = [counter];
      node._memberNames = [node.name];
      leaves.push(node);
      counter++;
      return;
    }
    for (const c of node.children) prep(c);
    node._lo = Math.min(...node.children.map((c) => c._lo));
    node._hi = Math.max(...node.children.map((c) => c._hi));
    node._members = [].concat(...node.children.map((c) => c._members));
    node._memberNames = [].concat(...node.children.map((c) => c._memberNames)).sort();
  }
  prep(tree);
  if (leaves.length < 2) return { ok: false };

  const nLeaves = leaves.length;

  if (layout === "rect") {
    const rowH = Math.max(26, Math.min(34, 700 / nLeaves));
    const H = Math.max(300, Math.min(780, 70 + nLeaves * rowH));
    const W = 980;
    const padL = 46, padT = 26;
    const labelRoom = 150;

    let totalRaw = 0;
    (function depths(node, acc) {
      node._raw = acc;
      if (node.isLeaf) { totalRaw = Math.max(totalRaw, acc); return; }
      for (const c of node.children) depths(c, acc + c._len);
    })(tree, 0);
    if (totalRaw <= 0) return { ok: false };

    const scale = (W - padL - labelRoom) / totalRaw;
    (function placeY(node) {
      if (node.isLeaf) { node._row = node._idx; return; }
      node.children.forEach(placeY);
      const rows = node.children.map((c) => c._row);
      node._row = (Math.min(...rows) + Math.max(...rows)) / 2;
    })(tree);

    const edges = [];
    const nodesList = [];
    (function collect(node) {
      nodesList.push(node);
      node._px = padL + node._raw * scale;
      node._py = padT + node._row * rowH;
      if (node.isLeaf) return;
      for (const c of node.children) {
        c._px = padL + c._raw * scale;
        c._py = padT + c._row * rowH;
        edges.push({
          id: `${node._id}-${c._id}`, parent: node, child: c,
          d: `M ${node._px} ${node._py} V ${c._py} H ${c._px}`,
        });
        collect(c);
      }
    })(tree);

    const steps = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];
    const scaleStep = steps.find((s) => s * scale >= 64) ?? steps[steps.length - 1];

    return { ok: true, W, H, padL, padT, scale, edges, nodes: nodesList, leaves, totalDepth: totalRaw, scaleStep };
  }

  // ---- radial ----
  const size = Math.min(860, Math.max(520, nLeaves * 44));
  const W = size, H = size;
  const cxp = W / 2, cyp = H / 2;
  const Rmax = size / 2 - 110;
  const k = nLeaves;

  let totalRaw = 0;
  (function depths(node, acc) {
    node._raw = acc;
    if (node.isLeaf) { totalRaw = Math.max(totalRaw, acc); return; }
    for (const c of node.children) depths(c, acc + c._len);
  })(tree, 0);
  if (totalRaw <= 0) return { ok: false };
  const scale = Rmax / totalRaw;

  (function placeAngle(node) {
    if (node.isLeaf) { node._ang = -Math.PI / 2 + (2 * Math.PI * (node._idx + 0.5)) / k; return; }
    node.children.forEach(placeAngle);
    let sin = 0, cos = 0;
    for (const c of node.children) {
      const w = c._members.length;
      sin += Math.sin(c._ang) * w;
      cos += Math.cos(c._ang) * w;
    }
    node._ang = Math.atan2(sin, cos);
  })(tree);

  const edges = [];
  const nodesList = [];
  (function collect(node) {
    nodesList.push(node);
    node._pr = node._raw * scale;
    node._px = cxp + node._pr * Math.cos(node._ang);
    node._py = cyp + node._pr * Math.sin(node._ang);
    if (node.isLeaf) {
      const lr = node._pr + 12;
      node._lx = cxp + lr * Math.cos(node._ang);
      node._ly = cyp + lr * Math.sin(node._ang);
      const c = Math.cos(node._ang);
      node._anchor = c >= 0.08 ? "start" : c <= -0.08 ? "end" : "middle";
      return;
    }
    for (const ch of node.children) {
      ch._pr = ch._raw * scale;
      ch._px = cxp + ch._pr * Math.cos(ch._ang);
      ch._py = cyp + ch._pr * Math.sin(ch._ang);
      const x1 = cxp + node._pr * Math.cos(node._ang);
      const y1 = cyp + node._pr * Math.sin(node._ang);
      edges.push({ id: `${node._id}-${ch._id}`, parent: node, child: ch, d: `M ${x1} ${y1} L ${ch._px} ${ch._py}` });
      collect(ch);
    }
  })(tree);

  const steps = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1];
  const scaleStep = steps.find((s) => s * scale >= 56) ?? steps[steps.length - 1];

  return { ok: true, W, H, padL: 0, padT: 0, scale, edges, nodes: nodesList, leaves, totalDepth: totalRaw, scaleStep };
}

function idxContains(lo, hi, idx) { return idx >= lo && idx <= hi; }
function rangeWithin(outer, inner) { return inner._lo >= outer._lo && inner._hi <= outer._hi; }

function supportColor(s) {
  if (!Number.isFinite(s)) return C.textFaint;
  if (s >= 90) return C.good;
  if (s >= 70) return C.qc;
  return C.bad;
}

function trimName(name) {
  const s = String(name ?? "");
  return s.length > 26 ? `${s.slice(0, 24)}…` : s;
}

function IconBtn({ children, onClick }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", padding: 6, borderRadius: 6,
      border: `1px solid ${C.border}`, color: C.textDim, display: "flex",
    }}>{children}</button>
  );
}
