import { FileText } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";

export function SourceCard({ title, sub, active, enabled, onBrowse, color }) {
  return (
    <button onClick={() => { if (enabled) onBrowse(); }} style={{
      all: "unset", cursor: enabled ? "pointer" : "default", flex: "1 1 200px", padding: "12px 14px", borderRadius: 8,
      background: active ? `${color}14` : "#05070a",
      border: `1px solid ${active ? color : C.border}`,
      opacity: enabled ? 1 : 0.45,
    }}>
      <div style={{ fontSize: 12.5, color: active ? color : C.textDim }}>{title}</div>
      <div style={{ fontSize: 11, color: C.textFaint, marginTop: 3 }}>{sub}</div>
    </button>
  );
}

/** Labeled radio-chip group with per-option hint text under the label. */
export function Segmented({ label, value, onChange, options }) {
  const activeHint = options.find((o) => o.v === value)?.hint;
  return (
    <div>
      <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 5, alignItems: "baseline" }}>
        <span style={{ fontSize: 12.5, color: C.text }}>{label}</span>
        <span style={{ fontSize: 10, color: C.textFaint }}>{activeHint}</span>
      </div>
      <div style={{ display: "flex", gap: 4, flexWrap: "wrap" }}>
        {options.map((o) => (
          <button key={o.v} onClick={() => onChange(o.v)} title={o.hint} style={{
            all: "unset", cursor: "pointer", fontSize: 11, padding: "4px 9px", borderRadius: 6,
            background: value === o.v ? `${C.phylo}22` : "transparent",
            border: `1px solid ${value === o.v ? C.phylo : C.border}`,
            color: value === o.v ? C.phylo : C.textDim,
          }}>{o.label}</button>
        ))}
      </div>
    </div>
  );
}

export function Slider({ label, value, min, max, step, fmt, onChange, color }) {
  return (
    <div>
      <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 6 }}>
        <span style={{ fontSize: 12.5, color: C.text }}>{label}</span>
        <span style={{ fontFamily: FONT_DISPLAY, color, fontSize: 13 }}>{fmt(value)}</span>
      </div>
      <input type="range" min={min} max={max} step={step} value={value}
        onChange={(e) => onChange(+e.target.value)} style={{ width: "100%" }} />
    </div>
  );
}

export function NumberField({ label, value, onChange, hint }) {
  return (
    <div>
      <div style={{ display: "flex", justifyContent: "space-between", marginBottom: 6 }}>
        <span style={{ fontSize: 12.5, color: C.text }}>{label}</span>
        <span style={{ fontSize: 10, color: C.textFaint }}>{hint}</span>
      </div>
      <input type="number" min={0} max={999999} value={value}
        onChange={(e) => onChange(Math.max(0, Math.floor(+e.target.value || 0)))}
        style={{
          width: 110, background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 6,
          padding: "5px 8px", color: C.text, fontFamily: FONT_DISPLAY, fontSize: 13,
        }} />
    </div>
  );
}

export function ParamNote({ children }) {
  return (
    <div style={{
      alignSelf: "start", fontSize: 11.5, lineHeight: 1.5, color: C.textFaint,
      background: "rgba(255,255,255,0.03)", border: `1px solid ${C.border}`,
      borderRadius: 8, padding: "8px 10px",
    }}>{children}</div>
  );
}

export function DownloadBtn({ label, onClick }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 6,
      fontSize: 12, color: C.good, border: `1px solid ${C.good}66`, borderRadius: 6, padding: "6px 12px",
    }}>
      <FileText size={12} /> {label}
    </button>
  );
}
