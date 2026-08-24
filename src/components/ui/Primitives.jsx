import React from "react";
import { Sparkles, AlertTriangle } from "lucide-react";
import { C, FONT_DISPLAY, FONT_HEAD } from "../../theme.js";

export function Panel({ children, style, className = "" }) {
  return (
    <div
      className={className}
      style={{
        background: `linear-gradient(180deg, ${C.bgPanel2}, ${C.bgPanel})`,
        border: `1px solid ${C.border}`, borderRadius: 12,
        backdropFilter: "blur(14px)",
        boxShadow: "0 1px 0 rgba(255,255,255,0.04) inset, 0 12px 28px -16px rgba(0,0,0,0.6)",
        ...style,
      }}
    >{children}</div>
  );
}

export function Eyebrow({ color, children }) {
  return (
    <div style={{
      fontFamily: FONT_DISPLAY, fontSize: 11, letterSpacing: "0.14em", textTransform: "uppercase",
      color, display: "flex", alignItems: "center", gap: 8, marginBottom: 6,
    }}>
      <span style={{ width: 6, height: 6, borderRadius: 999, background: color, boxShadow: `0 0 8px ${color}` }} />
      {children}
    </div>
  );
}

export function StatCard({ label, value, color, unit }) {
  return (
    <div style={{ background: C.bgPanel2, border: `1px solid ${C.border}`, borderRadius: 8, padding: "10px 14px", minWidth: 108 }}>
      <div style={{ fontSize: 10.5, color: C.textFaint, textTransform: "uppercase", letterSpacing: "0.08em" }}>{label}</div>
      <div style={{ fontFamily: FONT_DISPLAY, fontSize: 20, color: color || C.text, marginTop: 3 }}>
        {value}<span style={{ fontSize: 11, color: C.textDim, marginLeft: 3 }}>{unit}</span>
      </div>
    </div>
  );
}

export function ExplainBox({ explainMode, color, children }) {
  if (!explainMode) return null;
  return (
    <div style={{
      marginTop: 10, padding: "10px 12px", borderRadius: 8, fontSize: 13, lineHeight: 1.5,
      background: `${color}14`, border: `1px solid ${color}44`, color: C.text,
      display: "flex", gap: 8, alignItems: "flex-start",
    }}>
      <Sparkles size={14} color={color} style={{ marginTop: 2, flexShrink: 0 }} />
      <span>{children}</span>
    </div>
  );
}

export function LimitBanner({ children }) {
  return (
    <div style={{
      marginTop: 14, padding: "10px 12px", borderRadius: 8, fontSize: 12.5, lineHeight: 1.5,
      background: "rgba(230,104,95,0.08)", border: "1px solid rgba(230,104,95,0.3)", color: "#f2b3ad",
      display: "flex", gap: 8, alignItems: "flex-start",
    }}>
      <AlertTriangle size={14} style={{ marginTop: 2, flexShrink: 0, color: C.bad }} />
      <span>{children}</span>
    </div>
  );
}

export function SectionTitle({ icon: Icon, color, title, subtitle }) {
  return (
    <div style={{ display: "flex", alignItems: "flex-start", gap: 12, marginBottom: 18 }}>
      <div style={{
        width: 38, height: 38, borderRadius: 9, background: `${color}1c`, border: `1px solid ${color}55`,
        display: "flex", alignItems: "center", justifyContent: "center", flexShrink: 0,
      }}>
        <Icon size={18} color={color} />
      </div>
      <div>
        <div style={{ fontFamily: FONT_HEAD, fontWeight: 600, fontSize: 20, color: C.text, letterSpacing: "-0.01em" }}>{title}</div>
        <div style={{ fontSize: 13, color: C.textDim, marginTop: 2 }}>{subtitle}</div>
      </div>
    </div>
  );
}

export const tooltipStyle = { background: "#0e131a", border: `1px solid ${C.border}`, borderRadius: 8, fontSize: 12, color: C.text };
