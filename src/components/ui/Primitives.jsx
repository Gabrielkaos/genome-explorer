import { C, FONT } from "../../theme.js";

export function Panel({ children, style, className = "", title }) {
  return (
    <div
      className={`scanlines ${className}`}
      style={{
        background: C.bgPanel,
        border: `1px solid ${C.border}`,
        borderRadius: 2,
        position: "relative",
        overflow: "hidden",
        ...style,
      }}
    >
      {title && (
        <div style={{
          padding: "4px 10px",
          borderBottom: `1px solid ${C.border}`,
          fontSize: 11,
          color: C.prompt,
          fontFamily: FONT,
          display: "flex",
          alignItems: "center",
          gap: 6,
          background: "#0a0e16",
        }}>
          <span style={{ color: C.textFaint }}>$</span> {title}
          <span className="cursor-blink" />
        </div>
      )}
      {children}
    </div>
  );
}

export function Eyebrow({ color, children }) {
  return (
    <div style={{
      fontFamily: FONT, fontSize: 11, letterSpacing: "0.06em",
      textTransform: "uppercase", color,
      display: "flex", alignItems: "center", gap: 6, marginBottom: 6,
    }}>
      <span style={{ color: C.textFaint }}>›</span>
      {children}
    </div>
  );
}

export function StatCard({ label, value, color, unit }) {
  return (
    <div style={{
      background: C.bgPanel,
      border: `1px solid ${C.border}`,
      borderRadius: 2,
      padding: "8px 12px",
      minWidth: 100,
    }}>
      <div style={{
        fontSize: 10, color: C.textFaint,
        textTransform: "uppercase", letterSpacing: "0.08em",
        fontFamily: FONT,
      }}>{label}</div>
      <div style={{
        fontFamily: FONT, fontSize: 18, color: color || C.text,
        marginTop: 2, textShadow: color ? `0 0 8px ${color}44` : "none",
      }}>
        {value}<span style={{ fontSize: 10, color: C.textDim, marginLeft: 3 }}>{unit}</span>
      </div>
    </div>
  );
}

export function ExplainBox({ explainMode, color, children }) {
  if (!explainMode) return null;
  return (
    <div style={{
      marginTop: 10, padding: "8px 12px", borderRadius: 2,
      fontSize: 12, lineHeight: 1.55, fontFamily: FONT,
      background: `${color}08`, borderLeft: `2px solid ${color}`,
      color: C.text,
    }}>
      <span style={{ color, marginRight: 6 }}>#</span>
      {children}
    </div>
  );
}

export function LimitBanner({ children }) {
  return (
    <div style={{
      marginTop: 12, padding: "8px 12px", borderRadius: 2,
      fontSize: 11.5, lineHeight: 1.5, fontFamily: FONT,
      background: "rgba(230,104,95,0.06)",
      borderLeft: "2px solid rgba(230,104,95,0.6)",
      color: "#f2b3ad",
    }}>
      <span style={{ color: C.bad, marginRight: 6 }}>⚠</span>
      {children}
    </div>
  );
}

export function SectionTitle({ icon: Icon, color, title, subtitle }) {
  return (
    <div style={{ marginBottom: 16 }}>
      <div style={{
        fontFamily: FONT, fontSize: 13.5, color,
        display: "flex", alignItems: "center", gap: 8,
        textShadow: `0 0 10px ${color}33`,
      }}>
        <span style={{ color: C.prompt }}>genome@explorer:~$</span>
        <Icon size={14} color={color} />
        {title.toLowerCase().replace(/\s+/g, '_')}
      </div>
      {subtitle && (
        <div style={{
          fontSize: 11.5, color: C.textFaint, marginTop: 3,
          paddingLeft: 24, fontFamily: FONT,
          borderLeft: `1px solid ${C.border}`,
          marginLeft: 8,
        }}>
          # {subtitle}
        </div>
      )}
    </div>
  );
}

