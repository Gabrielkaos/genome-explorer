import { useState } from "react";
import { Check, Copy, Download, FileInput, X } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow } from "../ui/Primitives.jsx";
import { parseNewick } from "../../lib/phylo/newick.js";

/**
 * Newick exchange panel: shows the inferred tree's Newick string with
 * copy/download, and accepts a pasted Newick (e.g. an IQ-TREE .treefile)
 * for visualization in the viewer above.
 */
export default function NewickPanel({ newick, onImport, importActive, onClearImport }) {
  const [text, setText] = useState("");
  const [copied, setCopied] = useState(false);
  const [error, setError] = useState(null);

  async function handleCopy() {
    try {
      await navigator.clipboard.writeText(newick);
      setCopied(true);
      setTimeout(() => setCopied(false), 1400);
    } catch { /* clipboard unavailable */ }
  }

  function handleVisualize() {
    try {
      const tree = parseNewick(text);
      setError(null);
      onImport(tree);
    } catch (err) {
      setError(err.message);
      onImport(null);
    }
  }

  return (
    <Panel style={{ padding: 16 }}>
      <Eyebrow color={C.phylo}>Newick string &amp; import</Eyebrow>

      <textarea readOnly value={newick} rows={5} spellCheck={false} onFocus={(e) => e.target.select()}
        style={{
          width: "100%", marginTop: 8, background: "#05070a", border: `1px solid ${C.border}`,
          borderRadius: 2, padding: "10px 12px", color: C.textDim, fontFamily: FONT_DISPLAY,
          fontSize: 11.5, resize: "vertical",
        }} />
      <div style={{ display: "flex", gap: 8, marginTop: 8, flexWrap: "wrap" }}>
        <MiniBtn onClick={handleCopy}>{copied ? <Check size={12} color={C.good} /> : <Copy size={12} />}{copied ? " [ copied ]" : " [ copy ]"}</MiniBtn>
        <MiniBtn onClick={() => saveText(newick, "tree.nwk")}><Download size={12} /> [ .nwk ]</MiniBtn>
        <span style={{ fontSize: 10.5, color: C.textFaint, alignSelf: "center", fontFamily: FONT_DISPLAY }}>
          paste into iTOL / FigTree / MEGA, or keep as a record
        </span>
      </div>

      <div style={{ borderTop: `1px solid ${C.border}`, margin: "14px 0 12px" }} />
      <div style={{ fontSize: 12, color: C.text, marginBottom: 6, fontFamily: FONT_DISPLAY }}>
        <FileInput size={12} style={{ verticalAlign: "-2px" }} /> Visualize a tree you already have
      </div>
      <textarea value={text} onChange={(e) => setText(e.target.value)} rows={3} spellCheck={false}
        placeholder={"Paste any Newick here - e.g. IQ-TREE's .treefile or RAxML bestTree:\n(A:0.1,B:0.2,(C:0.15,D:0.05)97:0.1);"}
        style={{
          width: "100%", background: "#05070a", border: `1px solid ${error ? `${C.bad}88` : C.border}`,
          borderRadius: 2, padding: "10px 12px", color: C.text, fontFamily: FONT_DISPLAY,
          fontSize: 11.5, resize: "vertical",
        }} />
      <div style={{ display: "flex", gap: 8, marginTop: 8, flexWrap: "wrap", alignItems: "center" }}>
        <MiniBtn accent onClick={handleVisualize}><FileInput size={12} /> [ visualize_pasted_tree ]</MiniBtn>
        {importActive && (
          <>
            <span style={{ fontSize: 11, color: C.qc, fontFamily: FONT_DISPLAY }}>showing imported tree</span>
            <button onClick={() => { onClearImport(); setText(""); }}
              style={{ all: "unset", cursor: "pointer", fontSize: 11, color: C.textDim, display: "inline-flex", alignItems: "center", gap: 3, fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em" }}>
              <X size={11} /> [ clear ]
            </button>
          </>
        )}
      </div>
      {error && (
        <div style={{ fontSize: 11.5, color: "#f2b3ad", marginTop: 8, fontFamily: FONT_DISPLAY }}>{error}</div>
      )}
    </Panel>
  );
}

function MiniBtn({ children, onClick, accent }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 5,
      fontSize: 11.5, padding: "5px 10px", borderRadius: 2,
      border: `1px solid ${accent ? C.raw : C.border}66`,
      color: accent ? C.raw : C.textDim,
      fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
      background: "#05070a",
      textShadow: accent ? `0 0 6px ${C.raw}44` : "none",
    }}>{children}</button>
  );
}

function saveText(text, fileName) {
  const blob = new Blob([text], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = fileName;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}
