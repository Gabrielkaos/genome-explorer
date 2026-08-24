import { useMemo } from "react";
import { Copy, Download } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow } from "../ui/Primitives.jsx";
import { translate } from "../../lib/annotation/protein.js";
import { revcomp } from "../../lib/explorer/search.js";
import { buildRegionFasta } from "../../lib/explorer/exportExplorer.js";
import { downloadBlob } from "../../lib/fastq/exportFastq.js";

/**
 * Detail view for the currently selected feature or region: coordinates,
 * computed GC, extracted sequence(s), translation for CDS features, and
 * copy/download actions.
 */
export default function FeatureInspector({ contig, feature, selection }) {
  const seq = contig.seq;

  const data = useMemo(() => {
    if (!feature) return null;
    const sub = seq.slice(feature.start, feature.end);
    const featSeq = feature.strand === "-" ? revcomp(sub) : sub;
    const gc = (() => {
      let n = 0;
      for (const b of sub) if (b === "G" || b === "C") n++;
      return sub.length ? n / sub.length : 0;
    })();
    return {
      ...feature,
      lengthNt: sub.length,
      gc,
      dna: featSeq,
      protein: feature.type === "CDS" && featSeq.length >= 3 ? translate(featSeq) : null,
    };
  }, [seq, feature]);

  const selData = useMemo(() => {
    if (selection && selection.end - selection.start0 > 0) {
      return { start0: selection.start0, end: selection.end };
    }
    return null;
  }, [selection]);

  function copy(text) {
    navigator.clipboard?.writeText(text).catch(() => {});
  }

  function downloadRegion() {
    const r = selData ?? (data ? { start0: data.start, end: data.end } : null);
    if (!r) return;
    const strand = data && !selData ? data.strand : "+";
    const name = selData
      ? `${contig.id}_${r.start0 + 1}_${r.end}_region`
      : `${data.locusTag.replace(/\W+/g, "_")}`;
    const fasta = buildRegionFasta(seq, `${name} ${contig.id}:${r.start0 + 1}-${r.end}${strand === "-" ? " (revcomp)" : ""}`, r.start0, r.end, strand);
    downloadBlob(new Blob([fasta], { type: "text/plain;charset=utf-8" }), `${name}.fasta`);
  }

  return (
    <Panel style={{ padding: 16 }}>
      <Eyebrow color={C.annotation}>{selData ? "Selected region" : "Feature inspector"}</Eyebrow>

      {!data && !selData && (
        <div style={{ fontSize: 12.5, color: C.textFaint }}>
          Click a gene on the map to inspect it, or drag with the Select tool to define a region.
        </div>
      )}

      {selData && !data && (
        <>
          <Row k="Contig" v={contig.id} />
          <Row k="Coordinates" v={`${(selData.start0 + 1).toLocaleString()} – ${selData.end.toLocaleString()} (${(selData.end - selData.start0).toLocaleString()} bp)`} />
          <Actions onCopy={() => copy(buildRegionFasta(seq, `${contig.id}:${selData.start0 + 1}-${selData.end}`, selData.start0, selData.end))} onDownload={downloadRegion} />
        </>
      )}

      {data && (
        <>
          <div style={{ fontFamily: FONT_DISPLAY, fontSize: 15, color: C.text, marginBottom: 2 }}>{data.locusTag}</div>
          <div style={{ fontSize: 12.5, color: C.textDim, marginBottom: 10 }}>{data.product || "(no product)"}</div>
          <Row k="Type" v={data.type} />
          <Row k="Contig" v={`${data.contigId}`} />
          <Row k="Coordinates" v={`${(data.start + 1).toLocaleString()} – ${data.end.toLocaleString()} (${data.strand === "-" ? "reverse" : "forward"} strand)`} />
          <Row k="Length" v={`${data.lengthNt.toLocaleString()} nt${data.protein ? ` · ${data.protein.length.toLocaleString()} aa` : ""}`} />
          <Row k="GC content" v={`${(data.gc * 100).toFixed(1)}%`} />

          {data.protein && (
            <>
              <SubHead>Translation</SubHead>
              <MonoBlock text={wrap(data.protein)} />
            </>
          )}
          <SubHead>Nucleotide ({data.strand === "-" ? "reverse complement" : "forward"})</SubHead>
          <MonoBlock text={wrap(data.dna.slice(0, 900)) + (data.dna.length > 900 ? "\n…" : "")} />

          <Actions
            onCopy={() => copy(data.protein ? `>${data.locusTag}\n${wrap(data.protein)}` : `>${data.locusTag}\n${wrap(data.dna)}`)}
            onDownload={downloadRegion}
          />
        </>
      )}
    </Panel>
  );
}

const wrap = (s) => (s.match(/.{1,60}/g)?.join("\n") ?? "");

function Row({ k, v }) {
  return (
    <div style={{ display: "flex", justifyContent: "space-between", gap: 12, padding: "4px 0", borderBottom: `1px solid ${C.border}`, fontSize: 12.5 }}>
      <span style={{ color: C.textFaint }}>{k}</span>
      <span style={{ color: C.text, fontFamily: FONT_DISPLAY, textAlign: "right" }}>{v}</span>
    </div>
  );
}

function SubHead({ children }) {
  return (
    <div style={{ fontSize: 10.5, textTransform: "uppercase", letterSpacing: "0.08em", color: C.textFaint, margin: "12px 0 5px" }}>{children}</div>
  );
}

function MonoBlock({ text }) {
  return (
    <div style={{
      background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 6, padding: "8px 10px",
      fontFamily: FONT_DISPLAY, fontSize: 11, color: C.raw, whiteSpace: "pre-wrap", wordBreak: "break-all",
      maxHeight: 130, overflowY: "auto", lineHeight: 1.5,
    }}>{text}</div>
  );
}

function Actions({ onCopy, onDownload }) {
  return (
    <div style={{ display: "flex", gap: 8, marginTop: 12, flexWrap: "wrap" }}>
      <MiniBtn onClick={onCopy}><Copy size={11} /> Copy FASTA</MiniBtn>
      <MiniBtn onClick={onDownload}><Download size={11} /> Download FASTA</MiniBtn>
    </div>
  );
}

function MiniBtn({ children, onClick }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 5,
      fontSize: 11.5, color: C.textDim, border: `1px solid ${C.border}`, borderRadius: 6, padding: "5px 10px",
    }}>{children}</button>
  );
}
