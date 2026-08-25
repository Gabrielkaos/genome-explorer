import { useCallback, useMemo, useRef, useState } from "react";
import { CircuitBoard, Upload, FileText, MousePointer2, Hand, ZoomIn, ZoomOut, Maximize } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel, Eyebrow, SectionTitle, LimitBanner, ExplainBox } from "../ui/Primitives.jsx";
import { parseGff3, parseGenbank } from "../../lib/explorer/gff.js";
import { buildTracks, skewExtrema } from "../../lib/explorer/tracks.js";
import { generateExampleGenome } from "../../lib/sampleData/generateExampleGenome.js";
import { buildFeaturesGff3, buildFeaturesTsv } from "../../lib/explorer/exportExplorer.js";
import { downloadBlob } from "../../lib/fastq/exportFastq.js";
import GenomeCanvas from "./GenomeCanvas.jsx";
import FeatureInspector from "./FeatureInspector.jsx";
import SearchPanel from "./SearchPanel.jsx";

const TYPE_ORDER = { CDS: 0, tRNA: 1, rRNA: 2, tmRNA: 3, ncRNA: 4 };

function saveText(text, fileName) {
  downloadBlob(new Blob([text], { type: "text/plain;charset=utf-8" }), fileName);
}

function parseFastaText(text) {
  const records = [];
  let cur = null;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (cur) records.push(cur);
      const header = line.slice(1).trim();
      cur = { id: header.split(/\s+/)[0] || `seq_${records.length + 1}`, desc: header, seq: "" };
    } else if (cur) cur.seq += line.trim().toUpperCase();
  }
  if (cur) records.push(cur);
  return records.filter((r) => r.seq.length > 0);
}

function synthesizeContigRecords(featuresByContig, sequenceRegions = new Map()) {
  const records = [];
  const contigIds = new Set([...sequenceRegions.keys(), ...featuresByContig.keys()]);
  for (const id of contigIds) {
    const reg = sequenceRegions.get(id);
    const feats = featuresByContig.get(id) ?? [];
    const maxEnd = feats.reduce((m, f) => Math.max(m, f.end), 0);
    const len = reg?.length ?? Math.max(maxEnd, 1000);
    records.push({
      id,
      desc: `${id} (from GFF3 annotations)`,
      seq: "N".repeat(len),
      circular: reg?.circular ?? false,
      placeholderSeq: true,
    });
  }
  return records.length ? records : [{ id: "sequence", desc: "sequence", seq: "N".repeat(1000), circular: false, placeholderSeq: true }];
}

export default function ExplorerSection({ explainMode }) {
  const [source, setSource] = useState("example"); // upload | example
  const [upload, setUpload] = useState(null);      // {records:[{id,desc,seq,circular}], featuresByContig, name, placeholderSeq}
  const [contigIdx, setContigIdx] = useState(0);
  const [mode, setMode] = useState("pan");         // pan | select
  const [parseError, setParseError] = useState("");
  // Per-contig UI state (viewport / selection / selection tag / motif hits),
  // keyed by contig+source so switching data resets it without an effect.
  const [uiByContig, setUiByContig] = useState({});
  const genomeRef = useRef(null), annotRef = useRef(null);

  /* --------------------------- source data --------------------------- */

  const example = useMemo(() => {
    const g = generateExampleGenome(11);
    return {
      records: [{ id: g.id, desc: g.desc, seq: g.seq, circular: true }],
      featuresByContig: new Map([[g.id, g.features]]),
      name: "Example replicon",
    };
  }, []);

  const active = source === "upload" ? (upload ?? example) : example;
  const record = active.records && active.records.length > 0
    ? active.records[Math.min(Math.max(0, contigIdx), active.records.length - 1)]
    : null;

  const features = useMemo(() => {
    if (!record || !active.featuresByContig) return [];
    return (active.featuresByContig.get(record.id) ?? []).slice().sort(
      (a, b) => a.start - b.start || (TYPE_ORDER[a.type] ?? 9) - (TYPE_ORDER[b.type] ?? 9)
    );
  }, [active, record]);

  const tracks = useMemo(
    () => (!record?.seq ? { window: 100, step: 50, starts: [], gc: [], skew: [] } : buildTracks(record.seq)),
    [record]
  );
  const oriPos = useMemo(
    () => (record?.circular ? skewExtrema(tracks)?.originPos : undefined),
    [record?.circular, tracks]
  );

  // Per-contig UI state, keyed by contig+source so switching data resets it
  // during render (no effect needed).
  const uiKey = record ? `${source}:${record.id}` : `${source}:none`;
  const curUi = uiByContig[uiKey] ?? { view: null, selection: null, tag: null, hits: [] };
  const defaultView = useMemo(
    () => (!record?.seq ? { start: 0, end: 0 } : { start: 0, end: Math.min(record.seq.length, 16000) }),
    [record]
  );
  const view = curUi.view ?? defaultView;
  // Stable updater so child components can hold it in effects safely.
  const patchUi = useCallback((partial) => {
    setUiByContig((m) => {
      const base = m[uiKey] ?? { view: null, selection: null, tag: null, hits: [] };
      return { ...m, [uiKey]: { ...base, ...partial } };
    });
  }, [uiKey]);
  const setView = useCallback((v) => patchUi({ view: v }), [patchUi]);
  const setSelection = useCallback((s) => patchUi({ selection: s }), [patchUi]);
  const setSelectedTag = useCallback((t) => patchUi({ tag: t }), [patchUi]);
  const setMotifHits = useCallback((h) => patchUi({ hits: h }), [patchUi]);
  const selectedTag = curUi.tag;
  const motifHits = curUi.hits;
  const selection = curUi.selection;

  const selectedFeature = useMemo(
    () => features.find((f) => f.locusTag === selectedTag) ?? null,
    [features, selectedTag]
  );

  /* ------------------------------ uploads ------------------------------ */

  function handleGenomeFile(file) {
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const text = String(reader.result);
        let parsed;
        if (/^\s*LOCUS\s/m.test(text)) {
          const gb = parseGenbank(text);
          parsed = {
            records: gb.records.map((r) => ({ ...r, circular: /circular/i.test(text.slice(0, 400)) })),
            featuresByContig: gb.featuresByContig,
            name: file.name,
            placeholderSeq: false,
          };
          if (!parsed.records.some((r) => r.seq.length)) throw new Error("GenBank file carries no ORIGIN sequence.");
        } else if (/^\s*##gff-version/m.test(text) || (text.includes("\t") && /\t(CDS|gene|rRNA|tRNA|exon|region)\t/i.test(text))) {
          // GFF3 / GTF uploaded via genome button
          const gff = parseGff3(text, "sequence");
          let records = gff.records;
          let placeholderSeq = false;
          if (!records.length) {
            records = synthesizeContigRecords(gff.featuresByContig, gff.sequenceRegions);
            placeholderSeq = true;
          }
          parsed = {
            records,
            featuresByContig: gff.featuresByContig,
            name: file.name,
            placeholderSeq,
          };
        } else {
          const recs = parseFastaText(text);
          if (!recs.length) throw new Error("No sequences found.");
          const baseFeatures = upload?.featuresByContig ?? new Map();
          parsed = {
            records: recs.map((r) => ({ id: r.id, desc: r.desc ?? r.id, seq: r.seq.toUpperCase().replace(/[^ACGTN]/gi, "N"), circular: false })),
            featuresByContig: baseFeatures,
            name: file.name,
            placeholderSeq: false,
          };
        }
        setUpload(parsed);
        setSource("upload");
        setContigIdx(0);
        setParseError("");
      } catch (err) {
        setParseError(`Genome file: ${err.message}`);
      }
    };
    reader.readAsText(file);
  }

  function handleAnnotFile(file) {
    const reader = new FileReader();
    reader.onload = () => {
      try {
        const text = String(reader.result);
        const existingRecords = upload?.records?.length
          ? upload.records
          : [];

        const fallbackId = existingRecords[0]?.id ?? "sequence";

        if (/^\s*LOCUS\s/m.test(text)) {
          const gb = parseGenbank(text);
          const records = gb.records.length ? gb.records : existingRecords;
          const merged = new Map(upload?.featuresByContig ?? new Map());
          for (const [k, v] of gb.featuresByContig) merged.set(k, v);
          setUpload({
            records: records.length ? records : synthesizeContigRecords(merged),
            featuresByContig: merged,
            name: upload?.name ? `${upload.name} + ${file.name}` : file.name,
            placeholderSeq: records.length === 0,
          });
        } else {
          const gff = parseGff3(text, fallbackId);
          if (!gff.numFeatures) throw new Error("No usable features found in this GFF / annotation file.");

          let records = [];
          let placeholderSeq = false;

          if (gff.records.length > 0) {
            records = gff.records;
          } else if (existingRecords.length > 0) {
            records = existingRecords.slice();
            // Synthesize any contigs from GFF that are not in existingRecords
            for (const seqid of gff.featuresByContig.keys()) {
              if (!records.some((r) => r.id === seqid)) {
                const len = gff.sequenceRegions.get(seqid)?.length ?? Math.max(...(gff.featuresByContig.get(seqid) ?? []).map((f) => f.end), 1000);
                records.push({
                  id: seqid,
                  desc: `${seqid} (from GFF3 annotations)`,
                  seq: "N".repeat(len),
                  circular: gff.sequenceRegions.get(seqid)?.circular ?? false,
                  placeholderSeq: true,
                });
              }
            }
          } else {
            records = synthesizeContigRecords(gff.featuresByContig, gff.sequenceRegions);
            placeholderSeq = true;
          }

          // Remap single contig if FASTA has 1 record and GFF has 1 contig with different name
          if (records.length === 1 && gff.featuresByContig.size === 1 && !gff.featuresByContig.has(records[0].id)) {
            const onlyGffSeqId = [...gff.featuresByContig.keys()][0];
            const remappedFeatures = (gff.featuresByContig.get(onlyGffSeqId) ?? []).map((f) => ({
              ...f,
              contigId: records[0].id,
            }));
            gff.featuresByContig.set(records[0].id, remappedFeatures);
          }

          setUpload({
            records,
            featuresByContig: gff.featuresByContig,
            name: upload?.name ? `${upload.name} + ${file.name}` : file.name,
            placeholderSeq,
          });
        }
        setSource("upload");
        setContigIdx(0);
        setParseError("");
      } catch (err) {
        setParseError(`Annotation file: ${err.message}`);
      }
    };
    reader.readAsText(file);
  }

  /* ------------------------------- view ------------------------------- */

  const span = view ? view.end - view.start : 16000;
  function zoom(factor) {
    if (!view) return;
    const mid = (view.start + view.end) / 2;
    setView({
      start: Math.max(0, mid - (span * factor) / 2),
      end: Math.min(record.seq.length, mid + (span * factor) / 2),
    });
  }

  return (
    <div>
      <SectionTitle icon={CircuitBoard} color={C.raw}
        title="Unified Genome Explorer"
        subtitle="Interactive browser over your contigs: features, GC content/skew, IUPAC motif search, sequence extraction - runs entirely in your browser" />

      <LimitBanner>
        A lightweight Artemis-style viewer, not a database browser: it shows whatever you load - your own
        FASTA plus GFF3, or a self-contained GenBank file. No remote genomes, no BLAST, no variant overlays;
        products are exactly what the annotation carried.
      </LimitBanner>

      {/* ------------------------- INPUT ------------------------- */}
      <Panel style={{ padding: 18, marginTop: 16, marginBottom: 16 }}>
        <Eyebrow color={C.raw}>Load a genome</Eyebrow>
        <div style={{ display: "flex", gap: 10, marginTop: 8, flexWrap: "wrap" }}>
          <SrcCard title="Uploaded files" color={C.raw}
            enabled={!!upload}
            sub={upload ? `${upload.name} — ${upload.records.length} record(s), ${[...upload.featuresByContig.values()].reduce((a, v) => a + v.length, 0)} features` : "FASTA / GenBank (.gbk,.gbff) + optional GFF3"}
            onClick={() => { setSource("upload"); setContigIdx(0); }} />
          <SrcCard title="Example replicon" color={C.qc}
            sub="64 kb circular plasmid · ~50 features · planted motifs"
            onClick={() => { setSource("example"); setContigIdx(0); }} />
        </div>

        <div style={{ display: "flex", gap: 8, marginTop: 12, flexWrap: "wrap", alignItems: "center" }}>
          <button onClick={() => genomeRef.current.click()} style={upBtn}><Upload size={12} /> [ upload_fasta_gbk ]</button>
          <input ref={genomeRef} type="file" accept=".fasta,.fa,.fna,.fas,.gbk,.gbff,.gb,.genbank,.ffn,.txt" hidden onChange={(e) => e.target.files[0] && handleGenomeFile(e.target.files[0])} />
          <button onClick={() => annotRef.current.click()} style={upBtn}><Upload size={12} /> [ upload_gff3 ]</button>
          <input ref={annotRef} type="file" accept=".gff,.gff3,.txt" hidden onChange={(e) => e.target.files[0] && handleAnnotFile(e.target.files[0])} />
          <span style={{ fontSize: 11, color: C.textFaint }}>GFF seqids must match FASTA headers; GenBank needs its ORIGIN.</span>
        </div>

        {parseError && (
          <div style={{ marginTop: 10, fontSize: 12, color: "#f2b3ad" }}>⚠ {parseError}</div>
        )}
      </Panel>

      {/* ------------------------- BROWSER ------------------------- */}
      {view && record && (
        <>
          {record.placeholderSeq && (
            <div style={{
              padding: "9px 13px", marginBottom: 12, borderRadius: 2,
              background: "#0c131d", border: `1px solid ${C.raw}66`,
              color: C.raw, fontSize: 12, display: "flex", alignItems: "center", gap: 10,
              fontFamily: FONT_DISPLAY,
            }}>
              <span style={{ fontSize: 14 }}>ℹ</span>
              <span>
                <strong>GFF3 annotations loaded without sequence.</strong> Feature map, coordinates, and search are active. Upload your genome FASTA to enable GC curves, skew, and base-level sequence.
              </span>
            </div>
          )}

          <div style={{ display: "flex", justifyContent: "space-between", gap: 10, flexWrap: "wrap", marginBottom: 12 }}>
            <div style={{ display: "flex", gap: 8, flexWrap: "wrap" }}>
              {active.records.map((r, i) => (
                <Chip key={r.id} active={i === contigIdx} onClick={() => setContigIdx(i)}>{r.id}</Chip>
              ))}
            </div>
            <div style={{ display: "flex", gap: 6, alignItems: "center", flexWrap: "wrap" }}>
              <Chip active={mode === "pan"} onClick={() => setMode("pan")}><Hand size={11} style={{ verticalAlign: "-2px" }} /> Pan</Chip>
              <Chip active={mode === "select"} onClick={() => setMode("select")}><MousePointer2 size={11} style={{ verticalAlign: "-2px" }} /> Select region</Chip>
              <span style={{ width: 1, height: 16, background: C.border }} />
              <IconBtn title="Zoom out" onClick={() => zoom(2)}><ZoomOut size={13} /></IconBtn>
              <IconBtn title="Zoom in" onClick={() => zoom(0.5)}><ZoomIn size={13} /></IconBtn>
              <IconBtn title="Whole contig" onClick={() => setView({ start: 0, end: record.seq.length })}><Maximize size={13} /></IconBtn>
            </div>
          </div>

          <Panel style={{ padding: 14, marginBottom: 16 }}>
            <GenomeCanvas
              contig={record}
              features={features}
              tracks={tracks}
              oriPos={oriPos}
              view={view}
              onViewChange={setView}
              selection={selection}
              onSelectionChange={setSelection}
              mode={mode}
              selectedFeatureId={selectedTag}
              onSelectFeature={(tag) => { setSelectedTag(tag); }}
              motifHits={motifHits}
            />
            <div style={{ display: "flex", gap: 14, flexWrap: "wrap", marginTop: 10, fontSize: 11, color: C.textFaint }}>
              <LegendDot color="#ef7fa3" label="CDS" />
              <LegendDot color="#e8a95a" label="rRNA" />
              <LegendDot color="#e8c15a" label="tRNA" />
              <LegendDot color="#a08cf0" label="GC skew" />
              <LegendDot color="#5ec8d8" label="GC%" />
              <LegendDot color="#68c98f" label="ori? (skew min)" />
              <span>scroll = zoom · drag = pan · double-click gene = focus</span>
            </div>
          </Panel>

          <div style={{ display: "grid", gridTemplateColumns: "1fr 340px", gap: 16, alignItems: "start" }}>
            <FeatureList features={features} selectedTag={selectedTag} onSelect={(tag) => {
              setSelectedTag(tag);
              const f = features.find((x) => x.locusTag === tag);
              if (f) {
                const pad = Math.max(60, (f.end - f.start));
                setView({ start: Math.max(0, f.start - pad), end: Math.min(record.seq.length, f.end + pad) });
              }
            }} />

            <div style={{ display: "grid", gap: 16 }}>
              <SearchPanel contig={record} features={features} onJump={setView} onMarkHits={setMotifHits} />
              <FeatureInspector contig={record} feature={selectedFeature} selection={selection} />
            </div>
          </div>

          <div style={{ display: "flex", gap: 8, flexWrap: "wrap", margin: "16px 0 4px" }}>
            <span style={{ fontSize: 11, color: C.textFaint, alignSelf: "center" }}>Exports:</span>
            <ExpBtn label="Features TSV" onClick={() => saveText(buildFeaturesTsv(features), `${record.id}_features.tsv`)} />
            <ExpBtn label="Features GFF3" onClick={() => saveText(buildFeaturesGff3(features, [[record.id, record.seq.length]]), `${record.id}_features.gff3`)} />
            <ExpBtn label="Sequence FASTA" onClick={() => saveText(
              `>${record.id}${record.desc ? ` ${record.desc}` : ""}\n${(record.seq.match(/.{1,70}/g) ?? []).join("\n")}\n`,
              `${record.id}.fasta`)} />
          </div>

          <ExplainBox explainMode={explainMode} color={C.raw}>
            One screen, every representation: coordinates on the ruler, genes as arrows on their strands, composition
            underneath as GC% and skew curves, raw bases once you zoom close enough. The skew minimum on a circular
            replicon usually marks the replication origin - click around it and see which genes live nearby.
          </ExplainBox>
        </>
      )}

      {!view && (
        <Panel style={{ padding: 16, fontSize: 12.5, color: C.textFaint }}>
          Load a genome above to start browsing.
        </Panel>
      )}
    </div>
  );
}

/* ------------------------------- pieces ------------------------------- */

const upBtn = {
  all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 5,
  padding: "6px 11px", borderRadius: 2, border: `1px solid ${C.border}`, color: C.textDim, fontSize: 12,
  fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
  background: "#05070a",
};

function SrcCard({ title, sub, enabled, onClick, color }) {
  return (
    <button onClick={() => enabled && onClick()} style={{
      all: "unset", cursor: enabled ? "pointer" : "default", flex: "1 1 220px", padding: "12px 14px", borderRadius: 2,
      background: "#05070a", border: `1px solid ${C.border}`,
      opacity: enabled ? 1 : 0.5,
    }}>
      <div style={{ fontSize: 12.5, color: enabled ? (color || C.text) : C.textDim, fontFamily: FONT_DISPLAY, textShadow: enabled && color ? `0 0 6px ${color}33` : "none" }}>{title}</div>
      <div style={{ fontSize: 11, color: C.textFaint, marginTop: 3, fontFamily: FONT_DISPLAY }}>{sub}</div>
    </button>
  );
}

function Chip({ children, active, onClick }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", fontSize: 11.5, padding: "4px 10px", borderRadius: 2,
      background: active ? `${C.raw}22` : "#05070a",
      border: `1px solid ${active ? C.raw : C.border}`,
      color: active ? C.raw : C.textDim,
      fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
      textShadow: active ? `0 0 6px ${C.raw}44` : "none",
    }}>{children}</button>
  );
}

function IconBtn({ children, onClick, title }) {
  return (
    <button onClick={onClick} title={title} style={{
      all: "unset", cursor: "pointer", display: "inline-flex", padding: 5, borderRadius: 2,
      border: `1px solid ${C.border}`, color: C.textDim, background: "#05070a",
    }}>{children}</button>
  );
}

function LegendDot({ color, label }) {
  return (
    <span style={{ display: "inline-flex", alignItems: "center", gap: 4 }}>
      <span style={{ width: 7, height: 7, borderRadius: 999, background: color, display: "inline-block" }} />{label}
    </span>
  );
}

function ExpBtn({ label, onClick }) {
  return (
    <button onClick={onClick} style={{
      all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 5,
      fontSize: 12, color: C.good, border: `1px solid ${C.good}66`, borderRadius: 2, padding: "6px 12px",
      fontFamily: FONT_DISPLAY, textTransform: "uppercase", letterSpacing: "0.06em",
      textShadow: `0 0 6px ${C.good}44`,
    }}>
      <FileText size={12} /> [ {label} ]
    </button>
  );
}

const TYPE_COLORS = {
  CDS: "#ef7fa3", tRNA: "#e8c15a", rRNA: "#e8a95a", tmRNA: "#e8c15a", ncRNA: "#c9b3f5",
};

function FeatureList({ features, selectedTag, onSelect }) {
  const [query, setQuery] = useState("");
  const shown = useMemo(() => {
    const q = query.trim().toLowerCase();
    const base = q ? features.filter((f) => f.locusTag.toLowerCase().includes(q) || f.product?.toLowerCase().includes(q)) : features;
    return base.slice(0, 300);
  }, [features, query]);

  return (
    <Panel style={{ padding: 14 }}>
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "baseline", marginBottom: 8, gap: 10, flexWrap: "wrap" }}>
        <Eyebrow color={C.annotation}>Features ({features.length.toLocaleString()})</Eyebrow>
        <input value={query} onChange={(e) => setQuery(e.target.value)} placeholder="filter…" spellCheck={false}
          style={{
            width: 150, background: "#05070a", border: `1px solid ${C.border}`, borderRadius: 2,
            padding: "4px 9px", color: C.text, fontFamily: FONT_DISPLAY, fontSize: 11.5,
          }} />
      </div>
      <div style={{ maxHeight: 430, overflowY: "auto" }}>
        <table style={{ borderCollapse: "collapse", width: "100%", fontFamily: FONT_DISPLAY, fontSize: 11.5 }}>
          <tbody>
            {shown.map((f, i) => (
              <tr key={`${f.locusTag}_${f.start}_${f.end}_${i}`} onClick={() => onSelect(f.locusTag)}
                style={{ cursor: "pointer", borderBottom: `1px solid ${C.border}`,
                  background: selectedTag === f.locusTag ? `${C.raw}18` : "transparent" }}>
                <td style={{ padding: "4px 8px", whiteSpace: "nowrap" }}>
                  <span style={{ display: "inline-block", width: 8, height: 8, borderRadius: 2, background: TYPE_COLORS[f.type] ?? "#93a0ae", marginRight: 6 }} />
                  <span style={{ color: selectedTag === f.locusTag ? C.raw : C.text, textShadow: selectedTag === f.locusTag ? `0 0 6px ${C.raw}44` : "none" }}>{f.locusTag}</span>
                </td>
                <td style={{ padding: "4px 8px", color: C.textDim, textAlign: "right", whiteSpace: "nowrap" }}>
                  {(f.start + 1).toLocaleString()}–{f.end.toLocaleString()}<span style={{ color: f.strand === "+" ? C.good : C.bad }}>{f.strand === "+" ? " →" : " ←"}</span>
                </td>
                <td style={{ padding: "4px 8px", color: C.textFaint, maxWidth: 220, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>{f.product}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </Panel>
  );
}
