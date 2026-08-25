import React, { useState, Suspense, lazy } from "react";
import {
  Microscope, GitBranch, Network,
  Filter, ChevronRight, X,
  Layers, Table2, Waypoints, CircuitBoard,
  ScrollText,
} from "lucide-react";
import { C, FONT, FONT_IMPORT } from "./theme.js";
import { Panel, Eyebrow, ExplainBox, LimitBanner, SectionTitle } from "./components/ui/Primitives.jsx";
import { FastqDataProvider } from "./state/FastqDataContext.jsx";
import { ErrorBoundary } from "./components/ui/ErrorBoundary.jsx";

const FastqSection = lazy(() => import("./components/fastq/FastqSection.jsx"));
const QCDashboardSection = lazy(() => import("./components/fastq/QCDashboardSection.jsx"));
const AssemblySection = lazy(() => import("./components/assembly/AssemblySection.jsx"));
const AnnotationSection = lazy(() => import("./components/annotation/AnnotationSection.jsx"));
const AlignmentSection = lazy(() => import("./components/msa/MsaSection.jsx"));
const PhyloSection = lazy(() => import("./components/phylo/PhyloSection.jsx"));
const AssocSection = lazy(() => import("./components/assoc/AssocSection.jsx"));
const ExplorerSection = lazy(() => import("./components/explorer/ExplorerSection.jsx"));

/* =====================================================================================
   DESIGN TOKENS
   Moved to ./theme.js and ./components/ui/Primitives.jsx as the app becomes modular.
   Every section is now a real tool operating on user data - FASTQ (./components/fastq/),
   QC Dashboard, Assembly (./components/assembly/), Annotation
   (./components/annotation/), Multiple Sequence Alignment (./components/msa/), Phylogenetic Inference
   (./components/phylo/), Pan-GWAS Association (./components/assoc/) and the
   Unified Genome Explorer (./components/explorer/) - each with its own web
   worker where heavy compute is involved.
 ===================================================================================== */

/* =====================================================================================
   STAGE METADATA  (drives the pipeline overview + detail panels)
 ===================================================================================== */
const STAGES = [
  {
    id: "fastq", navId: "fastq", label: "Nanopore FASTQ", short: "FASTQ", color: C.raw, icon: ScrollText,
    tagline: "Raw long reads straight off the sequencer",
    input: "Electrical signal (\"squiggle\") already basecalled into nucleotide sequences with per-base quality scores.",
    does: "This is a storage format, not an algorithm: four lines per read - an identifier, the called sequence, a '+' separator, and a Phred-encoded quality string, one character per base.",
    output: "Thousands of reads of uneven length (here, ~200-3000 bp) and uneven accuracy, typical of nanopore's raw per-base error rate.",
    question: "What did the sequencer actually observe, before any cleanup?",
    limits: "FASTQ says nothing about correctness on its own - a read can be high-confidence and still wrong, or low-confidence and still useful in aggregate.",
  },
  {
    id: "qc", navId: "qc", label: "Quality Control", short: "QC", color: C.qc, icon: Microscope,
    tagline: "NanoPlot-style statistics and visualization - not filtering",
    input: "The full raw FASTQ read set.",
    does: "Computes summary statistics and renders diagnostic plots (length distribution, quality distribution, yield) so you can judge whether a run is usable. NanoPlot itself never removes or edits a single base - it only measures and visualizes.",
    output: "Read-length and quality histograms, a length-vs-quality scatterplot, cumulative yield, and headline numbers like N50 and mean quality.",
    question: "Was this sequencing run good enough to proceed with?",
    limits: "QC describes the dataset; it doesn't decide what to do about a bad one. A separate filtering step (next) makes that call.",
  },
  {
    id: "filter", label: "Filtering", short: "Filter", color: C.qc, icon: Filter,
    tagline: "Discarding reads below a chosen quality/length threshold",
    input: "Raw reads plus a chosen mean-quality (and optionally length) cutoff.",
    does: "A separate, explicit step (conceptually distinct from NanoPlot) that removes whole reads failing the threshold. It does not correct individual low-quality bases within a retained read.",
    output: "A smaller, higher-confidence read set that assembly will actually use.",
    question: "Which reads are trustworthy enough to build a genome from?",
    limits: "Aggressive filtering trades coverage for accuracy; too strict a cutoff can starve the assembler of coverage in some regions.",
  },
  {
    id: "flye", navId: "assembly", label: "Flye Assembly", short: "Flye", color: C.assembly, icon: Waypoints,
    tagline: "Overlap-based long-read genome assembly",
    input: "Filtered long reads.",
    does: "Finds overlaps between reads that must originate from the same genomic region, builds an assembly graph from those overlaps, and resolves it into contiguous consensus sequences.",
    output: "A small number of long contigs in a `assembly.fasta` file, ideally approaching the true chromosome/plasmid structure.",
    question: "What is the underlying genome sequence that produced these reads?",
    limits: "Repeats longer than read overlaps, uneven coverage, or a noisy read set can fragment the assembly into more contigs than true replicons, or introduce local errors.",
  },
  {
    id: "bakta", navId: "annotation", label: "Genome Annotation", short: "Bakta", color: C.annotation, icon: Layers,
    tagline: "De novo gene calling and protein-level analysis",
    input: "Assembled contig sequences (no reads involved anymore).",
    does: "Predicts protein-coding genes de novo: trains a hexamer coding-statistics model on long ORFs, scans all six frames for stop-to-stop open reading frames, chooses each gene's true start codon using ribosome-binding-site (Shine-Dalgarno) evidence plus coding potential, then computes translated-protein properties (size, pI, hydropathy, transmembrane helices, signal peptides).",
    output: "A structured feature table with coordinates, strand, start-codon evidence and computed protein properties — plus standards-compliant GFF3/GenBank/FASTA exports you can feed into server-side annotators or NCBI.",
    question: "Where are the protein-coding genes, and what can be computed about them without external databases?",
    limits: "Without similarity search against reference databases, every product stays honestly \"hypothetical\" — this tool cannot name a gene's function, and RNA features (rRNA/tRNA) need covariance-model tools like Barrnap/tRNAscan-SE. Gene prediction identifies coding capacity, never confirmed function or phenotype.",
  },
  {
    id: "msa", navId: "msa", label: "Multiple Sequence Alignment", short: "Alignment", color: C.phylo, icon: Table2,
    tagline: "Lining up homologous sequences to compare them",
    input: "Two or more homologous nucleotide sequences - the same gene (or locus) from several strains, or whole small genomes like virus/mitochondria/plasmid. Upload a multi-FASTA, paste raw sequence text, or align this session's assembled contigs.",
    does: "Runs a genuine progressive multiple sequence alignment in your browser, the architecture ClustalW made standard: all-against-all k-mer distance estimates build an UPGMA guide tree, then closest pairs are aligned first and merged into growing profiles via global dynamic programming with affine gap penalties and sum-of-pairs scoring.",
    output: "The aligned matrix itself (equal-length gapped sequences) plus column-level analysis: conserved / variable / parsimony-informative site classification, a consensus sequence, the pairwise identity matrix and a variant-site table - exportable as aligned FASTA, Clustal, NEXUS, PHYLIP, variant TSV or consensus FASTA.",
    question: "Which positions differ between these strains, and where?",
    limits: "Guide-tree order comes from approximate k-mer distances, and the result is single-pass progressive alignment with no iterative refinement (MAFFT/ClustalOmega polish rounds) or substitution-model selection. Quality degrades with high divergence or many large indels; N/ambiguity codes are treated as gaps; a poor alignment silently corrupts everything built on top of it.",
  },
  {
    id: "iqtree", navId: "tree", label: "Phylogenetic Inference", short: "Phylo", color: C.phylo, icon: GitBranch,
    tagline: "Distance-based tree inference with bootstrap support",
    input: "An ALIGNED sequence matrix - the Alignment section's output, an uploaded aligned FASTA, or pasted homologs (3+ sequences).",
    does: "Estimates every pairwise distance under a substitution model (raw p-distance, Jukes-Cantor, or Kimura 2-parameter, with pairwise/complete gap deletion), then builds the tree whose path lengths best fit that matrix via Neighbor-Joining (clock-free) or UPGMA (molecular-clock, rooted). Confidence per clade comes from genuine nonparametric bootstraps: alignment columns are resampled with replacement and the whole pipeline is rebuilt for every replicate.",
    output: "A phylogenetic tree with branch lengths in expected substitutions/site and a bootstrap support percentage on each internal split - rendered as rectangular or radial phylogram/cladogram, and exportable as Newick/NEXUS plus a TSV distance matrix. Existing trees (IQ-TREE .treefile etc.) can be pasted in as Newick for visualization.",
    question: "How are these strains related by descent - and which parts of that hypothesis are actually well supported?",
    limits: "Distance methods are fast but not maximum likelihood: no substitution-model selection, no ML branch optimization, and bootstrap percentages are not posterior probabilities. UPGMA's clock assumption is often violated in real data. Relatedness is also not causation - shared ancestry alone can produce shared traits.",
  },
  {
    id: "assoc", navId: "assoc", label: "Pan-GWAS Association", short: "Association", color: C.pheno, icon: Network,
    tagline: "Gene presence/absence vs phenotype, corrected for multiple testing and population structure",
    input: "A gene presence/absence matrix (ROARY/Panaroo .Rtab or gene_presence_absence.csv) and a phenotype metadata table (CSV/TSV) keyed by sample id.",
    does: "Tests every feature with Fisher's exact test on its 2x2 contingency table (Welch's t-test for quantitative traits), reports Haldane-corrected odds ratios with 95% confidence intervals, adjusts p-values with Benjamini-Hochberg FDR and Bonferroni, and - when clade labels or this session's phylogeny are supplied - reruns each comparison within lineages via the Cochran-Mantel-Haenszel stratified test so inherited lineage markers stop masquerading as associations.",
    output: "A volcano plot of effect size vs evidence, a fully sortable results table, and TSV exports carrying every statistic (counts, OR + CI, raw p, FDR q, stratified CMH p) for audit or downstream plotting.",
    question: "Which genes track with this trait more than chance allows - and which of those survive once shared ancestry is accounted for?",
    limits: "Bacterial GWAS on clonal genomes is inherently confounded: even CMH stratification cannot rescue genes perfectly locked to one lineage, rare phenotypes starve the test of power, and a surviving association is still a lead rather than a mechanism. Confirm hits by linkage analysis, context curation and experiment.",
  },
];

/* =====================================================================================
   PIPELINE OVERVIEW  (top-level animated diagram + slide-over detail)
   (Panel/Eyebrow/StatCard/ExplainBox/LimitBanner/SectionTitle imported from
   ./components/ui/Primitives.jsx)
===================================================================================== */
function PipelineOverview({ explainMode, onJump }) {
  const [openStage, setOpenStage] = useState(null);
  const stage = STAGES.find((s) => s.id === openStage);

  return (
    <div>
      <SectionTitle icon={Waypoints} color={C.raw} title="Pipeline Overview"
        subtitle="Nanopore reads → cleaned data → assembled genome → annotated genes → aligned genomes → evolutionary tree → genotype/phenotype analysis → interactive browser" />

      <Panel style={{ padding: 0, position: "relative", overflow: "hidden" }}>
        <div style={{ padding: "20px 20px 24px", overflowX: "auto" }}>
        <div style={{ display: "flex", alignItems: "center", minWidth: 1080 }}>
          {STAGES.map((s, i) => (
            <React.Fragment key={s.id}>
              <button
                onClick={() => setOpenStage(s.id)}
                style={{
                  all: "unset", cursor: "pointer", display: "flex", flexDirection: "column", alignItems: "center",
                  gap: 6, width: 108, flexShrink: 0, textAlign: "center",
                }}
              >
                <div style={{
                  width: 48, height: 48, borderRadius: 2, background: `${s.color}10`,
                  border: `1px solid ${s.color}44`, display: "flex", alignItems: "center", justifyContent: "center",
                  transition: "transform .15s, box-shadow .15s",
                }}
                  onMouseEnter={(e) => { e.currentTarget.style.transform = "translateY(-2px)"; e.currentTarget.style.boxShadow = `0 0 12px ${s.color}33`; }}
                  onMouseLeave={(e) => { e.currentTarget.style.transform = "translateY(0)"; e.currentTarget.style.boxShadow = "none"; }}
                >
                  <s.icon size={20} color={s.color} />
                </div>
                <div style={{ fontSize: 10.5, color: C.text, fontFamily: FONT, lineHeight: 1.25 }}>{s.short.toLowerCase()}</div>
              </button>
              {i < STAGES.length - 1 && (
                <div style={{ position: "relative", flex: 1, height: 1, background: `linear-gradient(90deg, ${s.color}44, ${STAGES[i + 1].color}44)`, minWidth: 26, marginBottom: 20 }}>
                  <div className="flow-dot" style={{ background: s.color, boxShadow: `0 0 8px ${s.color}` }} />
                </div>
              )}
            </React.Fragment>
          ))}
        </div>
        </div>
        <div style={{ position: "absolute", top: 0, bottom: 0, right: 0, width: 48, pointerEvents: "none", background: `linear-gradient(90deg, transparent, ${C.bgPanel})` }} />
        <div style={{ position: "absolute", top: 0, bottom: 0, left: 0, width: 20, pointerEvents: "none", background: `linear-gradient(90deg, ${C.bgPanel}, transparent)` }} />
      </Panel>
      <div style={{ fontSize: 10.5, color: C.textFaint, marginTop: 6, textAlign: "right" }}>← scroll for full pipeline →</div>

      <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(230px,1fr))", gap: 10, marginTop: 14 }}>
        {STAGES.map((s) => (
          <Panel key={s.id} style={{ padding: "14px 16px", cursor: "pointer", borderLeft: `2px solid ${s.color}44` }} className="hoverlift">
            <div onClick={() => setOpenStage(s.id)}>
              <Eyebrow color={s.color}>{s.short.toLowerCase()}</Eyebrow>
              <div style={{ fontSize: 12.5, color: C.text, lineHeight: 1.45, fontWeight: 400, fontFamily: FONT }}>{s.tagline}</div>
              <div style={{ marginTop: 8, display: "flex", justifyContent: "space-between", alignItems: "center" }}>
                <span style={{ fontSize: 11, color: s.color, display: "flex", alignItems: "center", gap: 4, fontFamily: FONT }}>
                  [ details ] <ChevronRight size={11} />
                </span>
                {onJump && s.navId && (
                  <button onClick={(e) => { e.stopPropagation(); onJump(s.navId); }}
                    style={{ all: "unset", cursor: "pointer", fontSize: 10.5, color: C.textDim, border: `1px solid ${C.border}`, borderRadius: 2, padding: "2px 8px", fontFamily: FONT }}>
                    [ open ]
                  </button>
                )}
              </div>
            </div>
          </Panel>
        ))}
      </div>

      {stage && (
        <div style={{ position: "fixed", inset: 0, zIndex: 50, display: "flex", justifyContent: "flex-end" }}>
          <div onClick={() => setOpenStage(null)} style={{ position: "absolute", inset: 0, background: "rgba(3,5,8,0.7)" }} />
          <div style={{
            position: "relative", width: "min(460px, 92vw)", height: "100%", background: "#0a0e16",
            borderLeft: `2px solid ${stage.color}44`, padding: "22px 20px", overflowY: "auto",
            boxShadow: "-20px 0 40px rgba(0,0,0,0.6)",
          }}>
            <button onClick={() => setOpenStage(null)} style={{ all: "unset", cursor: "pointer", position: "absolute", top: 18, right: 18, color: C.textDim }}>
              <X size={16} />
            </button>
            <Eyebrow color={stage.color}>{stage.short.toLowerCase()}</Eyebrow>
            <div style={{ fontFamily: FONT, fontWeight: 500, fontSize: 18, color: C.text, textShadow: `0 0 10px ${stage.color}22` }}>{stage.label}</div>
            <div style={{ fontSize: 12.5, color: C.textDim, marginTop: 3, marginBottom: 18, fontFamily: FONT }}>{stage.tagline}</div>

            {[["input", stage.input], ["process", stage.does], ["output", stage.output], ["question", stage.question]].map(([h, v]) => (
              <div key={h} style={{ marginBottom: 14 }}>
                <div style={{ fontSize: 10.5, textTransform: "uppercase", letterSpacing: "0.08em", color: stage.color, marginBottom: 3, fontFamily: FONT }}>{h}</div>
                <div style={{ fontSize: 12.5, color: C.text, lineHeight: 1.55, fontFamily: FONT }}>{v}</div>
              </div>
            ))}
            <LimitBanner>{stage.limits}</LimitBanner>
            <ExplainBox explainMode={explainMode} color={stage.color}>{EXPLAIN[stage.id]}</ExplainBox>
          </div>
        </div>
      )}

      <style>{`
        .flow-dot { position:absolute; top:-3px; width:7px; height:7px; border-radius:1px;
          animation: flow 2.4s linear infinite; }
        @keyframes flow { 0% { left:0%; opacity:0 } 8% { opacity:1 } 92% { opacity:1 } 100% { left:100%; opacity:0 } }
        .hoverlift:hover { border-color: ${C.borderStrong} !important; transform: translateY(-1px); transition: all .12s; }
      `}</style>
    </div>
  );
}

const EXPLAIN = {
  fastq: "Think of a FASTQ file like a stack of Polaroids with a confidence sticker on each one: here's what the machine thinks it saw, and here's how sure it was.",
  qc: "NanoPlot is a camera, not a filter - it just shows you the picture of your data's quality so you can decide what to do next.",
  filter: "This is the actual cleanup step: reads that don't meet your bar get set aside before anything is built from them.",
  flye: "Flye reconstructs longer DNA sequences by finding where reads overlap, like reassembling a shredded document from overlapping strips of text.",
  bakta: "Annotation reads the assembled genome and points at the parts: it finds open reading frames, decides where each gene really starts using ribosome-binding-site evidence, and computes what each protein would look like — without ever claiming to know what it does.",
  msa: "Alignment lines sequences up so 'the same position' actually means the same evolutionary spot across every strain being compared - here with a real progressive algorithm: guide tree first, then closest sequences merge into growing profiles.",
  iqtree: "A phylogeny is a family tree for sequences: branch points are common ancestors, and branch lengths measure how much change accumulated. Neighbor-Joining builds it from pairwise distances, and bootstrapping re-runs the whole inference on reshuffled data to see which groupings survive.",
  assoc: "Every gene becomes a 2×2 story: present or absent, sick or healthy. The pooled test asks 'do they correlate?', and the stratified test asks the harder follow-up - 'does that correlation survive inside families of related strains?' Only hits answering yes twice deserve bench time.",
};

/* =====================================================================================
   ALL EIGHT SECTIONS ARE REAL TOOLS NOW
   1. FASTQ viewer .................... ./components/fastq/
   2. QC Dashboard .................... ./components/fastq/ (read-only diagnostics)
   3. Assembly ........................ ./components/assembly/ (OLC assembler worker)
   4. Annotation ...................... ./components/annotation/ (Prodigal-style caller)
   5. MSA ............................. ./components/msa/ (progressive aligner worker)
   6. Phylogeny ....................... ./components/phylo/ (NJ/UPGMA + bootstrap)
   7. Pan-GWAS Association ............ ./components/assoc/ (Fisher + CMH + FDR worker)
   8. Unified Genome Explorer ......... ./components/explorer/ (canvas genome browser)
   Results flow through FastqDataContext: reads -> contigs -> genes -> alignment ->
   tree, and the tree feeds back into the association tool as population structure.
 ===================================================================================== */

/* =====================================================================================
   APP ROOT
 ===================================================================================== */
const NAV = [
  { id: "overview", label: "Overview", icon: Waypoints, color: C.raw },
  { id: "fastq", label: "FASTQ", icon: ScrollText, color: C.raw },
  { id: "qc", label: "QC Dashboard", icon: Microscope, color: C.qc },
  { id: "assembly", label: "Assembly", icon: Waypoints, color: C.assembly },
  { id: "annotation", label: "Annotation", icon: Layers, color: C.annotation },
  { id: "msa", label: "Alignment", icon: Table2, color: C.phylo },
  { id: "tree", label: "Phylogeny", icon: GitBranch, color: C.phylo },
  { id: "assoc", label: "Association", icon: Network, color: C.pheno },
  { id: "explorer", label: "Explorer", icon: CircuitBoard, color: C.raw },
];

export default function App() {
  return (
    <FastqDataProvider>
      <AppShell />
    </FastqDataProvider>
  );
}

function AppShell() {
  const [section, setSection] = useState("overview");
  const [explainMode, setExplainMode] = useState(true);

  const active = NAV.find((n) => n.id === section);

  return (
    <div style={{
      background: C.bg,
      minHeight: "100%", color: C.text, fontFamily: FONT, display: "flex",
    }}>
      <style>{`
        ${FONT_IMPORT}
        * { box-sizing: border-box; }
        input[type=range] { height: 3px; border-radius: 1px; background: ${C.border}; -webkit-appearance:none; }
        input[type=range]::-webkit-slider-thumb { -webkit-appearance:none; width:12px; height:12px; border-radius:2px; background:${C.raw}; box-shadow:0 0 6px rgba(0,255,65,0.3); cursor:pointer; margin-top:-5px; }
        ::-webkit-scrollbar { width: 7px; height: 7px; }
        ::-webkit-scrollbar-thumb { background: ${C.border}; border-radius: 1px; }
        ::-webkit-scrollbar-thumb:hover { background: ${C.borderStrong}; }
        button:focus-visible, input:focus-visible { outline: 1px solid ${C.raw}; outline-offset: 1px; }
      `}</style>

      {/* Sidebar */}
      <div style={{ width: 212, flexShrink: 0, borderRight: `1px solid ${C.border}`, padding: "18px 10px", display: "flex", flexDirection: "column", gap: 2, background: "#0a0e14" }}>
        <div style={{ display: "flex", alignItems: "center", gap: 8, padding: "0 8px 18px" }}>
          <span style={{ color: C.prompt, fontFamily: FONT, fontSize: 16, fontWeight: 700 }}>&gt;_</span>
          <div style={{ fontFamily: FONT, fontWeight: 600, fontSize: 13, color: C.prompt, lineHeight: 1.3, letterSpacing: "0.04em", textShadow: "0 0 10px rgba(0,255,65,0.3)" }}>GENOME<br /><span style={{ color: C.textFaint, fontSize: 9.5, fontWeight: 400, letterSpacing: "0.1em" }}>PIPELINE EXPLORER</span></div>
        </div>
        {NAV.map((n) => (
          <button key={n.id} onClick={() => setSection(n.id)}
            style={{
              all: "unset", cursor: "pointer", display: "flex", alignItems: "center", gap: 8, padding: "6px 10px",
              borderRadius: 2,
              background: section === n.id ? `${n.color}12` : "transparent",
              borderLeft: section === n.id ? `2px solid ${n.color}` : "2px solid transparent",
              color: section === n.id ? n.color : C.textDim, fontSize: 12, fontFamily: FONT,
              transition: "all .1s",
            }}
            onMouseEnter={(e) => { if (section !== n.id) e.currentTarget.style.background = "rgba(255,255,255,0.03)"; }}
            onMouseLeave={(e) => { if (section !== n.id) e.currentTarget.style.background = "transparent"; }}
          >
            <n.icon size={12} /> {n.label.toLowerCase().replace(/\s+/g, '_')}
          </button>
        ))}
        <div style={{ marginTop: "auto", padding: "12px 8px 0", borderTop: `1px solid ${C.border}` }}>
          <label style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer", fontSize: 11, color: C.textDim, fontFamily: FONT, marginBottom: 6 }}>
            <input type="checkbox" checked={explainMode} onChange={(e) => setExplainMode(e.target.checked)} style={{ accentColor: C.prompt }} />
            --explain
          </label>
          <div style={{ fontSize: 10, color: C.textFaint, lineHeight: 1.45, fontFamily: FONT }}>
            All computation runs locally in your browser. No uploads, no accounts — your data never leaves this machine.
          </div>
        </div>
      </div>

      {/* Main */}
      <div style={{ flex: 1, minWidth: 0, padding: "24px 36px 60px" }}>
        <div style={{
          fontSize: 11, color: C.textFaint, fontFamily: FONT,
          padding: "8px 0 10px",
          borderBottom: `1px solid ${C.border}`,
          marginBottom: 18, display: "flex", justifyContent: "space-between", alignItems: "center",
        }}>
          <div>
            <span style={{ color: C.good }}>●</span>{' '}
            <span style={{ color: C.prompt }}>local</span>{' '}
            <span style={{ color: C.textFaint }}>·</span>{' '}
            {active?.label?.toLowerCase().replace(/\s+/g, '_') ?? "overview"}
          </div>
          {explainMode && <span style={{ color: C.raw }}>--explain</span>}
        </div>

        <ErrorBoundary>
          <Suspense fallback={<div style={{ padding: 40, color: C.textDim, fontFamily: FONT }}>Loading section...</div>}>
            {section === "overview" && <PipelineOverview explainMode={explainMode} onJump={setSection} />}
            {section === "fastq" && <FastqSection explainMode={explainMode} />}
            {section === "qc" && <QCDashboardSection explainMode={explainMode} />}
            {section === "assembly" && <AssemblySection explainMode={explainMode} />}
            {section === "annotation" && <AnnotationSection explainMode={explainMode} />}
            {section === "msa" && <AlignmentSection explainMode={explainMode} />}
            {section === "tree" && <PhyloSection explainMode={explainMode} />}
            {section === "assoc" && <AssocSection explainMode={explainMode} />}
            {section === "explorer" && <ExplorerSection explainMode={explainMode} />}
          </Suspense>
        </ErrorBoundary>
      </div>
    </div>
  );
}
