import { Download, ExternalLink, ArrowRight, Globe, Home } from "lucide-react";
import { C, FONT_DISPLAY } from "../../theme.js";
import { Panel } from "../ui/Primitives.jsx";
import { useFastqData } from "../../state/FastqDataContext.jsx";
import { downloadBlob } from "../../lib/fastq/exportFastq.js";

function downloadAssemblyFasta(contigs) {
  const text = contigs
    .map((c) => `>${c.id} length=${c.length}${c.circular ? " circular=true" : ""}\n${(c.seq.match(/.{1,80}/g) || [c.seq]).join("\n")}`)
    .join("\n");
  downloadBlob(new Blob([text], { type: "text/plain" }), "assembly.fasta");
}

const DESTINATIONS = [
  {
    name: "NCBI BLAST",
    url: "https://blast.ncbi.nlm.nih.gov/",
    what: "Identify individual sequences",
    how: 'Paste one contig (choose blastn) or protein hits from a .faa file (blastp) and match them against NCBI\'s sequence databases. Best when you want to name a few specific genes rather than annotate everything at once.',
    returns: "Match descriptions with E-values, % identity and accession numbers you can look up in NCBI records.",
  },
  {
    name: "Bakta web server",
    url: "https://bakta.computational.bio/",
    what: "Full bacterial annotation with reference databases",
    how: "Upload your assembly.fasta as-is. Bakta runs gene calling plus similarity search against curated RefSeq/UniRef/IPG databases, and also finds rRNA/tRNA — the RNA features this app skips.",
    returns: "A GFF3, an annotated GenBank (.gbff), and protein/nucleotide FASTA files — every CDS carrying a real product name where a match existed.",
  },
  {
    name: "NCBI PGAP",
    url: "https://www.ncbi.nlm.nih.gov/genome/annotation_prok/",
    what: "NCBI's submission-grade prokaryotic pipeline",
    how: "The pipeline used for official RefSeq genomes. Heavier to run than Bakta but produces GenBank files acceptable for direct sequence submission.",
    returns: "Submission-ready annotated GenBank records with products backed by NCBI's protein databases.",
  },
];

export default function ExternalAnnotationGuide() {
  const { asm } = useFastqData();
  const asmAvailable = asm.status === "done" && asm.contigs?.length > 0;

  return (
    <Panel style={{ padding: 18, marginTop: 14 }}>
      <div style={{ fontSize: 12.5, color: C.textDim, lineHeight: 1.55, marginTop: 4 }}>
        Only similarity search against public reference databases can say <em>what genes do</em>. Browsers can't
        host those multi-GB databases — so use a free server together with this app. Three steps, five minutes:
      </div>

      {/* Step 1 */}
      <div style={{ display: "flex", gap: 12, marginTop: 16 }}>
        <StepNum n={1} />
        <div style={{ flex: 1 }}>
          <StepHead title="Export your assembled genome" />
          <div style={{ fontSize: 12.5, color: C.textDim, lineHeight: 1.55, marginTop: 3 }}>
            {asmAvailable ? (
              <>
                This session's Assembly produced <strong style={{ color: C.assembly }}>{asm.contigs.length} contig(s)</strong> — grab them right here:
              </>
            ) : (
              <>
                Run the <strong style={{ color: C.assembly }}>Assembly</strong> section first (sidebar), then come back — its
                contigs feed this page automatically, and its <em>Download FASTA</em> button exports <code style={codeStyle}>assembly.fasta</code>.
              </>
            )}
          </div>
          {asmAvailable && (
            <button onClick={() => downloadAssemblyFasta(asm.contigs)} style={{
              all: "unset", cursor: "pointer", display: "inline-flex", alignItems: "center", gap: 6,
              marginTop: 8, fontSize: 12, color: C.good, border: `1px solid ${C.good}66`, borderRadius: 6, padding: "6px 12px",
            }}>
              <Download size={12} /> Download assembly.fasta ({asm.contigs.length} contig{asm.contigs.length !== 1 ? "s" : ""})
            </button>
          )}
        </div>
      </div>

      {/* Step 2 */}
      <div style={{ display: "flex", gap: 12, marginTop: 18 }}>
        <StepNum n={2} />
        <div style={{ flex: 1, minWidth: 0 }}>
          <StepHead icon={Globe} title="Run it against a public reference database" />
          <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(230px, 1fr))", gap: 10, marginTop: 8 }}>
            {DESTINATIONS.map((d) => (
              <a key={d.name} href={d.url} target="_blank" rel="noreferrer" style={{
                textDecoration: "none", display: "block", background: "#05070a", border: `1px solid ${C.border}`,
                borderRadius: 8, padding: "12px 13px", transition: "border-color .15s", cursor: "pointer",
              }}
                onMouseEnter={(e) => (e.currentTarget.style.borderColor = `${C.raw}88`)}
                onMouseLeave={(e) => (e.currentTarget.style.borderColor = C.border)}>
                <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between" }}>
                  <span style={{ fontFamily: FONT_DISPLAY, fontSize: 12.5, color: C.raw }}>{d.name}</span>
                  <ExternalLink size={11} color={C.textFaint} />
                </div>
                <div style={{ fontSize: 11, color: C.annotation, marginTop: 5 }}>{d.what}</div>
                <div style={{ fontSize: 11.5, color: C.textDim, lineHeight: 1.5, marginTop: 5 }}>{d.how}</div>
                <div style={{ fontSize: 11, color: C.textFaint, lineHeight: 1.5, marginTop: 6 }}>
                  You get back: {d.returns}
                </div>
              </a>
            ))}
          </div>
          <div style={{ fontSize: 11, color: C.textFaint, marginTop: 8, lineHeight: 1.5 }}>
            Note: these are external websites — uploading there sends your sequences over the network, unlike everything else in this app which stays fully local.
          </div>
        </div>
      </div>

      {/* Step 3 */}
      <div style={{ display: "flex", gap: 12, marginTop: 18 }}>
        <StepNum n={3} />
        <div style={{ flex: 1 }}>
          <StepHead icon={Home} title="Come back and load the results here" />
          <div style={{ fontSize: 12.5, color: C.textDim, lineHeight: 1.55, marginTop: 3 }}>
            Open the <strong style={{ color: C.raw }}>Explorer</strong> section in the sidebar and pick{" "}
            <strong>Uploaded files</strong>. Either works:
          </div>
          <ul style={{ margin: "6px 0 0", paddingLeft: 18, fontSize: 12.5, color: C.textDim, lineHeight: 1.7 }}>
            <li>The self-contained GenBank file (<code style={codeStyle}>.gbff</code>/<code style={codeStyle}>.gbk</code>) from Bakta or PGAP — sequence + annotations in one file.</li>
            <li>Your <code style={codeStyle}>assembly.fasta</code> as the genome, then their <code style={codeStyle}>.gff3</code> as the annotation file.</li>
          </ul>
          <div style={{ fontSize: 12.5, color: C.textDim, lineHeight: 1.55, marginTop: 6 }}>
            Your genes render on this app's interactive canvas with real product names. One honest caveat: this
            round trip is a <strong>viewing</strong> branch only — it does not feed the rest of the pipeline. Alignment
            wants multiple homologous sequences (e.g. the same gene from several strains), Phylogeny consumes that
            alignment, and Association expects a ROARY/Panaroo-style presence/absence matrix plus a phenotype table
            across many samples — none of those read GFF3/GenBank annotation output.
          </div>
        </div>
      </div>

      <div style={{ display: "flex", alignItems: "center", gap: 8, marginTop: 16, fontSize: 11, color: C.textFaint }}>
        <ArrowRight size={11} /> Local gene calling here → public database names the functions → back into the Explorer browser (viewing only).
      </div>
    </Panel>
  );
}

function StepNum({ n }) {
  return (
    <div style={{
      width: 22, height: 22, flexShrink: 0, borderRadius: 999, background: `${C.raw}18`,
      border: `1px solid ${C.raw}55`, color: C.raw, fontFamily: FONT_DISPLAY,
      fontSize: 11.5, display: "flex", alignItems: "center", justifyContent: "center", marginTop: 1,
    }}>{n}</div>
  );
}

function StepHead({ icon: Icon, title }) {
  return (
    <div style={{ display: "flex", alignItems: "center", gap: 7 }}>
      {Icon && <Icon size={13} color={C.raw} />}
      <span style={{ fontSize: 13, color: C.text, fontWeight: 500 }}>{title}</span>
    </div>
  );
}

const codeStyle = {
  fontFamily: FONT_DISPLAY, fontSize: 11, background: "#05070a",
  border: `1px solid ${C.border}`, borderRadius: 4, padding: "1px 5px", color: C.textDim,
};
