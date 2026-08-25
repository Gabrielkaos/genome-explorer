import { Layers } from "lucide-react";
import { C } from "../../theme.js";
import { SectionTitle, LimitBanner } from "../ui/Primitives.jsx";
import ExternalAnnotationGuide from "./ExternalAnnotationGuide.jsx";

export default function AnnotationSection() {
  return (
    <div>
      <SectionTitle icon={Layers} color={C.annotation} title="Genome Annotation"
        subtitle="De novo bacterial gene calling with protein-level analysis — runs entirely in your browser" />

      <LimitBanner>
        This performs genuine unsupervised gene prediction (hexamer coding-model training, six-frame ORF scanning,
        Shine-Dalgarno-scored start selection, overlap resolution) plus computed protein properties on your real
        contigs — architecturally a lightweight Prodigal-style caller. It deliberately does <strong>not</strong> do
        similarity search against reference databases (so every product is honestly "hypothetical protein") and does
        <strong> not</strong> call RNA features (rRNA/tRNA require covariance models — use Barrnap or tRNAscan-SE for those).
        For submission-grade annotation, take the exported GenBank/GFF3 files to a full Bakta/PGAP server.
      </LimitBanner>

      <ExternalAnnotationGuide />
    </div>
  );
}
