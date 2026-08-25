/* eslint-disable react-refresh/only-export-components */
import { createContext, useCallback, useContext, useState } from "react";
import { useFastqParser } from "../hooks/useFastqParser.js";
import { useAssembly } from "../hooks/useAssembly.js";
import { useAnnotation } from "../hooks/useAnnotation.js";
import { useMsa } from "../hooks/useMsa.js";
import { usePhylo } from "../hooks/usePhylo.js";

const FastqDataContext = createContext(null);

/**
 * Wraps the app so any section can read the currently-loaded FASTQ dataset
 * without re-parsing or re-uploading. QC Dashboard is intentionally
 * read-only here: it visualizes whatever was loaded on the FASTQ page (or
 * loaded directly from its own uploader), but never filters anything - that
 * distinction is a deliberate, recurring point in this app.
 *
 * Assembly results also live here (not inside AssemblySection's local state)
 * so downstream stages - currently Annotation - can
 * consume the contigs without re-running the assembler when navigating
 * between sections. The progressive MSA result (`msa`)
 * is kept here so Phylogeny can consume the aligned matrix directly, and
 * the inferred phylogeny (`phyl`) survives navigation in turn. Annotation
 * (`ann`) is shared as well: the Genome Explorer renders its feature calls.
 */
export function FastqDataProvider({ children }) {
  const api = useFastqParser();
  const asm = useAssembly();
  const ann = useAnnotation();
  const msa = useMsa();
  const phyl = usePhylo();
  const [file, setFile] = useState(null);

  const parseFile = useCallback((f) => {
    setFile(f);
    api.parseFile(f);
  }, [api]);

  /** Batch-fetch full sequence text for many read indices (used by Assembly). */
  const getManyRecords = useCallback(async (indices, { concurrency = 80 } = {}) => {
    const out = new Array(indices.length);
    let cursor = 0;
    async function worker() {
      while (cursor < indices.length) {
        const myIdx = cursor++;
        out[myIdx] = await api.getRecord(indices[myIdx]);
      }
    }
    await Promise.all(Array.from({ length: Math.min(concurrency, indices.length) }, worker));
    return out;
  }, [api]);

  const value = { ...api, file, parseFile, getManyRecords, asm, ann, msa, phyl };
  return <FastqDataContext.Provider value={value}>{children}</FastqDataContext.Provider>;
}

export function useFastqData() {
  const ctx = useContext(FastqDataContext);
  if (!ctx) throw new Error("useFastqData must be used within a FastqDataProvider");
  return ctx;
}
