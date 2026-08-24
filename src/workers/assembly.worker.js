import { runAssembly } from "../lib/engine/assemblyEngine.js";

self.onmessage = async (e) => {
  if (e.data?.type !== "assemble") return;
  try {
    const { reads, params } = { reads: e.data.reads, params: e.data.params || {} };
    const { contigs, stats } = await runAssembly(reads, params, (m) =>
      self.postMessage({ type: "progress", ...m })
    );

    // strip raw sequence from members before sending back (UI re-derives what
    // it needs from contig.seq + reads list already on main thread)
    const contigsOut = contigs.map((c, idx) => ({
      id: `contig_${idx + 1}`,
      seq: c.seq,
      length: c.length,
      circular: c.circular,
      readCount: c.readCount,
      members: c.members.map((m) => ({ readIdx: m.readIdx, strand: m.strand, contigStart: m.contigStart, contigEnd: m.contigEnd })),
    }));

    self.postMessage({ type: "done", contigs: contigsOut, stats });
  } catch (err) {
    self.postMessage({ type: "error", message: String(err?.message || err) });
  }
};
