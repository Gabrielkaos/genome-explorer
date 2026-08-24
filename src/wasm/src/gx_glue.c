/*
 * Emscripten JS glue points (single TU - EM_JS must not be duplicated).
 */
#include "common.h"

#ifdef __EMSCRIPTEN__

/* JS side assigns globalThis.__gxWasmProgress(stageCode, pct) before a run. */
EM_JS(void, gx_progress, (int code, int pct), {
  var f = globalThis.__gxWasmProgress;
  if (f) f(code, pct);
});

/* Matches the JS engine's timing source so timings are comparable. */
EM_JS(double, gx_now, (), { return performance.now(); });

#endif
