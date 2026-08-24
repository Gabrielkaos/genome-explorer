/**
 * Incremental byte-level line scanner.
 *
 * Operates directly on raw Uint8Array chunks (rather than decoding to text
 * first) so we can track exact byte offsets - required to later re-slice
 * individual FASTQ records out of the original file on demand - and so we
 * avoid decoding megabytes of sequence/quality text we may never display.
 *
 * Usage: create one instance per parse job, call `push(chunk)` for each
 * incoming Uint8Array, and it invokes `onLine(bytes, startOffset, endOffset)`
 * for every complete line found (newline stripped). Call `flush()` at EOF to
 * emit any trailing line that wasn't newline-terminated.
 */
const NEWLINE = 0x0a; // '\n'
const CR = 0x0d; // '\r' (stripped if present before \n)

export function createLineScanner(onLine) {
  let carry = new Uint8Array(0);
  let absoluteOffset = 0; // total bytes consumed so far (start of `carry`)

  function concat(a, b) {
    const out = new Uint8Array(a.length + b.length);
    out.set(a, 0);
    out.set(b, a.length);
    return out;
  }

  function push(chunk) {
    let buf = carry.length ? concat(carry, chunk) : chunk;
    let searchStart = 0;
    while (true) {
      const nlIndex = buf.indexOf(NEWLINE, searchStart);
      if (nlIndex === -1) break;
      let end = nlIndex;
      if (end > 0 && buf[end - 1] === CR) end -= 1; // strip \r for CRLF files
      const lineBytes = buf.subarray(searchStart, end);
      const lineStartAbs = absoluteOffset + searchStart;
      const lineEndAbs = absoluteOffset + nlIndex + 1; // include the newline itself
      onLine(lineBytes, lineStartAbs, lineEndAbs);
      searchStart = nlIndex + 1;
    }
    // stash remainder as carry, rebase absoluteOffset
    carry = buf.subarray(searchStart);
    absoluteOffset += searchStart;
  }

  function flush() {
    if (carry.length > 0) {
      onLine(carry, absoluteOffset, absoluteOffset + carry.length);
      absoluteOffset += carry.length;
      carry = new Uint8Array(0);
    }
  }

  function bytesConsumed() {
    return absoluteOffset + carry.length;
  }

  return { push, flush, bytesConsumed };
}

const textDecoder = new TextDecoder("utf-8");
export function bytesToString(bytes) {
  return textDecoder.decode(bytes);
}
