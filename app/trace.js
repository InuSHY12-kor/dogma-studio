/* GeneCanvas — Sanger trace (.ab1 / .scf) reading and alignment to a reference. Exposes window.TRACE */
(function (root) {
  'use strict';
  const DNA = root.DNA;

  // ---------------------------------------------------------------- ABIF (.ab1) reader
  const ETYPE = { BYTE: 1, CHAR: 2, WORD: 3, SHORT: 4, LONG: 5, FLOAT: 7, DOUBLE: 8, DATE: 10, TIME: 11, PSTRING: 18, CSTRING: 19 };

  function readEntry(dv, off) {
    const name = String.fromCharCode(dv.getUint8(off), dv.getUint8(off + 1), dv.getUint8(off + 2), dv.getUint8(off + 3));
    return {
      name, number: dv.getInt32(off + 4), type: dv.getUint16(off + 8), size: dv.getUint16(off + 10),
      count: dv.getInt32(off + 12), dataSize: dv.getInt32(off + 16), dataOffset: dv.getInt32(off + 20),
      inlineAt: off + 20
    };
  }
  function entryData(dv, e) {
    // Values of four bytes or fewer are stored in the offset field itself.
    const at = e.dataSize <= 4 ? e.inlineAt : e.dataOffset;
    const out = [];
    switch (e.type) {
      case ETYPE.BYTE: for (let i = 0; i < e.count; i++) out.push(dv.getUint8(at + i)); return out;
      case ETYPE.CHAR: { let s = ''; for (let i = 0; i < e.count; i++) s += String.fromCharCode(dv.getUint8(at + i)); return s; }
      case ETYPE.WORD: for (let i = 0; i < e.count; i++) out.push(dv.getUint16(at + i * 2)); return out;
      case ETYPE.SHORT: for (let i = 0; i < e.count; i++) out.push(dv.getInt16(at + i * 2)); return out;
      case ETYPE.LONG: for (let i = 0; i < e.count; i++) out.push(dv.getInt32(at + i * 4)); return out;
      case ETYPE.FLOAT: for (let i = 0; i < e.count; i++) out.push(dv.getFloat32(at + i * 4)); return out;
      case ETYPE.PSTRING: { const L = dv.getUint8(at); let s = ''; for (let i = 1; i <= L; i++) s += String.fromCharCode(dv.getUint8(at + i)); return s; }
      case ETYPE.CSTRING: { let s = ''; for (let i = 0; i < e.count - 1; i++) s += String.fromCharCode(dv.getUint8(at + i)); return s; }
      default: return null;
    }
  }

  function parseAB1(buf, filename) {
    const dv = new DataView(buf);
    const magic = String.fromCharCode(dv.getUint8(0), dv.getUint8(1), dv.getUint8(2), dv.getUint8(3));
    if (magic !== 'ABIF') throw new Error('This is not an ABIF (.ab1) file.');
    const rootDir = readEntry(dv, 6);
    const n = rootDir.count, start = rootDir.dataOffset;
    const tags = {};
    for (let i = 0; i < n; i++) {
      const e = readEntry(dv, start + i * 28);
      const key = e.name + e.number;
      if (!(key in tags)) tags[key] = entryData(dv, e);
    }
    const bases = (tags.PBAS2 || tags.PBAS1 || '').toUpperCase().replace(/[^ACGTN]/g, 'N');
    if (!bases) throw new Error('The file has no base calls.');
    const quals = tags.PCON2 || tags.PCON1 || [];
    const peaks = tags.PLOC2 || tags.PLOC1 || [];
    const order = (tags.FWO_1 || 'GATC').toUpperCase();
    const traces = {};
    // DATA 9-12 hold the analysed traces, in the channel order given by FWO_
    const ch = [tags.DATA9, tags.DATA10, tags.DATA11, tags.DATA12];
    for (let i = 0; i < 4; i++) if (ch[i]) traces[order[i] || 'ACGT'[i]] = ch[i];
    return {
      kind: 'ab1',
      name: (tags.SMPL1 || (filename || 'trace').replace(/\.[^.]+$/, '')).trim() || 'trace',
      file: filename || '',
      seq: bases, quals, peaks, traces,
      hasTrace: Object.keys(traces).length === 4 && peaks.length === bases.length,
      machine: tags.MODL1 || '', ran: tags.RUND1 || null
    };
  }

  // A plain FASTA/text read, so alignment also works without a chromatogram.
  function parsePlainRead(text, filename) {
    const r = DNA.parseAny(text, filename);
    return { kind: 'seq', name: r.name, file: filename || '', seq: r.seq, quals: [], peaks: [], traces: {}, hasTrace: false };
  }
  function parseRead(data, filename) {
    if (data instanceof ArrayBuffer) {
      const dv = new DataView(data);
      const magic = String.fromCharCode(dv.getUint8(0), dv.getUint8(1), dv.getUint8(2), dv.getUint8(3));
      if (magic === 'ABIF') return parseAB1(data, filename);
      return parsePlainRead(new TextDecoder().decode(data), filename);
    }
    return parsePlainRead(data, filename);
  }

  // ---------------------------------------------------------------- quality trimming
  // Trims the low-quality run at each end using a sliding mean, as sequencing viewers do.
  function trimRange(quals, win = 20, cutoff = 20) {
    const m = quals.length;
    if (!m) return { start: 0, end: 0, trimmed: false };
    if (m <= win) return { start: 0, end: m, trimmed: false };
    const mean = (a, b) => { let s = 0; for (let i = a; i < b; i++) s += quals[i] || 0; return s / (b - a); };
    let start = 0;
    while (start + win <= m && mean(start, start + win) < cutoff) start++;
    let end = m;
    while (end - win > start && mean(end - win, end) < cutoff) end--;
    if (start >= end) return { start: 0, end: m, trimmed: false };
    return { start, end, trimmed: start > 0 || end < m };
  }

  // ---------------------------------------------------------------- alignment
  const K = 12;
  function seedOffsets(ref, read) {
    // Vote for the diagonal that the most shared 12-mers agree on.
    const idx = new Map();
    for (let i = 0; i + K <= ref.length; i++) {
      const k = ref.substr(i, K);
      let a = idx.get(k); if (!a) idx.set(k, a = []);
      if (a.length < 40) a.push(i);
    }
    const votes = new Map();
    for (let j = 0; j + K <= read.length; j += 3) {
      const hits = idx.get(read.substr(j, K)); if (!hits) continue;
      for (const i of hits) { const d = i - j; votes.set(d, (votes.get(d) || 0) + 1); }
    }
    if (!votes.size) return [];
    // merge nearby diagonals (indels shift them slightly)
    const sorted = [...votes.entries()].sort((a, b) => b[1] - a[1]);
    const picked = [];
    for (const [d, v] of sorted) {
      if (picked.some(p => Math.abs(p.d - d) < 60)) continue;
      picked.push({ d, votes: v });
      if (picked.length >= 3) break;
    }
    return picked;
  }

  const MATCH = 2, MISMATCH = -4, GAP = -6;
  /* Semi-global: the whole read is aligned, the reference may overhang at both ends.
   * Only a band around the seed diagonal is filled, since a sequencing read never
   * drifts far from it; `band` is how many bases of indel are allowed either way.
   * That turns the work from read length x window length into read length x band. */
  const NEG = -1e9;
  function alignWindow(refWin, read, band, diag) {
    const m = read.length, N = refWin.length;
    const W = N + 1;
    const wide = !band || band * 2 >= N;
    const score = new Int32Array((m + 1) * W);
    const ptr = new Uint8Array((m + 1) * W); // 0 diag, 1 up (gap in ref), 2 left (gap in read)
    if (!wide) score.fill(NEG);
    for (let j = 0; j <= N; j++) { score[j] = 0; ptr[j] = 2; }       // free leading reference
    for (let i = 1; i <= m; i++) { score[i * W] = GAP * i; ptr[i * W] = 1; }
    for (let i = 1; i <= m; i++) {
      const rowP = (i - 1) * W, row = i * W, rc = read.charCodeAt(i - 1);
      // the band this row covers, clamped to the window
      const centre = wide ? 0 : i + (diag || 0);
      const lo = wide ? 1 : Math.max(1, centre - band);
      const hi = wide ? N : Math.min(N, centre + band);
      for (let j = lo; j <= hi; j++) {
        const dj = score[rowP + j - 1];
        const d = dj <= NEG ? NEG : dj + (refWin.charCodeAt(j - 1) === rc ? MATCH : MISMATCH);
        const uj = score[rowP + j], u = uj <= NEG ? NEG : uj + GAP;
        const lj = score[row + j - 1], l = lj <= NEG ? NEG : lj + GAP;
        let best = d, p = 0;
        if (u > best) { best = u; p = 1; }
        if (l > best) { best = l; p = 2; }
        score[row + j] = best; ptr[row + j] = p;
      }
    }
    // best end position along the last row (free trailing reference)
    let endJ = 0, best = -Infinity;
    for (let j = 0; j <= N; j++) { const v = score[m * W + j]; if (v > NEG && v >= best) { best = v; endJ = j; } }
    if (best === -Infinity) return null;
    let i = m, j = endJ, A = '', B = '';
    let guard = (m + N) * 2 + 10;
    while (i > 0 && guard-- > 0) {
      const p = ptr[i * W + j];
      if (p === 0 && j > 0) { A += refWin[j - 1]; B += read[i - 1]; i--; j--; }
      else if (p === 1 || j === 0) { A += '-'; B += read[i - 1]; i--; }
      else { A += refWin[j - 1]; B += '-'; j--; }
    }
    return { refAln: [...A].reverse().join(''), readAln: [...B].reverse().join(''), startJ: j, endJ, score: best };
  }

  // Aligns one read against a document. Tries both strands and returns the better placement.
  function align(doc, read, opts = {}) {
    const refFull = doc.circular ? doc.seq + doc.seq : doc.seq;
    const n = doc.seq.length;
    const trim = opts.trim !== false && read.quals.length ? trimRange(read.quals) : { start: 0, end: read.seq.length, trimmed: false };
    const sub = read.seq.slice(trim.start, trim.end);
    if (sub.length < 25) throw new Error(`“${read.name}” has too little usable sequence to align.`);
    // every candidate placement from both strands, strongest seed first; the first one
    // that is clearly good enough ends the search, which is the usual case
    const cands = [];
    for (const strand of [1, -1]) {
      const q = strand === 1 ? sub : DNA.rc(sub);
      for (const seed of seedOffsets(refFull, q)) cands.push({ strand, q, seed });
    }
    cands.sort((a, b) => b.seed.votes - a.seed.votes);
    let best = null;
    const perfect = sub.length * MATCH;
    for (const c of cands.slice(0, 4)) {
      const q = c.q;
      const margin = Math.max(60, Math.round(q.length * 0.15));
      const from = Math.max(0, c.seed.d - margin), to = Math.min(refFull.length, c.seed.d + q.length + margin);
      if (to - from < q.length * 0.5) continue;
      const a = alignWindow(refFull.slice(from, to), q, margin, c.seed.d - from);
      if (!a) continue;
      if (!best || a.score > best.a.score) best = { a, strand: c.strand, from, q };
      // 90% of the best possible score means the placement is certain
      if (best.a.score > perfect * 0.9) break;
    }
    if (!best) throw new Error(`“${read.name}” does not match this sequence.`);
    const { a, strand, from, q } = best;
    const refStart = from + a.startJ, refEnd = from + a.endJ;
    // walk the alignment to collect differences in reference coordinates
    const diffs = [];
    let ri = refStart, qi = 0, match = 0;
    for (let k = 0; k < a.refAln.length; k++) {
      const r = a.refAln[k], s = a.readAln[k];
      if (r === '-') { diffs.push({ pos: ri % n, type: 'ins', ref: '', read: s, readIdx: qi }); qi++; }
      else if (s === '-') { diffs.push({ pos: ri % n, type: 'del', ref: r, read: '', readIdx: qi }); ri++; }
      else { if (r === s) match++; else diffs.push({ pos: ri % n, type: 'mis', ref: r, read: s, readIdx: qi }); ri++; qi++; }
    }
    const cols = a.refAln.length;
    // map each aligned read base back to its index in the original (untrimmed) read
    const toRaw = i => strand === 1 ? trim.start + i : trim.start + (q.length - 1 - i);
    return {
      read, strand, refStart: refStart % n, refEnd: ((refEnd - 1) % n) + 1, wraps: refEnd > n,
      refAln: a.refAln, readAln: a.readAln, diffs, identity: cols ? match / cols * 100 : 0,
      aligned: q.length, trim, score: a.score, toRaw,
      coverage: { start: refStart % n, end: ((refEnd - 1) % n) + 1, len: refEnd - refStart }
    };
  }

  // Consensus differences shared by several reads, which is what tells a real variant from a read error.
  function summarize(alns) {
    const byPos = new Map();
    for (const al of alns) for (const d of al.diffs) {
      const k = d.pos + ':' + d.type + ':' + d.read;
      let e = byPos.get(k); if (!e) byPos.set(k, e = { ...d, reads: [] });
      e.reads.push(al.read.name);
    }
    return [...byPos.values()].sort((a, b) => a.pos - b.pos);
  }

  root.TRACE = { parseAB1, parseRead, parsePlainRead, trimRange, align, summarize, alignWindow, seedOffsets };
})(typeof window !== 'undefined' ? window : globalThis);
