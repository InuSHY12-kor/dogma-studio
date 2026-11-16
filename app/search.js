/* GeneCanvas — sequence search. Exposes window.SEARCH
 *
 * The query syntax:
 *   GGTCT(AG)CTTC   a set of allowed letters at one position
 *   GGTCTRCTTC      an IUPAC ambiguity code, expanded to the letters it stands for
 *   !R  !(AG)       letters that must NOT appear at that position
 *   "GGTCTRCTTC"    quoted: a degenerate letter matches only itself, not what it stands for
 * and three kinds of search: exact, similar (mismatches, optionally with gaps), and a
 * protein search over translations of the sequence.
 */
(function (root) {
  'use strict';
  const DNA = root.DNA;

  const NT = { A: 'A', C: 'C', G: 'G', T: 'T', U: 'T', R: 'AG', Y: 'CT', S: 'GC', W: 'AT', K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT' };
  const AAX = { X: 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*', B: 'DN', Z: 'EQ' };
  const ALL_AA = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*';

  /* Parse a query into one token per position.
   * Each token is { allow } or { deny }, holding a set of letters. */
  function parseQuery(q, kind) {
    let s = String(q || '').trim();
    let literal = false;
    if (s.length > 1 && s[0] === '"' && s[s.length - 1] === '"') { literal = true; s = s.slice(1, -1); }
    const table = kind === 'protein' ? AAX : NT;
    const universe = kind === 'protein' ? ALL_AA : 'ACGT';
    const up = s.toUpperCase().replace(/\s+/g, '');
    const toks = [];
    let i = 0;
    const letters = group => {
      // the letters a group of query characters stands for. Unquoted, a degenerate code
      // matches the letters it stands for AND a base literally written as that code,
      // which is the behaviour the guide describes for R versus "R".
      let out = '';
      for (const ch of group) out += literal ? ch : ((table[ch] || ch) + ch);
      return new Set(out.split(''));
    };
    while (i < up.length) {
      let deny = false;
      if (up[i] === '!') { deny = true; i++; if (i >= up.length) break; }
      let group;
      if (up[i] === '(') {
        const close = up.indexOf(')', i);
        if (close < 0) return null;
        group = up.slice(i + 1, close);
        i = close + 1;
      } else { group = up[i]; i++; }
      if (!group) return null;
      const set = letters(group);
      // an excluded degenerate letter also excludes the letter itself, as the guide says
      if (deny) for (const ch of group) set.add(ch);
      if (deny) {
        const allow = new Set();
        for (const ch of universe) if (!set.has(ch)) allow.add(ch);
        toks.push({ allow, deny: set });
      } else {
        toks.push({ allow: set });
      }
    }
    if (!toks.length) return null;
    return { toks, literal, kind: kind || 'dna' };
  }
  const hit = (tok, ch) => tok.allow.has(ch);

  // The reverse complement of a parsed DNA query, for searching the bottom strand.
  function rcQuery(p) {
    const comp = { A: 'T', C: 'G', G: 'C', T: 'A' };
    const toks = p.toks.slice().reverse().map(t => {
      const allow = new Set();
      for (const ch of t.allow) allow.add(comp[ch] || ch);
      return { allow, deny: t.deny };
    });
    return { toks, literal: p.literal, kind: p.kind };
  }

  function matchAt(text, at, toks) {
    for (let j = 0; j < toks.length; j++) if (!hit(toks[j], text[at + j])) return false;
    return true;
  }

  /* An exact search, on both strands of a DNA sequence. A query that reads the same
   * on both strands is reported once. */
  function exact(seq, circular, query) {
    const p = parseQuery(query, 'dna');
    if (!p) return [];
    const n = seq.length, L = p.toks.length;
    if (!L || L > n) return [];
    const text = circular ? seq + seq.slice(0, L - 1) : seq;
    const res = [], seen = new Set();
    const push = (start, strand) => {
      if (start >= n) return;
      const k = start + ':' + strand;
      if (seen.has(k)) return;
      seen.add(k);
      let end = start + L; if (end > n) end -= n;
      res.push({ start, end, strand, errors: 0, mismatches: [] });
    };
    for (let i = 0; i + L <= text.length; i++) if (matchAt(text, i, p.toks)) push(i, 1);
    const rp = rcQuery(p);
    const same = rp.toks.length === p.toks.length && rp.toks.every((t, j) => {
      const a = [...t.allow].sort().join(''), b = [...p.toks[j].allow].sort().join('');
      return a === b;
    });
    if (!same) for (let i = 0; i + L <= text.length; i++) if (matchAt(text, i, rp.toks)) push(i, -1);
    return res.sort((a, b) => a.start - b.start);
  }

  /* A similar search: up to `maxErr` differences. With gaps allowed this is a banded
   * edit distance; without them it is a Hamming count, which also records which
   * positions differ so the sequence view can box them. */
  function similar(seq, circular, query, opts = {}) {
    const p = parseQuery(query, 'dna');
    if (!p) return [];
    const maxErr = Math.max(0, opts.maxErr == null ? 1 : opts.maxErr);
    const gaps = !!opts.gaps;
    const n = seq.length, L = p.toks.length;
    if (!L || L > n + maxErr) return [];
    const span = L + maxErr;
    const text = circular ? seq + seq.slice(0, Math.min(n, span)) : seq;
    const out = [];
    const run = (toks, strand) => {
      for (let i = 0; i < (circular ? n : Math.max(0, n - L + maxErr + 1)); i++) {
        if (!gaps) {
          if (i + L > text.length) break;
          let err = 0; const mm = [];
          for (let j = 0; j < L && err <= maxErr; j++) {
            if (!hit(toks[j], text[i + j])) { err++; mm.push(j); }
          }
          if (err <= maxErr) out.push({ start: i, len: L, strand, errors: err, mismatches: mm });
        } else {
          const r = bandedEdit(text, i, toks, maxErr);
          if (r) out.push({ start: i, len: r.len, strand, errors: r.err, mismatches: [] });
        }
      }
    };
    run(p.toks, 1);
    const rp = rcQuery(p);
    run(rp.toks, -1);
    // keep the best hit of each overlapping cluster, per strand
    out.sort((a, b) => a.start - b.start || a.errors - b.errors);
    const keep = [];
    for (const h of out) {
      const prev = keep.find(k => k.strand === h.strand && h.start < k.start + k.len && k.start < h.start + h.len);
      if (!prev) { keep.push(h); continue; }
      if (h.errors < prev.errors) Object.assign(prev, h);
    }
    return keep.filter(h => h.start < n).map(h => {
      let end = h.start + h.len; if (end > n) end -= n;
      return { start: h.start % n, end, strand: h.strand, errors: h.errors, mismatches: h.mismatches };
    }).sort((a, b) => a.start - b.start);
  }

  // Edit distance of the query against text starting at `at`, allowing at most k
  // differences, computed over a band of width 2k+1. Returns the shortest match found.
  function bandedEdit(text, at, toks, k) {
    const L = toks.length;
    let prev = new Array(L + 1);
    for (let j = 0; j <= L; j++) prev[j] = j;          // deleting query letters
    let best = null;
    for (let i = 1; i <= L + k; i++) {
      const ch = text[at + i - 1];
      if (ch === undefined) break;
      const cur = new Array(L + 1);
      cur[0] = i;
      const lo = Math.max(1, i - k), hi = Math.min(L, i + k);
      for (let j = 0; j < lo; j++) cur[j] = j === 0 ? i : Infinity;
      for (let j = lo; j <= hi; j++) {
        const sub = prev[j - 1] + (hit(toks[j - 1], ch) ? 0 : 1);
        const del = (prev[j] === undefined ? Infinity : prev[j]) + 1;
        const ins = (cur[j - 1] === undefined ? Infinity : cur[j - 1]) + 1;
        cur[j] = Math.min(sub, del, ins);
      }
      for (let j = hi + 1; j <= L; j++) cur[j] = Infinity;
      if (cur[L] <= k && (!best || cur[L] < best.err)) best = { err: cur[L], len: i };
      prev = cur;
    }
    return best;
  }

  /* A protein search over the translations of a DNA sequence.
   *   opts.frames  which frames to search: 'all6' (default), 'top3', 'bottom3' or a list
   *   opts.code    genetic code
   *   opts.regions limit the search to these ranges (visible translations) */
  function protein(seq, circular, query, opts = {}) {
    const p = parseQuery(query, 'protein');
    if (!p) return [];
    const L = p.toks.length, bp = L * 3;
    const n = seq.length;
    if (!L || bp > n) return [];
    const table = (DNA.GENETIC_CODES && DNA.GENETIC_CODES[opts.code]) || DNA.CODON;
    const ext = circular ? seq + seq.slice(0, Math.min(n, bp)) : seq;
    const res = [], seen = new Set();
    const push = (start, strand) => {
      if (start < 0 || start >= n) return;
      const k = start + ':' + strand;
      if (seen.has(k)) return;
      seen.add(k);
      let end = start + bp; if (end > n) end -= n;
      res.push({ start, end, strand, errors: 0, mismatches: [], protein: true });
    };
    const frames = opts.frames === 'top3' ? [1] : (opts.frames === 'bottom3' ? [-1] : [1, -1]);
    for (const strand of frames) {
      const src = strand === 1 ? ext : DNA.rc(ext);
      for (let f = 0; f < 3; f++) {
        // translate once per frame, then walk the residues
        let aa = '';
        for (let i = f; i + 3 <= src.length; i += 3) aa += table[src.substr(i, 3)] || 'X';
        for (let r = 0; r + L <= aa.length; r++) {
          let ok = true;
          for (let j = 0; j < L; j++) if (!hit(p.toks[j], aa[r + j])) { ok = false; break; }
          if (!ok) continue;
          const offInSrc = f + r * 3;
          if (strand === 1) push(offInSrc, 1);
          else push(src.length - (offInSrc + bp), -1);
        }
      }
    }
    let out = res.sort((a, b) => a.start - b.start);
    if (opts.regions && opts.regions.length) {
      out = out.filter(h => opts.regions.some(rg => DNA.overlaps(h, rg, n)));
    }
    return out;
  }

  /* A search by name over the enzymes, features and primers of a document.
   * `sources` is built by the caller, since only it knows what is currently visible. */
  function byName(sources, term) {
    const q = String(term || '').trim().toLowerCase();
    if (!q) return [];
    return sources.filter(s => (s.name || '').toLowerCase().includes(q))
      .sort((a, b) => a.start - b.start);
  }

  /* Feature metadata search: name, type, note and every qualifier. Returns the
   * matching features with the fields that matched, for the Features view. */
  function featureMeta(features, term) {
    const q = String(term || '').trim().toLowerCase();
    if (!q) return [];
    const out = [];
    for (const f of features) {
      const where = [];
      if ((f.name || '').toLowerCase().includes(q)) where.push('name');
      if ((f.type || '').toLowerCase().includes(q)) where.push('type');
      if ((f.notes || '').toLowerCase().includes(q)) where.push('/note');
      for (const [k, v] of Object.entries(f.quals || {})) {
        if (String(v).toLowerCase().includes(q) || k.toLowerCase().includes(q)) where.push('/' + k);
      }
      if (where.length) out.push({ feature: f, where });
    }
    return out;
  }

  root.SEARCH = { parseQuery, exact, similar, protein, byName, featureMeta };
})(window);
