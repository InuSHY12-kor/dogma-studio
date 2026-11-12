/* GeneCanvas — primer design. Exposes window.PRIMER */
(function (root) {
  'use strict';
  const DNA = root.DNA;

  const DEF = {
    minLen: 18, maxLen: 30, optLen: 20,
    minTm: 57, maxTm: 63, optTm: 60,
    minGC: 40, maxGC: 60,
    maxShift: 15,        // how far the outer end may move to find a better primer
    maxTmDiff: 3,        // between the two primers of a pair
    gcClamp: true,
    na: 0.05, conc: 250e-9
  };

  const comp = { A: 'T', T: 'A', G: 'C', C: 'G' };
  const pairs = (a, b) => comp[a] === b;

  // Longest run of complementary bases when two primers are slid past each other.
  function maxDimer(a, b) {
    let best = 0, where = null;
    const rb = [...b].reverse().join('');   // 3'→5' so that a 5'→3' lines up against it
    for (let off = -(a.length - 1); off < rb.length; off++) {
      let run = 0;
      for (let i = 0; i < a.length; i++) {
        const j = i + off;
        if (j < 0 || j >= rb.length) { run = 0; continue; }
        if (pairs(a[i], rb[j])) { run++; if (run > best) { best = run; where = { off, end: i }; } }
        else run = 0;
      }
    }
    return { len: best, at: where };
  }
  // Complementarity involving the last five bases, which is what actually primes a dimer.
  function endDimer(a, b) {
    const tail = a.slice(-5);
    let best = 0;
    const rb = [...b].reverse().join('');
    for (let off = 0; off <= rb.length - 1; off++) {
      let run = 0;
      for (let i = 0; i < tail.length; i++) {
        const j = i + off; if (j >= rb.length) break;
        if (pairs(tail[i], rb[j])) run++; else run = 0;
        if (run > best) best = run;
      }
    }
    return best;
  }
  // A hairpin needs a stem of paired bases with at least three unpaired bases in the loop.
  function hairpin(s) {
    let best = 0;
    for (let i = 0; i < s.length; i++) {
      for (let j = i + 7; j < s.length; j++) {   // >= 3 nt loop between the arms
        let run = 0, a = i, b = j;
        while (a >= 0 && b < s.length && pairs(s[a], s[b]) && (b - a) > 4) { run++; a--; b++; }
        if (run > best) best = run;
      }
    }
    return best;
  }
  function longestRun(s) { let b = 1, r = 1; for (let i = 1; i < s.length; i++) { r = s[i] === s[i - 1] ? r + 1 : 1; if (r > b) b = r; } return b; }

  function evaluate(seq, o) {
    const tm = DNA.tm(seq, o.na, o.conc), gc = DNA.gcPct(seq);
    const last = seq[seq.length - 1];
    const clamp = last === 'G' || last === 'C';
    const gcLast5 = (seq.slice(-5).match(/[GC]/g) || []).length;
    const hp = hairpin(seq), self = endDimer(seq, seq), run = longestRun(seq);
    let pen = 0;
    pen += Math.abs(tm - o.optTm) * 3;
    if (tm < o.minTm || tm > o.maxTm) pen += 18;
    if (gc < o.minGC || gc > o.maxGC) pen += 14;
    pen += Math.abs(gc - 50) * 0.25;
    pen += Math.abs(seq.length - o.optLen) * 0.7;
    if (o.gcClamp && !clamp) pen += 9;
    if (gcLast5 >= 4) pen += 6;            // too stable a 3' end causes mispriming
    if (hp >= 4) pen += (hp - 3) * 7;
    if (self >= 3) pen += (self - 2) * 7;
    if (run >= 4) pen += (run - 3) * 5;
    return { seq, tm, gc, len: seq.length, clamp, hairpin: hp, selfDimer: self, run, penalty: pen };
  }

  function warnings(p, o) {
    const w = [];
    if (p.tm < o.minTm) w.push(`Tm ${p.tm.toFixed(1)} °C is below ${o.minTm} °C`);
    if (p.tm > o.maxTm) w.push(`Tm ${p.tm.toFixed(1)} °C is above ${o.maxTm} °C`);
    if (p.gc < o.minGC || p.gc > o.maxGC) w.push(`GC ${p.gc.toFixed(0)}% is outside ${o.minGC}–${o.maxGC}%`);
    if (!p.clamp) w.push('does not end in G or C');
    if (p.hairpin >= 4) w.push(`hairpin stem of ${p.hairpin} bp`);
    if (p.selfDimer >= 3) w.push(`3′ self-dimer of ${p.selfDimer} bp`);
    if (p.run >= 4) w.push(`run of ${p.run} identical bases`);
    if (p.sites !== undefined && p.sites > 1) w.push(`3′ end matches the template in ${p.sites} places`);
    return w;
  }

  // Candidate primers whose outer (5') end sits at or outside `anchor`.
  // dir = 1 builds a forward primer starting at anchor; dir = -1 a reverse primer ending at anchor.
  function candidates(doc, anchor, dir, o) {
    const seq = doc.seq, n = seq.length, out = [];
    for (let shift = 0; shift <= o.maxShift; shift++) {
      for (let L = o.minLen; L <= o.maxLen; L++) {
        let s;
        if (dir === 1) {
          const p = ((anchor - shift) % n + n) % n;
          if (!doc.circular && (anchor - shift < 0 || p + L > n)) continue;
          s = DNA.sub(seq, p, (p + L) % n || n);
          if (s.length !== L) continue;
          out.push({ ...evaluate(s, o), start: p, end: (p + L) % n, strand: 1, shift });
        } else {
          const q = ((anchor + shift) % n + n) % n || (anchor + shift === n ? n : 0);
          const e = anchor + shift;
          if (!doc.circular && (e > n || e - L < 0)) continue;
          const st = ((e - L) % n + n) % n;
          s = DNA.rc(DNA.sub(seq, st, e % n || n));
          if (s.length !== L) continue;
          out.push({ ...evaluate(s, o), start: st, end: e % n || n, strand: -1, shift });
        }
      }
    }
    for (const c of out) c.penalty += c.shift * 1.2;   // prefer primers right at the requested boundary
    return out.sort((a, b) => a.penalty - b.penalty);
  }

  function countSites(doc, primerSeq) {
    const probe = primerSeq.slice(-15);
    return DNA.searchDNA(doc.seq, doc.circular, probe).length;
  }

  // Design a pair that amplifies [start, end) of the document.
  function designPair(doc, start, end, opts = {}) {
    const o = { ...DEF, ...opts };
    const n = doc.seq.length;
    const fwd = candidates(doc, start, 1, o).slice(0, 60);
    const rev = candidates(doc, end, -1, o).slice(0, 60);
    if (!fwd.length || !rev.length) throw new Error('No primer fits in this region. Try a longer target or a wider length range.');
    const tail5 = DNA.clean(o.fwdTail || ''), tail3 = DNA.clean(o.revTail || '');
    const out = [];
    for (const f of fwd.slice(0, 25)) for (const r of rev.slice(0, 25)) {
      const dTm = Math.abs(f.tm - r.tm);
      if (dTm > o.maxTmDiff + 3) continue;
      const cross = endDimer(f.seq, r.seq);
      let pen = f.penalty + r.penalty + dTm * 5 + (cross >= 3 ? (cross - 2) * 7 : 0);
      if (dTm > o.maxTmDiff) pen += 12;
      const amp = DNA.rlen({ start: f.start, end: r.end }, n);
      out.push({ fwd: f, rev: r, dTm, cross, penalty: pen, ampliconStart: f.start, ampliconEnd: r.end % n, amplicon: amp });
    }
    out.sort((a, b) => a.penalty - b.penalty);
    const top = out.slice(0, opts.limit || 5);
    for (const p of top) {
      p.fwd.sites = countSites(doc, p.fwd.seq); p.rev.sites = countSites(doc, p.rev.seq);
      p.fwd.warnings = warnings(p.fwd, o); p.rev.warnings = warnings(p.rev, o);
      p.fwd.full = tail5 + p.fwd.seq; p.rev.full = tail3 + p.rev.seq;
      p.pairWarnings = [];
      if (p.dTm > o.maxTmDiff) p.pairWarnings.push(`Tm difference of ${p.dTm.toFixed(1)} °C`);
      if (p.cross >= 3) p.pairWarnings.push(`3′ cross-dimer of ${p.cross} bp`);
    }
    if (!top.length) throw new Error('No primer pair met the Tm and dimer limits. Widen the Tm range.');
    return top;
  }

  // Primers spaced along a sequence so that overlapping reads cover all of it.
  function designWalk(doc, opts = {}) {
    const o = { ...DEF, ...opts };
    const n = doc.seq.length;
    const spacing = opts.spacing || 500;        // how much usable sequence one read gives
    const lead = opts.lead || 60;               // bases lost before a read becomes readable
    const both = opts.bothStrands !== false;
    const from = opts.start || 0, span = opts.span || n;
    const out = [];
    for (let off = 0; off < span; off += spacing) {
      const target = (from + off - lead + n) % n;
      if (!doc.circular && target + o.maxLen > n) break;
      const c = candidates(doc, target, 1, { ...o, maxShift: Math.min(60, o.maxShift * 4) })[0];
      if (!c) continue;
      c.sites = countSites(doc, c.seq); c.warnings = warnings(c, o);
      c.name = `${doc.name}-F${out.filter(x => x.strand === 1).length + 1}`;
      c.reads = { start: (c.end + lead) % n, end: (c.end + lead + spacing) % n };
      out.push(c);
    }
    if (both) {
      for (let off = 0; off < span; off += spacing) {
        const target = (from + off + spacing + lead) % n;
        if (!doc.circular && target > n) break;
        const c = candidates(doc, target, -1, { ...o, maxShift: Math.min(60, o.maxShift * 4) })[0];
        if (!c) continue;
        c.sites = countSites(doc, c.seq); c.warnings = warnings(c, o);
        c.name = `${doc.name}-R${out.filter(x => x.strand === -1).length + 1}`;
        out.push(c);
      }
    }
    return out;
  }

  // Primers that put a new sequence (a tag, a mutation) in place of [start, end).
  function designMutagenesis(doc, start, end, insert, opts = {}) {
    const o = { ...DEF, ...opts, optTm: opts.optTm || 62 };
    const n = doc.seq.length;
    const ins = DNA.clean(insert || '');
    const armLen = Math.max(12, Math.min(25, Math.round((o.optLen * 2 - ins.length) / 2)));
    const left = DNA.sub(doc.seq, ((start - armLen) % n + n) % n, start);
    const right = DNA.sub(doc.seq, end % n, (end + armLen) % n || armLen);
    const fwd = left + ins + right;
    const rev = DNA.rc(fwd);
    const anneal = evaluate(left.slice(-12) + ins + right.slice(0, 12), o);
    return {
      fwd: { ...evaluate(fwd, o), seq: fwd, name: 'Mutagenesis-F', strand: 1 },
      rev: { ...evaluate(rev, o), seq: rev, name: 'Mutagenesis-R', strand: -1 },
      insert: ins, removed: DNA.rlen({ start, end }, n), armLen, anneal
    };
  }

  root.PRIMER = { DEF, designPair, designWalk, designMutagenesis, evaluate, candidates, hairpin, maxDimer, endDimer, warnings, countSites };
})(typeof window !== 'undefined' ? window : globalThis);
