/* GeneCanvas — the cloning engine behind the Vector / Fragment / Product window.
 * Exposes window.ASSEMBLY
 *
 * A run is an ordered, circular list of pieces. The first piece is the vector backbone;
 * the rest are the inserts. Each piece says how it is prepared — cut with enzymes,
 * amplified by PCR, or used as it is — and that decides which side of a junction can
 * carry the overlap or the overhang that joins them.
 */
(function (root) {
  'use strict';
  const DNA = root.DNA;

  // ---------------------------------------------------------------- preparing a piece
  // Returns the stretch a piece contributes, with the coordinates it came from so the
  // features of the source document can be carried across.
  function prepare(doc, spec) {
    const n = doc.seq.length;
    if (!doc || !n) throw new Error('No sequence.');
    if (spec.prep === 'direct') {
      if (doc.circular) throw new Error(`“${doc.name}” is circular, so it cannot be used directly. Cut it or amplify it first.`);
      return { seq: doc.seq, from: 0, to: n, doc, label: 'used directly', ends: { left: { type: 'blunt', x: '' }, right: { type: 'blunt', x: '' } } };
    }
    if (spec.prep === 'pcr') {
      const from = spec.from | 0, to = spec.to | 0;
      const seq = DNA.sub(doc.seq, from, to);
      if (!seq.length) throw new Error('The region to amplify is empty.');
      return { seq, from, to, doc, label: `PCR ${from + 1}..${to}`, ends: { left: { type: 'blunt', x: '' }, right: { type: 'blunt', x: '' } } };
    }
    // cut with one or two enzymes
    const cuts = spec.cuts || [];
    if (!cuts.length) throw new Error('Choose at least one enzyme.');
    if (cuts.length === 1) {
      if (!doc.circular) throw new Error(`“${doc.name}” is linear, so one cut gives two pieces. Choose which piece, or cut with two enzymes.`);
      const c = cuts[0];
      return { seq: DNA.sub(doc.seq, c.cut, c.cut) || doc.seq, from: c.cut, to: c.cut, doc, label: `${c.enzyme} (${c.cut})`, ends: { left: endOf(doc.seq, c), right: endOf(doc.seq, c) } };
    }
    const [a, b] = cuts;
    // Two cuts give two arcs. The backbone is the longer one unless asked otherwise, and
    // the arc that is kept decides which cut end is at its start and which at its end.
    const arc = (from, to) => DNA.rlen({ start: from.cut, end: to.cut }, n) || n;
    const aFirst = arc(a, b), bFirst = arc(b, a);
    const takeAtoB = spec.keepInner ? (aFirst <= bFirst) : (aFirst >= bFirst);
    const startCut = takeAtoB ? a : b, endCut = takeAtoB ? b : a;
    const seq = DNA.sub(doc.seq, startCut.cut, endCut.cut) || doc.seq;
    if (!seq.length) throw new Error('The two cuts leave nothing.');
    return {
      seq, from: startCut.cut, to: endCut.cut, doc,
      label: `${startCut.enzyme} (${startCut.cut}) – ${endCut.enzyme} (${endCut.cut})`,
      removed: n - seq.length,
      ends: { left: endOf(doc.seq, startCut), right: endOf(doc.seq, endCut) }
    };
  }
  function endOf(seq, c) {
    if (c.cut === c.cut2) return { type: 'blunt', x: '', enzyme: c.enzyme };
    const lo = Math.min(c.cut, c.cut2), hi = Math.max(c.cut, c.cut2);
    return { type: c.cut2 > c.cut ? '5' : '3', x: seq.slice(lo, hi), enzyme: c.enzyme };
  }

  // Flips a piece so it goes into the product the other way round.
  function flip(piece) {
    return {
      ...piece, seq: DNA.rc(piece.seq), rev: !piece.rev,
      ends: { left: flipEnd(piece.ends.right), right: flipEnd(piece.ends.left) }
    };
  }
  const flipEnd = e => e ? { ...e, x: DNA.rc(e.x || '') } : e;

  // ---------------------------------------------------------------- overlaps for Gibson
  const tmOf = s => (s && s.length >= 2) ? DNA.tm(s) : 0;

  // Chooses the stretch of homology for one junction. `where` says which neighbour it is
  // copied from: the end of the piece on the left, or the start of the piece on the right.
  function chooseOverlap(seq, opts, where) {
    const take = L => where === 'start' ? seq.slice(0, L) : seq.slice(-L);
    const lo = Math.max(10, opts.minOverlap | 0), hi = Math.max(lo, opts.maxOverlap | 0);
    if (opts.fixed) {
      const L = Math.min(opts.fixedLength, seq.length);
      return { len: L, seq: take(L), tm: tmOf(take(L)) };
    }
    let best = null;
    for (let L = lo; L <= Math.min(hi, seq.length); L++) {
      const s = take(L), t = tmOf(s), d = Math.abs(t - opts.overlapTm);
      if (!best || d < best.d) best = { len: L, seq: s, tm: t, d };
    }
    if (!best) { const L = Math.min(lo, seq.length); best = { len: L, seq: take(L), tm: tmOf(take(L)) }; }
    return best;
  }
  // Grows the annealing part of a primer until it reaches the wanted melting temperature.
  function annealPart(seq, targetTm, minLen, maxLen) {
    let best = null;
    for (let L = minLen; L <= Math.min(maxLen, seq.length); L++) {
      const s = seq.slice(0, L), t = DNA.tm(s), d = Math.abs(t - targetTm);
      if (!best || d < best.d) best = { seq: s, tm: t, d, len: L };
    }
    return best || { seq: seq.slice(0, minLen), tm: 0, len: Math.min(minLen, seq.length) };
  }

  // Works out every junction of the circle and, where a piece is made by PCR, the primers
  // that give it the overlapping ends Gibson needs.
  function planGibson(pieces, opts) {
    const o = {
      minOverlap: 15, maxOverlap: 25, fixed: false, fixedLength: 40,
      overlapTm: 50, primerTm: 60, minAnneal: 16, maxAnneal: 36, ...opts
    };
    const N = pieces.length;
    const junctions = [], primers = pieces.map(() => null);
    for (let i = 0; i < N; i++) {
      const left = pieces[i], right = pieces[(i + 1) % N];
      const carrier = right.spec.prep === 'pcr' ? 'right' : (left.spec.prep === 'pcr' ? 'left' : null);
      if (!carrier) {
        // Neither side is amplified, so the two ends have to match already.
        const nat = naturalOverlap(left.seq, right.seq, o.minOverlap, 200);
        junctions.push({
          from: left.name, to: right.name, carrier: null,
          len: nat, seq: nat ? left.seq.slice(-nat) : '', tm: nat ? tmOf(left.seq.slice(-nat)) : 0,
          ok: nat >= o.minOverlap,
          why: nat >= o.minOverlap ? '' : `${left.name} and ${right.name} share no end of at least ${o.minOverlap} bp, and neither is amplified, so nothing can add one.`
        });
        continue;
      }
      // The tail copies the neighbour that is not being amplified, so the homology is
      // taken from the end of the left piece or the start of the right piece accordingly.
      const ov = carrier === 'right'
        ? chooseOverlap(left.seq, o, 'end')
        : chooseOverlap(right.seq, o, 'start');
      junctions.push({ from: left.name, to: right.name, carrier, len: ov.len, seq: ov.seq, tm: ov.tm, ok: true, why: '' });
    }
    // Build the two primers for each amplified piece from the junctions on either side.
    for (let i = 0; i < N; i++) {
      const p = pieces[i];
      if (p.spec.prep !== 'pcr') continue;
      const before = junctions[(i - 1 + N) % N];   // junction entering this piece
      const after = junctions[i];                  // junction leaving this piece
      const fwdAnneal = annealPart(p.seq, o.primerTm, o.minAnneal, o.maxAnneal);
      const revAnneal = annealPart(DNA.rc(p.seq), o.primerTm, o.minAnneal, o.maxAnneal);
      const fwdTail = before.carrier === 'right' ? before.seq : '';
      const revTail = after.carrier === 'left' ? DNA.rc(after.seq) : '';
      primers[i] = {
        fwd: { tail: fwdTail, anneal: fwdAnneal.seq, seq: fwdTail + fwdAnneal.seq, tm: fwdAnneal.tm, name: `${p.name}.FOR` },
        rev: { tail: revTail, anneal: revAnneal.seq, seq: revTail + revAnneal.seq, tm: revAnneal.tm, name: `${p.name}.REV` }
      };
    }
    return { junctions, primers, opts: o };
  }
  function naturalOverlap(a, b, min, max) {
    for (let L = Math.min(max, a.length, b.length); L >= min; L--) if (a.slice(-L) === b.slice(0, L)) return L;
    return 0;
  }

  // ---------------------------------------------------------------- restriction junctions
  function planRestriction(pieces) {
    const N = pieces.length, junctions = [];
    for (let i = 0; i < N; i++) {
      const left = pieces[i], right = pieces[(i + 1) % N];
      const a = left.ends.right, b = right.ends.left;
      const same = a && b && a.type === b.type && (a.x || '') === (b.x || '');
      junctions.push({
        from: left.name, to: right.name, carrier: null,
        len: (a && a.x ? a.x.length : 0), seq: a ? a.x : '', tm: 0, ok: !!same,
        overhang: a ? (a.type === 'blunt' ? 'blunt' : `${a.type}′-${a.x}`) : '',
        why: same ? '' : `${left.name} ends ${describeEnd(a)} but ${right.name} begins ${describeEnd(b)}, so they cannot be ligated.`
      });
    }
    return { junctions, primers: pieces.map(() => null), opts: {} };
  }
  const describeEnd = e => !e ? 'with an open end' : (e.type === 'blunt' ? 'blunt' : `with a ${e.type}′ overhang ${e.x}`);

  // ---------------------------------------------------------------- the product
  function build(pieces, plan, name) {
    let seq = '', features = [], primers = [];
    const bounds = [];
    for (const p of pieces) {
      const off = seq.length;
      bounds.push({ name: p.name, start: off, end: off + p.seq.length, piece: p });
      seq += p.seq;
      const n = p.doc ? p.doc.seq.length : 0;
      if (p.doc) {
        for (const f of p.doc.features) {
          const rs = ((f.start - p.from) % n + n) % n, L = DNA.rlen(f, n);
          if (rs + L > p.seq.length) continue;
          const start = p.rev ? off + (p.seq.length - rs - L) : off + rs;
          features.push({ ...f, id: undefined, start, end: start + L, strand: p.rev ? (-f.strand || 0) : f.strand });
        }
      }
    }
    for (let i = 0; i < pieces.length; i++) {
      const pr = plan.primers[i];
      if (pr) primers.push({ name: pr.fwd.name, seq: pr.fwd.seq }, { name: pr.rev.name, seq: pr.rev.seq });
    }
    return { name: name || 'Assembled', seq, circular: true, features, primers, bounds };
  }

  // ---------------------------------------------------------------- reading frames
  const STOPS = { TAA: 1, TAG: 1, TGA: 1 };
  // Reports coding features that the assembly has run together out of frame. Only pairs
  // that really would be read as one open reading frame are counted: same strand, close
  // together, and with no stop codon between them in the upstream frame.
  function frameIssues(product, bounds, maxGap) {
    const gapLimit = maxGap === undefined ? 60 : maxGap;
    const seq = product.seq, n = seq.length;
    const cds = product.features
      .filter(f => f.type === 'CDS' && f.strand !== 0 && DNA.rlen(f, n) >= 30)
      .sort((a, b) => a.start - b.start);
    const junctionAt = bounds.slice(1).map(b => b.start).concat([0]);
    const out = [];
    for (let i = 0; i < cds.length; i++) {
      for (let j = 0; j < cds.length; j++) {
        if (i === j) continue;
        const a = cds[i], b = cds[j];
        if (a.strand !== b.strand) continue;
        const gap = b.start - a.end;
        if (gap < 0 || gap > gapLimit) continue;
        // only when the two came from different pieces, which is what the assembly changed
        const crosses = junctionAt.some(x => x > a.start && x <= b.start);
        if (!crosses) continue;
        if (a.strand === 1) {
          if (stopBetween(seq, a.end - ((a.end - a.start) % 3), b.start)) continue;
          if ((b.start - a.start) % 3 === 0) continue;
          out.push({ a, b, shift: ((b.start - a.start) % 3 + 3) % 3, gap });
        } else {
          if ((b.end - a.end) % 3 === 0) continue;
          out.push({ a, b, shift: ((b.end - a.end) % 3 + 3) % 3, gap });
        }
      }
    }
    return out;
  }
  function stopBetween(seq, from, to) {
    for (let i = Math.max(0, from); i + 3 <= to; i += 3) if (STOPS[seq.substr(i, 3)]) return true;
    return false;
  }

  // Tries every orientation of the inserts and reports one that makes all the joins fit.
  // Only the inserts are turned; the vector sets the direction of the whole circle.
  function suggestFlips(pieces, planner) {
    const n = pieces.length - 1;
    if (n < 1 || n > 6) return null;
    for (let mask = 1; mask < (1 << n); mask++) {
      const trial = pieces.map((p, i) => {
        if (i === 0 || !((mask >> (i - 1)) & 1)) return p;
        const f = flip(p);
        return Object.assign(f, { name: p.name, spec: p.spec, role: p.role, index: p.index, doc: p.doc, from: p.from, to: p.to, label: p.label });
      });
      const plan = planner(trial);
      if (plan.junctions.every(j => j.ok)) {
        const which = [];
        for (let i = 0; i < n; i++) if ((mask >> i) & 1) which.push(pieces[i + 1].name);
        return { mask, flip: which };
      }
    }
    return null;
  }

  root.ASSEMBLY = { prepare, flip, planGibson, planRestriction, build, frameIssues, suggestFlips, chooseOverlap, annealPart, naturalOverlap, describeEnd, endOf };
})(typeof window !== 'undefined' ? window : globalThis);
