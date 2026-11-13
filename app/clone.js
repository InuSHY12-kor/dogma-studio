/* GeneCanvas — cloning simulation: restriction ligation, Gibson, Gateway. Exposes window.CLONE */
(function (root) {
  'use strict';
  const DNA = root.DNA;

  // ---------------------------------------------------------------- fragments & ends
  // An end is described by {type, x}: type '5' (5' overhang), '3' (3' overhang) or 'b' (blunt),
  // x = the single-stranded sequence read on the top strand. Two ends can be ligated when
  // their type and x are identical, because both were produced from the same duplex sequence.
  function endOf(seq, c1, c2) {
    const n = seq.length;
    if (c1 === c2) return { type: 'b', x: '' };
    const a = Math.min(c1, c2), b = Math.max(c1, c2);
    const x = b - a > n / 2 ? DNA.sub(seq, b, a) : seq.slice(a, b); // normal case: no wrap
    return { type: c2 > c1 ? '5' : '3', x };
  }
  function sameEnd(a, b) { return a && b && a.type === b.type && a.x === b.x; }
  function endLabel(e) { return e.type === 'b' ? 'blunt' : `${e.type}′-${e.x}`; }

  // All cuts made by the chosen enzymes, sorted along the top strand.
  function cutList(doc, sites, names) {
    const out = [];
    for (const nm of names) for (const s of (sites[nm] || [])) {
      if (s.cut == null || s.cut2 == null) continue;
      out.push({ c1: s.cut, c2: s.cut2, enzyme: nm });
    }
    out.sort((a, b) => a.c1 - b.c1);
    // drop duplicate cut positions (two enzymes cutting at the same place)
    return out.filter((c, i) => i === 0 || c.c1 !== out[i - 1].c1);
  }

  // Cut a document into fragments. Each fragment carries the top strand plus both end signatures.
  function cut(doc, sites, names) {
    const seq = doc.seq, n = seq.length;
    const cuts = cutList(doc, sites, names);
    const frags = [];
    if (!cuts.length) return frags;
    if (doc.circular) {
      for (let i = 0; i < cuts.length; i++) {
        const a = cuts[i], b = cuts[(i + 1) % cuts.length];
        frags.push(mkFrag(doc, a.c1, b.c1, endOf(seq, a.c1, a.c2), endOf(seq, b.c1, b.c2), a.enzyme, b.enzyme));
      }
    } else {
      const open = { type: 'o', x: '' }; // an original end of a linear molecule: cannot be ligated
      for (let i = 0; i <= cuts.length; i++) {
        const a = i === 0 ? null : cuts[i - 1], b = i === cuts.length ? null : cuts[i];
        const s = a ? a.c1 : 0, e = b ? b.c1 : n;
        if (e <= s && !(i === 0 && cuts.length === 0)) continue;
        frags.push(mkFrag(doc, s, e, a ? endOf(seq, a.c1, a.c2) : open, b ? endOf(seq, b.c1, b.c2) : open, a ? a.enzyme : 'end', b ? b.enzyme : 'end'));
      }
    }
    return frags;
  }
  function mkFrag(doc, start, end, left, right, lenz, renz) {
    const n = doc.seq.length;
    const size = end > start ? end - start : (end === start ? n : end + n - start);
    return { doc, start, end, size, left, right, lenz, renz, seq: DNA.sub(doc.seq, start, end) || doc.seq, label: `${lenz} – ${renz}` };
  }

  // ---------------------------------------------------------------- assembling products
  // Parts are stretches of source sequence; features are carried across automatically.
  function Assembly(name) { this.name = name; this.parts = []; this.seq = ''; this.notes = []; }
  Assembly.prototype.add = function (doc, start, end, rev) {
    let s = DNA.sub(doc.seq, start, end);
    if (start === end) s = doc.seq; // whole molecule
    if (rev) s = DNA.rc(s);
    this.parts.push({ doc, start, end, rev, off: this.seq.length, len: s.length });
    this.seq += s;
    return this;
  };
  Assembly.prototype.addSeq = function (s) { this.parts.push({ raw: true, off: this.seq.length, len: s.length }); this.seq += s; return this; };
  Assembly.prototype.features = function () {
    const out = [];
    for (const p of this.parts) {
      if (p.raw) continue;
      const d = p.doc, n = d.seq.length;
      for (const f of d.features) {
        const rs = ((f.start - p.start) % n + n) % n, L = DNA.rlen(f, n);
        if (rs + L > p.len) continue; // not fully inside this part
        const start = p.rev ? p.off + (p.len - rs - L) : p.off + rs;
        out.push({ name: f.name, type: f.type, start, end: start + L, strand: p.rev ? -f.strand || 0 : f.strand, color: f.color, notes: f.notes });
      }
    }
    return out;
  };
  Assembly.prototype.primers = function () {
    const seen = new Set(), out = [];
    for (const p of this.parts) if (!p.raw) for (const pr of p.doc.primers) if (!seen.has(pr.seq)) { seen.add(pr.seq); out.push({ name: pr.name, seq: pr.seq }); }
    return out;
  };
  Assembly.prototype.doc = function (circular) {
    return { name: this.name, seq: this.seq, circular: circular !== false, features: this.features(), primers: this.primers(), notes: this.notes };
  };

  function fragAssembly(name, frags) {
    const a = new Assembly(name);
    for (const f of frags) {
      if (f.rev) { const r = revFrag(f); a.add(r.doc, r.start, r.end, true); }
      else a.add(f.doc, f.start, f.end, false);
    }
    return a;
  }
  // A reversed fragment: its top strand is the reverse complement, so ends swap and flip.
  function revFrag(f) {
    const n = f.doc.seq.length;
    return { doc: f.doc, start: f.start, end: f.end, rev: true, size: f.size, seq: DNA.rc(f.seq), left: flipEnd(f.right), right: flipEnd(f.left), lenz: f.renz, renz: f.lenz, label: f.label + ' (flipped)' };
  }
  function flipEnd(e) { return e.type === 'o' ? e : { type: e.type, x: DNA.rc(e.x) }; }

  // ---------------------------------------------------------------- restriction ligation
  // Try every order and orientation of the chosen fragments and keep the circles that close.
  function ligate(frags, opts = {}) {
    const name = opts.name || 'Ligation product';
    const first = frags[0];                        // the first fragment anchors the orientation
    const rest = frags.slice(1);
    const out = [];
    const perms = permutations(rest.map((f, i) => i));
    for (const order of perms) {
      const variants = 1 << order.length;
      for (let v = 0; v < variants; v++) {
        const chain = [first];
        for (let i = 0; i < order.length; i++) {
          const f = rest[order[i]];
          chain.push((v >> i) & 1 ? revFrag(f) : f);
        }
        let ok = true;
        for (let i = 0; i < chain.length; i++) {
          const a = chain[i], b = chain[(i + 1) % chain.length];
          if (!sameEnd(a.right, b.left)) { ok = false; break; }
        }
        if (!ok) continue;
        const asm = fragAssembly(name, chain);
        const junctions = chain.map((f, i) => {
          const b = chain[(i + 1) % chain.length];
          return `${f.renz}/${b.lenz} (${endLabel(f.right)})`;
        });
        out.push({ doc: asm.doc(true), chain, junctions, size: asm.seq.length });
      }
    }
    // Deduplicate rotations/orientations of the same product
    const seen = new Set(), uniq = [];
    for (const p of out) {
      const key = canonicalCircle(p.doc.seq);
      if (seen.has(key)) continue;
      seen.add(key); uniq.push(p);
    }
    return uniq;
  }
  function canonicalCircle(s) {
    let best = null;
    for (const t of [s, DNA.rc(s)]) {
      const d = t + t;
      for (let i = 0; i < t.length; i++) { const r = d.slice(i, i + t.length); if (best === null || r < best) best = r; }
    }
    return best;
  }
  function permutations(a) {
    if (a.length <= 1) return [a];
    const out = [];
    for (let i = 0; i < a.length; i++) {
      const rest = a.slice(0, i).concat(a.slice(i + 1));
      for (const p of permutations(rest)) out.push([a[i], ...p]);
    }
    return out;
  }

  // ---------------------------------------------------------------- Gibson / overlap assembly
  function overlapLen(a, b, min, max) {
    const lim = Math.min(max, a.length, b.length);
    for (let L = lim; L >= min; L--) if (a.slice(a.length - L) === b.slice(0, L)) return L;
    return 0;
  }
  // Each input is {doc, seq, name}; linear pieces joined by their shared ends into a circle.
  function gibson(pieces, opts = {}) {
    const min = opts.minOverlap || 15, max = opts.maxOverlap || 120;
    const name = opts.name || 'Assembled construct';
    const first = pieces[0], rest = pieces.slice(1);
    const results = [];
    for (const order of permutations(rest.map((p, i) => i))) {
      const variants = 1 << order.length;
      for (let v = 0; v < variants; v++) {
        const chain = [{ p: first, rev: false }];
        for (let i = 0; i < order.length; i++) chain.push({ p: rest[order[i]], rev: !!((v >> i) & 1) });
        const seqs = chain.map(c => c.rev ? DNA.rc(c.p.seq) : c.p.seq);
        const ov = [];
        let ok = true;
        for (let i = 0; i < seqs.length; i++) {
          const L = overlapLen(seqs[i], seqs[(i + 1) % seqs.length], min, max);
          if (!L) { ok = false; break; }
          ov.push(L);
        }
        if (!ok) continue;
        // Build the circle: keep each overlap once.
        const asm = new Assembly(name);
        for (let i = 0; i < chain.length; i++) {
          const c = chain[i];
          const trimStart = i === 0 ? 0 : ov[i - 1];     // this overlap already came from the previous piece
          const trimEnd = i === chain.length - 1 ? ov[ov.length - 1] : 0; // closing overlap removed once
          const full = seqs[i];
          const piece = full.slice(trimStart, full.length - trimEnd);
          if (!piece.length) { ok = false; break; }
          if (c.p.doc) {
            const d = c.p.doc, n = d.seq.length;
            // map the kept stretch back onto the source document
            const s0 = c.p.start === undefined ? 0 : c.p.start;
            if (!c.rev) asm.add(d, (s0 + trimStart) % n, (s0 + trimStart + piece.length) % n || n, false);
            else {
              const e0 = c.p.end === undefined ? n : c.p.end;
              const st = ((e0 - trimStart - piece.length) % n + n) % n;
              asm.add(d, st, (st + piece.length) % n || n, true);
            }
          } else asm.addSeq(piece);
        }
        if (!ok) continue;
        results.push({
          doc: asm.doc(true), size: asm.seq.length,
          junctions: chain.map((c, i) => ({ from: c.p.name + (c.rev ? ' (flipped)' : ''), to: chain[(i + 1) % chain.length].p.name, len: ov[i], seq: seqs[i].slice(seqs[i].length - ov[i]) }))
        });
      }
    }
    const seen = new Set(), uniq = [];
    for (const r of results) { const k = canonicalCircle(r.doc.seq); if (seen.has(k)) continue; seen.add(k); uniq.push(r); }
    return uniq;
  }

  // ---------------------------------------------------------------- Golden Gate
  // Type IIS enzymes cut outside their own site, so the piece between two inward-facing
  // sites keeps no site of its own. That piece is what each part contributes.
  const TYPE2S = ['BsaI', 'BsmBI', 'BbsI', 'SapI'];

  function hasSite(seq, enzyme) {
    const e = DNA.ENZ_BY_NAME[enzyme];
    return DNA.searchDNA(seq, false, e.site).length > 0;
  }
  // The fragments of one part that carry no recognition site and have two cut ends.
  function ggFragments(doc, enzyme) {
    const sites = DNA.findSites(doc.seq, doc.circular);
    const frags = cut(doc, sites, [enzyme]);
    return frags.filter(f => f.left.type !== 'o' && f.right.type !== 'o' && !hasSite(f.seq, enzyme));
  }

  function goldenGate(parts, enzyme, opts = {}) {
    if (!DNA.ENZ_BY_NAME[enzyme]) throw new Error(`Unknown enzyme ${enzyme}.`);
    const chosen = [], notes = [], undomesticated = [];
    for (const p of parts) {
      const doc = p.doc || p;
      const n = (DNA.findSites(doc.seq, doc.circular)[enzyme] || []).length;
      const all = ggFragments(doc, enzyme);
      if (!all.length) {
        throw new Error(n < 2
          ? `“${doc.name}” has ${n} ${enzyme} site${n === 1 ? '' : 's'}. A Golden Gate part needs two, pointing inward.`
          : `“${doc.name}” has ${enzyme} sites, but none of its fragments is free of them. Check that the sites point inward.`);
      }
      // A part is meant to carry exactly the two flanking sites. Any more sit inside it
      // and will be cut during the reaction, so they have to be removed first.
      if (n > 2) {
        undomesticated.push({ name: doc.name, sites: n });
        notes.push(`“${doc.name}” has ${n} ${enzyme} sites, so ${n - 2} lie inside the part. Remove them with codon optimization before assembling.`);
      }
      const pick = p.fragment != null ? all[p.fragment] : all.slice().sort((a, b) => b.size - a.size)[0];
      if (!pick) throw new Error(`“${doc.name}” has no usable fragment.`);
      if (all.length > 1) notes.push(`“${doc.name}” gave ${all.length} site-free fragments; the ${fmtbp(pick.size)} one was used.`);
      chosen.push(pick);
    }
    // Every junction must use a different overhang, or the order is ambiguous.
    const ends = chosen.map(f => f.right.x);
    const dupes = ends.filter((x, i) => ends.indexOf(x) !== i);
    const products = ligate(chosen, { name: opts.name || 'Golden Gate product' });
    for (const pr of products) {
      pr.leftover = DNA.searchDNA(pr.doc.seq, true, DNA.ENZ_BY_NAME[enzyme].site).length;
      pr.overhangs = pr.chain.map(f => f.right.x);
    }
    return { products, parts: chosen, notes, duplicateOverhangs: [...new Set(dupes)], undomesticated, enzyme };
  }
  function fmtbp(n) { return n.toLocaleString('en-US') + ' bp'; }

  // ---------------------------------------------------------------- Gateway
  // Recombination happens inside a 7 bp core. att1 sites carry TTTGTAC, att2 sites carry the
  // same core in the opposite orientation, so both molecules are cut in the core and swapped.
  const ATT = {
    1: { core: 'TTTGTAC', ctx: 'TTTGTACAAAAAAGCAGGCT' },
    2: { core: 'TTTGTAC', ctx: 'TTTGTACAAGAAAGCTGGGT' }   // written in the attB2 reverse orientation
  };
  const ATT_B = { 1: 'ACAAGTTTGTACAAAAAAGCAGGCT', 2: 'ACCCAGCTTTCTTGTACAAAGTGGT' };
  const ARM = 'ATAATGCCAAC'; // lambda arm, present in attL/attP/attR but not in attB

  // Finds att1 / att2 sites with the position and orientation of their crossover core.
  function findAtt(seq, circular) {
    const out = [];
    for (const fam of [1, 2]) {
      const ctx = ATT[fam].ctx;
      for (const hit of DNA.searchDNA(seq, circular, ctx)) {
        const n = seq.length;
        // the core sits at the start of ctx on the forward strand, at its end on the reverse
        const corePos = hit.strand === 1 ? hit.start : (((hit.end - ATT[fam].core.length) % n) + n) % n;
        const near = DNA.sub(seq, ((corePos - 80) % n + n) % n, (corePos + 80) % n);
        // attB sites are short; attL, attR and attP all carry a lambda arm beside the core
        const hasArm = near.includes(ARM) || near.includes(DNA.rc(ARM));
        const kind = hasArm ? 'att' + fam : 'attB' + fam;
        out.push({ family: fam, pos: corePos, strand: hit.strand, hasArm, kind, label: kind + (hasArm ? ' (λ arm)' : '') });
      }
    }
    // keep one site per family: the clearest pair
    out.sort((a, b) => a.pos - b.pos);
    return out;
  }

  // The segment of `doc` running from its att1 core to its att2 core, and the rest of the molecule.
  function attSplit(doc) {
    const sites = findAtt(doc.seq, doc.circular);
    const s1 = sites.filter(s => s.family === 1), s2 = sites.filter(s => s.family === 2);
    if (!s1.length || !s2.length) return null;
    // choose the pair whose inner span is shortest, which is the cassette
    let best = null;
    for (const a of s1) for (const b of s2) {
      const L = DNA.rlen({ start: a.pos, end: b.pos }, doc.seq.length);
      if (!best || L < best.len) best = { a, b, len: L };
    }
    return { site1: best.a, site2: best.b, inner: { start: best.a.pos, end: best.b.pos }, outer: { start: best.b.pos, end: best.a.pos }, sites };
  }

  // donor contributes the stretch between its att sites; acceptor contributes the rest of itself.
  function gateway(donor, acceptor, opts = {}) {
    const A = attSplit(donor), B = attSplit(acceptor);
    if (!A) throw new Error(`No att1 and att2 pair found in “${donor.name}”.`);
    if (!B) throw new Error(`No att1 and att2 pair found in “${acceptor.name}”.`);
    const main = new Assembly(opts.name || 'Recombination product');
    main.add(donor, A.inner.start, A.inner.end, false);     // att1 core → att2 core, carries the insert
    main.add(acceptor, B.outer.start, B.outer.end, false);  // att2 core → att1 core, the backbone
    const by = new Assembly(opts.byName || 'By-product');
    by.add(acceptor, B.inner.start, B.inner.end, false);
    by.add(donor, A.outer.start, A.outer.end, false);
    const product = main.doc(true), byproduct = by.doc(true);
    // Report which sites the junctions actually spell out.
    const junction = fam => {
      const want = ATT_B[fam];
      return DNA.searchDNA(product.seq, true, want).length ? 'attB' + fam : 'att' + fam;
    };
    return {
      product, byproduct,
      donorSites: [A.site1, A.site2], acceptorSites: [B.site1, B.site2],
      junctions: [junction(1), junction(2)],
      insertSize: DNA.rlen(A.inner, donor.seq.length)
    };
  }

  root.CLONE = { cut, ligate, gibson, gateway, findAtt, attSplit, endLabel, sameEnd, Assembly, canonicalCircle, overlapLen, ATT_B, goldenGate, ggFragments, TYPE2S, hasSite };
})(typeof window !== 'undefined' ? window : globalThis);
