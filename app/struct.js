/* Dogma Studio — structure and thermodynamics. Exposes window.STRUCT
 *
 * Three separate things live here:
 *   1. Oligo structures (hairpin, self-dimer, cross-dimer) from nearest-neighbour
 *      free energies. These are derived from the same SantaLucia 1998 parameters the
 *      melting temperature uses, so the numbers are on a published footing.
 *   2. An RNA secondary structure fold. The recursion is a normal Zuker minimum-free-energy
 *      search, but the energy model here is a simplified one, not the full Turner set, so
 *      its free energies are approximate and the UI says so.
 *   3. Protein properties. Mass, extinction coefficient, hydropathy and the aliphatic
 *      index are exact formulas; the isoelectric point depends on which pKa set is used
 *      and the set in use is reported with the answer.
 */
(function (root) {
  'use strict';
  const DNA = root.DNA;

  const R = 1.9872e-3;        // kcal/(mol·K)
  const T37 = 310.15;

  // ---------------------------------------------------------------- DNA nearest neighbours
  // ΔH (kcal/mol) and ΔS (cal/(mol·K)) for each stack, read on the top strand 5'→3'.
  const NN = {
    AA: [-7.9, -22.2], TT: [-7.9, -22.2], AT: [-7.2, -20.4], TA: [-7.2, -21.3],
    CA: [-8.5, -22.7], TG: [-8.5, -22.7], GT: [-8.4, -22.4], AC: [-8.4, -22.4],
    CT: [-7.8, -21.0], AG: [-7.8, -21.0], GA: [-8.2, -22.2], TC: [-8.2, -22.2],
    CG: [-10.6, -27.2], GC: [-9.8, -24.4], GG: [-8.0, -19.9], CC: [-8.0, -19.9]
  };
  const dgOf = (h, s, T = T37) => h - T * s / 1000;
  // The tables are built at a temperature, so a structure can be scored at the
  // temperature it will actually be held at. 37 C is the default everywhere.
  const tableCache = new Map();
  function dgTable(tempC) {
    const T = 273.15 + (tempC == null ? 37 : tempC);
    const key = T.toFixed(2);
    if (tableCache.has(key)) return tableCache.get(key);
    const dg = {};
    for (const [k, [h, s2]] of Object.entries(NN)) dg[k] = dgOf(h, s2, T);
    const t = { T, dg, initGC: dgOf(0.1, -2.8, T), initAT: dgOf(2.3, 4.1, T) };
    tableCache.set(key, t);
    return t;
  }
  const DG = dgTable(37).dg;
  // Helix initiation, split between the two ends by which base closes them.
  const INIT_GC = dgTable(37).initGC, INIT_AT = dgTable(37).initAT;
  const initFor = b => (b === 'G' || b === 'C') ? INIT_GC : INIT_AT;

  const COMP = { A: 'T', T: 'A', G: 'C', C: 'G' };
  const pairs = (a, b) => COMP[a] === b;

  // ---------------------------------------------------------------- oligo duplexes
  // Slides two oligos past one another and scores every register. `b` is read 5'→3'
  // and laid antiparallel under `a`.
  function duplexes(a, b, opts = {}) {
    a = DNA.clean(a); b = DNA.clean(b);
    if (!a || !b) return [];
    const T = dgTable(opts.tempC);
    const DGx = T.dg, initAt = ch => (ch === 'G' || ch === 'C') ? T.initGC : T.initAT;
    const br = [...b].reverse().join('');       // b written 3'→5', which is how it faces a
    const out = [];
    for (let k = -(a.length - 1); k <= br.length - 1; k++) {
      const pat = [];
      let any = false;
      for (let i = 0; i < a.length; i++) {
        const j = i - k;
        const ok = j >= 0 && j < br.length && pairs(a[i], br[j]);
        pat.push(ok);
        if (ok) any = true;
      }
      if (!any) continue;
      // free energy from the stacks inside each run of consecutive pairs
      let dg = 0, runs = 0, best = 0, cur = 0, total = 0;
      for (let i = 0; i < a.length; i++) {
        if (!pat[i]) { cur = 0; continue; }
        total++;
        cur++; if (cur > best) best = cur;
        if (cur === 1) { runs++; dg += initAt(a[i]); }
        else dg += DGx[a[i - 1] + a[i]];
        // closing penalty when the run ends
        if (i === a.length - 1 || !pat[i + 1]) dg += initAt(a[i]);
      }
      if (runs > 1) dg += (runs - 1) * 1.6;     // crossing an internal mismatch costs energy
      // how many of the last five bases of each primer are paired, which is what primes
      const end3a = pat.slice(-5).filter(Boolean).length;
      let end3b = 0;
      for (let i = 0; i < a.length; i++) { const j = i - k; if (pat[i] && j >= br.length - 5) end3b++; }
      out.push({ offset: k, pattern: pat, dg, pairs: total, longest: best, runs, end3a, end3b });
    }
    out.sort((x, y) => x.dg - y.dg);
    const limit = opts.limit || 4;
    // keep distinct registers only
    const seen = new Set(), keep = [];
    for (const d of out) {
      const key = d.pattern.join('');
      if (seen.has(key)) continue;
      seen.add(key); keep.push(d);
      if (keep.length >= limit) break;
    }
    return keep;
  }

  // Renders a register as the three lines a lab notebook would draw.
  function duplexDiagram(a, b, d) {
    a = DNA.clean(a); b = DNA.clean(b);
    const br = [...b].reverse().join('');
    const k = d.offset;
    const padA = Math.max(0, -k), padB = Math.max(0, k);
    const top = ' '.repeat(padA) + a;
    const bot = ' '.repeat(padB) + br;
    let bars = '';
    for (let i = 0; i < Math.max(top.length, bot.length); i++) {
      const ai = i - padA, bi = i - padB;
      bars += (ai >= 0 && ai < a.length && d.pattern[ai] && bi >= 0 && bi < br.length) ? '|' : ' ';
    }
    const w = Math.max(top.length, bot.length);
    return {
      top: "5'-" + top.padEnd(w) + "-3'",
      bars: '   ' + bars.padEnd(w) + '   ',
      bot: "3'-" + bot.padEnd(w) + "-5'"
    };
  }

  // ---------------------------------------------------------------- hairpins
  // A stem of paired bases closing a loop of at least three unpaired bases.
  function hairpins(s, opts = {}) {
    s = DNA.clean(s);
    const TT = dgTable(opts.tempC), DGx = TT.dg;
    const n = s.length, minLoop = opts.minLoop || 3, minStem = opts.minStem || 3;
    const found = [];
    for (let i = 0; i < n; i++) {
      for (let j = i + minLoop + 1; j < n; j++) {
        if (!pairs(s[i], s[j])) continue;
        // grow the stem inward as far as the bases keep pairing
        let stem = 0, x = i, y = j;
        while (x < y && pairs(s[x], s[y]) && (y - x - 1) >= minLoop) { stem++; x++; y--; }
        if (stem < minStem) continue;
        // only report a stem that cannot be extended outward, so each one appears once
        if (i > 0 && j < n - 1 && pairs(s[i - 1], s[j + 1])) continue;
        const loop = y - x + 1;
        let dg = 0;
        for (let q = 1; q < stem; q++) dg += DGx[s[i + q - 1] + s[i + q]];
        dg += loopPenalty(loop);
        found.push({ start: i, end: j, stem, loop, loopStart: x, dg });
      }
    }
    found.sort((p, q) => p.dg - q.dg);
    const keep = [], used = [];
    for (const h of found) {
      if (used.some(u => !(h.end < u.start || h.start > u.end))) continue;
      used.push(h); keep.push(h);
      if (keep.length >= (opts.limit || 3)) break;
    }
    return keep;
  }
  // Loop cost grows with size; this is the usual logarithmic form, not a measured table.
  function loopPenalty(n) {
    const base = { 3: 3.5, 4: 3.5, 5: 3.3, 6: 3.4, 7: 3.5, 8: 3.6, 9: 3.7 };
    if (base[n] !== undefined) return base[n];
    if (n < 3) return 9;
    return 3.7 + 1.75 * R * T37 * Math.log(n / 9);
  }
  function hairpinDiagram(s, h) {
    s = DNA.clean(s);
    const left = s.slice(h.start, h.start + h.stem);
    const right = s.slice(h.end - h.stem + 1, h.end + 1);
    const loop = s.slice(h.loopStart, h.loopStart + h.loop);
    return {
      pre: s.slice(0, h.start), post: s.slice(h.end + 1),
      left, right: [...right].reverse().join(''), loop,
      bars: '|'.repeat(h.stem)
    };
  }

  // Everything a primer pair can do to itself and to its partner.
  function primerStructures(fwd, rev, opts = {}) {
    const o = { limit: 3, tempC: opts.tempC, minPairs: opts.minPairs };
    const r = { fwd: { seq: fwd }, tempC: opts.tempC == null ? 37 : opts.tempC };
    r.fwd.hairpins = hairpins(fwd, { tempC: opts.tempC, minStem: opts.minStem, minLoop: opts.minLoop });
    r.fwd.self = duplexes(fwd, fwd, o);
    if (rev) {
      r.rev = { seq: rev, hairpins: hairpins(rev, { tempC: opts.tempC, minStem: opts.minStem, minLoop: opts.minLoop }), self: duplexes(rev, rev, o) };
      r.cross = duplexes(fwd, rev, o);
    }
    return r;
  }

  // ---------------------------------------------------------------- RNA fold
  // Zuker-style minimum free energy search. The recursion is standard; the energy model
  // is deliberately simple and its numbers are approximate.
  const RNA_PAIR = { AU: 1, UA: 1, GC: 1, CG: 1, GU: 1, UG: 1 };
  // Stacking energies, kcal/mol, keyed by the closing pair then the pair stacked on it.
  const RNA_STACK = {
    CG: { CG: -3.3, GC: -2.4, GU: -2.1, UG: -1.4, AU: -2.1, UA: -2.1 },
    GC: { CG: -3.4, GC: -3.3, GU: -2.5, UG: -1.5, AU: -2.2, UA: -2.4 },
    GU: { CG: -2.5, GC: -2.1, GU: -1.3, UG: -0.5, AU: -1.4, UA: -1.3 },
    UG: { CG: -2.1, GC: -1.4, GU: 0.3, UG: -1.3, AU: -1.0, UA: -1.0 },
    AU: { CG: -2.4, GC: -2.1, GU: -1.3, UG: -1.0, AU: -0.9, UA: -1.1 },
    UA: { CG: -2.2, GC: -2.1, GU: -1.4, UG: -0.6, AU: -1.3, UA: -0.9 }
  };
  const HAIRPIN_INIT = { 3: 5.4, 4: 5.6, 5: 5.7, 6: 5.4, 7: 6.0, 8: 5.5, 9: 6.4 };
  const BULGE_INIT = { 1: 3.8, 2: 2.8, 3: 3.2, 4: 3.6, 5: 4.0, 6: 4.4 };
  const INTERNAL_INIT = { 2: 1.5, 3: 1.8, 4: 1.7, 5: 1.8, 6: 2.0, 7: 2.2, 8: 2.3 };
  const ML_A = 3.4, ML_B = 0.4, ML_C = 0.0;   // multibranch: closing, per helix, per unpaired
  const AU_END = 0.45;                        // penalty for a helix ending on A·U or G·U

  const logExtrap = (tab, n, k) => {
    const keys = Object.keys(tab).map(Number);
    const max = Math.max(...keys);
    if (tab[n] !== undefined) return tab[n];
    return tab[max] + k * R * T37 * Math.log(n / max);
  };
  const hairpinE = n => n < 3 ? Infinity : logExtrap(HAIRPIN_INIT, n, 1.75);
  const bulgeE = n => logExtrap(BULGE_INIT, n, 1.75);
  const internalE = n => logExtrap(INTERNAL_INIT, n, 1.08);
  const auPenalty = p => (p === 'AU' || p === 'UA' || p === 'GU' || p === 'UG') ? AU_END : 0;

  function toRNA(s) { return DNA.clean(s).replace(/T/g, 'U'); }

  function rnaFold(seq, opts = {}) {
    const s = toRNA(seq);
    const n = s.length;
    const MAXLOOP = 30, MINLOOP = 3;
    if (n < 5) return { seq: s, structure: '.'.repeat(n), dg: 0, pairs: [], approximate: true };
    const cap = opts.maxLength || 600;
    if (n > cap) throw new Error(`This sequence is ${n} nt. Folding here is limited to ${cap} nt, because the search grows with the cube of the length. Select a shorter region, or fold the whole thing with ViennaRNA RNAfold.`);

    const INF = 1e9;
    const V = new Float64Array(n * n).fill(INF);   // i pairs with j
    const W = new Float64Array(n * n).fill(0);     // best for the subsequence i..j
    const WM = new Float64Array(n * n).fill(INF);  // inside a multibranch loop
    const at = (i, j) => i * n + j;
    const pairOf = (i, j) => s[i] + s[j];
    const canPair = (i, j) => j - i - 1 >= MINLOOP && RNA_PAIR[pairOf(i, j)];
    // For each i, the j it is already known to pair with. Filled as the table grows, and
    // always ascending, so the inner-loop search can skip straight to the useful range.
    const vPartners = Array.from({ length: n }, () => []);

    for (let d = 1; d < n; d++) {
      for (let i = 0; i + d < n; i++) {
        const j = i + d;
        // ---- V: i pairs j
        if (canPair(i, j)) {
          const p = pairOf(i, j);
          let best = hairpinE(j - i - 1) + auPenalty(p);
          // stack, bulge or internal loop onto an inner pair
          for (let k = i + 1; k <= Math.min(i + MAXLOOP + 1, j - MINLOOP - 2); k++) {
            const list = vPartners[k];
            if (!list.length) continue;
            const lo = Math.max(k + MINLOOP + 1, j - 1 - (MAXLOOP - (k - i - 1)));
            let a0 = 0, a1 = list.length;
            while (a0 < a1) { const m = (a0 + a1) >> 1; if (list[m] < lo) a0 = m + 1; else a1 = m; }
            for (let t = a0; t < list.length; t++) {
              const l = list[t];
              if (l >= j) break;
              const u1 = k - i - 1, u2 = j - l - 1;
              if (u1 + u2 > MAXLOOP) continue;
              const q = pairOf(k, l);
              let e;
              if (u1 === 0 && u2 === 0) e = (RNA_STACK[p] && RNA_STACK[p][q] !== undefined) ? RNA_STACK[p][q] : -0.5;
              else if (u1 === 0 || u2 === 0) e = bulgeE(u1 + u2) + auPenalty(p) + auPenalty(q);
              else e = internalE(u1 + u2) + 0.6 * Math.abs(u1 - u2) + auPenalty(p) + auPenalty(q);
              const cand = e + V[at(k, l)];
              if (cand < best) best = cand;
            }
          }
          // close a multibranch loop
          const rowI1 = (i + 1) * n, jm1 = j - 1, mlBase = ML_A + ML_B + auPenalty(p);
          for (let k = i + 2; k < jm1; k++) {
            const a = WM[rowI1 + k];
            if (a >= INF) continue;
            const b = WM[(k + 1) * n + jm1];
            if (b >= INF) continue;
            const cand = a + b + mlBase;
            if (cand < best) best = cand;
          }
          V[at(i, j)] = best;
          if (best < INF) vPartners[i].push(j);
        }
        // ---- WM: a piece of a multibranch loop
        const rowI = i * n;
        let wm = INF;
        const vij = V[rowI + j];
        if (vij < INF) wm = vij + ML_B + auPenalty(pairOf(i, j));
        const up = WM[rowI + n + j];
        if (up < INF && up + ML_C < wm) wm = up + ML_C;
        const left = WM[rowI + j - 1];
        if (left < INF && left + ML_C < wm) wm = left + ML_C;
        for (let k = i + 1; k < j; k++) {
          const a = WM[rowI + k];
          if (a >= INF) continue;
          const b = WM[k * n + n + j];
          if (b >= INF) continue;
          const c = a + b;
          if (c < wm) wm = c;
        }
        WM[rowI + j] = wm;
        // ---- W: best overall for i..j
        let w = W[at(i, j - 1)];                       // j unpaired
        if (V[at(i, j)] < INF) w = Math.min(w, V[at(i, j)] + auPenalty(pairOf(i, j)));
        for (let k = i; k < j; k++) {
          if (V[at(k + 1, j)] >= INF) continue;
          const cand = W[at(i, k)] + V[at(k + 1, j)] + auPenalty(pairOf(k + 1, j));
          if (cand < w) w = cand;
        }
        W[at(i, j)] = Math.min(w, 0);
      }
    }

    // ---- traceback
    const pairsOut = [];
    const struct = new Array(n).fill('.');
    const stack = [[0, n - 1, 'W']];
    while (stack.length) {
      const [i, j, kind] = stack.pop();
      if (i >= j) continue;
      if (kind === 'W') {
        if (Math.abs(W[at(i, j)] - W[at(i, j - 1)]) < 1e-9) { stack.push([i, j - 1, 'W']); continue; }
        if (V[at(i, j)] < INF && Math.abs(W[at(i, j)] - (V[at(i, j)] + auPenalty(pairOf(i, j)))) < 1e-9) { stack.push([i, j, 'V']); continue; }
        let done = false;
        for (let k = i; k < j && !done; k++) {
          if (V[at(k + 1, j)] >= INF) continue;
          if (Math.abs(W[at(i, j)] - (W[at(i, k)] + V[at(k + 1, j)] + auPenalty(pairOf(k + 1, j)))) < 1e-9) {
            stack.push([i, k, 'W']); stack.push([k + 1, j, 'V']); done = true;
          }
        }
        if (!done && Math.abs(W[at(i, j)]) < 1e-9) continue;
        continue;
      }
      if (kind === 'WM') {
        if (V[at(i, j)] < INF && Math.abs(WM[at(i, j)] - (V[at(i, j)] + ML_B + auPenalty(pairOf(i, j)))) < 1e-9) { stack.push([i, j, 'V']); continue; }
        if (WM[at(i + 1, j)] < INF && Math.abs(WM[at(i, j)] - (WM[at(i + 1, j)] + ML_C)) < 1e-9) { stack.push([i + 1, j, 'WM']); continue; }
        if (WM[at(i, j - 1)] < INF && Math.abs(WM[at(i, j)] - (WM[at(i, j - 1)] + ML_C)) < 1e-9) { stack.push([i, j - 1, 'WM']); continue; }
        let done = false;
        for (let k = i + 1; k < j && !done; k++) {
          if (WM[at(i, k)] < INF && WM[at(k + 1, j)] < INF && Math.abs(WM[at(i, j)] - (WM[at(i, k)] + WM[at(k + 1, j)])) < 1e-9) {
            stack.push([i, k, 'WM']); stack.push([k + 1, j, 'WM']); done = true;
          }
        }
        continue;
      }
      // kind === 'V'
      if (!canPair(i, j) || V[at(i, j)] >= INF) continue;
      struct[i] = '('; struct[j] = ')'; pairsOut.push([i, j]);
      const p = pairOf(i, j), target = V[at(i, j)];
      if (Math.abs(target - (hairpinE(j - i - 1) + auPenalty(p))) < 1e-9) continue;
      let done = false;
      for (let k = i + 1; k <= Math.min(i + MAXLOOP + 1, j - MINLOOP - 2) && !done; k++) {
        for (let l = Math.max(k + MINLOOP + 1, j - 1 - (MAXLOOP - (k - i - 1))); l < j && !done; l++) {
          if (V[at(k, l)] >= INF) continue;
          const u1 = k - i - 1, u2 = j - l - 1;
          if (u1 + u2 > MAXLOOP) continue;
          const q = pairOf(k, l);
          let e;
          if (u1 === 0 && u2 === 0) e = (RNA_STACK[p] && RNA_STACK[p][q] !== undefined) ? RNA_STACK[p][q] : -0.5;
          else if (u1 === 0 || u2 === 0) e = bulgeE(u1 + u2) + auPenalty(p) + auPenalty(q);
          else e = internalE(u1 + u2) + 0.6 * Math.abs(u1 - u2) + auPenalty(p) + auPenalty(q);
          if (Math.abs(target - (e + V[at(k, l)])) < 1e-9) { stack.push([k, l, 'V']); done = true; }
        }
      }
      if (done) continue;
      for (let k = i + 2; k < j - 1 && !done; k++) {
        if (WM[at(i + 1, k)] >= INF || WM[at(k + 1, j - 1)] >= INF) continue;
        if (Math.abs(target - (WM[at(i + 1, k)] + WM[at(k + 1, j - 1)] + ML_A + ML_B + auPenalty(p))) < 1e-9) {
          stack.push([i + 1, k, 'WM']); stack.push([k + 1, j - 1, 'WM']); done = true;
        }
      }
    }
    return { seq: s, structure: struct.join(''), dg: W[at(0, n - 1)], pairs: pairsOut.sort((a, b) => a[0] - b[0]), approximate: true };
  }

  // Groups a dot-bracket string into the hairpins it contains, for a readable summary.
  function stemList(structure) {
    const stack = [], stems = [];
    for (let i = 0; i < structure.length; i++) {
      if (structure[i] === '(') stack.push(i);
      else if (structure[i] === ')') {
        const j = stack.pop();
        const last = stems[stems.length - 1];
        if (last && last.i === j + 1 && last.j === i - 1) { last.i = j; last.j = i; last.len++; }
        else stems.push({ i: j, j: i, len: 1 });
      }
    }
    return stems.filter(s => s.len >= 2).sort((a, b) => a.i - b.i);
  }

  // ---------------------------------------------------------------- protein
  const AA_MONO = { A: 71.0788, R: 156.1875, N: 114.1038, D: 115.0886, C: 103.1388, E: 129.1155, Q: 128.1307, G: 57.0519, H: 137.1411, I: 113.1594, L: 113.1594, K: 128.1741, M: 131.1926, F: 147.1766, P: 97.1167, S: 87.0782, T: 101.1051, W: 186.2132, Y: 163.1760, V: 99.1326 };
  const KD = { I: 4.5, V: 4.2, L: 3.8, F: 2.8, C: 2.5, M: 1.9, A: 1.8, G: -0.4, T: -0.7, S: -0.8, W: -0.9, Y: -1.3, P: -1.6, H: -3.2, E: -3.5, Q: -3.5, D: -3.5, N: -3.5, K: -3.9, R: -4.5 };
  const AA_NAME = { A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', E: 'Glu', Q: 'Gln', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val' };
  // EMBOSS pKa set. Different tools use different sets and differ by a few tenths of a pH unit.
  const PKA = { cTerm: 3.6, nTerm: 8.6, D: 3.9, E: 4.1, C: 8.5, Y: 10.1, H: 6.5, K: 10.8, R: 12.5 };
  const PKA_SET = 'EMBOSS';

  function netCharge(counts, pH) {
    const pos = (pk, n) => n / (1 + Math.pow(10, pH - pk));
    const neg = (pk, n) => -n / (1 + Math.pow(10, pk - pH));
    return pos(PKA.nTerm, 1) + pos(PKA.K, counts.K || 0) + pos(PKA.R, counts.R || 0) + pos(PKA.H, counts.H || 0)
      + neg(PKA.cTerm, 1) + neg(PKA.D, counts.D || 0) + neg(PKA.E, counts.E || 0) + neg(PKA.C, counts.C || 0) + neg(PKA.Y, counts.Y || 0);
  }
  function proteinStats(protein) {
    const p = String(protein).toUpperCase().replace(/[^A-Z*]/g, '').replace(/\*+$/, '');
    if (!p) throw new Error('No protein sequence.');
    const unknown = [...new Set([...p].filter(c => !AA_MONO[c]))];
    const counts = {};
    for (const c of p) counts[c] = (counts[c] || 0) + 1;
    let mw = 18.0153;
    for (const c of p) mw += AA_MONO[c] || 0;
    let lo = 0, hi = 14;
    for (let i = 0; i < 100; i++) { const mid = (lo + hi) / 2; if (netCharge(counts, mid) > 0) lo = mid; else hi = mid; }
    const pI = (lo + hi) / 2;
    const W = counts.W || 0, Y = counts.Y || 0, C = counts.C || 0;
    const extRed = W * 5500 + Y * 1490;
    const extOx = extRed + Math.floor(C / 2) * 125;
    const gravy = p.length ? [...p].reduce((a, c) => a + (KD[c] || 0), 0) / p.length : 0;
    const aliphatic = p.length ? 100 * ((counts.A || 0) + 2.9 * (counts.V || 0) + 3.9 * ((counts.I || 0) + (counts.L || 0))) / p.length : 0;
    return {
      protein: p, length: p.length, mw, pI, pkaSet: PKA_SET, unknown,
      counts, gravy, aliphatic,
      extinction: { reduced: extRed, cystine: extOx, trp: W, tyr: Y, cys: C },
      a280: { reduced: extRed / mw, cystine: extOx / mw },
      charge7: netCharge(counts, 7.0),
      composition: Object.keys(AA_MONO).map(c => ({ aa: c, name: AA_NAME[c], n: counts[c] || 0, pct: p.length ? (counts[c] || 0) / p.length * 100 : 0 })).sort((a, b) => b.n - a.n)
    };
  }
  // Mean hydropathy over a sliding window, as Kyte and Doolittle described it.
  function hydropathy(protein, window = 9) {
    const p = String(protein).toUpperCase().replace(/[^A-Z]/g, '');
    const half = Math.floor(window / 2), out = [];
    for (let i = 0; i < p.length; i++) {
      let sum = 0, n = 0;
      for (let j = Math.max(0, i - half); j <= Math.min(p.length - 1, i + half); j++) { sum += KD[p[j]] || 0; n++; }
      out.push(sum / n);
    }
    return out;
  }
  function chargeCurve(protein, step = 0.25) {
    const p = String(protein).toUpperCase().replace(/[^A-Z]/g, '');
    const counts = {};
    for (const c of p) counts[c] = (counts[c] || 0) + 1;
    const out = [];
    for (let pH = 0; pH <= 14 + 1e-9; pH += step) out.push({ pH, charge: netCharge(counts, pH) });
    return out;
  }

  root.STRUCT = {
    DG, dgTable, duplexes, duplexDiagram, hairpins, hairpinDiagram, primerStructures, loopPenalty,
    rnaFold, stemList, toRNA,
    proteinStats, hydropathy, chargeCurve, KD, AA_NAME, PKA_SET
  };
})(typeof window !== 'undefined' ? window : globalThis);
