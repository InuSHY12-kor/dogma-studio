/* GeneCanvas — codon optimisation. Exposes window.CODON
 *
 * Works on a coding sequence and rewrites synonymous codons for a chosen organism,
 * while respecting regions the user wants left alone, removing unwanted restriction
 * sites and, where the protein allows it, writing wanted ones in silently.
 */
(function (root) {
  'use strict';
  const DNA = root.DNA;

  // ---------------------------------------------------------------- genetic code helpers
  const SYN = {};                      // amino acid -> its codons
  for (const [codon, aa] of Object.entries(DNA.CODON)) (SYN[aa] = SYN[aa] || []).push(codon);
  const IUPAC = { A: 'A', C: 'C', G: 'G', T: 'T', R: 'AG', Y: 'CT', S: 'GC', W: 'AT', K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT' };
  // For each amino acid, which bases can appear at each of the three codon positions.
  const AA_BASES = {};
  for (const [aa, list] of Object.entries(SYN)) {
    AA_BASES[aa] = [new Set(), new Set(), new Set()];
    for (const c of list) for (let i = 0; i < 3; i++) AA_BASES[aa][i].add(c[i]);
  }

  function translate(s) { return DNA.translate(s); }
  function codonsOf(s) { const out = []; for (let i = 0; i + 3 <= s.length; i += 3) out.push(s.substr(i, 3)); return out; }

  // ---------------------------------------------------------------- usage tables
  function tables() {
    const T = (root.CODON_TABLES || {}).tables || {};
    return Object.entries(T).map(([key, t]) => ({ key, ...t }));
  }
  const cache = new Map();
  // Per-codon frequency within its amino acid, and relative adaptiveness w used by CAI.
  function usage(key) {
    if (cache.has(key)) return cache.get(key);
    const src = (root.CODON_TABLES || {}).tables || {};
    const list = (root.CODON_TABLES || {}).codons || [];
    const t = src[key];
    if (!t) throw new Error(`No codon usage table for “${key}”.`);
    const count = {};
    list.forEach((c, i) => { count[c] = t.counts[i] || 0; });
    const frac = {}, w = {}, perK = {};
    const total = Object.values(count).reduce((a, b) => a + b, 0) || 1;
    for (const [aa, codons] of Object.entries(SYN)) {
      const sum = codons.reduce((a, c) => a + count[c], 0);
      const max = Math.max(...codons.map(c => count[c]));
      for (const c of codons) {
        frac[c] = sum ? count[c] / sum : 1 / codons.length;
        w[c] = max ? Math.max(count[c] / max, 1e-4) : 1;   // a floor keeps CAI finite
        perK[c] = count[c] / total * 1000;
      }
    }
    const u = { key, name: t.name, short: t.short, source: t.source, genes: t.genes, count, frac, w, perK };
    cache.set(key, u);
    return u;
  }

  // ---------------------------------------------------------------- measures
  function cai(seq, key) {
    const u = typeof key === 'string' ? usage(key) : key;
    let sum = 0, n = 0;
    for (const c of codonsOf(seq)) {
      const aa = DNA.CODON[c];
      if (!aa || aa === '*' || SYN[aa].length < 2) continue;   // no choice means no information
      sum += Math.log(u.w[c] || 1e-4); n++;
    }
    return n ? Math.exp(sum / n) : 0;
  }
  function gc3(seq) {
    let g = 0, n = 0;
    for (const c of codonsOf(seq)) { if (c[2] === 'G' || c[2] === 'C') g++; n++; }
    return n ? g / n * 100 : 0;
  }
  function longestRun(s) { let b = 0, r = 0, p = ''; for (const ch of s) { r = ch === p ? r + 1 : 1; p = ch; if (r > b) b = r; } return b; }

  // ---------------------------------------------------------------- validation
  // Checks that a stretch really is a reading frame before anything is changed.
  function inspectORF(seq) {
    const s = DNA.clean(seq);
    const issues = [];
    if (s.length < 6) issues.push('shorter than two codons');
    if (s.length % 3 !== 0) issues.push(`length ${s.length} is not a multiple of 3`);
    const starts = /^(ATG|GTG|TTG)/.test(s);
    if (!starts) issues.push('does not begin with a start codon (ATG, GTG or TTG)');
    const prot = translate(s);
    const stopAtEnd = prot.endsWith('*');
    const internal = prot.slice(0, -1).indexOf('*');
    if (internal >= 0) issues.push(`stop codon inside the reading frame at codon ${internal + 1}`);
    if (!stopAtEnd) issues.push('does not end with a stop codon');
    return {
      ok: issues.length === 0, issues, length: s.length, codons: Math.floor(s.length / 3),
      protein: prot, aa: stopAtEnd ? prot.length - 1 : prot.length, startsWith: s.slice(0, 3), stopsWith: s.slice(-3), hasStop: stopAtEnd, internalStop: internal >= 0
    };
  }

  // ---------------------------------------------------------------- site matching
  function siteRegex(site) { return new RegExp(site.split('').map(c => `[${IUPAC[c] || c}]`).join(''), 'g'); }
  // Every position where an enzyme cuts, on either strand.
  function findSite(seq, site) {
    const out = [];
    for (const [pat, strand] of [[site, 1], [DNA.rc(site), -1]]) {
      if (strand === -1 && DNA.rc(site) === site) continue;
      const re = siteRegex(pat);
      let m;
      while ((m = re.exec(seq))) { out.push({ start: m.index, end: m.index + pat.length, strand }); re.lastIndex = m.index + 1; }
    }
    return out.sort((a, b) => a.start - b.start);
  }
  function enzymeSite(name) {
    const e = DNA.ENZ_BY_NAME[name];
    if (!e) throw new Error(`Unknown enzyme ${name}.`);
    return e.site;
  }

  // ---------------------------------------------------------------- picking codons
  function makeRng(seed) {
    let a = (seed >>> 0) || 1;
    return () => { a ^= a << 13; a >>>= 0; a ^= a >> 17; a ^= a << 5; a >>>= 0; return a / 4294967296; };
  }
  function pickBest(aa, u) {
    let best = null, bv = -1;
    for (const c of SYN[aa]) if (u.frac[c] > bv) { bv = u.frac[c]; best = c; }
    return best;
  }
  function pickWeighted(aa, u, rng, floor) {
    const codons = SYN[aa];
    const ws = codons.map(c => (floor && u.frac[c] < floor ? 0 : u.frac[c]));
    let total = ws.reduce((a, b) => a + b, 0);
    if (total <= 0) return pickBest(aa, u);
    let r = rng() * total;
    for (let i = 0; i < codons.length; i++) { r -= ws[i]; if (r <= 0) return codons[i]; }
    return codons[codons.length - 1];
  }
  // Keeps the source organism's rank: a codon that was rare there stays rare here.
  function pickHarmonised(aa, original, src, u) {
    const codons = SYN[aa];
    if (!src) return pickBest(aa, u);
    const target = src.frac[original] !== undefined ? src.frac[original] : 1 / codons.length;
    let best = codons[0], bd = Infinity;
    for (const c of codons) { const d = Math.abs(u.frac[c] - target); if (d < bd) { bd = d; best = c; } }
    return best;
  }

  const ALGORITHMS = [
    { key: 'best', name: 'Most frequent codon', note: 'Always uses the organism’s commonest codon for each amino acid. Highest CAI, but it repeats the same codons and can strain single tRNA pools.' },
    { key: 'threshold', name: 'Frequency sampling, rare codons dropped', note: 'Picks codons in proportion to their natural frequency after discarding those below the cutoff. A good default for expression.' },
    { key: 'weighted', name: 'Frequency sampling', note: 'Picks codons in proportion to their natural frequency in the host, keeping a natural-looking mix.' },
    { key: 'harmonised', name: 'Harmonised to the source organism', note: 'Gives each codon a frequency in the host close to what it had in the source, preserving translation pauses that can matter for folding.' }
  ];

  // ---------------------------------------------------------------- main
  function optimize(cds, opts = {}) {
    const seq = DNA.clean(cds);
    const u = usage(opts.organism);
    const src = opts.sourceOrganism ? usage(opts.sourceOrganism) : null;
    const algo = opts.algorithm || 'threshold';
    const floor = opts.rareCutoff !== undefined ? opts.rareCutoff : 0.10;
    const rng = makeRng(opts.seed || 20260101);
    const keepStart = opts.keepStart !== false;
    const keepStop = opts.keepStop !== false;
    const n = seq.length;
    const nCodons = Math.floor(n / 3);

    // which codons may not be touched
    const locked = new Array(nCodons).fill(false);
    for (const r of (opts.preserve || [])) {
      const a = Math.max(0, Math.floor(r.start / 3)), b = Math.min(nCodons - 1, Math.floor((r.end - 1) / 3));
      for (let i = a; i <= b; i++) locked[i] = true;
    }
    const original = codonsOf(seq);
    if (keepStart && nCodons) locked[0] = true;
    if (keepStop && nCodons && DNA.CODON[original[nCodons - 1]] === '*') locked[nCodons - 1] = true;

    const out = original.slice();
    for (let i = 0; i < nCodons; i++) {
      if (locked[i]) continue;
      const aa = DNA.CODON[original[i]];
      if (!aa || aa === 'X') continue;
      if (aa === '*') { out[i] = opts.preferredStop ? pickBest('*', u) : original[i]; continue; }
      if (SYN[aa].length < 2) continue;
      out[i] = algo === 'best' ? pickBest(aa, u)
        : algo === 'weighted' ? pickWeighted(aa, u, rng, 0)
          : algo === 'harmonised' ? pickHarmonised(aa, original[i], src, u)
            : pickWeighted(aa, u, rng, floor);
    }

    const log = [];
    let cur = out.join('');
    const protein = translate(seq);

    // --- avoid unwanted sites --------------------------------------------------
    const avoid = (opts.avoidEnzymes || []).map(name => ({ name, site: enzymeSite(name) }));
    const extraAvoid = (opts.avoidSequences || []).filter(Boolean).map(s => ({ name: s, site: DNA.clean(s) }));
    const allAvoid = [...avoid, ...extraAvoid];
    const removal = removeSites(cur, allAvoid, locked, u, protein, log);
    cur = removal.seq;

    // --- write wanted sites in silently ----------------------------------------
    const added = [];
    for (const name of (opts.addEnzymes || [])) {
      const site = enzymeSite(name);
      const r = addSite(cur, site, locked, u, protein, opts.addNear);
      if (r) {
        cur = r.seq;
        added.push({ name, site, position: r.position, codonsChanged: r.codonsChanged });
        log.push(`Wrote a ${name} site at ${r.position + 1} by changing ${r.codonsChanged} codon${r.codonsChanged === 1 ? '' : 's'}.`);
        // the new site must not be undone later, so lock the codons it spans
        for (let i = Math.floor(r.position / 3); i <= Math.floor((r.position + site.length - 1) / 3); i++) locked[i] = true;
      } else {
        added.push({ name, site, position: null, codonsChanged: 0 });
        log.push(`Could not write a ${name} site anywhere without changing the protein.`);
      }
    }

    // --- break up long single-base runs ----------------------------------------
    const maxRun = opts.maxRun || 0;
    if (maxRun >= 4) cur = breakRuns(cur, maxRun, locked, u, protein, log);

    if (translate(cur) !== protein) throw new Error('Internal error: the optimised sequence does not encode the same protein.');

    const changes = [];
    const after = codonsOf(cur);
    for (let i = 0; i < nCodons; i++) if (after[i] !== original[i]) changes.push({ codon: i, aa: DNA.CODON[original[i]], from: original[i], to: after[i], fromFreq: u.frac[original[i]], toFreq: u.frac[after[i]] });

    const remaining = [];
    for (const a of allAvoid) { const hits = findSite(cur, a.site); if (hits.length) remaining.push({ name: a.name, count: hits.length, at: hits.map(h => h.start + 1) }); }

    return {
      seq: cur, protein, organism: u, algorithm: algo, log, changes, added, remaining,
      locked: locked.map((v, i) => v ? i : -1).filter(i => i >= 0),
      before: { cai: cai(seq, u), gc: DNA.gcPct(seq), gc3: gc3(seq), run: longestRun(seq), sites: countSites(seq, allAvoid) },
      after: { cai: cai(cur, u), gc: DNA.gcPct(cur), gc3: gc3(cur), run: longestRun(cur), sites: countSites(cur, allAvoid) },
      changedCodons: changes.length, totalCodons: nCodons
    };
  }
  function countSites(seq, list) { let n = 0; for (const a of list) n += findSite(seq, a.site).length; return n; }

  // Replaces codons until an unwanted site is gone, preferring the smallest loss of CAI.
  function removeSites(seq, list, locked, u, protein, log) {
    let cur = seq;
    for (const a of list) {
      if (!a.site) continue;
      for (let guard = 0; guard < 200; guard++) {
        const hits = findSite(cur, a.site);
        if (!hits.length) break;
        const h = hits[0];
        const fixed = killSiteAt(cur, h.start, h.end, locked, u, protein);
        if (!fixed) {
          log.push(`Could not remove the ${a.name} site at ${h.start + 1}: every synonymous change there is blocked or would alter the protein.`);
          // skip past this one by masking it out of further attempts
          break;
        }
        cur = fixed.seq;
        log.push(`Removed a ${a.name} site at ${h.start + 1} by changing codon ${fixed.codon + 1} (${fixed.from} → ${fixed.to}).`);
      }
    }
    return { seq: cur };
  }
  function killSiteAt(seq, start, end, locked, u, protein) {
    const first = Math.floor(start / 3), last = Math.floor((end - 1) / 3);
    let best = null;
    for (let i = first; i <= last; i++) {
      if (locked[i]) continue;
      const cur = seq.substr(i * 3, 3);
      const aa = DNA.CODON[cur];
      if (!aa || SYN[aa].length < 2) continue;
      for (const alt of SYN[aa]) {
        if (alt === cur) continue;
        const cand = seq.slice(0, i * 3) + alt + seq.slice(i * 3 + 3);
        if (cand.slice(start, end) === seq.slice(start, end)) continue;   // must actually change the site
        const loss = (u.frac[cur] || 0) - (u.frac[alt] || 0);
        if (!best || loss < best.loss) best = { seq: cand, codon: i, from: cur, to: alt, loss };
      }
    }
    if (!best) return null;
    if (translate(best.seq) !== protein) return null;
    return best;
  }

  // Finds a place where synonymous codons can spell a site, then writes it there.
  function addSite(seq, site, locked, u, protein, near) {
    const L = site.length, n = seq.length;
    const allowed = site.split('').map(c => new Set((IUPAC[c] || c).split('')));
    let best = null;
    for (let p = 0; p + L <= n; p++) {
      const first = Math.floor(p / 3), last = Math.floor((p + L - 1) / 3);
      let blocked = false;
      for (let i = first; i <= last; i++) if (locked[i]) { blocked = true; break; }
      if (blocked) continue;
      // quick test: can the amino acids here produce the required base at every position?
      let feasible = true;
      for (let k = 0; k < L && feasible; k++) {
        const abs = p + k, ci = Math.floor(abs / 3), off = abs % 3;
        const aa = DNA.CODON[seq.substr(ci * 3, 3)];
        if (!aa || !AA_BASES[aa]) { feasible = false; break; }
        let any = false;
        for (const b of AA_BASES[aa][off]) if (allowed[k].has(b)) { any = true; break; }
        if (!any) feasible = false;
      }
      if (!feasible) continue;
      // enumerate the synonymous combinations across the few codons involved
      const idxs = [];
      for (let i = first; i <= last; i++) idxs.push(i);
      const options = idxs.map(i => SYN[DNA.CODON[seq.substr(i * 3, 3)]] || []);
      const combos = options.reduce((a, o) => a * o.length, 1);
      if (combos > 20000) continue;
      const cand = tryCombos(seq, idxs, options, p, L, site);
      if (!cand) continue;
      let changed = 0, loss = 0;
      for (const i of idxs) {
        const a = seq.substr(i * 3, 3), b = cand.substr(i * 3, 3);
        if (a !== b) { changed++; loss += (u.frac[a] || 0) - (u.frac[b] || 0); }
      }
      const distance = near === undefined ? 0 : Math.abs(p - near);
      const score = loss + changed * 0.05 + distance * 1e-5;
      if (!best || score < best.score) best = { seq: cand, position: p, codonsChanged: changed, score };
    }
    if (!best) return null;
    if (translate(best.seq) !== protein) return null;
    return best;
  }
  function tryCombos(seq, idxs, options, p, L, site) {
    const re = new RegExp('^' + site.split('').map(c => `[${IUPAC[c] || c}]`).join('') + '$');
    const pick = new Array(idxs.length).fill(0);
    const total = options.reduce((a, o) => a * o.length, 1);
    for (let t = 0; t < total; t++) {
      let r = t;
      for (let i = 0; i < idxs.length; i++) { pick[i] = r % options[i].length; r = Math.floor(r / options[i].length); }
      let cand = seq;
      for (let i = 0; i < idxs.length; i++) cand = cand.slice(0, idxs[i] * 3) + options[i][pick[i]] + cand.slice(idxs[i] * 3 + 3);
      if (re.test(cand.substr(p, L))) return cand;
    }
    return null;
  }

  // Shortens runs such as AAAAAA that can cause slippage during synthesis.
  function breakRuns(seq, maxRun, locked, u, protein, log) {
    let cur = seq;
    for (let guard = 0; guard < 300; guard++) {
      const m = cur.match(new RegExp(`(.)\\1{${maxRun},}`));
      if (!m) break;
      const at = m.index, len = m[0].length;
      const first = Math.floor(at / 3), last = Math.floor((at + len - 1) / 3);
      let done = null;
      for (let i = first; i <= last && !done; i++) {
        if (locked[i]) continue;
        const c = cur.substr(i * 3, 3), aa = DNA.CODON[c];
        if (!aa || SYN[aa].length < 2) continue;
        let bestAlt = null, bestLoss = Infinity;
        for (const alt of SYN[aa]) {
          if (alt === c) continue;
          const cand = cur.slice(0, i * 3) + alt + cur.slice(i * 3 + 3);
          if (longestRun(cand.slice(Math.max(0, at - 3), at + len + 3)) >= len) continue;
          const loss = (u.frac[c] || 0) - (u.frac[alt] || 0);
          if (loss < bestLoss) { bestLoss = loss; bestAlt = cand; }
        }
        if (bestAlt) done = bestAlt;
      }
      if (!done) { log.push(`Could not shorten the run of ${len} ${m[1]}s at ${at + 1}.`); break; }
      cur = done;
    }
    if (translate(cur) !== protein) return seq;
    return cur;
  }

  // Builds a coding sequence from a protein, then optimises it.
  function fromProtein(protein, opts = {}) {
    const p = String(protein).toUpperCase().replace(/[^A-Z*]/g, '');
    if (!p) throw new Error('Enter a protein sequence.');
    const bad = [...p].find(c => !SYN[c]);
    if (bad) throw new Error(`“${bad}” is not an amino acid letter.`);
    const u = usage(opts.organism);
    let cds = '';
    for (const aa of p) cds += pickBest(aa, u);
    if (!p.endsWith('*') && opts.addStop !== false) cds += pickBest('*', u);
    if (!cds.startsWith('ATG') && opts.addMet !== false && p[0] !== 'M') cds = 'ATG' + cds;
    return optimize(cds, opts);
  }

  root.CODON = { tables, usage, cai, gc3, optimize, fromProtein, inspectORF, findSite, enzymeSite, SYN, ALGORITHMS, longestRun, codonsOf };
})(typeof window !== 'undefined' ? window : globalThis);
