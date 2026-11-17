/* Dogma Studio — silent mutagenesis of restriction sites. Exposes window.SILENT
 *
 * Two jobs, both done without changing a single amino acid:
 *   addSites()    finds every place inside a coding sequence where an enzyme site
 *                 can be written in using alternative synonymous codons.
 *   removeSite()  lists the synonymous codon swaps that destroy one chosen site.
 *
 * All of the work happens on the coding strand of the feature, so a reverse-strand
 * CDS is handled by reverse-complementing first and mapping the answer back at the end.
 */
(function (root) {
  'use strict';
  const DNA = root.DNA;

  const SYN = {};                        // amino acid -> its codons
  for (const [codon, aa] of Object.entries(DNA.CODON)) (SYN[aa] = SYN[aa] || []).push(codon);
  const IUPAC = { A: 'A', C: 'C', G: 'G', T: 'T', R: 'AG', Y: 'CT', S: 'GC', W: 'AT', K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT' };
  const matches = (base, code) => (IUPAC[code] || code).indexOf(base) >= 0;
  const matchesAt = (seq, at, pat) => {
    for (let i = 0; i < pat.length; i++) if (!matches(seq[at + i], pat[i])) return false;
    return true;
  };
  const diff = (a, b) => { let n = 0; for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) n++; return n; };

  // Patterns an enzyme can be recognised by: the site itself, plus its reverse
  // complement when the site is not palindromic.
  function patterns(site) {
    const r = DNA.rc(site);
    return r === site ? [site] : [site, r];
  }

  // ---------------------------------------------------------------- writing a site in
  // Can the window [i, i+L) of the coding sequence be made to read `pat` using only
  // synonymous codons? Returns the cheapest replacement codons, or null.
  function fitPattern(cds, i, pat) {
    const L = pat.length, c0 = Math.floor(i / 3), c1 = Math.floor((i + L - 1) / 3);
    if ((c1 + 1) * 3 > cds.length) return null;
    const origs = [], options = [];
    for (let c = c0; c <= c1; c++) {
      const codon = cds.substr(c * 3, 3), aa = DNA.CODON[codon];
      if (!aa || aa === '*' || aa === 'X') return null;   // stop codons and ambiguity are left alone
      origs.push(codon);
      options.push(SYN[aa]);
    }
    let best = null, bestCost = Infinity;
    const pick = [];
    (function dfs(k, cost) {
      if (cost >= bestCost) return;
      if (k === options.length) { best = pick.slice(); bestCost = cost; return; }
      for (const cand of options[k]) {
        // every base of this codon that falls inside the window must match the pattern
        let ok = true;
        const base = (c0 + k) * 3;
        for (let j = 0; j < 3; j++) {
          const p = base + j - i;
          if (p >= 0 && p < L && !matches(cand[j], pat[p])) { ok = false; break; }
        }
        if (!ok) continue;
        pick[k] = cand;
        dfs(k + 1, cost + diff(cand, origs[k]));
      }
    })(0, 0);
    if (!best) return null;
    return { firstCodon: c0, codons: best, origs, changes: bestCost };
  }

  function applyCodons(cds, firstCodon, codons) {
    return cds.slice(0, firstCodon * 3) + codons.join('') + cds.slice((firstCodon + codons.length) * 3);
  }

  /* Find new enzyme sites that can be added silently.
   *   cds       the coding sequence, in frame, on its own coding strand
   *   enzymes   array of enzyme records ({name, site, core, ...} from DNA.ENZYMES)
   *   opts      maxChanges, minSite, maxSite, existing ('any'|'none'|'noneInCDS'),
   *             fullSeq (the whole molecule, for counting sites already present),
   *             limit (stop after this many candidates)
   */
  function addSites(cds, enzymes, opts = {}) {
    const maxChanges = opts.maxChanges == null ? 3 : opts.maxChanges;
    const limit = opts.limit || 4000;
    const whole = opts.fullSeq || cds;
    const out = [];
    for (const e of enzymes) {
      const L = e.site.length;
      if (L < (opts.minSite || 0) || L > (opts.maxSite || 99)) continue;
      const pats = patterns(e.site);
      const inWhole = countPat(whole, pats, !!opts.circular);
      const inCDS = countPat(cds, pats, false);
      if (opts.existing === 'none' && inWhole) continue;
      if (opts.existing === 'noneInCDS' && inCDS) continue;
      for (let i = 0; i + L <= cds.length; i++) {
        let hit = null;
        for (const pat of pats) {
          const f = fitPattern(cds, i, pat);
          if (f && f.changes > 0 && f.changes <= maxChanges && (!hit || f.changes < hit.changes)) hit = Object.assign({ pat }, f);
        }
        if (!hit) continue;
        const mutated = applyCodons(cds, hit.firstCodon, hit.codons);
        out.push({
          enzyme: e.name, site: e.site, core: e.core, start: i, end: i + L,
          changes: hit.changes, firstCodon: hit.firstCodon,
          codons: hit.codons, origs: hit.origs, mutatedCDS: mutated,
          existingInSequence: inWhole, existingInCDS: inCDS
        });
        if (out.length >= limit) return out;
      }
    }
    return out;
  }

  function countPat(seq, pats, circular) {
    const text = circular ? seq + seq.slice(0, 10) : seq;
    let n = 0;
    for (const pat of pats) {
      const re = new RegExp(pat.split('').map(c => `[${IUPAC[c] || c}]`).join(''), 'g');
      let m;
      while ((m = re.exec(text))) { if (m.index < seq.length) n++; re.lastIndex = m.index + 1; }
    }
    return n;
  }

  // ---------------------------------------------------------------- removing a site
  /* For a site occupying [start, end) of the coding sequence, list each overlapping
   * codon and every synonymous alternative, marked with whether it kills the site. */
  function removeSite(cds, start, end, site) {
    const pats = patterns(site);
    const c0 = Math.floor(start / 3), c1 = Math.floor((end - 1) / 3), out = [];
    for (let c = c0; c <= c1; c++) {
      const codon = cds.substr(c * 3, 3), aa = DNA.CODON[codon];
      if (!aa || aa === '*') continue;
      const alts = (SYN[aa] || []).filter(x => x !== codon).map(x => {
        const mutated = applyCodons(cds, c, [x]);
        const still = pats.some(p => matchesAt(mutated, start, p));
        return { codon: x, changes: diff(x, codon), kills: !still, mutatedCDS: mutated };
      });
      out.push({ codonIndex: c, codon, aa, aa3: DNA.AA3[aa] || aa, alts, canKill: alts.some(a => a.kills) });
    }
    return out;
  }

  root.SILENT = { addSites, removeSite, patterns, applyCodons, fitPattern, countPat, SYN };
})(window);
