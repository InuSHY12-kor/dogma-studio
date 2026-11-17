/* Dogma Studio — comparing whole sequences to one another. Exposes window.MSA
 *
 * Plasmids being compared are usually near-identical, so the alignment is anchored on
 * exact unique matches and only the stretches between anchors go through dynamic
 * programming. That keeps a multi-kb comparison fast without giving up accuracy.
 */
(function (root) {
  'use strict';
  const DNA = root.DNA;

  const MAX_CELLS = 4e6;      // largest DP block to attempt before anchoring again
  const KS = [24, 16, 11, 8]; // anchor sizes, tried from coarse to fine

  // Positions of k-mers that occur exactly once in the text, which makes an anchor trustworthy.
  function uniqueIndex(text, k) {
    const idx = new Map();
    for (let i = 0; i + k <= text.length; i++) {
      const s = text.substr(i, k);
      idx.set(s, idx.has(s) ? -1 : i);
    }
    return idx;
  }

  // Anchor points shared by both sequences, kept in increasing order on both.
  function chain(ref, qry, k) {
    if (ref.length < k || qry.length < k) return [];
    const idx = uniqueIndex(ref, k);
    const seen = new Map();
    const pts = [];
    for (let j = 0; j + k <= qry.length; j++) {
      const s = qry.substr(j, k);
      if (seen.has(s)) continue;          // only unique on the query side too
      seen.set(s, 1);
      const r = idx.get(s);
      if (r !== undefined && r >= 0) pts.push([r, j]);
    }
    if (!pts.length) return [];
    pts.sort((a, b) => a[1] - b[1] || a[0] - b[0]);
    // longest increasing subsequence on the reference coordinate
    const tails = [], tailIdx = [], prev = new Array(pts.length).fill(-1);
    for (let i = 0; i < pts.length; i++) {
      const r = pts[i][0];
      let lo = 0, hi = tails.length;
      while (lo < hi) { const mid = (lo + hi) >> 1; if (tails[mid] < r) lo = mid + 1; else hi = mid; }
      tails[lo] = r; tailIdx[lo] = i;
      prev[i] = lo > 0 ? tailIdx[lo - 1] : -1;
    }
    const out = [];
    for (let i = tailIdx.length ? tailIdx[tailIdx.length - 1] : -1; i >= 0; i = prev[i]) out.push(pts[i]);
    out.reverse();
    // merge points that sit on the same diagonal into one run
    const runs = [];
    for (const [r, q] of out) {
      const last = runs[runs.length - 1];
      if (last && r - last.r === q - last.q && r <= last.r + last.len) {
        const end = Math.max(last.r + last.len, r + k);
        last.len = end - last.r;
      } else runs.push({ r, q, len: k });
    }
    // trim any run that would overlap the previous one on either sequence
    const clean = [];
    for (const run of runs) {
      const last = clean[clean.length - 1];
      if (last && (run.r < last.r + last.len || run.q < last.q + last.len)) continue;
      clean.push(run);
    }
    return clean;
  }

  const MATCH = 2, MISMATCH = -3, GAP = -5;
  // Full global alignment, used only on the short stretches between anchors.
  function nw(a, b) {
    const m = a.length, n = b.length;
    if (!m && !n) return ['', ''];
    if (!m) return ['-'.repeat(n), b];
    if (!n) return [a, '-'.repeat(m)];
    const W = n + 1;
    const score = new Int32Array((m + 1) * W);
    const ptr = new Uint8Array((m + 1) * W);
    for (let j = 1; j <= n; j++) { score[j] = GAP * j; ptr[j] = 2; }
    for (let i = 1; i <= m; i++) { score[i * W] = GAP * i; ptr[i * W] = 1; }
    for (let i = 1; i <= m; i++) {
      const row = i * W, rowP = (i - 1) * W, ac = a.charCodeAt(i - 1);
      for (let j = 1; j <= n; j++) {
        const d = score[rowP + j - 1] + (b.charCodeAt(j - 1) === ac ? MATCH : MISMATCH);
        const u = score[rowP + j] + GAP;
        const l = score[row + j - 1] + GAP;
        let best = d, p = 0;
        if (u > best) { best = u; p = 1; }
        if (l > best) { best = l; p = 2; }
        score[row + j] = best; ptr[row + j] = p;
      }
    }
    let i = m, j = n, A = '', B = '';
    while (i > 0 || j > 0) {
      const p = i === 0 ? 2 : (j === 0 ? 1 : ptr[i * W + j]);
      if (p === 0) { A += a[i - 1]; B += b[j - 1]; i--; j--; }
      else if (p === 1) { A += a[i - 1]; B += '-'; i--; }
      else { A += '-'; B += b[j - 1]; j--; }
    }
    return [[...A].reverse().join(''), [...B].reverse().join('')];
  }

  // Aligns two stretches, anchoring again on finer k-mers when the block is too big for DP.
  function alignBlock(a, b, depth) {
    if (!a.length || !b.length || a.length * b.length <= MAX_CELLS) {
      if (a.length * b.length <= MAX_CELLS) return nw(a, b);
    }
    const k = KS[depth];
    if (k === undefined) {
      // Too divergent to align sensibly: show the two stretches one after the other.
      return [a + '-'.repeat(b.length), '-'.repeat(a.length) + b];
    }
    const runs = chain(a, b, k);
    if (!runs.length) return alignBlock(a, b, depth + 1);
    let A = '', B = '', ra = 0, qb = 0;
    for (const run of runs) {
      const [x, y] = alignBlock(a.slice(ra, run.r), b.slice(qb, run.q), depth + 1);
      A += x; B += y;
      A += a.substr(run.r, run.len); B += b.substr(run.q, run.len);
      ra = run.r + run.len; qb = run.q + run.len;
    }
    const [x, y] = alignBlock(a.slice(ra), b.slice(qb), depth + 1);
    return [A + x, B + y];
  }

  // The rotation of a circular query that lines it up best with the reference.
  function bestRotation(ref, qry) {
    const k = 16;
    if (qry.length < k * 2) return 0;
    const idx = uniqueIndex(ref, k);
    const votes = new Map();
    for (let j = 0; j + k <= qry.length; j += 2) {
      const r = idx.get(qry.substr(j, k));
      if (r === undefined || r < 0) continue;
      const d = ((r - j) % qry.length + qry.length) % qry.length;
      votes.set(d, (votes.get(d) || 0) + 1);
    }
    if (!votes.size) return 0;
    let best = 0, bv = 0;
    for (const [d, v] of votes) if (v > bv) { bv = v; best = d; }
    return best;
  }
  const rotate = (s, k) => { k = ((k % s.length) + s.length) % s.length; return s.slice(k) + s.slice(0, k); };

  // Aligns one sequence against the reference, trying both strands and, for circles, rotation.
  function pairAlign(ref, qry, opts = {}) {
    const tries = [];
    const add = (seq, strand, rotation) => tries.push({ seq, strand, rotation });
    if (opts.circular) {
      // The vote can land a few bases off when indels split the diagonal, so the
      // unrotated sequence stays in the running and the better alignment wins.
      const rc = DNA.rc(qry);
      const d1 = bestRotation(ref, qry), d2 = bestRotation(ref, rc);
      add(qry, 1, 0);
      if (d1) add(rotate(qry, -d1), 1, d1);
      add(rc, -1, 0);
      if (d2) add(rotate(rc, -d2), -1, d2);
    } else {
      add(qry, 1, 0);
      if (opts.bothStrands !== false) add(DNA.rc(qry), -1, 0);
    }
    let best = null;
    for (const t of tries) {
      const [refAln, qryAln] = alignBlock(ref, t.seq, 0);
      let match = 0;
      for (let i = 0; i < refAln.length; i++) if (refAln[i] === qryAln[i] && refAln[i] !== '-') match++;
      // prefer more matches, and among equals the alignment with fewer columns (fewer gaps)
      if (!best || match > best.match || (match === best.match && refAln.length < best.refAln.length)) best = { ...t, refAln, qryAln, match };
    }
    const cols = best.refAln.length;
    return { refAln: best.refAln, qryAln: best.qryAln, strand: best.strand, rotation: best.rotation, identity: cols ? best.match / cols * 100 : 0, matches: best.match, columns: cols };
  }

  // Builds one column layout holding the reference and every other sequence.
  function compare(seqs, opts = {}) {
    if (seqs.length < 2) throw new Error('Choose at least two sequences to compare.');
    const ref = seqs[0], n = ref.seq.length;
    const pairs = seqs.slice(1).map(s => {
      const a = pairAlign(ref.seq, s.seq, { circular: ref.circular && s.circular, bothStrands: opts.bothStrands });
      // per reference position: the aligned base, plus anything inserted before it
      const at = new Array(n).fill('-');
      const ins = new Array(n + 1).fill('');
      let ri = 0;
      for (let k = 0; k < a.refAln.length; k++) {
        const r = a.refAln[k], q = a.qryAln[k];
        if (r === '-') ins[Math.min(ri, n)] += q;
        else { if (ri < n) at[ri] = q; ri++; }
      }
      return { seq: s, aln: a, at, ins };
    });

    const insWidth = new Array(n + 1).fill(0);
    for (const p of pairs) for (let i = 0; i <= n; i++) insWidth[i] = Math.max(insWidth[i], p.ins[i].length);

    const refRow = [], refPos = [], rows = pairs.map(() => []);
    const push = (chars, pos) => { refRow.push(chars[0]); refPos.push(pos); for (let r = 0; r < rows.length; r++) rows[r].push(chars[r + 1]); };
    for (let i = 0; i <= n; i++) {
      for (let w = 0; w < insWidth[i]; w++) {
        const chars = ['-'];
        for (const p of pairs) chars.push(p.ins[i][w] || '-');
        push(chars, null);
      }
      if (i < n) {
        const chars = [ref.seq[i]];
        for (const p of pairs) chars.push(p.at[i]);
        push(chars, i);
      }
    }

    const stats = pairs.map((p, r) => {
      let sub = 0, del = 0, insC = 0, same = 0;
      const diffs = [];
      for (let c = 0; c < refRow.length; c++) {
        const a = refRow[c], b = rows[r][c];
        if (a === '-' && b === '-') continue;
        if (a === '-') { insC++; diffs.push({ col: c, pos: nearestPos(refPos, c), type: 'ins', ref: '', alt: b }); }
        else if (b === '-') { del++; diffs.push({ col: c, pos: refPos[c], type: 'del', ref: a, alt: '' }); }
        else if (a !== b) { sub++; diffs.push({ col: c, pos: refPos[c], type: 'sub', ref: a, alt: b }); }
        else same++;
      }
      return {
        name: p.seq.name, strand: p.aln.strand, rotation: p.aln.rotation,
        identity: p.aln.identity, substitutions: sub, deletions: del, insertions: insC,
        matches: same, diffs, length: p.seq.seq.length
      };
    });
    return { refName: ref.name, refLength: n, columns: refRow.length, refRow: refRow.join(''), refPos, rows: rows.map(r => r.join('')), names: pairs.map(p => p.seq.name), stats };
  }
  function nearestPos(refPos, c) {
    for (let i = c; i < refPos.length; i++) if (refPos[i] !== null) return refPos[i];
    for (let i = c; i >= 0; i--) if (refPos[i] !== null) return refPos[i];
    return 0;
  }

  root.MSA = { compare, pairAlign, chain, nw, bestRotation, rotate };
})(typeof window !== 'undefined' ? window : globalThis);
