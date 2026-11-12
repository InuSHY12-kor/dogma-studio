/* GeneCanvas — DNA logic (no DOM). Exposes window.DNA */
(function (root) {
  'use strict';

  const COMP = { A: 'T', T: 'A', G: 'C', C: 'G', U: 'A', N: 'N', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', V: 'B', D: 'H', H: 'D' };
  const VALID = /[^ACGTURYSWKMBDHVN]/g;

  function clean(s) { return String(s || '').toUpperCase().replace(VALID, '').replace(/U/g, 'T'); }
  function rc(s) { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] || 'N'; return o; }
  function comp(s) { let o = ''; for (let i = 0; i < s.length; i++) o += COMP[s[i]] || 'N'; return o; }

  // ---------- Genetic code ----------
  const BASES = 'TCAG', AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const CODON = {};
  for (let i = 0; i < 64; i++) CODON[BASES[i >> 4] + BASES[(i >> 2) & 3] + BASES[i & 3]] = AAS[i];
  function translate(s) { let o = ''; for (let i = 0; i + 2 < s.length; i += 3) o += CODON[s.substr(i, 3)] || 'X'; return o; }
  const AA3 = { A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Stop', X: 'Xaa' };
  const AA_MW = { A: 71.08, R: 156.19, N: 114.1, D: 115.09, C: 103.14, Q: 128.13, E: 129.12, G: 57.05, H: 137.14, I: 113.16, L: 113.16, K: 128.17, M: 131.19, F: 147.18, P: 97.12, S: 87.08, T: 101.1, W: 186.21, Y: 163.18, V: 99.13 };
  function proteinMW(p) { let m = 18.02; for (const c of p) m += AA_MW[c] || 0; return m; }

  function gcPct(s) { if (!s.length) return 0; let g = 0; for (const c of s) if (c === 'G' || c === 'C') g++; return g / s.length * 100; }
  function dsDnaMW(n) { return n * 617.96 + 36.04; }

  // ---------- Melting temperature (SantaLucia 1998 nearest neighbour) ----------
  const NN = { AA: [-7.9, -22.2], TT: [-7.9, -22.2], AT: [-7.2, -20.4], TA: [-7.2, -21.3], CA: [-8.5, -22.7], TG: [-8.5, -22.7], GT: [-8.4, -22.4], AC: [-8.4, -22.4], CT: [-7.8, -21.0], AG: [-7.8, -21.0], GA: [-8.2, -22.2], TC: [-8.2, -22.2], CG: [-10.6, -27.2], GC: [-9.8, -24.4], GG: [-8.0, -19.9], CC: [-8.0, -19.9] };
  function tm(s, na = 0.05, ct = 250e-9) {
    s = String(s).toUpperCase().replace(/[^ACGT]/g, '');
    const N = s.length;
    if (N < 2) return 0;
    if (N < 14) { const g = (s.match(/[GC]/g) || []).length; return 2 * (N - g) + 4 * g; }
    let dH = 0, dS = 0;
    for (let i = 0; i < N - 1; i++) { const p = NN[s.substr(i, 2)]; dH += p[0]; dS += p[1]; }
    for (const e of [s[0], s[N - 1]]) { if (e === 'G' || e === 'C') { dH += 0.1; dS += -2.8; } else { dH += 2.3; dS += 4.1; } }
    dS += 0.368 * (N - 1) * Math.log(na);
    return dH * 1000 / (dS + 1.987 * Math.log(ct / 4)) - 273.15;
  }

  // ---------- Restriction enzymes ----------
  // name site topCut bottomCut  (cut offsets measured from the left end of the site on the top strand)
  const ENZ_RAW = `AatII GACGTC 5 1
AccI GTMKAC 2 4
AflII CTTAAG 1 5
AgeI ACCGGT 1 5
AhdI GACNNNNNGTC 6 5
AleI CACNNNNGTG 5 5
AluI AGCT 2 2
ApaI GGGCCC 5 1
ApaLI GTGCAC 1 5
AscI GGCGCGCC 2 6
AvrII CCTAGG 1 5
BamHI GGATCC 1 5
BanII GRGCYC 5 1
BbsI GAAGAC 8 12
BclI TGATCA 1 5
BglII AGATCT 1 5
BmtI GCTAGC 5 1
BsaI GGTCTC 7 11
BsiWI CGTACG 1 5
BsmBI CGTCTC 7 11
BspEI TCCGGA 1 5
BsrGI TGTACA 1 5
BstBI TTCGAA 2 4
BstEII GGTNACC 1 6
BstXI CCANNNNNNTGG 8 4
ClaI ATCGAT 2 4
DpnII GATC 0 4
DraI TTTAAA 3 3
EagI CGGCCG 1 5
EcoNI CCTNNNNNAGG 5 6
EcoRI GAATTC 1 5
EcoRV GATATC 3 3
FseI GGCCGGCC 6 2
HaeIII GGCC 2 2
HhaI GCGC 3 1
HincII GTYRAC 3 3
HindIII AAGCTT 1 5
HpaI GTTAAC 3 3
KasI GGCGCC 1 5
KpnI GGTACC 5 1
MfeI CAATTG 1 5
MluI ACGCGT 1 5
MspI CCGG 1 3
NaeI GCCGGC 3 3
NcoI CCATGG 1 5
NdeI CATATG 2 4
NheI GCTAGC 1 5
NotI GCGGCCGC 2 6
NruI TCGCGA 3 3
NsiI ATGCAT 5 1
PacI TTAATTAA 5 3
PciI ACATGT 1 5
PmeI GTTTAAAC 4 4
PstI CTGCAG 5 1
PvuI CGATCG 4 2
PvuII CAGCTG 3 3
RsaI GTAC 2 2
SacI GAGCTC 5 1
SacII CCGCGG 4 2
SalI GTCGAC 1 5
SapI GCTCTTC 8 11
SbfI CCTGCAGG 6 2
ScaI AGTACT 3 3
SmaI CCCGGG 3 3
SnaBI TACGTA 3 3
SpeI ACTAGT 1 5
SphI GCATGC 5 1
SspI AATATT 3 3
StuI AGGCCT 3 3
SwaI ATTTAAAT 4 4
TaqI TCGA 1 3
XbaI TCTAGA 1 5
XhoI CTCGAG 1 5
XmaI CCCGGG 1 5
XmnI GAANNNNTTC 5 5
ZraI GACGTC 3 3`;
  const IUPAC_RE = { A: 'A', C: 'C', G: 'G', T: 'T', R: '[AG]', Y: '[CT]', S: '[GC]', W: '[AT]', K: '[GT]', M: '[AC]', B: '[CGT]', D: '[AGT]', H: '[ACT]', V: '[ACG]', N: '[ACGT]' };
  function iupacPattern(site) { return site.split('').map(c => IUPAC_RE[c] || c).join(''); }
  const ENZYMES = ENZ_RAW.split('\n').map(line => {
    const [name, site, top, bot] = line.trim().split(/\s+/);
    const pal = rc(site) === site;
    return {
      name, site, top: +top, bot: +bot, pal,
      core: site.replace(/N/g, '').length,
      re: new RegExp('(?=' + iupacPattern(site) + ')', 'g'),
      rcre: pal ? null : new RegExp('(?=' + iupacPattern(rc(site)) + ')', 'g'),
      overhang: (+bot) - (+top) // >0 5' overhang, <0 3' overhang, 0 blunt
    };
  });
  const ENZ_BY_NAME = Object.fromEntries(ENZYMES.map(e => [e.name, e]));

  function scan(re, text, fn) { re.lastIndex = 0; let m; while ((m = re.exec(text))) { fn(m.index); re.lastIndex = m.index + 1; } }

  function findSites(seq, circular) {
    const n = seq.length, out = {};
    for (const e of ENZYMES) {
      const L = e.site.length;
      if (n < L) continue;
      const text = circular ? seq + seq.slice(0, L - 1) : seq;
      const res = [];
      const norm = c => {
        if (circular) return ((c % n) + n) % n;
        return (c <= 0 || c >= n) ? null : c;
      };
      scan(e.re, text, p => { if (p < n) res.push({ pos: p, len: L, strand: 1, cut: norm(p + e.top), cut2: norm(p + e.bot) }); });
      if (e.rcre) scan(e.rcre, text, p => { if (p < n) res.push({ pos: p, len: L, strand: -1, cut: norm(p + L - e.bot), cut2: norm(p + L - e.top) }); });
      if (res.length) out[e.name] = res.sort((a, b) => a.pos - b.pos);
    }
    return out;
  }

  // ---------- Ranges (0-based, end exclusive; end <= start means it wraps the origin) ----------
  function rlen(r, n) { return r.end > r.start ? r.end - r.start : (r.end === r.start ? 0 : r.end + n - r.start); }
  function segs(r, n) { return r.end > r.start ? [[r.start, r.end]] : (r.end === r.start ? [] : [[r.start, n], [0, r.end]]); }
  function sub(seq, start, end) { return end > start ? seq.slice(start, end) : (end === start ? '' : seq.slice(start) + seq.slice(0, end)); }
  function overlaps(a, b, n) {
    for (const [a1, a2] of segs(a, n)) for (const [b1, b2] of segs(b, n)) if (a1 < b2 && b1 < a2) return true;
    return false;
  }

  // ---------- ORFs ----------
  function findORFs(seq, circular, minAA) {
    const n = seq.length; if (n < 6) return [];
    const s = circular ? seq + seq : seq, m = s.length;
    const seen = new Set(), out = [];
    for (const strand of [1, -1]) {
      const t = strand === 1 ? s : rc(s);
      for (let f = 0; f < 3; f++) {
        let st = -1;
        for (let i = f; i + 3 <= m; i += 3) {
          const c = t.substr(i, 3);
          if (st < 0 && c === 'ATG') st = i;
          else if (st >= 0 && (c === 'TAA' || c === 'TAG' || c === 'TGA')) {
            const aa = (i - st) / 3;
            const len = i + 3 - st;
            if (aa >= minAA && len <= n) {
              let a, b;
              if (strand === 1) { a = st; b = i + 3; } else { a = m - (i + 3); b = m - st; }
              const start = a % n, end = ((b - 1) % n) + 1;
              const key = start + ':' + end + ':' + strand;
              if (!seen.has(key) && (circular || b <= n)) { seen.add(key); out.push({ start, end, strand, aa }); }
            }
            st = -1;
          }
        }
      }
    }
    // remove ORFs nested inside a longer ORF on the same strand that shares the stop codon
    out.sort((x, y) => y.aa - x.aa);
    const keep = [];
    for (const o of out) {
      const stop = o.strand === 1 ? o.end : o.start;
      if (keep.some(k => k.strand === o.strand && (k.strand === 1 ? k.end : k.start) === stop)) continue;
      keep.push(o);
    }
    return keep.sort((x, y) => x.start - y.start);
  }

  // ---------- Search (IUPAC DNA or protein) ----------
  function searchDNA(seq, circular, q) {
    q = clean(q); const n = seq.length; if (!q || q.length > n) return [];
    const text = circular ? seq + seq.slice(0, q.length - 1) : seq;
    const res = [], seen = new Set();
    const add = (p, strand) => { if (p >= n) return; const k = p + ':' + strand; if (seen.has(k)) return; seen.add(k); let e = p + q.length; if (e > n) e -= n; res.push({ start: p, end: e, strand }); };
    scan(new RegExp('(?=' + iupacPattern(q) + ')', 'g'), text, p => add(p, 1));
    const r = rc(q); if (r !== q) scan(new RegExp('(?=' + iupacPattern(r) + ')', 'g'), text, p => add(p, -1));
    return res.sort((a, b) => a.start - b.start);
  }
  function searchProtein(seq, circular, q) {
    q = String(q).toUpperCase().replace(/[^A-Z*]/g, ''); const n = seq.length; if (!q) return [];
    const L = q.length * 3; if (L > n) return [];
    const ext = circular ? seq + seq.slice(0, L - 1) : seq, m = ext.length;
    const res = [], seen = new Set();
    const add = (p, strand) => { if (p < 0 || p >= n) return; const k = p + ':' + strand; if (seen.has(k)) return; seen.add(k); let e = p + L; if (e > n) e -= n; res.push({ start: p, end: e, strand }); };
    const rcx = rc(ext);
    for (let f = 0; f < 3; f++) {
      const pf = translate(ext.slice(f)); let i = pf.indexOf(q);
      while (i >= 0) { add(f + 3 * i, 1); i = pf.indexOf(q, i + 1); }
      const pr = translate(rcx.slice(f)); i = pr.indexOf(q);
      while (i >= 0) { add(m - (f + 3 * i + L), -1); i = pr.indexOf(q, i + 1); }
    }
    return res.sort((a, b) => a.start - b.start);
  }

  // ---------- Digest ----------
  function digest(seq, circular, sites, names) {
    const n = seq.length; const cuts = new Map();
    for (const nm of names) for (const s of (sites[nm] || [])) if (s.cut != null) {
      if (!cuts.has(s.cut)) cuts.set(s.cut, []);
      cuts.get(s.cut).push(nm);
    }
    const pos = [...cuts.keys()].sort((a, b) => a - b);
    const frags = [];
    if (circular) {
      if (!pos.length) return [{ size: n, start: 0, end: 0, left: '', right: '', uncut: true }];
      for (let i = 0; i < pos.length; i++) {
        const a = pos[i], b = pos[(i + 1) % pos.length];
        frags.push({ start: a, end: b, size: ((b - a + n) % n) || n, left: cuts.get(a).join('/'), right: cuts.get(b).join('/') });
      }
    } else {
      const all = [0, ...pos, n];
      for (let i = 0; i < all.length - 1; i++) frags.push({ start: all[i], end: all[i + 1], size: all[i + 1] - all[i], left: i === 0 ? 'end' : cuts.get(all[i]).join('/'), right: i === all.length - 2 ? 'end' : cuts.get(all[i + 1]).join('/') });
    }
    return frags.sort((a, b) => b.size - a.size);
  }

  // ---------- Common feature library (for "Detect Common Features") ----------
  const COMMON = [
    ['T7 promoter', 'promoter', 'TAATACGACTCACTATAGG'],
    ['SP6 promoter', 'promoter', 'ATTTAGGTGACACTATAG'],
    ['T3 promoter', 'promoter', 'AATTAACCCTCACTAAAGG'],
    ['lac operator', 'protein_bind', 'TTGTGAGCGGATAACAA'],
    ['lac promoter', 'promoter', 'TTTACACTTTATGCTTCCGGCTCGTATGTTG'],
    ['T7 terminator', 'terminator', 'CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG'],
    ['loxP', 'protein_bind', 'ATAACTTCGTATAGCATACATTATACGAAGTTAT'],
    ['lox2272', 'protein_bind', 'ATAACTTCGTATAAAGTATCCTATACGAAGTTAT'],
    ['FRT', 'protein_bind', 'GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC'],
    ['attB1', 'protein_bind', 'ACAAGTTTGTACAAAAAAGCAGGCT'],
    ['attB2', 'protein_bind', 'ACCCAGCTTTCTTGTACAAAGTGGT'],
    ['6xHis', 'CDS', 'CATCATCATCATCATCAT'],
    ['6xHis', 'CDS', 'CACCACCACCACCACCAC'],
    ['FLAG', 'CDS', 'GATTACAAGGATGACGACGATAAG'],
    ['HA', 'CDS', 'TACCCATACGATGTTCCAGATTACGCT'],
    ['Myc', 'CDS', 'GAACAAAAACTCATCTCAGAAGAGGATCTG'],
    ['V5 tag', 'CDS', 'GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG'],
    ['SV40 NLS', 'CDS', 'CCAAAGAAGAAGCGGAAGGTC'],
    ['TEV site', 'CDS', 'GAGAATCTTTATTTTCAGGGC'],
    ['thrombin site', 'CDS', 'CTGGTTCCGCGTGGATCC'],
    ['Kozak sequence', 'regulatory', 'GCCACCATGG'],
    ['M13 fwd', 'primer_bind', 'GTAAAACGACGGCCAGT'],
    ['M13 rev', 'primer_bind', 'CAGGAAACAGCTATGAC'],
    ['bGH poly(A) signal', 'polyA_signal', 'CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCG']
  ];
  function detectCommon(seq, circular, existing) {
    const n = seq.length, found = [];
    for (const [name, type, s] of COMMON) {
      for (const h of searchDNA(seq, circular, s)) {
        if (existing.some(f => f.start === h.start && f.end === h.end)) continue;
        if (found.some(f => f.start === h.start && f.end === h.end)) continue;
        found.push({ name, type, start: h.start, end: h.end, strand: h.strand });
      }
    }
    return found;
  }

  // ---------- File formats ----------
  function parseLocation(loc, n, circular) {
    const strand = /complement/.test(loc) ? -1 : 1;
    const ranges = [];
    const re = /<?(\d+)(?:\.\.>?(\d+)|\^(\d+))?/g; let m;
    while ((m = re.exec(loc))) { const a = +m[1], b = m[2] ? +m[2] : (m[3] ? +m[3] : +m[1]); ranges.push([a, b]); }
    if (!ranges.length) return null;
    if (strand === -1 && /complement\(join/.test(loc) === false && /join\(complement/.test(loc)) ranges.reverse();
    const first = ranges[0], last = ranges[ranges.length - 1];
    let start, end;
    if (circular && ranges.length > 1 && /join|order/.test(loc)) {
      const sorted = ranges.slice();
      const wraps = sorted.some((r, i) => i > 0 && r[0] < sorted[i - 1][1]);
      if (wraps) {
        const ordered = strand === -1 && /join\(complement/.test(loc) ? ranges.slice().reverse() : ranges;
        start = ordered[0][0] - 1; end = ordered[ordered.length - 1][1];
        return { start, end: end % n === 0 ? n : end % n, strand };
      }
    }
    start = Math.min(...ranges.map(r => r[0])) - 1;
    end = Math.max(...ranges.map(r => r[1]));
    if (start > end - 1 && circular) { start = first[0] - 1; end = last[1]; }
    return { start, end, strand };
  }

  function parseGenBank(txt) {
    const lines = txt.split(/\r?\n/);
    let name = 'Untitled', circular = false, mode = '', seq = '', cur = null, lastQ = null, definition = '';
    const raw = [];
    for (const ln of lines) {
      if (ln.startsWith('LOCUS')) { const t = ln.split(/\s+/); name = t[1] || name; circular = /\bcircular\b/i.test(ln); continue; }
      if (ln.startsWith('DEFINITION')) { definition = ln.slice(12).trim(); continue; }
      if (ln.startsWith('FEATURES')) { mode = 'F'; continue; }
      if (ln.startsWith('ORIGIN')) { mode = 'O'; continue; }
      if (ln.startsWith('//')) break;
      if (mode === 'F') {
        if (/^ {5}\S/.test(ln)) { cur = { type: ln.slice(5, 21).trim(), loc: ln.slice(21).trim(), q: {} }; raw.push(cur); lastQ = null; }
        else if (/^ {21}\//.test(ln) && cur) {
          const m = ln.trim().match(/^\/([^=]+)(?:=(.*))?$/);
          if (m) { lastQ = m[1]; cur.q[lastQ] = (m[2] || '').replace(/^"/, ''); }
        } else if (/^ {21}\S/.test(ln) && cur) {
          if (lastQ) cur.q[lastQ] += (lastQ === 'translation' ? '' : ' ') + ln.trim();
          else cur.loc += ln.trim();
        }
      } else if (mode === 'O') seq += ln.replace(/[^a-zA-Z]/g, '');
      else if (/^\s+\S/.test(ln) && !mode && definition) { /* continued definition */ }
    }
    seq = clean(seq);
    const n = seq.length;
    const features = [], primers = [];
    for (const f of raw) {
      for (const k in f.q) f.q[k] = f.q[k].replace(/"$/, '').trim();
      if (f.type === 'source') continue;
      const loc = parseLocation(f.loc, n, circular); if (!loc) continue;
      const fname = f.q.label || f.q.gene || f.q.product || f.q.standard_name || f.q.note || f.type;
      let color = f.q.ApEinfo_fwdcolor || null;
      const cm = (f.q.note || '').match(/color:\s*(#[0-9a-fA-F]{6})/); if (cm) color = cm[1];
      if (color && !/^#[0-9a-fA-F]{6}$/.test(color)) color = null;
      features.push({ name: fname.slice(0, 80), type: f.type, start: loc.start, end: loc.end, strand: loc.strand, color, notes: f.q.note || '' });
    }
    return { name, seq, circular, features, primers, definition };
  }

  function parseFASTA(txt) {
    const lines = txt.split(/\r?\n/);
    let name = 'Untitled', seq = '', got = false;
    for (const l of lines) {
      if (l.startsWith('>')) { if (got) break; got = true; name = l.slice(1).trim().split(/\s+/)[0] || name; }
      else if (got) seq += l;
    }
    return { name, seq: clean(seq), circular: false, features: [], primers: [] };
  }

  function parseAny(txt, filename) {
    const t = txt.replace(/^﻿/, '').trimStart();
    const base = (filename || '').replace(/\.[^.]+$/, '');
    let r;
    if (t.startsWith('LOCUS')) r = parseGenBank(t);
    else if (t.startsWith('>')) r = parseFASTA(t);
    else r = { name: base || 'Untitled', seq: clean(t), circular: false, features: [], primers: [] };
    if (base && (r.name === 'Untitled' || r.name === 'Exported')) r.name = base;
    if (!r.seq.length) throw new Error('No DNA sequence found in this file.');
    return r;
  }

  // .dna binary (block stream: type byte, uint32 BE length, payload)
  function parseDnaBinary(buf, filename) {
    const dv = new DataView(buf); let o = 0;
    let seq = '', circular = false; const features = [], primers = [];
    const dec = new TextDecoder();
    let sawCookie = false;
    while (o + 5 <= buf.byteLength) {
      const type = dv.getUint8(o), len = dv.getUint32(o + 1), start = o + 5;
      o = start + len;
      if (o > buf.byteLength) break;
      if (type === 9) sawCookie = true;
      if (type === 0) {
        circular = !!(dv.getUint8(start) & 1);
        seq = dec.decode(new Uint8Array(buf, start + 1, len - 1));
      } else if (type === 10 || type === 5) {
        const xml = dec.decode(new Uint8Array(buf, start, len));
        let doc;
        try { doc = new DOMParser().parseFromString(xml, 'application/xml'); } catch (e) { continue; }
        if (type === 10) {
          for (const fe of doc.querySelectorAll('Feature')) {
            const segsEl = [...fe.querySelectorAll('Segment')];
            if (!segsEl.length) continue;
            const ranges = segsEl.map(s => (s.getAttribute('range') || '').split('-').map(Number)).filter(r => r.length === 2 && !isNaN(r[0]));
            if (!ranges.length) continue;
            const dir = +(fe.getAttribute('directionality') || 0);
            let qname = fe.getAttribute('name') || fe.getAttribute('type') || 'feature';
            let st = ranges[0][0] - 1, en = ranges[ranges.length - 1][1];
            if (ranges.length > 1 && dir === 2) { st = Math.min(...ranges.map(r => r[0])) - 1; en = Math.max(...ranges.map(r => r[1])); }
            features.push({ name: qname, type: fe.getAttribute('type') || 'misc_feature', start: st, end: en, strand: dir === 2 ? -1 : (dir === 1 ? 1 : 0), color: segsEl[0].getAttribute('color') || null });
          }
        } else {
          for (const p of doc.querySelectorAll('Primer')) {
            const s = p.getAttribute('sequence'); if (!s) continue;
            primers.push({ name: p.getAttribute('name') || 'Primer', seq: clean(s) });
          }
        }
      }
    }
    if (!sawCookie && !seq) throw new Error('This is not a .dna file.');
    seq = clean(seq);
    if (!seq) throw new Error('The .dna file has no DNA sequence block.');
    const n = seq.length;
    for (const f of features) { if (f.end > n) f.end = n; if (f.start >= n) f.start = 0; if (f.color && !/^#[0-9a-fA-F]{6}$/.test(f.color)) f.color = null; }
    return { name: (filename || 'Untitled').replace(/\.[^.]+$/, ''), seq, circular, features, primers };
  }

  function gbLoc(f, n) {
    let s;
    if (f.end > f.start) s = (f.end - f.start === 1) ? `${f.start + 1}` : `${f.start + 1}..${f.end}`;
    else s = `join(${f.start + 1}..${n},1..${f.end})`;
    return f.strand === -1 ? `complement(${s})` : s;
  }
  function wrapQual(key, val) {
    const full = `/${key}="${String(val).replace(/"/g, "'")}"`;
    const out = []; let rest = full;
    while (rest.length > 58) { out.push(rest.slice(0, 58)); rest = rest.slice(58); }
    out.push(rest);
    return out.map(l => ' '.repeat(21) + l).join('\n');
  }
  function toGenBank(d, primerSitesFn) {
    const n = d.seq.length;
    const months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'];
    const now = new Date();
    const date = `${String(now.getDate()).padStart(2, '0')}-${months[now.getMonth()]}-${now.getFullYear()}`;
    const nm = d.name.replace(/\s+/g, '_').slice(0, 16);
    let s = `LOCUS       ${nm.padEnd(16)} ${String(n).padStart(11)} bp    DNA     ${d.circular ? 'circular' : 'linear  '} SYN ${date}\n`;
    s += `DEFINITION  ${d.name}.\nACCESSION   .\nVERSION     .\nKEYWORDS    GeneCanvas\nSOURCE      synthetic DNA construct\n  ORGANISM  synthetic DNA construct\n`;
    s += `FEATURES             Location/Qualifiers\n`;
    s += `     source          1..${n}\n${wrapQual('mol_type', 'other DNA')}\n`;
    const feats = [...d.features];
    if (primerSitesFn) for (const p of d.primers) for (const b of primerSitesFn(p)) feats.push({ type: 'primer_bind', name: p.name, start: b.start, end: b.end, strand: b.strand, notes: 'sequence: ' + p.seq });
    feats.sort((a, b) => a.start - b.start);
    for (const f of feats) {
      s += `     ${f.type.slice(0, 15).padEnd(16)}${gbLoc(f, n)}\n`;
      s += wrapQual('label', f.name) + '\n';
      if (f.notes) s += wrapQual('note', f.notes) + '\n';
      if (f.color) s += wrapQual('ApEinfo_fwdcolor', f.color) + '\n' + wrapQual('ApEinfo_revcolor', f.color) + '\n';
    }
    s += 'ORIGIN\n';
    const lc = d.seq.toLowerCase();
    for (let i = 0; i < n; i += 60) {
      let line = String(i + 1).padStart(9);
      for (let j = i; j < Math.min(i + 60, n); j += 10) line += ' ' + lc.slice(j, Math.min(j + 10, n));
      s += line + '\n';
    }
    return s + '//\n';
  }
  function toFASTA(name, seq) {
    let s = '>' + name.replace(/\s+/g, '_') + '\n';
    for (let i = 0; i < seq.length; i += 70) s += seq.slice(i, i + 70) + '\n';
    return s;
  }

  // ---------- Demo plasmid (synthetic sequence, deterministic) ----------
  function buildDemo() {
    let a = 20231117 >>> 0;
    const rnd = () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
    const rs = k => { let s = ''; for (let i = 0; i < k; i++) s += 'ACGT'[Math.floor(rnd() * 4)]; return s; };
    const sense = Object.keys(CODON).filter(c => CODON[c] !== '*');
    const rcod = k => { let s = ''; for (let i = 0; i < k; i++) s += sense[Math.floor(rnd() * sense.length)]; return s; };
    let seq = ''; const features = [];
    const add = (s, f) => { const st = seq.length; seq += s; if (f) features.push(Object.assign({ start: st, end: seq.length }, f)); };
    add(rs(60));
    add('TTTACACTTTATGCTTCCGGCTCGTATGTTG', { name: 'lac promoter', type: 'promoter', strand: 1 });
    add('TG');
    add('TGGAATTGTGAGCGGATAACAATT', { name: 'lac operator', type: 'protein_bind', strand: 1 });
    add('TCACA');
    add('CAGGAAACAGCT');
    const lacZ = 'ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTAA';
    const lzStart = seq.length;
    add(lacZ, { name: 'lacZα', type: 'CDS', strand: 1 });
    const mcsS = lzStart + lacZ.indexOf('AAGCTT'), mcsE = lzStart + lacZ.indexOf('GAATTC') + 6;
    features.push({ name: 'MCS', type: 'misc_feature', start: mcsS, end: mcsE, strand: 0 });
    add(rs(140));
    add(rs(87), { name: 'rrnB T1 terminator', type: 'terminator', strand: 1 });
    add(rs(160));
    add(rs(589), { name: 'ori', type: 'rep_origin', strand: -1 });
    add(rs(140));
    const ampCDS = 'ATG' + rcod(285) + 'TAA';
    const ampS = seq.length;
    add(rc(ampCDS), { name: 'AmpR', type: 'CDS', strand: -1 });
    add(rs(105), { name: 'AmpR promoter', type: 'promoter', strand: -1 });
    add(rs(270));
    // keep MCS enzymes unique
    const protect = [[lzStart - 60, lzStart + lacZ.length]];
    const ampOK = s => { const p = translate(rc(s.slice(ampS, ampS + ampCDS.length))); return p.indexOf('*') === p.length - 1; };
    const uniq = ['HindIII', 'SphI', 'PstI', 'SbfI', 'SalI', 'XbaI', 'BamHI', 'SmaI', 'XmaI', 'KpnI', 'SacI', 'EcoRI'];
    for (let iter = 0; iter < 400; iter++) {
      const sites = findSites(seq, true);
      let fixed = false;
      for (const nm of uniq) {
        const list = (sites[nm] || []).filter(s => !(s.pos >= mcsS - 2 && s.pos < mcsE));
        if (!list.length) continue;
        const s = list[0];
        for (const off of [2, 3, 1, 4]) {
          const q = s.pos + off;
          if (protect.some(([x, y]) => q >= x && q < y)) continue;
          let done = false;
          for (const b of 'ACGT') {
            if (b === seq[q]) continue;
            const t = seq.slice(0, q) + b + seq.slice(q + 1);
            if (!ampOK(t)) continue;
            seq = t; done = true; break;
          }
          if (done) { fixed = true; break; }
        }
        if (fixed) break;
      }
      if (!fixed) break;
    }
    const primers = [
      { name: 'M13 fwd', seq: 'GTAAAACGACGGCCAGT' },
      { name: 'M13 rev', seq: 'CAGGAAACAGCTATGAC' },
      { name: 'AmpR-R', seq: rc(seq.slice(ampS + 300, ampS + 321)) }
    ];
    return { name: 'pGC-Demo', seq, circular: true, features, primers, demo: true };
  }

  root.DNA = {
    COMP, clean, rc, comp, CODON, translate, AA3, proteinMW, gcPct, dsDnaMW, tm,
    ENZYMES, ENZ_BY_NAME, findSites, rlen, segs, sub, overlaps, findORFs, searchDNA, searchProtein,
    digest, detectCommon, COMMON, parseGenBank, parseFASTA, parseAny, parseDnaBinary, toGenBank, toFASTA, buildDemo
  };
})(typeof window !== 'undefined' ? window : globalThis);
