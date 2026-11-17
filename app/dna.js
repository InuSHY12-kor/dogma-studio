/* Dogma Studio — DNA logic (no DOM). Exposes window.DNA */
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

  /* ---------- Nucleotide masses and composition ----------
   * Average monoisotopic-free (average) masses of the monophosphate nucleotides, less
   * the water lost at each phosphodiester bond. An ambiguity code is counted as the
   * mean of the bases it stands for, which is the usual convention.
   * A 5'-phosphate adds HPO3 (79.98); an unphosphorylated 5' end has none. */
  const NT_MW = { A: 331.22, C: 307.20, G: 347.22, T: 322.21, U: 324.18 };
  const H2O = 18.02, HPO3 = 79.98;
  const IUPAC_SETS = { R: 'AG', Y: 'CT', S: 'GC', W: 'AT', K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT' };
  function ntMW(ch, rna) {
    if (NT_MW[ch] !== undefined) return rna && ch === 'T' ? NT_MW.U : NT_MW[ch];
    const set = IUPAC_SETS[ch];
    if (!set) return 0;
    let t = 0;
    for (const c of set) t += rna && c === 'T' ? NT_MW.U : NT_MW[c];
    return t / set.length;
  }
  /* Molecular weight of a strand.
   *   ends: 'phosphorylated' | 'unphosphorylated' | 'circular'
   * A circular molecule loses one more water than a linear one, and carries no free end. */
  function strandMW(seq, { rna = false, ends = 'unphosphorylated' } = {}) {
    let t = 0;
    for (const ch of seq) t += ntMW(ch, rna);
    // each residue above is a monophosphate; joining n of them loses (n-1) waters
    t -= H2O * Math.max(0, seq.length - 1);
    if (ends === 'circular') t -= H2O;
    else if (ends === 'unphosphorylated') t -= HPO3;      // the 5' phosphate is not there
    return t;
  }
  function composition(seq) {
    const out = { A: 0, C: 0, G: 0, T: 0, other: 0 };
    for (const ch of seq) { if (out[ch] !== undefined) out[ch]++; else out.other++; }
    return out;
  }
  /* %GC over a sliding window, which is what the GC content plot draws.
   * A degenerate position counts as the fraction of its bases that are G or C. */
  function gcCurve(seq, win = 25, circular = false) {
    const n = seq.length;
    if (!n) return [];
    const frac = ch => {
      if (ch === 'G' || ch === 'C' || ch === 'S') return 1;
      if (ch === 'A' || ch === 'T' || ch === 'U' || ch === 'W') return 0;
      const set = IUPAC_SETS[ch];
      if (!set) return 0;
      let g = 0; for (const c of set) if (c === 'G' || c === 'C') g++;
      return g / set.length;
    };
    const w = Math.max(3, Math.min(win, n));
    const half = Math.floor(w / 2);
    const out = new Array(n);
    let sum = 0;
    const at = i => frac(seq[circular ? ((i % n) + n) % n : Math.max(0, Math.min(n - 1, i))]);
    for (let i = -half; i < w - half; i++) sum += at(i);
    for (let i = 0; i < n; i++) {
      out[i] = sum / w * 100;
      sum += at(i + w - half) - at(i - half);
    }
    return out;
  }

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
AciI CCGC 1 3
AclI AACGTT 2 4
AflII CTTAAG 1 5
AflIII ACRYGT 1 5
AgeI ACCGGT 1 5
AhdI GACNNNNNGTC 6 5
AleI CACNNNNGTG 5 5
AluI AGCT 2 2
AlwNI CAGNNNCTG 6 3
ApaI GGGCCC 5 1
ApaLI GTGCAC 1 5
ApeKI GCWGC 1 4
AscI GGCGCGCC 2 6
AseI ATTAAT 2 4
AvaI CYCGRG 1 5
AvaII GGWCC 1 4
AvrII CCTAGG 1 5
BamHI GGATCC 1 5
BanI GGYRCC 1 5
BanII GRGCYC 5 1
BbsI GAAGAC 8 12
BbvCI CCTCAGC 2 5
BccI CCATC 9 10
BciVI GTATCC 18 16
BclI TGATCA 1 5
BfaI CTAG 1 3
BglI GCCNNNNNGGC 7 5
BglII AGATCT 1 5
BlpI GCTNAGC 2 5
BmgBI CACGTC 3 3
BmrI ACTGGG 11 10
BmtI GCTAGC 5 1
BpmI CTGGAG 22 20
Bpu10I CCTNAGC 2 5
BpuEI CTTGAG 22 20
BsaAI YACGTR 3 3
BsaBI GATNNNNATC 5 5
BsaHI GRCGYC 2 4
BsaI GGTCTC 7 11
BsaJI CCNNGG 1 5
BsaWI WCCGGW 1 5
BseRI GAGGAG 16 14
BseYI CCCAGC 1 5
BsgI GTGCAG 22 20
BsiEI CGRYCG 4 2
BsiHKAI GWGCWC 5 1
BsiWI CGTACG 1 5
BslI CCNNNNNNNGG 7 5
BsmAI GTCTC 6 10
BsmBI CGTCTC 7 11
BsmFI GGGAC 15 19
BsmI GAATGC 7 5
BsoBI CYCGRG 1 5
Bsp1286I GDGCHC 5 1
BspCNI CTCAG 14 12
BspDI ATCGAT 2 4
BspEI TCCGGA 1 5
BspHI TCATGA 1 5
BspMI ACCTGC 10 14
BspQI GCTCTTC 8 11
BsrBI CCGCTC 3 3
BsrDI GCAATG 8 6
BsrFI RCCGGY 1 5
BsrGI TGTACA 1 5
BsrI ACTGG 6 4
BssHII GCGCGC 1 5
BssSI CACGAG 1 5
BstAPI GCANNNNNTGC 7 5
BstBI TTCGAA 2 4
BstEII GGTNACC 1 6
BstNI CCWGG 2 3
BstUI CGCG 2 2
BstXI CCANNNNNNTGG 8 4
BstYI RGATCY 1 5
BstZ17I GTATAC 3 3
Bsu36I CCTNAGG 2 5
BtgI CCRYGG 1 5
BtgZI GCGATG 16 20
BtsI GCAGTG 8 6
Cac8I GCNNGC 3 3
ClaI ATCGAT 2 4
CviAII CATG 1 3
CviKI-1 RGCY 2 2
CviQI GTAC 1 3
DdeI CTNAG 1 4
DpnII GATC 0 4
DraI TTTAAA 3 3
DraIII CACNNNGTG 6 3
DrdI GACNNNNNNGTC 7 5
EaeI YGGCCR 1 5
EagI CGGCCG 1 5
EarI CTCTTC 7 10
EciI GGCGGA 17 15
Eco53kI GAGCTC 3 3
EcoNI CCTNNNNNAGG 5 6
EcoO109I RGGNCCY 2 5
EcoRI GAATTC 1 5
EcoRV GATATC 3 3
FatI CATG 0 4
FauI CCCGC 9 11
Fnu4HI GCNGC 2 3
FokI GGATG 14 18
FseI GGCCGGCC 6 2
FspI TGCGCA 3 3
HaeII RGCGCY 5 1
HaeIII GGCC 2 2
HgaI GACGC 10 15
HhaI GCGC 3 1
HincII GTYRAC 3 3
HindIII AAGCTT 1 5
HinfI GANTC 1 4
HinP1I GCGC 1 3
HpaI GTTAAC 3 3
HpaII CCGG 1 3
HphI GGTGA 13 12
Hpy166II GTNNAC 3 3
Hpy188I TCNGA 3 2
HpyAV CCTTC 11 10
HpyCH4III ACNGT 3 2
HpyCH4IV ACGT 1 3
HpyCH4V TGCA 2 2
KasI GGCGCC 1 5
KpnI GGTACC 5 1
MboI GATC 0 4
MboII GAAGA 13 12
MfeI CAATTG 1 5
MluCI AATT 0 4
MluI ACGCGT 1 5
MlyI GAGTC 10 10
MmeI TCCRAC 26 24
MnlI CCTC 11 10
MscI TGGCCA 3 3
MseI TTAA 1 3
MslI CAYNNNNRTG 5 5
MspA1I CMGCKG 3 3
MspI CCGG 1 3
MwoI GCNNNNNNNGC 7 5
NaeI GCCGGC 3 3
NciI CCSGG 2 3
NcoI CCATGG 1 5
NdeI CATATG 2 4
NgoMIV GCCGGC 1 5
NheI GCTAGC 1 5
NlaIII CATG 4 0
NlaIV GGNNCC 3 3
NmeAIII GCCGAG 27 25
NotI GCGGCCGC 2 6
NruI TCGCGA 3 3
NsiI ATGCAT 5 1
NspI RCATGY 5 1
PacI TTAATTAA 5 3
PaeR7I CTCGAG 1 5
PciI ACATGT 1 5
PflFI GACNNNGTC 4 5
PflMI CCANNNNNTGG 7 5
PleI GAGTC 9 10
PluTI GGCGCC 5 1
PmeI GTTTAAAC 4 4
PpuMI RGGWCCY 2 5
PshAI GACNNNNGTC 5 5
PspGI CCWGG 0 5
PspOMI GGGCCC 1 5
PspXI VCTCGAGB 2 6
PstI CTGCAG 5 1
PvuI CGATCG 4 2
PvuII CAGCTG 3 3
RsaI GTAC 2 2
RsrII CGGWCCG 2 5
SacI GAGCTC 5 1
SacII CCGCGG 4 2
SalI GTCGAC 1 5
SapI GCTCTTC 8 11
SbfI CCTGCAGG 6 2
ScaI AGTACT 3 3
SexAI ACCWGGT 1 6
SfaNI GCATC 10 14
SfcI CTRYAG 1 5
SfiI GGCCNNNNNGGCC 8 5
SgrAI CRCCGGYG 2 6
SmaI CCCGGG 3 3
SmlI CTYRAG 1 5
SnaBI TACGTA 3 3
SpeI ACTAGT 1 5
SphI GCATGC 5 1
SrfI GCCCGGGC 4 4
SspI AATATT 3 3
StuI AGGCCT 3 3
StyD4I CCNGG 0 5
StyI CCWWGG 1 5
SwaI ATTTAAAT 4 4
TaqI TCGA 1 3
TfiI GAWTC 1 4
Tsp45I GTSAC 0 5
TspMI CCCGGG 1 5
Tth111I GACNNNGTC 4 5
XbaI TCTAGA 1 5
XcmI CCANNNNNNNNNTGG 8 7
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

  /* Suppliers sell the same enzyme under their own names — FastDigest® EcoRI, EcoRI-HF®,
   * Anza™ EcoRI — and those are the names on the tube in the freezer. They are not added
   * as separate enzymes, which would draw the same site twice on every map; they hang off
   * the enzyme they are, so searching for one finds it. A branded name is only accepted
   * when the supplier lists the same recognition sequence, so a coincidence of naming
   * cannot give an enzyme a cut position that is not its own. */
  const BRANDS = /^(FastDigest®?|Anza™|QuickCut™|EZ-CleanCut™|Fast\s*Digest|FuniCut™|SureCut™|WarmStart®?)\s+/i;
  function baseEnzymeName(name) {
    return String(name || '').trim()
      .replace(BRANDS, '')
      .replace(/\s*[-\u2013]?\s*HF®?\b/i, '')
      .replace(/\s*-?v\d+$/i, '')
      .replace(/[®™]/g, '')
      .trim();
  }
  (function supplierInfo() {
    const db = root.ENZYME_DB;
    if (!db) return;
    for (const s of db.suppliers || []) {
      for (const e of s.enzymes || []) {
        const canon = ENZ_BY_NAME[baseEnzymeName(e.name)];
        if (!canon || canon.site !== e.site) continue;
        if (!canon.suppliers) { canon.suppliers = []; canon.aliases = []; }
        if (!canon.suppliers.includes(s.name)) canon.suppliers.push(s.name);
        if (e.name !== canon.name && !canon.aliases.includes(e.name)) canon.aliases.push(e.name);
      }
    }
    for (const e of ENZYMES) {
      if (e.suppliers) { e.suppliers.sort(); e.aliases.sort(); }
      else { e.suppliers = []; e.aliases = []; }
    }
  })();

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
  function findORFs(seq, circular, minAA, opts) {
    opts = opts || {};
    const table = GENETIC_CODES[opts.code] || CODON;
    const starts = opts.starts && opts.starts.length ? opts.starts : ['ATG'];
    const anyStart = !!opts.anyStart;      // start at any sense codon after a stop
    const strands = opts.strands === 'top' ? [1] : (opts.strands === 'bottom' ? [-1] : [1, -1]);
    const n = seq.length; if (n < 6) return [];
    const s = circular ? seq + seq : seq, m = s.length;
    const seen = new Set(), out = [];
    for (const strand of strands) {
      const t = strand === 1 ? s : rc(s);
      for (let f = 0; f < 3; f++) {
        let st = -1;
        for (let i = f; i + 3 <= m; i += 3) {
          const c = t.substr(i, 3);
          const isStop = (table[c] || 'X') === '*';
          if (st < 0 && !isStop && (anyStart || starts.includes(c))) st = i;
          else if (st >= 0 && isStop) {
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

  /* ---------- Genetic codes ----------
   * NCBI translation tables, written as their differences from the standard code so the
   * table numbers stay recognisable. Only the codes in real use are listed. */
  const CODE_DIFFS = {
    1: { name: 'Standard', diff: {} },
    2: { name: 'Vertebrate Mitochondrial', diff: { AGA: '*', AGG: '*', ATA: 'M', TGA: 'W' } },
    3: { name: 'Yeast Mitochondrial', diff: { ATA: 'M', CTT: 'T', CTC: 'T', CTA: 'T', CTG: 'T', TGA: 'W' } },
    4: { name: 'Mold, Protozoan & Coelenterate Mitochondrial; Mycoplasma', diff: { TGA: 'W' } },
    5: { name: 'Invertebrate Mitochondrial', diff: { AGA: 'S', AGG: 'S', ATA: 'M', TGA: 'W' } },
    6: { name: 'Ciliate, Dasycladacean & Hexamita Nuclear', diff: { TAA: 'Q', TAG: 'Q' } },
    9: { name: 'Echinoderm & Flatworm Mitochondrial', diff: { AAA: 'N', AGA: 'S', AGG: 'S', TGA: 'W' } },
    10: { name: 'Euplotid Nuclear', diff: { TGA: 'C' } },
    11: { name: 'Bacterial, Archaeal & Plant Plastid', diff: {} },
    12: { name: 'Alternative Yeast Nuclear', diff: { CTG: 'S' } },
    13: { name: 'Ascidian Mitochondrial', diff: { AGA: 'G', AGG: 'G', ATA: 'M', TGA: 'W' } },
    14: { name: 'Alternative Flatworm Mitochondrial', diff: { AAA: 'N', AGA: 'S', AGG: 'S', TAA: 'Y', TGA: 'W' } },
    16: { name: 'Chlorophycean Mitochondrial', diff: { TAG: 'L' } },
    21: { name: 'Trematode Mitochondrial', diff: { AAA: 'N', AGA: 'S', AGG: 'S', ATA: 'M', TGA: 'W' } }
  };
  const GENETIC_CODES = {};
  const CODE_LIST = [];
  for (const [num, spec] of Object.entries(CODE_DIFFS)) {
    const t = Object.assign({}, CODON, spec.diff);
    GENETIC_CODES[num] = t;
    GENETIC_CODES[spec.name] = t;
    CODE_LIST.push({ id: num, name: spec.name, changed: Object.keys(spec.diff) });
  }

  // ---------- Feature segments ----------
  // A feature is one span unless it carries `segs`, in which case those are its exons
  // and the gaps between them are introns.
  function featSegs(f, n) {
    if (f.segs && f.segs.length) return f.segs.slice().sort((a, b) => a.start - b.start);
    return [{ start: f.start, end: f.end }];
  }
  function featLen(f, n) {
    return featSegs(f, n).reduce((t, g) => t + rlen(g, n), 0);
  }
  // The spliced coding sequence, plus the molecule position of every base in it, so a
  // translation can be drawn over exons and numbered like a protein.
  function featCoding(seq, f, n) {
    const segList = featSegs(f, n);
    const order = f.strand === -1 ? segList.slice().reverse() : segList;
    let cds = '', pos = [];
    for (const g of order) {
      const piece = sub(seq, g.start, g.end);
      if (f.strand === -1) {
        cds += rc(piece);
        for (let i = rlen(g, n) - 1; i >= 0; i--) pos.push((g.start + i) % n);
      } else {
        cds += piece;
        for (let i = 0; i < rlen(g, n); i++) pos.push((g.start + i) % n);
      }
    }
    return { cds, pos };
  }
  /* Translate a feature across its exons.
   *   opts.code      a genetic code name from GENETIC_CODES
   *   opts.startMet  translate the first codon as methionine
   *   opts.frame     0, 1 or 2
   *   opts.startAt   the residue number of the first amino acid (default 1)
   * Returns one entry per codon: its letter, residue number and the position of its
   * middle base in the molecule, which is where the letter is drawn. */
  function featTranslate(seq, f, n, opts = {}) {
    const { cds, pos } = featCoding(seq, f, n);
    const fr = opts.frame || 0;
    const table = GENETIC_CODES[opts.code] || null;
    const out = [];
    let res = opts.startAt || 1;
    for (let i = fr; i + 3 <= cds.length; i += 3) {
      const codon = cds.substr(i, 3);
      let aa = (table ? table[codon] : CODON[codon]) || 'X';
      const forced = opts.startMet && i === fr && aa !== 'M';
      if (forced) aa = 'M';
      out.push({ aa, res: res++, mid: pos[i + 1], a: pos[i], c: pos[i + 2], forced, codon });
    }
    return out;
  }

  // ---------- Common feature library (for "Detect Common Features") ----------
  // Common features. `kind` says how a match is found: "dna" compares the sequence
  // directly on both strands, "protein" translates all six frames first, which finds a
  // tag or peptide whatever codons were used for it.
  const COMMON = [
    // --- promoters and transcription signals
    { name: 'T7 promoter', type: 'promoter', kind: 'dna', seq: 'TAATACGACTCACTATAGG' },
    { name: 'SP6 promoter', type: 'promoter', kind: 'dna', seq: 'ATTTAGGTGACACTATAG' },
    { name: 'T3 promoter', type: 'promoter', kind: 'dna', seq: 'AATTAACCCTCACTAAAGG' },
    { name: 'lac promoter', type: 'promoter', kind: 'dna', seq: 'TTTACACTTTATGCTTCCGGCTCGTATGTTG' },
    { name: 'tac promoter', type: 'promoter', kind: 'dna', seq: 'TTGACAATTAATCATCGGCTCGTATAATG' },
    { name: 'trc promoter', type: 'promoter', kind: 'dna', seq: 'TTGACAATTAATCATCCGGCTCGTATAATG' },
    { name: 'T7 terminator', type: 'terminator', kind: 'dna', seq: 'CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG' },
    { name: 'lac operator', type: 'protein_bind', kind: 'dna', seq: 'TTGTGAGCGGATAACAA' },
    { name: 'lac operator (full)', type: 'protein_bind', kind: 'dna', seq: 'GGAATTGTGAGCGGATAACAATT' },
    { name: 'Kozak sequence', type: 'regulatory', kind: 'dna', seq: 'GCCACCATGG' },
    { name: 'poly(A) signal', type: 'polyA_signal', kind: 'dna', seq: 'AATAAA', short: true },
    { name: 'bGH poly(A) signal', type: 'polyA_signal', kind: 'dna', seq: 'CTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCCTGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAGGAAATTGCATCG' },
    { name: 'Shine-Dalgarno', type: 'RBS', kind: 'dna', seq: 'AGGAGG', short: true },

    // --- recombination and integration sites
    { name: 'loxP', type: 'protein_bind', kind: 'dna', seq: 'ATAACTTCGTATAGCATACATTATACGAAGTTAT' },
    { name: 'lox2272', type: 'protein_bind', kind: 'dna', seq: 'ATAACTTCGTATAAAGTATCCTATACGAAGTTAT' },
    { name: 'lox511', type: 'protein_bind', kind: 'dna', seq: 'ATAACTTCGTATAATGTATACTATACGAAGTTAT' },
    { name: 'FRT', type: 'protein_bind', kind: 'dna', seq: 'GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC' },
    { name: 'attB1', type: 'protein_bind', kind: 'dna', seq: 'ACAAGTTTGTACAAAAAAGCAGGCT' },
    { name: 'attB2', type: 'protein_bind', kind: 'dna', seq: 'ACCCAGCTTTCTTGTACAAAGTGGT' },
    { name: 'I-SceI site', type: 'protein_bind', kind: 'dna', seq: 'TAGGGATAACAGGGTAAT' },
    { name: 'Tn5 mosaic end', type: 'misc_feature', kind: 'dna', seq: 'AGATGTGTATAAGAGACAG' },

    // --- sequencing primers
    { name: 'M13 fwd', type: 'primer_bind', kind: 'dna', seq: 'GTAAAACGACGGCCAGT' },
    { name: 'M13 rev', type: 'primer_bind', kind: 'dna', seq: 'CAGGAAACAGCTATGAC' },

    // --- affinity and epitope tags, matched on the protein
    { name: '6xHis', type: 'CDS', kind: 'protein', seq: 'HHHHHH' },
    { name: '8xHis', type: 'CDS', kind: 'protein', seq: 'HHHHHHHH' },
    { name: '10xHis', type: 'CDS', kind: 'protein', seq: 'HHHHHHHHHH' },
    { name: 'FLAG', type: 'CDS', kind: 'protein', seq: 'DYKDDDDK' },
    { name: '3xFLAG', type: 'CDS', kind: 'protein', seq: 'DYKDHDGDYKDHDIDYKDDDDK' },
    { name: 'HA', type: 'CDS', kind: 'protein', seq: 'YPYDVPDYA' },
    { name: 'Myc', type: 'CDS', kind: 'protein', seq: 'EQKLISEEDL' },
    { name: 'V5 tag', type: 'CDS', kind: 'protein', seq: 'GKPIPNPLLGLDST' },
    { name: 'Strep-tag II', type: 'CDS', kind: 'protein', seq: 'WSHPQFEK' },
    { name: 'S-tag', type: 'CDS', kind: 'protein', seq: 'KETAAAKFERQHMDS' },
    { name: 'T7 tag', type: 'CDS', kind: 'protein', seq: 'MASMTGGQQMG' },
    { name: 'Avi-tag', type: 'CDS', kind: 'protein', seq: 'GLNDIFEAQKIEWHE' },
    { name: 'Sortase A site', type: 'CDS', kind: 'protein', seq: 'LPETG', short: true },

    // --- protease cleavage sites
    { name: 'TEV site', type: 'CDS', kind: 'protein', seq: 'ENLYFQ', short: true },
    { name: 'HRV 3C / PreScission site', type: 'CDS', kind: 'protein', seq: 'LEVLFQGP' },
    { name: 'Thrombin site', type: 'CDS', kind: 'protein', seq: 'LVPRGS', short: true },
    { name: 'Factor Xa site', type: 'CDS', kind: 'protein', seq: 'IEGR', short: true },
    { name: 'Enterokinase site', type: 'CDS', kind: 'protein', seq: 'DDDDK', short: true },

    // --- localisation signals
    { name: 'SV40 NLS', type: 'CDS', kind: 'protein', seq: 'PKKKRKV' },
    { name: 'Nucleoplasmin NLS', type: 'CDS', kind: 'protein', seq: 'KRPAATKKAGQAKKKK' },

    // --- self-cleaving 2A peptides
    { name: 'P2A', type: 'CDS', kind: 'protein', seq: 'ATNFSLLKQAGDVEENPGP' },
    { name: 'T2A', type: 'CDS', kind: 'protein', seq: 'EGRGSLLTCGDVEENPGP' },
    { name: 'E2A', type: 'CDS', kind: 'protein', seq: 'QCTNYALLKLAGDVESNPGP' },
    { name: 'F2A', type: 'CDS', kind: 'protein', seq: 'VKQTLNFDLLKLAGDVESNPGP' },

    // --- linkers
    { name: '(G4S)3 linker', type: 'CDS', kind: 'protein', seq: 'GGGGSGGGGSGGGGS' },
    { name: 'G4S linker', type: 'CDS', kind: 'protein', seq: 'GGGGS', short: true }
  ];

  // `short` entries are left out by default: a motif of a few residues turns up by chance.
  /* Detection walks a long list — well over a thousand standard features — so each one
   * gets a cheap look before an expensive one. The sequence is laid out once, forwards,
   * reverse complemented and translated in all six frames, and a feature only reaches
   * the real search if a plain substring of it turns up somewhere in that. Nothing is
   * decided by the quick look: it only rules candidates out, so what is found is exactly
   * what the full search would have found on its own. */
  const PROBE = 14;
  function prepScan(seq, circular, list) {
    let longest = 0;
    for (const f of list) if (f && f.seq) longest = Math.max(longest, f.kind === 'protein' ? f.seq.length * 3 : f.seq.length);
    const ext = circular && longest > 1 ? seq + seq.slice(0, Math.min(seq.length, longest - 1)) : seq;
    const p = { dna: ext, rna: null, rc: rc(ext), aa: [] };
    for (let f = 0; f < 3; f++) p.aa.push(translate(ext.slice(f)), translate(p.rc.slice(f)));
    return p;
  }
  // the first stretch of a feature with no ambiguity code in it, or '' when it has none
  function probeOf(f) {
    if (f.kind === 'protein') return f.seq.slice(0, Math.min(10, f.seq.length)).replace(/\*/g, '');
    const s = String(f.seq).toUpperCase();
    if (s.length < PROBE) return /^[ACGT]+$/.test(s) ? s : '';
    for (let i = 0; i + PROBE <= s.length; i++) {
      const w = s.slice(i, i + PROBE);
      if (/^[ACGT]+$/.test(w)) return w;
    }
    return '';
  }
  function mayBeThere(f, p, probe) {
    if (!probe) return true;                       // nothing cheap to test: let it through
    if (f.kind === 'protein') return p.aa.some(fr => fr.indexOf(probe) >= 0);
    return p.dna.indexOf(probe) >= 0 || p.rc.indexOf(probe) >= 0;
  }

  /* The core list above is what the program knows on its own. features-db.js, when it
   * ships, adds the rest of the standard library beside it; anything already named here
   * with the same sequence is left alone, so the core stays the one definition. */
  (function mergeLibrary() {
    const lib = root.FEATURE_DB;
    if (!Array.isArray(lib)) return;
    const have = new Set(COMMON.map(f => (f.kind || 'dna') + ':' + f.seq));
    for (const f of lib) {
      if (!f || !f.seq) continue;
      const k = (f.kind || 'dna') + ':' + f.seq;
      if (have.has(k)) continue;
      have.add(k);
      COMMON.push(f);
    }
  })();

  function detectCommon(seq, circular, existing, opts) {
    opts = opts || {};
    const found = [];
    const list = COMMON.concat(opts.extra || []);
    const prep = prepScan(seq, circular, list);
    for (const f of list) {
      if (!f || !f.seq) continue;
      if (f.short && !opts.includeShort) continue;
      if (!mayBeThere(f, prep, probeOf(f))) continue;
      const hits = f.kind === 'protein' ? searchProtein(seq, circular, f.seq) : searchDNA(seq, circular, f.seq);
      for (const h of hits) {
        if (existing.some(x => x.start === h.start && x.end === h.end)) continue;
        if (found.some(x => x.start === h.start && x.end === h.end)) continue;
        found.push({
          name: f.name, type: f.type || 'misc_feature', start: h.start, end: h.end, strand: h.strand,
          via: f.kind, molecule: f.molecule || (f.kind === 'protein' ? 'protein' : 'DNA'),
          matchProtein: f.kind === 'protein', color: f.color || null, mine: !!f.mine
        });
      }
    }
    return found;
  }

  // ---------- File formats ----------
  const HIST_TAG = 'DogmaStudio-History-v1';
  // files written before the program was renamed carry the old tag; they still open
  const HIST_TAGS = [HIST_TAG, 'GeneCanvas-History-v1'];
  function b64encode(str) {
    const bytes = new TextEncoder().encode(str);
    let bin = '';
    for (const b of bytes) bin += String.fromCharCode(b);
    return typeof btoa === 'function' ? btoa(bin) : Buffer.from(bytes).toString('base64');
  }
  function b64decode(b64) {
    if (typeof atob === 'function') {
      const bin = atob(b64);
      const bytes = new Uint8Array(bin.length);
      for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
      return new TextDecoder().decode(bytes);
    }
    return Buffer.from(b64, 'base64').toString('utf8');
  }

  function parseLocation(loc, n, circular) {
    const strand = /complement/.test(loc) ? -1 : 1;
    const parts = [];
    { const re2 = /<?(\d+)\.\.>?(\d+)/g; let mm; while ((mm = re2.exec(loc))) parts.push([+mm[1] - 1, +mm[2]]); }
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
        return { start, end: end % n === 0 ? n : end % n, strand, parts };
      }
    }
    start = Math.min(...ranges.map(r => r[0])) - 1;
    end = Math.max(...ranges.map(r => r[1]));
    if (start > end - 1 && circular) { start = first[0] - 1; end = last[1]; }
    return { start, end, strand, parts };
  }

  function parseGenBank(txt) {
    const lines = txt.split(/\r?\n/);
    let name = 'Untitled', circular = false, mode = '', seq = '', cur = null, lastQ = null, definition = '';
    let histLines = null;
    // The header fields the Info panel shows, kept as they are read.
    const header = { accession: '', organism: '', source: '', division: '', comments: [], references: [] };
    let refKey = null, comment = null;
    const raw = [];
    for (const ln of lines) {
      if (ln.startsWith('LOCUS')) {
        const t = ln.split(/\s+/); name = t[1] || name; circular = /\bcircular\b/i.test(ln);
        const dv = ln.slice(64, 67).trim(); if (/^[A-Z]{3}$/.test(dv)) header.division = dv;
        continue;
      }
      if (ln.startsWith('DEFINITION')) { definition = ln.slice(12).trim(); mode = 'D'; continue; }
      if (mode === 'D' && /^\s{5,}\S/.test(ln)) { definition += ' ' + ln.trim(); continue; }
      if (ln.startsWith('ACCESSION')) { header.accession = ln.slice(12).trim().split(/\s+/)[0] || ''; mode = ''; continue; }
      if (ln.startsWith('SOURCE')) { header.source = ln.slice(12).trim(); mode = 'S'; continue; }
      if (/^\s{2}ORGANISM/.test(ln)) { header.organism = ln.slice(12).trim(); mode = 'S'; continue; }
      if (ln.startsWith('REFERENCE')) { header.references.push({ authors: '', title: '', journal: '', pmid: '' }); refKey = null; mode = 'R'; continue; }
      if (mode === 'R') {
        const cur2 = header.references[header.references.length - 1];
        const m2 = ln.match(/^\s{2}(AUTHORS|TITLE|JOURNAL)\s+(.*)$/);
        if (m2) { refKey = m2[1].toLowerCase(); cur2[refKey] = m2[2].trim(); continue; }
        const mp = ln.match(/^\s{3}PUBMED\s+(\d+)/);
        if (mp) { cur2.pmid = mp[1]; refKey = null; continue; }
        if (/^\s{5,}\S/.test(ln) && refKey) { cur2[refKey] += ' ' + ln.trim(); continue; }
        if (/^\S/.test(ln)) { mode = ''; refKey = null; }
      }
      if (ln.startsWith('COMMENT')) {
        mode = 'C';
        if (HIST_TAGS.includes(ln.slice(12).trim())) { histLines = []; comment = null; }
        else { comment = [ln.slice(12).trim()]; header.comments.push(comment); }
        continue;
      }
      if (mode === 'C') {
        if (/^\s/.test(ln)) {
          if (histLines) histLines.push(ln.trim());
          else if (comment) comment.push(ln.trim());
          continue;
        }
        mode = ''; comment = null; // a new top-level key ends the comment
      }
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
      const nf = { name: fname.slice(0, 80), type: f.type, start: loc.start, end: loc.end, strand: loc.strand, color, notes: f.q.note || '' };
      // A join() over several ranges is a segmented feature: exons with introns between them.
      if (loc.parts && loc.parts.length > 1 && !(loc.parts.length === 2 && loc.parts[0][1] === n && loc.parts[1][0] === 0))
        nf.segs = loc.parts.map(([a, b]) => ({ start: a, end: b }));
      features.push(nf);
    }
    let provenance = null, info = null;
    if (histLines && histLines.length) {
      try {
        const obj = JSON.parse(b64decode(histLines.join('')));
        // the block used to hold the provenance alone; now it is a wrapper
        if (obj && obj.__gc) { provenance = obj.provenance || null; info = obj.info || null; }
        else provenance = obj;
      } catch (e) { provenance = null; }
    }
    if (!info) {
      info = {
        molType: /synthetic/i.test(header.source || '') ? 'synthetic' : (header.source ? 'natural' : ''),
        organism: header.organism || '', accession: header.accession || '',
        description: definition.replace(/\.$/, ''), comments: (header.comments || []).map(c => c.join(' ')).join('\n\n'),
        references: header.references || [], author: '', code: '', host: '', strain: '', methylation: '',
        seqClass: header.division || '', confirmed: false, files: []
      };
    }
    return { name, seq, circular, features, primers, definition, provenance, info };
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

  /* The .dna binary layout: a stream of blocks, each a type byte, a big-endian
   * uint32 length and that many bytes of payload. Read only; nothing is written in
   * this format. */
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
    if (!sawCookie && !seq) throw new Error('This file is not in the .dna format.');
    seq = clean(seq);
    if (!seq) throw new Error('The .dna file has no DNA sequence block.');
    const n = seq.length;
    for (const f of features) { if (f.end > n) f.end = n; if (f.start >= n) f.start = 0; if (f.color && !/^#[0-9a-fA-F]{6}$/.test(f.color)) f.color = null; }
    return { name: (filename || 'Untitled').replace(/\.[^.]+$/, ''), seq, circular, features, primers };
  }

  function gbLoc(f, n) {
    let s;
    if (f.segs && f.segs.length > 1) {
      // exons in ascending order, which is how GenBank writes a join
      const parts = f.segs.slice().sort((a, b) => a.start - b.start)
        .map(g => g.end > g.start ? `${g.start + 1}..${g.end}` : `${g.start + 1}..${n},1..${g.end}`);
      s = `join(${parts.join(',')})`;
    }
    else if (f.end > f.start) s = (f.end - f.start === 1) ? `${f.start + 1}` : `${f.start + 1}..${f.end}`;
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
    const info = d.info || {};
    const natural = info.molType === 'natural';
    const org = (natural && info.organism) ? info.organism : 'synthetic DNA construct';
    const division = natural ? (info.seqClass || 'UNA').slice(0, 3).toUpperCase() : 'SYN';
    const wrapField = (key, text, width = 67) => {
      const words = String(text).replace(/\s+/g, ' ').trim().split(' ');
      const lines = []; let line = '';
      for (const w of words) { if ((line + ' ' + w).trim().length > width) { lines.push(line.trim()); line = w; } else line += ' ' + w; }
      if (line.trim()) lines.push(line.trim());
      if (!lines.length) lines.push('');
      return lines.map((l, i) => (i === 0 ? key.padEnd(12) : ' '.repeat(12)) + l).join('\n') + '\n';
    };
    let s = `LOCUS       ${nm.padEnd(16)} ${String(n).padStart(11)} bp    DNA     ${d.circular ? 'circular' : 'linear  '} ${division} ${date}\n`;
    s += wrapField('DEFINITION', (info.description || d.name).replace(/\.$/, '') + '.');
    s += `ACCESSION   ${info.accession || '.'}\nVERSION     ${info.accession || '.'}\n`;
    s += `KEYWORDS    ${[info.code, 'Dogma Studio'].filter(Boolean).join('; ')}\n`;
    s += wrapField('SOURCE', org);
    s += wrapField('  ORGANISM', org);
    (info.references || []).forEach((r, i) => {
      s += `REFERENCE   ${i + 1}${r.from ? `  (bases ${r.from} to ${r.to || n})` : ''}\n`;
      if (r.authors) s += wrapField('  AUTHORS', r.authors);
      if (r.title) s += wrapField('  TITLE', r.title);
      if (r.journal) s += wrapField('  JOURNAL', r.journal);
      if (r.pmid) s += `   PUBMED   ${r.pmid}\n`;
    });
    if (info.author) { s += `REFERENCE   ${(info.references || []).length + 1}  (bases 1 to ${n})\n`; s += wrapField('  AUTHORS', info.author); s += wrapField('  TITLE', 'Direct Submission'); }
    const notes = [];
    if (info.comments) notes.push(info.comments);
    if (info.molType === 'synthetic') notes.push(`Synthetic DNA${info.confirmed ? ', confirmed experimentally' : ''}${info.host ? `; laboratory host ${info.host}` : ''}${info.strain ? `; transformation strain ${info.strain}` : ''}${info.methylation ? `; methylation ${info.methylation}` : ''}.`);
    for (const t of notes) s += wrapField('COMMENT', t);
    // Anything GenBank has no field for rides in a base64 block any reader will skip.
    const extra = {};
    if (d.provenance) extra.provenance = d.provenance;
    if (d.info) extra.info = d.info;
    if (Object.keys(extra).length) {
      const enc = b64encode(JSON.stringify({ __gc: 1, ...extra }));
      s += 'COMMENT     ' + HIST_TAG + '\n';
      for (let i = 0; i < enc.length; i += 68) s += '            ' + enc.slice(i, i + 68) + '\n';
    }
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
    return { name: 'pDogma-Demo', seq, circular: true, features, primers, demo: true };
  }

  root.DNA = {
    COMP, clean, rc, comp, CODON, translate, AA3, proteinMW, gcPct, dsDnaMW, tm,
    GENETIC_CODES, CODE_LIST, featSegs, featLen, featCoding, featTranslate,
    NT_MW, ntMW, strandMW, composition, gcCurve,
    ENZYMES, ENZ_BY_NAME, baseEnzymeName, findSites, rlen, segs, sub, overlaps, findORFs, searchDNA, searchProtein,
    digest, detectCommon, COMMON, parseGenBank, parseFASTA, parseAny, parseDnaBinary, toGenBank, toFASTA, buildDemo
  };
})(typeof window !== 'undefined' ? window : globalThis);
