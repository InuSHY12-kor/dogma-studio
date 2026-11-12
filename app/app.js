/* GeneCanvas — application UI */
(function () {
  'use strict';
  const DNA = window.DNA, CLONE = window.CLONE, TRACE = window.TRACE, PRIMER = window.PRIMER, MSA = window.MSA, CODON = window.CODON, STRUCT = window.STRUCT, GEL = window.GEL;
  const host = window.gcHost || null; // provided by the Electron preload; null in a plain browser
  const $ = s => document.querySelector(s);
  const esc = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const fmt = n => Number(n).toLocaleString('en-US');
  const TAU = Math.PI * 2;

  const TYPE_COL = {
    CDS: '#f4c95d', gene: '#f4c95d', promoter: '#7cc37a', terminator: '#e06b63', rep_origin: '#f0e27a',
    misc_feature: '#c9ced6', primer_bind: '#8fb3e8', protein_bind: '#69c4c6', RBS: '#c69be0', polyA_signal: '#e8a06a',
    enhancer: '#a3d56f', sig_peptide: '#f2a3c9', LTR: '#a8a8de', regulatory: '#9fd6b4', intron: '#d8d8d8', exon: '#f6d98a',
    "5'UTR": '#cfe3a5', "3'UTR": '#cfe3a5', mobile_element: '#d9b38c', misc_RNA: '#b4d4f0', ncRNA: '#b4d4f0', ORF: '#b9bfc8'
  };
  // Rules for the map SVG. Kept here so exported and printed copies can carry them;
  // a stylesheet loaded from a file cannot be read back at runtime.
  const SVG_CSS = `
svg .bb { fill: none; stroke: var(--bb); stroke-width: 1.2; }
svg .tick { stroke: var(--faint); stroke-width: 1; }
svg .tickl { fill: var(--faint); font: 10px var(--f-ui); }
svg .feat { stroke: rgba(0, 0, 0, .38); stroke-width: .8; cursor: pointer; }
svg .feat:hover { stroke: #000; stroke-width: 1.6; }
svg .ftxt { fill: #111; font: 600 11px var(--f-ui); pointer-events: none; }
svg .lbl { font: 11px var(--f-ui); fill: var(--ink); cursor: pointer; }
svg .lbl.e { fill: #3a3f4b; }
svg .lbl.e.u { font-weight: 600; }
svg .lbl.f { font-weight: 600; }
svg .lbl:hover { fill: var(--accent); text-decoration: underline; }
svg .lead { fill: none; stroke: #9aa2af; stroke-width: .8; }
svg .ecut { stroke: #48505e; stroke-width: 1; }
svg .ctr1 { font: 700 17px var(--f-ui); fill: var(--ink); }
svg .ctr2 { font: 13px var(--f-ui); fill: var(--muted); }
svg .selarea { fill: var(--sel); stroke: var(--selline); stroke-width: .8; pointer-events: none; }
svg .caret { stroke: var(--selline); stroke-width: 1.5; pointer-events: none; }
svg .orf { stroke: none; cursor: pointer; }
svg .prim { cursor: pointer; }
`;
  const FEATURE_TYPES = Object.keys(TYPE_COL).filter(t => t !== 'ORF');
  const CDS_PALETTE = ['#f4c95d', '#8ec3f2', '#f29e8e', '#a6dd8b', '#d6a8f0', '#f6b26b', '#7fd1c7'];

  // ---------------------------------------------------------------- state
  const S = { docs: [], active: -1, view: 'map', zoom: 1, sideHidden: false, seq: null, mapGeom: null, find: { q: '', mode: 'dna', res: [], i: -1 } };
  let UID = 1; const uid = () => UID++;
  const D = () => S.docs[S.active] || null;

  function defaultOpts(circular) {
    return {
      mapMode: circular ? 'circular' : 'linear', showFeat: true, showPrimers: true, showORFs: false,
      orfMin: 100, showTrans: true, enzMode: 'unique6', showComp: true,
      labelMode: 'auto',      // auto | inside | outside
      featOnLine: false,      // draw every feature on the backbone instead of in its own lane
      showRuler: true,
      baseColors: { enabled: false, A: '#2e9e4f', C: '#2563c9', G: '#1c2230', T: '#c43c3c' }
    };
  }

  function makeDoc(src) {
    const n = src.seq.length;
    let cdsI = 0;
    const d = {
      id: uid(), name: src.name || 'Untitled', seq: DNA.clean(src.seq), circular: !!src.circular,
      features: (src.features || []).map(f => {
        let color = f.color;
        if (!color) color = f.type === 'CDS' ? CDS_PALETTE[(cdsI++) % CDS_PALETTE.length] : (TYPE_COL[f.type] || '#c9ced6');
        return { id: uid(), name: f.name || f.type, type: f.type || 'misc_feature', start: Math.max(0, Math.min(n, f.start)), end: Math.max(0, Math.min(n, f.end)), strand: f.strand === undefined ? 1 : f.strand, color, notes: f.notes || '' };
      }),
      primers: (src.primers || []).map(p => ({ id: uid(), name: p.name, seq: DNA.clean(p.seq), notes: p.notes || '' })),
      provenance: src.provenance || null,
      opts: Object.assign(defaultOpts(!!src.circular), SET.mapView || {}, SET.baseColors ? { baseColors: { ...SET.baseColors } } : {}), chosen: [], sel: { start: 0, end: 0 }, anchor: 0, alignments: [], designed: [],
      hist: [], hi: -1, ver: 0, cache: null, dirty: false, path: src.path || null
    };
    commit(d, src.histLabel || 'Opened ' + d.name, { initial: true });
    return d;
  }

  function snapshot(d) { return JSON.stringify({ name: d.name, seq: d.seq, circular: d.circular, features: d.features, primers: d.primers }); }
  function commit(d, label, o = {}) {
    d.hist = d.hist.slice(0, d.hi + 1);
    const last = d.hist[d.hist.length - 1];
    const entry = { label, t: new Date(), snap: snapshot(d), len: d.seq.length };
    if (o.merge && last && last.mergeKey === o.merge && (Date.now() - last.t.getTime()) < 2500) { entry.mergeKey = o.merge; d.hist[d.hist.length - 1] = entry; }
    else { if (o.merge) entry.mergeKey = o.merge; d.hist.push(entry); }
    if (d.hist.length > 300) d.hist.shift();
    d.hi = d.hist.length - 1;
    d.ver++; d.cache = null;
    if (!o.initial) d.dirty = true;
  }
  function restore(d, i) {
    const s = JSON.parse(d.hist[i].snap);
    d.name = s.name; d.seq = s.seq; d.circular = s.circular; d.features = s.features; d.primers = s.primers;
    d.hi = i; d.ver++; d.cache = null; d.dirty = true;
    const n = d.seq.length;
    d.sel = { start: Math.min(d.sel.start, n), end: Math.min(d.sel.end, n) };
    if (!d.circular && d.sel.end < d.sel.start) d.sel = { start: d.sel.end, end: d.sel.start };
  }
  function undo() { const d = D(); if (d && d.hi > 0) { restore(d, d.hi - 1); toast('Undo: ' + d.hist[d.hi + 1].label); renderAll(); } }
  function redo() { const d = D(); if (d && d.hi < d.hist.length - 1) { restore(d, d.hi + 1); toast('Redo: ' + d.hist[d.hi].label); renderAll(); } }

  // ---------------------------------------------------------------- analysis
  function primerBind(d, p) {
    const n = d.seq.length, s = p.seq; if (!s || s.length > n || s.length < 6) return [];
    const out = [];
    const push = (h, full) => out.push({ start: h.start, end: h.end, strand: h.strand, full });
    const hits = DNA.searchDNA(d.seq, d.circular, s);
    for (const h of hits) push(h, true);
    if (!hits.length && s.length > 15) {
      const tail = s.slice(-15);
      for (const h of DNA.searchDNA(d.seq, d.circular, tail)) push(h, false);
    }
    return out;
  }
  function analyze(d) {
    if (d.cache && d.cache.ver === d.ver && d.cache.orfMin === d.opts.orfMin) return d.cache;
    const t = { ver: d.ver, orfMin: d.opts.orfMin };
    t.sites = DNA.findSites(d.seq, d.circular);
    t.orfs = DNA.findORFs(d.seq, d.circular, d.opts.orfMin);
    t.psites = {};
    for (const p of d.primers) t.psites[p.id] = primerBind(d, p);
    d.cache = t; return t;
  }
  function shownEnz(d, A) {
    const m = d.opts.enzMode, out = [];
    if (m === 'none') return out;
    for (const e of DNA.ENZYMES) {
      const s = A.sites[e.name]; if (!s) continue;
      let ok = false;
      if (m === 'unique6') ok = s.length === 1 && e.core >= 6;
      else if (m === 'unique') ok = s.length === 1;
      else if (m === 'dual6') ok = s.length <= 2 && e.core >= 6;
      else if (m === 'chosen') ok = d.chosen.includes(e.name);
      if (ok) out.push({ name: e.name, sites: s, unique: s.length === 1 });
    }
    return out;
  }
  function collectItems(d, A) {
    const items = [], o = d.opts;
    if (o.showFeat) for (const f of d.features) items.push({ kind: 'f', key: 'f:' + f.id, ref: f, start: f.start, end: f.end, strand: f.strand, name: f.name, color: f.color, type: f.type });
    if (o.showPrimers) for (const p of d.primers) (A.psites[p.id] || []).forEach((b, i) => items.push({ kind: 'p', key: 'p:' + p.id + ':' + i, ref: p, start: b.start, end: b.end, strand: b.strand, name: p.name, color: b.full ? '#2f6fd0' : '#7d9fd6', type: 'primer' }));
    if (o.showORFs) A.orfs.forEach((r, i) => items.push({ kind: 'o', key: 'o:' + i, ref: r, start: r.start, end: r.end, strand: r.strand, name: 'ORF: ' + r.aa + ' aa', color: '#9aa2ae', type: 'ORF' }));
    return items;
  }
  function assignLanes(list, n) {
    const lanes = [];
    const sorted = list.slice().sort((a, b) => DNA.rlen(b, n) - DNA.rlen(a, n));
    for (const it of sorted) {
      let k = 0;
      while (lanes[k] && lanes[k].some(o => DNA.overlaps(o, it, n))) k++;
      (lanes[k] = lanes[k] || []).push(it); it.lane = k;
    }
    return lanes.length;
  }
  function findItem(d, key) {
    const A = analyze(d);
    return collectItems(d, A).find(i => i.key === key) || null;
  }

  // ---------------------------------------------------------------- selection helpers
  function selLen(d) { return DNA.rlen(d.sel, d.seq.length); }
  function selSeq(d) { return DNA.sub(d.seq, d.sel.start, d.sel.end); }
  function setSel(start, end, o = {}) {
    const d = D(); if (!d) return;
    const n = d.seq.length;
    start = Math.max(0, Math.min(n, start)); end = Math.max(0, Math.min(n, end));
    if (!d.circular && end < start) [start, end] = [end, start];
    d.sel = { start, end };
    if (o.anchor !== undefined) d.anchor = o.anchor;
    refreshSel(o.scroll);
  }
  function selectRange(r, o = {}) { setSel(r.start, r.end, Object.assign({ anchor: r.start, scroll: true }, o)); }

  // ---------------------------------------------------------------- editing
  function editSeq(d, pos, del, ins) {
    d.seq = d.seq.slice(0, pos) + ins + d.seq.slice(pos + del);
    const L = ins.length, nn = d.seq.length;
    const mS = x => x < pos ? x : (x >= pos + del ? x - del + L : pos + L);
    const mE = x => x <= pos ? x : (x >= pos + del ? x - del + L : pos);
    d.features = d.features.filter(f => {
      const wrap = f.end <= f.start;
      let s = mS(f.start), e = mE(f.end);
      if (!wrap && del > 0 && f.start >= pos && f.end <= pos + del) return false;
      if (!wrap && e <= s) return false;
      if (wrap && s === e) return false;
      if (s >= nn) s = nn - 1; if (e > nn) e = nn;
      f.start = s; f.end = e; return true;
    });
  }
  function deleteRange(d, r) {
    const n = d.seq.length;
    if (r.end > r.start) editSeq(d, r.start, r.end - r.start, '');
    else if (r.end < r.start) { editSeq(d, r.start, n - r.start, ''); editSeq(d, 0, r.end, ''); }
  }
  function replaceSel(text, label, merge) {
    const d = D(); if (!d) return;
    text = DNA.clean(text);
    const L = selLen(d);
    let pos = d.sel.start;
    if (L) {
      if (d.sel.end < d.sel.start) { deleteRange(d, d.sel); pos = d.seq.length; }
      else deleteRange(d, d.sel);
    }
    if (!text && !L) return;
    if (!d.seq.length && !text) { toast('A sequence cannot be empty.'); restore(d, d.hi); renderAll(); return; }
    editSeq(d, pos, 0, text);
    d.sel = { start: pos + text.length, end: pos + text.length }; d.anchor = d.sel.start;
    commit(d, label, merge ? { merge } : {});
    renderAll(true);
  }
  function deleteSel() {
    const d = D(); if (!d || !selLen(d)) return;
    if (selLen(d) >= d.seq.length) { toast('Cannot delete the whole sequence.'); return; }
    const L = selLen(d);
    const pos = d.sel.end < d.sel.start ? 0 : d.sel.start;
    deleteRange(d, d.sel);
    const p = Math.min(pos, d.seq.length);
    d.sel = { start: p, end: p }; d.anchor = p;
    commit(d, `Deleted ${fmt(L)} bp`);
    renderAll(true);
  }

  function setOrigin() {
    const d = D(); if (!d) return;
    if (!d.circular) { toast('Set Origin needs a circular sequence.'); return; }
    const p = d.sel.start, n = d.seq.length; if (!p) { toast('The origin is already at the cursor.'); return; }
    d.seq = d.seq.slice(p) + d.seq.slice(0, p);
    for (const f of d.features) { f.start = (f.start - p + n) % n; let e = (f.end - p + n) % n; if (e === 0) e = n; f.end = e; }
    d.sel = { start: 0, end: 0 };
    commit(d, `Set origin at ${fmt(p + 1)}`); renderAll(); toast('New origin set at the former position ' + fmt(p + 1));
  }
  function revCompAll() {
    const d = D(); if (!d) return; const n = d.seq.length;
    d.seq = DNA.rc(d.seq);
    for (const f of d.features) { const s = n - f.end, e = n - f.start; f.start = s; f.end = e; f.strand = -f.strand || 0; }
    d.sel = { start: 0, end: 0 };
    commit(d, 'Reverse-complemented entire sequence'); renderAll();
  }
  function revCompSel() {
    const d = D(); if (!d) return; const L = selLen(d);
    if (!L) { revCompAll(); return; }
    if (d.sel.end < d.sel.start) { toast('Select a region that does not span the origin.'); return; }
    const s = d.sel.start, e = d.sel.end;
    const r = DNA.rc(d.seq.slice(s, e));
    d.seq = d.seq.slice(0, s) + r + d.seq.slice(e);
    for (const f of d.features) if (f.start >= s && f.end <= e && f.end > f.start) { const a = s + (e - f.end), b = s + (e - f.start); f.start = a; f.end = b; f.strand = -f.strand || 0; }
    commit(d, `Reverse-complemented ${fmt(L)} bp`); renderAll();
  }
  function setTopology(circ) {
    const d = D(); if (!d || d.circular === circ) return;
    d.circular = circ; const n = d.seq.length;
    if (!circ) for (const f of d.features) if (f.end <= f.start) f.end = n;
    if (!circ && d.sel.end < d.sel.start) d.sel = { start: d.sel.start, end: n };
    d.opts.mapMode = circ ? 'circular' : 'linear';
    commit(d, 'Changed topology to ' + (circ ? 'circular' : 'linear')); renderAll();
  }

  // ---------------------------------------------------------------- render orchestration
  function renderAll(keepScroll) {
    renderMenuState(); renderToolbar(); renderDocTabs(); renderViewTabs(); renderSide(); renderView(keepScroll); renderStatus();
    document.title = D() ? `${D().name}${D().dirty ? ' •' : ''} — GeneCanvas` : 'GeneCanvas';
  }
  function refreshSel(scroll) {
    if (S.view === 'map') drawSelLayer();
    else if (S.view === 'seq') { updateSeqSel(); if (scroll) scrollSeqTo(D().sel.start); }
    else if (S.view === 'enz' || S.view === 'feat' || S.view === 'primer') { /* tables don't show selection */ }
    renderSide(); renderStatus();
  }

  // ---------------------------------------------------------------- icons
  const IC = {
    new: '<path d="M6 3h8l4 4v14H6z"/><path d="M14 3v4h4"/><path d="M12 11v6M9 14h6"/>',
    open: '<path d="M3 7h6l2 2h10v10H3z"/><path d="M3 7V5h6l2 2"/>',
    save: '<path d="M5 3h12l3 3v15H5z"/><path d="M8 3v6h8V3"/><rect x="8" y="13" width="8" height="6"/>',
    undo: '<path d="M9 7 4 12l5 5"/><path d="M4 12h11a5 5 0 0 1 0 10h-3"/>',
    redo: '<path d="m15 7 5 5-5 5"/><path d="M20 12H9a5 5 0 0 0 0 10h3"/>',
    find: '<circle cx="10.5" cy="10.5" r="6"/><path d="m15 15 6 6"/>',
    circ: '<circle cx="12" cy="12" r="8"/><path d="M12 4a8 8 0 0 1 8 8" stroke-width="3.2"/>',
    lin: '<path d="M3 12h18"/><path d="M6 9v6M18 9v6"/><rect x="8" y="10" width="7" height="4" fill="currentColor" stroke="none"/>',
    feat: '<path d="M3 9h13l4 3-4 3H3z"/>',
    primer: '<path d="M3 12h14"/><path d="m13 8 4 4-4 4"/><path d="M3 18h7"/>',
    enz: '<path d="M6 4l12 16M18 4 6 20"/><circle cx="6" cy="4" r="1.6"/><circle cx="18" cy="4" r="1.6"/>',
    gel: '<rect x="4" y="3" width="16" height="18" rx="1"/><path d="M7 7h3M14 7h3M7 11h3M14 13h3M7 16h3M14 17h3"/>',
    origin: '<circle cx="12" cy="12" r="8"/><path d="M12 2v6"/><circle cx="12" cy="12" r="1.5" fill="currentColor"/>',
    rc: '<path d="M4 8h14l-3-3M20 16H6l3 3"/>',
    tr: '<path d="M4 6h16M4 12h10M4 18h6"/><path d="m16 15 3 3-3 3"/>',
    side: '<rect x="3" y="4" width="18" height="16" rx="1"/><path d="M9 4v16"/>',
    map: '<circle cx="12" cy="12" r="7"/><path d="M12 5a7 7 0 0 1 7 7" stroke-width="3"/>',
    seq: '<path d="M4 7h16M4 12h16M4 17h10"/>',
    tbl: '<rect x="3" y="5" width="18" height="14"/><path d="M3 10h18M9 5v14"/>',
    hist: '<circle cx="12" cy="12" r="8"/><path d="M12 7v5l3 2"/>',
    align: '<path d="M4 7h16M4 12h11M4 17h14"/><path d="M19 10v4"/>',
    compare: '<path d="M4 6h7M4 12h7M4 18h7"/><path d="M15 6h5M15 12h5M15 18h5"/><path d="M12 4v16" stroke-dasharray="2 2"/>',
    codon: '<path d="M5 8h14M5 12h14M5 16h14"/><circle cx="9" cy="8" r="1.6" fill="currentColor" stroke="none"/><circle cx="14" cy="16" r="1.6" fill="currentColor" stroke="none"/>',
    ncbi: '<path d="M12 3v12"/><path d="m7 11 5 5 5-5"/><path d="M4 20h16"/>',
    clone: '<circle cx="7" cy="8" r="3.4"/><circle cx="17" cy="16" r="3.4"/><path d="M9.6 10.4 14.4 14"/>',
    design: '<path d="M4 16h8"/><path d="m9 13 3 3-3 3"/><path d="M4 7h16"/><circle cx="17" cy="7" r="2"/>',
    print: '<path d="M7 9V4h10v5"/><rect x="4" y="9" width="16" height="7" rx="1"/><path d="M7 14h10v6H7z"/>',
    zoomin: '<circle cx="10.5" cy="10.5" r="6"/><path d="m15 15 6 6M8 10.5h5M10.5 8v5"/>',
    zoomout: '<circle cx="10.5" cy="10.5" r="6"/><path d="m15 15 6 6M8 10.5h5"/>'
  };
  const icon = k => `<svg viewBox="0 0 24 24">${IC[k]}</svg>`;

  // ---------------------------------------------------------------- menus
  function menus() {
    const d = D(), has = !!d && !isGel(d), sl = has && selLen(d) > 0;
    return [
      { label: 'File', items: [
        { label: 'New DNA File…', key: 'Ctrl+N', act: dlgNew },
        { label: 'Open…', key: 'Ctrl+O', act: openFiles },
        { label: 'Open from NCBI…', key: 'Ctrl+D', act: dlgNCBI },
        { label: 'Open Sequencing Trace (.ab1)…', dis: !has, act: openTraces },
        { label: 'Open Example Plasmid', act: () => addDoc(DNA.buildDemo()) },
        '-',
        { label: 'Save', key: 'Ctrl+S', dis: !has, act: () => saveDoc(false) },
        { label: 'Save As…', key: 'Ctrl+Shift+S', dis: !has, act: () => saveDoc(true) },
        { label: 'Export FASTA…', dis: !has, act: exportFasta },
        { label: 'Export Map Image (SVG)…', dis: !has, act: exportSVG },
        '-',
        { label: 'Print…', key: 'Ctrl+P', dis: !has, act: dlgPrint },
        { label: 'Export PDF…', dis: !has, act: () => dlgPrint(true) },
        '-',
        { label: 'Close File', key: 'Ctrl+W', dis: !has, act: () => closeDoc(S.active) },
        ...(host ? ['-', { label: 'Exit', act: () => host.quit() }] : [])
      ] },
      { label: 'Edit', items: [
        { label: 'Undo', key: 'Ctrl+Z', dis: !has || d.hi <= 0, act: undo },
        { label: 'Redo', key: 'Ctrl+Y', dis: !has || d.hi >= d.hist.length - 1, act: redo },
        '-',
        { label: 'Copy', key: 'Ctrl+C', dis: !sl, act: () => copyText(selSeq(D()), 'Copied ' + fmt(selLen(D())) + ' bp') },
        { label: 'Copy Reverse Complement', dis: !sl, act: () => copyText(DNA.rc(selSeq(D())), 'Copied reverse complement') },
        { label: 'Copy Protein Translation', dis: !sl, act: () => copyText(DNA.translate(selSeq(D())), 'Copied translation') },
        { label: 'Paste', key: 'Ctrl+V', dis: !has, act: pasteFromClipboard },
        { label: 'Delete Selection', key: 'Del', dis: !sl, act: deleteSel },
        '-',
        { label: 'Insert Bases…', dis: !has, act: dlgInsert },
        { label: 'Select All', key: 'Ctrl+A', dis: !has, act: () => setSel(0, D().seq.length, { anchor: 0 }) },
        { label: 'Select Range…', dis: !has, act: dlgGoto },
        { label: 'Go To Position…', key: 'Ctrl+G', dis: !has, act: dlgGoto },
        '-',
        { label: 'Find…', key: 'Ctrl+F', dis: !has, act: openFind },
        { label: 'Rename Sequence…', dis: !has, act: dlgRename }
      ] },
      { label: 'View', items: [
        { label: 'Map', key: 'Ctrl+1', chk: S.view === 'map', dis: !has, act: () => setView('map') },
        { label: 'Sequence', key: 'Ctrl+2', chk: S.view === 'seq', dis: !has, act: () => setView('seq') },
        { label: 'Enzymes', key: 'Ctrl+3', chk: S.view === 'enz', dis: !has, act: () => setView('enz') },
        { label: 'Features', key: 'Ctrl+4', chk: S.view === 'feat', dis: !has, act: () => setView('feat') },
        { label: 'Primers', key: 'Ctrl+5', chk: S.view === 'primer', dis: !has, act: () => setView('primer') },
        { label: 'History', key: 'Ctrl+6', chk: S.view === 'hist', dis: !has, act: () => setView('hist') },
        '-',
        { label: 'Circular Map', chk: has && d.opts.mapMode === 'circular', dis: !has, act: () => setOpt('mapMode', 'circular') },
        { label: 'Linear Map', chk: has && d.opts.mapMode === 'linear', dis: !has, act: () => setOpt('mapMode', 'linear') },
        '-',
        { label: 'Show Features', chk: has && d.opts.showFeat, dis: !has, act: () => setOpt('showFeat', !d.opts.showFeat) },
        { label: 'Show Primers', chk: has && d.opts.showPrimers, dis: !has, act: () => setOpt('showPrimers', !d.opts.showPrimers) },
        { label: 'Show ORFs', chk: has && d.opts.showORFs, dis: !has, act: () => setOpt('showORFs', !d.opts.showORFs) },
        { label: 'Show Translations', chk: has && d.opts.showTrans, dis: !has, act: () => setOpt('showTrans', !d.opts.showTrans) },
        { label: 'Show Complementary Strand', chk: has && d.opts.showComp, dis: !has, act: () => setOpt('showComp', !d.opts.showComp) },
        '-',
        { label: 'Show Side Panel', chk: !S.sideHidden, act: () => { S.sideHidden = !S.sideHidden; renderAll(); } },
        { label: 'Zoom In', key: 'Ctrl++', dis: !has, act: () => zoom(1.25) },
        { label: 'Zoom Out', key: 'Ctrl+-', dis: !has, act: () => zoom(0.8) },
        { label: 'Reset Zoom', key: 'Ctrl+0', dis: !has, act: () => { S.zoom = 1; renderView(); } }
      ] },
      { label: 'Features', items: [
        { label: 'Add Feature…', key: 'Ctrl+K', dis: !has, act: () => dlgFeature(null) },
        { label: 'Detect Common Features', dis: !has, act: detectCommon },
        { label: 'Show Features Table', dis: !has, act: () => setView('feat') },
        '-',
        { label: 'Add ORF as Feature', dis: !sl, act: selToCDS }
      ] },
      { label: 'Primers', items: [
        { label: 'Add Primer…', dis: !has, act: () => dlgPrimer(null) },
        { label: 'Primer from Selection (Top Strand)', dis: !sl, act: () => primerFromSel(1) },
        { label: 'Primer from Selection (Bottom Strand)', dis: !sl, act: () => primerFromSel(-1) },
        '-',
        { label: 'Show Primers Table', dis: !has, act: () => setView('primer') }
      ] },
      { label: 'Cloning', items: [
        { label: 'Restriction Cloning…', dis: S.docs.length < 1, act: dlgRestrictionClone },
        { label: 'Golden Gate Assembly…', dis: S.docs.length < 1, act: dlgGoldenGate },
        { label: 'Gibson / Overlap Assembly…', dis: S.docs.length < 1, act: dlgGibson },
        { label: 'Gateway BP / LR…', dis: S.docs.length < 2, act: dlgGateway },
        '-',
        { label: 'Insert Fragment from Clipboard…', dis: !has, act: dlgInsert },
        { label: 'Show att Sites in This Sequence', dis: !has, act: showAttSites }
      ] },
      { label: 'Tools', items: [
        { label: 'Codon Optimisation…', dis: !has, act: dlgCodon },
        { label: 'Compare Sequences…', dis: S.docs.length < 2, act: dlgCompare },
        '-',
        { label: 'Design PCR Primers…', dis: !has, act: dlgDesignPCR },
        { label: 'Design Sequencing Primers…', dis: !has, act: dlgDesignWalk },
        { label: 'Design Mutagenesis Primers…', dis: !has || !sl, act: dlgDesignMut },
        '-',
        { label: 'Align Sequencing Traces…', dis: !has, act: openTraces },
        { label: 'Show Alignments', dis: !has, act: () => setView('align') },
        { label: 'Remove All Alignments', dis: !has || !d.alignments.length, act: () => { d.alignments = []; toast('Alignments removed'); renderAll(); } }
      ] },
      { label: 'Analyze', items: [
        { label: 'Choose Enzymes…', dis: !has, act: () => dlgEnzymeChooser() },
        { label: 'Enzyme Database…', act: dlgEnzymeDB },
        { label: 'Feature Database…', act: dlgFeatureDB },
        '-',
        { label: 'Primer Structures…', dis: !has || !d.primers.length, act: () => dlgPrimerStruct(null, null) },
        { label: 'RNA Secondary Structure…', dis: !has, act: dlgRNA },
        { label: 'Protein Properties…', dis: !has, act: dlgProtein },
        '-',
        { label: 'Restriction Digest & Gel…', dis: !has, act: dlgDigest },
        { label: 'New Agarose Gel…', dis: !S.docs.some(x => !isGel(x)), act: () => newGelDoc({ doc: has ? d : null }) },
        { label: 'Translate Selection…', dis: !sl, act: dlgTranslate }
      ] },
      { label: 'Actions', items: [
        { label: 'Restriction Digest & Gel…', dis: !has, act: dlgDigest },
        { label: 'Translate Selection…', dis: !sl, act: dlgTranslate },
        '-',
        { label: 'Reverse Complement Selection', dis: !has, act: revCompSel },
        { label: 'Reverse Complement Entire Sequence', dis: !has, act: revCompAll },
        { label: 'Set Origin at Cursor', dis: !has || !d.circular, act: setOrigin },
        '-',
        { label: 'Make Circular', chk: has && d.circular, dis: !has, act: () => setTopology(true) },
        { label: 'Make Linear', chk: has && !d.circular, dis: !has, act: () => setTopology(false) }
      ] },
      { label: 'Help', items: [
        { label: 'Keyboard Shortcuts', act: dlgShortcuts },
        { label: 'About GeneCanvas', act: dlgAbout }
      ] }
    ];
  }
  let openMenu = -1;
  function renderMenuBar() {
    const mb = $('#menubar');
    mb.innerHTML = menus().map((m, i) => `<div class="mi" data-i="${i}">${m.label}</div>`).join('') + `<div class="brand"><span class="brand-dot"></span>GeneCanvas</div>`;
    mb.onmousedown = e => { const it = e.target.closest('.mi'); if (!it) return; e.preventDefault(); const i = +it.dataset.i; if (openMenu === i) closeMenu(); else showMenu(i); };
    mb.onmouseover = e => { const it = e.target.closest('.mi'); if (it && openMenu >= 0 && +it.dataset.i !== openMenu) showMenu(+it.dataset.i); };
  }
  function renderMenuState() { /* menu items are generated when opened */ }
  function showMenu(i) {
    closeCtx();
    const m = menus()[i]; openMenu = i;
    document.querySelectorAll('.menubar .mi').forEach(el => el.classList.toggle('open', +el.dataset.i === i));
    const anchor = document.querySelector(`.menubar .mi[data-i="${i}"]`).getBoundingClientRect();
    const pop = $('#menupop');
    pop.innerHTML = popItems(m.items);
    pop.hidden = false;
    pop.style.left = Math.min(anchor.left, window.innerWidth - 270) + 'px'; pop.style.top = anchor.bottom + 'px';
    pop.onmousedown = e => e.preventDefault();
    pop.onclick = e => { const it = e.target.closest('.it'); if (!it || it.classList.contains('dis')) return; const a = m.items[+it.dataset.i]; closeMenu(); a.act(); };
  }
  function popItems(items) {
    return items.map((a, j) => a === '-' ? '<div class="sep"></div>' : `<div class="it${a.dis ? ' dis' : ''}${a.chk ? ' chk' : ''}" data-i="${j}">${esc(a.label)}${a.key ? `<span class="k">${a.key}</span>` : ''}</div>`).join('');
  }
  function closeMenu() { openMenu = -1; $('#menupop').hidden = true; document.querySelectorAll('.menubar .mi.open').forEach(el => el.classList.remove('open')); }

  // context menu
  function showCtx(x, y, items) {
    closeMenu();
    const c = $('#ctx'); c.innerHTML = items.map((a, j) => a === '-' ? '<div class="sep"></div>' : `<div class="it${a.dis ? ' dis' : ''}" data-i="${j}">${esc(a.label)}</div>`).join('');
    c.hidden = false;
    const w = 220, h = items.length * 26;
    c.style.left = Math.min(x, window.innerWidth - w - 8) + 'px'; c.style.top = Math.min(y, window.innerHeight - h - 8) + 'px';
    c.onclick = e => { const it = e.target.closest('.it'); if (!it) return; closeCtx(); items[+it.dataset.i].act(); };
  }
  function closeCtx() { $('#ctx').hidden = true; }
  function selCtxItems(extra = []) {
    const d = D(), sl = selLen(d) > 0;
    return [...extra,
      { label: 'Copy', dis: !sl, act: () => copyText(selSeq(D()), 'Copied ' + fmt(selLen(D())) + ' bp') },
      { label: 'Copy Reverse Complement', dis: !sl, act: () => copyText(DNA.rc(selSeq(D())), 'Copied reverse complement') },
      '-',
      { label: 'Add Feature…', dis: !sl, act: () => dlgFeature(null) },
      { label: 'Primer from Selection (Top)', dis: !sl, act: () => primerFromSel(1) },
      { label: 'Primer from Selection (Bottom)', dis: !sl, act: () => primerFromSel(-1) },
      { label: 'Translate Selection…', dis: !sl, act: dlgTranslate },
      '-',
      { label: 'Insert Bases…', act: dlgInsert },
      { label: 'Delete Selection', dis: !sl, act: deleteSel },
      { label: 'Set Origin Here', dis: !d.circular, act: setOrigin },
      { label: 'Select All', act: () => setSel(0, D().seq.length, { anchor: 0 }) }
    ];
  }
  function itemCtx(d, key) {
    const it = findItem(d, key); if (!it) return [];
    if (it.kind === 'f') return [
      { label: `Edit Feature “${it.name}”…`, act: () => dlgFeature(it.ref) },
      { label: 'Translate Feature…', act: () => { selectRange(it); dlgTranslate(); } },
      { label: 'Delete Feature', act: () => removeFeature(it.ref.id) }, '-'];
    if (it.kind === 'p') return [
      { label: `Edit Primer “${it.name}”…`, act: () => dlgPrimer(it.ref) },
      { label: 'Delete Primer', act: () => removePrimer(it.ref.id) }, '-'];
    if (it.kind === 'o') return [{ label: 'Add ORF as CDS Feature', act: () => { selectRange(it); addFeature({ name: 'ORF ' + it.ref.aa + ' aa', type: 'CDS', start: it.start, end: it.end, strand: it.strand }); } }, '-'];
    return [];
  }

  // ---------------------------------------------------------------- toolbar / tabs
  function renderToolbar() {
    const d = D(), has = !!d && !isGel(d);
    const b = (k, label, act, o = {}) => `<button class="tb${o.on ? ' on' : ''}" data-a="${act}" title="${esc(o.title || label)}" ${o.dis ? 'disabled' : ''}>${icon(k)}<span>${label}</span></button>`;
    $('#toolbar').innerHTML = [
      b('new', 'New', 'new'), b('open', 'Open', 'open'), b('ncbi', 'NCBI', 'ncbi', { title: 'Download a sequence from NCBI' }), b('save', 'Save', 'save', { dis: !has }), b('print', 'Print', 'print', { dis: !has }),
      '<div class="tbsep"></div>',
      b('undo', 'Undo', 'undo', { dis: !has || d.hi <= 0 }), b('redo', 'Redo', 'redo', { dis: !has || d.hi >= d.hist.length - 1 }),
      '<div class="tbsep"></div>',
      b('find', 'Find', 'find', { dis: !has }),
      b('circ', 'Circular', 'circ', { dis: !has, on: has && d.opts.mapMode === 'circular', title: 'Circular map' }),
      b('lin', 'Linear', 'lin', { dis: !has, on: has && d.opts.mapMode === 'linear', title: 'Linear map' }),
      '<div class="tbsep"></div>',
      b('feat', 'Feature', 'feat', { dis: !has, title: 'Add feature from selection' }),
      b('primer', 'Primer', 'primer', { dis: !has, title: 'Add primer' }),
      b('enz', 'Enzymes', 'enz', { dis: !has, title: 'Choose enzymes' }),
      b('gel', 'Digest', 'gel', { dis: !has, title: 'Restriction digest and simulated gel' }),
      '<div class="tbsep"></div>',
      b('clone', 'Clone', 'clone', { dis: !has, title: 'Restriction cloning' }),
      b('design', 'Primers', 'design', { dis: !has, title: 'Design PCR primers' }),
      b('align', 'Trace', 'align', { dis: !has, title: 'Align a sequencing trace' }),
      '<div class="tbsep"></div>',
      b('tr', 'Translate', 'tr', { dis: !has || !selLen(d) }),
      b('rc', 'Rev Comp', 'rc', { dis: !has }),
      b('origin', 'Origin', 'origin', { dis: !has || !d.circular, title: 'Set origin at cursor' }),
      '<div class="tbsep"></div>',
      b('side', 'Panel', 'side', { on: !S.sideHidden, title: 'Show or hide the side panel' })
    ].join('');
  }
  const TB_ACT = {
    new: () => dlgNew(), open: () => openFiles(), save: () => saveDoc(false), undo, redo, find: () => openFind(),
    circ: () => setOpt('mapMode', 'circular', true), lin: () => setOpt('mapMode', 'linear', true),
    feat: () => dlgFeature(null), primer: () => (selLen(D()) ? primerFromSel(1) : dlgPrimer(null)), enz: () => setView('enz'), gel: () => dlgDigest(),
    tr: () => dlgTranslate(), rc: () => revCompSel(), origin: () => setOrigin(), side: () => { S.sideHidden = !S.sideHidden; renderAll(); },
    ncbi: () => dlgNCBI(), print: () => dlgPrint(), clone: () => dlgRestrictionClone(), design: () => dlgDesignPCR(), align: () => openTraces()
  };

  function renderDocTabs() {
    $('#doctabs').innerHTML = S.docs.map((d, i) => `<div class="dtab${i === S.active ? ' on' : ''}${isGel(d) ? ' gel' : ''}" data-i="${i}" title="${esc(d.name)}${isGel(d) ? ' \u2014 agarose gel' : ' \u2014 ' + fmt(d.seq.length) + ' bp'}"><span class="nm">${d.dirty ? '<span class="dirty">● </span>' : ''}${esc(d.name)}</span><button class="x" data-x="${i}" title="Close">×</button></div>`).join('');
  }
  const VIEWS = [['map', 'Map', 'map'], ['seq', 'Sequence', 'seq'], ['enz', 'Enzymes', 'enz'], ['feat', 'Features', 'feat'], ['primer', 'Primers', 'primer'], ['align', 'Alignment', 'align'], ['compare', 'Compare', 'compare'], ['hist', 'History', 'hist']];
  function renderViewTabs() {
    const vt = $('#viewtabs');
    vt.classList.toggle('noside', S.sideHidden || !(S.view === 'map' || S.view === 'seq'));
    vt.innerHTML = D() && !isGel(D()) ? VIEWS.map(([k, l, ic]) => `<button class="vt${S.view === k ? ' on' : ''}" data-v="${k}">${icon(ic)}${l}</button>`).join('') : '';
  }
  function setView(v) { if (!D()) return; S.view = v; renderAll(); if (v === 'seq') setTimeout(() => { scrollSeqTo(D().sel.start); $('#seqscroll') && $('#seqscroll').focus(); }, 0); }
  function setOpt(k, v, toggleView) {
    const d = D(); if (!d) return; d.opts[k] = v;
    if (k === 'orfMin') d.cache = null;
    if (toggleView && S.view !== 'map') S.view = 'map';
    renderAll(true);
  }
  function zoom(f) { S.zoom = Math.max(0.5, Math.min(6, S.zoom * f)); if (S.view === 'map') renderView(true); }

  // ---------------------------------------------------------------- side panel
  function renderSide() {
    const side = $('#side'), d = D();
    const show = d && !isGel(d) && !S.sideHidden && (S.view === 'map' || S.view === 'seq');
    side.classList.toggle('hidden', !show);
    if (!show) return;
    const o = d.opts, n = d.seq.length, L = selLen(d), A = analyze(d);
    const chk = (k, label) => `<label><input type="checkbox" data-o="${k}" ${o[k] ? 'checked' : ''}> ${label}</label>`;
    let h = `<h4>Map</h4><div class="seg"><button data-mm="circular" class="${o.mapMode === 'circular' ? 'on' : ''}">Circular</button><button data-mm="linear" class="${o.mapMode === 'linear' ? 'on' : ''}">Linear</button></div>`;
    h += `<h4>Enzymes</h4><select id="side-enz">
      <option value="unique6" ${o.enzMode === 'unique6' ? 'selected' : ''}>Unique 6+ cutters</option>
      <option value="unique" ${o.enzMode === 'unique' ? 'selected' : ''}>All unique cutters</option>
      <option value="dual6" ${o.enzMode === 'dual6' ? 'selected' : ''}>Unique &amp; dual 6+ cutters</option>
      <option value="chosen" ${o.enzMode === 'chosen' ? 'selected' : ''}>Chosen enzymes (${d.chosen.length})</option>
      <option value="none" ${o.enzMode === 'none' ? 'selected' : ''}>None</option></select>
      <button class="btn" id="side-pick" style="width:100%;margin-top:4px">Choose Enzymes…</button>`;
    h += `<h4>Show</h4>${chk('showFeat', 'Features')}${chk('showPrimers', 'Primers')}${chk('showORFs', 'ORFs')}
      <div class="row2" style="padding-left:22px">min <input type="number" id="side-orfmin" min="20" max="2000" step="10" value="${o.orfMin}"> aa</div>
      ${chk('showTrans', 'Translations')}${chk('showComp', 'Complementary strand')}`;
    h += `<h4>Selection</h4>`;
    if (L) {
      const s = selSeq(d), a = d.sel.start + 1, b = d.sel.end === 0 ? n : d.sel.end;
      h += `<dl class="kv"><dt>Range</dt><dd>${fmt(a)} – ${fmt(b)}</dd><dt>Length</dt><dd>${fmt(L)} bp</dd><dt>GC</dt><dd>${DNA.gcPct(s).toFixed(1)} %</dd>`;
      if (L <= 80) h += `<dt>Tm</dt><dd>${DNA.tm(s).toFixed(1)} °C</dd>`;
      h += `<dt>MW (ds)</dt><dd>${fmt(Math.round(DNA.dsDnaMW(L)))} Da</dd></dl>`;
      if (L >= 3 && L <= 3000) h += `<div style="margin-top:6px;color:var(--muted);font-size:11.5px">Translation (frame +1)</div><div class="aa-box">${esc(DNA.translate(s))}</div>`;
    } else {
      h += `<dl class="kv"><dt>Cursor</dt><dd>after ${fmt(d.sel.start)}</dd></dl><div style="color:var(--faint);font-size:11.5px;margin-top:4px">Drag on the map or sequence to select. Click a feature or enzyme to select it.</div>`;
    }
    h += `<h4>Sequence</h4><dl class="kv"><dt>Length</dt><dd>${fmt(n)} bp</dd><dt>Topology</dt><dd>${d.circular ? 'circular' : 'linear'}</dd><dt>GC</dt><dd>${DNA.gcPct(d.seq).toFixed(1)} %</dd><dt>Features</dt><dd>${d.features.length}</dd><dt>Primers</dt><dd>${d.primers.length}</dd><dt>Enzymes cutting</dt><dd>${Object.keys(A.sites).length}</dd><dt>ORFs ≥ ${o.orfMin} aa</dt><dd>${A.orfs.length}</dd></dl>`;
    side.innerHTML = h;
  }

  function renderStatus() {
    const d = D(), st = $('#status');
    if (!d) { st.innerHTML = '<span>No file open</span>'; return; }
    if (isGel(d)) {
      const o = d.gel, r = GEL.rangeFor(o.percent);
      st.innerHTML = `<span><b>${esc(d.name)}</b></span><span>${o.percent.toFixed(1)}% agarose, 1\u00d7 ${esc(o.buffer)}, run ${Math.round(o.run * 100)}%</span><span>${o.laneCount} lanes</span><span class="r">resolves about ${fmt(Math.round(r[1]))}\u2013${fmt(Math.round(r[0]))} bp</span>`;
      return;
    }
    const n = d.seq.length, L = selLen(d);
    let sel;
    if (L) {
      const s = selSeq(d);
      sel = `Selected <b>${fmt(d.sel.start + 1)} – ${fmt(d.sel.end === 0 ? n : d.sel.end)}</b> = <b>${fmt(L)} bp</b> &nbsp; GC ${DNA.gcPct(s).toFixed(0)}%` + (L <= 80 ? ` &nbsp; Tm ${DNA.tm(s).toFixed(1)} °C` : '');
    } else sel = `Cursor after base <b>${fmt(d.sel.start)}</b>`;
    st.innerHTML = `<span><b>${esc(d.name)}</b> &nbsp; ${fmt(n)} bp &nbsp; ${d.circular ? 'circular' : 'linear'}</span><span>${sel}</span><span class="r">${S.view === 'seq' ? 'Type A/C/G/T to edit · Shift+arrows extend' : ''}</span>`;
  }

  // ---------------------------------------------------------------- main view dispatcher
  function renderView(keepScroll) {
    const v = $('#view'), d = D();
    if (isGel(d)) { v.onclick = null; renderGel(); return; }
    if (!d) {
      v.innerHTML = `<div class="empty"><h2>GeneCanvas</h2><p>Open a GenBank, FASTA or .dna file, or start a new sequence. You can also drop a file onto this window.</p><div class="btns"><button class="btn pri" data-a="open">Open File…</button><button class="btn" data-a="new">New DNA File…</button><button class="btn" data-a="demo">Open Example Plasmid</button></div></div>`;
      v.onclick = e => { const b = e.target.closest('[data-a]'); if (!b) return; ({ open: openFiles, new: dlgNew, demo: () => addDoc(DNA.buildDemo()) })[b.dataset.a](); };
      return;
    }
    v.onclick = null;
    if (S.view === 'map') renderMap(keepScroll);
    else if (S.view === 'seq') renderSeq(keepScroll);
    else if (S.view === 'enz') renderEnzymes();
    else if (S.view === 'feat') renderFeatures();
    else if (S.view === 'primer') renderPrimers();
    else if (S.view === 'align') renderAlign();
    else if (S.view === 'compare') renderCompare();
    else if (S.view === 'hist') renderHistory();
  }

  // ================================================================ MAP VIEW
  // Where the middle of the window sits, as a fraction of the whole drawing, so the same
  // place stays in view when the drawing is redrawn at another size.
  function scrollCentre(el) {
    const mx = el.scrollWidth - el.clientWidth, my = el.scrollHeight - el.clientHeight;
    return { x: mx > 0 ? el.scrollLeft / mx : 0.5, y: my > 0 ? el.scrollTop / my : 0.5 };
  }
  function applyCentre(el, c) {
    const mx = el.scrollWidth - el.clientWidth, my = el.scrollHeight - el.clientHeight;
    el.scrollLeft = mx > 0 ? Math.round(mx * c.x) : 0;
    el.scrollTop = my > 0 ? Math.round(my * c.y) : 0;
  }
  function renderMap(keepScroll) {
    const v = $('#view');
    let wrap = $('#mapwrap');
    const prev = wrap ? scrollCentre(wrap) : null;
    v.innerHTML = `<div class="mapwrap" id="mapwrap"></div>${mapSideHTML()}<div class="map-zoom"><button data-z="out" title="Zoom out">−</button><button data-z="reset" title="Fit">⤢</button><button data-z="in" title="Zoom in">+</button></div>`;
    wrap = $('#mapwrap');
    const d = D(), A = analyze(d);
    const W = Math.max(420, wrap.clientWidth), H = Math.max(380, wrap.clientHeight);
    const svg = d.opts.mapMode === 'circular' ? drawCircular(d, A, W * S.zoom, H * S.zoom) : drawLinear(d, A, W * S.zoom);
    wrap.innerHTML = svg;
    // The drawing is normally wider or taller than the window, so centre it rather than
    // leaving it pinned to the top left corner.
    applyCentre(wrap, keepScroll && prev ? prev : { x: 0.5, y: 0.5 });
    drawSelLayer();
    bindMap(wrap.querySelector('svg'));
    v.querySelector('.map-zoom').onclick = e => { const b = e.target.closest('button'); if (!b) return; if (b.dataset.z === 'in') zoom(1.25); else if (b.dataset.z === 'out') zoom(0.8); else { S.zoom = 1; renderView(true); } };
    bindMapSide();
  }

  // ---------------------------------------------------------------- map view panel
  const MS_ICONS = {
    feat: '<path d="M3 9h12l4 3-4 3H3z"/>',
    primer: '<path d="M3 12h13"/><path d="m12 8 4 4-4 4"/>',
    enz: '<path d="M6 4l12 16M18 4 6 20"/>',
    orf: '<path d="M4 12h12"/><path d="M16 9l4 3-4 3"/><circle cx="6" cy="12" r="2"/>',
    label: '<rect x="3" y="7" width="18" height="10" rx="2"/><path d="M7 12h10"/>',
    online: '<circle cx="12" cy="12" r="8"/><path d="M4 12h16" stroke-dasharray="2 2"/>',
    ruler: '<rect x="3" y="9" width="18" height="6" rx="1"/><path d="M7 9v3M11 9v3M15 9v3M19 9v3"/>',
    paint: '<path d="M6 3h9v6H6z"/><path d="M12 9v4a2 2 0 0 1-2 2H9v6"/>'
  };
  function mapSideHTML() {
    const d = D(); if (!d) return '';
    const o = d.opts;
    const b = (k, on, title, act) => `<button class="msb${on ? ' on' : ''}" data-ms="${act}" title="${esc(title)}"><svg viewBox="0 0 24 24">${MS_ICONS[k]}</svg></button>`;
    return `<div class="mapside" id="mapside">
      ${b('feat', o.showFeat, 'Show features', 'showFeat')}
      ${b('primer', o.showPrimers, 'Show primers', 'showPrimers')}
      ${b('enz', o.enzMode !== 'none', 'Show enzyme sites', 'enz')}
      ${b('orf', o.showORFs, 'Show open reading frames', 'showORFs')}
      <div class="mssep"></div>
      ${b('label', o.labelMode !== 'auto', `Feature labels: ${o.labelMode}`, 'label')}
      ${b('online', o.featOnLine, 'Draw features on the backbone', 'featOnLine')}
      ${b('ruler', o.showRuler, 'Show the position ruler', 'showRuler')}
      <div class="mssep"></div>
      ${b('paint', o.baseColors.enabled, 'Colour the bases', 'colors')}
    </div>`;
  }
  function bindMapSide() {
    const el = $('#mapside'); if (!el) return;
    el.onclick = e => {
      const btn = e.target.closest('[data-ms]'); if (!btn) return;
      const d = D(), o = d.opts, a = btn.dataset.ms;
      if (a === 'enz') o.enzMode = o.enzMode === 'none' ? 'unique6' : 'none';
      else if (a === 'label') o.labelMode = o.labelMode === 'auto' ? 'inside' : (o.labelMode === 'inside' ? 'outside' : 'auto');
      else if (a === 'colors') { dlgBaseColors(); return; }
      else o[a] = !o[a];
      d.cache = null;
      renderAll(true);
      if (a === 'label') toast('Feature labels: ' + d.opts.labelMode);
    };
  }
  function dlgBaseColors() {
    const d = D(), c = d.opts.baseColors;
    dialog({
      title: 'Base Colours',
      body: `<div class="form">
        <label></label><label><input type="checkbox" id="bcOn" ${c.enabled ? 'checked' : ''}> Colour each base in the sequence view</label>
        ${['A', 'C', 'G', 'T'].map(b => `<label for="bc${b}">${b}</label><input type="color" id="bc${b}" value="${c[b]}">`).join('')}
        <div class="full hint">Colours apply to both strands in the Sequence view. They are saved with the program, so every file you open uses them.</div></div>`,
      buttons: [{ label: 'Reset', act: () => { Object.assign(d.opts.baseColors, { A: '#2e9e4f', C: '#2563c9', G: '#1c2230', T: '#c43c3c' }); SET.baseColors = { ...d.opts.baseColors }; saveSettings(); renderAll(true); } },
      { label: 'Cancel' }, {
        label: 'Apply', primary: true, act: () => {
          c.enabled = $('#bcOn').checked;
          for (const b of ['A', 'C', 'G', 'T']) c[b] = val('bc' + b);
          SET.baseColors = { ...c };
          saveSettings();
          renderAll(true);
        }
      }]
    });
  }

  function arcPath(cx, cy, r, a1, a2, sweep) {
    const span = Math.abs(a2 - a1);
    const p1 = [cx + r * Math.cos(a1), cy + r * Math.sin(a1)], p2 = [cx + r * Math.cos(a2), cy + r * Math.sin(a2)];
    return `M${p1[0].toFixed(2)},${p1[1].toFixed(2)} A${r},${r} 0 ${span > Math.PI ? 1 : 0} ${sweep} ${p2[0].toFixed(2)},${p2[1].toFixed(2)}`;
  }

  function drawCircular(d, A, W, H) {
    const n = d.seq.length, o = d.opts;
    const items = collectItems(d, A);
    const feats = items.filter(i => i.kind === 'f'), prims = items.filter(i => i.kind === 'p'), orfs = items.filter(i => i.kind === 'o');
    let nF = assignLanes(feats, n), nP = assignLanes(prims, n), nO = assignLanes(orfs, n);
    if (o.featOnLine) { for (const f of feats) f.lane = 0; nF = 1; }
    const enz = shownEnz(d, A);
    // label count estimation decides how much room labels need
    const labelGuess = enz.reduce((s, e) => s + e.sites.length, 0) + feats.length;
    const needH = Math.ceil(labelGuess / 2) * 14 + 60;
    H = Math.max(H, Math.min(needH, 2400));
    const cx = W / 2, cy = H / 2;
    const laneW = 14, laneGap = 5;
    const pBand = nP ? nP * 9 + 8 : 0;
    const R = Math.max(70, Math.min(W / 2 - 175, H / 2 - 60) - pBand);
    const ang = bp => bp / n * TAU - Math.PI / 2;
    const P = (r, a) => [cx + r * Math.cos(a), cy + r * Math.sin(a)];
    const fmtP = p => p[0].toFixed(1) + ',' + p[1].toFixed(1);
    const fR = k => R - k * (laneW + laneGap);
    const innerF = fR(Math.max(nF, 1) - 1) - laneW / 2;
    const oR = k => innerF - 7 - k * 8;
    const innerAll = nO ? oR(nO - 1) - 4 : innerF - 2;
    const pR = k => R + laneW / 2 + 7 + k * 9;
    const outer = nP ? pR(nP - 1) + 4 : R + laneW / 2 + 2;
    const Rl = outer + 30;
    const g = [];
    let defs = '';

    // ticks inside
    const step = niceStep(n / 10);
    const tr = Math.max(20, innerAll - 3);
    for (let t = step; o.showRuler && t < n; t += step) {
      const a = ang(t);
      g.push(`<line class="tick" x1="${fmtP(P(tr, a)).replace(',', '" y1="')}" x2="${fmtP(P(tr - 5, a)).replace(',', '" y2="')}"/>`);
      if (tr > 70) { const q = P(tr - 15, a); g.push(`<text class="tickl" x="${q[0].toFixed(1)}" y="${q[1].toFixed(1)}" text-anchor="middle" dy=".35em">${fmt(t)}</text>`); }
    }
    // backbone
    g.push(`<circle class="bb" cx="${cx}" cy="${cy}" r="${R - 2}"/><circle class="bb" cx="${cx}" cy="${cy}" r="${R + 2}"/>`);
    // origin mark
    g.push(`<line class="ecut" x1="${cx}" y1="${cy - R - 6}" x2="${cx}" y2="${cy - R + 6}" stroke-width="1.6"/>`);

    const labels = [];
    // arrows
    const arrow = (it, r, w, headPx) => {
      const a1 = ang(it.start);
      let span = DNA.rlen(it, n) / n * TAU; if (span >= TAU) span = TAU - 0.0005;
      const a2 = a1 + span;
      const ro = r + w / 2, ri = r - w / 2;
      const h = it.strand ? Math.min(headPx / r, span * 0.6) : 0;
      const lg = s => s > Math.PI ? 1 : 0;
      let path;
      if (it.strand === 1) {
        const b = a2 - h;
        path = `M${fmtP(P(ro, a1))} A${ro},${ro} 0 ${lg(b - a1)} 1 ${fmtP(P(ro, b))} L${fmtP(P(r, a2))} L${fmtP(P(ri, b))} A${ri},${ri} 0 ${lg(b - a1)} 0 ${fmtP(P(ri, a1))} Z`;
      } else if (it.strand === -1) {
        const b = a1 + h;
        path = `M${fmtP(P(r, a1))} L${fmtP(P(ro, b))} A${ro},${ro} 0 ${lg(a2 - b)} 1 ${fmtP(P(ro, a2))} L${fmtP(P(ri, a2))} A${ri},${ri} 0 ${lg(a2 - b)} 0 ${fmtP(P(ri, b))} Z`;
      } else {
        path = `M${fmtP(P(ro, a1))} A${ro},${ro} 0 ${lg(span)} 1 ${fmtP(P(ro, a2))} L${fmtP(P(ri, a2))} A${ri},${ri} 0 ${lg(span)} 0 ${fmtP(P(ri, a1))} Z`;
      }
      return { path, a1, a2, h, span };
    };
    let tpi = 0;
    for (const it of feats) {
      const r = fR(it.lane);
      const ar = arrow(it, r, laneW, 11);
      g.push(`<path class="feat" d="${ar.path}" fill="${it.color}" data-k="${it.key}"><title>${esc(it.name)} (${esc(it.type)}) ${fmt(it.start + 1)}..${fmt(it.end)} · ${fmt(DNA.rlen(it, n))} bp</title></path>`);
      const textW = it.name.length * 6.4 + 10;
      const avail = (ar.span - ar.h) * r;
      const fits = avail > textW;
      const inside = o.labelMode === 'inside' ? avail > it.name.length * 4 : (o.labelMode === 'outside' ? false : fits);
      if (inside) {
        const mid = ar.a1 + ar.span / 2;
        const b1 = it.strand === -1 ? ar.a1 + ar.h : ar.a1, b2 = it.strand === 1 ? ar.a2 - ar.h : ar.a2;
        const id = 'tp' + (tpi++);
        const bottom = Math.sin(mid) > 0.05;
        defs += `<path id="${id}" d="${bottom ? arcPath(cx, cy, r, b2, b1, 0) : arcPath(cx, cy, r, b1, b2, 1)}"/>`;
        g.push(`<text class="ftxt" dy=".35em"><textPath href="#${id}" startOffset="50%" text-anchor="middle">${esc(it.name)}</textPath></text>`);
      } else {
        labels.push({ a: ang(it.start) + DNA.rlen(it, n) / n * Math.PI, text: it.name, cls: 'f', k: it.key, r0: r + laneW / 2 });
      }
    }
    for (const it of orfs) {
      const ar = arrow(it, oR(it.lane), 4, 7);
      g.push(`<path class="orf" d="${ar.path}" fill="${it.color}" data-k="${it.key}"><title>${esc(it.name)} ${fmt(it.start + 1)}..${fmt(it.end)} (${it.strand > 0 ? 'forward' : 'reverse'})</title></path>`);
    }
    for (const it of prims) {
      const ar = arrow(it, pR(it.lane), 5, 7);
      g.push(`<path class="prim" d="${ar.path}" fill="${it.color}" data-k="${it.key}"><title>Primer ${esc(it.name)} ${fmt(it.start + 1)}..${fmt(it.end)}</title></path>`);
      labels.push({ a: ar.a1 + ar.span / 2, text: it.name, cls: 'p', k: it.key, r0: pR(it.lane) + 3 });
    }
    // enzyme cut marks
    for (const e of enz) e.sites.forEach((s, i) => {
      if (s.cut == null) return;
      const a = ang(s.cut);
      g.push(`<line class="ecut" x1="${fmtP(P(R - 7, a)).replace(',', '" y1="')}" x2="${fmtP(P(R + 7, a)).replace(',', '" y2="')}"/>`);
      labels.push({ a, text: `${e.name} (${fmt(s.cut)})`, cls: 'e' + (e.unique ? ' u' : ''), k: `e:${e.name}:${i}`, r0: Math.max(R + 7, outer) });
    });
    // label layout
    const sp = 14;
    const sides = [[], []];
    for (const l of labels) { l.side = Math.cos(l.a) >= -1e-6 ? 0 : 1; l.ty = cy + Rl * Math.sin(l.a); sides[l.side].push(l); }
    for (const arr of sides) {
      arr.sort((a, b) => a.ty - b.ty);
      for (let i = 1; i < arr.length; i++) if (arr[i].ty < arr[i - 1].ty + sp) arr[i].ty = arr[i - 1].ty + sp;
      const maxY = H - 10;
      if (arr.length && arr[arr.length - 1].ty > maxY) {
        arr[arr.length - 1].ty = maxY;
        for (let i = arr.length - 2; i >= 0; i--) if (arr[i].ty > arr[i + 1].ty - sp) arr[i].ty = arr[i + 1].ty - sp;
      }
      if (arr.length && arr[0].ty < 12) { arr[0].ty = 12; for (let i = 1; i < arr.length; i++) if (arr[i].ty < arr[i - 1].ty + sp) arr[i].ty = arr[i - 1].ty + sp; }
    }
    for (const l of labels) {
      const dy = l.ty - cy;
      let dx = Math.abs(dy) < Rl ? Math.sqrt(Rl * Rl - dy * dy) : 0;
      dx = Math.max(dx, 28) + 6;
      const tx = l.side === 0 ? cx + dx : cx - dx;
      const e1 = P(l.r0 + 1, l.a), e2 = P(Rl - 10, l.a);
      g.push(`<polyline class="lead" points="${fmtP(e1)} ${fmtP(e2)} ${(l.side === 0 ? tx - 3 : tx + 3).toFixed(1)},${l.ty.toFixed(1)}"/>`);
      g.push(`<text class="lbl ${l.cls}" x="${tx.toFixed(1)}" y="${l.ty.toFixed(1)}" dy=".35em" text-anchor="${l.side === 0 ? 'start' : 'end'}" data-k="${l.k}">${esc(l.text)}</text>`);
    }
    // center
    g.push(`<text class="ctr1" x="${cx}" y="${cy - 4}" text-anchor="middle">${esc(d.name)}</text><text class="ctr2" x="${cx}" y="${cy + 16}" text-anchor="middle">${fmt(n)} bp</text>`);
    S.mapGeom = { type: 'c', cx, cy, n, r1: Math.max(8, innerAll - 24), r2: Rl - 14, R, ang };
    return `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="Segoe UI, Arial, sans-serif"><defs>${defs}</defs><rect width="${W}" height="${H}" style="fill:var(--canvas)" data-bg="1"/><g id="sellayer"></g>${g.join('')}</svg>`;
  }

  function drawLinear(d, A, W) {
    const n = d.seq.length;
    const items = collectItems(d, A);
    const feats = items.filter(i => i.kind === 'f'), prims = items.filter(i => i.kind === 'p'), orfs = items.filter(i => i.kind === 'o');
    const m = 46, w = W - 2 * m;
    const X = bp => m + bp / n * w;
    // lanes by pixel extent, so labels drawn under short bars never collide
    const pxLanes = (list, withLabel) => {
      const ext = it => {
        let lo = Infinity, hi = -Infinity;
        for (const [a, b] of DNA.segs(it, n)) {
          let x1 = X(a), x2 = X(b);
          const tw = it.name.length * 6.3 + 10;
          if (withLabel && x2 - x1 < tw) { const c = (x1 + x2) / 2; x1 = Math.min(x1, c - tw / 2); x2 = Math.max(x2, c + tw / 2); }
          lo = Math.min(lo, x1); hi = Math.max(hi, x2);
        }
        return [lo, hi];
      };
      const lanes = [];
      for (const it of list.slice().sort((p, q) => p.start - q.start)) {
        const [lo, hi] = ext(it);
        let k = 0; while (lanes[k] !== undefined && lanes[k] > lo - 4) k++;
        lanes[k] = hi; it.lane = k;
      }
      return lanes.length;
    };
    const o = d.opts;
    let nF = pxLanes(feats, true), nP = pxLanes(prims, true), nO = pxLanes(orfs, false);
    if (o.featOnLine) { for (const f of feats) f.lane = 0; nF = 1; }
    const enz = shownEnz(d, A);
    // enzyme label rows
    const elabels = [];
    for (const e of enz) e.sites.forEach((s, i) => { if (s.cut != null) elabels.push({ x: X(s.cut), text: `${e.name} (${fmt(s.cut)})`, k: `e:${e.name}:${i}`, u: e.unique }); });
    elabels.sort((a, b) => a.x - b.x);
    const rowEnd = [];
    for (const l of elabels) {
      const wd = l.text.length * 6.1 + 6;
      let r = 0; while (rowEnd[r] !== undefined && rowEnd[r] > l.x - 2) r++;
      rowEnd[r] = l.x + wd; l.row = r;
    }
    const nRows = rowEnd.length;
    const top = 26;
    const baseY = top + nRows * 14 + 22;
    const g = [];
    for (const l of elabels) {
      const y = top + (nRows - 1 - l.row) * 14;
      g.push(`<line class="lead" x1="${l.x.toFixed(1)}" y1="${baseY - 7}" x2="${l.x.toFixed(1)}" y2="${y + 4}"/>`);
      g.push(`<text class="lbl e${l.u ? ' u' : ''}" x="${(l.x + 2).toFixed(1)}" y="${y}" data-k="${l.k}">${esc(l.text)}</text>`);
      g.push(`<line class="ecut" x1="${l.x.toFixed(1)}" y1="${baseY - 6}" x2="${l.x.toFixed(1)}" y2="${baseY + 6}"/>`);
    }
    g.push(`<line class="bb" x1="${m}" y1="${baseY - 2}" x2="${m + w}" y2="${baseY - 2}"/><line class="bb" x1="${m}" y1="${baseY + 2}" x2="${m + w}" y2="${baseY + 2}"/>`);
    const step = niceStep(n / Math.max(4, Math.floor(w / 90)));
    for (let t = 0; o.showRuler && t <= n; t += step) {
      const x = X(t);
      g.push(`<line class="tick" x1="${x.toFixed(1)}" y1="${baseY + 4}" x2="${x.toFixed(1)}" y2="${baseY + 9}"/><text class="tickl" x="${x.toFixed(1)}" y="${baseY + 20}" text-anchor="middle">${fmt(t || 1)}</text>`);
    }
    let y = baseY + 32;
    const laneH = 30;
    const bar = (it, yy, h, labelBelow) => {
      for (const [a, b] of DNA.segs(it, n)) {
        const x1 = X(a), x2 = X(b), hd = it.strand ? Math.min(9, (x2 - x1) * 0.5) : 0;
        const tipR = it.strand === 1 && b === it.end, tipL = it.strand === -1 && a === it.start;
        const y1 = yy, y2 = yy + h, ym = yy + h / 2;
        let pts;
        if (tipR) pts = `${x1},${y1} ${x2 - hd},${y1} ${x2},${ym} ${x2 - hd},${y2} ${x1},${y2}`;
        else if (tipL) pts = `${x1},${ym} ${x1 + hd},${y1} ${x2},${y1} ${x2},${y2} ${x1 + hd},${y2}`;
        else pts = `${x1},${y1} ${x2},${y1} ${x2},${y2} ${x1},${y2}`;
        const cls = it.kind === 'f' ? 'feat' : (it.kind === 'p' ? 'prim' : 'orf');
        g.push(`<polygon class="${cls}" points="${pts}" fill="${it.color}" data-k="${it.key}"><title>${esc(it.name)} ${fmt(it.start + 1)}..${fmt(it.end)}</title></polygon>`);
        if (labelBelow !== false) {
          const tw = it.name.length * 6.3 + 8;
          const fits = x2 - x1 - hd > tw;
          const inside = o.labelMode === 'inside' ? (x2 - x1 > it.name.length * 4) : (o.labelMode === 'outside' ? false : fits);
          if (it.kind === 'f' && inside) g.push(`<text class="ftxt" x="${((x1 + x2) / 2).toFixed(1)}" y="${ym}" dy=".35em" text-anchor="middle">${esc(it.name)}</text>`);
          else g.push(`<text class="lbl ${it.kind === 'f' ? 'f' : 'p'}" x="${((x1 + x2) / 2).toFixed(1)}" y="${y2 + 11}" text-anchor="middle" data-k="${it.key}">${esc(it.name)}</text>`);
        }
      }
    };
    for (const it of feats) bar(it, y + it.lane * laneH, 15);
    y += Math.max(nF, 0) * laneH;
    for (const it of prims) bar(it, y + it.lane * 24, 7);
    y += nP * 24;
    for (const it of orfs) bar(it, y + it.lane * 10, 5, false);
    y += nO * 10;
    const H = Math.max(y + 20, 200);
    g.push(`<text class="ctr1" x="${m}" y="16" style="font-size:14px">${esc(d.name)}</text><text class="ctr2" x="${m + d.name.length * 8.6 + 10}" y="16" style="font-size:12px">${fmt(n)} bp · ${d.circular ? 'circular (shown linear)' : 'linear'}</text>`);
    S.mapGeom = { type: 'l', m, w, n, X, y1: baseY - 12, y2: H - 6 };
    return `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="Segoe UI, Arial, sans-serif"><rect width="${W}" height="${H}" style="fill:var(--canvas)" data-bg="1"/><g id="sellayer"></g>${g.join('')}</svg>`;
  }

  function niceStep(x) {
    const p = Math.pow(10, Math.floor(Math.log10(Math.max(1, x))));
    for (const k of [1, 2, 2.5, 5, 10]) if (k * p >= x) return k * p;
    return 10 * p;
  }

  function drawSelLayer() {
    const layer = document.getElementById('sellayer'), G = S.mapGeom, d = D();
    if (!layer || !G || !d) return;
    const n = d.seq.length, L = selLen(d);
    let h = '';
    if (G.type === 'c') {
      const { cx, cy, r1, r2, ang } = G;
      const P = (r, a) => (cx + r * Math.cos(a)).toFixed(1) + ',' + (cy + r * Math.sin(a)).toFixed(1);
      if (L) {
        if (L >= n) h = `<path class="selarea" fill-rule="evenodd" d="M${cx - r2},${cy} a${r2},${r2} 0 1 0 ${2 * r2},0 a${r2},${r2} 0 1 0 ${-2 * r2},0 M${cx - r1},${cy} a${r1},${r1} 0 1 0 ${2 * r1},0 a${r1},${r1} 0 1 0 ${-2 * r1},0"/>`;
        else {
          const a1 = ang(d.sel.start), a2 = a1 + L / n * TAU, lg = (a2 - a1) > Math.PI ? 1 : 0;
          h = `<path class="selarea" d="M${P(r1, a1)} L${P(r2, a1)} A${r2},${r2} 0 ${lg} 1 ${P(r2, a2)} L${P(r1, a2)} A${r1},${r1} 0 ${lg} 0 ${P(r1, a1)} Z"/>`;
        }
      } else {
        const a = ang(d.sel.start);
        h = `<line class="caret" x1="${P(r1, a).replace(',', '" y1="')}" x2="${P(r2, a).replace(',', '" y2="')}"/>`;
      }
    } else {
      const { X, y1, y2 } = G;
      if (L) for (const [a, b] of DNA.segs(d.sel, n)) h += `<rect class="selarea" x="${X(a).toFixed(1)}" y="${y1}" width="${Math.max(1, X(b) - X(a)).toFixed(1)}" height="${y2 - y1}"/>`;
      else { const x = X(d.sel.start).toFixed(1); h = `<line class="caret" x1="${x}" y1="${y1}" x2="${x}" y2="${y2}"/>`; }
    }
    layer.innerHTML = h;
  }

  function bindMap(svg) {
    if (!svg) return;
    const d = D();
    const toBp = e => {
      const pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY;
      const p = pt.matrixTransform(svg.getScreenCTM().inverse());
      const G = S.mapGeom, n = d.seq.length;
      if (G.type === 'c') {
        let a = Math.atan2(p.y - G.cy, p.x - G.cx) + Math.PI / 2; if (a < 0) a += TAU;
        return Math.round(a / TAU * n) % n;
      }
      return Math.max(0, Math.min(n, Math.round((p.x - G.m) / G.w * n)));
    };
    let drag = null;
    svg.addEventListener('mousedown', e => {
      if (e.button !== 0) return;
      const k = e.target.getAttribute && e.target.getAttribute('data-k');
      if (k) { clickItem(k, e); return; }
      e.preventDefault();
      const bp = toBp(e);
      if (e.shiftKey) { extendTo(bp); return; }
      drag = { anchor: bp, last: bp, cum: 0 };
      setSel(bp, bp, { anchor: bp });
    });
    svg.addEventListener('dblclick', e => {
      const k = e.target.getAttribute && e.target.getAttribute('data-k');
      if (k && k.startsWith('f:')) { const it = findItem(D(), k); if (it) dlgFeature(it.ref); }
      else if (k && k.startsWith('p:')) { const it = findItem(D(), k); if (it) dlgPrimer(it.ref); }
    });
    svg.addEventListener('contextmenu', e => {
      e.preventDefault();
      const k = e.target.getAttribute && e.target.getAttribute('data-k');
      if (k && !selLen(D())) clickItem(k, e);
      showCtx(e.clientX, e.clientY, selCtxItems(k ? itemCtx(D(), k) : []));
    });
    const move = e => {
      if (!drag) return;
      const bp = toBp(e), n = d.seq.length;
      if (S.mapGeom.type === 'c') {
        let delta = bp - drag.last; if (delta > n / 2) delta -= n; if (delta < -n / 2) delta += n;
        drag.cum += delta; drag.last = bp;
        if (Math.abs(drag.cum) >= n) { setSel(drag.anchor, drag.anchor === 0 ? n : drag.anchor); return; }
        if (!d.circular) { setSel(Math.min(drag.anchor, bp), Math.max(drag.anchor, bp)); return; }
        if (drag.cum >= 0) setSel(drag.anchor, bp); else setSel(bp, drag.anchor);
      } else setSel(Math.min(drag.anchor, bp), Math.max(drag.anchor, bp));
    };
    const up = () => { drag = null; window.removeEventListener('mousemove', move); };
    svg.addEventListener('mousedown', () => { if (drag) { window.addEventListener('mousemove', move); window.addEventListener('mouseup', up, { once: true }); } });
  }

  function extendTo(bp) {
    const d = D(); const a = d.anchor;
    setSel(Math.min(a, bp), Math.max(a, bp));
  }

  function clickItem(k, e) {
    const d = D();
    if (k.startsWith('e:')) {
      const [, name, idx] = k.split(':'); const s = analyze(d).sites[name][+idx];
      let end = s.pos + s.len; if (end > d.seq.length) end -= d.seq.length;
      selectRange({ start: s.pos, end });
      toast(`${name}: ${DNA.ENZ_BY_NAME[name].site} — cuts after ${fmt(s.cut)}`);
      return;
    }
    const it = findItem(d, k); if (!it) return;
    selectRange(it);
  }

  // ================================================================ SEQUENCE VIEW
  let CW = 7.8;
  function measureCW() {
    const sp = document.createElement('span');
    sp.style.cssText = 'position:absolute;visibility:hidden;white-space:pre;font:13px/17px ' + getComputedStyle(document.documentElement).getPropertyValue('--f-mono');
    sp.textContent = 'ACGT'.repeat(50); document.body.appendChild(sp);
    CW = sp.getBoundingClientRect().width / 200; sp.remove();
  }

  function renderSeq(keepScroll) {
    const v = $('#view');
    const old = $('#seqscroll'); const prevTop = old ? old.scrollTop : 0;
    const hadFocus = !!old && document.activeElement === old;
    const f = S.find;
    v.innerHTML = `<div class="findbar"><span>Find</span><input type="text" id="findq" placeholder="DNA (IUPAC ok) or protein" value="${esc(f.q)}" spellcheck="false">
      <select id="findmode"><option value="dna" ${f.mode === 'dna' ? 'selected' : ''}>DNA</option><option value="prot" ${f.mode === 'prot' ? 'selected' : ''}>Protein</option></select>
      <button id="findprev" title="Previous (Shift+Enter)">▲</button><button id="findnext" title="Next (Enter)">▼</button><span class="cnt" id="findcnt"></span><span class="sp"></span>
      <button id="gotobtn">Go To…</button></div><div class="seqscroll" id="seqscroll" tabindex="0"><div id="seqrows"></div></div>`;
    buildSeqRows();
    const sc = $('#seqscroll');
    if (keepScroll) sc.scrollTop = prevTop;
    if (hadFocus) sc.focus({ preventScroll: true });
    bindSeq(sc);
    bindFind();
    updateFindCount();
  }

  function buildSeqRows() {
    const d = D(), A = analyze(d), n = d.seq.length, o = d.opts;
    const sc = $('#seqscroll'), host = $('#seqrows');
    const G = 78;
    const bpr = Math.max(10, Math.floor((sc.clientWidth - G - 30) / CW / 10) * 10);
    // enzyme cuts
    const cuts = [];
    for (const e of shownEnz(d, A)) e.sites.forEach((s, i) => { if (s.cut != null) cuts.push({ pos: s.cut, name: e.name, k: `e:${e.name}:${i}`, u: e.unique }); });
    cuts.sort((a, b) => a.pos - b.pos);
    // feature segments
    const items = collectItems(d, A);
    const segsAll = [];
    for (const it of items) for (const [a, b] of DNA.segs(it, n)) segsAll.push({ a, b, it, tipR: it.strand === 1 && b === it.end, tipL: it.strand === -1 && a === it.start });
    // translations, carrying the residue number so each one can be numbered like a protein
    const aaRows = {};
    if (o.showTrans) for (const it of items) {
      if (!((it.kind === 'f' && it.type === 'CDS') || it.kind === 'o')) continue;
      const L = DNA.rlen(it, n); const fs = DNA.sub(d.seq, it.start, it.end);
      const s = it.strand === -1 ? DNA.rc(fs) : fs;
      for (let j = 0; j + 3 <= L; j += 3) {
        const aa = DNA.CODON[s.substr(j, 3)] || 'X';
        const mid = it.strand === -1 ? ((it.end - j - 2) % n + n) % n : (it.start + j + 1) % n;
        const r = Math.floor(mid / bpr);
        (aaRows[r] = aaRows[r] || []).push({ mid, aa, key: it.key, res: j / 3 + 1 });
      }
    }
    const rows = [], html = [];
    const comp = o.showComp;
    const NT_H = 17, RULER_H = 15, NUM_H = 13, AA_H = 15;
    for (let rs = 0, ri = 0; rs < n; rs += bpr, ri++) {
      const re = Math.min(n, rs + bpr), len = re - rs;
      const parts = [];
      // --- enzyme labels stacked above the strands
      const lineEnd = [];
      const myCuts = cuts.filter(c => c.pos >= rs && (c.pos < re || (re === n && c.pos === n)));
      for (const c of myCuts) {
        const x = (c.pos - rs) * CW, w = c.name.length * 6.6 + 6;
        let l = 0; while (lineEnd[l] !== undefined && lineEnd[l] > x - 1) l++;
        lineEnd[l] = x + w; c.line = l;
      }
      const nL = lineEnd.length, enzH = nL ? nL * 14 + 6 : 3;
      const y0 = 4 + enzH;
      for (const c of myCuts) {
        const x = G + (c.pos - rs) * CW, y = 4 + (nL - 1 - c.line) * 14;
        parts.push(`<div class="enz${c.u ? ' u' : ''}" data-k="${c.k}" style="left:${x - 1}px;top:${y}px">${esc(c.name)}</div><div class="etick" style="left:${x}px;top:${y + 13}px;height:${y0 - y - 12}px"></div>`);
      }
      // --- top strand, tick ruler, bottom strand
      parts.push(`<div class="num" style="top:${y0 + 1}px">${fmt(rs + 1)}</div>`);
      parts.push(`<div class="nt" style="left:${G}px;top:${y0}px">${paintBases(d.seq.slice(rs, re), o)}</div>`);
      const rulerTop = y0 + NT_H;
      let ruler = '';
      for (let p = rs; p < re; p++) ruler += ((p + 1) % 10 === 0) ? '+' : '-';
      parts.push(`<div class="tick" style="left:${G}px;top:${rulerTop}px">${ruler}</div>`);
      let y = rulerTop + RULER_H;
      if (comp) { parts.push(`<div class="nt b" style="left:${G}px;top:${y}px">${paintBases(DNA.comp(d.seq.slice(rs, re)), o)}</div>`); y += NT_H; }
      const strandBottom = y;
      y += 2;
      // --- features, each translated one carrying its own residue numbering above it
      const rsegs = segsAll.filter(s => s.a < re && s.b > rs).map(s => ({ ...s, x1: Math.max(s.a, rs), x2: Math.min(s.b, re) })).sort((p, q) => p.x1 - q.x1 || (q.x2 - q.x1) - (p.x2 - p.x1));
      const laneEnd = [], laneTrans = [];
      for (const s of rsegs) {
        let l = 0; while (laneEnd[l] !== undefined && laneEnd[l] > s.x1) l++;
        laneEnd[l] = s.x2; s.lane = l;
        if (o.showTrans && ((s.it.kind === 'f' && s.it.type === 'CDS') || s.it.kind === 'o')) laneTrans[l] = true;
      }
      const laneBar = [], laneAA = [], laneNum = [];
      for (let l = 0; l < laneEnd.length; l++) {
        if (laneTrans[l]) { laneNum[l] = y; laneAA[l] = y + NUM_H; laneBar[l] = y + NUM_H + AA_H; y += NUM_H + AA_H + 16; }
        else { laneBar[l] = y; y += 17; }
      }
      const keyLane = {};
      for (const s of rsegs) {
        const x = G + (s.x1 - rs) * CW, w = (s.x2 - s.x1) * CW;
        let clip = '';
        if (s.tipR && s.x2 === s.b) clip = 'clip-path:polygon(0 0,calc(100% - 8px) 0,100% 50%,calc(100% - 8px) 100%,0 100%);';
        else if (s.tipL && s.x1 === s.a) clip = 'clip-path:polygon(8px 0,100% 0,100% 100%,8px 100%,0 50%);';
        const label = w > s.it.name.length * 6.4 + 16 ? esc(s.it.name) : '';
        const h = s.it.kind === 'o' ? 8 : 14;
        parts.push(`<div class="fbar${s.it.kind === 'o' ? ' orf' : ''}" data-k="${s.it.key}" title="${esc(s.it.name)} ${fmt(s.it.start + 1)}..${fmt(s.it.end)}" style="left:${x}px;top:${laneBar[s.lane]}px;width:${w}px;height:${h}px;background:${s.it.color};${clip}${s.tipL && s.x1 === s.a ? 'padding-left:10px;' : ''}">${s.it.kind === 'o' ? '' : label}</div>`);
        keyLane[s.it.key] = s.lane;
      }
      for (const a of (aaRows[ri] || [])) {
        const l = keyLane[a.key]; if (l === undefined || laneAA[l] === undefined) continue;
        const x = G + (a.mid - rs) * CW - CW;
        parts.push(`<div class="aa${a.aa === '*' ? ' st' : (a.aa === 'M' ? ' m' : '')}" style="left:${x}px;top:${laneAA[l]}px;width:${CW * 3}px">${a.aa}</div>`);
        if (a.res % 5 === 0) parts.push(`<div class="aanum" style="left:${x - CW * 1.5}px;top:${laneNum[l]}px;width:${CW * 6}px">${a.res}</div>`);
      }
      const hgt = y + 5;
      rows.push({ rs, re, y0, y1: strandBottom });
      html.push(`<div class="srow" data-i="${ri}" style="height:${hgt}px">${parts.join('')}<div class="selov" hidden style="top:${y0 - 1}px;height:${strandBottom - y0 + 2}px"></div><div class="caret-el" hidden style="top:${y0 - 1}px;height:${strandBottom - y0 + 2}px"></div></div>`);
    }
    host.innerHTML = html.join('');
    S.seq = { rows, G, bpr, els: [...host.querySelectorAll('.srow')] };
    updateSeqSel();
  }

  // Wraps each base in a coloured span, but only while base colouring is switched on,
  // since a span per base is wasted work otherwise.
  function paintBases(str, o) {
    if (!o.baseColors || !o.baseColors.enabled) return str;
    let out = '', run = '', cur = '';
    for (const ch of str) {
      const c = o.baseColors[ch] || 'inherit';
      if (c === cur) { run += ch; continue; }
      if (run) out += `<span style="color:${cur}">${run}</span>`;
      cur = c; run = ch;
    }
    if (run) out += `<span style="color:${cur}">${run}</span>`;
    return out;
  }
  function updateSeqSel() {
    const d = D(), V = S.seq; if (!V || !d) return;
    const n = d.seq.length, L = selLen(d);
    const ranges = L ? DNA.segs(d.sel, n) : [];
    const p = d.sel.start;
    V.rows.forEach((r, i) => {
      const el = V.els[i]; if (!el) return;
      const ov = el.querySelector('.selov'), ca = el.querySelector('.caret-el');
      let a = Infinity, b = -Infinity;
      for (const [s, e] of ranges) { const x1 = Math.max(s, r.rs), x2 = Math.min(e, r.re); if (x2 > x1) { a = Math.min(a, x1); b = Math.max(b, x2); } }
      if (b > a) { ov.hidden = false; ov.style.left = (V.G + (a - r.rs) * CW) + 'px'; ov.style.width = ((b - a) * CW) + 'px'; } else ov.hidden = true;
      const hasCaret = !L && ((p >= r.rs && p < r.re) || (p === n && r.re === n));
      ca.hidden = !hasCaret;
      if (hasCaret) ca.style.left = (V.G + (p - r.rs) * CW - 1) + 'px';
    });
  }
  function scrollSeqTo(pos) {
    const V = S.seq, sc = $('#seqscroll'); if (!V || !sc) return;
    const i = Math.min(V.rows.length - 1, Math.floor(pos / V.bpr)); const el = V.els[i]; if (!el) return;
    const top = el.offsetTop, h = el.offsetHeight;
    if (top < sc.scrollTop || top + h > sc.scrollTop + sc.clientHeight) sc.scrollTop = Math.max(0, top - sc.clientHeight / 3);
  }

  function bindSeq(sc) {
    const d = D();
    const posFrom = (e, rowEl) => {
      const r = S.seq.rows[+rowEl.dataset.i], rect = rowEl.getBoundingClientRect();
      const col = Math.round((e.clientX - rect.left - S.seq.G) / CW);
      return r.rs + Math.max(0, Math.min(r.re - r.rs, col));
    };
    let drag = false, auto = null, lastEv = null;
    sc.addEventListener('mousedown', e => {
      if (e.button !== 0) return;
      const k = e.target.closest('[data-k]');
      if (k) { clickItem(k.dataset.k, e); sc.focus(); e.preventDefault(); return; }
      const row = e.target.closest('.srow'); if (!row) return;
      e.preventDefault(); sc.focus();
      const p = posFrom(e, row);
      if (e.shiftKey) { extendTo(p); return; }
      drag = true; setSel(p, p, { anchor: p });
    });
    sc.addEventListener('dblclick', e => {
      const k = e.target.closest('[data-k]'); if (!k) return;
      const it = findItem(D(), k.dataset.k);
      if (it && it.kind === 'f') dlgFeature(it.ref); else if (it && it.kind === 'p') dlgPrimer(it.ref);
    });
    sc.addEventListener('contextmenu', e => {
      e.preventDefault();
      const k = e.target.closest('[data-k]');
      if (k && !selLen(D())) clickItem(k.dataset.k, e);
      showCtx(e.clientX, e.clientY, selCtxItems(k ? itemCtx(D(), k.dataset.k) : []));
    });
    const move = e => {
      if (!drag) return; lastEv = e;
      const rect = sc.getBoundingClientRect();
      const y = Math.max(rect.top + 2, Math.min(rect.bottom - 2, e.clientY));
      const el = document.elementFromPoint(Math.max(rect.left + 2, Math.min(rect.right - 2, e.clientX)), y);
      const row = el && el.closest && el.closest('.srow');
      if (row) { const p = posFrom(e, row); setSel(Math.min(D().anchor, p), Math.max(D().anchor, p)); }
      clearInterval(auto); auto = null;
      if (e.clientY > rect.bottom - 10 || e.clientY < rect.top + 10) auto = setInterval(() => { sc.scrollTop += (lastEv.clientY > rect.bottom - 10 ? 18 : -18); move(lastEv); }, 40);
    };
    window.addEventListener('mousemove', move);
    window.addEventListener('mouseup', () => { drag = false; clearInterval(auto); auto = null; });
    sc.addEventListener('keydown', e => {
      const d = D(); if (!d || e.ctrlKey || e.metaKey || e.altKey) return;
      const n = d.seq.length, k = e.key;
      if (k.length === 1 && /[ACGTNRYSWKMBDHVacgtnryswkmbdhv]/.test(k)) { e.preventDefault(); replaceSel(k.toUpperCase(), 'Typed bases', 'type'); return; }
      if (k === 'Backspace') {
        e.preventDefault();
        if (selLen(d)) deleteSel();
        else if (d.sel.start > 0) { const p = d.sel.start - 1; editSeq(d, p, 1, ''); d.sel = { start: p, end: p }; d.anchor = p; commit(d, 'Deleted bases', { merge: 'bs' }); renderAll(true); }
        return;
      }
      if (k === 'Delete') {
        e.preventDefault();
        if (selLen(d)) deleteSel();
        else if (d.sel.start < n && n > 1) { const p = d.sel.start; editSeq(d, p, 1, ''); commit(d, 'Deleted bases', { merge: 'del' }); renderAll(true); }
        return;
      }
      const bpr = S.seq.bpr;
      const mv = { ArrowLeft: -1, ArrowRight: 1, ArrowUp: -bpr, ArrowDown: bpr, Home: 'home', End: 'end' }[k];
      if (mv === undefined) return;
      e.preventDefault();
      const L = selLen(d);
      let cur = e.shiftKey ? (d.sel.start === d.anchor ? d.sel.end : d.sel.start) : (L ? (mv < 0 ? d.sel.start : d.sel.end) : d.sel.start);
      if (mv === 'home') cur = e.shiftKey ? 0 : Math.floor(cur / bpr) * bpr;
      else if (mv === 'end') cur = e.shiftKey ? n : Math.min(n, Math.floor(cur / bpr) * bpr + bpr);
      else if (!L || e.shiftKey) cur += mv;
      cur = Math.max(0, Math.min(n, cur));
      if (e.shiftKey) setSel(Math.min(d.anchor, cur), Math.max(d.anchor, cur), { scroll: true });
      else setSel(cur, cur, { anchor: cur, scroll: true });
    });
    sc.addEventListener('paste', e => {
      const t = (e.clipboardData || window.clipboardData).getData('text');
      e.preventDefault(); pasteText(t);
    });
  }
  function pasteText(t) {
    const s = DNA.clean(t);
    if (!s) { toast('The clipboard has no DNA bases to paste.'); return; }
    replaceSel(s, `Pasted ${fmt(s.length)} bp`);
  }
  async function pasteFromClipboard() {
    try { const t = await navigator.clipboard.readText(); pasteText(t); }
    catch (err) { toast('Click in the sequence view and press Ctrl+V to paste.'); }
  }

  // ---------------------------------------------------------------- find
  function bindFind() {
    const q = $('#findq'), m = $('#findmode');
    q.addEventListener('input', () => { S.find.q = q.value; runFind(false); });
    m.addEventListener('change', () => { S.find.mode = m.value; runFind(false); });
    q.addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); stepFind(e.shiftKey ? -1 : 1); } if (e.key === 'Escape') { $('#seqscroll').focus(); } });
    $('#findnext').onclick = () => stepFind(1);
    $('#findprev').onclick = () => stepFind(-1);
    $('#gotobtn').onclick = dlgGoto;
  }
  function runFind(jump) {
    const d = D(), f = S.find;
    f.res = !f.q ? [] : (f.mode === 'dna' ? DNA.searchDNA(d.seq, d.circular, f.q) : DNA.searchProtein(d.seq, d.circular, f.q));
    f.i = -1;
    if (f.res.length) { const p = d.sel.start; f.i = f.res.findIndex(r => r.start >= p); if (f.i < 0) f.i = 0; selectRange(f.res[f.i]); }
    updateFindCount();
  }
  function stepFind(dir) {
    const f = S.find; if (!f.res.length) { runFind(true); return; }
    f.i = (f.i + dir + f.res.length) % f.res.length;
    selectRange(f.res[f.i]); updateFindCount();
  }
  function updateFindCount() {
    const c = $('#findcnt'); if (!c) return;
    const f = S.find;
    c.textContent = !f.q ? '' : (f.res.length ? `${f.i + 1} of ${f.res.length}` + (f.res[f.i] && f.res[f.i].strand < 0 ? ' (−)' : '') : 'No matches');
  }
  function openFind() { if (!D()) return; if (S.view !== 'seq') setView('seq'); setTimeout(() => { const q = $('#findq'); if (q) { q.focus(); q.select(); } }, 10); }

  // ================================================================ ENZYMES TABLE
  let enzFilter = '', enzSet = 'cutters', enzSort = 'name';
  function renderEnzymes() {
    const d = D(), A = analyze(d), n = d.seq.length;
    const v = $('#view');
    let list = DNA.ENZYMES.map(e => ({ e, s: A.sites[e.name] || [] }));
    if (enzSet === 'cutters') list = list.filter(x => x.s.length);
    else if (enzSet === 'unique') list = list.filter(x => x.s.length === 1);
    else if (enzSet === 'noncutters') list = list.filter(x => !x.s.length);
    if (enzFilter) { const q = enzFilter.toUpperCase(); list = list.filter(x => x.e.name.toUpperCase().includes(q) || x.e.site.includes(q)); }
    list.sort((a, b) => enzSort === 'count' ? a.s.length - b.s.length || a.e.name.localeCompare(b.e.name) : a.e.name.localeCompare(b.e.name));
    const siteHtml = e => { const s = e.site; const t = e.top; if (t >= 0 && t <= s.length) return esc(s.slice(0, t)) + '<i>^</i>' + esc(s.slice(t)); return esc(s) + ` (${e.top}/${e.bot})`; };
    const ends = e => e.overhang > 0 ? `5′ (${e.overhang})` : (e.overhang < 0 ? `3′ (${-e.overhang})` : 'blunt');
    v.innerHTML = `<div class="panel"><div class="ptools">
      <input type="text" id="enzq" placeholder="Filter by name or site" value="${esc(enzFilter)}">
      <select id="enzset"><option value="cutters" ${enzSet === 'cutters' ? 'selected' : ''}>Enzymes that cut</option><option value="unique" ${enzSet === 'unique' ? 'selected' : ''}>Unique cutters</option><option value="noncutters" ${enzSet === 'noncutters' ? 'selected' : ''}>Non-cutters</option><option value="all" ${enzSet === 'all' ? 'selected' : ''}>All enzymes</option></select>
      <button class="btn pri" id="enzpick">Choose Enzymes…</button><button class="btn" id="enzshow">Show chosen on map</button><button class="btn" id="enzclear">Clear chosen</button><button class="btn pri" id="enzdig" ${d.chosen.length ? '' : 'disabled'}>Digest chosen (${d.chosen.length})…</button>
      <span class="note">${list.length} enzymes</span></div>
      <div class="tscroll"><table class="grid"><thead><tr><th style="width:30px">✓</th><th class="sort" data-s="name">Enzyme</th><th>Recognition site</th><th class="sort" data-s="count"># Sites</th><th>Cut positions</th><th>Ends</th></tr></thead><tbody>
      ${list.map(x => `<tr data-e="${x.e.name}"><td><input type="checkbox" data-c="${x.e.name}" ${d.chosen.includes(x.e.name) ? 'checked' : ''}></td><td><b>${x.e.name}</b></td><td class="cut-site">${siteHtml(x.e)}</td><td>${x.s.length}</td><td class="wrap">${x.s.map(s => s.cut == null ? '—' : fmt(s.cut)).join(', ')}</td><td>${ends(x.e)}</td></tr>`).join('')}
      </tbody></table></div></div>`;
    const q = $('#enzq');
    q.oninput = () => { enzFilter = q.value; const p = q.selectionStart; renderEnzymes(); const q2 = $('#enzq'); q2.focus(); q2.setSelectionRange(p, p); };
    $('#enzset').onchange = e => { enzSet = e.target.value; renderEnzymes(); };
    v.querySelectorAll('th.sort').forEach(th => th.onclick = () => { enzSort = th.dataset.s; renderEnzymes(); });
    v.querySelector('tbody').onclick = e => {
      const c = e.target.closest('[data-c]');
      if (c) { const nm = c.dataset.c; d.chosen = c.checked ? [...new Set([...d.chosen, nm])] : d.chosen.filter(x => x !== nm); const b = $('#enzdig'); b.disabled = !d.chosen.length; b.textContent = `Digest chosen (${d.chosen.length})…`; return; }
      const tr = e.target.closest('tr[data-e]'); if (!tr) return;
      const s = (A.sites[tr.dataset.e] || [])[0]; if (!s) return;
      let end = s.pos + s.len; if (end > n) end -= n;
      selectRange({ start: s.pos, end });
      v.querySelectorAll('tr.on').forEach(r => r.classList.remove('on')); tr.classList.add('on');
    };
    v.querySelector('tbody').ondblclick = e => { const tr = e.target.closest('tr[data-e]'); if (tr && A.sites[tr.dataset.e]) { const nm = tr.dataset.e; if (!d.chosen.includes(nm)) d.chosen.push(nm); d.opts.enzMode = 'chosen'; S.view = 'map'; renderAll(); } };
    $('#enzpick').onclick = () => dlgEnzymeChooser();
    $('#enzshow').onclick = () => { d.opts.enzMode = 'chosen'; S.view = 'map'; renderAll(); };
    $('#enzclear').onclick = () => { d.chosen = []; renderEnzymes(); };
    $('#enzdig').onclick = () => dlgDigest();
  }

  // ================================================================ FEATURES TABLE
  function renderFeatures() {
    const d = D(), n = d.seq.length, v = $('#view');
    const fs = d.features.slice().sort((a, b) => a.start - b.start);
    v.innerHTML = `<div class="panel"><div class="ptools"><button class="btn pri" id="fadd">Add Feature…</button><button class="btn" id="fdet">Detect Common Features</button><span class="note">${fs.length} features · double-click a row to edit</span></div>
      <div class="tscroll">${fs.length ? `<table class="grid"><thead><tr><th></th><th>Name</th><th>Type</th><th>Location</th><th>Length</th><th>Direction</th><th>Notes</th><th></th></tr></thead><tbody>
      ${fs.map(f => `<tr data-f="${f.id}"><td><span class="sw" style="background:${f.color}"></span></td><td><b>${esc(f.name)}</b></td><td>${esc(f.type)}</td><td>${fmt(f.start + 1)} .. ${fmt(f.end)}${f.end <= f.start ? ' (spans origin)' : ''}</td><td>${fmt(DNA.rlen(f, n))} bp${f.type === 'CDS' ? ` · ${Math.floor(DNA.rlen(f, n) / 3) - 1} aa` : ''}</td><td>${f.strand === 1 ? '→ forward' : f.strand === -1 ? '← reverse' : '— none'}</td><td class="wrap">${esc(f.notes || '')}</td><td><button class="ib" data-ed="${f.id}">Edit</button><button class="ib dng" data-del="${f.id}">Delete</button></td></tr>`).join('')}
      </tbody></table>` : `<div class="empty"><h2>No features yet</h2><p>Select a region in the Map or Sequence view and choose Add Feature, or let GeneCanvas detect common features such as promoters, tags and recombination sites.</p></div>`}</div></div>`;
    $('#fadd').onclick = () => dlgFeature(null);
    $('#fdet').onclick = detectCommon;
    const tb = v.querySelector('tbody'); if (!tb) return;
    tb.onclick = e => {
      const ed = e.target.closest('[data-ed]'); if (ed) { dlgFeature(d.features.find(f => f.id === +ed.dataset.ed)); return; }
      const del = e.target.closest('[data-del]'); if (del) { removeFeature(+del.dataset.del); return; }
      const tr = e.target.closest('tr[data-f]'); if (!tr) return;
      const f = d.features.find(x => x.id === +tr.dataset.f); selectRange(f);
      tb.querySelectorAll('tr.on').forEach(r => r.classList.remove('on')); tr.classList.add('on');
    };
    tb.ondblclick = e => { const tr = e.target.closest('tr[data-f]'); if (tr) dlgFeature(d.features.find(x => x.id === +tr.dataset.f)); };
  }

  // ================================================================ PRIMERS TABLE
  function renderPrimers() {
    const d = D(), A = analyze(d), v = $('#view');
    v.innerHTML = `<div class="panel"><div class="ptools"><button class="btn pri" id="padd">Add Primer…</button><button class="btn" id="ptop" ${selLen(d) ? '' : 'disabled'}>From Selection (Top)</button><button class="btn" id="pbot" ${selLen(d) ? '' : 'disabled'}>From Selection (Bottom)</button><button class="btn" id="pstruct" ${d.primers.length ? '' : 'disabled'}>Structures…</button><span class="note">Tm: nearest-neighbour, 50 mM Na⁺, 250 nM primer</span></div>
      <div class="tscroll">${d.primers.length ? `<table class="grid"><thead><tr><th>Name</th><th>Sequence (5′→3′)</th><th>Length</th><th>Tm</th><th>GC</th><th>Binding sites</th><th></th></tr></thead><tbody>
      ${d.primers.map(p => { const b = A.psites[p.id] || []; return `<tr data-p="${p.id}"><td><b>${esc(p.name)}</b></td><td class="mono wrap">${esc(p.seq)}</td><td>${p.seq.length}</td><td>${DNA.tm(p.seq).toFixed(1)} °C</td><td>${DNA.gcPct(p.seq).toFixed(0)}%</td><td>${b.length ? b.map(x => `${fmt(x.start + 1)}..${fmt(x.end)} ${x.strand > 0 ? '→' : '←'}${x.full ? '' : ' (3′ 15 nt)'}`).join('<br>') : '<span style="color:var(--danger)">does not bind</span>'}</td><td><button class="ib" data-struct="${p.id}">Structure</button><button class="ib" data-ed="${p.id}">Edit</button><button class="ib dng" data-del="${p.id}">Delete</button></td></tr>`; }).join('')}
      </tbody></table>` : `<div class="empty"><h2>No primers yet</h2><p>Select bases and choose Primer from Selection, or add a primer by sequence. Primers with 5′ tails are matched by their 3′ end.</p></div>`}</div></div>`;
    $('#padd').onclick = () => dlgPrimer(null);
    $('#ptop').onclick = () => primerFromSel(1);
    $('#pbot').onclick = () => primerFromSel(-1);
    $('#pstruct').onclick = () => dlgPrimerStruct(null, null);
    const tb = v.querySelector('tbody'); if (!tb) return;
    tb.onclick = e => {
      const st = e.target.closest('[data-struct]'); if (st) { dlgPrimerStruct(d.primers.find(p => p.id === +st.dataset.struct), null); return; }
      const ed = e.target.closest('[data-ed]'); if (ed) { dlgPrimer(d.primers.find(p => p.id === +ed.dataset.ed)); return; }
      const del = e.target.closest('[data-del]'); if (del) { removePrimer(+del.dataset.del); return; }
      const tr = e.target.closest('tr[data-p]'); if (!tr) return;
      const b = (A.psites[+tr.dataset.p] || [])[0]; if (b) selectRange(b);
      tb.querySelectorAll('tr.on').forEach(r => r.classList.remove('on')); tr.classList.add('on');
    };
  }

  // ================================================================ HISTORY
  // ---------------------------------------------------------------- construction history
  // A small picture of a sequence, used in the tree so each piece is recognisable.
  function miniMap(src, W, H) {
    const n = src.seq ? src.seq.length : (src.length || 0);
    if (!n) return '';
    const feats = (src.features || []).filter(f => f.end !== f.start);
    if (src.circular) {
      const cx = W / 2, cy = H / 2, R = Math.min(W, H) / 2 - 10;
      const P = (r, a) => `${(cx + r * Math.cos(a - Math.PI / 2)).toFixed(1)},${(cy + r * Math.sin(a - Math.PI / 2)).toFixed(1)}`;
      const arcs = feats.map(f => {
        const a1 = f.start / n * TAU, span = Math.max(0.03, DNA.rlen(f, n) / n * TAU), a2 = a1 + span;
        const lg = span > Math.PI ? 1 : 0;
        return `<path d="M${P(R + 4, a1)} A${R + 4},${R + 4} 0 ${lg} 1 ${P(R + 4, a2)} L${P(R - 4, a2)} A${R - 4},${R - 4} 0 ${lg} 0 ${P(R - 4, a1)} Z" fill="${f.color || '#c9ced6'}" stroke="rgba(0,0,0,.25)" stroke-width=".6"><title>${esc(f.name)}</title></path>`;
      }).join('');
      return `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}"><circle cx="${cx}" cy="${cy}" r="${R}" fill="none" stroke="#2f3644" stroke-width="1.6"/>${arcs}</svg>`;
    }
    const m = 8, w = W - 2 * m, y = H / 2;
    const boxes = feats.map(f => {
      const x1 = m + f.start / n * w, x2 = m + Math.min(n, f.end) / n * w;
      return `<rect x="${x1.toFixed(1)}" y="${y - 5}" width="${Math.max(1.5, x2 - x1).toFixed(1)}" height="10" fill="${f.color || '#c9ced6'}" stroke="rgba(0,0,0,.25)" stroke-width=".6"><title>${esc(f.name)}</title></rect>`;
    }).join('');
    return `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}"><line x1="${m}" y1="${y}" x2="${m + w}" y2="${y}" stroke="#2f3644" stroke-width="1.6"/>${boxes}</svg>`;
  }

  function provCard(node, depth) {
    const sizeTxt = `${fmt((node.seq || '').length)} bp`;
    const kids = node.source ? [node.source] : [];
    return `<div class="pvnode">
      <div class="pvcard" data-pv="${esc(node.uid)}">
        <div class="pvmap">${miniMap(node, 150, node.circular ? 86 : 44)}</div>
        <div class="pvname">${esc(node.name)}</div>
        <div class="pvsize">${sizeTxt}${node.circular ? ' · circular' : ''}</div>
        ${node.features && node.features.length ? `<div class="pvfeat">${node.features.slice(0, 4).map(f => esc(f.name)).join(', ')}${node.features.length > 4 ? ` +${node.features.length - 4}` : ''}</div>` : ''}
        <button class="btn pvopen" data-open="${esc(node.uid)}">Open</button>
      </div>
      ${kids.map(k => `<div class="pvstep"><span class="pvmethod">${esc(k.method)}</span>${k.detail ? `<span class="pvdetail">${esc(k.detail)}</span>` : ''}${k.primers && k.primers.length ? `<span class="pvprimers">${k.primers.map(esc).join(' · ')}</span>` : ''}</div>
        <div class="pvrow">${(k.inputs || []).map(c => provCard(c, depth + 1)).join('')}</div>`).join('')}
    </div>`;
  }

  // Walks the stored tree and gives every node an id so its Open button can find it.
  function tagProv(node, seen) {
    node.uid = 'p' + (seen.n++);
    seen.map[node.uid] = node;
    if (node.source) for (const c of (node.source.inputs || [])) tagProv(c, seen);
    return node;
  }

  function renderHistory() {
    const d = D();
    if (S.histTab === undefined) S.histTab = d.provenance ? 'tree' : 'edits';
    const v = $('#view');
    const tabBar = `<div class="ptools">
      <div class="seg" style="width:260px"><button data-ht="tree" class="${S.histTab === 'tree' ? 'on' : ''}">How it was made</button><button data-ht="edits" class="${S.histTab === 'edits' ? 'on' : ''}">Edits (${d.hist.length})</button></div>
      <span class="note">${d.provenance ? esc(d.provenance.method) : 'No construction record — this sequence was opened rather than assembled here.'}</span></div>`;
    if (S.histTab === 'tree') {
      if (!d.provenance) {
        v.innerHTML = `<div class="panel">${tabBar}<div class="empty"><h2>No construction record</h2><p>When a sequence is built here by restriction cloning, Gibson or Golden Gate assembly, Gateway recombination or codon optimisation, this view shows what it was made from: the method, the pieces, the primers, and the sequence of each piece.</p></div></div>`;
      } else {
        const seen = { n: 0, map: {} };
        const root = tagProv(JSON.parse(JSON.stringify({ name: d.name, seq: d.seq, circular: d.circular, features: d.features, source: d.provenance })), seen);
        v.innerHTML = `<div class="panel">${tabBar}<div class="tscroll"><div class="pvtree">${provCard(root, 0)}</div></div></div>`;
        v.querySelector('.pvtree').onclick = e => {
          const b = e.target.closest('[data-open]'); if (!b) return;
          const node = seen.map[b.dataset.open];
          if (!node || !node.seq) return;
          addDoc({ name: node.name, seq: node.seq, circular: !!node.circular, features: node.features || [], primers: [], histLabel: 'Opened from the construction history' });
        };
      }
    } else {
      const t = dt => new Date(dt).toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' });
      v.innerHTML = `<div class="panel">${tabBar}<div class="tscroll">
        ${d.hist.map((h, i) => ({ h, i })).reverse().map(({ h, i }) => `<div class="hist-item${i === d.hi ? ' on' : ''}${i > d.hi ? ' future' : ''}" data-h="${i}"><span class="t">${t(h.t)}</span><span>${esc(h.label)}</span><span class="len">${fmt(h.len)} bp</span></div>`).join('')}</div></div>`;
      v.querySelector('.tscroll').onclick = e => { const it = e.target.closest('[data-h]'); if (!it) return; restore(d, +it.dataset.h); renderAll(); };
    }
    v.querySelector('.seg').onclick = e => { const b = e.target.closest('[data-ht]'); if (!b) return; S.histTab = b.dataset.ht; renderView(); };
  }

  function oldRenderHistory() {
    const d = D(), v = $('#view');
    const t = dt => dt.toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' });
    v.innerHTML = `<div class="panel"><div class="ptools"><span class="note">Click any step to return the sequence to that state. Later steps stay available until you make a new change.</span></div><div class="tscroll">
      ${d.hist.map((h, i) => ({ h, i })).reverse().map(({ h, i }) => `<div class="hist-item${i === d.hi ? ' on' : ''}${i > d.hi ? ' future' : ''}" data-h="${i}"><span class="t">${t(h.t)}</span><span>${esc(h.label)}</span><span class="len">${fmt(h.len)} bp</span></div>`).join('')}</div></div>`;
    v.querySelector('.tscroll').onclick = e => { const it = e.target.closest('[data-h]'); if (!it) return; restore(d, +it.dataset.h); renderAll(); };
  }

  // ================================================================ feature / primer ops
  function addFeature(f) {
    const d = D();
    const nf = { id: uid(), name: f.name, type: f.type, start: f.start, end: f.end, strand: f.strand, color: f.color || (f.type === 'CDS' ? CDS_PALETTE[d.features.filter(x => x.type === 'CDS').length % CDS_PALETTE.length] : (TYPE_COL[f.type] || '#c9ced6')), notes: f.notes || '' };
    d.features.push(nf); commit(d, `Added feature ${nf.name}`); renderAll(true);
    return nf;
  }
  function removeFeature(id) {
    const d = D(); const f = d.features.find(x => x.id === id); if (!f) return;
    d.features = d.features.filter(x => x.id !== id); commit(d, `Deleted feature ${f.name}`); renderAll(true);
  }
  function removePrimer(id) {
    const d = D(); const p = d.primers.find(x => x.id === id); if (!p) return;
    d.primers = d.primers.filter(x => x.id !== id); commit(d, `Deleted primer ${p.name}`); renderAll(true);
  }
  function primerFromSel(strand) {
    const d = D(); const L = selLen(d);
    if (!L) { toast('Select the primer binding region first.'); return; }
    if (L > 120) { toast('Primers are limited to 120 nt. Select a shorter region.'); return; }
    const s = selSeq(d);
    dlgPrimer(null, { name: `Primer ${d.primers.length + 1}`, seq: strand === 1 ? s : DNA.rc(s) });
  }
  function selToCDS() {
    const d = D(); if (!selLen(d)) return;
    dlgFeature(null, { type: 'CDS' });
  }
  function detectCommon() {
    const d = D();
    const found = DNA.detectCommon(d.seq, d.circular, d.features);
    if (!found.length) { toast('No additional common features found.'); return; }
    for (const f of found) d.features.push({ id: uid(), name: f.name, type: f.type, start: f.start, end: f.end, strand: f.strand, color: TYPE_COL[f.type] || '#c9ced6', notes: 'Detected automatically' });
    commit(d, `Detected ${found.length} common feature${found.length > 1 ? 's' : ''}`); renderAll(true);
    toast(`Added ${found.length}: ${[...new Set(found.map(f => f.name))].join(', ')}`);
  }

  // ================================================================ dialogs
  let dlgKeyHandler = null;
  function dialog({ title, body, buttons, wide, onOpen }) {
    closeMenu(); closeCtx();
    const m = $('#modal');
    // Start from a fresh node so listeners added by a previous dialog cannot fire on this one.
    const prev = $('#dlg'), box = prev.cloneNode(false);
    prev.parentNode.replaceChild(box, prev);
    box.className = 'dlg' + (wide ? ' wide' : '');
    box.innerHTML = `<div class="hd">${esc(title)}</div><div class="bd">${body}</div><div class="ft">${buttons.map((b, i) => `<button class="btn${b.primary ? ' pri' : ''}${b.danger ? ' dng' : ''}" data-b="${i}"${b.id ? ` id="${b.id}"` : ''}${b.disabled ? ' disabled' : ''}>${esc(b.label)}</button>`).join('')}</div>`;
    m.hidden = false;
    box.querySelector('.ft').onclick = e => {
      const b = e.target.closest('[data-b]'); if (!b) return;
      const r = buttons[+b.dataset.b].act ? buttons[+b.dataset.b].act() : true;
      // An async handler keeps the dialog open and closes it itself when it is done.
      if (r && typeof r.then === 'function') return;
      if (r !== false) closeDialog();
    };
    dlgKeyHandler = e => {
      if (e.key === 'Escape') { e.preventDefault(); closeDialog(); }
      else if (e.key === 'Enter' && e.target.tagName !== 'TEXTAREA' && !e.target.closest('.ft')) {
        const i = buttons.findIndex(b => b.primary); if (i >= 0) { e.preventDefault(); box.querySelector(`[data-b="${i}"]`).click(); }
      }
    };
    box.addEventListener('keydown', dlgKeyHandler);
    if (onOpen) onOpen(box);
    // Focus something inside the dialog, otherwise Escape and Enter never reach it.
    box.tabIndex = -1;
    const first = box.querySelector('.bd input, .bd textarea, .bd select');
    setTimeout(() => {
      if (first) { first.focus(); if (first.select) first.select(); }
      else box.focus();
    }, 0);
  }
  function closeDialog() { const box = $('#dlg'); if (dlgKeyHandler) box.removeEventListener('keydown', dlgKeyHandler); $('#modal').hidden = true; box.innerHTML = ''; }
  const val = id => document.getElementById(id).value;
  const showErr = (msg) => { const e = $('#dlg .err'); if (e) e.textContent = msg; return false; };

  function dlgNew() {
    dialog({
      title: 'New DNA File',
      body: `<div class="form"><label for="nn">Name</label><input type="text" id="nn" value="New DNA">
        <label for="nt">Topology</label><select id="nt"><option value="0">Linear</option><option value="1">Circular</option></select>
        <div class="full"><textarea id="ns" rows="9" placeholder="Paste DNA bases, FASTA or a full GenBank record"></textarea><div class="hint">Non-DNA characters such as numbers and spaces are removed. A pasted GenBank record keeps its features and topology.</div><div class="err"></div></div></div>`,
      buttons: [{ label: 'Cancel' }, {
        label: 'Create', primary: true, act: () => {
          const raw = val('ns');
          let src;
          try {
            if (/^\s*(LOCUS|>)/.test(raw)) { src = DNA.parseAny(raw); if (val('nn') && val('nn') !== 'New DNA') src.name = val('nn'); }
            else { const s = DNA.clean(raw); if (!s) return showErr('Enter at least one base.'); src = { name: val('nn') || 'New DNA', seq: s, circular: val('nt') === '1', features: [], primers: [] }; }
          } catch (err) { return showErr(err.message); }
          addDoc(src);
        }
      }]
    });
  }
  function dlgInsert() {
    const d = D();
    dialog({
      title: 'Insert Bases',
      body: `<div class="form"><label>Position</label><div>${selLen(d) ? `Replace selection ${fmt(d.sel.start + 1)}–${fmt(d.sel.end)}` : `After base ${fmt(d.sel.start)}`}</div>
        <div class="full"><textarea id="ins" rows="7" placeholder="Bases to insert"></textarea><div class="err"></div></div></div>`,
      buttons: [{ label: 'Cancel' }, { label: 'Insert', primary: true, act: () => { const s = DNA.clean(val('ins')); if (!s) return showErr('Enter at least one base.'); replaceSel(s, `Inserted ${fmt(s.length)} bp`); } }]
    });
  }
  function dlgGoto() {
    const d = D(), n = d.seq.length;
    dialog({
      title: 'Go To / Select Range',
      body: `<div class="form"><label for="g1">From</label><input type="number" id="g1" min="1" max="${n}" value="${d.sel.start + 1}"><label for="g2">To (optional)</label><input type="number" id="g2" min="1" max="${n}" value="${selLen(d) ? (d.sel.end || n) : ''}"><div class="full hint">Leave “To” empty to place the cursor before the “From” base. Positions are 1-based.</div><div class="full err"></div></div>`,
      buttons: [{ label: 'Cancel' }, {
        label: 'Go', primary: true, act: () => {
          const a = parseInt(val('g1'), 10), b = parseInt(val('g2'), 10);
          if (!(a >= 1 && a <= n)) return showErr(`From must be between 1 and ${fmt(n)}.`);
          if (val('g2') && !(b >= 1 && b <= n)) return showErr(`To must be between 1 and ${fmt(n)}.`);
          if (!val('g2')) setSel(a - 1, a - 1, { anchor: a - 1, scroll: true });
          else if (b >= a) setSel(a - 1, b, { anchor: a - 1, scroll: true });
          else if (d.circular) setSel(a - 1, b, { anchor: a - 1, scroll: true });
          else return showErr('“To” must not be smaller than “From” on a linear sequence.');
        }
      }]
    });
  }
  function dlgRename() {
    const d = D();
    dialog({
      title: 'Rename Sequence', body: `<div class="form"><label for="rn">Name</label><input type="text" id="rn" value="${esc(d.name)}"><div class="full err"></div></div>`,
      buttons: [{ label: 'Cancel' }, { label: 'Rename', primary: true, act: () => { const v = val('rn').trim(); if (!v) return showErr('Enter a name.'); d.name = v; commit(d, 'Renamed to ' + v); renderAll(true); } }]
    });
  }

  function dlgFeature(f, preset = {}) {
    const d = D(), n = d.seq.length, isNew = !f;
    if (isNew && !selLen(d)) { toast('Select a region first, then add a feature.'); return; }
    const cur = f || { name: preset.name || 'New feature', type: preset.type || 'misc_feature', start: d.sel.start, end: d.sel.end, strand: preset.strand !== undefined ? preset.strand : 1, color: '', notes: '' };
    const color = cur.color || TYPE_COL[cur.type] || '#c9ced6';
    dialog({
      title: isNew ? 'Add Feature' : `Edit Feature “${cur.name}”`,
      body: `<div class="form">
        <label for="fn">Name</label><input type="text" id="fn" value="${esc(cur.name)}">
        <label for="ft">Type</label><select id="ft">${FEATURE_TYPES.map(t => `<option ${t === cur.type ? 'selected' : ''}>${t}</option>`).join('')}</select>
        <label for="fs">Start</label><input type="number" id="fs" min="1" max="${n}" value="${cur.start + 1}">
        <label for="fe">End</label><input type="number" id="fe" min="1" max="${n}" value="${cur.end === 0 ? n : cur.end}">
        <label for="fd">Direction</label><select id="fd"><option value="1" ${cur.strand === 1 ? 'selected' : ''}>→ Forward</option><option value="-1" ${cur.strand === -1 ? 'selected' : ''}>← Reverse</option><option value="0" ${!cur.strand ? 'selected' : ''}>— None</option></select>
        <label for="fc">Color</label><input type="color" id="fc" value="${color}">
        <label for="fno">Notes</label><input type="text" id="fno" value="${esc(cur.notes || '')}">
        <div class="full hint" id="finfo"></div><div class="full err"></div></div>`,
      buttons: [...(isNew ? [] : [{ label: 'Delete', danger: true, act: () => removeFeature(f.id) }]), { label: 'Cancel' }, {
        label: isNew ? 'Add' : 'Save', primary: true, act: () => {
          const name = val('fn').trim(), a = parseInt(val('fs'), 10), b = parseInt(val('fe'), 10);
          if (!name) return showErr('Enter a feature name.');
          if (!(a >= 1 && a <= n && b >= 1 && b <= n)) return showErr(`Start and end must be between 1 and ${fmt(n)}.`);
          if (b < a && !d.circular) return showErr('End must not be before start on a linear sequence.');
          const nf = { name, type: val('ft'), start: a - 1, end: b, strand: +val('fd'), color: val('fc'), notes: val('fno').trim() };
          if (b < a) nf.end = b; // spans origin on a circular sequence
          if (isNew) addFeature(nf);
          else { Object.assign(f, nf); commit(d, `Edited feature ${name}`); renderAll(true); }
        }
      }],
      onOpen: box => {
        const upd = () => {
          const a = parseInt(val('fs'), 10), b = parseInt(val('fe'), 10);
          if (!(a >= 1 && b >= 1)) return;
          const r = { start: a - 1, end: b }; const L = DNA.rlen(r, n);
          let s = `${fmt(L)} bp`;
          if (val('ft') === 'CDS' && L >= 3) { const seq = DNA.sub(d.seq, r.start, r.end); const p = DNA.translate(+val('fd') === -1 ? DNA.rc(seq) : seq); s += ` · ${p.length} codons · ${p.startsWith('M') ? 'starts with Met' : 'no ATG start'} · ${p.indexOf('*') === p.length - 1 ? 'ends with a stop codon' : (p.includes('*') ? 'internal stop codon' : 'no stop codon')}`; }
          $('#finfo').textContent = s;
        };
        box.querySelectorAll('input,select').forEach(el => el.addEventListener('input', upd));
        $('#ft').addEventListener('change', () => { if (isNew) $('#fc').value = TYPE_COL[val('ft')] || '#c9ced6'; upd(); });
        upd();
      }
    });
  }

  function dlgPrimer(p, preset) {
    const d = D(), isNew = !p;
    const cur = p || preset || { name: `Primer ${d.primers.length + 1}`, seq: '' };
    dialog({
      title: isNew ? 'Add Primer' : `Edit Primer “${cur.name}”`,
      body: `<div class="form"><label for="pn">Name</label><input type="text" id="pn" value="${esc(cur.name)}">
        <label for="ps">Sequence 5′→3′</label><textarea id="ps" rows="3" style="min-height:60px">${esc(cur.seq)}</textarea>
        <label>Properties</label><div id="pinfo" class="hint" style="margin:0"></div>
        <div class="full err"></div></div>`,
      buttons: [...(isNew ? [] : [{ label: 'Delete', danger: true, act: () => removePrimer(p.id) }]), { label: 'Cancel' }, {
        label: isNew ? 'Add' : 'Save', primary: true, act: () => {
          const name = val('pn').trim(), s = DNA.clean(val('ps'));
          if (!name) return showErr('Enter a primer name.');
          if (s.length < 6) return showErr('A primer needs at least 6 bases.');
          if (isNew) { d.primers.push({ id: uid(), name, seq: s, notes: '' }); commit(d, `Added primer ${name}`); }
          else { p.name = name; p.seq = s; commit(d, `Edited primer ${name}`); }
          renderAll(true);
        }
      }],
      onOpen: () => {
        const upd = () => {
          const s = DNA.clean(val('ps'));
          const b = s.length >= 6 ? primerBind(d, { seq: s }) : [];
          $('#pinfo').innerHTML = s.length ? `${s.length} nt · Tm ${DNA.tm(s).toFixed(1)} °C · GC ${DNA.gcPct(s).toFixed(0)}% · ${b.length ? `binds at ${b.map(x => `${fmt(x.start + 1)}..${fmt(x.end)} ${x.strand > 0 ? '→' : '←'}`).join(', ')}` : 'no binding site in this sequence'}` : 'Enter a sequence.';
        };
        $('#ps').addEventListener('input', upd); upd();
      }
    });
  }

  function dlgTranslate() {
    const d = D(); if (!selLen(d)) { toast('Select bases to translate.'); return; }
    const s = selSeq(d);
    const frames = [1, 2, 3, -1, -2, -3].map(f => { const src = f > 0 ? s.slice(f - 1) : DNA.rc(s).slice(-f - 1); return { f, p: DNA.translate(src) }; });
    const p1 = frames[0].p;
    const clean = p1.replace(/\*$/, '');
    dialog({
      title: `Translation of ${fmt(selLen(d))} bp`, wide: true,
      body: `<table class="grid"><thead><tr><th>Frame</th><th>Protein</th><th>Length</th></tr></thead><tbody>${frames.map(x => `<tr><td>${x.f > 0 ? '+' : ''}${x.f}</td><td class="mono wrap">${esc(x.p).replace(/\*/g, '<b style="color:var(--danger)">*</b>')}</td><td>${x.p.length} aa</td></tr>`).join('')}</tbody></table>
        <div class="hint">Frame +1 protein: ${clean.length} aa · ${(DNA.proteinMW(clean.replace(/\*/g, '')) / 1000).toFixed(2)} kDa</div>`,
      buttons: [{ label: 'Copy Frame +1', act: () => { copyText(p1, 'Copied protein sequence'); return false; } }, { label: 'Close', primary: true }]
    });
  }

  // ---------- Digest & gel
  const LADDER = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250];
  function gelSVG(lanes) {
    const W = 120 + lanes.length * 88, H = 480, top = 56, bot = 452;
    const lmax = Math.log10(12000), lmin = Math.log10(150);
    const y = s => top + 8 + (lmax - Math.log10(Math.min(Math.max(s, 150), 12000))) / (lmax - lmin) * (bot - top - 16);
    let g = `<rect width="${W}" height="${H}" fill="#121419" rx="4"/><defs><filter id="glow" x="-20%" y="-200%" width="140%" height="500%"><feGaussianBlur stdDeviation="1.3"/></filter></defs>`;
    lanes.forEach((ln, i) => {
      const x = 70 + i * 88;
      g += `<rect x="${x}" y="${top - 12}" width="64" height="7" fill="#2a2e37" rx="1"/>`;
      g += `<text x="${x + 32}" y="22" fill="#d4d9e3" font-size="11" text-anchor="middle" font-weight="600">${esc(ln.label)}</text>`;
      if (ln.sub) g += `<text x="${x + 32}" y="36" fill="#8e96a5" font-size="10" text-anchor="middle">${esc(ln.sub)}</text>`;
      const maxAmt = Math.max(...ln.bands.map(b => b.amt), 1);
      for (const b of ln.bands) {
        const op = 0.35 + 0.65 * Math.min(1, b.amt / maxAmt);
        const yy = y(b.size);
        g += `<rect x="${x + 2}" y="${(yy - 2.2).toFixed(1)}" width="60" height="4.4" fill="#f5f8ff" opacity="${op.toFixed(2)}" filter="url(#glow)" rx="1"><title>${fmt(b.size)} bp${b.note ? ' · ' + b.note : ''}</title></rect>`;
        if (ln.ladder) g += `<text x="${x - 6}" y="${(yy + 3.5).toFixed(1)}" fill="#a8b0be" font-size="10" text-anchor="end">${b.size >= 1000 ? (b.size / 1000) + ' kb' : b.size + ' bp'}</text>`;
      }
    });
    return `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="Segoe UI, Arial, sans-serif">${g}</svg>`;
  }
  function dlgDigest() {
    const d = D(), A = analyze(d), n = d.seq.length;
    const cutters = DNA.ENZYMES.filter(e => A.sites[e.name] && A.sites[e.name].some(s => s.cut != null)).map(e => e.name);
    let chosen = d.chosen.filter(x => cutters.includes(x));
    if (!chosen.length) {
      const common = ['EcoRI', 'HindIII', 'BamHI', 'XhoI', 'NotI', 'PstI', 'NdeI', 'XbaI', 'SalI', 'KpnI'];
      chosen = common.filter(nm => cutters.includes(nm)).slice(0, 2);
      if (!chosen.length) chosen = cutters.filter(nm => A.sites[nm].length === 1 && DNA.ENZ_BY_NAME[nm].core >= 6).slice(0, 2);
    }
    const render = () => {
      const fr = DNA.digest(d.seq, d.circular, A.sites, chosen);
      const lanes = [
        { label: '1 kb ladder', ladder: true, bands: LADDER.map(s => ({ size: s, amt: s === 3000 ? 2.2 : 1 })) },
        { label: 'Uncut', sub: d.circular ? 'supercoiled' : '', bands: [{ size: d.circular ? Math.round(n * 0.7) : n, amt: 1, note: d.circular ? 'supercoiled plasmid runs faster than its size' : '' }] },
        { label: chosen.length ? chosen.join(' + ') : 'No enzyme', sub: `${fr.length} fragment${fr.length > 1 ? 's' : ''}`, bands: fr.map(f => ({ size: f.uncut && d.circular ? Math.round(n * 0.7) : f.size, amt: f.size, note: `${f.left} → ${f.right}` })) }
      ];
      $('#gel').innerHTML = gelSVG(lanes);
      $('#frags').innerHTML = `<table class="grid"><thead><tr><th>#</th><th>Size</th><th>Coordinates</th><th>Ends</th></tr></thead><tbody>${fr.map((f, i) => `<tr><td>${i + 1}</td><td><b>${fmt(f.size)} bp</b></td><td>${f.uncut ? '—' : `${fmt(f.start + 1)} .. ${fmt(f.end)}`}</td><td>${f.uncut ? 'uncut circle' : `${esc(f.left)} / ${esc(f.right)}`}</td></tr>`).join('')}</tbody></table>`;
    };
    dialog({
      title: 'Restriction Digest & Simulated Agarose Gel', wide: true,
      body: `<div class="enz-pick" id="epick">${cutters.map(nm => `<label><input type="checkbox" value="${nm}" ${chosen.includes(nm) ? 'checked' : ''}> ${nm} <span style="color:var(--faint)">(${A.sites[nm].length})</span></label>`).join('')}</div>
        <div class="gel-wrap"><div id="gel"></div><div id="frags" style="flex:1;min-width:260px"></div></div>
        <div class="hint">Simulated 1% agarose gel, 1 kb ladder. Band brightness is proportional to fragment mass.</div>`,
      buttons: [{ label: 'Use as Chosen Enzymes', act: () => { d.chosen = chosen.slice(); d.opts.enzMode = 'chosen'; renderAll(true); } }, { label: 'Close', primary: true }],
      onOpen: () => {
        $('#epick').addEventListener('change', () => { chosen = [...document.querySelectorAll('#epick input:checked')].map(i => i.value); render(); });
        render();
      }
    });
  }

  function dlgShortcuts() {
    const rows = [['Ctrl+N', 'New DNA file'], ['Ctrl+O', 'Open file'], ['Ctrl+S / Ctrl+Shift+S', 'Save / Save As'], ['Ctrl+W', 'Close file'], ['Ctrl+Z / Ctrl+Y', 'Undo / Redo'], ['Ctrl+C / Ctrl+V', 'Copy / Paste bases'], ['Ctrl+A', 'Select all'], ['Ctrl+F', 'Find DNA or protein'], ['Ctrl+G', 'Go to position'], ['Ctrl+K', 'Add feature from selection'], ['Ctrl+1 … Ctrl+7', 'Map, Sequence, Enzymes, Features, Primers, Alignment, History'], ['Ctrl+D', 'Open from NCBI'], ['Ctrl+P', 'Print'], ['Ctrl + / Ctrl − / Ctrl+0', 'Zoom the map'], ['A C G T N', 'Type bases in the Sequence view'], ['Shift+click, Shift+arrows', 'Extend the selection'], ['Double-click a feature', 'Edit the feature'], ['Right-click', 'Context menu']];
    dialog({ title: 'Keyboard Shortcuts', body: `<table class="grid"><tbody>${rows.map(r => `<tr><td><b>${r[0]}</b></td><td>${r[1]}</td></tr>`).join('')}</tbody></table>`, buttons: [{ label: 'Close', primary: true }] });
  }
  function dlgAbout() {
    dialog({
      title: 'About GeneCanvas',
      body: `<p><b>GeneCanvas</b> ${host && host.version ? 'version ' + esc(host.version) : ''}<br>A desktop editor for plasmids and DNA sequences: maps, sequence editing, restriction analysis, features, primers and simulated gels.</p><p class="hint">Opens GenBank (.gb, .gbk), FASTA and .dna files. Saves GenBank and FASTA. The example plasmid is a synthetic sequence made for demonstration.</p>`,
      buttons: [{ label: 'Close', primary: true }]
    });
  }

  // ================================================================ NCBI download
  // Electron wraps a main-process error; show only the part written for the reader.
  function hostErr(e) { return String((e && e.message) || e).replace(/^Error invoking remote method '[^']*':\s*(Error:\s*)?/, ''); }
  async function ncbiSearch(db, term) {
    if (host && host.ncbiSearch) return host.ncbiSearch({ db, term, retmax: 25 }).catch(e => { throw new Error(hostErr(e)); });
    const base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    const r = await fetch(`${base}esearch.fcgi?db=${db}&term=${encodeURIComponent(term)}&retmax=25&retmode=json&tool=GeneCanvas`);
    const ids = ((await r.json()).esearchresult || {}).idlist || [];
    if (!ids.length) return { ids: [], items: [] };
    const s = await (await fetch(`${base}esummary.fcgi?db=${db}&id=${ids.join(',')}&retmode=json&tool=GeneCanvas`)).json();
    const res = s.result || {};
    return { ids, items: (res.uids || []).map(u => ({ uid: u, accession: res[u].caption || u, title: res[u].title || '', slen: res[u].slen || 0, topology: res[u].topology || '', organism: res[u].organism || '' })) };
  }
  async function ncbiFetch(db, id) {
    if (host && host.ncbiFetch) return host.ncbiFetch({ db, id }).catch(e => { throw new Error(hostErr(e)); });
    const r = await fetch(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=${db}&id=${encodeURIComponent(id)}&rettype=gbwithparts&retmode=text&tool=GeneCanvas`);
    if (!r.ok) throw new Error(`NCBI returned HTTP ${r.status}.`);
    return r.text();
  }
  function dlgNCBI() {
    let items = [], busy = false;
    const setStatus = (msg, isErr) => { const e = $('#ncstat'); if (e) { e.textContent = msg; e.className = isErr ? 'err' : 'hint'; } };
    const download = async id => {
      if (busy) return; busy = true;
      const auto = $('#ncauto') && $('#ncauto').checked;
      setStatus(`Downloading ${id} from NCBI…`);
      try {
        const text = await ncbiFetch(val('ncdb'), id);
        const src = DNA.parseAny(text, id + '.gb');
        src.histLabel = `Downloaded ${id} from NCBI`;
        let detected = 0;
        // Many database records carry no annotation beyond the source line.
        if (auto) {
          for (const f of DNA.detectCommon(src.seq, src.circular, src.features)) {
            src.features.push({ ...f, color: TYPE_COL[f.type] || '#c9ced6', notes: 'Detected automatically' });
            detected++;
          }
        }
        closeDialog();
        addDoc(src);
        toast(`${src.name} — ${fmt(src.seq.length)} bp, ${src.features.length} features${detected ? ` (${detected} detected)` : ''}`);
      } catch (err) { setStatus(err.message, true); }
      busy = false;
    };
    const renderList = () => {
      $('#ncres').innerHTML = !items.length ? '' : `<table class="grid"><thead><tr><th>Accession</th><th>Description</th><th>Length</th><th>Form</th><th></th></tr></thead><tbody>${items.map(i => `<tr data-id="${esc(i.accession)}"><td><b>${esc(i.accession)}</b></td><td class="wrap">${esc(i.title)}</td><td>${fmt(i.slen)} bp</td><td>${esc(i.topology || '')}</td><td><button class="ib" data-get="${esc(i.accession)}">Open</button></td></tr>`).join('')}</tbody></table>`;
    };
    dialog({
      title: 'Open from NCBI', wide: true,
      body: `<div class="form">
        <label for="ncdb">Database</label><select id="ncdb"><option value="nuccore">Nucleotide (GenBank)</option><option value="nucleotide">Nucleotide (alias)</option></select>
        <label for="ncq">Accession</label><input type="text" id="ncq" placeholder="e.g. L09137, NM_000546, U49845, or a search like: pUC19 complete" spellcheck="false">
        <label></label><label><input type="checkbox" id="ncauto" checked> Detect common features after downloading</label>
        <div class="full hint" id="ncstat">Type an accession number and choose Download, or type words and choose Search.</div>
        <div class="full" id="ncres"></div></div>`,
      buttons: [{ label: 'Close' }, {
        label: 'Search', act: async () => {
          const term = val('ncq').trim(); if (!term) return setStatus('Enter an accession number or search words.', true);
          setStatus('Searching NCBI…'); items = [];
          try { const r = await ncbiSearch(val('ncdb'), term); items = r.items; renderList(); setStatus(items.length ? `${items.length} records found. Choose Open on one of them.` : 'Nothing matched that search.', !items.length); }
          catch (err) { setStatus(err.message, true); }
          return false;
        }
      }, {
        label: 'Download', primary: true, act: () => {
          const term = val('ncq').trim(); if (!term) { setStatus('Enter an accession number.', true); return false; }
          download(term); return false;
        }
      }],
      onOpen: box => {
        box.querySelector('#ncres').addEventListener('click', e => { const b = e.target.closest('[data-get]'); if (b) download(b.dataset.get); });
      }
    });
  }

  // ================================================================ sequencing traces
  function openTraces() {
    const d = D(); if (!d) { toast('Open the reference sequence first.'); return; }
    if (host && host.openFiles) { host.openFiles().then(files => { for (const f of files || []) { const bin = Uint8Array.from(atob(f.data), c => c.charCodeAt(0)); addTrace(f.name, bin.buffer); } }); return; }
    const inp = document.createElement('input');
    inp.type = 'file'; inp.multiple = true; inp.accept = '.ab1,.fa,.fasta,.seq,.txt';
    inp.onchange = () => { for (const f of inp.files) { const r = new FileReader(); r.onload = () => addTrace(f.name, r.result); r.readAsArrayBuffer(f); } };
    inp.click();
  }
  function addTrace(name, buf) {
    const d = D(); if (!d) return;
    try {
      const read = TRACE.parseRead(buf, name);
      const al = TRACE.align(d, read);
      al.id = uid();
      al.colQi = [];
      let qi = 0;
      for (let k = 0; k < al.readAln.length; k++) { al.colQi.push(al.readAln[k] === '-' ? -1 : qi++); }
      d.alignments.push(al);
      S.alignSel = al.id;
      setView('align');
      toast(`${read.name}: ${al.identity.toFixed(1)}% identity over ${fmt(al.refAln.length)} bases, ${al.diffs.length} difference${al.diffs.length === 1 ? '' : 's'}`);
    } catch (err) {
      dialog({ title: 'Could Not Align Trace', body: `<p>${esc(name)}</p><p class="err">${esc(err.message)}</p><p class="hint">Traces are matched to the open sequence on both strands after trimming low-quality ends. A read from a different construct will not align.</p>`, buttons: [{ label: 'OK', primary: true }] });
    }
  }

  const TRACE_COL = { A: '#1a9850', C: '#2563c9', G: '#1c2230', T: '#c43c3c', N: '#8d95a3' };
  function chromaSVG(al, c0, c1, h) {
    const r = al.read;
    if (!r.hasTrace) return '';
    const W = (c1 - c0) * CW;
    // only use peaks the trace arrays actually reach, so missing data draws nothing
    let tlen = Infinity;
    for (const b of 'ACGT') tlen = Math.min(tlen, (r.traces[b] || []).length);
    const cols = [];
    for (let k = c0; k < c1; k++) {
      const qi = al.colQi[k];
      if (qi < 0) { cols.push(null); continue; }
      const raw = al.toRaw(qi);
      cols.push(raw >= 0 && raw < r.peaks.length && r.peaks[raw] < tlen ? raw : null);
    }
    const present = cols.filter(c => c !== null);
    if (present.length < 2) return '';
    // average spacing between peaks sets how much trace one column shows
    let step = 12;
    for (let i = 1; i < present.length; i++) { const dd = Math.abs(r.peaks[present[i]] - r.peaks[present[i - 1]]); if (dd > 0 && dd < 60) { step = dd; break; } }
    let max = 1;
    for (const raw of present) for (const b of 'ACGT') { const t = r.traces[b]; if (t && t[r.peaks[raw]] > max) max = t[r.peaks[raw]]; }
    max *= 1.15;
    const paths = {};
    for (const b of 'ACGT') paths[b] = [];
    const SAMP = 6;
    for (let i = 0; i < cols.length; i++) {
      const raw = cols[i]; if (raw === null) continue;
      const peak = r.peaks[raw];
      for (const b of 'ACGT') {
        const chan = al.strand === 1 ? b : DNA.COMP[b];
        const t = r.traces[chan]; if (!t) continue;
        for (let s = 0; s <= SAMP; s++) {
          const frac = s / SAMP - 0.5;
          const pos = Math.round(peak + (al.strand === 1 ? frac : -frac) * step);
          const v = pos >= 0 && pos < t.length ? t[pos] : 0;
          const x = (i + 0.5 + frac) * CW;
          const y = h - Math.max(0, Math.min(1, v / max)) * (h - 2) - 1;
          paths[b].push(`${x.toFixed(1)},${y.toFixed(1)}`);
        }
      }
    }
    const lines = Object.entries(paths).filter(([, p]) => p.length > 1)
      .map(([b, p]) => `<polyline points="${p.join(' ')}" fill="none" stroke="${TRACE_COL[b]}" stroke-width="1" opacity=".85"/>`).join('');
    return `<svg class="chroma" width="${W}" height="${h}" viewBox="0 0 ${W} ${h}">${lines}</svg>`;
  }

  function renderAlign() {
    const d = D(), v = $('#view');
    if (!d.alignments.length) {
      v.innerHTML = `<div class="empty"><h2>No sequencing traces</h2><p>Open one or more .ab1 files to align them against <b>${esc(d.name)}</b>. Low-quality ends are trimmed, both strands are tried, and differences from the reference are listed.</p><div class="btns"><button class="btn pri" data-a="addtrace">Add Trace Files…</button></div></div>`;
      v.querySelector('[data-a]').onclick = openTraces;
      return;
    }
    if (!d.alignments.some(a => a.id === S.alignSel)) S.alignSel = d.alignments[0].id;
    const al = d.alignments.find(a => a.id === S.alignSel);
    const sum = TRACE.summarize(d.alignments);
    v.innerHTML = `<div class="panel">
      <div class="ptools">
        <select id="alsel">${d.alignments.map(a => `<option value="${a.id}" ${a.id === S.alignSel ? 'selected' : ''}>${esc(a.read.name)} (${a.identity.toFixed(1)}%)</option>`).join('')}</select>
        <button class="btn" id="aladd">Add Trace…</button><button class="btn dng" id="aldel">Remove</button>
        <span class="note">${esc(al.read.name)} · ${al.strand > 0 ? 'forward' : 'reverse'} strand · covers ${fmt(al.refStart + 1)}–${fmt(al.refEnd)} · ${al.identity.toFixed(1)}% identity · ${al.diffs.length} difference${al.diffs.length === 1 ? '' : 's'}${al.trim.trimmed ? ` · trimmed ${al.trim.start} + ${al.read.seq.length - al.trim.end} low-quality bases` : ''}${al.read.hasTrace ? '' : ' · no chromatogram in this file'}</span>
      </div>
      ${sum.length ? `<div class="ptools" style="background:var(--win)"><span class="note">Differences: ${sum.map(s => `<a href="#" data-go="${s.pos}" class="difflink">${fmt(s.pos + 1)} ${s.type === 'mis' ? esc(s.ref + '→' + s.read) : (s.type === 'del' ? 'deleted ' + esc(s.ref) : 'inserted ' + esc(s.read))}${s.reads.length > 1 ? ` (${s.reads.length} reads)` : ''}</a>`).join(' · ')}</span></div>` : `<div class="ptools" style="background:var(--win)"><span class="note">No differences from the reference in any aligned read.</span></div>`}
      <div class="tscroll" id="alrows"></div></div>`;
    buildAlignRows(al);
    $('#alsel').onchange = e => { S.alignSel = +e.target.value; renderView(); };
    $('#aladd').onclick = openTraces;
    $('#aldel').onclick = () => { d.alignments = d.alignments.filter(a => a.id !== S.alignSel); renderAll(); };
    v.querySelectorAll('.difflink').forEach(a => a.onclick = e => {
      e.preventDefault();
      const p = +a.dataset.go;
      setSel(p, p + 1, { anchor: p });
      const el = $('#alrows').querySelector(`[data-ref="${p}"]`);
      if (el) el.scrollIntoView({ block: 'center' });
    });
  }

  function buildAlignRows(al) {
    const d = D(), n = d.seq.length, host2 = $('#alrows');
    const G = 94, per = Math.max(20, Math.floor((host2.clientWidth - G - 26) / CW / 10) * 10);
    const hasT = al.read.hasTrace;
    const chH = hasT ? 54 : 0;
    const out = [];
    // reference coordinate at each alignment column
    const refPos = new Array(al.refAln.length);
    let rp = al.refStart;
    for (let k = 0; k < al.refAln.length; k++) { refPos[k] = rp % n; if (al.refAln[k] !== '-') rp++; }
    for (let c0 = 0; c0 < al.refAln.length; c0 += per) {
      const c1 = Math.min(al.refAln.length, c0 + per);
      const ref = al.refAln.slice(c0, c1), read = al.readAln.slice(c0, c1);
      let readHtml = '', qualHtml = '';
      for (let k = c0; k < c1; k++) {
        const r = al.refAln[k], s = al.readAln[k];
        const cls = s === '-' ? 'd' : (r === '-' ? 'i' : (r === s ? '' : 'm'));
        readHtml += cls ? `<span class="${cls}">${s}</span>` : s;
        const qi = al.colQi[k];
        if (qi >= 0 && al.read.quals.length) {
          const qv = al.read.quals[al.toRaw(qi)] || 0;
          const o = Math.max(0.08, Math.min(1, qv / 45));
          qualHtml += `<div class="qb" style="left:${(k - c0) * CW}px;width:${CW - 0.5}px;height:${(2 + o * 6).toFixed(1)}px;opacity:${(0.25 + o * 0.6).toFixed(2)}"></div>`;
        }
      }
      const y0 = 4;
      out.push(`<div class="arow" style="height:${y0 + 54 + chH}px">
        <div class="num" style="top:${y0 + 1}px">${fmt(refPos[c0] + 1)}</div>
        <div class="lbl2" style="top:${y0}px">ref</div>
        <div class="nt" style="left:${G}px;top:${y0}px">${ref}</div>
        <div class="lbl2" style="top:${y0 + 18}px">read</div>
        <div class="nt rd" style="left:${G}px;top:${y0 + 18}px">${readHtml}</div>
        <div class="quals" style="left:${G}px;top:${y0 + 36}px">${qualHtml}</div>
        ${hasT ? `<div class="chw" style="left:${G}px;top:${y0 + 46}px">${chromaSVG(al, c0, c1, chH)}</div>` : ''}
        ${Array.from({ length: c1 - c0 }, (_, i) => al.refAln[c0 + i] !== '-' && al.readAln[c0 + i] !== al.refAln[c0 + i] ? `<div class="mk" data-ref="${refPos[c0 + i]}" style="left:${G + i * CW}px;top:${y0}px"></div>` : '').join('')}
      </div>`);
    }
    host2.innerHTML = out.join('');
    host2.onclick = e => {
      const row = e.target.closest('.arow'); if (!row) return;
      const rect = row.getBoundingClientRect();
      const i = Math.floor((e.clientX - rect.left - G) / CW);
      const idx = [...host2.children].indexOf(row) * per + i;
      if (idx >= 0 && idx < refPos.length) { const p = refPos[idx]; setSel(p, p + 1, { anchor: p }); }
    };
  }

  // A record of one piece that went into a construct, kept with its own sequence so the
  // history stands on its own even after the parent file is closed.
  function pieceOfFrag(f, detail) {
    return {
      name: `${f.doc.name} ${f.lenz}–${f.renz} fragment`, seq: f.seq, circular: false,
      features: f.doc.features.filter(x => {
        const n = f.doc.seq.length, rs = ((x.start - f.start) % n + n) % n;
        return rs + DNA.rlen(x, n) <= f.size;
      }).map(x => { const n = f.doc.seq.length; const rs = ((x.start - f.start) % n + n) % n; return { ...x, start: rs, end: rs + DNA.rlen(x, n) }; }),
      detail,
      source: { method: 'Cut from ' + f.doc.name, detail: `${f.lenz} (${fmt(f.start + 1)}) – ${f.renz} (${fmt(f.end)})`, inputs: [{ name: f.doc.name, seq: f.doc.seq, circular: f.doc.circular, features: f.doc.features }] }
    };
  }

  // ================================================================ cloning dialogs
  function docOptions(sel) { return S.docs.map((x, i) => `<option value="${i}" ${i === sel ? 'selected' : ''}>${esc(x.name)} (${fmt(x.seq.length)} bp, ${x.circular ? 'circular' : 'linear'})</option>`).join(''); }
  function cutters(doc) {
    const sites = DNA.findSites(doc.seq, doc.circular);
    return DNA.ENZYMES.filter(e => (sites[e.name] || []).some(s => s.cut != null)).map(e => ({ name: e.name, n: sites[e.name].length })).sort((a, b) => a.n - b.n || a.name.localeCompare(b.name));
  }
  function fragLabel(f) { return `${fmt(f.size)} bp · ${esc(f.lenz)} / ${esc(f.renz)} · ends ${esc(CLONE.endLabel(f.left))} … ${esc(CLONE.endLabel(f.right))}`; }

  function dlgRestrictionClone() {
    if (!S.docs.length) return;
    let vi = S.active, ii = S.docs.length > 1 ? (S.active + 1) % S.docs.length : S.active;
    let vEnz = [], iEnz = [], products = [];
    const enzBox = (id, doc, chosen) => `<div class="enz-pick" id="${id}" style="max-height:110px">${cutters(doc).map(e => `<label><input type="checkbox" value="${e.name}" ${chosen.includes(e.name) ? 'checked' : ''}> ${e.name} <span style="color:var(--faint)">(${e.n})</span></label>`).join('') || '<span class="hint">No enzyme from the list cuts this sequence.</span>'}</div>`;
    const guess = doc => { const c = cutters(doc).filter(e => e.n === 1 && DNA.ENZ_BY_NAME[e.name].core >= 6); return c.slice(0, 2).map(e => e.name); };
    const fragsOf = (doc, enz) => enz.length ? CLONE.cut(doc, DNA.findSites(doc.seq, doc.circular), enz) : [];
    // The enzyme lists are rebuilt only when the document changes, so ticking a box keeps its place.
    const refreshEnz = () => {
      $('#vbox').innerHTML = enzBox('venz', S.docs[vi], vEnz);
      $('#ibox').innerHTML = enzBox('ienz', S.docs[ii], iEnz);
    };
    const refresh = () => {
      const V = S.docs[vi], I = S.docs[ii];
      const vf = fragsOf(V, vEnz), inf = fragsOf(I, iEnz);
      const usable = f => f.left.type !== 'o' && f.right.type !== 'o';
      // Largest first: the vector backbone is almost always the big fragment.
      const bySize = (a, b) => b.size - a.size;
      const vlist = vf.filter(usable).sort(bySize), ilist = inf.filter(usable).sort(bySize);
      $('#vfr').innerHTML = vlist.length ? vlist.map((f, k) => `<label><input type="radio" name="vf" value="${k}" ${k === 0 ? 'checked' : ''}> ${fragLabel(f)}</label>`).join('') : `<span class="hint">${vEnz.length ? 'These enzymes do not produce a fragment with two cut ends.' : 'Choose the enzymes that cut the vector.'}</span>`;
      $('#ifr').innerHTML = ilist.length ? ilist.map((f, k) => `<label><input type="radio" name="if" value="${k}" ${k === 0 ? 'checked' : ''}> ${fragLabel(f)}</label>`).join('') : `<span class="hint">${iEnz.length ? 'These enzymes do not produce a fragment with two cut ends.' : 'Choose the enzymes that cut the insert.'}</span>`;
      $('#rcres').innerHTML = '';
      products = [];
      refresh.vlist = vlist; refresh.ilist = ilist;
    };
    const simulate = () => {
      const vk = +(document.querySelector('input[name=vf]:checked') || {}).value;
      const ik = +(document.querySelector('input[name=if]:checked') || {}).value;
      const V = refresh.vlist[vk], I = refresh.ilist[ik];
      if (!V || !I) { $('#rcres').innerHTML = '<div class="err">Choose one vector fragment and one insert fragment.</div>'; return; }
      products = CLONE.ligate([V, I], { name: `${S.docs[vi].name}-${S.docs[ii].name}` });
      if (!products.length) {
        $('#rcres').innerHTML = `<div class="err">These two fragments cannot be joined. The vector ends are ${esc(CLONE.endLabel(V.left))} and ${esc(CLONE.endLabel(V.right))}; the insert ends are ${esc(CLONE.endLabel(I.left))} and ${esc(CLONE.endLabel(I.right))}. Ligation needs matching overhangs at both joins.</div>`;
        return;
      }
      $('#rcres').innerHTML = `<table class="grid"><thead><tr><th>Product</th><th>Size</th><th>Orientation</th><th>Junctions</th><th></th></tr></thead><tbody>${products.map((p, k) => `<tr><td>${k + 1}</td><td><b>${fmt(p.size)} bp</b></td><td>${p.chain[1].rev ? 'insert reversed' : 'insert forward'}</td><td class="wrap">${esc(p.junctions.join(' · '))}</td><td><button class="ib" data-mk="${k}">Create</button></td></tr>`).join('')}</tbody></table><div class="hint">Features on both fragments are carried into the product.</div>`;
    };
    vEnz = guess(S.docs[vi]); iEnz = guess(S.docs[ii]);
    dialog({
      title: 'Restriction Cloning', wide: true,
      body: `<div class="form">
        <label for="vdoc">Vector</label><select id="vdoc">${docOptions(vi)}</select>
        <label>Cut with</label><div id="vbox"></div>
        <label>Fragment</label><div class="fraglist" id="vfr"></div>
        <div class="full" style="border-top:1px solid var(--border-soft);margin:4px 0"></div>
        <label for="idoc">Insert from</label><select id="idoc">${docOptions(ii)}</select>
        <label>Cut with</label><div id="ibox"></div>
        <label>Fragment</label><div class="fraglist" id="ifr"></div>
        <div class="full" id="rcres"></div></div>`,
      buttons: [{ label: 'Close' }, { label: 'Simulate Ligation', primary: true, act: () => { simulate(); return false; } }],
      onOpen: box => {
        refreshEnz(); refresh();
        box.addEventListener('change', e => {
          if (e.target.id === 'vdoc') { vi = +e.target.value; vEnz = guess(S.docs[vi]); refreshEnz(); refresh(); }
          else if (e.target.id === 'idoc') { ii = +e.target.value; iEnz = guess(S.docs[ii]); refreshEnz(); refresh(); }
          else if (e.target.closest('#venz')) { vEnz = [...box.querySelectorAll('#venz input:checked')].map(i => i.value); refresh(); }
          else if (e.target.closest('#ienz')) { iEnz = [...box.querySelectorAll('#ienz input:checked')].map(i => i.value); refresh(); }
        });
        box.addEventListener('click', e => {
          const b = e.target.closest('[data-mk]'); if (!b) return;
          const p = products[+b.dataset.mk];
          closeDialog();
          addDoc({
            ...p.doc, histLabel: `Ligated ${S.docs[vi].name} with ${S.docs[ii].name}`,
            provenance: {
              method: 'Restriction ligation',
              detail: p.junctions.join('  ·  '),
              inputs: p.chain.map((f, k) => pieceOfFrag(f.rev ? { ...f, seq: f.seq } : f, k === 0 ? 'Backbone' : 'Insert' + (f.rev ? ', flipped' : '')))
            }
          });
          toast(`Created ${fmt(p.size)} bp construct`);
        });
      }
    });
  }

  function dlgGibson() {
    let picked = S.docs.map((_, i) => i === S.active), results = [], idxUsed = [];
    const pieceOf = i => {
      const d = S.docs[i];
      const useSel = d.circular && selLen(d) > 0 && i === S.active;
      const seq = useSel ? selSeq(d) : (d.circular ? d.seq : d.seq);
      return { name: d.name + (useSel ? ' (selection)' : ''), seq, doc: d, start: useSel ? d.sel.start : 0, end: useSel ? d.sel.end : d.seq.length };
    };
    const simulate = () => {
      const idx = picked.map((p, i) => p ? i : -1).filter(i => i >= 0);
      if (idx.length < 2) { $('#gbres').innerHTML = '<div class="err">Choose at least two pieces.</div>'; return; }
      if (idx.length > 6) { $('#gbres').innerHTML = '<div class="err">Choose at most six pieces.</div>'; return; }
      const min = Math.max(10, Math.min(80, parseInt(val('gbov'), 10) || 20));
      const pieces = idx.map(pieceOf);
      idxUsed = idx;
      results = CLONE.gibson(pieces, { minOverlap: min, name: 'Assembled construct' });
      if (!results.length) {
        $('#gbres').innerHTML = `<div class="err">These pieces do not share ends of at least ${min} bp, so they cannot be assembled. Gibson assembly needs each piece to end in the same sequence that the next piece begins with — usually added as primer tails.</div>`;
        return;
      }
      $('#gbres').innerHTML = `<table class="grid"><thead><tr><th>Product</th><th>Size</th><th>Junction overlaps</th><th></th></tr></thead><tbody>${results.map((r, k) => `<tr><td>${k + 1}</td><td><b>${fmt(r.size)} bp</b></td><td class="wrap">${r.junctions.map(j => `${esc(j.from)} → ${esc(j.to)}: ${j.len} bp`).join(' · ')}</td><td><button class="ib" data-mk="${k}">Create</button></td></tr>`).join('')}</tbody></table>`;
    };
    dialog({
      title: 'Gibson / Overlap Assembly', wide: true,
      body: `<div class="form">
        <label>Pieces</label><div class="enz-pick" id="gbp" style="columns:1;max-height:150px">${S.docs.map((d, i) => `<label><input type="checkbox" value="${i}" ${picked[i] ? 'checked' : ''}> ${esc(d.name)} — ${fmt(d.seq.length)} bp${d.circular ? ' (circular; the current selection is used if there is one)' : ''}</label>`).join('')}</div>
        <label for="gbov">Minimum overlap</label><div class="row2"><input type="number" id="gbov" value="20" min="10" max="80" style="width:80px"> bp</div>
        <div class="full hint">Each piece must end in the same sequence that the next piece starts with. Pieces are tried in every order and orientation.</div>
        <div class="full" id="gbres"></div></div>`,
      buttons: [{ label: 'Close' }, { label: 'Simulate Assembly', primary: true, act: () => { simulate(); return false; } }],
      onOpen: box => {
        box.querySelector('#gbp').addEventListener('change', () => { picked = S.docs.map((_, i) => !!box.querySelector(`#gbp input[value="${i}"]`).checked); });
        box.addEventListener('click', e => {
          const b = e.target.closest('[data-mk]'); if (!b) return;
          const r = results[+b.dataset.mk];
          closeDialog();
          addDoc({
            ...r.doc, histLabel: 'Assembled from ' + r.junctions.length + ' pieces',
            provenance: {
              method: 'Gibson / overlap assembly',
              detail: r.junctions.map(j => `${j.from} → ${j.to}: ${j.len} bp overlap`).join('  ·  '),
              inputs: idxUsed.map(i => {
                const pc = pieceOf(i);
                return { name: pc.name, seq: pc.seq, circular: false, features: [], detail: 'Overlap and insert 1..' + pc.seq.length,
                  source: { method: 'From ' + pc.doc.name, detail: `${fmt(pc.start + 1)}..${fmt(pc.end)}`, inputs: [{ name: pc.doc.name, seq: pc.doc.seq, circular: pc.doc.circular, features: pc.doc.features }] } };
              })
            }
          });
          toast(`Created ${fmt(r.size)} bp construct`);
        });
      }
    });
  }

  function dlgGateway() {
    if (S.docs.length < 2) { toast('Open both the entry clone and the destination vector first.'); return; }
    const sizes = S.docs.map(d => d.seq.length);
    let di = sizes.indexOf(Math.min(...sizes)), ai = S.docs.findIndex((_, i) => i !== di);
    let res = null;
    const describe = i => {
      const d = S.docs[i], sites = CLONE.findAtt(d.seq, d.circular);
      return sites.length ? sites.map(s => `${esc(s.label)} at ${fmt(s.pos + 1)}`).join(', ') : 'no att sites found';
    };
    const refresh = () => {
      $('#gwd').innerHTML = describe(di);
      $('#gwa').innerHTML = describe(ai);
      $('#gwres').innerHTML = '';
      res = null;
    };
    const simulate = () => {
      if (di === ai) { $('#gwres').innerHTML = '<div class="err">Choose two different sequences.</div>'; return; }
      try {
        res = CLONE.gateway(S.docs[di], S.docs[ai], { name: `${S.docs[ai].name}-${S.docs[di].name}`, byName: 'By-product' });
        const warn = res.junctions.filter(j => !j.startsWith('attB')).length;
        $('#gwres').innerHTML = `<table class="grid"><thead><tr><th>Product</th><th>Size</th><th>Features</th><th></th></tr></thead><tbody>
          <tr><td><b>Recombination product</b></td><td>${fmt(res.product.seq.length)} bp</td><td>${res.product.features.length}</td><td><button class="ib" data-mk="p">Create</button></td></tr>
          <tr><td>By-product</td><td>${fmt(res.byproduct.seq.length)} bp</td><td>${res.byproduct.features.length}</td><td><button class="ib" data-mk="b">Create</button></td></tr>
          </tbody></table><div class="hint">Transferred segment: ${fmt(res.insertSize)} bp. Junctions read as ${esc(res.junctions.join(' and '))}.${warn ? ' A junction that is not a full attB site usually means the att sites are not the standard ones.' : ''}</div>`;
      } catch (err) { $('#gwres').innerHTML = `<div class="err">${esc(err.message)}</div>`; }
    };
    dialog({
      title: 'Gateway Recombination', wide: true,
      body: `<div class="form">
        <label for="gwdon">Insert donor</label><select id="gwdon">${docOptions(di)}</select>
        <label></label><div class="hint" id="gwd" style="margin:0"></div>
        <label for="gwacc">Backbone</label><select id="gwacc">${docOptions(ai)}</select>
        <label></label><div class="hint" id="gwa" style="margin:0"></div>
        <div class="full hint">The stretch between the donor’s att1 and att2 sites replaces the stretch between the backbone’s own att sites. For an LR reaction the donor is the entry clone and the backbone is the destination vector; for BP, the donor carries the attB-flanked insert and the backbone is the donor vector.</div>
        <div class="full" id="gwres"></div></div>`,
      buttons: [{ label: 'Close' }, { label: 'Simulate Reaction', primary: true, act: () => { simulate(); return false; } }],
      onOpen: box => {
        refresh();
        box.addEventListener('change', e => {
          if (e.target.id === 'gwdon') di = +e.target.value;
          else if (e.target.id === 'gwacc') ai = +e.target.value;
          else return;
          refresh();
        });
        box.addEventListener('click', e => {
          const b = e.target.closest('[data-mk]'); if (!b || !res) return;
          const p = b.dataset.mk === 'p' ? res.product : res.byproduct;
          const don = S.docs[di], acc = S.docs[ai];
          closeDialog();
          addDoc({
            ...p, histLabel: 'Gateway recombination',
            provenance: {
              method: 'Gateway recombination',
              detail: `Junctions read as ${res.junctions.join(' and ')}; ${fmt(res.insertSize)} bp transferred`,
              inputs: [
                { name: don.name, seq: don.seq, circular: don.circular, features: don.features, detail: 'Insert donor' },
                { name: acc.name, seq: acc.seq, circular: acc.circular, features: acc.features, detail: 'Backbone' }
              ]
            }
          });
          toast(`Created ${fmt(p.seq.length)} bp construct`);
        });
      }
    });
  }

  function showAttSites() {
    const d = D(), sites = CLONE.findAtt(d.seq, d.circular);
    if (!sites.length) { toast('No att sites found in this sequence.'); return; }
    for (const s of sites) {
      if (d.features.some(f => f.start === s.pos && f.name.startsWith('att'))) continue;
      d.features.push({ id: uid(), name: s.kind, type: 'protein_bind', start: s.pos, end: (s.pos + 7) % d.seq.length, strand: s.strand, color: '#69c4c6', notes: 'Recombination core detected by GeneCanvas' });
    }
    commit(d, `Marked ${sites.length} att site${sites.length > 1 ? 's' : ''}`); renderAll(true);
    toast(sites.map(s => s.label).join(', '));
  }

  // ================================================================ Golden Gate
  function dlgGoldenGate() {
    let enzyme = 'BsaI';
    let picked = S.docs.map((_, i) => i === S.active);
    let result = null;
    const describe = (d, enz) => {
      const n = (DNA.findSites(d.seq, d.circular)[enz] || []).length;
      let frags = [];
      try { frags = CLONE.ggFragments(d, enz); } catch (e) { frags = []; }
      if (!n) return `<span style="color:var(--danger)">no ${enz} site</span>`;
      if (!frags.length) return `<span style="color:var(--danger)">${n} ${enz} site${n > 1 ? 's' : ''}, but no site-free fragment</span>`;
      const f = frags.slice().sort((a, b) => b.size - a.size)[0];
      const inner = n > 2 ? ` <span style="color:var(--danger)">— ${n - 2} site${n - 2 > 1 ? 's' : ''} inside the part</span>` : '';
      return `${n} ${enz} site${n > 1 ? 's' : ''} → part of ${fmt(f.size)} bp, overhangs <span class="mono">${esc(f.left.x)}</span> / <span class="mono">${esc(f.right.x)}</span>${inner}`;
    };
    const refreshList = () => {
      $('#ggp').innerHTML = S.docs.map((d, i) => `<label><input type="checkbox" value="${i}" ${picked[i] ? 'checked' : ''}> <b>${esc(d.name)}</b> — ${describe(d, enzyme)}</label>`).join('');
    };
    const simulate = () => {
      const idx = picked.map((p, i) => p ? i : -1).filter(i => i >= 0);
      if (idx.length < 2) { $('#ggres').innerHTML = '<div class="err">Choose at least two parts.</div>'; return; }
      try {
        result = CLONE.goldenGate(idx.map(i => ({ doc: S.docs[i] })), enzyme, { name: 'Golden Gate product' });
      } catch (err) { $('#ggres').innerHTML = `<div class="err">${esc(err.message)}</div>`; return; }
      const warn = [];
      if (result.duplicateOverhangs.length) warn.push(`Two parts end in the same overhang (${esc(result.duplicateOverhangs.join(', '))}), so the order is ambiguous.`);
      warn.push(...result.notes);
      if (!result.products.length) {
        $('#ggres').innerHTML = `<div class="err">The parts do not form a circle. Their overhangs are ${result.parts.map(f => esc(f.left.x + '…' + f.right.x)).join(', ')}. Each part's right overhang must match the next part's left overhang.</div>`;
        return;
      }
      $('#ggres').innerHTML = `<table class="grid"><thead><tr><th>Product</th><th>Size</th><th>Junction overhangs</th><th>${esc(enzyme)} sites left</th><th></th></tr></thead><tbody>${result.products.map((p, k) => `<tr><td>${k + 1}</td><td><b>${fmt(p.size)} bp</b></td><td class="mono">${p.overhangs.map(esc).join(' · ')}</td><td>${p.leftover === 0 ? '<span style="color:#1a7f37">none</span>' : `<span style="color:var(--danger)">${p.leftover}</span>`}</td><td><button class="ib" data-mk="${k}">Create</button></td></tr>`).join('')}</tbody></table>
        ${warn.length ? `<div class="hint" style="color:var(--danger)">${warn.map(esc).join('<br>')}</div>` : ''}
        <div class="hint">A finished Golden Gate product keeps no ${esc(enzyme)} site, because both sites are carried away on the discarded ends.</div>`;
    };
    dialog({
      title: 'Golden Gate Assembly', wide: true,
      body: `<div class="form">
        <label for="ggenz">Enzyme</label><select id="ggenz">${CLONE.TYPE2S.map(e => `<option ${e === enzyme ? 'selected' : ''}>${e}</option>`).join('')}</select>
        <label>Parts</label><div class="enz-pick" id="ggp" style="columns:1;max-height:180px"></div>
        <div class="full hint">Each part needs two inward-facing ${esc(enzyme)} sites. The piece between them carries no site of its own and is what goes into the product; the overhangs decide the order.</div>
        <div class="full" id="ggres"></div></div>`,
      buttons: [{ label: 'Close' }, { label: 'Simulate Assembly', primary: true, act: () => { simulate(); return false; } }],
      onOpen: box => {
        refreshList();
        box.addEventListener('change', e => {
          if (e.target.id === 'ggenz') { enzyme = e.target.value; refreshList(); $('#ggres').innerHTML = ''; }
          else if (e.target.closest('#ggp')) picked = S.docs.map((_, i) => { const el = box.querySelector(`#ggp input[value="${i}"]`); return !!(el && el.checked); });
        });
        box.addEventListener('click', e => {
          const b = e.target.closest('[data-mk]'); if (!b || !result) return;
          const p = result.products[+b.dataset.mk];
          closeDialog();
          addDoc({
            ...p.doc, histLabel: `Golden Gate with ${enzyme}`,
            provenance: {
              method: `Golden Gate assembly (${enzyme})`,
              detail: 'Junction overhangs ' + p.overhangs.join(' · '),
              inputs: p.chain.map(f => pieceOfFrag(f, `Overhangs ${f.left.x} … ${f.right.x}`))
            }
          });
          toast(`Created ${fmt(p.size)} bp construct`);
        });
      }
    });
  }

  // ================================================================ compare sequences
  function dlgCompare() {
    let picked = S.docs.map((_, i) => i === S.active);
    dialog({
      title: 'Compare Sequences', wide: true,
      body: `<div class="form">
        <label>Sequences</label><div class="enz-pick" id="cmpp" style="columns:1;max-height:200px">${S.docs.map((d, i) => `<label><input type="checkbox" value="${i}" ${picked[i] ? 'checked' : ''}> ${esc(d.name)} — ${fmt(d.seq.length)} bp, ${d.circular ? 'circular' : 'linear'}</label>`).join('')}</div>
        <div class="full hint">The first sequence you tick is the reference; the rest are aligned to it. Circular sequences are rotated and, where needed, reverse-complemented to line up.</div>
        <div class="full err"></div></div>`,
      buttons: [{ label: 'Cancel' }, {
        label: 'Compare', primary: true, act: () => {
          const idx = picked.map((p, i) => p ? i : -1).filter(i => i >= 0);
          if (idx.length < 2) return showErr('Choose at least two sequences.');
          if (idx.length > 8) return showErr('Choose at most eight sequences.');
          const total = idx.reduce((a, i) => a + S.docs[i].seq.length, 0);
          if (total > 600000) return showErr('These sequences are too long to compare in one go.');
          try {
            S.compare = MSA.compare(idx.map(i => ({ name: S.docs[i].name, seq: S.docs[i].seq, circular: S.docs[i].circular })));
            S.compareDiffOnly = false;
          } catch (err) { return showErr(err.message); }
          closeDialog();
          S.view = 'compare';
          renderAll();
        }
      }],
      onOpen: box => {
        box.querySelector('#cmpp').addEventListener('change', () => {
          picked = S.docs.map((_, i) => { const el = box.querySelector(`#cmpp input[value="${i}"]`); return !!(el && el.checked); });
        });
      }
    });
  }

  function renderCompare() {
    const v = $('#view'), C = S.compare;
    if (!C) {
      v.innerHTML = `<div class="empty"><h2>No comparison yet</h2><p>Compare two or more open sequences to see every substitution, insertion and deletion between them.</p><div class="btns"><button class="btn pri" data-a="cmp">Compare Sequences…</button></div></div>`;
      v.querySelector('[data-a]').onclick = dlgCompare;
      return;
    }
    v.innerHTML = `<div class="panel">
      <div class="ptools">
        <button class="btn" id="cmpnew">Compare Other Sequences…</button>
        <label><input type="checkbox" id="cmpdiff" ${S.compareDiffOnly ? 'checked' : ''}> Show only rows with differences</label>
        <span class="note">Reference <b>${esc(C.refName)}</b> · ${fmt(C.refLength)} bp · ${fmt(C.columns)} alignment columns</span>
      </div>
      <div class="ptools" style="background:var(--win);display:block">
        <table class="grid" style="width:auto"><thead><tr><th>Sequence</th><th>Length</th><th>Strand</th><th>Identity</th><th>Substitutions</th><th>Insertions</th><th>Deletions</th></tr></thead><tbody>
        ${C.stats.map(s => `<tr><td><b>${esc(s.name)}</b></td><td>${fmt(s.length)} bp</td><td>${s.strand > 0 ? 'same' : 'reverse complement'}${s.rotation ? `, rotated ${fmt(s.rotation)}` : ''}</td><td>${s.identity.toFixed(2)}%</td><td>${fmt(s.substitutions)}</td><td>${fmt(s.insertions)}</td><td>${fmt(s.deletions)}</td></tr>`).join('')}
        </tbody></table></div>
      <div class="tscroll" id="cmprows"></div></div>`;
    buildCompareRows();
    $('#cmpnew').onclick = dlgCompare;
    $('#cmpdiff').onchange = e => { S.compareDiffOnly = e.target.checked; buildCompareRows(); };
  }

  function buildCompareRows() {
    const C = S.compare, host2 = $('#cmprows');
    const labelW = Math.min(170, 24 + Math.max(...[C.refName, ...C.names].map(n => n.length)) * 6.6);
    const G = labelW + 62;
    const per = Math.max(20, Math.floor((host2.clientWidth - G - 26) / CW / 10) * 10);
    const rows = [C.refRow, ...C.rows], names = [C.refName, ...C.names];
    const out = [];
    for (let c0 = 0; c0 < C.columns; c0 += per) {
      const c1 = Math.min(C.columns, c0 + per);
      let anyDiff = false;
      for (let c = c0; c < c1 && !anyDiff; c++) for (let r = 1; r < rows.length; r++) if (rows[r][c] !== rows[0][c]) anyDiff = true;
      if (S.compareDiffOnly && !anyDiff) continue;
      let pos = null;
      for (let c = c0; c < c1 && pos === null; c++) if (C.refPos[c] !== null) pos = C.refPos[c];
      const lines = rows.map((row, r) => {
        let html = '';
        for (let c = c0; c < c1; c++) {
          const ch = row[c], ref = rows[0][c];
          if (r === 0) { html += ch; continue; }
          if (ch === ref) html += '<span class="s">.</span>';
          else if (ch === '-') html += '<span class="d">-</span>';
          else if (ref === '-') html += `<span class="i">${ch}</span>`;
          else html += `<span class="m">${ch}</span>`;
        }
        return `<div class="cname" style="top:${4 + r * 17}px;width:${labelW}px">${esc(names[r])}</div><div class="nt${r ? ' alt' : ''}" style="left:${G}px;top:${4 + r * 17}px">${html}</div>`;
      }).join('');
      out.push(`<div class="crow" style="height:${8 + rows.length * 17}px"><div class="num" style="left:${labelW}px;top:5px">${pos === null ? '' : fmt(pos + 1)}</div>${lines}</div>`);
    }
    host2.innerHTML = out.length ? out.join('') : '<div class="empty"><p>No differences between these sequences.</p></div>';
  }

  // ================================================================ codon optimisation
  function dlgCodon() {
    const d = D(), n = d.seq.length;
    const orfs = analyze(d).orfs;
    const cdsFeatures = d.features.filter(f => f.type === 'CDS' || f.type === 'gene');
    let mode = 'region';
    let preserve = [];
    let avoid = new Set(), add = new Set();
    let result = null, lastCds = '';

    const regionSeq = () => {
      const a = parseInt(val('coFrom'), 10) - 1, b = parseInt(val('coTo'), 10);
      if (!(a >= 0 && b >= 1 && a < n && b <= n)) return '';
      let s = DNA.sub(d.seq, a, b);
      if ($('#coRev').checked) s = DNA.rc(s);
      return s;
    };
    const currentCds = () => mode === 'protein' ? '' : regionSeq();

    const checkORF = () => {
      const box = $('#coOrf'); if (!box) return;
      if (mode === 'protein') {
        const p = val('coProt').toUpperCase().replace(/[^A-Z*]/g, '');
        box.innerHTML = p ? `<span style="color:#1a7f37">✓</span> ${fmt(p.length)} amino acids. A coding sequence will be built for the chosen organism.` : 'Paste a protein sequence.';
        return;
      }
      const s = regionSeq();
      lastCds = s;
      if (!s) { box.innerHTML = '<span style="color:var(--danger)">Enter a range inside the sequence.</span>'; return; }
      const r = CODON.inspectORF(s);
      if (r.ok) box.innerHTML = `<span style="color:#1a7f37">✓ Valid reading frame</span> — ${fmt(r.codons)} codons, ${fmt(r.aa)} amino acids, ${esc(r.startsWith)} … ${esc(r.stopsWith)}`;
      else box.innerHTML = `<span style="color:var(--danger)">Not a clean reading frame:</span> ${esc(r.issues.join('; '))}.<br><span style="color:var(--muted)">${fmt(r.codons)} codons read as ${esc(r.protein.slice(0, 40))}${r.protein.length > 40 ? '…' : ''}</span>`;
    };
    const setRange = (a, b, rev) => { $('#coFrom').value = a + 1; $('#coTo').value = b; $('#coRev').checked = !!rev; checkORF(); refreshEnzymeHints(); };

    const refreshPreserve = () => {
      $('#coPres').innerHTML = preserve.length
        ? preserve.map((p, i) => `<label><span class="mono">${fmt(p.start + 1)}–${fmt(p.end)}</span> <span style="color:var(--muted)">${esc(p.name || '')}</span> <button class="ib dng" data-rmp="${i}">remove</button></label>`).join('')
        : '<span class="hint">Nothing preserved. Codons in a preserved range are left exactly as they are.</span>';
    };
    const refreshEnzymeHints = () => {
      const s = currentCds();
      const el = $('#coPresent'); if (!el) return;
      if (!s) { el.textContent = ''; return; }
      const present = DNA.ENZYMES.filter(e => e.core >= 6 && CODON.findSite(s, e.site).length).map(e => e.name);
      el.innerHTML = present.length ? `Sites already inside this region: ${present.map(x => `<a href="#" data-avoid="${esc(x)}">${esc(x)}</a>`).join(', ')}` : 'No common restriction site inside this region.';
    };

    const run = () => {
      const organism = val('coOrg');
      const opts = {
        organism, algorithm: val('coAlgo'),
        sourceOrganism: val('coAlgo') === 'harmonised' ? val('coSrc') : null,
        rareCutoff: Math.max(0, Math.min(0.5, (parseFloat(val('coRare')) || 10) / 100)),
        avoidEnzymes: [...avoid], addEnzymes: [...add],
        maxRun: parseInt(val('coRun'), 10) || 0,
        keepStart: $('#coKeepStart').checked, keepStop: $('#coKeepStop').checked,
        seed: parseInt(val('coSeed'), 10) || 1
      };
      try {
        if (mode === 'protein') {
          result = CODON.fromProtein(val('coProt'), opts);
        } else {
          const s = regionSeq();
          if (!s) throw new Error('Enter a range inside the sequence.');
          const r = CODON.inspectORF(s);
          if (!r.ok && !$('#coForce').checked) throw new Error(`This is not a clean reading frame (${r.issues.join('; ')}). Fix the range, or tick “optimise anyway”.`);
          const off = parseInt(val('coFrom'), 10) - 1;
          opts.preserve = preserve.map(p => $('#coRev').checked
            ? { start: s.length - (p.end - off), end: s.length - (p.start - off) }
            : { start: p.start - off, end: p.end - off });
          result = CODON.optimize(s, opts);
        }
      } catch (err) { $('#coRes').innerHTML = `<div class="err">${esc(err.message)}</div>`; result = null; return; }
      showResult();
    };
    const showResult = () => {
      const r = result, b = r.before, a = r.after;
      const row = (label, x, y, fmtFn, better) => {
        const f = fmtFn || (v => v.toFixed(2));
        return `<tr><td>${label}</td><td>${f(x)}</td><td><b>${f(y)}</b></td><td>${better ? better(x, y) : ''}</td></tr>`;
      };
      const arrow = (x, y) => y > x ? '<span style="color:#1a7f37">↑</span>' : (y < x ? '<span style="color:var(--danger)">↓</span>' : '–');
      const seqHtml = [];
      for (let i = 0; i < r.seq.length; i += 60) seqHtml.push(esc(String(i + 1).padStart(6) + '  ' + r.seq.slice(i, i + 60)));
      $('#coRes').innerHTML = `
        <table class="grid" style="width:auto"><thead><tr><th>Measure</th><th>Before</th><th>After</th><th></th></tr></thead><tbody>
          ${row('CAI (codon adaptation index)', b.cai, a.cai, v => v.toFixed(3), arrow)}
          ${row('GC content', b.gc, a.gc, v => v.toFixed(1) + '%')}
          ${row('GC at third positions', b.gc3, a.gc3, v => v.toFixed(1) + '%')}
          ${row('Longest single-base run', b.run, a.run, v => String(v))}
          ${row('Unwanted sites', b.sites, a.sites, v => String(v))}
        </tbody></table>
        <div class="hint">${fmt(r.changedCodons)} of ${fmt(r.totalCodons)} codons changed (${(r.changedCodons / r.totalCodons * 100).toFixed(0)}%) for <b>${esc(r.organism.name)}</b>. CAI 1.00 means every codon is the host’s commonest choice.${a.gc3 > 70 ? ' <span style="color:var(--danger)">GC at third positions is very high, which can make the gene hard to synthesise.</span>' : ''}</div>
        ${r.remaining.length ? `<div class="err">Still present: ${r.remaining.map(x => `${esc(x.name)} ×${x.count}`).join(', ')}</div>` : ''}
        ${r.log.length ? `<div class="hint">${r.log.map(esc).join('<br>')}</div>` : ''}
        <div style="margin-top:8px"><b>Optimised sequence</b></div>
        <div class="aa-box" style="max-height:170px">${seqHtml.join('<br>')}</div>`;
      $('#coApply').disabled = false;
      $('#coOpen').disabled = false;
      $('#coCopy').disabled = false;
    };

    const algoNote = () => {
      const a = CODON.ALGORITHMS.find(x => x.key === val('coAlgo'));
      $('#coAlgoNote').textContent = a ? a.note : '';
      $('#coSrcRow').style.display = val('coAlgo') === 'harmonised' ? '' : 'none';
      $('#coRareRow').style.display = val('coAlgo') === 'threshold' ? '' : 'none';
    };

    const orgOptions = CODON.tables().map(t => `<option value="${t.key}" ${t.key === 'ecoli' ? 'selected' : ''}>${esc(t.name)}</option>`).join('');
    const enzChoices = DNA.ENZYMES.filter(e => e.core >= 6).map(e => e.name);
    const sel = selLen(d);

    dialog({
      title: 'Codon Optimisation', wide: true,
      body: `<div class="form">
        <label>Optimise</label><div class="row2">
          <label><input type="radio" name="comode" value="region" checked> a region of ${esc(d.name)}</label>
          <label><input type="radio" name="comode" value="protein"> a protein sequence</label></div>

        <div class="full" id="coRegionBox"><div class="form" style="margin:0">
          <label for="coFrom">Coding region</label><div class="row2">
            <input type="number" id="coFrom" min="1" max="${n}" value="${sel ? d.sel.start + 1 : 1}"> to
            <input type="number" id="coTo" min="1" max="${n}" value="${sel ? (d.sel.end || n) : Math.min(n, 300)}">
            <label><input type="checkbox" id="coRev"> reverse strand</label></div>
          <label>Use</label><div class="row2">
            <select id="coPick"><option value="">choose a feature or ORF…</option>
              ${cdsFeatures.map((f, i) => `<option value="f${i}">${esc(f.name)} (${esc(f.type)}, ${fmt(DNA.rlen(f, n))} bp)</option>`).join('')}
              ${orfs.map((o, i) => `<option value="o${i}">ORF ${fmt(o.start + 1)}–${fmt(o.end)} (${o.aa} aa, ${o.strand > 0 ? 'forward' : 'reverse'})</option>`).join('')}
            </select></div>
          <label>Reading frame</label><div class="full" id="coOrf" style="grid-column:2/-1"></div>
          <label></label><label><input type="checkbox" id="coForce"> Optimise anyway, even if it is not a clean reading frame</label>
        </div></div>

        <div class="full" id="coProteinBox" style="display:none"><div class="form" style="margin:0">
          <label for="coProt">Protein</label><textarea id="coProt" rows="4" placeholder="Paste a protein sequence in one-letter code"></textarea>
        </div></div>

        <label for="coOrg">Target organism</label><select id="coOrg">${orgOptions}</select>
        <label></label><div class="hint" id="coOrgNote" style="margin:0"></div>

        <label for="coAlgo">Algorithm</label><select id="coAlgo">${CODON.ALGORITHMS.map(a => `<option value="${a.key}" ${a.key === 'threshold' ? 'selected' : ''}>${esc(a.name)}</option>`).join('')}</select>
        <label></label><div class="hint" id="coAlgoNote" style="margin:0"></div>
        <div class="full" id="coRareRow"><div class="form" style="margin:0"><label for="coRare">Drop codons below</label><div class="row2"><input type="number" id="coRare" value="10" min="0" max="50"> % of their amino acid</div></div></div>
        <div class="full" id="coSrcRow" style="display:none"><div class="form" style="margin:0"><label for="coSrc">Source organism</label><select id="coSrc">${CODON.tables().map(t => `<option value="${t.key}" ${t.key === 'human' ? 'selected' : ''}>${esc(t.name)}</option>`).join('')}</select></div></div>

        <label>Remove sites</label><div><input type="text" id="coAvoidQ" placeholder="filter enzymes" style="width:150px"><div class="enz-pick" id="coAvoid" style="max-height:92px"></div><div class="hint" id="coPresent" style="margin:2px 0 0"></div></div>
        <label>Add sites</label><div><input type="text" id="coAddQ" placeholder="filter enzymes" style="width:150px"><div class="enz-pick" id="coAdd" style="max-height:92px"></div><div class="hint" style="margin:2px 0 0">A site is written in only where synonymous codons allow it, so the protein never changes.</div></div>

        <label>Preserve</label><div><div class="fraglist" id="coPres"></div>
          <div class="row2" style="margin-top:4px"><button class="btn" id="coAddSel" type="button">Add current selection</button><button class="btn" id="coAddFeat" type="button">Add a feature…</button></div></div>

        <label for="coRun">Also</label><div class="row2">
          break runs longer than <input type="number" id="coRun" value="6" min="0" max="12"> bases
          <label><input type="checkbox" id="coKeepStart" checked> keep start codon</label>
          <label><input type="checkbox" id="coKeepStop" checked> keep stop codon</label>
          <label>seed <input type="number" id="coSeed" value="1" min="1" max="99999" style="width:70px"></label></div>

        <div class="full" id="coRes"></div></div>`,
      buttons: [
        { label: 'Close' },
        { label: 'Copy', id: 'coCopy', disabled: true, act: () => { if (result) copyText(result.seq, 'Copied optimised sequence'); return false; } },
        { label: 'Open as New Sequence', id: 'coOpen', disabled: true, act: () => {
          if (!result) return false;
          const srcSeq = mode === 'protein' ? '' : regionSeq();
          closeDialog();
          addDoc({
            name: (mode === 'protein' ? 'optimised' : d.name + '-opt'), seq: result.seq, circular: false,
            features: [{ name: 'optimised CDS', type: 'CDS', start: 0, end: result.seq.length, strand: 1 }], primers: [],
            histLabel: `Codon-optimised for ${result.organism.short}`,
            provenance: {
              method: `Codon optimisation for ${result.organism.name}`,
              detail: `${result.algorithm} · ${fmt(result.changedCodons)} of ${fmt(result.totalCodons)} codons changed · CAI ${result.before.cai.toFixed(3)} → ${result.after.cai.toFixed(3)}`,
              inputs: srcSeq ? [{ name: d.name + ' (source region)', seq: srcSeq, circular: false, features: [], detail: `${val('coFrom')}..${val('coTo')}` }] : []
            }
          });
        } },
        { label: 'Replace in Sequence', id: 'coApply', disabled: true, act: () => { applyCodon(); return false; } },
        { label: 'Optimise', primary: true, act: () => { run(); return false; } }
      ],
      onOpen: box => {
        const mkList = (id, chosen, q) => {
          const query = (q || '').toUpperCase();
          box.querySelector('#' + id).innerHTML = enzChoices.filter(x => !query || x.toUpperCase().includes(query))
            .map(x => `<label><input type="checkbox" value="${x}" ${chosen.has(x) ? 'checked' : ''}> ${x}</label>`).join('');
        };
        mkList('coAvoid', avoid, ''); mkList('coAdd', add, '');
        refreshPreserve(); checkORF(); algoNote(); refreshEnzymeHints();
        const orgNote = () => { const t = CODON.tables().find(x => x.key === val('coOrg')); $('#coOrgNote').textContent = t ? `${t.source}, ${fmt(t.codons)} codons counted.` : ''; };
        orgNote();
        box.addEventListener('change', e => {
          const t = e.target;
          if (t.name === 'comode') {
            mode = t.value;
            $('#coRegionBox').style.display = mode === 'region' ? '' : 'none';
            $('#coProteinBox').style.display = mode === 'protein' ? '' : 'none';
            checkORF();
          } else if (t.id === 'coFrom' || t.id === 'coTo' || t.id === 'coRev') { checkORF(); refreshEnzymeHints(); }
          else if (t.id === 'coPick') {
            const v2 = t.value; if (!v2) return;
            if (v2[0] === 'f') { const f = cdsFeatures[+v2.slice(1)]; setRange(f.start, f.end, f.strand === -1); }
            else { const o = orfs[+v2.slice(1)]; setRange(o.start, o.end, o.strand === -1); }
          }
          else if (t.id === 'coAlgo') algoNote();
          else if (t.id === 'coOrg') orgNote();
          else if (t.closest('#coAvoid')) { avoid = new Set([...box.querySelectorAll('#coAvoid input:checked')].map(i => i.value)); }
          else if (t.closest('#coAdd')) { add = new Set([...box.querySelectorAll('#coAdd input:checked')].map(i => i.value)); }
        });
        box.addEventListener('input', e => {
          if (e.target.id === 'coAvoidQ') { avoid = new Set([...box.querySelectorAll('#coAvoid input:checked')].map(i => i.value).concat([...avoid])); mkList('coAvoid', avoid, e.target.value); }
          else if (e.target.id === 'coAddQ') { add = new Set([...box.querySelectorAll('#coAdd input:checked')].map(i => i.value).concat([...add])); mkList('coAdd', add, e.target.value); }
          else if (e.target.id === 'coProt') checkORF();
        });
        box.addEventListener('click', e => {
          const av = e.target.closest('[data-avoid]');
          if (av) { e.preventDefault(); avoid.add(av.dataset.avoid); mkList('coAvoid', avoid, val('coAvoidQ')); toast(`${av.dataset.avoid} added to the sites to remove`); return; }
          const rm = e.target.closest('[data-rmp]');
          if (rm) { preserve.splice(+rm.dataset.rmp, 1); refreshPreserve(); return; }
          if (e.target.id === 'coAddSel') {
            if (!selLen(d)) { toast('Select a region in the sequence first.'); return; }
            preserve.push({ start: d.sel.start, end: d.sel.end || n, name: 'selection' }); refreshPreserve(); return;
          }
          if (e.target.id === 'coAddFeat') {
            const fs = d.features;
            if (!fs.length) { toast('This sequence has no features.'); return; }
            // a small inline chooser avoids opening a second modal on top of this one
            const cur = $('#coPresPick');
            if (cur) { cur.remove(); return; }
            const div = document.createElement('div');
            div.id = 'coPresPick'; div.className = 'fraglist'; div.style.marginTop = '4px';
            div.innerHTML = fs.map((f, i) => `<label><button class="ib" data-pf="${i}">${esc(f.name)}</button> <span style="color:var(--muted)">${fmt(f.start + 1)}–${fmt(f.end)}</span></label>`).join('');
            e.target.parentElement.after(div);
            div.addEventListener('click', ev => {
              const b = ev.target.closest('[data-pf]'); if (!b) return;
              const f = fs[+b.dataset.pf];
              preserve.push({ start: f.start, end: f.end, name: f.name });
              div.remove(); refreshPreserve();
            });
          }
        });
      }
    });

    function applyCodon() {
      if (!result) return;
      if (mode === 'protein') { toast('Use “Open as New Sequence” for a pasted protein.'); return; }
      const a = parseInt(val('coFrom'), 10) - 1, b = parseInt(val('coTo'), 10);
      const rev = $('#coRev').checked;
      const repl = rev ? DNA.rc(result.seq) : result.seq;
      if (b <= a) { toast('A region that wraps the origin cannot be replaced in place.'); return; }
      const dd = D();
      editSeq(dd, a, b - a, repl);
      dd.features.push({ id: uid(), name: `optimised for ${result.organism.short}`, type: 'CDS', start: a, end: a + repl.length, strand: rev ? -1 : 1, color: '#9fd6b4', notes: `Codon-optimised, CAI ${result.after.cai.toFixed(3)}` });
      commit(dd, `Codon-optimised ${fmt(result.changedCodons)} codons for ${result.organism.short}`);
      closeDialog(); renderAll(true);
      toast(`Replaced ${fmt(b - a)} bp with the optimised sequence`);
    }
  }

  // ================================================================ primer design dialogs
  function primerRow(p, label) {
    return `<tr><td><b>${label}</b></td><td class="mono wrap">${esc(p.full || p.seq)}</td><td>${(p.full || p.seq).length}</td><td>${p.tm.toFixed(1)} °C</td><td>${p.gc.toFixed(0)}%</td><td class="wrap">${p.warnings && p.warnings.length ? `<span style="color:var(--danger)">${esc(p.warnings.join('; '))}</span>` : 'no warnings'}</td></tr>`;
  }
  function addDesigned(names) {
    const d = D();
    for (const p of names) d.primers.push({ id: uid(), name: p.name, seq: p.seq, notes: p.notes || '' });
    commit(d, `Added ${names.length} designed primer${names.length > 1 ? 's' : ''}`);
    renderAll(true); toast(`Added ${names.length} primer${names.length > 1 ? 's' : ''}`);
  }
  function dlgDesignPCR() {
    const d = D(), n = d.seq.length;
    const hasSel = selLen(d) > 0;
    let pairs = [];
    const run = () => {
      const a = parseInt(val('pdfrom'), 10) - 1, b = parseInt(val('pdto'), 10);
      if (!(a >= 0 && a < n && b >= 1 && b <= n)) { $('#pdres').innerHTML = '<div class="err">Enter a target range inside the sequence.</div>'; return; }
      try {
        pairs = PRIMER.designPair(d, a, b, {
          optTm: +val('pdtm'), minTm: +val('pdtm') - 3, maxTm: +val('pdtm') + 3,
          minLen: +val('pdmin'), maxLen: +val('pdmax'),
          fwdTail: val('pdftail'), revTail: val('pdrtail'), limit: 5
        });
      } catch (err) { $('#pdres').innerHTML = `<div class="err">${esc(err.message)}</div>`; return; }
      $('#pdres').innerHTML = pairs.map((p, k) => `<div class="pcard"><div class="pch"><b>Pair ${k + 1}</b> — product ${fmt(p.amplicon)} bp (${fmt(p.ampliconStart + 1)}–${fmt(p.ampliconEnd)}), ΔTm ${p.dTm.toFixed(1)} °C${p.pairWarnings.length ? ` · <span style="color:var(--danger)">${esc(p.pairWarnings.join('; '))}</span>` : ''}<button class="btn" data-add="${k}">Add to Primers</button></div>
        <table class="grid"><thead><tr><th></th><th>Sequence 5′→3′</th><th>Length</th><th>Tm</th><th>GC</th><th>Notes</th></tr></thead><tbody>${primerRow(p.fwd, 'Forward')}${primerRow(p.rev, 'Reverse')}</tbody></table></div>`).join('');
    };
    dialog({
      title: 'Design PCR Primers', wide: true,
      body: `<div class="form">
        <label for="pdfrom">Amplify from</label><div class="row2"><input type="number" id="pdfrom" value="${hasSel ? d.sel.start + 1 : 1}" min="1" max="${n}"> to <input type="number" id="pdto" value="${hasSel ? (d.sel.end || n) : Math.min(n, 1000)}" min="1" max="${n}"> (1-based)</div>
        <label for="pdtm">Target Tm</label><div class="row2"><input type="number" id="pdtm" value="60" min="45" max="75" style="width:70px"> °C &nbsp; length <input type="number" id="pdmin" value="18" min="15" max="40" style="width:60px"> – <input type="number" id="pdmax" value="30" min="15" max="45" style="width:60px"> nt</div>
        <label for="pdftail">Forward 5′ tail</label><input type="text" id="pdftail" placeholder="optional, e.g. GCGAATTC for an EcoRI site" spellcheck="false">
        <label for="pdrtail">Reverse 5′ tail</label><input type="text" id="pdrtail" placeholder="optional" spellcheck="false">
        <div class="full hint">Tm uses the nearest-neighbour model at 50 mM Na⁺ and 250 nM primer, and is calculated on the annealing part only, not on any tail. Each primer’s 3′ end is checked for a second match in the template.</div>
        <div class="full" id="pdres"></div></div>`,
      buttons: [{ label: 'Close' }, { label: 'Design', primary: true, act: () => { run(); return false; } }],
      onOpen: box => {
        run();
        box.addEventListener('click', e => {
          const b = e.target.closest('[data-add]'); if (!b) return;
          const p = pairs[+b.dataset.add];
          closeDialog();
          addDesigned([
            { name: `${d.name}-F`, seq: p.fwd.full, notes: `Tm ${p.fwd.tm.toFixed(1)} °C, designed by GeneCanvas` },
            { name: `${d.name}-R`, seq: p.rev.full, notes: `Tm ${p.rev.tm.toFixed(1)} °C, designed by GeneCanvas` }
          ]);
        });
      }
    });
  }
  function dlgDesignWalk() {
    const d = D();
    let list = [];
    const run = () => {
      list = PRIMER.designWalk(d, { spacing: Math.max(150, +val('pwsp') || 500), bothStrands: $('#pwboth').checked });
      $('#pwres').innerHTML = list.length ? `<table class="grid"><thead><tr><th>Name</th><th>Sequence 5′→3′</th><th>Position</th><th>Strand</th><th>Tm</th><th>Notes</th></tr></thead><tbody>${list.map(p => `<tr><td>${esc(p.name)}</td><td class="mono">${esc(p.seq)}</td><td>${fmt(p.start + 1)}</td><td>${p.strand > 0 ? '→' : '←'}</td><td>${p.tm.toFixed(1)} °C</td><td class="wrap">${p.warnings.length ? `<span style="color:var(--danger)">${esc(p.warnings.join('; '))}</span>` : 'no warnings'}</td></tr>`).join('')}</tbody></table>` : '<div class="err">No primers could be placed. Try a larger spacing.</div>';
    };
    dialog({
      title: 'Design Sequencing Primers', wide: true,
      body: `<div class="form">
        <label for="pwsp">Read length</label><div class="row2"><input type="number" id="pwsp" value="500" min="150" max="1200" style="width:80px"> bp of usable sequence per read</div>
        <label></label><label><input type="checkbox" id="pwboth" checked> Also design reverse-strand primers</label>
        <div class="full hint">Primers are spaced so that consecutive reads overlap, with the first 60 bases of each read treated as unreadable.</div>
        <div class="full" id="pwres"></div></div>`,
      buttons: [{ label: 'Close' }, { label: 'Add All to Primers', act: () => { if (list.length) addDesigned(list.map(p => ({ name: p.name, seq: p.seq, notes: `Sequencing primer, Tm ${p.tm.toFixed(1)} °C` }))); } }, { label: 'Design', primary: true, act: () => { run(); return false; } }],
      onOpen: () => run()
    });
  }
  function dlgDesignMut() {
    const d = D(); if (!selLen(d)) { toast('Select the bases to change first.'); return; }
    const a = d.sel.start, b = d.sel.end;
    let r = null;
    const run = () => {
      r = PRIMER.designMutagenesis(d, a, b, val('mtnew'));
      const old = DNA.sub(d.seq, a, b);
      $('#mtres').innerHTML = `<table class="grid"><thead><tr><th></th><th>Sequence 5′→3′</th><th>Length</th><th>Tm</th><th>GC</th></tr></thead><tbody>
        <tr><td><b>Forward</b></td><td class="mono wrap">${esc(r.fwd.seq)}</td><td>${r.fwd.len}</td><td>${r.fwd.tm.toFixed(1)} °C</td><td>${r.fwd.gc.toFixed(0)}%</td></tr>
        <tr><td><b>Reverse</b></td><td class="mono wrap">${esc(r.rev.seq)}</td><td>${r.rev.len}</td><td>${r.rev.tm.toFixed(1)} °C</td><td>${r.rev.gc.toFixed(0)}%</td></tr>
        </tbody></table><div class="hint">Replaces ${fmt(r.removed)} bp (${esc(old.slice(0, 30))}${old.length > 30 ? '…' : ''}) with ${r.insert.length ? esc(r.insert) : 'nothing'}, with ${r.armLen} bp of matching sequence on each side.</div>`;
    };
    dialog({
      title: 'Design Mutagenesis Primers', wide: true,
      body: `<div class="form">
        <label>Replace</label><div>${fmt(a + 1)}–${fmt(b)} (${fmt(selLen(d))} bp)</div>
        <label for="mtnew">With</label><input type="text" id="mtnew" placeholder="new bases, or leave empty to delete" spellcheck="false">
        <div class="full" id="mtres"></div></div>`,
      buttons: [{ label: 'Close' }, { label: 'Add to Primers', act: () => { if (r) addDesigned([{ name: 'Mut-F', seq: r.fwd.seq, notes: 'Mutagenesis primer' }, { name: 'Mut-R', seq: r.rev.seq, notes: 'Mutagenesis primer' }]); } }, { label: 'Design', primary: true, act: () => { run(); return false; } }],
      onOpen: box => { run(); $('#mtnew').addEventListener('input', run); }
    });
  }

  // ================================================================ agarose gel document
  const isGel = x => !!(x && x.gel);
  function newGelDoc(opts) {
    const dnaDocs = S.docs.filter(x => !isGel(x));
    const first = (opts && opts.doc) || dnaDocs[0];
    const g = {
      id: uid(), name: 'Agarose Gel', gel: {
        percent: 1.0, buffer: 'TAE', run: 0.75, laneCount: 6, ladder: 'neb-1kb',
        lanes: []
      }, dirty: false
    };
    g.gel.lanes = [{ src: 'ladder:neb-1kb' }];
    if (first) g.gel.lanes.push({ src: 'doc:' + first.id, enzymes: (opts && opts.enzymes) || first.chosen.slice(0, 2) });
    while (g.gel.lanes.length < g.gel.laneCount) g.gel.lanes.push({ src: '' });
    S.docs.push(g);
    S.active = S.docs.length - 1;
    renderAll();
    return g;
  }

  function laneBands(lane) {
    if (!lane || !lane.src) return { bands: [], label: '', sub: '' };
    if (lane.src.startsWith('ladder:')) {
      const l = GEL.ladderByKey(lane.src.slice(7), SET.ladders);
      if (!l) return { bands: [], label: 'missing ladder', sub: '' };
      const bright = new Set(l.bright || []);
      return { bands: l.bands.map(size => ({ size, mass: bright.has(size) ? size * 2.4 : size, ref: bright.has(size) })), label: l.name, sub: l.brand, ladder: true };
    }
    if (lane.src.startsWith('doc:')) {
      const doc = S.docs.find(x => String(x.id) === lane.src.slice(4) && !isGel(x));
      if (!doc) return { bands: [], label: 'missing sequence', sub: '' };
      const enz = (lane.enzymes || []).filter(nm => DNA.ENZ_BY_NAME[nm]);
      if (!enz.length) {
        const n = doc.seq.length;
        // A supercoiled circle is more compact than a linear fragment and runs ahead of it.
        return { bands: [{ size: doc.circular ? Math.round(n * 0.65) : n, mass: n, real: n, sc: doc.circular }], label: doc.name, sub: doc.circular ? 'uncut, supercoiled' : 'uncut' };
      }
      const sites = DNA.findSites(doc.seq, doc.circular);
      const frags = DNA.digest(doc.seq, doc.circular, sites, enz);
      return { bands: frags.map(f => ({ size: f.size, mass: f.size })), label: doc.name, sub: enz.join(' + ') };
    }
    if (lane.src === 'custom') {
      const list = (lane.bands || []);
      return { bands: list.map(size => ({ size, mass: size })), label: lane.label || 'Custom', sub: list.length + ' bands' };
    }
    return { bands: [], label: '', sub: '' };
  }

  function gelSVGDoc(g, width) {
    const o = g.gel;
    const lanes = o.lanes.slice(0, o.laneCount);
    const avail = Math.max(360, width - 120);
    const laneW = Math.max(54, Math.min(170, Math.floor(avail / Math.max(1, lanes.length))));
    const left = 86, top = 58, gelH = 520, wellY = top + 10, runH = gelH - 54;
    const W = left + lanes.length * laneW + 26, H = top + gelH + 34;
    const data = lanes.map(laneBands);
    const parts = [];
    parts.push(`<rect width="${W}" height="${H}" fill="#0f1116" rx="5"/>`);
    parts.push(`<rect x="${left - 10}" y="${top}" width="${lanes.length * laneW + 20}" height="${gelH}" fill="#15181f" rx="3"/>`);
    // ladder size labels down the left, taken from the first ladder lane
    const refIdx = data.findIndex(x => x.ladder);
    if (refIdx >= 0) {
      // Bands that pile up near the well cannot all be labelled, so skip the ones that
      // would print on top of the label above them.
      let lastY = -99;
      for (const b of data[refIdx].bands.slice().sort((p, q) => q.size - p.size)) {
        const y = wellY + GEL.mobility(b.size, o) * runH;
        if (y > top + gelH - 6) continue;
        if (y - lastY < 11) continue;
        lastY = y;
        parts.push(`<text x="${left - 16}" y="${(y + 3.5).toFixed(1)}" fill="#9aa3b2" font-size="10" text-anchor="end">${b.size >= 1000 ? (b.size / 1000) + ' kb' : b.size + ' bp'}</text>`);
        parts.push(`<line x1="${left - 13}" y1="${y.toFixed(1)}" x2="${left - 9}" y2="${y.toFixed(1)}" stroke="#4a5262"/>`);
      }
    }
    lanes.forEach((lane, i) => {
      const x = left + i * laneW, cx = x + laneW / 2, bw = Math.min(laneW - 12, 120);
      parts.push(`<rect x="${(cx - bw / 2).toFixed(1)}" y="${wellY - 9}" width="${bw}" height="7" fill="#2b3040" rx="1"/>`);
      const d2 = data[i];
      const head = d2.label || '—';
      parts.push(`<text x="${cx.toFixed(1)}" y="26" fill="#d6dbe5" font-size="11" text-anchor="middle" font-weight="600">${esc(head.length > Math.floor(laneW / 6) ? head.slice(0, Math.floor(laneW / 6) - 1) + '…' : head)}</text>`);
      if (d2.sub) parts.push(`<text x="${cx.toFixed(1)}" y="40" fill="#8e97a8" font-size="9.5" text-anchor="middle">${esc(d2.sub.length > Math.floor(laneW / 5) ? d2.sub.slice(0, Math.floor(laneW / 5) - 1) + '…' : d2.sub)}</text>`);
      parts.push(`<text x="${cx.toFixed(1)}" y="${H - 10}" fill="#6f7889" font-size="10" text-anchor="middle">${i + 1}</text>`);
      if (!d2.bands.length) return;
      // bands close enough to overlap are drawn as one, which is what a real gel shows
      const placed = [];
      const maxMass = Math.max(...d2.bands.map(b => b.mass), 1);
      for (const b of d2.bands.slice().sort((p, q) => q.size - p.size)) {
        const y = wellY + GEL.mobility(b.size, o) * runH;
        if (y > top + gelH - 4) continue;
        const near = placed.find(p => Math.abs(p.y - y) < 2.2);
        if (near) { near.mass += b.mass; near.sizes.push(b.size); continue; }
        placed.push({ y, mass: b.mass, sizes: [b.size], ref: b.ref, sc: b.sc, real: b.real });
      }
      for (const p of placed) {
        const op = 0.3 + 0.7 * Math.min(1, p.mass / maxMass);
        const h = 3 + 1.6 * Math.min(1, p.mass / maxMass);
        const title = p.sizes.map(s => fmt(s) + ' bp').join(' + ') + (p.sc ? ` (supercoiled, really ${fmt(p.real)} bp)` : '');
        parts.push(`<rect x="${(cx - bw / 2).toFixed(1)}" y="${(p.y - h / 2).toFixed(1)}" width="${bw}" height="${h.toFixed(1)}" fill="#f2f6ff" opacity="${op.toFixed(2)}" rx="1"><title>${esc(title)}</title></rect>`);
      }
    });
    return `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="Segoe UI, Arial, sans-serif">${parts.join('')}</svg>`;
  }

  function renderGel() {
    const g = D(), o = g.gel, v = $('#view');
    const dnaDocs = S.docs.filter(x => !isGel(x));
    const ladders = GEL.ladderList(SET.ladders);
    const byBrand = {};
    for (const l of ladders) (byBrand[l.brand] = byBrand[l.brand] || []).push(l);
    const srcOptions = lane => {
      const cur = lane.src || '';
      let h = `<option value="" ${cur === '' ? 'selected' : ''}>— empty —</option>`;
      h += Object.entries(byBrand).map(([brand, list]) => `<optgroup label="${esc(brand)}">${list.map(l => `<option value="ladder:${esc(l.key)}" ${cur === 'ladder:' + l.key ? 'selected' : ''}>${esc(l.name)}</option>`).join('')}</optgroup>`).join('');
      if (dnaDocs.length) h += `<optgroup label="Open sequences">${dnaDocs.map(x => `<option value="doc:${x.id}" ${cur === 'doc:' + x.id ? 'selected' : ''}>${esc(x.name)}</option>`).join('')}</optgroup>`;
      return h;
    };
    const lanes = o.lanes.slice(0, o.laneCount);
    v.innerHTML = `<div class="panel">
      <div class="ptools">
        <label>Agarose <select id="glPct">${[0.5, 0.7, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0].map(p => `<option value="${p}" ${p === o.percent ? 'selected' : ''}>${p.toFixed(1)} %</option>`).join('')}</select></label>
        <label>Buffer <select id="glBuf">${['TAE', 'TBE'].map(x => `<option ${x === o.buffer ? 'selected' : ''}>${x}</option>`).join('')}</select></label>
        <label>Run <input type="range" id="glRun" min="20" max="100" step="5" value="${Math.round(o.run * 100)}" style="width:110px"> <span id="glRunV">${Math.round(o.run * 100)}%</span></label>
        <label>Lanes <select id="glN">${Array.from({ length: 20 }, (_, i) => i + 1).map(k => `<option value="${k}" ${k === o.laneCount ? 'selected' : ''}>${k}</option>`).join('')}</select></label>
        <button class="btn" id="glLadder">Add a ladder…</button>
        <button class="btn" id="glPNG">Export SVG…</button>
        <span class="note">${(() => { const r = GEL.rangeFor(o.percent); return `resolves about ${fmt(Math.round(r[1]))}–${fmt(Math.round(r[0]))} bp`; })()}</span>
      </div>
      <div class="tscroll" style="background:var(--chrome2)"><div id="glImg" style="padding:14px"></div>
        <div style="padding:0 14px 14px">
          <table class="grid" style="width:auto;background:var(--win)"><thead><tr><th>Lane</th><th>Contents</th><th>Cut with</th><th>Bands</th></tr></thead><tbody>
          ${lanes.map((lane, i) => {
            const d2 = laneBands(lane);
            const isDoc = (lane.src || '').startsWith('doc:');
            return `<tr><td>${i + 1}</td>
              <td><select data-lane="${i}" class="glSrc">${srcOptions(lane)}</select></td>
              <td>${isDoc ? `<button class="btn" data-enz="${i}">${(lane.enzymes || []).length ? esc((lane.enzymes || []).join(' + ')) : 'uncut — choose…'}</button>` : '<span class="note">—</span>'}</td>
              <td class="wrap">${d2.bands.length ? d2.bands.slice().sort((a, b) => b.size - a.size).map(b => fmt(b.size)).join(', ') + ' bp' : '<span class="note">—</span>'}</td></tr>`;
          }).join('')}
          </tbody></table>
        </div>
      </div>
      <div class="ptools"><span class="note">Migration is a model of a 1× ${esc(o.buffer)} gel, not a measurement. Ladder band sizes follow the suppliers' published values — check them against the datasheet in your box.</span></div>
    </div>`;
    const img = $('#glImg');
    img.innerHTML = gelSVGDoc(g, img.clientWidth || 900);
    const mark = () => { g.dirty = true; renderAll(true); };
    $('#glPct').onchange = e => { o.percent = parseFloat(e.target.value); mark(); };
    $('#glBuf').onchange = e => { o.buffer = e.target.value; mark(); };
    $('#glRun').oninput = e => { o.run = +e.target.value / 100; $('#glRunV').textContent = e.target.value + '%'; img.innerHTML = gelSVGDoc(g, img.clientWidth || 900); };
    $('#glRun').onchange = () => { g.dirty = true; };
    $('#glN').onchange = e => {
      o.laneCount = +e.target.value;
      while (o.lanes.length < o.laneCount) o.lanes.push({ src: '' });
      mark();
    };
    $('#glLadder').onclick = dlgCustomLadder;
    $('#glPNG').onclick = () => {
      const svg = img.querySelector('svg');
      const content = svg.outerHTML;
      if (host && host.saveFile) host.saveFile({ path: null, defaultName: 'gel.svg', content, filters: [{ name: 'SVG image', extensions: ['svg'] }] }).then(p => { if (p) toast('Exported ' + p.split(/[\\/]/).pop()); });
      else downloadText('gel.svg', content);
    };
    v.querySelectorAll('.glSrc').forEach(sel => sel.onchange = e => {
      const i = +e.target.dataset.lane;
      o.lanes[i] = { src: e.target.value };
      if (e.target.value.startsWith('doc:')) {
        const doc = S.docs.find(x => String(x.id) === e.target.value.slice(4));
        o.lanes[i].enzymes = doc ? doc.chosen.slice(0, 2) : [];
      }
      mark();
    });
    v.querySelectorAll('[data-enz]').forEach(btn => btn.onclick = () => {
      const i = +btn.dataset.enz;
      const doc = S.docs.find(x => String(x.id) === o.lanes[i].src.slice(4));
      if (!doc) return;
      dlgEnzymeChooser({ doc, initial: o.lanes[i].enzymes || [], title: `Cut ${doc.name} with`, onApply: names => { o.lanes[i].enzymes = names; mark(); } });
    });
  }

  function dlgCustomLadder() {
    dialog({
      title: 'Add a Ladder',
      body: `<div class="form">
        <label for="clN">Name</label><input type="text" id="clN" placeholder="e.g. Bioneer 1 kb">
        <label for="clB">Band sizes</label><textarea id="clB" rows="5" placeholder="10000, 8000, 6000, 5000, 4000, 3000*, 2000, 1500, 1000, 500"></textarea>
        <div class="full hint">Sizes in base pairs, separated by commas, spaces or line breaks. Write <b>3000*</b> to mark a brighter reference band. A number under 50 is read as kilobases.</div>
        <div class="full" id="clPrev"></div><div class="full err"></div></div>`,
      buttons: [{ label: 'Cancel' }, {
        label: 'Save', primary: true, act: () => {
          const name = val('clN').trim();
          if (!name) return showErr('Give the ladder a name.');
          let parsed;
          try { parsed = GEL.parseLadder(val('clB')); } catch (e) { return showErr(e.message); }
          SET.ladders[name] = parsed;
          saveSettings().then(() => { toast(`Saved the ladder “${name}” with ${parsed.bands.length} bands`); renderAll(true); });
        }
      }],
      onOpen: box => {
        const upd = () => {
          try { const r = GEL.parseLadder(val('clB')); $('#clPrev').innerHTML = `<div class="hint">${r.bands.length} bands: ${r.bands.map(b => fmt(b)).join(', ')}${r.bright.length ? ` · reference: ${r.bright.map(fmt).join(', ')}` : ''}</div>`; }
          catch (e) { $('#clPrev').innerHTML = ''; }
        };
        box.querySelector('#clB').addEventListener('input', upd);
      }
    });
  }

  // ================================================================ stored settings
  // Preferences shared by every file: saved enzyme sets, map and sequence view defaults.
  const SETTINGS_KEY = 'genecanvas.settings';
  const DEFAULT_SETTINGS = { enzymeSets: {}, lastEnzymeSet: '', mapView: null, seqView: null, ladders: {} };
  let SET = { ...DEFAULT_SETTINGS };

  async function loadSettings() {
    let raw = null;
    try {
      if (host && host.loadSettings) raw = await host.loadSettings();
      else { const t = localStorage.getItem(SETTINGS_KEY); raw = t ? JSON.parse(t) : null; }
    } catch (e) { raw = null; }
    SET = { ...DEFAULT_SETTINGS, ...(raw || {}) };
    SET.enzymeSets = SET.enzymeSets || {};
    SET.ladders = SET.ladders || {};
    return SET;
  }
  async function saveSettings() {
    try {
      if (host && host.saveSettings) await host.saveSettings(SET);
      else localStorage.setItem(SETTINGS_KEY, JSON.stringify(SET));
      return true;
    } catch (e) { toast('Could not save settings: ' + e.message); return false; }
  }

  // Sets that are always offered, computed from the sequence rather than stored.
  function builtinSets(d, A) {
    const cutters = DNA.ENZYMES.filter(e => (A.sites[e.name] || []).some(s => s.cut != null));
    return {
      'Unique cutters, 6+ bases': cutters.filter(e => A.sites[e.name].length === 1 && e.core >= 6).map(e => e.name),
      'Unique cutters': cutters.filter(e => A.sites[e.name].length === 1).map(e => e.name),
      'Cuts once or twice, 6+ bases': cutters.filter(e => A.sites[e.name].length <= 2 && e.core >= 6).map(e => e.name),
      'All enzymes that cut': cutters.map(e => e.name),
      'Type IIS (BsaI, BsmBI, BbsI, SapI)': CLONE.TYPE2S.slice(),
      'Blunt cutters': cutters.filter(e => e.overhang === 0).map(e => e.name)
    };
  }

  // ================================================================ enzyme chooser
  function dlgEnzymeChooser(opts) {
    opts = typeof opts === 'function' ? { onApply: opts } : (opts || {});
    const d = opts.doc || D(), A = analyze(d);
    const chosen = new Set(opts.initial || d.chosen);
    const built = builtinSets(d, A);
    let filter = '', show = 'cutters', sort = 'name';
    const countOf = nm => (A.sites[nm] || []).length;

    const listFor = () => {
      let list = DNA.ENZYMES.slice();
      if (show === 'cutters') list = list.filter(e => countOf(e.name) > 0);
      else if (show === 'unique') list = list.filter(e => countOf(e.name) === 1);
      else if (show === 'noncut') list = list.filter(e => countOf(e.name) === 0);
      else if (show === 'chosen') list = list.filter(e => chosen.has(e.name));
      else if (show === 'six') list = list.filter(e => e.core >= 6);
      if (filter) {
        const q = filter.toUpperCase();
        list = list.filter(e => e.name.toUpperCase().includes(q) || e.site.includes(q));
      }
      if (sort === 'count') list.sort((a, b) => countOf(a.name) - countOf(b.name) || a.name.localeCompare(b.name));
      else if (sort === 'site') list.sort((a, b) => a.site.localeCompare(b.site) || a.name.localeCompare(b.name));
      else list.sort((a, b) => a.name.localeCompare(b.name));
      return list;
    };
    const renderGrid = () => {
      const list = listFor();
      $('#ecGrid').innerHTML = list.length ? list.map(e => {
        const n = countOf(e.name);
        return `<label class="ecitem${n ? '' : ' nocut'}" title="${esc(e.name)} ${esc(e.site)} — ${n ? n + ' site' + (n > 1 ? 's' : '') : 'does not cut'}">
          <input type="checkbox" value="${esc(e.name)}" ${chosen.has(e.name) ? 'checked' : ''}>
          <span class="ecn">${esc(e.name)}</span><span class="ecc">${n || ''}</span></label>`;
      }).join('') : '<div class="hint">Nothing matches that filter.</div>';
      $('#ecCount').textContent = `${list.length} shown · ${chosen.size} chosen`;
    };
    const renderSets = () => {
      const saved = Object.keys(SET.enzymeSets).sort();
      $('#ecSet').innerHTML = `<option value="">— choose a set —</option>`
        + `<optgroup label="From this sequence">${Object.keys(built).map(k => `<option value="b:${esc(k)}">${esc(k)} (${built[k].length})</option>`).join('')}</optgroup>`
        + (saved.length ? `<optgroup label="My sets">${saved.map(k => `<option value="s:${esc(k)}">${esc(k)} (${SET.enzymeSets[k].length})</option>`).join('')}</optgroup>` : '');
      $('#ecDel').disabled = !saved.length;
    };
    dialog({
      title: opts.title || 'Choose Enzymes', wide: true,
      body: `<div class="ecbar">
          <input type="text" id="ecQ" placeholder="Filter by name or site" style="width:160px">
          <select id="ecShow">
            <option value="cutters">Enzymes that cut</option>
            <option value="unique">Unique cutters</option>
            <option value="six">6+ base recognition</option>
            <option value="chosen">Currently chosen</option>
            <option value="noncut">Non-cutters</option>
            <option value="all">All ${DNA.ENZYMES.length} enzymes</option>
          </select>
          <select id="ecSort"><option value="name">Sort by name</option><option value="count">Sort by site count</option><option value="site">Sort by site</option></select>
          <span class="sp" style="flex:1"></span>
          <button class="btn" id="ecAll">Check shown</button><button class="btn" id="ecNone">Uncheck shown</button>
        </div>
        <div class="ecbar">
          <span>Set</span><select id="ecSet" style="min-width:220px"></select>
          <button class="btn" id="ecLoad">Load</button><button class="btn" id="ecAdd">Add to selection</button>
          <span class="sp" style="flex:1"></span>
          <button class="btn" id="ecSave">Save as my set…</button><button class="btn dng" id="ecDel">Delete my set</button>
        </div>
        <div class="ecgrid" id="ecGrid"></div>
        <div class="hint" id="ecCount"></div>
        <div class="hint">A saved set is stored with the program, not with the file, so it is there for every sequence you open.</div>`,
      buttons: [{ label: 'Cancel' }, {
        label: 'Apply', primary: true, act: () => {
          const names = [...chosen].sort();
          if (opts.onApply) { opts.onApply(names); return; }
          d.chosen = names;
          d.opts.enzMode = 'chosen';
          renderAll(true);
          toast(`${names.length} enzyme${names.length === 1 ? '' : 's'} chosen`);
        }
      }],
      onOpen: box => {
        renderSets(); renderGrid();
        box.querySelector('#ecQ').addEventListener('input', e => { filter = e.target.value; renderGrid(); });
        box.querySelector('#ecShow').addEventListener('change', e => { show = e.target.value; renderGrid(); });
        box.querySelector('#ecSort').addEventListener('change', e => { sort = e.target.value; renderGrid(); });
        box.querySelector('#ecGrid').addEventListener('change', e => {
          const i = e.target.closest('input'); if (!i) return;
          if (i.checked) chosen.add(i.value); else chosen.delete(i.value);
          $('#ecCount').textContent = `${listFor().length} shown · ${chosen.size} chosen`;
        });
        box.addEventListener('click', e => {
          const id = e.target.id;
          if (id === 'ecAll') { for (const el of box.querySelectorAll('#ecGrid input')) { el.checked = true; chosen.add(el.value); } renderGrid(); }
          else if (id === 'ecNone') { for (const el of box.querySelectorAll('#ecGrid input')) { el.checked = false; chosen.delete(el.value); } renderGrid(); }
          else if (id === 'ecLoad' || id === 'ecAdd') {
            const v = val('ecSet'); if (!v) { toast('Choose a set first.'); return; }
            const names = v.startsWith('b:') ? built[v.slice(2)] : SET.enzymeSets[v.slice(2)];
            if (!names) return;
            if (id === 'ecLoad') chosen.clear();
            for (const nm of names) chosen.add(nm);
            renderGrid();
          } else if (id === 'ecSave') {
            if (!chosen.size) { toast('Check some enzymes first.'); return; }
            promptName('Save enzyme set', 'Name', '', async name => {
              SET.enzymeSets[name] = [...chosen].sort();
              SET.lastEnzymeSet = name;
              await saveSettings();
              renderSets();
              $('#ecSet').value = 's:' + name;
              toast(`Saved the set “${name}” with ${chosen.size} enzymes`);
            });
          } else if (id === 'ecDel') {
            const v = val('ecSet');
            if (!v.startsWith('s:')) { toast('Choose one of your own sets to delete.'); return; }
            const name = v.slice(2);
            delete SET.enzymeSets[name];
            saveSettings().then(() => { renderSets(); toast(`Deleted the set “${name}”`); });
          }
        });
      }
    });
  }
  // A one-field prompt that can sit on top of another dialog.
  function promptName(title, label, initial, done) {
    const bg = document.createElement('div');
    bg.className = 'modal-bg';
    bg.style.zIndex = 120;
    bg.innerHTML = `<div class="dlg" style="width:360px"><div class="hd">${esc(title)}</div>
      <div class="bd"><div class="form"><label for="pnV">${esc(label)}</label><input type="text" id="pnV" value="${esc(initial)}"><div class="full err" id="pnE"></div></div></div>
      <div class="ft"><button class="btn" data-x="1">Cancel</button><button class="btn pri" data-ok="1">OK</button></div></div>`;
    document.body.appendChild(bg);
    const input = bg.querySelector('#pnV');
    const close = () => bg.remove();
    const ok = () => {
      const v = input.value.trim();
      if (!v) { bg.querySelector('#pnE').textContent = 'Enter a name.'; return; }
      close(); done(v);
    };
    bg.addEventListener('click', e => { if (e.target.dataset.x) close(); else if (e.target.dataset.ok) ok(); });
    bg.addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); ok(); } else if (e.key === 'Escape') { e.preventDefault(); close(); } });
    setTimeout(() => { input.focus(); input.select(); }, 0);
  }

  // ================================================================ databases
  function ends(e) { return e.overhang > 0 ? `5′ overhang (${e.overhang})` : (e.overhang < 0 ? `3′ overhang (${-e.overhang})` : 'blunt'); }
  function cutSiteHtml(e) {
    const s = e.site;
    if (e.top >= 0 && e.top <= s.length) return esc(s.slice(0, e.top)) + '<i>▼</i>' + esc(s.slice(e.top));
    const off = e.top - s.length;
    return esc(s) + `<span class="note">(${off}/${e.bot - s.length})</span>`;
  }
  function dlgEnzymeDB() {
    const d = D();
    const A = d ? analyze(d) : { sites: {} };
    let q = '', only = 'all';
    const render = () => {
      let list = DNA.ENZYMES.slice();
      if (only === 'blunt') list = list.filter(e => e.overhang === 0);
      else if (only === '5') list = list.filter(e => e.overhang > 0);
      else if (only === '3') list = list.filter(e => e.overhang < 0);
      else if (only === 'iis') list = list.filter(e => e.top > e.site.length);
      else if (only === 'six') list = list.filter(e => e.core >= 6);
      if (q) { const u = q.toUpperCase(); list = list.filter(e => e.name.toUpperCase().includes(u) || e.site.includes(u)); }
      list.sort((a, b) => a.name.localeCompare(b.name));
      $('#edbT').innerHTML = `<table class="grid"><thead><tr><th>Enzyme</th><th>Recognition site</th><th>Length</th><th>Ends</th><th>Cuts (top/bottom)</th>${d ? `<th>In ${esc(d.name)}</th>` : ''}</tr></thead><tbody>
        ${list.map(e => `<tr><td><b>${esc(e.name)}</b></td><td class="cut-site">${cutSiteHtml(e)}</td><td>${e.core}</td><td>${ends(e)}</td><td>${e.top} / ${e.bot}</td>${d ? `<td>${(A.sites[e.name] || []).length || '—'}</td>` : ''}</tr>`).join('')}
      </tbody></table>`;
      $('#edbN').textContent = `${list.length} of ${DNA.ENZYMES.length} enzymes`;
    };
    dialog({
      title: 'Enzyme Database', wide: true,
      body: `<div class="ecbar"><input type="text" id="edbQ" placeholder="Filter by name or site" style="width:180px">
        <select id="edbF"><option value="all">All enzymes</option><option value="six">6+ base recognition</option><option value="5">5′ overhang</option><option value="3">3′ overhang</option><option value="blunt">Blunt</option><option value="iis">Type IIS (cuts outside the site)</option></select>
        <span class="sp" style="flex:1"></span><span class="note" id="edbN"></span></div>
        <div class="tscroll" style="max-height:52vh" id="edbT"></div>
        <div class="hint">▼ marks where the top strand is cut. For a Type IIS enzyme the cut lies outside the recognition site, so the two offsets are given instead. Cut positions are counted from the left end of the site on the top strand.</div>`,
      buttons: [{ label: 'Export as CSV', act: () => { exportEnzymeCSV(); return false; } }, { label: 'Close', primary: true }],
      onOpen: box => {
        render();
        box.querySelector('#edbQ').addEventListener('input', e => { q = e.target.value; render(); });
        box.querySelector('#edbF').addEventListener('change', e => { only = e.target.value; render(); });
      }
    });
  }
  function exportEnzymeCSV() {
    const rows = [['Enzyme', 'Site', 'CoreLength', 'TopCut', 'BottomCut', 'Ends']];
    for (const e of DNA.ENZYMES) rows.push([e.name, e.site, e.core, e.top, e.bot, ends(e).replace(/[′]/g, "'")]);
    const csv = rows.map(r => r.join(',')).join('\n');
    if (host && host.saveFile) host.saveFile({ path: null, defaultName: 'genecanvas-enzymes.csv', content: csv, filters: [{ name: 'CSV', extensions: ['csv'] }] }).then(p => { if (p) toast('Exported ' + p.split(/[\\/]/).pop()); });
    else downloadText('genecanvas-enzymes.csv', csv);
  }
  function dlgFeatureDB() {
    let q = '';
    const render = () => {
      let list = DNA.COMMON.slice();
      if (q) { const u = q.toUpperCase(); list = list.filter(f => f.name.toUpperCase().includes(u) || f.type.toUpperCase().includes(u) || f.seq.includes(u)); }
      $('#fdbT').innerHTML = `<table class="grid"><thead><tr><th>Name</th><th>Type</th><th>Matched on</th><th>Length</th><th>Sequence</th></tr></thead><tbody>
        ${list.map(f => `<tr><td><b>${esc(f.name)}</b></td><td>${esc(f.type)}</td><td>${f.kind === 'protein' ? 'protein' : 'DNA'}${f.short ? ' <span class="note">(short — off by default)</span>' : ''}</td><td>${f.seq.length} ${f.kind === 'protein' ? 'aa' : 'bp'}</td><td class="mono wrap">${esc(f.seq)}</td></tr>`).join('')}
      </tbody></table>`;
      $('#fdbN').textContent = `${list.length} of ${DNA.COMMON.length} features`;
    };
    dialog({
      title: 'Feature Database', wide: true,
      body: `<div class="ecbar"><input type="text" id="fdbQ" placeholder="Filter by name, type or sequence" style="width:220px"><span class="sp" style="flex:1"></span><span class="note" id="fdbN"></span></div>
        <div class="tscroll" style="max-height:52vh" id="fdbT"></div>
        <div class="hint">Features matched <b>on protein</b> are found whatever codons were used, which is what makes a tag turn up in a gene that was codon-optimised. Entries marked short are skipped by Detect Common Features unless you ask for them, because a motif of a few residues appears by chance.</div>`,
      buttons: [{ label: 'Export as CSV', act: () => { exportFeatureCSV(); return false; } }, { label: 'Close', primary: true }],
      onOpen: box => { render(); box.querySelector('#fdbQ').addEventListener('input', e => { q = e.target.value; render(); }); }
    });
  }
  function exportFeatureCSV() {
    const rows = [['Name', 'Type', 'MatchedOn', 'Length', 'Sequence']];
    for (const f of DNA.COMMON) rows.push([JSON.stringify(f.name), f.type, f.kind, f.seq.length, f.seq]);
    const csv = rows.map(r => r.join(',')).join('\n');
    if (host && host.saveFile) host.saveFile({ path: null, defaultName: 'genecanvas-features.csv', content: csv, filters: [{ name: 'CSV', extensions: ['csv'] }] }).then(p => { if (p) toast('Exported ' + p.split(/[\\/]/).pop()); });
    else downloadText('genecanvas-features.csv', csv);
  }

  // ================================================================ structure analysis
  // Links open in the default browser; the Electron main process handles that.
  function extLink(url, label) { return `<a href="${esc(url)}" target="_blank" rel="noreferrer">${esc(label)} ↗</a>`; }
  function toolLinks(kind, seq) {
    const rows = kind === 'rna' ? [
      ['ViennaRNA RNAfold', 'http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi', 'Minimum free energy and base-pair probabilities with the full Turner parameters.'],
      ['mfold / UNAFold', 'http://www.unafold.org/mfold/applications/rna-folding-form.php', 'The classic Zuker server, with suboptimal structures.'],
      ['RNAcentral', 'https://rnacentral.org/sequence-search', 'Find out whether the sequence is a known non-coding RNA.']
    ] : [
      ['AlphaFold DB', 'https://alphafold.ebi.ac.uk/', 'Predicted 3D structures for most known proteins. Search by name or UniProt accession.'],
      ['ESMFold (ESM Atlas)', 'https://esmatlas.com/resources?action=fold', 'Folds a pasted sequence directly, with no database entry needed.'],
      ['SWISS-MODEL', 'https://swissmodel.expasy.org/interactive', 'Homology modelling against solved structures.'],
      ['InterPro', 'https://www.ebi.ac.uk/interpro/search/sequence/', 'Domains, families and functional sites.'],
      ['NCBI BLASTp', 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch', 'Find related sequences.']
    ];
    return `<div class="hint" style="margin-top:8px"><b>Open in an external tool.</b> These sites take a pasted sequence, so copy it first.</div>
      <div class="row2" style="margin:4px 0"><button class="btn" id="copySeqBtn">Copy sequence</button><span class="note" id="copySeqNote"></span></div>
      <table class="grid"><tbody>${rows.map(([n2, u, why]) => `<tr><td style="white-space:nowrap">${extLink(u, n2)}</td><td class="wrap">${esc(why)}</td></tr>`).join('')}</tbody></table>`;
  }
  function bindCopySeq(seq, what) {
    const b = $('#copySeqBtn');
    if (b) b.onclick = () => { copyText(seq, `Copied ${what}`); const n2 = $('#copySeqNote'); if (n2) n2.textContent = 'Copied. Paste it into the tool you open.'; };
  }

  // ---------------------------------------------------------------- primer structures
  function dgClass(dg) { return dg <= -9 ? 'bad' : (dg <= -6 ? 'warn' : 'okdg'); }
  function dgNote(dg) { return dg <= -9 ? 'strong — likely to interfere' : (dg <= -6 ? 'moderate' : 'weak'); }
  function dlgPrimerStruct(primer, partner) {
    const d = D();
    const list = d.primers;
    let aId = primer ? primer.id : (list[0] && list[0].id);
    let bId = partner ? partner.id : (list.find(p => p.id !== aId) || {}).id || '';
    const render = () => {
      const a = list.find(p => p.id === aId);
      if (!a) { $('#psRes').innerHTML = '<div class="hint">This sequence has no primers yet.</div>'; return; }
      const b = bId ? list.find(p => p.id === bId) : null;
      const r = STRUCT.primerStructures(a.seq, b ? b.seq : null);
      const dupBlock = (title, seqA, seqB, dups, note) => {
        if (!dups.length) return `<h4 class="psh">${esc(title)}</h4><div class="hint">None found.</div>`;
        return `<h4 class="psh">${esc(title)}</h4>` + dups.map(x => {
          const dg = STRUCT.duplexDiagram(seqA, seqB, x);
          return `<div class="psbox"><div class="pshead"><span class="dg ${dgClass(x.dg)}">ΔG ${x.dg.toFixed(2)} kcal/mol</span> <span class="note">${x.pairs} pairs, longest run ${x.longest}${x.end3a || x.end3b ? ` · ${Math.max(x.end3a, x.end3b)} of the last 5 bases paired` : ''} · ${dgNote(x.dg)}</span></div><pre class="psdiag">${esc(dg.top)}\n${esc(dg.bars)}\n${esc(dg.bot)}</pre></div>`;
        }).join('') + (note ? `<div class="hint">${esc(note)}</div>` : '');
      };
      const hpBlock = (title, seq, hps) => {
        if (!hps.length) return `<h4 class="psh">${esc(title)}</h4><div class="hint">None found.</div>`;
        return `<h4 class="psh">${esc(title)}</h4>` + hps.map(h => {
          const g = STRUCT.hairpinDiagram(seq, h);
          // The stem is drawn as a ladder and the loop hangs off its closed end.
          const postRev = [...g.post].reverse().join('');
          const lead = Math.max(g.pre.length, postRev.length);
          const l1 = "5'-" + g.pre.padStart(lead) + g.left + '\\';
          const l2 = '   ' + ' '.repeat(lead) + g.bars + ' ' + g.loop;
          const l3 = "3'-" + postRev.padStart(lead) + g.right + '/';
          return `<div class="psbox"><div class="pshead"><span class="dg ${dgClass(h.dg)}">ΔG ${h.dg.toFixed(2)} kcal/mol</span> <span class="note">stem ${h.stem} bp, loop ${h.loop} nt, at ${fmt(h.start + 1)}–${fmt(h.end + 1)} · ${dgNote(h.dg)}</span></div><pre class="psdiag">${esc(l1)}\n${esc(l2)}\n${esc(l3)}</pre></div>`;
        }).join('');
      };
      $('#psRes').innerHTML = `
        <div class="note" style="margin-bottom:6px"><b>${esc(a.name)}</b> <span class="mono">${esc(a.seq)}</span> · ${a.seq.length} nt · Tm ${DNA.tm(a.seq).toFixed(1)} °C</div>
        ${hpBlock('Hairpins in ' + a.name, a.seq, r.fwd.hairpins)}
        ${dupBlock('Self-dimers of ' + a.name, a.seq, a.seq, r.fwd.self)}
        ${b ? `<div class="note" style="margin:10px 0 6px"><b>${esc(b.name)}</b> <span class="mono">${esc(b.seq)}</span> · ${b.seq.length} nt · Tm ${DNA.tm(b.seq).toFixed(1)} °C</div>
        ${hpBlock('Hairpins in ' + b.name, b.seq, r.rev.hairpins)}
        ${dupBlock('Self-dimers of ' + b.name, b.seq, b.seq, r.rev.self)}
        ${dupBlock('Cross-dimers between the two', a.seq, b.seq, r.cross, 'A cross-dimer that pairs the 3′ ends is the one that produces primer-dimer product.')}` : ''}
        <div class="hint">Free energies use the SantaLucia 1998 nearest-neighbour parameters at 37 °C. A structure weaker than about −6 kcal/mol rarely matters; below about −9 it usually does, especially when the 3′ end is involved.</div>`;
    };
    dialog({
      title: 'Primer Structures', wide: true,
      body: `<div class="form">
        <label for="psA">Primer</label><select id="psA">${list.map(p => `<option value="${p.id}" ${p.id === aId ? 'selected' : ''}>${esc(p.name)}</option>`).join('')}</select>
        <label for="psB">Compare with</label><select id="psB"><option value="">— none —</option>${list.map(p => `<option value="${p.id}" ${String(p.id) === String(bId) ? 'selected' : ''}>${esc(p.name)}</option>`).join('')}</select>
        <div class="full" id="psRes"></div></div>`,
      buttons: [{ label: 'Close', primary: true }],
      onOpen: box => {
        render();
        box.addEventListener('change', e => {
          if (e.target.id === 'psA') aId = +e.target.value;
          else if (e.target.id === 'psB') bId = e.target.value ? +e.target.value : '';
          else return;
          render();
        });
      }
    });
  }

  // ---------------------------------------------------------------- RNA structure
  function rnaDiagram(res, W) {
    // An arc diagram: the sequence on a line, each base pair drawn as an arc above it.
    const n = res.seq.length;
    const m = 16, w = Math.max(200, W - 2 * m);
    const x = i => m + (n <= 1 ? 0 : i / (n - 1) * w);
    const maxSpan = res.pairs.reduce((a, p) => Math.max(a, p[1] - p[0]), 1);
    const H = Math.min(190, 40 + maxSpan / n * 240);
    const base = H - 22;
    const arcs = res.pairs.map(([i, j]) => {
      const x1 = x(i), x2 = x(j), r = (x2 - x1) / 2;
      const h = Math.min(base - 6, 10 + (j - i) / maxSpan * (base - 20));
      return `<path d="M${x1.toFixed(1)},${base} C${x1.toFixed(1)},${(base - h).toFixed(1)} ${x2.toFixed(1)},${(base - h).toFixed(1)} ${x2.toFixed(1)},${base}" fill="none" stroke="#2563c9" stroke-width="1" opacity=".55"/>`;
    }).join('');
    const ticks = [];
    const step = niceStep(n / 8);
    for (let t = 0; t < n; t += step) ticks.push(`<line x1="${x(t).toFixed(1)}" y1="${base}" x2="${x(t).toFixed(1)}" y2="${base + 4}" stroke="#8d95a3"/><text x="${x(t).toFixed(1)}" y="${base + 15}" font-size="10" fill="#8d95a3" text-anchor="middle">${t + 1}</text>`);
    const paired = new Array(n).fill(false);
    for (const [i, j] of res.pairs) { paired[i] = true; paired[j] = true; }
    const dots = n <= 400 ? res.seq.split('').map((c, i) => `<circle cx="${x(i).toFixed(1)}" cy="${base}" r="1.6" fill="${paired[i] ? '#2563c9' : '#c3c9d3'}"/>`).join('') : '';
    return `<svg width="${w + 2 * m}" height="${H}" viewBox="0 0 ${w + 2 * m} ${H}"><line x1="${m}" y1="${base}" x2="${m + w}" y2="${base}" stroke="#c6cbd3"/>${arcs}${dots}${ticks.join('')}</svg>`;
  }
  function dlgRNA() {
    const d = D(), n = d.seq.length;
    const sel = selLen(d);
    let res = null;
    const run = () => {
      const a = parseInt(val('rnFrom'), 10) - 1, b = parseInt(val('rnTo'), 10);
      if (!(a >= 0 && b >= 1 && a < n && b <= n)) { $('#rnRes').innerHTML = '<div class="err">Enter a range inside the sequence.</div>'; return; }
      let s = DNA.sub(d.seq, a, b);
      if ($('#rnRev').checked) s = DNA.rc(s);
      $('#rnRes').innerHTML = '<div class="hint">Folding…</div>';
      // let the message paint before the search starts
      setTimeout(() => {
        const t0 = Date.now();
        try { res = STRUCT.rnaFold(s); } catch (err) { $('#rnRes').innerHTML = `<div class="err">${esc(err.message)}</div>${toolLinks('rna', s)}`; bindCopySeq(STRUCT.toRNA(s), 'RNA sequence'); res = null; return; }
        const ms = Date.now() - t0;
        const stems = STRUCT.stemList(res.structure);
        const rows = [];
        for (let i = 0; i < res.seq.length; i += 60) rows.push(esc(String(a + i + 1).padStart(6) + '  ' + res.seq.slice(i, i + 60)) + '<br>' + esc(' '.repeat(8) + res.structure.slice(i, i + 60)));
        $('#rnRes').innerHTML = `
          <div class="pshead"><span class="dg ${res.dg <= -20 ? 'bad' : res.dg <= -8 ? 'warn' : 'okdg'}">ΔG ≈ ${res.dg.toFixed(2)} kcal/mol</span>
            <span class="note">${fmt(res.seq.length)} nt · ${res.pairs.length} base pairs · ${stems.length} stem${stems.length === 1 ? '' : 's'} · folded in ${ms} ms</span></div>
          <div id="rnDiag"></div>
          ${stems.length ? `<table class="grid"><thead><tr><th>Stem</th><th>Pairs</th><th>Positions</th><th>Loop</th></tr></thead><tbody>${stems.slice(0, 20).map((s2, i) => `<tr><td>${i + 1}</td><td>${s2.len} bp</td><td>${fmt(a + s2.i + 1)}–${fmt(a + s2.j + 1)}</td><td>${fmt(s2.j - s2.i - 1 - 2 * (s2.len - 1))} nt</td></tr>`).join('')}</tbody></table>` : '<div class="hint">No stems found.</div>'}
          <div style="margin-top:8px"><b>Dot-bracket</b></div>
          <div class="aa-box" style="max-height:150px">${rows.join('<br>')}</div>
          <div class="hint"><b>These free energies are approximate.</b> The search is a standard minimum-free-energy fold, but it uses a simplified energy model rather than the full Turner parameter set, so treat ΔG as a rough guide and the structure as a starting point. For a number you would publish, fold the sequence with one of the tools below.</div>
          ${toolLinks('rna', res.seq)}`;
        $('#rnDiag').innerHTML = rnaDiagram(res, Math.min(760, $('#rnRes').clientWidth - 10));
        bindCopySeq(res.seq, 'RNA sequence');
      }, 30);
    };
    dialog({
      title: 'RNA Secondary Structure', wide: true,
      body: `<div class="form">
        <label for="rnFrom">Fold</label><div class="row2">
          <input type="number" id="rnFrom" min="1" max="${n}" value="${sel ? d.sel.start + 1 : 1}"> to
          <input type="number" id="rnTo" min="1" max="${n}" value="${sel ? (d.sel.end || n) : Math.min(n, 200)}">
          <label><input type="checkbox" id="rnRev"> reverse strand</label></div>
        <div class="full hint">T is read as U. Folding is limited to 600 nt here; longer sequences go to an external server.</div>
        <div class="full" id="rnRes"></div></div>`,
      buttons: [{ label: 'Close' }, { label: 'Fold', primary: true, act: () => { run(); return false; } }],
      onOpen: () => run()
    });
  }

  // ---------------------------------------------------------------- protein
  function hydroSVG(vals, W, win) {
    const m = 34, w = Math.max(200, W - m - 16), H = 150, top = 10, bot = H - 22;
    const x = i => m + (vals.length <= 1 ? 0 : i / (vals.length - 1) * w);
    const y = v => bot - (v + 4.5) / 9 * (bot - top);
    const zero = y(0);
    const pts = vals.map((v, i) => `${x(i).toFixed(1)},${y(v).toFixed(1)}`).join(' ');
    const area = `M${m},${zero} L${pts.split(' ').join(' L')} L${(m + w).toFixed(1)},${zero} Z`;
    const ticks = [-4, -2, 0, 2, 4].map(v => `<line x1="${m}" y1="${y(v).toFixed(1)}" x2="${(m + w).toFixed(1)}" y2="${y(v).toFixed(1)}" stroke="${v === 0 ? '#9aa2ae' : '#e8ebef'}"/><text x="${m - 5}" y="${(y(v) + 3.5).toFixed(1)}" font-size="10" fill="#8d95a3" text-anchor="end">${v}</text>`).join('');
    const step = niceStep(vals.length / 7);
    const xt = [];
    for (let t = 0; t < vals.length; t += step) xt.push(`<text x="${x(t).toFixed(1)}" y="${H - 6}" font-size="10" fill="#8d95a3" text-anchor="middle">${t + 1}</text>`);
    return `<svg width="${m + w + 16}" height="${H}" viewBox="0 0 ${m + w + 16} ${H}">${ticks}
      <clipPath id="hcUp"><rect x="${m}" y="${top}" width="${w}" height="${(zero - top).toFixed(1)}"/></clipPath>
      <clipPath id="hcDn"><rect x="${m}" y="${zero.toFixed(1)}" width="${w}" height="${(bot - zero).toFixed(1)}"/></clipPath>
      <path d="${area}" fill="#e8a06a" opacity=".55" clip-path="url(#hcUp)"/>
      <path d="${area}" fill="#8fb3e8" opacity=".55" clip-path="url(#hcDn)"/>
      <polyline points="${pts}" fill="none" stroke="#2f3644" stroke-width="1"/>${xt.join('')}
      <text x="${m}" y="${top - 1}" font-size="10" fill="#8d95a3">hydrophobic ↑ · Kyte–Doolittle, window ${win}</text></svg>`;
  }
  function chargeSVG(curve, pI, W) {
    const m = 34, w = Math.max(200, W - m - 16), H = 150, top = 12, bot = H - 22;
    const lo = Math.min(...curve.map(c => c.charge)), hi = Math.max(...curve.map(c => c.charge));
    const span = Math.max(2, hi - lo);
    const x = pH => m + pH / 14 * w;
    const y = c => bot - (c - lo) / span * (bot - top);
    const pts = curve.map(c => `${x(c.pH).toFixed(1)},${y(c.charge).toFixed(1)}`).join(' ');
    const zero = y(0);
    const xt = [0, 2, 4, 6, 8, 10, 12, 14].map(p => `<text x="${x(p).toFixed(1)}" y="${H - 6}" font-size="10" fill="#8d95a3" text-anchor="middle">${p}</text>`).join('');
    return `<svg width="${m + w + 16}" height="${H}" viewBox="0 0 ${m + w + 16} ${H}">
      <line x1="${m}" y1="${zero.toFixed(1)}" x2="${(m + w).toFixed(1)}" y2="${zero.toFixed(1)}" stroke="#9aa2ae"/>
      <text x="${m - 5}" y="${(zero + 3.5).toFixed(1)}" font-size="10" fill="#8d95a3" text-anchor="end">0</text>
      <line x1="${x(pI).toFixed(1)}" y1="${top}" x2="${x(pI).toFixed(1)}" y2="${bot}" stroke="#c43c3c" stroke-dasharray="3 2"/>
      <text x="${(x(pI) + 4).toFixed(1)}" y="${top + 9}" font-size="10" fill="#c43c3c">pI ${pI.toFixed(2)}</text>
      <polyline points="${pts}" fill="none" stroke="#2563c9" stroke-width="1.4"/>${xt}
      <text x="${m}" y="${top - 2}" font-size="10" fill="#8d95a3">net charge vs pH</text></svg>`;
  }
  function dlgProtein() {
    const d = D(), n = d.seq.length;
    const orfs = analyze(d).orfs;
    const cdsFeatures = d.features.filter(f => f.type === 'CDS' || f.type === 'gene');
    const sel = selLen(d);
    let mode = 'region', stats = null;
    const proteinOf = () => {
      if (mode === 'paste') return val('prPaste').toUpperCase().replace(/[^A-Z*]/g, '');
      const a = parseInt(val('prFrom'), 10) - 1, b = parseInt(val('prTo'), 10);
      if (!(a >= 0 && b >= 1 && a < n && b <= n)) return '';
      let s = DNA.sub(d.seq, a, b);
      if ($('#prRev').checked) s = DNA.rc(s);
      return DNA.translate(s);
    };
    const run = () => {
      const p = proteinOf();
      if (!p) { $('#prRes').innerHTML = '<div class="err">Choose a coding region, or paste a protein sequence.</div>'; return; }
      let st;
      try { st = STRUCT.proteinStats(p); } catch (err) { $('#prRes').innerHTML = `<div class="err">${esc(err.message)}</div>`; return; }
      stats = st;
      const internal = p.slice(0, -1).includes('*');
      const comp = st.composition.filter(c => c.n > 0);
      const seqRows = [];
      for (let i = 0; i < st.protein.length; i += 60) seqRows.push(esc(String(i + 1).padStart(5) + '  ' + st.protein.slice(i, i + 60)));
      $('#prRes').innerHTML = `
        ${internal ? '<div class="err">This translation contains a stop codon in the middle, so the reading frame is probably wrong.</div>' : ''}
        ${st.unknown.length ? `<div class="err">Unrecognised letters ignored: ${esc(st.unknown.join(', '))}</div>` : ''}
        <table class="grid" style="width:auto"><tbody>
          <tr><td>Length</td><td><b>${fmt(st.length)}</b> amino acids</td></tr>
          <tr><td>Molecular weight</td><td><b>${(st.mw / 1000).toFixed(2)} kDa</b> (${fmt(Math.round(st.mw))} Da)</td></tr>
          <tr><td>Isoelectric point</td><td><b>${st.pI.toFixed(2)}</b> <span class="note">(${esc(st.pkaSet)} pKa set)</span></td></tr>
          <tr><td>Net charge at pH 7</td><td><b>${st.charge7 >= 0 ? '+' : ''}${st.charge7.toFixed(2)}</b></td></tr>
          <tr><td>Extinction at 280 nm</td><td><b>${fmt(st.extinction.reduced)}</b> M⁻¹cm⁻¹ with free cysteines${st.extinction.cystine !== st.extinction.reduced ? `, <b>${fmt(st.extinction.cystine)}</b> with all cysteines in cystines` : ''}<br><span class="note">A280 of a 1 mg/mL solution: ${st.a280.reduced.toFixed(3)}${st.extinction.cystine !== st.extinction.reduced ? ` / ${st.a280.cystine.toFixed(3)}` : ''} · ${st.extinction.trp} Trp, ${st.extinction.tyr} Tyr, ${st.extinction.cys} Cys</span></td></tr>
          <tr><td>GRAVY</td><td><b>${st.gravy.toFixed(3)}</b> <span class="note">${st.gravy > 0 ? 'overall hydrophobic' : 'overall hydrophilic'}</span></td></tr>
          <tr><td>Aliphatic index</td><td><b>${st.aliphatic.toFixed(1)}</b></td></tr>
        </tbody></table>
        <div style="margin-top:10px"><b>Hydropathy</b></div><div id="prHydro"></div>
        <div><b>Charge</b></div><div id="prCharge"></div>
        <div style="margin-top:8px"><b>Composition</b></div>
        <table class="grid"><thead><tr><th>Residue</th><th>Count</th><th>%</th><th>Residue</th><th>Count</th><th>%</th></tr></thead><tbody>
          ${Array.from({ length: Math.ceil(comp.length / 2) }, (_, i) => {
            const L = comp[i], Rr = comp[i + Math.ceil(comp.length / 2)];
            const cell = c => c ? `<td>${esc(c.name)} (${c.aa})</td><td>${c.n}</td><td>${c.pct.toFixed(1)}%</td>` : '<td></td><td></td><td></td>';
            return `<tr>${cell(L)}${cell(Rr)}</tr>`;
          }).join('')}
        </tbody></table>
        <div style="margin-top:8px"><b>Sequence</b></div>
        <div class="aa-box" style="max-height:130px">${seqRows.join('<br>')}</div>
        <div class="hint">Mass, extinction coefficient, GRAVY and the aliphatic index are exact calculations. The isoelectric point depends on which pKa set is used, and different tools will differ by a few tenths of a unit. Secondary and three-dimensional structure are not predicted here — use one of the tools below.</div>
        ${toolLinks('protein', st.protein)}`;
      $('#prHydro').innerHTML = hydroSVG(STRUCT.hydropathy(st.protein, parseInt(val('prWin'), 10) || 9), Math.min(740, $('#prRes').clientWidth - 10), parseInt(val('prWin'), 10) || 9);
      $('#prCharge').innerHTML = chargeSVG(STRUCT.chargeCurve(st.protein), st.pI, Math.min(740, $('#prRes').clientWidth - 10));
      bindCopySeq(st.protein, 'protein sequence');
    };
    dialog({
      title: 'Protein Properties', wide: true,
      body: `<div class="form">
        <label>Analyse</label><div class="row2">
          <label><input type="radio" name="prmode" value="region" checked> a coding region</label>
          <label><input type="radio" name="prmode" value="paste"> a pasted protein</label></div>
        <div class="full" id="prRegionBox"><div class="form" style="margin:0">
          <label for="prFrom">Region</label><div class="row2">
            <input type="number" id="prFrom" min="1" max="${n}" value="${sel ? d.sel.start + 1 : 1}"> to
            <input type="number" id="prTo" min="1" max="${n}" value="${sel ? (d.sel.end || n) : Math.min(n, 300)}">
            <label><input type="checkbox" id="prRev"> reverse strand</label></div>
          <label>Use</label><select id="prPick"><option value="">choose a feature or ORF…</option>
            ${cdsFeatures.map((f, i) => `<option value="f${i}">${esc(f.name)} (${esc(f.type)}, ${fmt(DNA.rlen(f, n))} bp)</option>`).join('')}
            ${orfs.map((o, i) => `<option value="o${i}">ORF ${fmt(o.start + 1)}–${fmt(o.end)} (${o.aa} aa, ${o.strand > 0 ? 'forward' : 'reverse'})</option>`).join('')}</select>
        </div></div>
        <div class="full" id="prPasteBox" style="display:none"><div class="form" style="margin:0">
          <label for="prPaste">Protein</label><textarea id="prPaste" rows="4" placeholder="Paste a protein sequence in one-letter code"></textarea>
        </div></div>
        <label for="prWin">Hydropathy window</label><div class="row2"><input type="number" id="prWin" value="9" min="3" max="31"> residues</div>
        <div class="full" id="prRes"></div></div>`,
      buttons: [{ label: 'Close' }, { label: 'Analyse', primary: true, act: () => { run(); return false; } }],
      onOpen: box => {
        run();
        box.addEventListener('change', e => {
          const t = e.target;
          if (t.name === 'prmode') {
            mode = t.value;
            $('#prRegionBox').style.display = mode === 'region' ? '' : 'none';
            $('#prPasteBox').style.display = mode === 'paste' ? '' : 'none';
          } else if (t.id === 'prPick') {
            const v2 = t.value; if (!v2) return;
            if (v2[0] === 'f') { const f = cdsFeatures[+v2.slice(1)]; $('#prFrom').value = f.start + 1; $('#prTo').value = f.end; $('#prRev').checked = f.strand === -1; }
            else { const o = orfs[+v2.slice(1)]; $('#prFrom').value = o.start + 1; $('#prTo').value = o.end; $('#prRev').checked = o.strand === -1; }
            run();
          } else if (/^(prFrom|prTo|prRev|prWin)$/.test(t.id)) run();
        });
      }
    });
  }

  // ================================================================ printing
  function buildPrintHTML(d, opts) {
    const A = analyze(d), n = d.seq.length;
    const parts = [];
    const css = `@page{margin:14mm}
      body{font:11pt/1.4 "Segoe UI",Arial,sans-serif;color:#111;margin:0}
      h1{font-size:17pt;margin:0 0 2pt}h2{font-size:12pt;margin:16pt 0 4pt;border-bottom:1px solid #bbb;padding-bottom:2pt}
      .sub{color:#555;font-size:9.5pt;margin-bottom:10pt}
      table{border-collapse:collapse;width:100%;font-size:9pt;margin-top:4pt}
      th{text-align:left;border-bottom:1px solid #999;padding:2pt 4pt}
      td{border-bottom:1px solid #e3e3e3;padding:2pt 4pt;vertical-align:top}
      .mono{font-family:Consolas,"DejaVu Sans Mono",monospace}
      .seqblock{font-family:Consolas,"DejaVu Sans Mono",monospace;font-size:9.2pt;line-height:1.5;white-space:pre}
      svg{max-width:100%;height:auto}
      .page{page-break-after:always}
      .page:last-child{page-break-after:auto}
      ${opts.map ? mapCSS() : ''}`;
    parts.push(`<h1>${esc(d.name)}</h1><div class="sub">${fmt(n)} bp · ${d.circular ? 'circular' : 'linear'} · GC ${DNA.gcPct(d.seq).toFixed(1)}% · ${d.features.length} features · printed ${new Date().toLocaleDateString()}</div>`);
    if (opts.map) {
      const svg = $('#mapwrap svg');
      let markup = '';
      if (svg) markup = svg.outerHTML.replace(/<g id="sellayer">[\s\S]*?<\/g>/, '');
      else markup = (d.opts.mapMode === 'circular' ? drawCircular(d, A, 760, 760) : drawLinear(d, A, 760)).replace(/<g id="sellayer">[\s\S]*?<\/g>/, '');
      markup = markup.replace(/style="fill:var\(--canvas\)"/g, 'fill="#ffffff"');
      parts.push(`<div class="page"><h2>Map</h2>${markup}</div>`);
    }
    if (opts.features && d.features.length) {
      const fs = d.features.slice().sort((a, b) => a.start - b.start);
      parts.push(`<h2>Features</h2><table><thead><tr><th>Name</th><th>Type</th><th>Location</th><th>Length</th><th>Direction</th></tr></thead><tbody>${fs.map(f => `<tr><td>${esc(f.name)}</td><td>${esc(f.type)}</td><td>${fmt(f.start + 1)}..${fmt(f.end)}</td><td>${fmt(DNA.rlen(f, n))} bp</td><td>${f.strand === 1 ? 'forward' : f.strand === -1 ? 'reverse' : '—'}</td></tr>`).join('')}</tbody></table>`);
    }
    if (opts.primers && d.primers.length) {
      parts.push(`<h2>Primers</h2><table><thead><tr><th>Name</th><th>Sequence 5′→3′</th><th>Length</th><th>Tm</th><th>Binds</th></tr></thead><tbody>${d.primers.map(p => { const b = A.psites[p.id] || []; return `<tr><td>${esc(p.name)}</td><td class="mono">${esc(p.seq)}</td><td>${p.seq.length}</td><td>${DNA.tm(p.seq).toFixed(1)} °C</td><td>${b.length ? b.map(x => `${fmt(x.start + 1)}..${fmt(x.end)} ${x.strand > 0 ? '→' : '←'}`).join(', ') : 'no site'}</td></tr>`; }).join('')}</tbody></table>`);
    }
    if (opts.enzymes) {
      const list = shownEnz(d, A);
      if (list.length) parts.push(`<h2>Restriction sites</h2><table><thead><tr><th>Enzyme</th><th>Site</th><th># Sites</th><th>Cuts after</th></tr></thead><tbody>${list.map(e => `<tr><td>${esc(e.name)}</td><td class="mono">${esc(DNA.ENZ_BY_NAME[e.name].site)}</td><td>${e.sites.length}</td><td>${e.sites.map(s => fmt(s.cut)).join(', ')}</td></tr>`).join('')}</tbody></table>`);
    }
    if (opts.sequence) {
      const per = 60, lines = [];
      for (let i = 0; i < n; i += per) {
        let row = String(i + 1).padStart(7) + ' ';
        for (let j = i; j < Math.min(i + per, n); j += 10) row += d.seq.slice(j, Math.min(j + 10, n)) + ' ';
        lines.push(esc(row.trimEnd()));
      }
      parts.push(`<h2>Sequence</h2><div class="seqblock">${lines.join('\n')}</div>`);
    }
    if (opts.alignments && d.alignments.length) {
      parts.push(`<h2>Sequencing alignments</h2><table><thead><tr><th>Read</th><th>Strand</th><th>Covers</th><th>Identity</th><th>Differences</th></tr></thead><tbody>${d.alignments.map(a => `<tr><td>${esc(a.read.name)}</td><td>${a.strand > 0 ? 'forward' : 'reverse'}</td><td>${fmt(a.refStart + 1)}..${fmt(a.refEnd)}</td><td>${a.identity.toFixed(1)}%</td><td>${a.diffs.length ? esc(a.diffs.map(x => `${fmt(x.pos + 1)} ${x.type === 'mis' ? x.ref + '→' + x.read : x.type}`).join(', ')) : 'none'}</td></tr>`).join('')}</tbody></table>`);
    }
    return `<!doctype html><html><head><meta charset="utf-8"><title>${esc(d.name)}</title><style>${css}</style></head><body>${parts.join('')}</body></html>`;
  }
  function dlgPrint(pdf) {
    const d = D();
    const cb = (id, label, on) => `<label><input type="checkbox" id="${id}" ${on ? 'checked' : ''}> ${label}</label>`;
    dialog({
      title: pdf ? 'Export PDF' : 'Print',
      body: `<div class="form"><label>Include</label><div>
        ${cb('prmap', 'Map', true)}${cb('prfeat', 'Feature table', true)}${cb('prprim', 'Primer table', d.primers.length > 0)}
        ${cb('prenz', 'Restriction sites', false)}${cb('prseq', 'Full sequence', false)}${cb('pralign', 'Sequencing alignments', d.alignments.length > 0)}</div>
        <div class="full hint">The map is printed as it is currently shown, so switch between the circular and linear map before printing.</div>
        <div class="full err"></div></div>`,
      buttons: [{ label: 'Cancel' }, {
        label: pdf ? 'Save PDF…' : 'Print…', primary: true, act: () => {
          const opts = { map: $('#prmap').checked, features: $('#prfeat').checked, primers: $('#prprim').checked, enzymes: $('#prenz').checked, sequence: $('#prseq').checked, alignments: $('#pralign').checked };
          if (!Object.values(opts).some(Boolean)) return showErr('Choose at least one section.');
          const html = buildPrintHTML(d, opts);
          closeDialog();
          doPrint(html, pdf, d);
          return false;
        }
      }]
    });
  }
  async function doPrint(html, pdf, d) {
    try {
      if (host && host.printHtml) {
        if (pdf) {
          const p = await host.pdfHtml({ html, defaultName: d.name.replace(/[\\/:*?"<>|]/g, '_') + '.pdf', landscape: false });
          if (p) toast('Saved ' + p.split(/[\\/]/).pop());
        } else {
          const r = await host.printHtml(html);
          if (r && !r.ok && r.reason && r.reason !== 'cancelled') toast('Printing failed: ' + r.reason);
        }
        return;
      }
      const w = window.open('', '_blank');
      if (!w) { toast('Allow pop-up windows to print.'); return; }
      w.document.write(html); w.document.close();
      w.onload = () => { w.focus(); w.print(); };
    } catch (err) { toast('Printing failed: ' + err.message); }
  }

  // ================================================================ files
  function addDoc(src) {
    const d = makeDoc(src);
    S.docs.push(d); S.active = S.docs.length - 1;
    if (S.view === 'hist') S.view = 'map';
    renderAll();
    return d;
  }
  function closeDoc(i, force) {
    const d = S.docs[i]; if (!d) return;
    if (isGel(d)) { S.docs.splice(i, 1); if (S.active >= S.docs.length) S.active = S.docs.length - 1; renderAll(); return; }
    if (d.dirty && !force) {
      dialog({
        title: 'Unsaved Changes', body: `<p>Save changes to “${esc(d.name)}” before closing?</p>`,
        buttons: [{ label: 'Don’t Save', danger: true, act: () => { setTimeout(() => closeDoc(i, true), 0); } }, { label: 'Cancel' }, { label: 'Save', primary: true, act: () => { S.active = i; saveDoc(false).then(ok => { if (ok) closeDoc(i, true); }); } }]
      });
      return;
    }
    S.docs.splice(i, 1);
    if (S.active >= S.docs.length) S.active = S.docs.length - 1;
    renderAll();
  }
  function isTraceFile(name, data) {
    if (/\.ab1$/i.test(name)) return true;
    if (data instanceof ArrayBuffer && data.byteLength > 4) {
      const dv = new DataView(data);
      return String.fromCharCode(dv.getUint8(0), dv.getUint8(1), dv.getUint8(2), dv.getUint8(3)) === 'ABIF';
    }
    return false;
  }
  function loadContent(name, data, path) {
    if (isTraceFile(name, data)) {
      if (!D()) { dialog({ title: 'Open the Reference First', body: `<p>${esc(name)} is a sequencing trace. Open the sequence it should be compared against before adding traces.</p>`, buttons: [{ label: 'OK', primary: true }] }); return; }
      addTrace(name, data);
      return;
    }
    try {
      let src;
      if (/\.dna$/i.test(name)) src = DNA.parseDnaBinary(data instanceof ArrayBuffer ? data : new TextEncoder().encode(data).buffer, name);
      else src = DNA.parseAny(typeof data === 'string' ? data : new TextDecoder().decode(data), name);
      src.path = path && /\.(gb|gbk|genbank|gbff)$/i.test(path) ? path : null;
      if (path && !src.path) src.histLabel = 'Imported ' + name;
      addDoc(src);
      toast(`Opened ${name} — ${fmt(src.seq.length)} bp, ${src.features.length} features`);
    } catch (err) {
      dialog({ title: 'Could Not Open File', body: `<p>${esc(name)}</p><p class="err">${esc(err.message)}</p><p class="hint">GeneCanvas opens GenBank (.gb, .gbk), FASTA (.fa, .fasta), .dna and plain-text sequences.</p>`, buttons: [{ label: 'OK', primary: true }] });
    }
  }
  async function openFiles() {
    if (host && host.openFiles) {
      const files = await host.openFiles();
      for (const f of files || []) {
        const bin = Uint8Array.from(atob(f.data), c => c.charCodeAt(0));
        loadContent(f.name, /\.(dna|ab1)$/i.test(f.name) ? bin.buffer : new TextDecoder().decode(bin), f.path);
      }
      return;
    }
    $('#fileinput').click();
  }
  $('#fileinput').addEventListener('change', e => {
    for (const f of e.target.files) readFileObj(f);
    e.target.value = '';
  });
  function readFileObj(f) {
    const r = new FileReader();
    r.onload = () => loadContent(f.name, r.result, f.path || null);
    if (/\.(dna|ab1)$/i.test(f.name)) r.readAsArrayBuffer(f); else r.readAsText(f);
  }
  function primerSitesFor(d) { const A = analyze(d); return p => A.psites[p.id] || []; }
  async function saveDoc(asNew) {
    const d = D(); if (!d || isGel(d)) return false;
    const content = DNA.toGenBank(d, primerSitesFor(d));
    const fname = d.name.replace(/[\\/:*?"<>|]/g, '_') + '.gb';
    if (host && host.saveFile) {
      const res = await host.saveFile({ path: asNew ? null : d.path, defaultName: fname, content, filters: [{ name: 'GenBank', extensions: ['gb', 'gbk'] }] });
      if (!res) return false;
      d.path = res; d.dirty = false; renderAll(true); toast('Saved ' + res.split(/[\\/]/).pop()); return true;
    }
    downloadText(fname, content); d.dirty = false; renderAll(true); return true;
  }
  async function exportFasta() {
    const d = D(); const content = DNA.toFASTA(d.name, d.seq), fname = d.name.replace(/[\\/:*?"<>|]/g, '_') + '.fasta';
    if (host && host.saveFile) { const r = await host.saveFile({ path: null, defaultName: fname, content, filters: [{ name: 'FASTA', extensions: ['fasta', 'fa'] }] }); if (r) toast('Exported ' + r.split(/[\\/]/).pop()); }
    else downloadText(fname, content);
  }
  // The map's SVG is styled by class, so any copy of it has to carry those rules with it.
  function mapCSS() {
    return `:root{--canvas:#fff;--bb:#2f3644;--faint:#8d95a3;--ink:#1c2230;--muted:#5e6776;--accent:#2563c9;--sel:rgba(0,0,0,0);--selline:transparent}\n${SVG_CSS}`;
  }
  async function exportSVG() {
    const d = D();
    if (S.view !== 'map') { S.view = 'map'; renderAll(); }
    const svg = $('#mapwrap svg'); if (!svg) return;
    const css = `<style>${mapCSS()}</style>`;
    const content = svg.outerHTML.replace(/<g id="sellayer">[\s\S]*?<\/g>/, '').replace('<defs>', css + '<defs>');
    const fname = d.name.replace(/[\\/:*?"<>|]/g, '_') + '-map.svg';
    if (host && host.saveFile) { const r = await host.saveFile({ path: null, defaultName: fname, content, filters: [{ name: 'SVG image', extensions: ['svg'] }] }); if (r) toast('Exported ' + r.split(/[\\/]/).pop()); }
    else downloadText(fname, content);
  }
  function downloadText(name, text) {
    const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([text], { type: 'text/plain' })); a.download = name; document.body.appendChild(a); a.click(); setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
  }
  function copyText(t, msg) {
    const done = () => toast(msg || 'Copied');
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(t).then(done, () => fallbackCopy(t, done));
    else fallbackCopy(t, done);
  }
  function fallbackCopy(t, done) { const ta = document.createElement('textarea'); ta.value = t; document.body.appendChild(ta); ta.select(); try { document.execCommand('copy'); done(); } catch (e) { toast('Copy failed'); } ta.remove(); }

  let toastT = null;
  function toast(msg) { const t = $('#toast'); t.textContent = msg; t.hidden = false; clearTimeout(toastT); toastT = setTimeout(() => t.hidden = true, 2600); }

  // ================================================================ global events
  document.addEventListener('mousedown', e => {
    if (!e.target.closest('#menupop') && !e.target.closest('#menubar')) closeMenu();
    if (!e.target.closest('#ctx')) closeCtx();
  });
  $('#toolbar').addEventListener('click', e => { const b = e.target.closest('[data-a]'); if (b && !b.disabled) TB_ACT[b.dataset.a](); });
  $('#doctabs').addEventListener('mousedown', e => {
    const x = e.target.closest('[data-x]'); if (x) { e.stopPropagation(); closeDoc(+x.dataset.x); return; }
    const t = e.target.closest('.dtab'); if (t) { S.active = +t.dataset.i; renderAll(); }
  });
  $('#viewtabs').addEventListener('click', e => { const b = e.target.closest('[data-v]'); if (b) setView(b.dataset.v); });
  $('#side').addEventListener('change', e => {
    const d = D(); if (!d) return;
    const t = e.target;
    if (t.dataset.o) setOpt(t.dataset.o, t.checked);
    else if (t.id === 'side-enz') { if (t.value === 'chosen') { setOpt('enzMode', 'chosen'); dlgEnzymeChooser(); } else setOpt('enzMode', t.value); }
    else if (t.id === 'side-orfmin') { const v = Math.max(20, Math.min(2000, parseInt(t.value, 10) || 100)); setOpt('orfMin', v); }
  });
  $('#side').addEventListener('click', e => {
    const b = e.target.closest('[data-mm]'); if (b) { setOpt('mapMode', b.dataset.mm); return; }
    if (e.target.id === 'side-pick') dlgEnzymeChooser();
  });

  document.addEventListener('keydown', e => {
    if (!$('#modal').hidden) return;
    if (e.key === 'Escape') { closeMenu(); closeCtx(); return; }
    const mod = e.ctrlKey || e.metaKey;
    const inField = /^(INPUT|TEXTAREA|SELECT)$/.test(e.target.tagName);
    if (!mod) { if (!inField && (e.key === 'Delete' || e.key === 'Backspace') && D() && S.view === 'map' && selLen(D())) { e.preventDefault(); deleteSel(); } return; }
    const k = e.key.toLowerCase();
    const map = {
      n: dlgNew, o: openFiles, s: () => saveDoc(e.shiftKey), w: () => D() && closeDoc(S.active), f: openFind, g: () => D() && dlgGoto(), k: () => D() && dlgFeature(null),
      d: dlgNCBI, p: () => D() && dlgPrint(),
      '1': () => setView('map'), '2': () => setView('seq'), '3': () => setView('enz'), '4': () => setView('feat'), '5': () => setView('primer'), '6': () => setView('align'), '7': () => setView('hist'),
      '=': () => zoom(1.25), '+': () => zoom(1.25), '-': () => zoom(0.8), '0': () => { S.zoom = 1; renderView(); }
    };
    if (!inField) {
      map.z = () => (e.shiftKey ? redo() : undo()); map.y = redo;
      map.a = () => D() && setSel(0, D().seq.length, { anchor: 0 });
    }
    if (map[k]) { e.preventDefault(); map[k](); }
  });
  document.addEventListener('copy', e => {
    if (/^(INPUT|TEXTAREA)$/.test(document.activeElement.tagName)) return;
    const d = D(); if (!d || !selLen(d)) return;
    e.clipboardData.setData('text/plain', selSeq(d)); e.preventDefault(); toast('Copied ' + fmt(selLen(d)) + ' bp');
  });
  document.addEventListener('cut', e => {
    if (/^(INPUT|TEXTAREA)$/.test(document.activeElement.tagName)) return;
    const d = D(); if (!d || !selLen(d) || S.view !== 'seq') return;
    e.clipboardData.setData('text/plain', selSeq(d)); e.preventDefault(); deleteSel();
  });
  // drag & drop
  let dz = null;
  window.addEventListener('dragover', e => { e.preventDefault(); if (!dz) { dz = document.createElement('div'); dz.className = 'dropzone'; dz.textContent = 'Drop to open'; document.body.appendChild(dz); } });
  window.addEventListener('dragleave', e => { if (e.relatedTarget === null && dz) { dz.remove(); dz = null; } });
  window.addEventListener('drop', e => { e.preventDefault(); if (dz) { dz.remove(); dz = null; } for (const f of e.dataTransfer.files) readFileObj(f); });
  let rT = null;
  window.addEventListener('resize', () => { clearTimeout(rT); rT = setTimeout(() => { if (S.view === 'map' || S.view === 'seq') renderView(true); }, 120); });

  // host integration (Electron)
  if (host) {
    if (host.onOpenFile) host.onOpenFile(f => { const bin = Uint8Array.from(atob(f.data), c => c.charCodeAt(0)); loadContent(f.name, /\.(dna|ab1)$/i.test(f.name) ? bin.buffer : new TextDecoder().decode(bin), f.path); });
    if (host.onCloseRequest) host.onCloseRequest(() => {
      const dirty = S.docs.filter(d => d.dirty);
      if (!dirty.length) { host.confirmClose(); return; }
      dialog({ title: 'Quit GeneCanvas', body: `<p>${dirty.length} file${dirty.length > 1 ? 's have' : ' has'} unsaved changes: ${dirty.map(d => '“' + esc(d.name) + '”').join(', ')}.</p><p>Quit without saving?</p>`, buttons: [{ label: 'Cancel', primary: true }, { label: 'Quit Without Saving', danger: true, act: () => host.confirmClose() }] });
    });
  }

  // ================================================================ start
  async function start() {
    await loadSettings();
    const st = document.createElement('style');
    st.textContent = SVG_CSS;
    document.head.appendChild(st);
    measureCW();
    renderMenuBar();
    addDoc(DNA.buildDemo());
    S.docs[0].dirty = false;
    renderAll();
    if (host && host.ready) host.ready();
  }
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(start); else start();
  // Test and automation hook: run a menu command by its labels.
  function menuAct(menuLabel, itemLabel) {
    const m = menus().find(x => x.label === menuLabel);
    if (!m) throw new Error('No menu ' + menuLabel);
    const it = m.items.find(x => x !== '-' && x.label === itemLabel);
    if (!it) throw new Error('No item ' + itemLabel + ' in ' + menuLabel);
    if (it.dis) throw new Error('Disabled: ' + itemLabel);
    return it.act();
  }
  window.__gc = { S, DNA, addDoc, setView, setSel, renderAll, menuAct, addTrace, D, buildPrintHTML, get SET() { return SET; }, saveSettings };
})();
