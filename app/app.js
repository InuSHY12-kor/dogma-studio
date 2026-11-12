/* GeneCanvas — application UI */
(function () {
  'use strict';
  const DNA = window.DNA;
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
  const FEATURE_TYPES = Object.keys(TYPE_COL).filter(t => t !== 'ORF');
  const CDS_PALETTE = ['#f4c95d', '#8ec3f2', '#f29e8e', '#a6dd8b', '#d6a8f0', '#f6b26b', '#7fd1c7'];

  // ---------------------------------------------------------------- state
  const S = { docs: [], active: -1, view: 'map', zoom: 1, sideHidden: false, seq: null, mapGeom: null, find: { q: '', mode: 'dna', res: [], i: -1 } };
  let UID = 1; const uid = () => UID++;
  const D = () => S.docs[S.active] || null;

  function defaultOpts(circular) {
    return { mapMode: circular ? 'circular' : 'linear', showFeat: true, showPrimers: true, showORFs: false, orfMin: 100, showTrans: true, enzMode: 'unique6', showComp: true };
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
      opts: defaultOpts(!!src.circular), chosen: [], sel: { start: 0, end: 0 }, anchor: 0,
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
    zoomin: '<circle cx="10.5" cy="10.5" r="6"/><path d="m15 15 6 6M8 10.5h5M10.5 8v5"/>',
    zoomout: '<circle cx="10.5" cy="10.5" r="6"/><path d="m15 15 6 6M8 10.5h5"/>'
  };
  const icon = k => `<svg viewBox="0 0 24 24">${IC[k]}</svg>`;

  // ---------------------------------------------------------------- menus
  function menus() {
    const d = D(), has = !!d, sl = has && selLen(d) > 0;
    return [
      { label: 'File', items: [
        { label: 'New DNA File…', key: 'Ctrl+N', act: dlgNew },
        { label: 'Open…', key: 'Ctrl+O', act: openFiles },
        { label: 'Open Example Plasmid', act: () => addDoc(DNA.buildDemo()) },
        '-',
        { label: 'Save', key: 'Ctrl+S', dis: !has, act: () => saveDoc(false) },
        { label: 'Save As…', key: 'Ctrl+Shift+S', dis: !has, act: () => saveDoc(true) },
        { label: 'Export FASTA…', dis: !has, act: exportFasta },
        { label: 'Export Map Image (SVG)…', dis: !has, act: exportSVG },
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
    const d = D(), has = !!d;
    const b = (k, label, act, o = {}) => `<button class="tb${o.on ? ' on' : ''}" data-a="${act}" title="${esc(o.title || label)}" ${o.dis ? 'disabled' : ''}>${icon(k)}<span>${label}</span></button>`;
    $('#toolbar').innerHTML = [
      b('new', 'New', 'new'), b('open', 'Open', 'open'), b('save', 'Save', 'save', { dis: !has }),
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
    tr: () => dlgTranslate(), rc: () => revCompSel(), origin: () => setOrigin(), side: () => { S.sideHidden = !S.sideHidden; renderAll(); }
  };

  function renderDocTabs() {
    $('#doctabs').innerHTML = S.docs.map((d, i) => `<div class="dtab${i === S.active ? ' on' : ''}" data-i="${i}" title="${esc(d.name)} — ${fmt(d.seq.length)} bp"><span class="nm">${d.dirty ? '<span class="dirty">● </span>' : ''}${esc(d.name)}</span><button class="x" data-x="${i}" title="Close">×</button></div>`).join('');
  }
  const VIEWS = [['map', 'Map', 'map'], ['seq', 'Sequence', 'seq'], ['enz', 'Enzymes', 'enz'], ['feat', 'Features', 'feat'], ['primer', 'Primers', 'primer'], ['hist', 'History', 'hist']];
  function renderViewTabs() {
    const vt = $('#viewtabs');
    vt.classList.toggle('noside', S.sideHidden || !(S.view === 'map' || S.view === 'seq'));
    vt.innerHTML = D() ? VIEWS.map(([k, l, ic]) => `<button class="vt${S.view === k ? ' on' : ''}" data-v="${k}">${icon(ic)}${l}</button>`).join('') : '';
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
    const show = d && !S.sideHidden && (S.view === 'map' || S.view === 'seq');
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
      <option value="none" ${o.enzMode === 'none' ? 'selected' : ''}>None</option></select>`;
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
    else if (S.view === 'hist') renderHistory();
  }

  // ================================================================ MAP VIEW
  function renderMap(keepScroll) {
    const v = $('#view');
    let wrap = $('#mapwrap');
    const prev = wrap ? { t: wrap.scrollTop, l: wrap.scrollLeft } : null;
    v.innerHTML = `<div class="mapwrap" id="mapwrap"></div><div class="map-zoom"><button data-z="out" title="Zoom out">−</button><button data-z="reset" title="Fit">⤢</button><button data-z="in" title="Zoom in">+</button></div>`;
    wrap = $('#mapwrap');
    const d = D(), A = analyze(d);
    const W = Math.max(420, wrap.clientWidth), H = Math.max(380, wrap.clientHeight);
    const svg = d.opts.mapMode === 'circular' ? drawCircular(d, A, W * S.zoom, H * S.zoom) : drawLinear(d, A, W * S.zoom);
    wrap.innerHTML = svg;
    if (keepScroll && prev) { wrap.scrollTop = prev.t; wrap.scrollLeft = prev.l; }
    drawSelLayer();
    bindMap(wrap.querySelector('svg'));
    v.querySelector('.map-zoom').onclick = e => { const b = e.target.closest('button'); if (!b) return; if (b.dataset.z === 'in') zoom(1.25); else if (b.dataset.z === 'out') zoom(0.8); else { S.zoom = 1; renderView(); } };
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
    const nF = assignLanes(feats, n), nP = assignLanes(prims, n), nO = assignLanes(orfs, n);
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
    for (let t = step; t < n; t += step) {
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
      if (avail > textW) {
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
    const nF = pxLanes(feats, true), nP = pxLanes(prims, true), nO = pxLanes(orfs, false);
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
    for (let t = 0; t <= n; t += step) {
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
          if (it.kind === 'f' && x2 - x1 - hd > tw) g.push(`<text class="ftxt" x="${((x1 + x2) / 2).toFixed(1)}" y="${ym}" dy=".35em" text-anchor="middle">${esc(it.name)}</text>`);
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
    const G = 72;
    const bpr = Math.max(10, Math.floor((sc.clientWidth - G - 26) / CW / 10) * 10);
    // enzyme cuts
    const cuts = [];
    for (const e of shownEnz(d, A)) e.sites.forEach((s, i) => { if (s.cut != null) cuts.push({ pos: s.cut, name: e.name, k: `e:${e.name}:${i}`, u: e.unique }); });
    cuts.sort((a, b) => a.pos - b.pos);
    // feature segments
    const items = collectItems(d, A);
    const segsAll = [];
    for (const it of items) for (const [a, b] of DNA.segs(it, n)) segsAll.push({ a, b, it, tipR: it.strand === 1 && b === it.end, tipL: it.strand === -1 && a === it.start });
    // translations
    const aaRows = {};
    if (o.showTrans) for (const it of items) {
      if (!((it.kind === 'f' && it.type === 'CDS') || it.kind === 'o')) continue;
      const L = DNA.rlen(it, n); const fs = DNA.sub(d.seq, it.start, it.end);
      const s = it.strand === -1 ? DNA.rc(fs) : fs;
      for (let j = 0; j + 3 <= L; j += 3) {
        const aa = DNA.CODON[s.substr(j, 3)] || 'X';
        const mid = it.strand === -1 ? ((it.end - j - 2) % n + n) % n : (it.start + j + 1) % n;
        const r = Math.floor(mid / bpr);
        (aaRows[r] = aaRows[r] || []).push({ mid, aa, key: it.key });
      }
    }
    const rows = [], html = [];
    const comp = o.showComp;
    for (let rs = 0, ri = 0; rs < n; rs += bpr, ri++) {
      const re = Math.min(n, rs + bpr), len = re - rs;
      const parts = [];
      // enzymes
      const lineEnd = [];
      const myCuts = cuts.filter(c => c.pos >= rs && (c.pos < re || (re === n && c.pos === n)));
      for (const c of myCuts) {
        const x = (c.pos - rs) * CW, w = c.name.length * 6.4 + 6;
        let l = 0; while (lineEnd[l] !== undefined && lineEnd[l] > x - 1) l++;
        lineEnd[l] = x + w; c.line = l;
      }
      const nL = lineEnd.length, enzH = nL ? nL * 13 + 5 : 2;
      const y0 = 4 + enzH;
      for (const c of myCuts) {
        const x = G + (c.pos - rs) * CW, y = 4 + (nL - 1 - c.line) * 13;
        parts.push(`<div class="enz${c.u ? ' u' : ''}" data-k="${c.k}" style="left:${x - 1}px;top:${y}px">${esc(c.name)}</div><div class="etick" style="left:${x}px;top:${y + 12}px;height:${y0 - y - 11}px"></div>`);
      }
      parts.push(`<div class="num" style="top:${y0 + 1}px">${fmt(rs + 1)}</div>`);
      parts.push(`<div class="nt" style="left:${G}px;top:${y0}px">${d.seq.slice(rs, re)}</div>`);
      const rulerTop = y0 + 17;
      parts.push(`<div class="ruler" style="left:${G}px;top:${rulerTop}px;width:${len * CW}px;background:repeating-linear-gradient(to right,transparent 0,transparent ${(CW * 10 - 1).toFixed(2)}px,#9aa2ae ${(CW * 10 - 1).toFixed(2)}px,#9aa2ae ${(CW * 10).toFixed(2)}px);background-position:${(-CW * 0.5).toFixed(2)}px 0"></div>`);
      let y = rulerTop + 6;
      if (comp) { parts.push(`<div class="nt b" style="left:${G}px;top:${y}px">${DNA.comp(d.seq.slice(rs, re))}</div>`); y += 17; }
      const strandBottom = y;
      // features in this row
      const rsegs = segsAll.filter(s => s.a < re && s.b > rs).map(s => ({ ...s, x1: Math.max(s.a, rs), x2: Math.min(s.b, re) })).sort((p, q) => p.x1 - q.x1 || (q.x2 - q.x1) - (p.x2 - p.x1));
      const laneEnd = [], laneTrans = [];
      for (const s of rsegs) {
        let l = 0; while (laneEnd[l] !== undefined && laneEnd[l] > s.x1) l++;
        laneEnd[l] = s.x2; s.lane = l;
        if (o.showTrans && ((s.it.kind === 'f' && s.it.type === 'CDS') || s.it.kind === 'o')) laneTrans[l] = true;
      }
      y += 5;
      const laneY = [];
      for (let l = 0; l < laneEnd.length; l++) { laneY[l] = y; y += 16 + (laneTrans[l] ? 14 : 0); }
      const keyLane = {};
      for (const s of rsegs) {
        const x = G + (s.x1 - rs) * CW, w = (s.x2 - s.x1) * CW;
        let clip = '';
        if (s.tipR && s.x2 === s.b) clip = 'clip-path:polygon(0 0,calc(100% - 7px) 0,100% 50%,calc(100% - 7px) 100%,0 100%);';
        else if (s.tipL && s.x1 === s.a) clip = 'clip-path:polygon(7px 0,100% 0,100% 100%,7px 100%,0 50%);';
        const label = w > s.it.name.length * 6 + 14 ? esc(s.it.name) : '';
        const h = s.it.kind === 'o' ? 7 : 12;
        parts.push(`<div class="fbar" data-k="${s.it.key}" title="${esc(s.it.name)} ${fmt(s.it.start + 1)}..${fmt(s.it.end)}" style="left:${x}px;top:${laneY[s.lane]}px;width:${w}px;height:${h}px;background:${s.it.color};${clip}${s.tipL && s.x1 === s.a ? 'padding-left:9px;' : ''}">${s.it.kind === 'o' ? '' : label}</div>`);
        keyLane[s.it.key] = s.lane;
      }
      for (const a of (aaRows[ri] || [])) {
        const l = keyLane[a.key]; if (l === undefined) continue;
        const x = G + (a.mid - rs) * CW - CW;
        parts.push(`<div class="aa${a.aa === '*' ? ' st' : (a.aa === 'M' ? ' m' : '')}" style="left:${x}px;top:${laneY[l] + 14}px;width:${CW * 3}px">${a.aa}</div>`);
      }
      const hgt = y + 6;
      rows.push({ rs, re, y0, y1: strandBottom });
      html.push(`<div class="srow" data-i="${ri}" style="height:${hgt}px">${parts.join('')}<div class="selov" hidden style="top:${y0 - 1}px;height:${strandBottom - y0 + 2}px"></div><div class="caret-el" hidden style="top:${y0 - 1}px;height:${strandBottom - y0 + 2}px"></div></div>`);
    }
    host.innerHTML = html.join('');
    S.seq = { rows, G, bpr, els: [...host.querySelectorAll('.srow')] };
    updateSeqSel();
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
      <button class="btn" id="enzshow">Show chosen on map</button><button class="btn" id="enzclear">Clear chosen</button><button class="btn pri" id="enzdig" ${d.chosen.length ? '' : 'disabled'}>Digest chosen (${d.chosen.length})…</button>
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
    v.innerHTML = `<div class="panel"><div class="ptools"><button class="btn pri" id="padd">Add Primer…</button><button class="btn" id="ptop" ${selLen(d) ? '' : 'disabled'}>From Selection (Top)</button><button class="btn" id="pbot" ${selLen(d) ? '' : 'disabled'}>From Selection (Bottom)</button><span class="note">Tm: nearest-neighbour, 50 mM Na⁺, 250 nM primer</span></div>
      <div class="tscroll">${d.primers.length ? `<table class="grid"><thead><tr><th>Name</th><th>Sequence (5′→3′)</th><th>Length</th><th>Tm</th><th>GC</th><th>Binding sites</th><th></th></tr></thead><tbody>
      ${d.primers.map(p => { const b = A.psites[p.id] || []; return `<tr data-p="${p.id}"><td><b>${esc(p.name)}</b></td><td class="mono wrap">${esc(p.seq)}</td><td>${p.seq.length}</td><td>${DNA.tm(p.seq).toFixed(1)} °C</td><td>${DNA.gcPct(p.seq).toFixed(0)}%</td><td>${b.length ? b.map(x => `${fmt(x.start + 1)}..${fmt(x.end)} ${x.strand > 0 ? '→' : '←'}${x.full ? '' : ' (3′ 15 nt)'}`).join('<br>') : '<span style="color:var(--danger)">does not bind</span>'}</td><td><button class="ib" data-ed="${p.id}">Edit</button><button class="ib dng" data-del="${p.id}">Delete</button></td></tr>`; }).join('')}
      </tbody></table>` : `<div class="empty"><h2>No primers yet</h2><p>Select bases and choose Primer from Selection, or add a primer by sequence. Primers with 5′ tails are matched by their 3′ end.</p></div>`}</div></div>`;
    $('#padd').onclick = () => dlgPrimer(null);
    $('#ptop').onclick = () => primerFromSel(1);
    $('#pbot').onclick = () => primerFromSel(-1);
    const tb = v.querySelector('tbody'); if (!tb) return;
    tb.onclick = e => {
      const ed = e.target.closest('[data-ed]'); if (ed) { dlgPrimer(d.primers.find(p => p.id === +ed.dataset.ed)); return; }
      const del = e.target.closest('[data-del]'); if (del) { removePrimer(+del.dataset.del); return; }
      const tr = e.target.closest('tr[data-p]'); if (!tr) return;
      const b = (A.psites[+tr.dataset.p] || [])[0]; if (b) selectRange(b);
      tb.querySelectorAll('tr.on').forEach(r => r.classList.remove('on')); tr.classList.add('on');
    };
  }

  // ================================================================ HISTORY
  function renderHistory() {
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
    const m = $('#modal'), box = $('#dlg');
    box.className = 'dlg' + (wide ? ' wide' : '');
    box.innerHTML = `<div class="hd">${esc(title)}</div><div class="bd">${body}</div><div class="ft">${buttons.map((b, i) => `<button class="btn${b.primary ? ' pri' : ''}${b.danger ? ' dng' : ''}" data-b="${i}">${esc(b.label)}</button>`).join('')}</div>`;
    m.hidden = false;
    box.querySelector('.ft').onclick = e => {
      const b = e.target.closest('[data-b]'); if (!b) return;
      const r = buttons[+b.dataset.b].act ? buttons[+b.dataset.b].act() : true;
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
    const first = box.querySelector('.bd input, .bd textarea, .bd select'); if (first) setTimeout(() => { first.focus(); if (first.select) first.select(); }, 0);
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
    const rows = [['Ctrl+N', 'New DNA file'], ['Ctrl+O', 'Open file'], ['Ctrl+S / Ctrl+Shift+S', 'Save / Save As'], ['Ctrl+W', 'Close file'], ['Ctrl+Z / Ctrl+Y', 'Undo / Redo'], ['Ctrl+C / Ctrl+V', 'Copy / Paste bases'], ['Ctrl+A', 'Select all'], ['Ctrl+F', 'Find DNA or protein'], ['Ctrl+G', 'Go to position'], ['Ctrl+K', 'Add feature from selection'], ['Ctrl+1 … Ctrl+6', 'Map, Sequence, Enzymes, Features, Primers, History'], ['Ctrl + / Ctrl − / Ctrl+0', 'Zoom the map'], ['A C G T N', 'Type bases in the Sequence view'], ['Shift+click, Shift+arrows', 'Extend the selection'], ['Double-click a feature', 'Edit the feature'], ['Right-click', 'Context menu']];
    dialog({ title: 'Keyboard Shortcuts', body: `<table class="grid"><tbody>${rows.map(r => `<tr><td><b>${r[0]}</b></td><td>${r[1]}</td></tr>`).join('')}</tbody></table>`, buttons: [{ label: 'Close', primary: true }] });
  }
  function dlgAbout() {
    dialog({
      title: 'About GeneCanvas',
      body: `<p><b>GeneCanvas</b> ${host && host.version ? 'version ' + esc(host.version) : ''}<br>A desktop editor for plasmids and DNA sequences: maps, sequence editing, restriction analysis, features, primers and simulated gels.</p><p class="hint">Opens GenBank (.gb, .gbk), FASTA and .dna files. Saves GenBank and FASTA. The example plasmid is a synthetic sequence made for demonstration.</p>`,
      buttons: [{ label: 'Close', primary: true }]
    });
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
  function loadContent(name, data, path) {
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
        loadContent(f.name, /\.dna$/i.test(f.name) ? bin.buffer : new TextDecoder().decode(bin), f.path);
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
    if (/\.dna$/i.test(f.name)) r.readAsArrayBuffer(f); else r.readAsText(f);
  }
  function primerSitesFor(d) { const A = analyze(d); return p => A.psites[p.id] || []; }
  async function saveDoc(asNew) {
    const d = D(); if (!d) return false;
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
  async function exportSVG() {
    const d = D();
    if (S.view !== 'map') { S.view = 'map'; renderAll(); }
    const svg = $('#mapwrap svg'); if (!svg) return;
    const css = `<style>${[...document.styleSheets].flatMap(s => { try { return [...s.cssRules].map(r => r.cssText).filter(t => t.startsWith('svg')); } catch (e) { return []; } }).join('\n')} :root{--canvas:#fff;--bb:#2f3644;--faint:#8d95a3;--ink:#1c2230;--muted:#5e6776;--accent:#2563c9;--sel:rgba(0,0,0,0);--selline:transparent}</style>`;
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
    else if (t.id === 'side-enz') setOpt('enzMode', t.value);
    else if (t.id === 'side-orfmin') { const v = Math.max(20, Math.min(2000, parseInt(t.value, 10) || 100)); setOpt('orfMin', v); }
  });
  $('#side').addEventListener('click', e => { const b = e.target.closest('[data-mm]'); if (b) setOpt('mapMode', b.dataset.mm); });

  document.addEventListener('keydown', e => {
    if (!$('#modal').hidden) return;
    if (e.key === 'Escape') { closeMenu(); closeCtx(); return; }
    const mod = e.ctrlKey || e.metaKey;
    const inField = /^(INPUT|TEXTAREA|SELECT)$/.test(e.target.tagName);
    if (!mod) { if (!inField && (e.key === 'Delete' || e.key === 'Backspace') && D() && S.view === 'map' && selLen(D())) { e.preventDefault(); deleteSel(); } return; }
    const k = e.key.toLowerCase();
    const map = {
      n: dlgNew, o: openFiles, s: () => saveDoc(e.shiftKey), w: () => D() && closeDoc(S.active), f: openFind, g: () => D() && dlgGoto(), k: () => D() && dlgFeature(null),
      '1': () => setView('map'), '2': () => setView('seq'), '3': () => setView('enz'), '4': () => setView('feat'), '5': () => setView('primer'), '6': () => setView('hist'),
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
    if (host.onOpenFile) host.onOpenFile(f => { const bin = Uint8Array.from(atob(f.data), c => c.charCodeAt(0)); loadContent(f.name, /\.dna$/i.test(f.name) ? bin.buffer : new TextDecoder().decode(bin), f.path); });
    if (host.onCloseRequest) host.onCloseRequest(() => {
      const dirty = S.docs.filter(d => d.dirty);
      if (!dirty.length) { host.confirmClose(); return; }
      dialog({ title: 'Quit GeneCanvas', body: `<p>${dirty.length} file${dirty.length > 1 ? 's have' : ' has'} unsaved changes: ${dirty.map(d => '“' + esc(d.name) + '”').join(', ')}.</p><p>Quit without saving?</p>`, buttons: [{ label: 'Cancel', primary: true }, { label: 'Quit Without Saving', danger: true, act: () => host.confirmClose() }] });
    });
  }

  // ================================================================ start
  function start() {
    measureCW();
    renderMenuBar();
    addDoc(DNA.buildDemo());
    S.docs[0].dirty = false;
    renderAll();
    if (host && host.ready) host.ready();
  }
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(start); else start();
  window.__gc = { S, DNA, addDoc, setView, setSel, renderAll };
})();
