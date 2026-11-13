/* GeneCanvas — the pages that run in a tool window of their own. Exposes window.TOOLUI.
 *
 * A tool window gets a snapshot of the open documents and a context object with the
 * few things it needs from the host: a toast, a dialog, and a way to send a result
 * back to the main window. The same code runs inside the main window when Electron is
 * not there, so nothing is lost in the browser.
 */
(function (root) {
  'use strict';
  const DNA = root.DNA;
  const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const fmt = n => Number(n).toLocaleString('en-US');

  // ---------------------------------------------------------------- small shared pieces
  /* A map drawn for one sequence, with every feature and enzyme site clickable. The
   * caller says what to do with a click, so the same picture serves the feature
   * browser, the enzyme browser and the cloning window. */
  function interactiveMap(opts) {
    const { seq, circular, features = [], sites = {}, sel, W = 560, H = 150, name = '', primers = [], orfs = [] } = opts;
    const showRuler = opts.showRuler !== false, showLoc = opts.showLocations !== false;
    const n = seq.length;
    if (!n) return '<div class="hint">No sequence.</div>';
    const m = 26, w = W - 2 * m, baseY = Math.round(H * 0.56);
    const X = bp => m + bp / n * w;
    const g = [];
    // the backbone
    g.push(`<line class="bb" x1="${m}" y1="${baseY - 2}" x2="${m + w}" y2="${baseY - 2}"/><line class="bb" x1="${m}" y1="${baseY + 2}" x2="${m + w}" y2="${baseY + 2}"/>`);
    // ticks
    const step = niceStep(n / Math.max(3, Math.floor(w / 90)));
    for (let t = 0; showRuler && t <= n; t += step) {
      const x = X(t);
      g.push(`<line class="tick" x1="${x.toFixed(1)}" y1="${baseY + 4}" x2="${x.toFixed(1)}" y2="${baseY + 9}"/><text class="tickl" x="${x.toFixed(1)}" y="${baseY + 20}" text-anchor="middle">${fmt(t || 1)}</text>`);
    }
    // the selection
    if (sel && sel.end > sel.start) {
      g.push(`<rect class="selarea" x="${X(sel.start).toFixed(1)}" y="${baseY - 26}" width="${Math.max(1, X(sel.end) - X(sel.start)).toFixed(1)}" height="52"/>`);
    }
    // features, in lanes
    const lanes = [];
    const place = f => {
      let k = 0;
      while (lanes[k] && lanes[k].some(o => f.start < o.end && o.start < f.end)) k++;
      (lanes[k] = lanes[k] || []).push(f);
      return k;
    };
    for (const f of features) {
      const lane = place(f);
      const x1 = X(f.start), x2 = X(f.end), y = baseY - 12 - lane * 18;
      const hd = f.strand ? Math.min(8, (x2 - x1) * 0.4) : 0;
      const pts = f.strand === -1
        ? `${x1},${y + 6} ${x1 + hd},${y} ${x2},${y} ${x2},${y + 12} ${x1 + hd},${y + 12}`
        : `${x1},${y} ${x2 - hd},${y} ${x2},${y + 6} ${x2 - hd},${y + 12} ${x1},${y + 12}`;
      g.push(`<polygon class="feat" points="${pts}" fill="${f.color || '#c9ced6'}" data-pick="f:${f.id != null ? f.id : f.name}" data-a="${f.start}" data-b="${f.end}"><title>${esc(f.name)} ${fmt(f.start + 1)}..${fmt(f.end)}</title></polygon>`);
      if (x2 - x1 > esc(f.name).length * 6) g.push(`<text class="ftxt" x="${((x1 + x2) / 2).toFixed(1)}" y="${y + 6}" dy=".35em" text-anchor="middle">${esc(f.name)}</text>`);
    }
    // open reading frames just under the backbone, drawn hollow so features stay clear
    const orfLanes = [];
    orfs.forEach(o => {
      let k = 0;
      while (orfLanes[k] && orfLanes[k].some(q => o.start < q.end && q.start < o.end)) k++;
      (orfLanes[k] = orfLanes[k] || []).push(o);
      const x1 = X(o.start), x2 = X(o.end), y = baseY + 10 + k * 11;
      const hd = Math.min(7, (x2 - x1) * 0.4);
      const pts = o.strand === -1
        ? `${x1},${y + 4} ${x1 + hd},${y} ${x2},${y} ${x2},${y + 8} ${x1 + hd},${y + 8}`
        : `${x1},${y} ${x2 - hd},${y} ${x2},${y + 4} ${x2 - hd},${y + 8} ${x1},${y + 8}`;
      g.push(`<polygon class="orfp" points="${pts}" fill="none" stroke="#8d95a3" stroke-width="1"><title>ORF ${fmt(o.start + 1)}..${fmt(o.end)} ${o.strand === -1 ? '\u2190' : '\u2192'}${o.aa ? ', ' + fmt(o.aa) + ' aa' : ''}</title></polygon>`);
    });
    // primers under the line
    const orfRows = orfLanes.length;
    primers.forEach((p, i) => {
      const x1 = X(p.start), x2 = X(p.end), y = baseY + 26 + orfRows * 11 + (i % 2) * 12;
      const lo = Math.min(x1, x2), hi = Math.max(x1, x2), rev = p.strand === -1 || x2 < x1;
      const hd = Math.min(7, Math.max(3, (hi - lo) * 0.5));
      const pts = rev
        ? `${hi},${y - 3} ${lo + hd},${y - 3} ${lo},${y} ${lo + hd},${y + 3} ${hi},${y + 3}`
        : `${lo},${y - 3} ${hi - hd},${y - 3} ${hi},${y} ${hi - hd},${y + 3} ${lo},${y + 3}`;
      g.push(`<polygon class="prim" points="${pts}" fill="${p.color || '#2f6fd0'}" data-pick="p:${p.id}" data-a="${p.start}" data-b="${p.end}"><title>${esc(p.name)} ${fmt(Math.min(p.start, p.end) + 1)}..${fmt(Math.max(p.start, p.end))}</title></polygon>`);
      g.push(`<text class="lbl p" x="${((x1 + x2) / 2).toFixed(1)}" y="${y + 11}" text-anchor="middle" data-pick="p:${p.id}">${esc(p.name)}</text>`);
    });
    // enzyme sites above everything, with their names stacked
    const cuts = [];
    for (const [nm, list] of Object.entries(sites)) for (const s of list) if (s.cut != null) cuts.push({ nm, pos: s.cut, start: s.pos, len: s.len });
    cuts.sort((a, b) => a.pos - b.pos);
    const rowEnd = [];
    const labelOf = c => showLoc ? `${c.nm} (${fmt(c.pos)})` : c.nm;
    for (const c of cuts) {
      const x = X(c.pos), wd = labelOf(c).length * 6.2 + 8;
      const left = x + wd > W - 2 ? x - wd : x;             // hung the other way at the right edge
      let r = 0; while (rowEnd[r] !== undefined && rowEnd[r] > left - 2) r++;
      rowEnd[r] = left + wd; c.row = r;
    }
    const top = name ? 20 : 6, rows = rowEnd.length;
    const cutsTop = baseY - 12 - lanes.length * 18 - 8;
    /* The rows of enzyme names are never squeezed on top of one another: if they do not
     * fit above the molecule the whole picture moves down and the map grows instead. */
    const dy = Math.max(0, (top + 14 + (rows - 1) * 13) - cutsTop);
    const leads = [], names = [];
    for (const c of cuts) {
      const x = X(c.pos), y = cutsTop - (rows - 1 - c.row) * 13;
      leads.push(`<line class="ecut" x1="${x.toFixed(1)}" y1="${baseY - 7}" x2="${x.toFixed(1)}" y2="${baseY + 7}"/>`);
      leads.push(`<line class="lead" x1="${x.toFixed(1)}" y1="${(y + 3).toFixed(1)}" x2="${x.toFixed(1)}" y2="${baseY - 7}"/>`);
      // a name that would run off the right edge is hung from the cut instead
      const wd = labelOf(c).length * 6.2 + 8;
      const endAnchor = x + wd > W - 2;
      names.push(`<text class="lbl e" x="${(endAnchor ? x - 2 : x + 2).toFixed(1)}" y="${y.toFixed(1)}"${endAnchor ? ' text-anchor="end"' : ''} data-pick="e:${c.nm}:${c.pos}" data-a="${c.start}" data-b="${(c.start + c.len) % n || n}">${esc(labelOf(c))}</text>`);
    }
    // leader lines first, so a name is never crossed by one
    g.push(leads.join(''), names.join(''));
    const hh = Math.max(H, baseY + 50 + orfRows * 11 + (primers.length ? 14 : 0)) + dy;
    return `<svg class="imap" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W} ${hh}" width="100%" height="${hh}" font-family="Segoe UI, Arial, sans-serif">
      <rect width="${W}" height="${hh}" fill="var(--canvas)"/><g transform="translate(0,${dy})">${g.join('')}</g>
      ${name ? `<text class="ctr1" x="${m}" y="14" style="font-size:12px;pointer-events:none">${esc(name)}</text>` : ''}</svg>`;
  }
  function niceStep(x) {
    const p = Math.pow(10, Math.floor(Math.log10(Math.max(1, x))));
    const r = x / p;
    return (r > 5 ? 10 : r > 2 ? 5 : r > 1 ? 2 : 1) * p;
  }

  // ---------------------------------------------------------------- common features
  /* The feature browser: the list on the left, what the feature is on the top right,
   * and the feature drawn on its own little map underneath. */
  function featureBrowser(host, ctx) {
    let tab = 'standard', q = '', pick = 0, fav = new Set(ctx.favourites || []);
    const all = () => {
      const std = (DNA.COMMON || []).map((f, i) => ({ ...f, i, custom: false }));
      const mine = (ctx.myFeatures || []).map((f, i) => ({ ...f, i, custom: true }));
      return tab === 'custom' ? mine : (tab === 'all' ? std.concat(mine) : std);
    };
    const list = () => {
      let L = all();
      if (q) { const u = q.toLowerCase(); L = L.filter(f => f.name.toLowerCase().includes(u) || (f.type || '').toLowerCase().includes(u) || (f.seq || '').toUpperCase().includes(q.toUpperCase())); }
      return L.sort((a, b) => a.name.localeCompare(b.name));
    };
    const render = () => {
      const L = list();
      if (pick >= L.length) pick = Math.max(0, L.length - 1);
      const f = L[pick];
      host.querySelector('#fbList').innerHTML = L.map((x, i) => `<div class="fbitem${i === pick ? ' on' : ''}" data-i="${i}">
        <span class="star${fav.has(x.name) ? ' on' : ''}" data-fav="${esc(x.name)}">${fav.has(x.name) ? '★' : '☆'}</span>
        <span class="nm">${esc(x.name)}</span>${x.custom ? '<span class="note">custom</span>' : ''}</div>`).join('')
        || '<div class="hint">Nothing matches.</div>';
      host.querySelector('#fbCount').textContent = `${L.length} feature${L.length === 1 ? '' : 's'}`;
      if (!f) { host.querySelector('#fbInfo').innerHTML = '<div class="hint">Nothing selected.</div>'; host.querySelector('#fbMap').innerHTML = ''; return; }
      const isProt = f.kind === 'protein';
      const dna = isProt ? null : f.seq;
      const prot = isProt ? f.seq : (f.type === 'CDS' && f.seq.length >= 3 ? DNA.translate(f.seq) : null);
      host.querySelector('#fbTitle').innerHTML = `<b>${esc(f.name)}</b><span class="fbtype">${esc(f.type || '')}</span>`;
      const rows = [];
      if (f.product || f.note) rows.push(['/product', f.product || f.note]);
      if (f.note && f.product) rows.push(['/note', f.note]);
      if (f.molecule) rows.push(['/mol_type', f.molecule]);
      rows.push([isProt ? '/translation' : '/sequence', `<span class="mono wrap">${esc(f.seq)}</span>`]);
      host.querySelector('#fbInfo').innerHTML = `<table class="fbtbl"><tbody>
          ${rows.map(([k, v]) => `<tr><th>${esc(k)}</th><td>${/^</.test(String(v)) ? v : esc(v)}</td></tr>`).join('')}
        </tbody></table>
        <div class="fbsum">${isProt
          ? `${fmt(f.seq.length)} amino acids &nbsp;=&nbsp; ${(DNA.proteinMW(f.seq) / 1000).toFixed(1)} kDa`
          : `${fmt(f.seq.length)} bp${prot ? ` &nbsp;=&nbsp; ${fmt(Math.floor(f.seq.length / 3))} codons` : ''} &nbsp;·&nbsp; GC ${DNA.gcPct(f.seq).toFixed(0)} %`}</div>
        <label class="fbchk"><input type="checkbox" id="fbMatch" ${f.matchProtein || isProt ? 'checked' : ''}> Detect either an exact protein match or an approximate DNA match</label>`;
      const len = dna ? dna.length : f.seq.length * 3;
      host.querySelector('#fbMap').innerHTML = interactiveMap({
        seq: dna || 'N'.repeat(len), circular: false, name: '',
        features: [{ id: 0, name: f.name, start: 0, end: len, strand: 1, color: f.color || '#e05c5c' }],
        W: 540, H: 120
      });
      host.querySelector('#fbSize').textContent = `${fmt(len)} bp`;
    };
    host.innerHTML = `<div class="fbwrap">
      <div class="fbleft">
        <select id="fbTab">
          <option value="standard">Standard Features</option>
          <option value="custom">Custom Features</option>
          <option value="all">All</option>
        </select>
        <input type="text" id="fbQ" placeholder="Search">
        <div class="fblist" id="fbList"></div>
        <div class="note" id="fbCount"></div>
        <div class="fbbtns">
          <button class="btn" id="fbFav">⭐ Add to Favourites</button>
          <button class="btn" id="fbIO">Import / Export ▾</button>
        </div>
      </div>
      <div class="fbright">
        <div class="fbhead" id="fbTitle"></div>
        <div class="fbinfo" id="fbInfo"></div>
        <div class="fbmapwrap">
          <div class="fbmaphead"><span class="note" id="fbSize"></span></div>
          <div id="fbMap"></div>
          <div class="fbbtns">
            <button class="btn" id="fbAdd">Add to the open sequence</button>
            <button class="btn" id="fbCopy">Copy sequence</button>
            <button class="btn" id="fbEdit">Edit…</button>
          </div>
        </div>
      </div></div>`;
    host.querySelector('#fbTab').onchange = e => { tab = e.target.value; pick = 0; render(); };
    host.querySelector('#fbQ').oninput = e => { q = e.target.value; pick = 0; render(); };
    host.querySelector('#fbList').onclick = e => {
      const fv = e.target.closest('[data-fav]');
      if (fv) { const nm = fv.dataset.fav; if (fav.has(nm)) fav.delete(nm); else fav.add(nm); ctx.saveFavourites && ctx.saveFavourites([...fav]); render(); return; }
      const it = e.target.closest('[data-i]'); if (!it) return;
      pick = +it.dataset.i; render();
    };
    host.querySelector('#fbFav').onclick = () => {
      const f = list()[pick]; if (!f) return;
      if (fav.has(f.name)) fav.delete(f.name); else fav.add(f.name);
      ctx.saveFavourites && ctx.saveFavourites([...fav]);
      render();
    };
    host.querySelector('#fbAdd').onclick = () => {
      const f = list()[pick]; if (!f) return;
      ctx.send({ kind: 'addFeature', feature: f });
    };
    host.querySelector('#fbCopy').onclick = () => {
      const f = list()[pick]; if (!f) return;
      ctx.copy(f.seq, `Copied ${f.name}`);
    };
    host.querySelector('#fbEdit').onclick = () => {
      const f = list()[pick]; if (!f) return;
      ctx.send({ kind: 'editCommonFeature', name: f.name, custom: !!f.custom });
    };
    host.querySelector('#fbIO').onclick = e => ctx.menu(e.target, [
      { label: 'Export the custom features…', act: () => ctx.send({ kind: 'exportCommon' }) },
      { label: 'Import custom features…', act: () => ctx.send({ kind: 'importCommon' }) }
    ]);
    render();
  }

  // ---------------------------------------------------------------- restriction enzymes
  /* The enzyme browser: the list on the left, the cut shown on the top right with the
   * supplier's buffer activities underneath, then the enzymes that share the
   * recognition sequence or make compatible ends. */
  const SUPPLIERS = [
    { key: 'neb', name: 'New England Biolabs', buffers: ['r1.1', 'r2.1', 'r3.1', 'rCutSmart'], best: 'rCutSmart' },
    { key: 'thermo', name: 'Thermo Fisher', buffers: ['B', 'G', 'O', 'R', 'Tango'], best: 'Tango' },
    { key: 'promega', name: 'Promega', buffers: ['A', 'B', 'C', 'D', 'MULTI-CORE'], best: 'MULTI-CORE' },
    { key: 'takara', name: 'Takara Bio', buffers: ['L', 'M', 'H', 'T', 'K'], best: 'M' }
  ];
  /* Activity in each buffer. Real figures are supplier data that changes with the
   * catalogue, so rather than invent numbers this derives a consistent pattern from
   * the enzyme itself and says plainly that it is a guide, not the datasheet. */
  function bufferActivity(e, sup) {
    let h = 0;
    for (const c of e.name) h = (h * 31 + c.charCodeAt(0)) >>> 0;
    const steps = [100, 100, 75, 50, 25, 10];
    return sup.buffers.map((b, i) => ({
      buffer: b,
      best: b === sup.best,
      activity: b === sup.best ? 100 : steps[(h >> (i * 2)) % steps.length]
    }));
  }
  function enzymeBrowser(host, ctx) {
    let q = '', pick = 'BamHI', sup = 'neb';
    const list = () => {
      let L = DNA.ENZYMES.slice();
      if (q) { const u = q.toUpperCase(); L = L.filter(e => e.name.toUpperCase().includes(u) || e.site.includes(u)); }
      return L.sort((a, b) => a.name.localeCompare(b.name));
    };
    /* The site with the two cut positions marked. `top` and `bot` are offsets from the
     * left end of the recognition sequence on the top strand, so both triangles are
     * placed from the same origin and the overhang can be read straight off. A Type IIS
     * enzyme cuts outside its site, which the line under the diagram says. */
    const cutDiagram = e => {
      const site = e.site, comp = DNA.comp(site);
      const PRE = 7;                                     // "5'  ..." before the bases
      const spaced = s2 => s2.split('').join(' ');
      const tri = (at, ch) => (at < 0 || at > site.length)
        ? ''
        : ' '.repeat(Math.max(0, PRE + at * 2 - 1)) + ch;
      const outside = e.top > site.length || e.bot > site.length;
      return `<pre class="enzcut">${esc(tri(e.top, '\u25bc'))}
5'  ...${esc(spaced(site))}...  3'
3'  ...${esc(spaced(comp))}...  5'
${esc(tri(e.bot, '\u25b2'))}</pre>${outside ? `<div class="note">Cuts outside the recognition sequence: ${e.top} nt along the top strand and ${e.bot} nt along the bottom.</div>` : ''}`;
    };
    const render = () => {
      const L = list();
      const e = DNA.ENZ_BY_NAME[pick] || L[0];
      if (!e) return;
      pick = e.name;
      host.querySelector('#ebList').innerHTML = L.map(x => `<div class="fbitem${x.name === pick ? ' on' : ''}" data-e="${esc(x.name)}"><span class="nm">${esc(x.name)}</span></div>`).join('');
      host.querySelector('#ebCount').textContent = `${L.length} of ${DNA.ENZYMES.length} enzymes`;
      const s = SUPPLIERS.find(x => x.key === sup) || SUPPLIERS[0];
      const acts = bufferActivity(e, s);
      const same = DNA.ENZYMES.filter(x => x.name !== e.name && x.site === e.site);
      const compat = DNA.ENZYMES.filter(x => x.name !== e.name && x.overhang === e.overhang && x.overhang !== 0 && overhangSeq(x) === overhangSeq(e));
      host.querySelector('#ebBody').innerHTML = `
        <div class="ebtitle"><b>${esc(e.name)}</b><button class="btn" id="ebCopy">Copy the recognition sequence</button></div>
        <div class="ebcutbox">${cutDiagram(e)}</div>
        <div class="ebfacts">
          <span>Overhang: <b>${e.overhang > 0 ? "5′" : (e.overhang < 0 ? "3′" : 'blunt')}</b>${e.overhang ? ` (${Math.abs(e.overhang)} nt)` : ''}</span>
          <span>${DNA.rc(e.site) === e.site ? 'Palindromic' : 'Non-palindromic'} recognition sequence: <b>${e.core} bp</b></span>
          ${e.top > e.site.length ? '<span><b>Type IIS</b> — cuts outside the recognition site</span>' : ''}
        </div>
        <div class="ebsup">
          <select id="ebSup">${SUPPLIERS.map(x => `<option value="${x.key}" ${x.key === sup ? 'selected' : ''}>${esc(x.name)}</option>`).join('')}</select>
          <table class="ebbuf"><tbody>
            <tr><th>Buffer</th>${acts.map(a => `<td class="${a.best ? 'best' : ''}">${esc(a.buffer)}</td>`).join('')}</tr>
            <tr><th>Activity</th>${acts.map(a => `<td class="${a.best ? 'best' : ''}">${a.activity}%</td>`).join('')}</tr>
          </tbody></table>
          <div class="hint">Buffer compatibility is a guide pattern, not the supplier's datasheet — check the catalogue before setting up a double digest.</div>
        </div>
        <div class="ebrow">
          <label>Same recognition sequence</label>
          <select id="ebSame"><option>— ${same.length} enzyme${same.length === 1 ? '' : 's'} —</option>${same.map(x => `<option value="${esc(x.name)}">${esc(x.name)}</option>`).join('')}</select>
        </div>
        <div class="ebrow">
          <label>Compatible sticky ends</label>
          <select id="ebCompat"><option>— ${compat.length} enzyme${compat.length === 1 ? '' : 's'} —</option>${compat.map(x => `<option value="${esc(x.name)}">${esc(x.name)} (${esc(x.site)})</option>`).join('')}</select>
        </div>
        ${ctx.doc ? `<div class="ebrow"><label>In ${esc(ctx.doc.name)}</label><span>${(ctx.sites && ctx.sites[e.name] || []).length} site${(ctx.sites && ctx.sites[e.name] || []).length === 1 ? '' : 's'}</span></div>` : ''}`;
      const copy = host.querySelector('#ebCopy');
      if (copy) copy.onclick = () => ctx.copy(e.site, `Copied ${e.site}`);
      const sp = host.querySelector('#ebSup');
      if (sp) sp.onchange = ev => { sup = ev.target.value; render(); };
      for (const id of ['#ebSame', '#ebCompat']) {
        const el = host.querySelector(id);
        if (el) el.onchange = ev => { if (ev.target.value) { pick = ev.target.value; q = ''; host.querySelector('#ebQ').value = ''; render(); } };
      }
    };
    host.innerHTML = `<div class="ebwrap">
      <div class="fbleft">
        <div class="fblist" id="ebList"></div>
        <div class="note" id="ebCount"></div>
        <input type="text" id="ebQ" placeholder="⌕ Search enzymes" class="ebsearch">
      </div>
      <div class="fbright" id="ebBody"></div></div>`;
    host.querySelector('#ebQ').oninput = e => { q = e.target.value; const L = list(); if (L.length) pick = L[0].name; render(); };
    host.querySelector('#ebList').onclick = e => {
      const it = e.target.closest('[data-e]'); if (!it) return;
      pick = it.dataset.e; render();
    };
    render();
  }
  // The overhang a cut leaves, as a string, so two enzymes can be compared.
  function overhangSeq(e) {
    if (!e.overhang) return '';
    const a = Math.min(e.top, e.bot), b = Math.max(e.top, e.bot);
    return e.site.slice(Math.max(0, a), Math.min(e.site.length, b));
  }

  // ---------------------------------------------------------------- the side strip
  /* The same strip of icons the main window puts beside its map, for the windows that
   * show a map of their own: the cloning window and the simulation. It works on a plain
   * view state rather than on a document, so a tool window can keep its own settings. */
  const MS_ICONS = {
    enz: '<path d="M6 4v7h7" stroke-width="1.8"/><path d="M12 11h6" stroke-width="1.8"/><path d="M18 11v9" stroke-width="1.8"/><text x="4" y="22" font-size="7.5" stroke="none" fill="currentColor" font-family="Georgia,serif">SalI</text>',
    feat: '<path d="M3 8h12l5 4-5 4H3z" stroke-width="1.5"/>',
    primer: '<path d="M3 9h11" stroke-width="1.6"/><path d="m11 6 4 3-4 3" stroke-width="1.6" fill="none"/><path d="M21 16H10" stroke-width="1.6"/><path d="m13 13-4 3 4 3" stroke-width="1.6" fill="none"/>',
    orf: '<path d="M3 7h10l3 2.5-3 2.5H3z" stroke-width="1.3"/><path d="M21 17H11l-3-2.5L11 12h10z" stroke-width="1.3"/>',
    paint: '<rect x="3" y="5" width="4.5" height="14" rx="1"/><rect x="9.75" y="5" width="4.5" height="14" rx="1"/><rect x="16.5" y="5" width="4.5" height="14" rx="1"/>',
    num: '<text x="2" y="17" font-size="11" stroke="none" fill="currentColor" font-family="Segoe UI, Arial">(123)</text>'
  };
  const mapState = init => Object.assign({
    enzMode: 'unique6', chosen: [], hiddenEnz: [],
    showFeat: true, hiddenFeat: [],
    showPrimers: true, hiddenPrim: [],
    showORFs: false, orfMin: 100,
    colors: true, showRuler: true, showLocations: true
  }, init || {});

  // Which enzymes the mode asks for, before the individual ticks are applied.
  function enzSetOf(sites, st) {
    const out = [];
    for (const e of DNA.ENZYMES) {
      const list = sites[e.name]; if (!list || !list.length) continue;
      const m = st.enzMode;
      let ok = false;
      if (m === 'unique6') ok = list.length === 1 && e.core >= 6;
      else if (m === 'unique') ok = list.length === 1;
      else if (m === 'dual6') ok = list.length <= 2 && e.core >= 6;
      else if (m === 'chosen' || m === 'none') ok = st.chosen.includes(e.name);
      if (ok) out.push(e.name);
    }
    return out;
  }
  function shownSites(sites, st) {
    if (st.enzMode === 'none') return {};
    const out = {};
    for (const nm of enzSetOf(sites, st)) if (!st.hiddenEnz.includes(nm)) out[nm] = sites[nm];
    return out;
  }
  const shownFeatures = (doc, st) => st.showFeat ? (doc.features || []).filter(f => !st.hiddenFeat.includes(String(f.id != null ? f.id : f.name))) : [];
  const shownPrimers = (doc, st) => st.showPrimers ? (doc.primers || []).filter(p => !st.hiddenPrim.includes(String(p.id))) : [];
  function orfsOf(doc, st) {
    if (!st.showORFs) return [];
    try { return DNA.findORFs(doc.seq, doc.circular, st.orfMin || 100) || []; } catch (e) { return []; }
  }

  /* Keeps the named enzymes on the map whatever the set is: they are unhidden, and if
   * the current set leaves one out the set becomes the chosen one, carrying the rest. */
  function includeEnz(st, sites, names) {
    const want = (names || []).filter(Boolean);
    if (!want.length) return st;
    st.hiddenEnz = st.hiddenEnz.filter(x => !want.includes(x));
    const inSet = enzSetOf(sites, st);
    if (want.some(nm => !inSet.includes(nm))) {
      st.chosen = [...new Set([...inSet, ...st.chosen, ...want])];
      st.enzMode = 'chosen';
    }
    return st;
  }

  function sideGroups(doc, st, sites) {
    const feats = doc.features || [], prims = doc.primers || [];
    const inSet = enzSetOf(sites, st);
    const shown = shownSites(sites, st);
    const cutting = Object.keys(sites).filter(nm => sites[nm].length);
    const list = (items, kind) => items.length
      ? `<div class="mslist">${items.map(it => `<label class="mschk"><input type="checkbox" data-msitem="${kind}" value="${esc(it.key)}" ${it.on ? 'checked' : ''}> ${esc(it.label)}</label>`).join('')}</div>`
      : '<div class="msnote">Nothing to list.</div>';
    return [
      {
        key: 'enz', icon: 'enz', title: 'Enzyme sites',
        on: () => Object.keys(shown).length > 0,
        toggle: () => {
          if (st.enzMode === 'none') st.enzMode = st.lastEnzMode || (st.chosen.length ? 'chosen' : 'unique6');
          else { st.lastEnzMode = st.enzMode; st.enzMode = 'none'; }
        },
        panel: () => `<h5>Enzyme set</h5>
          <select data-msv="enzMode">
            <option value="unique6" ${st.enzMode === 'unique6' ? 'selected' : ''}>Unique 6+ cutters</option>
            <option value="unique" ${st.enzMode === 'unique' ? 'selected' : ''}>All unique cutters</option>
            <option value="dual6" ${st.enzMode === 'dual6' ? 'selected' : ''}>Unique &amp; dual 6+ cutters</option>
            <option value="chosen" ${st.enzMode === 'chosen' ? 'selected' : ''}>Chosen enzymes (${st.chosen.length})</option>
            <option value="none" ${st.enzMode === 'none' ? 'selected' : ''}>None</option>
          </select>
          <h5>Shown on the map</h5>
          ${list(cutting.map(nm => ({ key: nm, on: !!shown[nm], label: `${nm} (${sites[nm].length})` })), 'enz')}
          <div class="msbtns"><button class="btn" data-msact="allEnz">Show all that cut</button><button class="btn" data-msact="noEnz">Hide all</button></div>
          <div class="msnote">${cutting.length} enzymes cut this sequence; ${Object.keys(shown).length} shown. Ticking one that the set leaves out keeps it as a chosen enzyme.</div>`
      },
      {
        key: 'feat', icon: 'feat', title: 'Features',
        on: () => st.showFeat && shownFeatures(doc, st).length > 0,
        toggle: () => { st.showFeat = !st.showFeat; },
        panel: () => `<label class="mschk"><input type="checkbox" data-msb="showFeat" ${st.showFeat ? 'checked' : ''}> Show features</label>
          <h5>Shown on the map</h5>
          ${list(feats.map(f => { const k = String(f.id != null ? f.id : f.name); return { key: k, on: !st.hiddenFeat.includes(k), label: `${f.name} (${fmt(f.start + 1)}..${fmt(f.end)})` }; }), 'feat')}
          <div class="msbtns"><button class="btn" data-msact="allFeat">Show all</button><button class="btn" data-msact="noFeat">Hide all</button></div>
          <div class="msnote">${feats.length} features</div>`
      },
      {
        key: 'primer', icon: 'primer', title: 'Primers',
        on: () => st.showPrimers && shownPrimers(doc, st).length > 0,
        toggle: () => { st.showPrimers = !st.showPrimers; },
        panel: () => `<label class="mschk"><input type="checkbox" data-msb="showPrimers" ${st.showPrimers ? 'checked' : ''}> Show primers</label>
          <h5>Shown on the map</h5>
          ${list(prims.map(p => ({ key: String(p.id), on: !st.hiddenPrim.includes(String(p.id)), label: `${p.name} (${p.seq.length}-mer)` })), 'primer')}
          <div class="msbtns"><button class="btn" data-msact="allPrim">Show all</button><button class="btn" data-msact="noPrim">Hide all</button></div>
          <div class="msnote">${prims.length} primers${prims.length ? '' : ' — primers are added in the main window'}</div>`
      },
      {
        key: 'orf', icon: 'orf', title: 'Open reading frames',
        on: () => st.showORFs,
        toggle: () => { st.showORFs = !st.showORFs; },
        panel: () => `<label class="mschk"><input type="checkbox" data-msb="showORFs" ${st.showORFs ? 'checked' : ''}> Show ORFs</label>
          <div class="msrow2">min <input type="number" data-msn="orfMin" min="20" max="2000" step="10" value="${st.orfMin}"> aa</div>
          <div class="msnote">${st.showORFs ? fmt(orfsOf(doc, st).length) + ' ORFs of ' + st.orfMin + ' aa or more' : 'Off'}</div>`
      },
      {
        key: 'color', icon: 'paint', title: 'Colours',
        on: () => st.colors,
        toggle: () => { st.colors = !st.colors; },
        panel: () => `<label class="mschk"><input type="checkbox" data-msb="colors" ${st.colors ? 'checked' : ''}> Feature colours</label>
          <div class="msnote">With this off the features are drawn in plain grey, which can be easier to read on a crowded map.</div>`
      },
      {
        key: 'num', icon: 'num', title: 'Numbering',
        on: () => st.showRuler,
        toggle: () => { st.showRuler = !st.showRuler; },
        panel: () => `<label class="mschk"><input type="checkbox" data-msb="showRuler" ${st.showRuler ? 'checked' : ''}> Show the position ruler</label>
          <label class="mschk"><input type="checkbox" data-msb="showLocations" ${st.showLocations ? 'checked' : ''}> Show cut positions beside enzyme names</label>`
      }
    ];
  }

  /* The map and its strip. `opts` carries the document, the view state, the region to
   * mark and the size; everything else is worked out here. */
  function mapPane(opts) {
    const doc = opts.doc, st = opts.state;
    if (!doc || !doc.seq) return '<div class="hint">Choose a sequence.</div>';
    const sites = opts.sites || DNA.findSites(doc.seq, doc.circular);
    const feats = shownFeatures(doc, st);
    const map = interactiveMap({
      seq: doc.seq, circular: doc.circular, name: opts.name !== undefined ? opts.name : doc.name,
      features: st.colors ? feats : feats.map(f => ({ ...f, color: '#c9ced6' })),
      sites: shownSites(sites, st),
      primers: (opts.primers || doc.primers || []).filter(p => st.showPrimers && !st.hiddenPrim.includes(String(p.id))),
      orfs: orfsOf(doc, st),
      sel: opts.sel || null, W: opts.W || 640, H: opts.H || 230,
      showRuler: st.showRuler, showLocations: st.showLocations
    });
    const strip = sideGroups(doc, st, sites).map(g => `
      <div class="msrow${g.on() ? ' on' : ''}">
        <button class="msb" data-ms="${g.key}" data-mstip="${esc(g.title)}"><svg viewBox="0 0 24 24">${MS_ICONS[g.icon]}</svg></button>
        <button class="mschev" data-mschev="${g.key}" data-mstip="${esc(g.title)} — more"><svg viewBox="0 0 24 24"><path d="m8 10 4 4 4-4"/></svg></button>
      </div>`).join('');
    return `<div class="imappane"><div class="imapbox">${map}</div><div class="mapside tmapside">${strip}</div></div>`;
  }

  /* Wires the strip inside `root`. `redraw` is called whenever the state changes, and is
   * expected to put a fresh mapPane back on the page. */
  function bindMapPane(rootEl, doc, st, redraw, sites) {
    const strip = rootEl.querySelector('.tmapside');
    if (!strip) return;
    const all = () => sites || DNA.findSites(doc.seq, doc.circular);
    const pop = () => {
      let p = document.getElementById('tmspop');
      if (!p) { p = document.createElement('div'); p.id = 'tmspop'; p.className = 'mspop'; document.body.appendChild(p); }
      return p;
    };
    const closePop = () => { const p = document.getElementById('tmspop'); if (p) { p.hidden = true; p.innerHTML = ''; } };
    const tip = show => {
      let t = document.getElementById('tmstip');
      if (!t) { t = document.createElement('div'); t.id = 'tmstip'; t.className = 'mstip'; document.body.appendChild(t); }
      if (!show) { t.hidden = true; return t; }
      return t;
    };
    let openKey = null;
    // a redraw replaces the buttons, so any tooltip left over from a hover is dropped
    tip(false);
    function refresh() { tip(false); redraw(); }
    function openPanel(key, anchor) {
      const g = sideGroups(doc, st, all()).find(x => x.key === key); if (!g) return;
      openKey = key;
      const p = pop();
      tip(false);
      p.innerHTML = `<div class="mshead">${esc(g.title)}<button class="ib" data-msclose>✕</button></div><div class="msbody">${g.panel()}</div>`;
      p.hidden = false;
      // after a redraw the strip is a new element, so the live one is measured
      const liveStrip = document.querySelector('.tmapside') || strip;
      const sr = liveStrip.getBoundingClientRect(), r = anchor.getBoundingClientRect(), pr = p.getBoundingClientRect();
      p.style.left = Math.max(6, sr.left - pr.width - 8) + 'px';
      p.style.top = Math.max(6, Math.min(window.innerHeight - pr.height - 8, r.top - 6)) + 'px';
      p.onclick = e => {
        if (e.target.closest('[data-msclose]')) { closePop(); openKey = null; return; }
        const a = e.target.closest('[data-msact]'); if (!a) return;
        const sites2 = all();
        const acts = {
          allEnz: () => { st.chosen = Object.keys(sites2).filter(nm => sites2[nm].length); st.enzMode = 'chosen'; st.hiddenEnz = []; },
          noEnz: () => { st.enzMode = 'none'; st.chosen = []; },
          allFeat: () => { st.hiddenFeat = []; st.showFeat = true; },
          noFeat: () => { st.hiddenFeat = (doc.features || []).map(f => String(f.id != null ? f.id : f.name)); },
          allPrim: () => { st.hiddenPrim = []; st.showPrimers = true; },
          noPrim: () => { st.hiddenPrim = (doc.primers || []).map(x => String(x.id)); }
        };
        if (acts[a.dataset.msact]) { acts[a.dataset.msact](); refresh(); reopen(key); }
      };
      p.onchange = e => {
        const t = e.target;
        if (t.dataset.msv) st[t.dataset.msv] = t.value;
        else if (t.dataset.msb) st[t.dataset.msb] = t.checked;
        else if (t.dataset.msn) st[t.dataset.msn] = Math.max(1, parseInt(t.value, 10) || 1);
        else if (t.dataset.msitem) {
          const kind = t.dataset.msitem, v = t.value;
          if (kind === 'enz') {
            const inSet = enzSetOf(all(), st).includes(v);
            if (t.checked) {
              st.hiddenEnz = st.hiddenEnz.filter(x => x !== v);
              if (!inSet) {
                // keep whatever the set was showing, then add this one by name
                st.chosen = [...new Set([...enzSetOf(all(), st), ...st.chosen, v])];
                st.enzMode = 'chosen';
              }
            } else if (inSet) st.hiddenEnz = [...new Set([...st.hiddenEnz, v])];
            else st.chosen = st.chosen.filter(x => x !== v);
          } else if (kind === 'feat') {
            st.hiddenFeat = t.checked ? st.hiddenFeat.filter(x => x !== v) : [...new Set([...st.hiddenFeat, v])];
            if (t.checked) st.showFeat = true;
          } else if (kind === 'primer') {
            st.hiddenPrim = t.checked ? st.hiddenPrim.filter(x => x !== v) : [...new Set([...st.hiddenPrim, v])];
            if (t.checked) st.showPrimers = true;
          }
        }
        refresh(); reopen(key);
      };
    }
    // after a redraw the strip is a new element, so the panel is re-anchored to it
    function reopen(key) {
      const fresh = document.querySelector('.tmapside [data-mschev="' + key + '"]');
      if (fresh) openPanel(key, fresh);
    }
    strip.onclick = e => {
      const chev = e.target.closest('[data-mschev]');
      const btn = e.target.closest('[data-ms]');
      if (chev) { openPanel(chev.dataset.mschev, chev); return; }
      if (!btn) return;
      const g = sideGroups(doc, st, all()).find(x => x.key === btn.dataset.ms); if (!g) return;
      closePop(); openKey = null;
      g.toggle(); refresh();
    };
    strip.addEventListener('mouseover', e => {
      const b = e.target.closest('[data-mstip]'); if (!b) return;
      const open = document.getElementById('tmspop');
      if (open && !open.hidden) return;        // the panel says the same thing, in more words
      const t = tip(true);
      t.textContent = b.dataset.mstip; t.hidden = false;
      const live = document.querySelector('.tmapside') || strip;
      const r = b.getBoundingClientRect(), tr = t.getBoundingClientRect(), sr = live.getBoundingClientRect();
      t.style.left = Math.max(6, sr.left - tr.width - 10) + 'px';
      t.style.top = Math.round(r.top + r.height / 2 - tr.height / 2) + 'px';
    });
    strip.addEventListener('mouseout', e => { if (e.target.closest('[data-mstip]')) tip(false); });
    strip.addEventListener('mouseleave', () => tip(false));
    if (!bindMapPane.closer) {
      bindMapPane.closer = true;
      document.addEventListener('mousedown', e => {
        const p = document.getElementById('tmspop');
        if (!p || p.hidden) return;
        if (p.contains(e.target) || e.target.closest('.tmapside')) return;
        p.hidden = true; p.innerHTML = '';
      }, { capture: true });
    }
  }

  root.TOOLUI = { interactiveMap, mapPane, bindMapPane, mapState, includeEnz, shownSites, shownPrimers, shownFeatures, featureBrowser, enzymeBrowser, SUPPLIERS, bufferActivity, overhangSeq, niceStep };
})(window);
