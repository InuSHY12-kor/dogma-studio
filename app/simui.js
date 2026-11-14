/* GeneCanvas — the Simulation window. Exposes window.SIMUI.
 *
 * A gel on the left and the sequence on the right. Each lane holds a ladder or a
 * sequence; a sequence lane is either cut with enzymes or amplified with a pair of
 * primers. Enzymes and primers are picked from the map itself: click one to choose it,
 * shift-click a second to make the pair or the double digest.
 */
(function (root) {
  'use strict';
  const DNA = root.DNA;
  const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const fmt = n => Number(n).toLocaleString('en-US');

  // ---------------------------------------------------------------- the gel picture
  /* Migration follows the resolving range of the gel: a band at the top of the range
   * barely moves, one at the bottom runs near the front. The model is the same one the
   * agarose document uses, so the two pictures agree. */
  const RANGES = [
    [0.5, 1000, 30000], [0.7, 800, 12000], [1.0, 500, 10000],
    [1.2, 400, 7000], [1.5, 200, 3000], [2.0, 50, 2000], [3.0, 20, 1000]
  ];
  function rangeFor(pct) {
    let best = RANGES[0];
    for (const r of RANGES) if (Math.abs(r[0] - pct) < Math.abs(best[0] - pct)) best = r;
    return [best[2], best[1]];
  }
  function migration(size, pct, run) {
    const [hi, lo] = rangeFor(pct);
    const t = (Math.log10(Math.max(1, size)) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo));
    const y = 1 - Math.max(0, Math.min(1, t));          // 0 at the well, 1 at the front
    return 0.06 + y * 0.88 * run;
  }
  const LADDERS = {
    '50 bp DNA Ladder': [50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1200, 1350],
    '100 bp DNA Ladder': [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517],
    '1 kb DNA Ladder': [500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000],
    '1 kb Plus DNA Ladder': [100, 200, 300, 400, 500, 700, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000],
    'Lambda HindIII': [564, 2027, 2322, 4361, 6557, 9416, 23130]
  };

  function open(hostEl, ctx) {
    const docs = (ctx.docs || []).filter(d => d.seq && d.seq.length);
    const state = {
      pct: 1.0, run: 0.85, laneCount: 10, ladder: '50 bp DNA Ladder',
      lanes: [], active: 1, tab: 'map'
    };
    for (let i = 0; i < 12; i++) state.lanes.push({ docId: docs[0] ? docs[0].id : null, mode: 'cut', enzymes: [], primers: [] });
    const docOf = lane => docs.find(d => d.id === lane.docId) || null;
    /* The strip beside the map keeps its settings per sequence, so switching lanes back
     * and forth does not lose what was hidden. */
    const viewStates = new Map();
    const viewState = doc => {
      const k = doc ? doc.id : 0;
      if (!viewStates.has(k)) viewStates.set(k, TOOLUI.mapState({ enzMode: 'dual6' }));
      return viewStates.get(k);
    };
    const sitesCache = new Map();
    const sitesOf = doc => {
      if (!doc) return {};
      if (!sitesCache.has(doc.id)) sitesCache.set(doc.id, DNA.findSites(doc.seq, doc.circular));
      return sitesCache.get(doc.id);
    };

    /* What a lane runs. A cut lane is a digest; a PCR lane is the amplicon between the
     * two chosen primers; with nothing chosen the whole molecule runs uncut. */
    function bandsFor(lane) {
      const doc = docOf(lane);
      if (!doc) return { bands: [], label: '', sub: '' };
      const n = doc.seq.length;
      if (lane.mode === 'pcr') {
        const ps = lane.primers.map(id => (doc.primers || []).find(p => p.id === id)).filter(Boolean);
        if (ps.length < 2) return { bands: [], label: doc.name, sub: 'pick two primers' };
        const hits = ps.map(p => ({ p, b: bindOf(doc, p) })).filter(x => x.b);
        if (hits.length < 2) return { bands: [], label: doc.name, sub: 'a primer does not bind' };
        const fwd = hits.find(x => x.b.strand === 1), rev = hits.find(x => x.b.strand === -1);
        if (!fwd || !rev) return { bands: [], label: doc.name, sub: 'both primers face the same way' };
        let len = rev.b.end - fwd.b.start;
        if (len <= 0) len = doc.circular ? len + n : 0;
        if (len <= 0) return { bands: [], label: doc.name, sub: 'the primers do not face each other' };
        return { bands: [{ size: len, amt: len }], label: doc.name, sub: `PCR ${esc(fwd.p.name)} + ${esc(rev.p.name)}`, amplicon: { start: fwd.b.start, end: rev.b.end, len } };
      }
      if (!lane.enzymes.length) {
        return { bands: [{ size: doc.circular ? Math.round(n * 0.7) : n, amt: n, sc: doc.circular }], label: doc.name, sub: doc.circular ? 'uncut, supercoiled' : 'uncut' };
      }
      const fr = DNA.digest(doc.seq, doc.circular, sitesOf(doc), lane.enzymes);
      return { bands: fr.map(f => ({ size: f.size, amt: f.size })), label: doc.name, sub: lane.enzymes.join(' + ') };
    }
    function bindOf(doc, p) {
      const hits = DNA.searchDNA(doc.seq, doc.circular, p.seq);
      if (hits.length) return hits[0];
      if (p.seq.length > 15) {
        const tail = p.seq.slice(-15);
        const h = DNA.searchDNA(doc.seq, doc.circular, tail);
        if (h.length) return h[0];
      }
      return null;
    }

    function gelSVG() {
      const W = 300, H = 440, left = 44, lw = (W - left - 10) / (state.laneCount + 1);
      const g = [];
      g.push(`<rect x="${left - 4}" y="18" width="${W - left + 2}" height="${H - 34}" fill="#161a20" rx="2"/>`);
      // wells
      for (let i = 0; i <= state.laneCount; i++) {
        const x = left + i * lw + 2;
        g.push(`<rect x="${x}" y="24" width="${lw - 5}" height="5" fill="#d8dde4"/>`);
        g.push(`<text x="${x + (lw - 5) / 2}" y="15" font-size="9" fill="#5e6776" text-anchor="middle">${i === 0 ? 'MW' : i}</text>`);
      }
      const bandY = frac => 30 + frac * (H - 70);
      // the ladder
      const lad = LADDERS[state.ladder] || [];
      for (const size of lad) {
        const y = bandY(migration(size, state.pct, state.run));
        g.push(`<rect x="${left + 2}" y="${y.toFixed(1)}" width="${lw - 5}" height="2.4" fill="#e8ecf2" opacity=".85"/>`);
        g.push(`<text x="${left - 7}" y="${(y + 3).toFixed(1)}" font-size="8.5" fill="#5e6776" text-anchor="end">${size}</text>`);
      }
      // the sample lanes
      for (let i = 1; i <= state.laneCount; i++) {
        const lane = state.lanes[i - 1];
        const r = bandsFor(lane);
        const x = left + i * lw + 2;
        if (i === state.active) g.push(`<rect x="${x - 2}" y="20" width="${lw - 1}" height="${H - 46}" fill="none" stroke="#5a8de0" stroke-width="1"/>`);
        for (const b of r.bands) {
          const y = bandY(migration(b.size, state.pct, state.run));
          const bright = Math.max(0.35, Math.min(1, Math.log10(Math.max(50, b.amt)) / 4));
          g.push(`<rect x="${x}" y="${y.toFixed(1)}" width="${lw - 5}" height="3" fill="#f2f6fb" opacity="${bright.toFixed(2)}"><title>${fmt(b.size)} bp</title></rect>`);
        }
      }
      return `<svg viewBox="0 0 ${W} ${H}" width="100%" height="${H}" font-family="Segoe UI, Arial, sans-serif">${g.join('')}</svg>`;
    }

    function laneListHTML() {
      const rows = [`<div class="simlane"><b>MW:</b> ${esc(state.ladder)}</div>`];
      for (let i = 1; i <= state.laneCount; i++) {
        const lane = state.lanes[i - 1], doc = docOf(lane), r = bandsFor(lane);
        rows.push(`<div class="simlane${i === state.active ? ' on' : ''}" data-lane="${i}">
          <b>${i}:</b> ${doc ? esc(doc.name) : '<span class="note">empty</span>'}
          ${doc ? `<div class="simsub">${esc(r.sub || '')}</div>
            <div class="simbands">${r.bands.map(b => `<span>${fmt(b.size)} bp</span>`).join('')}</div>` : ''}
        </div>`);
      }
      return rows.join('');
    }

    function rightHTML() {
      const lane = state.lanes[state.active - 1] || state.lanes[0];
      const doc = docOf(lane);
      const sites = sitesOf(doc);
      const r = bandsFor(lane);
      const chosen = lane.mode === 'cut'
        ? (lane.enzymes.length ? lane.enzymes.join(' + ') : 'nothing chosen')
        : (lane.primers.length ? lane.primers.map(id => (doc.primers || []).find(p => p.id === id)).filter(Boolean).map(p => p.name).join(' + ') : 'no primers chosen');
      return `<div class="simtop">
          <label>Lane ${state.active}:</label>
          <select id="simDoc">${docs.map(d => `<option value="${d.id}" ${doc && d.id === doc.id ? 'selected' : ''}>${esc(d.name)}</option>`).join('')}</select>
          <select id="simMode">
            <option value="cut" ${lane.mode === 'cut' ? 'selected' : ''}>cut with</option>
            <option value="pcr" ${lane.mode === 'pcr' ? 'selected' : ''}>amplify with</option>
          </select>
          <span class="note">${esc(chosen)}</span>
          <button class="btn" id="simClear">Clear</button>
          <button class="btn" id="simAll">Apply to All Lanes</button>
        </div>
        <div class="simsel">
          ${lane.mode === 'cut'
          ? `Click an enzyme on the map to cut with it, shift-click another for a double digest.`
          : `Click a primer on the map, shift-click the one facing it to amplify between them.`}
          ${r.amplicon ? ` &nbsp;·&nbsp; <b>${fmt(r.amplicon.len)} bp</b> from ${fmt(r.amplicon.start + 1)} to ${fmt(r.amplicon.end)}` : ''}
          ${lane.mode === 'cut' && r.bands.length ? ` &nbsp;·&nbsp; ${r.bands.length} fragment${r.bands.length === 1 ? '' : 's'}` : ''}
        </div>
        <div class="simtabs">
          ${['map', 'sequence', 'enzymes', 'primers'].map(t => `<button class="${state.tab === t ? 'on' : ''}" data-simtab="${t}">${t[0].toUpperCase() + t.slice(1)}</button>`).join('')}
        </div>
        <div class="simbody" id="simBody">${doc ? bodyHTML(doc, lane, sites, r) : '<div class="hint">No sequence in this lane.</div>'}</div>`;
    }
    function bodyHTML(doc, lane, sites, r) {
      if (state.tab === 'map') {
        const st = viewState(doc);
        // the enzymes this lane is cut with stay on the map whichever set is chosen
        TOOLUI.includeEnz(st, sites, lane.mode === 'cut' ? lane.enzymes : []);
        const primers = (doc.primers || []).map(p => {
          const b = bindOf(doc, p);
          return b ? { ...p, start: b.start, end: b.end, strand: b.strand, color: lane.primers.includes(p.id) ? '#c0392b' : '#7b3fb0' } : null;
        }).filter(Boolean);
        return TOOLUI.mapPane({
          doc, state: st, sites, primers,
          sel: r.amplicon || null, W: 700, H: 230
        });
      }
      if (state.tab === 'sequence') {
        const per = 60, out = [];
        for (let i = 0; i < Math.min(doc.seq.length, 6000); i += per) out.push(String(i + 1).padStart(7) + '  ' + doc.seq.slice(i, i + per));
        return `<pre class="simseq">${esc(out.join('\n'))}${doc.seq.length > 6000 ? '\n…' : ''}</pre>`;
      }
      if (state.tab === 'enzymes') {
        const rows = Object.entries(sites).sort((a, b) => a[0].localeCompare(b[0]));
        return `<table class="grid"><thead><tr><th></th><th>Enzyme</th><th>Site</th><th>Cuts</th><th>Positions</th></tr></thead><tbody>
          ${rows.map(([nm, list]) => `<tr data-enz="${esc(nm)}"${lane.enzymes.includes(nm) ? ' class="on"' : ''}>
            <td><input type="checkbox" data-enzck="${esc(nm)}" ${lane.enzymes.includes(nm) ? 'checked' : ''}></td>
            <td><b>${esc(nm)}</b></td><td class="mono">${esc((DNA.ENZ_BY_NAME[nm] || {}).site || '')}</td>
            <td>${list.length}</td><td class="wrap">${list.map(x => fmt(x.cut)).join(', ')}</td></tr>`).join('')}
        </tbody></table>`;
      }
      const ps = (doc.primers || []);
      return ps.length ? `<table class="grid"><thead><tr><th></th><th>Name</th><th>Sequence</th><th>Binds</th></tr></thead><tbody>
          ${ps.map(p => { const b = bindOf(doc, p); return `<tr data-prim="${p.id}"${lane.primers.includes(p.id) ? ' class="on"' : ''}>
            <td><input type="checkbox" data-primck="${p.id}" ${lane.primers.includes(p.id) ? 'checked' : ''}></td>
            <td><b>${esc(p.name)}</b></td><td class="mono wrap">${esc(p.seq)}</td>
            <td>${b ? `${fmt(b.start + 1)}..${fmt(b.end)} ${b.strand > 0 ? '→' : '←'}` : '<span style="color:var(--danger)">no site</span>'}</td></tr>`; }).join('')}
        </tbody></table>` : '<div class="hint">This sequence has no primers. Add some in the main window first.</div>';
    }

    function render() {
      hostEl.innerHTML = `<div class="simwrap">
        <div class="simleft">
          <div class="simgel" id="simGel"></div>
          <div class="simopts">
            <label>lanes <input type="number" id="simLanes" min="1" max="12" value="${state.laneCount}"></label>
            <label><input type="number" id="simPct" min="0.5" max="3" step="0.1" value="${state.pct}"> % agarose</label>
          </div>
          <div class="simopts">
            <label>run <input type="range" id="simRun" min="0.3" max="1" step="0.05" value="${state.run}"></label>
            <select id="simLadder">${Object.keys(LADDERS).map(k => `<option ${k === state.ladder ? 'selected' : ''}>${esc(k)}</option>`).join('')}</select>
          </div>
          <div class="simlanes" id="simLanes2">${laneListHTML()}</div>
          <div class="simopts"><button class="btn" id="simToGel">Open as an Agarose document</button></div>
        </div>
        <div class="simright" id="simRight">${rightHTML()}</div>
      </div>`;
      hostEl.querySelector('#simGel').innerHTML = gelSVG();
      bind();
    }
    function bindSide() {
      const lane = state.lanes[state.active - 1], doc = docOf(lane);
      if (!doc || state.tab !== 'map') return;
      TOOLUI.bindMapPane(hostEl, doc, viewState(doc), () => refresh(), sitesOf(doc), nm => ctx.send({ kind: 'openEnzyme', name: nm }));
    }
    function refresh() {
      hostEl.querySelector('#simGel').innerHTML = gelSVG();
      hostEl.querySelector('#simLanes2').innerHTML = laneListHTML();
      hostEl.querySelector('#simRight').innerHTML = rightHTML();
      bind();
    }
    function bind() {
      const q = s => hostEl.querySelector(s);
      const lane = state.lanes[state.active - 1];
      const doc = docOf(lane);
      const ln = q('#simLanes'); if (ln) ln.onchange = e => { state.laneCount = Math.max(1, Math.min(12, +e.target.value || 10)); if (state.active > state.laneCount) state.active = state.laneCount; render(); };
      const pc = q('#simPct'); if (pc) pc.onchange = e => { state.pct = Math.max(0.5, Math.min(3, parseFloat(e.target.value) || 1)); refresh(); };
      const rn = q('#simRun'); if (rn) rn.oninput = e => { state.run = parseFloat(e.target.value); q('#simGel').innerHTML = gelSVG(); };
      const ld = q('#simLadder'); if (ld) ld.onchange = e => { state.ladder = e.target.value; refresh(); };
      const ll = q('#simLanes2'); if (ll) ll.onclick = e => { const it = e.target.closest('[data-lane]'); if (!it) return; state.active = +it.dataset.lane; refresh(); };
      const dc = q('#simDoc'); if (dc) dc.onchange = e => { lane.docId = +e.target.value; lane.enzymes = []; lane.primers = []; refresh(); };
      const md = q('#simMode'); if (md) md.onchange = e => { lane.mode = e.target.value; refresh(); };
      const cl = q('#simClear'); if (cl) cl.onclick = () => { lane.enzymes = []; lane.primers = []; refresh(); };
      const al = q('#simAll'); if (al) al.onclick = () => {
        for (const other of state.lanes) { other.docId = lane.docId; other.mode = lane.mode; other.enzymes = lane.enzymes.slice(); other.primers = lane.primers.slice(); }
        refresh();
      };
      const tg = q('#simToGel'); if (tg) tg.onclick = () => {
        ctx.send({
          kind: 'gel', opts: {
            percent: state.pct, run: state.run, laneCount: state.laneCount, ladder: state.ladder,
            lanes: state.lanes.slice(0, state.laneCount).map(l => ({ docId: l.docId, enzymes: l.enzymes.slice() }))
          }
        });
      };
      hostEl.querySelectorAll('[data-simtab]').forEach(b => b.onclick = () => { state.tab = b.dataset.simtab; refresh(); });
      const body = q('#simBody');
      if (!body || !doc) return;
      bindSide();
      // picking on the map: click chooses, shift-click adds the second one
      body.addEventListener('click', e => {
        const pick = e.target.closest('[data-pick]');
        if (pick) {
          const id = pick.dataset.pick;
          if (id.startsWith('e:') && lane.mode === 'cut') {
            const nm = id.split(':')[1];
            lane.enzymes = e.shiftKey
              ? [...new Set([...lane.enzymes, nm])].slice(-4)
              : [nm];
            refresh(); return;
          }
          if (id.startsWith('p:') && lane.mode === 'pcr') {
            const pid = +id.split(':')[1];
            lane.primers = e.shiftKey ? [...new Set([...lane.primers, pid])].slice(-2) : [pid];
            refresh(); return;
          }
          if (id.startsWith('f:')) {
            const a = +pick.dataset.a, b = +pick.dataset.b;
            ctx.send({ kind: 'select', docId: doc.id, start: a, end: b });
            return;
          }
        }
        const ck = e.target.closest('[data-enzck]');
        if (ck) { const nm = ck.dataset.enzck; lane.enzymes = ck.checked ? [...new Set([...lane.enzymes, nm])] : lane.enzymes.filter(x => x !== nm); refresh(); return; }
        const pk = e.target.closest('[data-primck]');
        if (pk) { const pid = +pk.dataset.primck; lane.primers = pk.checked ? [...new Set([...lane.primers, pid])].slice(-2) : lane.primers.filter(x => x !== pid); refresh(); return; }
        const row = e.target.closest('[data-enz]');
        if (row && !e.target.closest('input')) { const nm = row.dataset.enz; lane.enzymes = e.shiftKey ? [...new Set([...lane.enzymes, nm])] : [nm]; refresh(); return; }
        const prow = e.target.closest('[data-prim]');
        if (prow && !e.target.closest('input')) { const pid = +prow.dataset.prim; lane.primers = e.shiftKey ? [...new Set([...lane.primers, pid])].slice(-2) : [pid]; refresh(); }
      });
    }
    render();
  }

  root.SIMUI = { open, migration, rangeFor, LADDERS };
})(window);
