/* GeneCanvas — bootstrap for a tool window. Works out which page it is, asks the main
 * process for its payload, and gives the page the small host it expects. */
(function () {
  'use strict';
  const host = window.gcHost;
  const $ = s => document.querySelector(s);
  const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  let toastT = null;
  const toast = (msg, ms) => {
    const t = $('#toast'); if (!t) return;
    t.textContent = msg; t.hidden = false;
    clearTimeout(toastT); toastT = setTimeout(() => { t.hidden = true; }, ms || 2600);
  };
  const copy = (text, msg) => {
    navigator.clipboard.writeText(text).then(() => toast(msg || 'Copied')).catch(() => {
      const ta = document.createElement('textarea');
      ta.value = text; document.body.appendChild(ta); ta.select();
      try { document.execCommand('copy'); toast(msg || 'Copied'); } catch (e) { toast('Could not copy'); }
      ta.remove();
    });
  };
  // a very small context menu, enough for the Import / Export button
  function menu(anchor, items) {
    const c = $('#ctx');
    c.innerHTML = items.map((a, i) => `<div class="it" data-i="${i}">${esc(a.label)}</div>`).join('');
    c.hidden = false;
    const r = anchor.getBoundingClientRect();
    c.style.left = Math.min(r.left, window.innerWidth - 240) + 'px';
    c.style.top = r.bottom + 'px';
    c.onclick = e => { const it = e.target.closest('.it'); if (!it) return; c.hidden = true; items[+it.dataset.i].act(); };
  }
  document.addEventListener('mousedown', e => { if (!e.target.closest('#ctx')) $('#ctx').hidden = true; });

  const send = msg => { if (host && host.toolSend) host.toolSend(msg); else toast('This window is not connected to a sequence window.'); };

  let state = { kind: '', payload: {} };
  function build() {
    const root = $('#toolroot');
    const p = state.payload || {};
    const doc = p.docs && p.docs.length ? p.docs.find(d => d.id === p.current) || p.docs[0] : null;
    const ctx = {
      send, copy, toast, menu,
      docs: p.docs || [],
      doc,
      sites: doc ? DNA.findSites(doc.seq, doc.circular) : null,
      myFeatures: p.myFeatures || [],
      favourites: p.favourites || [],
      payload: p,
      saveFavourites: list => send({ kind: 'favourites', list })
    };
    document.title = ({ featuredb: 'Common Features', enzymedb: 'Restriction Enzymes', simulator: 'Simulation', cloning: 'Cloning' })[state.kind] || 'GeneCanvas';
    if (state.kind === 'featuredb') TOOLUI.featureBrowser(root, ctx);
    else if (state.kind === 'enzymedb') TOOLUI.enzymeBrowser(root, ctx);
    else if (state.kind === 'simulator') SIMUI.open(root, ctx);
    else if (state.kind === 'cloning') CLONEUI.open(root, ctx);
    else root.innerHTML = '<div class="hint">Nothing to show.</div>';
  }

  const kind = new URLSearchParams(location.search).get('kind') || '';
  if (host && host.toolInit) {
    host.toolInit(kind).then(r => { state = { kind: r.kind || kind, payload: r.payload || {} }; build(); });
    host.onToolReload(r => { state = { kind: r.kind || state.kind, payload: r.payload || {} }; build(); });
    host.onToolDocs(docs => { state.payload.docs = docs; build(); });
  } else {
    state.kind = kind;
    build();
  }
})();
