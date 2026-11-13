/* GeneCanvas — bootstrap for the separate cloning window.
 *
 * This window has no documents of its own. It is handed a snapshot of the open
 * sequences when it opens, is told again whenever they change in the main window, and
 * sends the finished construct back so the main window can open it as a new tab.
 */
(function () {
  'use strict';
  const host = window.gcHost || null;
  const $ = s => document.querySelector(s);
  const esc = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

  let DOCS = [];

  // ---- a toast and a modal, the same ones the main window uses
  let toastT = null;
  function toast(msg) {
    const t = $('#toast');
    t.textContent = msg; t.hidden = false;
    clearTimeout(toastT); toastT = setTimeout(() => { t.hidden = true; }, 2600);
  }
  let dlgKey = null;
  function dialog({ title, body, buttons, wide, onOpen }) {
    const m = $('#modal'), prev = $('#dlg'), box = prev.cloneNode(false);
    prev.parentNode.replaceChild(box, prev);
    box.className = 'dlg' + (wide ? ' wide' : '');
    box.innerHTML = `<div class="hd">${esc(title)}</div><div class="bd">${body}</div><div class="ft">${buttons.map((b, i) => `<button class="btn${b.primary ? ' pri' : ''}${b.danger ? ' dng' : ''}" data-b="${i}"${b.id ? ` id="${b.id}"` : ''}${b.disabled ? ' disabled' : ''}>${esc(b.label)}</button>`).join('')}</div>`;
    m.hidden = false;
    box.querySelector('.ft').onclick = e => {
      const b = e.target.closest('[data-b]'); if (!b) return;
      const r = buttons[+b.dataset.b].act ? buttons[+b.dataset.b].act() : true;
      if (r && typeof r.then === 'function') return;
      if (r !== false) closeDialog();
    };
    dlgKey = e => {
      if (e.key === 'Escape') { e.preventDefault(); closeDialog(); }
      else if (e.key === 'Enter' && e.target.tagName !== 'TEXTAREA' && !e.target.closest('.ft')) {
        const i = buttons.findIndex(b => b.primary);
        if (i >= 0) { e.preventDefault(); box.querySelector(`[data-b="${i}"]`).click(); }
      }
    };
    box.addEventListener('keydown', dlgKey);
    if (onOpen) onOpen(box);
    box.tabIndex = -1;
    const first = box.querySelector('.bd input, .bd textarea, .bd select');
    setTimeout(() => { if (first) { first.focus(); if (first.select) first.select(); } else box.focus(); }, 0);
  }
  function closeDialog() {
    const box = $('#dlg');
    if (dlgKey) box.removeEventListener('keydown', dlgKey);
    $('#modal').hidden = true; box.innerHTML = '';
  }
  window.showErr = msg => { const e = $('#dlg .err'); if (e) e.textContent = msg; return false; };

  function mount(payload) {
    DOCS = payload.docs || [];
    document.title = (payload.method === 'restriction' ? 'Restriction Cloning' : (payload.method === 'infusion' ? 'In-Fusion Cloning' : 'Gibson Assembly')) + ' — GeneCanvas';
    window.ASMUI.open(payload.method || 'gibson', {
      standalone: true,
      docs: () => DOCS,
      current: payload.current,
      addDoc: src => {
        if (host && host.assemblyCreate) host.assemblyCreate(src);
        toast(`Sent ${Number(src.seq.length).toLocaleString('en-US')} bp to the main window`);
      },
      toast, dialog,
      onClose: () => { if (host && host.closeSelf) host.closeSelf(); }
    });
  }

  if (host && host.assemblyInit) {
    host.assemblyInit().then(mount);
    if (host.onAssemblyDocs) host.onAssemblyDocs(docs => {
      DOCS = docs || [];
      // A fresh snapshot only matters once the user looks at the lists again, so the
      // window is left as it is rather than being redrawn under their hands.
      toast('The open sequences changed in the main window');
    });
  } else {
    document.body.innerHTML = '<div class="empty"><h2>Cloning</h2><p>This window is opened by GeneCanvas.</p></div>';
  }
})();
