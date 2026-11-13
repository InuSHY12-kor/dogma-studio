/* GeneCanvas — the Vector / Fragment / Product cloning window.
 * Exposes window.ASMUI.open(method, ctx)
 *
 * The same window runs inside the main window and, on the desktop, in a window of its own.
 * Everything it needs from its host arrives in `ctx`, so it does not reach into the
 * application state directly:
 *   ctx.docs()      the nucleic-acid documents it may use
 *   ctx.current     the id of the document to start from, if any
 *   ctx.addDoc(src) make the finished construct
 *   ctx.toast(msg)  a short message
 *   ctx.dialog(o)   a modal, used by Choose Overlapping PCR Primers
 *   ctx.onClose()   called when the window is dismissed
 */
(function (root) {
  'use strict';
  const DNA = root.DNA, ASSEMBLY = root.ASSEMBLY;

  // ================================================================ cloning window
  // One window with Vector, Fragment and Product tabs, shared by every assembly method.
  const val = id => (document.getElementById(id) || {}).value;
  const $ = s => document.querySelector(s);
  const esc = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const fmt = n => Number(n).toLocaleString('en-US');
  const TAU = Math.PI * 2;
  const ASM_METHODS = {
    gibson: { title: 'Gibson Assembly: Insert Fragment', ready: 'Ready for Gibson Assembly', joins: 'overlap' },
    infusion: { title: 'In-Fusion Cloning: Insert Fragment', ready: 'Ready for In-Fusion Cloning', joins: 'overlap' },
    restriction: { title: 'Restriction Cloning: Insert Fragment', ready: 'Ready for Ligation', joins: 'ends' }
  };
  const FRAG_COLORS = ['#2f6fd0', '#1a9e5c', '#c9731b', '#8a4ec4', '#c43c3c', '#1b9aa8'];

  function openAssembly(method, ctx) {
    const M = ASM_METHODS[method] || ASM_METHODS.gibson;
    const docs = () => ctx.docs();
    if (docs().length < 1) { ctx.toast('Open a sequence first.'); return; }
    const first = ctx.current && docs().some(x => String(x.id) === String(ctx.current)) ? docs().find(x => String(x.id) === String(ctx.current)) : docs()[0];
    const A = {
      method, tab: 'vector', activeFrag: 0,
      vector: { docId: first.id, prep: first.circular ? 'enzymes' : 'direct', cuts: [{}, {}], keepInner: false, from: 0, to: first.seq.length, orient: 1 },
      fragments: [{ docId: (docs()[1] || first).id, prep: 'pcr', cuts: [{}, {}], from: 0, to: 0, orient: 1 }],
      opts: { minOverlap: 15, maxOverlap: 25, fixed: false, fixedLength: 40, overlapTm: 50, primerTm: 60, regenUp: false, regenDown: false },
      product: { name: 'Assembled', close: true, strain: '' },
      state: null
    };
    const f0 = docs()[1] || first;
    const firstCds = f0.features.find(x => x.type === 'CDS');
    if (firstCds) { A.fragments[0].from = firstCds.start; A.fragments[0].to = firstCds.end; }
    else { A.fragments[0].from = 0; A.fragments[0].to = Math.min(f0.seq.length, 800); }

    const docOf = id => docs().find(x => String(x.id) === String(id));

    // ---- build the ordered pieces from the current settings
    function compute() {
      const errs = [];
      let pieces = null, plan = null, built = null, frames = [], suggest = null;
      try {
        const vdoc = docOf(A.vector.docId);
        if (!vdoc) throw new Error('Choose a vector.');
        const V = ASSEMBLY.prepare(vdoc, specOf(A.vector, vdoc));
        V.name = 'Linearized Vector'; V.spec = specOf(A.vector, vdoc); V.role = 'vector';
        pieces = [A.vector.orient === -1 ? Object.assign(ASSEMBLY.flip(V), { name: V.name, spec: V.spec, role: 'vector', doc: V.doc, from: V.from, to: V.to, label: V.label }) : V];
        A.fragments.forEach((f, i) => {
          const fd = docOf(f.docId);
          if (!fd) throw new Error(`Choose a source for fragment ${i + 1}.`);
          const P = ASSEMBLY.prepare(fd, specOf(f, fd));
          P.name = A.fragments.length > 1 ? `Fragment ${i + 1}` : 'Fragment';
          P.spec = specOf(f, fd); P.role = 'fragment'; P.index = i;
          pieces.push(f.orient === -1 ? Object.assign(ASSEMBLY.flip(P), { name: P.name, spec: P.spec, role: 'fragment', index: i, doc: P.doc, from: P.from, to: P.to, label: P.label }) : P);
        });
        const planner = ps => M.joins === 'overlap' ? ASSEMBLY.planGibson(ps, A.opts) : ASSEMBLY.planRestriction(ps);
        plan = planner(pieces);
        built = ASSEMBLY.build(pieces, plan, A.product.name);
        frames = ASSEMBLY.frameIssues(built, built.bounds);
        if (plan.junctions.some(j => !j.ok)) suggest = ASSEMBLY.suggestFlips(pieces, planner);
      } catch (e) { errs.push(e.message); }
      A.state = { pieces, plan, built, frames, errs, suggest };
      return A.state;
    }
    function specOf(sel, doc) {
      if (sel.prep === 'enzymes') {
        const cuts = (sel.cuts || []).filter(c => c && c.enzyme).map(c => {
          const list = (DNA.findSites(doc.seq, doc.circular)[c.enzyme] || []).filter(x => x.cut != null);
          const s = list[Math.min(c.site || 0, list.length - 1)];
          return s ? { enzyme: c.enzyme, cut: s.cut, cut2: s.cut2 } : null;
        }).filter(Boolean);
        return { prep: 'enzymes', cuts, keepInner: !!sel.keepInner };
      }
      if (sel.prep === 'pcr') return { prep: 'pcr', from: sel.from | 0, to: sel.to | 0 };
      return { prep: 'direct' };
    }

    // ---- a small map with the chosen region picked out
    function asmMap(doc, region, w, h) {
      if (!doc) return '';
      const n = doc.seq.length;
      const feats = doc.features;
      if (doc.circular) {
        const cx = w / 2, cy = h / 2, R = Math.min(w, h) / 2 - 26;
        const P = (r, a) => `${(cx + r * Math.cos(a - Math.PI / 2)).toFixed(1)},${(cy + r * Math.sin(a - Math.PI / 2)).toFixed(1)}`;
        const arc = (f, rr, ww, col) => {
          const a1 = f.start / n * TAU, span = Math.max(0.02, DNA.rlen(f, n) / n * TAU), a2 = a1 + span, lg = span > Math.PI ? 1 : 0;
          return `<path d="M${P(rr + ww, a1)} A${rr + ww},${rr + ww} 0 ${lg} 1 ${P(rr + ww, a2)} L${P(rr - ww, a2)} A${rr - ww},${rr - ww} 0 ${lg} 0 ${P(rr - ww, a1)} Z" fill="${col}" stroke="rgba(0,0,0,.25)" stroke-width=".7"><title>${esc(f.name)}</title></path>`;
        };
        let g = `<circle cx="${cx}" cy="${cy}" r="${R}" fill="none" stroke="#2f3644" stroke-width="2"/>`;
        if (region && region.mark !== undefined) {
          const a = region.mark / n * TAU;
          g += `<line x1="${P(R - 12, a).replace(',', '" y1="')}" x2="${P(R + 12, a).replace(',', '" y2="')}" stroke="#c43c3c" stroke-width="2.4"/>`;
          const q = P(R + 20, a).split(',');
          g += `<text x="${q[0]}" y="${q[1]}" font-size="10.5" fill="#c43c3c" text-anchor="middle" font-weight="700">${esc(region.label || '')}</text>`;
        }
        if (region && region.start !== undefined) g += arc({ start: region.start, end: region.end, name: 'selected' }, R, 9, 'rgba(37,99,201,.35)');
        g += feats.map(f => arc(f, R, 5, f.color || '#c9ced6')).join('');
        g += `<text x="${cx}" y="${cy - 2}" text-anchor="middle" font-size="12.5" font-weight="700" fill="#1c2230">${esc(doc.name)}</text>`;
        g += `<text x="${cx}" y="${cy + 14}" text-anchor="middle" font-size="11" fill="#5e6776">${fmt(n)} bp</text>`;
        return `<svg width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" font-family="Segoe UI, Arial, sans-serif">${g}</svg>`;
      }
      const m = 18, ww = w - 2 * m, y = h / 2;
      const X = p => m + p / n * ww;
      let g = `<line x1="${m}" y1="${y}" x2="${m + ww}" y2="${y}" stroke="#2f3644" stroke-width="2"/>`;
      if (region && region.start !== undefined) g += `<rect x="${X(region.start).toFixed(1)}" y="${y - 12}" width="${Math.max(2, X(region.end) - X(region.start)).toFixed(1)}" height="24" fill="rgba(37,99,201,.3)"/>`;
      g += feats.map(f => `<rect x="${X(f.start).toFixed(1)}" y="${y - 6}" width="${Math.max(1.5, X(f.end) - X(f.start)).toFixed(1)}" height="12" fill="${f.color || '#c9ced6'}" stroke="rgba(0,0,0,.25)" stroke-width=".7"><title>${esc(f.name)}</title></rect>`).join('');
      g += `<text x="${w / 2}" y="18" text-anchor="middle" font-size="12.5" font-weight="700">${esc(doc.name)}</text>`;
      g += `<text x="${w / 2}" y="${h - 8}" text-anchor="middle" font-size="11" fill="#5e6776">${fmt(n)} bp</text>`;
      return `<svg width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" font-family="Segoe UI, Arial, sans-serif">${g}</svg>`;
    }

    // ---- the picture of the finished circle along the bottom of the Product tab
    function junctionDiagram(st, w) {
      if (!st.pieces || !st.plan) return '';
      const pieces = st.pieces, plan = st.plan;
      const total = pieces.reduce((a, p) => a + p.seq.length, 0) || 1;
      const m = 20, bw = w - 2 * m, barY = 76, barH = 26;
      let x = m, g = '';
      const centres = [];
      pieces.forEach((p, i) => {
        const pw = Math.max(34, p.seq.length / total * bw);
        const col = p.role === 'vector' ? '#1c2230' : FRAG_COLORS[(p.index || 0) % FRAG_COLORS.length];
        g += `<rect x="${x.toFixed(1)}" y="${barY}" width="${pw.toFixed(1)}" height="${barH}" rx="3" fill="${col}"/>`;
        const label = `${p.name} · ${fmt(p.seq.length)} bp`;
        if (pw > label.length * 5.6) g += `<text x="${(x + pw / 2).toFixed(1)}" y="${barY + 17}" text-anchor="middle" font-size="11" fill="#fff" font-weight="600">${esc(label)}</text>`;
        else g += `<text x="${(x + pw / 2).toFixed(1)}" y="${barY + barH + 13}" text-anchor="middle" font-size="10" fill="#5e6776">${esc(p.name)}</text>`;
        centres.push({ x: x + pw / 2, w: pw, piece: p, i });
        x += pw;
      });
      // the vector wraps the whole circle, drawn as a bracket under everything
      g += `<path d="M${m},${barY - 46} L${m},${barY - 54} L${(m + bw).toFixed(1)},${barY - 54} L${(m + bw).toFixed(1)},${barY - 46}" fill="none" stroke="#1c2230" stroke-width="1.6"/>`;
      g += `<text x="${(m + bw / 2).toFixed(1)}" y="${barY - 58}" text-anchor="middle" font-size="10.5" fill="#1c2230">circular product · ${fmt(total)} bp</text>`;
      // Junction labels sit on the boundaries, dropping to a second line where two of
      // them would otherwise print on top of one another.
      let jx = m;
      const taken = [];
      pieces.forEach((p, i) => {
        jx += Math.max(34, p.seq.length / total * bw);
        const j = plan.junctions[i];
        if (!j) return;
        const atEnd = i === pieces.length - 1;
        const px = atEnd ? m + bw : jx;
        g += `<line x1="${px.toFixed(1)}" y1="${barY - 6}" x2="${px.toFixed(1)}" y2="${barY + barH + 6}" stroke="${j.ok ? '#1a7f37' : '#c43c3c'}" stroke-width="1.4" stroke-dasharray="${j.ok ? '' : '3 2'}"/>`;
        const txt = M.joins === 'overlap'
          ? (j.ok ? `${j.len} bp · Tm ${j.tm.toFixed(0)}°C` : 'no overlap')
          : (j.ok ? (j.overhang || 'blunt') : 'ends differ');
        const wTxt = txt.length * 5.4;
        const x0 = atEnd ? px - wTxt : px - wTxt / 2, x1 = atEnd ? px : px + wTxt / 2;
        let row = 0;
        while (taken.some(t => t.row === row && x0 < t.x1 + 6 && t.x0 < x1 + 6)) row++;
        taken.push({ row, x0, x1 });
        g += `<text x="${px.toFixed(1)}" y="${(barY - 10 - row * 12).toFixed(1)}" text-anchor="${atEnd ? 'end' : 'middle'}" font-size="10" fill="${j.ok ? '#1a7f37' : '#c43c3c'}" font-weight="600">${esc(txt)}</text>`;
      });
      // an arrow under each fragment showing which way round it goes
      for (const c of centres) {
        if (c.piece.role !== 'fragment') continue;
        const i = c.piece.index;
        const rev = A.fragments[i].orient === -1;
        const ax = c.x, ay = barY + barH + 22;
        const half = Math.min(28, c.w / 2 - 4);
        if (half < 8) continue;
        g += rev
          ? `<path d="M${(ax + half).toFixed(1)},${ay} L${(ax - half + 7).toFixed(1)},${ay} M${(ax - half).toFixed(1)},${ay} l7,-4 v8 z" stroke="#5e6776" fill="#5e6776" stroke-width="1.3"/>`
          : `<path d="M${(ax - half).toFixed(1)},${ay} L${(ax + half - 7).toFixed(1)},${ay} M${(ax + half).toFixed(1)},${ay} l-7,-4 v8 z" stroke="#5e6776" fill="#5e6776" stroke-width="1.3"/>`;
        g += `<text x="${ax.toFixed(1)}" y="${ay + 15}" text-anchor="middle" font-size="9.5" fill="#8d95a3">${rev ? 'reverse' : 'forward'}</text>`;
      }
      return `<svg width="${w}" height="${barY + barH + 54}" viewBox="0 0 ${w} ${barY + barH + 54}" font-family="Segoe UI, Arial, sans-serif">${g}</svg>`;
    }

    // ---- status box, bottom right
    function statusHTML(st) {
      if (st.errs.length) return `<div class="asmstat warn">${st.errs.map(esc).join('<br>')}</div>`;
      const bad = st.plan.junctions.filter(j => !j.ok);
      if (bad.length) {
        const hint = st.suggest && st.suggest.flip.length
          ? `<br><b>Try turning ${esc(st.suggest.flip.join(' and '))} round</b> — the ends fit that way.`
          : '';
        return `<div class="asmstat warn"><b>Not ready yet.</b><br>${bad.map(j => esc(j.why)).join('<br>')}${hint}</div>`;
      }
      const need = st.pieces.filter((p, i) => p.spec.prep === 'pcr' && !st.plan.primers[i]);
      if (need.length) return `<div class="asmstat warn">Please specify primers for ${need.map(p => esc(p.name)).join(', ')}.</div>`;
      const fr = st.frames;
      const warn = fr.length
        ? `<div class="asmwarn" title="${esc(fr.map(f => `${f.a.name} and ${f.b.name} run together ${f.shift} base${f.shift > 1 ? 's' : ''} out of frame`).join('\n'))}">⚠ ${fr.length} reading-frame problem${fr.length > 1 ? 's' : ''} — hover to see which</div>`
        : '';
      return `<div class="asmstat ok"><b>${esc(M.ready)}</b><br>Product: ${fmt(st.built.seq.length)} bp${warn}</div>`;
    }

    // ---- primers panel, bottom left of the Vector and Fragment tabs
    function primerPanel(st, pieceIndex) {
      const p = st.pieces && st.pieces[pieceIndex];
      const pr = st.plan && st.plan.primers[pieceIndex];
      if (!p) return '<div class="asmprimers"><div class="hint">Choose a source first.</div></div>';
      if (p.spec.prep !== 'pcr') return `<div class="asmprimers"><div class="hint">${esc(p.name)} is ${p.spec.prep === 'direct' ? 'used as it is' : 'cut with enzymes'}, so it needs no primers.</div></div>`;
      const row = (n, q) => `<div class="asmprow">
          <label>Primer ${n}:</label>
          <input type="text" class="asmpname" value="${q ? esc(q.name) : ''}" data-pname="${pieceIndex}:${n}" placeholder="&lt;${esc(p.name)}.${n === 1 ? 'FOR' : 'REV'}&gt;">
          <span class="note">${q ? `${q.seq.length}-mer / Tm = ${q.tm.toFixed(0)}°C` : ''}</span>
        </div>
        <div class="asmpseq"><span class="e5">5′</span><span class="asmpbox">${q ? `<span class="tail">${esc(q.tail)}</span>${esc(q.anneal)}` : ''}</span><span class="e5">3′</span></div>`;
      return `<div class="asmprimers">${row(1, pr && pr.fwd)}${row(2, pr && pr.rev)}
        <div class="hint">The red part is the tail that creates the overlap; the rest anneals to the template.</div></div>`;
    }

    // ---- right-hand controls per tab
    function vectorPanel() {
      const v = A.vector, doc = docOf(v.docId);
      const cutters = doc ? DNA.ENZYMES.filter(e => (DNA.findSites(doc.seq, doc.circular)[e.name] || []).some(x => x.cut != null)) : [];
      const enzSel = (idx) => {
        const cur = (v.cuts[idx] || {}).enzyme || '';
        const list = cur && doc ? (DNA.findSites(doc.seq, doc.circular)[cur] || []).filter(x => x.cut != null) : [];
        const siteSel = `<select data-vsite="${idx}" ${list.length > 1 ? '' : 'disabled'}>${list.map((x, k) => `<option value="${k}" ${k === ((v.cuts[idx] || {}).site || 0) ? 'selected' : ''}>${fmt(x.cut)}</option>`).join('') || '<option>—</option>'}</select>`;
        return `<div class="asmrow"><span>Cut</span>${siteSel}<span>with</span>
          <select data-venz="${idx}"><option value="">—</option>${cutters.map(e => `<option ${e.name === cur ? 'selected' : ''}>${e.name}</option>`).join('')}</select></div>`;
      };
      const st = A.state;
      const sites = st && st.pieces && st.pieces[0] ? [st.pieces[0].label] : [];
      return `<div class="asmside">
        <label class="asmlbl">Vector:</label>
        <select id="asmVdoc">${docs().map(x => `<option value="${x.id}" ${String(x.id) === String(v.docId) ? 'selected' : ''}>${esc(x.name)}</option>`).join('')}</select>
        <label class="asmradio"><input type="radio" name="vprep" value="pcr" ${v.prep === 'pcr' ? 'checked' : ''}> Linearize by PCR</label>
        <label class="asmradio"><input type="radio" name="vprep" value="enzymes" ${v.prep === 'enzymes' ? 'checked' : ''} ${doc && doc.circular ? '' : ''}> Linearize with restriction enzymes</label>
        <label class="asmradio${doc && doc.circular ? ' dis' : ''}"><input type="radio" name="vprep" value="direct" ${v.prep === 'direct' ? 'checked' : ''} ${doc && doc.circular ? 'disabled' : ''}> Use directly as the linearized Vector</label>
        ${v.prep === 'enzymes' ? enzSel(0) + enzSel(1) : ''}
        ${v.prep === 'enzymes' && (v.cuts[0] || {}).enzyme && (v.cuts[1] || {}).enzyme ? `<label class="asmradio"><input type="checkbox" id="asmKeepInner" ${v.keepInner ? 'checked' : ''}> Keep the stretch between the cuts instead</label>` : ''}
        ${v.prep === 'pcr' ? `<div class="asmrow"><span>Keep</span><input type="number" id="asmVfrom" value="${v.from + 1}"><span>to</span><input type="number" id="asmVto" value="${v.to}"></div>` : ''}
        <label class="asmlbl">Site of insertion:</label>
        <div class="asmlist">${sites.length ? sites.map(sx => `<div class="asmlistitem on">${esc(sx)}</div>`).join('') : '<div class="hint" style="padding:6px">Choose how to linearize the vector.</div>'}</div>
        <div class="asmrow"><span>Orientation of Vector:</span>
          <button class="asmarrow${v.orient === 1 ? ' on' : ''}" data-vori="1">→</button>
          <button class="asmarrow${v.orient === -1 ? ' on' : ''}" data-vori="-1">←</button></div>
      </div>`;
    }

    function fragmentPanel() {
      const i = A.activeFrag, f = A.fragments[i], doc = docOf(f.docId);
      const many = A.fragments.length > 1;
      const cutters = doc ? DNA.ENZYMES.filter(e => (DNA.findSites(doc.seq, doc.circular)[e.name] || []).some(x => x.cut != null)) : [];
      const enzSel = idx => {
        const cur = (f.cuts[idx] || {}).enzyme || '';
        const list = cur && doc ? (DNA.findSites(doc.seq, doc.circular)[cur] || []).filter(x => x.cut != null) : [];
        return `<div class="asmrow"><span>Cut</span>
          <select data-fsite="${idx}" ${list.length > 1 ? '' : 'disabled'}>${list.map((x, k) => `<option value="${k}" ${k === ((f.cuts[idx] || {}).site || 0) ? 'selected' : ''}>${fmt(x.cut)}</option>`).join('') || '<option>—</option>'}</select>
          <span>with</span><select data-fenz="${idx}"><option value="">—</option>${cutters.map(e => `<option ${e.name === cur ? 'selected' : ''}>${e.name}</option>`).join('')}</select></div>`;
      };
      const feats = doc ? doc.features : [];
      return `<div class="asmside">
        <div class="asmrow"><label class="asmlbl" style="flex:1">Source of Fragment${many ? ' ' + (i + 1) : ''}:</label>
          <button class="asmmini" id="asmFadd" title="Add a fragment">+</button>
          <button class="asmmini" id="asmFdel" title="Remove this fragment" ${many ? '' : 'disabled'}>−</button></div>
        ${many ? `<div class="asmtabsmall">${A.fragments.map((_, k) => `<button class="${k === i ? 'on' : ''}" data-fsel="${k}">${k + 1}</button>`).join('')}</div>` : ''}
        <select id="asmFdoc">${docs().map(x => `<option value="${x.id}" ${String(x.id) === String(f.docId) ? 'selected' : ''}>${esc(x.name)}</option>`).join('')}</select>
        <label class="asmradio"><input type="radio" name="fprep" value="pcr" ${f.prep === 'pcr' ? 'checked' : ''}> Use as a template for PCR</label>
        <label class="asmradio"><input type="radio" name="fprep" value="enzymes" ${f.prep === 'enzymes' ? 'checked' : ''}> Cut with restriction enzymes</label>
        <label class="asmradio${doc && doc.circular ? ' dis' : ''}"><input type="radio" name="fprep" value="direct" ${f.prep === 'direct' ? 'checked' : ''} ${doc && doc.circular ? 'disabled' : ''}> Use directly as the Fragment</label>
        ${f.prep === 'enzymes' ? enzSel(0) + enzSel(1) : ''}
        ${f.prep === 'pcr' ? `<label class="asmlbl">Amplified region:</label>
          <div class="asmrow"><input type="number" id="asmFfrom" value="${f.from + 1}"><span>to</span><input type="number" id="asmFto" value="${f.to}">
            <span class="note">${fmt(DNA.rlen({ start: f.from, end: f.to }, doc ? doc.seq.length : 1))} bp</span></div>
          ${feats.length ? `<select id="asmFfeat"><option value="">use a feature…</option>${feats.map((x, k) => `<option value="${k}">${esc(x.name)} (${fmt(x.start + 1)}..${fmt(x.end)})</option>`).join('')}</select>` : ''}` : ''}
        <div class="asmrow"><span>Orientation of Fragment:</span>
          <button class="asmarrow${f.orient === 1 ? ' on' : ''}" data-fori="1">→</button>
          <button class="asmarrow${f.orient === -1 ? ' on' : ''}" data-fori="-1">←</button></div>
      </div>`;
    }

    function productPanel(st) {
      const rows = [];
      if (st.pieces) for (const p of st.pieces) rows.push(`<div class="asmsrc"><b>${esc(p.name)}</b> ${fmt(p.seq.length)} bp<br><span class="note">Source: ${esc(p.doc ? p.doc.name : '—')} · ${esc(p.label)}</span></div>`);
      return `<div class="asmside">
        ${M.joins === 'overlap' ? '<button class="btn" id="asmChoose" style="width:100%">Choose Overlapping PCR Primers…</button>' : ''}
        <label class="asmlbl">Bacterial transformation strain</label>
        <select id="asmStrain"><option value="">Unspecified</option><option>DH5α</option><option>TOP10</option><option>Stbl3</option><option>NEB Stable</option><option>ccdB Survival</option></select>
        <div class="asmsrcs">${rows.join('')}</div>
      </div>`;
    }

    // ---- the window
    function render() {
      const st = compute();
      const vdoc = docOf(A.vector.docId);
      const fsel = A.fragments[A.activeFrag];
      const fdoc = docOf(fsel.docId);
      let mapDoc = vdoc, region = null, selText = '';
      if (A.tab === 'vector') {
        const c0 = (A.vector.cuts[0] || {});
        if (A.vector.prep === 'enzymes' && st.pieces && st.pieces[0]) { region = { mark: st.pieces[0].to, label: c0.enzyme || '' }; selText = st.pieces[0].label + ' — ' + fmt(st.pieces[0].seq.length) + ' bp'; }
        else if (A.vector.prep === 'pcr') { region = { start: A.vector.from, end: A.vector.to }; selText = `${fmt(A.vector.from + 1)} .. ${fmt(A.vector.to)}`; }
        else selText = vdoc ? `${fmt(vdoc.seq.length)} bp, ${vdoc.circular ? 'circular' : 'linear'}` : '';
      } else if (A.tab === 'fragment') {
        mapDoc = fdoc;
        if (fsel.prep === 'pcr') { region = { start: fsel.from, end: fsel.to }; selText = `${fmt(fsel.from + 1)} .. ${fmt(fsel.to)} = ${fmt(DNA.rlen({ start: fsel.from, end: fsel.to }, fdoc ? fdoc.seq.length : 1))} bp`; }
        else if (st.pieces && st.pieces[A.activeFrag + 1]) selText = st.pieces[A.activeFrag + 1].label;
      }
      const tabBtn = (k, l) => `<button class="asmtab${A.tab === k ? ' on' : ''}" data-tab="${k}">${l}</button>`;
      const body = `
        <div class="asmtabs">${tabBtn('vector', 'Vector')}${tabBtn('fragment', 'Fragment')}${tabBtn('product', 'Product')}</div>
        <div class="asmmain">
          <div class="asmleft">
            ${A.tab === 'product'
          ? `<div class="asmprod">${junctionDiagram(st, 560)}</div>`
          : `<div class="asmmapbox" id="asmMapBox">${asmMap(mapDoc, region, 540, 300)}</div>
                 <div class="asmsel">Selected: <b>${esc(selText)}</b></div>`}
          </div>
          <div class="asmright">${A.tab === 'vector' ? vectorPanel() : (A.tab === 'fragment' ? fragmentPanel() : productPanel(st))}</div>
        </div>
        <div class="asmbottom">
          <div class="asmbl">${A.tab === 'product' ? `<div class="asmprimers">${(st.plan && st.plan.junctions || []).map(j => `<div class="asmjrow"><span class="asmjdot ${j.ok ? 'ok' : 'bad'}"></span>${esc(j.from)} → ${esc(j.to)}: ${j.ok ? esc(M.joins === 'overlap' ? `${j.len} bp overlap, Tm ${j.tm.toFixed(1)} °C` : (j.overhang || 'blunt')) : `<span style="color:var(--danger)">${esc(j.why)}</span>`}</div>`).join('')}</div>`
          : primerPanel(st, A.tab === 'vector' ? 0 : A.activeFrag + 1)}</div>
          <div class="asmbr">
            ${statusHTML(st)}
            <div class="asmcreate">
              <label><input type="checkbox" id="asmClose" ${A.product.close ? 'checked' : ''}> Create product and close this window</label>
              <input type="text" id="asmName" value="${esc(A.product.name)}">
            </div>
            <div class="asmbuttons"><button class="btn" id="asmCancel">Cancel</button><button class="btn pri" id="asmGo" ${st.errs.length || st.plan.junctions.some(j => !j.ok) ? 'disabled' : ''}>Assemble</button></div>
          </div>
        </div>`;
      let box = $('#asmWin');
      if (!box) {
        const bg = document.createElement('div');
        // In a window of its own there is nothing behind to dim, so the backdrop is dropped.
        bg.className = ctx.standalone ? 'asmstage' : 'modal-bg';
        bg.id = 'asmBg';
        bg.innerHTML = `<div class="asmdlg${ctx.standalone ? ' standalone' : ''}" id="asmWin">${ctx.standalone ? '' : `<div class="hd">${esc(M.title)}</div>`}<div class="asmbody"></div></div>`;
        document.body.appendChild(bg);
        box = $('#asmWin');
      }
      box.querySelector('.asmbody').innerHTML = body;
      bind();
    }
    const close = () => { const bg = $('#asmBg'); if (bg) bg.remove(); if (ctx.onClose) ctx.onClose(); };

    function bind() {
      const box = $('#asmWin');
      box.onclick = e => {
        const t = e.target;
        const tab = t.closest('[data-tab]'); if (tab) { A.tab = tab.dataset.tab; render(); return; }
        if (t.id === 'asmCancel') { close(); return; }
        if (t.id === 'asmGo') { finish(); return; }
        if (t.id === 'asmFadd') { const last = A.fragments[A.fragments.length - 1]; A.fragments.push({ ...last, cuts: [{}, {}] }); A.activeFrag = A.fragments.length - 1; render(); return; }
        if (t.id === 'asmFdel') { if (A.fragments.length > 1) { A.fragments.splice(A.activeFrag, 1); A.activeFrag = Math.max(0, A.activeFrag - 1); render(); } return; }
        if (t.id === 'asmChoose') { chooseOverlapDialog(); return; }
        const fs = t.closest('[data-fsel]'); if (fs) { A.activeFrag = +fs.dataset.fsel; render(); return; }
        const vo = t.closest('[data-vori]'); if (vo) { A.vector.orient = +vo.dataset.vori; render(); return; }
        const fo = t.closest('[data-fori]'); if (fo) { A.fragments[A.activeFrag].orient = +fo.dataset.fori; render(); return; }
      };
      box.onchange = e => {
        const t = e.target;
        if (t.name === 'vprep') { A.vector.prep = t.value; render(); return; }
        if (t.name === 'fprep') { A.fragments[A.activeFrag].prep = t.value; render(); return; }
        if (t.id === 'asmVdoc') { A.vector.docId = t.value; const d2 = docOf(t.value); A.vector.prep = d2 && d2.circular ? 'enzymes' : 'direct'; A.vector.cuts = [{}, {}]; A.vector.from = 0; A.vector.to = d2 ? d2.seq.length : 0; render(); return; }
        if (t.id === 'asmFdoc') { const f = A.fragments[A.activeFrag]; f.docId = t.value; const d2 = docOf(t.value); f.cuts = [{}, {}]; const c = d2 && d2.features.find(x => x.type === 'CDS'); f.from = c ? c.start : 0; f.to = c ? c.end : Math.min(d2 ? d2.seq.length : 0, 800); render(); return; }
        if (t.id === 'asmKeepInner') { A.vector.keepInner = t.checked; render(); return; }
        if (t.id === 'asmClose') { A.product.close = t.checked; return; }
        if (t.id === 'asmStrain') { A.product.strain = t.value; return; }
        if (t.id === 'asmFfeat') { const d2 = docOf(A.fragments[A.activeFrag].docId); const f2 = d2.features[+t.value]; if (f2) { A.fragments[A.activeFrag].from = f2.start; A.fragments[A.activeFrag].to = f2.end; } render(); return; }
        const ve = t.closest('[data-venz]'); if (ve) { const k = +ve.dataset.venz; A.vector.cuts[k] = { enzyme: t.value, site: 0 }; render(); return; }
        const vs = t.closest('[data-vsite]'); if (vs) { const k = +vs.dataset.vsite; A.vector.cuts[k] = { ...(A.vector.cuts[k] || {}), site: +t.value }; render(); return; }
        const fe = t.closest('[data-fenz]'); if (fe) { const k = +fe.dataset.fenz; A.fragments[A.activeFrag].cuts[k] = { enzyme: t.value, site: 0 }; render(); return; }
        const fsi = t.closest('[data-fsite]'); if (fsi) { const k = +fsi.dataset.fsite; const c = A.fragments[A.activeFrag].cuts[k] || {}; A.fragments[A.activeFrag].cuts[k] = { ...c, site: +t.value }; render(); return; }
      };
      box.oninput = e => {
        const t = e.target;
        if (t.id === 'asmName') { A.product.name = t.value; return; }
        if (t.id === 'asmVfrom') { A.vector.from = Math.max(0, (+t.value || 1) - 1); refreshSoft(); return; }
        if (t.id === 'asmVto') { A.vector.to = +t.value || 0; refreshSoft(); return; }
        if (t.id === 'asmFfrom') { A.fragments[A.activeFrag].from = Math.max(0, (+t.value || 1) - 1); refreshSoft(); return; }
        if (t.id === 'asmFto') { A.fragments[A.activeFrag].to = +t.value || 0; refreshSoft(); return; }
      };
      let soft = null;
      function refreshSoft() { clearTimeout(soft); soft = setTimeout(render, 420); }
      box.addEventListener('keydown', e => { if (e.key === 'Escape') { e.preventDefault(); close(); } });
    }

    function chooseOverlapDialog() {
      const o = A.opts;
      const vdoc = docOf(A.vector.docId);
      const enz = (A.vector.cuts[0] || {}).enzyme;
      ctx.dialog({
        title: 'Choose Overlapping PCR Primers',
        body: `<div class="form">
          <div class="full hint">Regions to amplify: ${A.fragments.map((f, i) => `Fragment${A.fragments.length > 1 ? ' ' + (i + 1) : ''}: ${fmt(f.from + 1)} .. ${fmt(f.to)} = ${fmt(DNA.rlen({ start: f.from, end: f.to }, (docOf(f.docId) || { seq: '' }).seq.length || 1))} bp`).join('<br>')}</div>
          <label for="ovTm">Target Tm for PCR primers</label><div class="row2"><input type="number" id="ovTm" value="${o.primerTm}" min="45" max="75"> °C</div>
          <label>Overlapping ends</label><div>
            <label><input type="radio" name="ovmode" value="fixed" ${o.fixed ? 'checked' : ''}> <input type="number" id="ovFixed" value="${o.fixedLength}" min="12" max="80" style="width:64px"> overlapping bases</label><br>
            <label><input type="radio" name="ovmode" value="range" ${o.fixed ? '' : 'checked'}> <input type="number" id="ovMin" value="${o.minOverlap}" min="10" max="60" style="width:56px"> to <input type="number" id="ovMax" value="${o.maxOverlap}" min="12" max="80" style="width:56px"> overlapping bases</label><br>
            <label style="margin-left:22px">with a target Tm for overlap of <input type="number" id="ovOvTm" value="${o.overlapTm}" min="35" max="70" style="width:56px"> °C</label></div>
          ${enz ? `<label></label><div><label><input type="checkbox" id="ovUp" ${o.regenUp ? 'checked' : ''}> Regenerate the upstream ${esc(enz)} site from the Vector</label><br>
            <label><input type="checkbox" id="ovDown" ${o.regenDown ? 'checked' : ''}> Regenerate the downstream ${esc(enz)} site from the Vector</label></div>` : ''}
          <div class="full err"></div></div>`,
        buttons: [{ label: 'Cancel' }, {
          label: 'Choose Primers', primary: true, act: () => {
            o.primerTm = +val('ovTm') || 60;
            o.fixed = (document.querySelector('input[name=ovmode]:checked') || {}).value === 'fixed';
            o.fixedLength = +val('ovFixed') || 40;
            o.minOverlap = +val('ovMin') || 15;
            o.maxOverlap = Math.max(o.minOverlap, +val('ovMax') || 25);
            o.overlapTm = +val('ovOvTm') || 50;
            if ($('#ovUp')) o.regenUp = $('#ovUp').checked;
            if ($('#ovDown')) o.regenDown = $('#ovDown').checked;
            render();
            ctx.toast('Primers chosen');
          }
        }]
      });
    }

    function finish() {
      const st = compute();
      if (st.errs.length || st.plan.junctions.some(j => !j.ok)) { ctx.toast('The pieces do not join yet.'); return; }
      const prod = st.built;
      const inputs = st.pieces.map((p, i) => ({
        name: p.name, seq: p.seq, circular: false,
        features: [], detail: p.label + (st.plan.primers[i] ? ` · primers ${st.plan.primers[i].fwd.name}, ${st.plan.primers[i].rev.name}` : ''),
        source: p.doc ? { method: 'From ' + p.doc.name, detail: p.label, primers: st.plan.primers[i] ? [st.plan.primers[i].fwd.name, st.plan.primers[i].rev.name] : [], inputs: [{ name: p.doc.name, seq: p.doc.seq, circular: p.doc.circular, features: p.doc.features }] } : null
      }));
      ctx.addDoc({
        name: A.product.name || 'Assembled', seq: prod.seq, circular: true,
        features: prod.features, primers: prod.primers,
        histLabel: `${M.ready.replace('Ready for ', '')} of ${st.pieces.length} pieces`,
        provenance: {
          method: M.title.split(':')[0],
          detail: st.plan.junctions.map(j => `${j.from} → ${j.to}: ${M.joins === 'overlap' ? `${j.len} bp overlap` : (j.overhang || 'blunt')}`).join('  ·  '),
          inputs
        }
      });
      if (st.frames.length) ctx.toast(`Created ${fmt(prod.seq.length)} bp — but ${st.frames.length} reading-frame problem${st.frames.length > 1 ? 's' : ''} remain`);
      else ctx.toast(`Created ${fmt(prod.seq.length)} bp construct`);
      if (A.product.close) close();
      else render();
    }

    render();
  }


  root.ASMUI = { open: openAssembly };
})(typeof window !== 'undefined' ? window : globalThis);
