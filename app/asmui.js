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
  const DNA = root.DNA, ASSEMBLY = root.ASSEMBLY, CLONE = root.CLONE;
  // toolui.js loads after this file, so the shared map is looked up when it is used
  const TOOLUI = () => root.TOOLUI;

  // ================================================================ cloning window
  // One window with Vector, Fragment and Product tabs, shared by every assembly method.
  const val = id => (document.getElementById(id) || {}).value;
  const $ = s => document.querySelector(s);
  const esc = s => String(s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const fmt = n => Number(n).toLocaleString('en-US');
  /* One window serves every method. The dropdown at its top right changes which
   * algorithm runs: how each piece is prepared, and how the junctions are judged. */
  const ASM_METHODS = {
    restriction: { label: 'Restriction Cloning', title: 'Restriction Cloning: Insert Fragment', ready: 'Ready for Ligation', joins: 'ends', preps: ['enzymes', 'pcr', 'direct'] },
    goldengate: { label: 'Golden Gate Assembly', title: 'Golden Gate Assembly: Insert Parts', ready: 'Ready for Golden Gate Assembly', joins: 'ends', preps: ['typeIIs', 'enzymes', 'direct'] },
    gibson: { label: 'Gibson Assembly', title: 'Gibson Assembly: Insert Fragment', ready: 'Ready for Gibson Assembly', joins: 'overlap', preps: ['pcr', 'enzymes', 'direct'] },
    infusion: { label: 'In-Fusion Cloning', title: 'In-Fusion Cloning: Insert Fragment', ready: 'Ready for In-Fusion Cloning', joins: 'overlap', preps: ['pcr', 'enzymes', 'direct'] },
    gateway: { label: 'Gateway BP / LR', title: 'Gateway Recombination: Insert Cassette', ready: 'Ready for Recombination', joins: 'att', preps: ['att'] }
  };
  const PREP_LABEL = {
    pcr: { v: 'Linearize by PCR', f: 'Use as a template for PCR' },
    enzymes: { v: 'Linearize with restriction enzymes', f: 'Cut with restriction enzymes' },
    typeIIs: { v: 'Cut out with a Type IIS enzyme', f: 'Cut out with a Type IIS enzyme' },
    direct: { v: 'Use directly as the linearized Vector', f: 'Use directly as the Fragment' },
    att: { v: 'Recombine at its att sites (destination)', f: 'Recombine at its att sites (entry clone)' }
  };
  const FRAG_COLORS = ['#2f6fd0', '#1a9e5c', '#c9731b', '#8a4ec4', '#c43c3c', '#1b9aa8'];

  function openAssembly(method, ctx) {
    const docs = () => ctx.docs();
    if (docs().length < 1) { ctx.toast('Open a sequence first.'); return; }
    const first = ctx.current && docs().some(x => String(x.id) === String(ctx.current)) ? docs().find(x => String(x.id) === String(ctx.current)) : docs()[0];
    const A = {
      method, tab: 'vector', activeFrag: 0,
      vector: { docId: first.id, prep: first.circular ? 'enzymes' : 'direct', cuts: [], keepInner: false, from: 0, to: first.seq.length, orient: 1, enzyme: 'BsaI', fragment: 0 },
      fragments: [{ docId: (docs()[1] || first).id, prep: 'pcr', cuts: [], from: 0, to: 0, orient: 1, enzyme: 'BsaI', fragment: 0 }],
      opts: { minOverlap: 15, maxOverlap: 25, fixed: false, fixedLength: 40, overlapTm: 50, primerTm: 60, regenUp: false, regenDown: false },
      product: { name: 'Assembled', close: true, strain: '' },
      state: null
    };
    const f0 = docs()[1] || first;
    const firstCds = f0.features.find(x => x.type === 'CDS');
    if (firstCds) { A.fragments[0].from = firstCds.start; A.fragments[0].to = firstCds.end; }
    else { A.fragments[0].from = 0; A.fragments[0].to = Math.min(f0.seq.length, 800); }

    const docOf = id => docs().find(x => String(x.id) === String(id));
    const meth = () => ASM_METHODS[A.method] || ASM_METHODS.restriction;
    /* Each method allows its own ways of preparing a piece, so changing the method moves
     * any piece whose preparation the new method does not offer onto its first one. */
    const prepsFor = () => meth().preps;
    const fixPreps = () => {
      const ok = prepsFor();
      // Golden Gate and Gateway have a way of their own to cut a piece out, so that is
      // where they start; the other methods leave a piece as it was set up.
      const force = ok[0] === 'typeIIs' || ok[0] === 'att';
      const fit = sel => { if (force || !ok.includes(sel.prep)) sel.prep = ok[0]; };
      fit(A.vector); A.fragments.forEach(fit);
    };
    function setMethod(m) {
      if (!ASM_METHODS[m]) return;
      A.method = m;
      fixPreps();
      render();
    }
    fixPreps();

    // ---- build the ordered pieces from the current settings
    function compute() {
      const errs = [];
      let pieces = null, plan = null, built = null, frames = [], suggest = null;
      try {
        const vdoc = docOf(A.vector.docId);
        if (!vdoc) throw new Error('Choose a vector.');
        const V = ASSEMBLY.prepare(vdoc, specOf(A.vector, vdoc, 'vector'));
        V.name = 'Linearized Vector'; V.spec = specOf(A.vector, vdoc, 'vector'); V.role = 'vector';
        pieces = [A.vector.orient === -1 ? Object.assign(ASSEMBLY.flip(V), { name: V.name, spec: V.spec, role: 'vector', doc: V.doc, from: V.from, to: V.to, label: V.label }) : V];
        A.fragments.forEach((f, i) => {
          const fd = docOf(f.docId);
          if (!fd) throw new Error(`Choose a source for fragment ${i + 1}.`);
          const P = ASSEMBLY.prepare(fd, specOf(f, fd, 'fragment'));
          P.name = A.fragments.length > 1 ? `Fragment ${i + 1}` : 'Fragment';
          P.spec = specOf(f, fd, 'fragment'); P.role = 'fragment'; P.index = i;
          pieces.push(f.orient === -1 ? Object.assign(ASSEMBLY.flip(P), { name: P.name, spec: P.spec, role: 'fragment', index: i, doc: P.doc, from: P.from, to: P.to, label: P.label }) : P);
        });
        const planner = ps => meth().joins === 'overlap' ? ASSEMBLY.planGibson(ps, A.opts)
          : meth().joins === 'att' ? ASSEMBLY.planAtt(ps) : ASSEMBLY.planRestriction(ps);
        plan = planner(pieces);
        built = ASSEMBLY.build(pieces, plan, A.product.name);
        frames = ASSEMBLY.frameIssues(built, built.bounds);
        if (plan.junctions.some(j => !j.ok)) suggest = ASSEMBLY.suggestFlips(pieces, planner);
      } catch (e) { errs.push(e.message); }
      A.state = { pieces, plan, built, frames, errs, suggest };
      return A.state;
    }
    function specOf(sel, doc, role) {
      if (sel.prep === 'enzymes') {
        const cuts = (sel.cuts || []).filter(c => c && c.enzyme).map(c => {
          const list = (DNA.findSites(doc.seq, doc.circular)[c.enzyme] || []).filter(x => x.cut != null);
          const s = list[Math.min(c.site || 0, list.length - 1)];
          return s ? { enzyme: c.enzyme, cut: s.cut, cut2: s.cut2 } : null;
        }).filter(Boolean);
        return { prep: 'enzymes', cuts, keepInner: !!sel.keepInner };
      }
      if (sel.prep === 'pcr') return { prep: 'pcr', from: sel.from | 0, to: sel.to | 0 };
      if (sel.prep === 'typeIIs') return { prep: 'typeIIs', enzyme: sel.enzyme, fragment: sel.fragment | 0 };
      if (sel.prep === 'att') return { prep: 'att', role: role === 'vector' ? 'acceptor' : 'donor' };
      return { prep: 'direct' };
    }

    /* The Vector and Fragment maps are the shared interactive map, so an enzyme, a
     * feature or a primer can be chosen on the picture instead of in the lists. */
    /* Each tab keeps its own view state for the strip beside the map, so hiding a
     * primer on the Fragment map does not change what the Vector map shows. */
    const mapStates = {};
    function viewState(key) {
      if (!mapStates[key]) mapStates[key] = TOOLUI().mapState({ showPrimers: true });
      return mapStates[key];
    }
    function pieceMap(doc, region, w) {
      if (!doc) return '<div class="hint">Choose a sequence.</div>';
      const sel = A.tab === 'vector' ? A.vector : A.fragments[A.activeFrag];
      const chosen = (sel.cuts || []).map(c => c.enzyme).filter(Boolean);
      const key = A.tab === 'vector' ? 'vector' : 'fragment' + A.activeFrag;
      const st = viewState(key);
      const sites = DNA.findSites(doc.seq, doc.circular);
      // whatever this piece is cut with stays on the map, whichever set is chosen
      TOOLUI().includeEnz(st, sites, chosen.concat(sel.prep === 'typeIIs' ? [sel.enzyme] : []));
      const primers = (doc.primers || []).map(p => {
        const h = DNA.searchDNA(doc.seq, doc.circular, p.seq)[0];
        return h ? { ...p, start: h.start, end: h.end, strand: h.strand, color: (A.pick && A.pick.primers || []).includes(p.id) ? '#c0392b' : '#2f6fd0' } : null;
      }).filter(Boolean);
      return TOOLUI().mapPane({
        doc, state: st, sites, primers,
        sel: region && region.start !== undefined ? region : null,
        W: w, H: 230
      }) + `<div class="asmmaphint">${sel.prep === 'enzymes'
        ? `Click an enzyme to cut with it, shift-click to add another${chosen.length ? ` \u2014 now ${esc(chosen.join(' + '))}` : ''}.`
        : sel.prep === 'pcr' ? 'Click a feature to amplify it, or a primer to set one end (shift-click the other).'
          : 'Click a feature to look at it.'}</div>`;
    }
    function productMap(st, w) {
      if (!st.built) return '<div class="hint">The pieces do not make a product yet.</div>';
      const b = st.built;
      const doc = { name: A.product.name || 'Product', seq: b.seq, circular: true, features: b.features || [], primers: b.primers || [] };
      return TOOLUI().mapPane({
        doc, state: viewState('product'),
        sel: A.pick && A.pick.tab === 'product' ? A.pick : null,
        W: w, H: 230
      }) + `<div class="asmmaphint">Click a feature on the product to select it, then copy it.</div>`;
    }
    // the document a tab's map is drawn from, for the strip's own lists
    function mapDocOf() {
      if (A.tab === 'product') {
        const b = A.state && A.state.built;
        return b ? { name: A.product.name || 'Product', seq: b.seq, circular: true, features: b.features || [], primers: b.primers || [] } : null;
      }
      return A.tab === 'vector' ? docOf(A.vector.docId) : docOf(A.fragments[A.activeFrag].docId);
    }
    function mapStateKey() {
      return A.tab === 'product' ? 'product' : (A.tab === 'vector' ? 'vector' : 'fragment' + A.activeFrag);
    }

    function selLine(pick, st, bare) {
      if (!pick) return 'Nothing selected on the map.';
      const len = pick.end - pick.start;
      return `${bare ? 'on the map:' : 'Selected'} <b>${esc(pick.name || '')}</b> ${fmt(pick.start + 1)}..${fmt(pick.end)} = <b>${fmt(len)} bp</b>
        <button class="asmmini" id="asmCopySel" title="Copy this stretch">Copy</button>
        ${A.tab === 'product' && st.built ? `<button class="asmmini" id="asmCopyAll" title="Copy the whole product">Copy all ${fmt(st.built.seq.length)} bp</button>` : ''}`;
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
        else {
          // a narrow piece at either edge would print its name off the picture
          const lx = Math.min(Math.max(x + pw / 2, m + 4), m + bw - 4);
          const anchor = lx <= m + 30 ? 'start' : (lx >= m + bw - 30 ? 'end' : 'middle');
          g += `<text x="${lx.toFixed(1)}" y="${barY + barH + 13}" text-anchor="${anchor}" font-size="10" fill="#5e6776">${esc(p.name)}</text>`;
        }
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
        const txt = meth().joins === 'overlap'
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
        const lx2 = Math.min(Math.max(ax, m + 4), m + bw - 4);
        const an2 = lx2 <= m + 24 ? 'start' : (lx2 >= m + bw - 24 ? 'end' : 'middle');
        g += `<text x="${lx2.toFixed(1)}" y="${ay + 15}" text-anchor="${an2}" font-size="9.5" fill="#8d95a3">${rev ? 'reverse' : 'forward'}</text>`;
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
      return `<div class="asmstat ok"><b>${esc(meth().ready)}</b><br>Product: ${fmt(st.built.seq.length)} bp${warn}${methodNotes(st)}</div>`;
    }

    /* What each method has to say about a product that is otherwise ready: Golden Gate
     * warns about repeated overhangs and sites left inside a part, Gateway reports the
     * by-product the reaction also makes. */
    function methodNotes(st) {
      const out = [];
      if (A.method === 'goldengate' && st.pieces) {
        const ends = st.pieces.map(p => p.ends.right && p.ends.right.x).filter(Boolean);
        const dupes = [...new Set(ends.filter((x, i) => ends.indexOf(x) !== i))];
        if (dupes.length) out.push(`\u26a0 the overhang ${dupes.map(esc).join(' and ')} is used more than once, so the order is ambiguous`);
        for (const p of st.pieces) {
          if (!p.doc || p.spec.prep !== 'typeIIs') continue;
          const k = (DNA.findSites(p.doc.seq, p.doc.circular)[p.spec.enzyme] || []).length;
          if (k > 2) out.push(`\u26a0 ${esc(p.doc.name)} has ${k} ${esc(p.spec.enzyme)} sites, so ${k - 2} lie inside the part \u2014 remove them first`);
        }
        const left = st.built ? DNA.searchDNA(st.built.seq, true, (DNA.ENZ_BY_NAME[st.pieces[0].spec.enzyme] || { site: 'X' }).site).length : 0;
        if (left) out.push(`\u26a0 the product still carries ${left} ${esc(st.pieces[0].spec.enzyme)} site${left === 1 ? '' : 's'}`);
      }
      if (A.method === 'gateway' && st.pieces && st.pieces.length === 2) {
        const v = st.pieces[0], f = st.pieces[1];
        if (v.att && f.att) {
          const by = DNA.rlen(v.att.inner, v.doc.seq.length) + DNA.rlen(f.att.outer, f.doc.seq.length);
          out.push(`The reaction also makes a by-product of about ${fmt(by)} bp.`);
        }
      }
      return out.length ? `<br><span class="note">${out.join('<br>')}</span>` : '';
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
    const cuttersOf = doc => doc ? DNA.ENZYMES.filter(e => (DNA.findSites(doc.seq, doc.circular)[e.name] || []).some(x => x.cut != null)) : [];
    /* One row per chosen enzyme, plus an empty row at the end: choosing an enzyme there
     * opens the next one, so a piece can be cut with as many enzymes as wanted. */
    function enzChain(sel, doc, who) {
      const cutters = cuttersOf(doc);
      const chosen = (sel.cuts || []).filter(c => c && c.enzyme);
      const rows = chosen.concat([{}]).map((c, idx) => {
        const cur = c.enzyme || '';
        const list = cur && doc ? (DNA.findSites(doc.seq, doc.circular)[cur] || []).filter(x => x.cut != null) : [];
        const site = list.length > 1
          ? `<select data-${who}site="${idx}">${list.map((x, k) => `<option value="${k}" ${k === (c.site || 0) ? 'selected' : ''}>at ${fmt(x.cut)}</option>`).join('')}</select>`
          : (list.length === 1 ? `<span class="note">at ${fmt(list[0].cut)}</span>` : '');
        return `<div class="asmrow"><span>${idx === 0 ? 'Cut with' : 'and'}</span>
          <select data-${who}enz="${idx}"><option value="">${idx === 0 ? '\u2014' : 'add an enzyme\u2026'}</option>${cutters.map(e => `<option ${e.name === cur ? 'selected' : ''}>${e.name}</option>`).join('')}</select>
          ${site}${cur ? `<button class="asmmini" data-${who}del="${idx}" title="Remove this enzyme">\u2212</button>` : ''}</div>`;
      });
      return rows.join('');
    }
    function typeIIsBox(sel, doc, who) {
      const list = CLONE && CLONE.TYPE2S ? CLONE.TYPE2S : ['BsaI', 'BsmBI', 'BbsI', 'SapI'];
      const n = doc ? (DNA.findSites(doc.seq, doc.circular)[sel.enzyme] || []).length : 0;
      let parts = 0;
      try { parts = doc ? CLONE.ggFragments(doc, sel.enzyme).length : 0; } catch (e) { parts = 0; }
      return `<div class="asmrow"><span>Cut out with</span>
          <select data-${who}t2s><option value="">\u2014</option>${list.map(nm => `<option ${nm === sel.enzyme ? 'selected' : ''}>${nm}</option>`).join('')}</select>
          <span class="note">${n} site${n === 1 ? '' : 's'}</span></div>
        ${parts > 1 ? `<div class="asmrow"><span>Part</span><select data-${who}t2f>${Array.from({ length: parts }, (_, k) => `<option value="${k}" ${k === (sel.fragment | 0) ? 'selected' : ''}>${k + 1} of ${parts}</option>`).join('')}</select><span class="note">largest first</span></div>` : ''}`;
    }
    function prepRadios(sel, who) {
      const doc = docOf(sel.docId);
      const kind = who === 'v' ? 'v' : 'f';
      return prepsFor().map(k => {
        const dis = k === 'direct' && doc && doc.circular;
        return `<label class="asmradio${dis ? ' dis' : ''}"><input type="radio" name="${kind}prep" value="${k}" ${sel.prep === k ? 'checked' : ''} ${dis ? 'disabled' : ''}> ${esc(PREP_LABEL[k][who === 'v' ? 'v' : 'f'])}</label>`;
      }).join('');
    }
    function prepBody(sel, who) {
      const doc = docOf(sel.docId);
      if (sel.prep === 'enzymes') {
        const many = (sel.cuts || []).filter(c => c && c.enzyme).length > 1;
        return enzChain(sel, doc, who) +
          (many ? `<label class="asmradio"><input type="checkbox" data-${who}inner ${sel.keepInner ? 'checked' : ''}> Keep the stretch between the cuts instead</label>` : '');
      }
      if (sel.prep === 'typeIIs') return typeIIsBox(sel, doc, who);
      if (sel.prep === 'att') {
        const sp = doc ? CLONE.attSplit(doc) : null;
        return `<div class="hint" style="padding:4px 0">${sp
          ? `${esc(sp.site1.kind)} at ${fmt(sp.site1.pos + 1)} and ${esc(sp.site2.kind)} at ${fmt(sp.site2.pos + 1)}.`
          : 'No att1 / att2 pair was found in this sequence.'}</div>`;
      }
      if (sel.prep === 'pcr') {
        const n = doc ? doc.seq.length : 1;
        return `<div class="asmrow"><span>Keep</span><input type="number" id="asm${who === 'v' ? 'V' : 'F'}from" value="${sel.from + 1}"><span>to</span><input type="number" id="asm${who === 'v' ? 'V' : 'F'}to" value="${sel.to}">
          <span class="note">${fmt(DNA.rlen({ start: sel.from, end: sel.to }, n))} bp</span></div>`;
      }
      return '';
    }

    function vectorPanel() {
      const v = A.vector, doc = docOf(v.docId);
      const st = A.state;
      const sites = st && st.pieces && st.pieces[0] ? [st.pieces[0].label] : [];
      return `<div class="asmside">
        <label class="asmlbl">${meth().joins === 'att' ? 'Destination vector:' : 'Vector:'}</label>
        <select id="asmVdoc">${docs().map(x => `<option value="${x.id}" ${String(x.id) === String(v.docId) ? 'selected' : ''}>${esc(x.name)}</option>`).join('')}</select>
        ${prepRadios(v, 'v')}
        ${prepBody(v, 'v')}
        <label class="asmlbl">Site of insertion:</label>
        <div class="asmlist">${sites.length ? sites.map(sx => `<div class="asmlistitem on">${esc(sx)}</div>`).join('') : '<div class="hint" style="padding:6px">Choose how to linearize the vector.</div>'}</div>
        <div class="asmrow"><span>Orientation of Vector:</span>
          <button class="asmarrow${v.orient === 1 ? ' on' : ''}" data-vori="1">\u2192</button>
          <button class="asmarrow${v.orient === -1 ? ' on' : ''}" data-vori="-1">\u2190</button></div>
      </div>`;
    }

    function fragmentPanel() {
      const i = A.activeFrag, f = A.fragments[i], doc = docOf(f.docId);
      const many = A.fragments.length > 1;
      const feats = doc ? doc.features : [];
      return `<div class="asmside">
        <div class="asmrow"><label class="asmlbl" style="flex:1">${meth().joins === 'att' ? 'Entry clone' : 'Source of Fragment'}${many ? ' ' + (i + 1) : ''}:</label>
          <button class="asmmini" id="asmFadd" title="Add a fragment">+</button>
          <button class="asmmini" id="asmFdel" title="Remove this fragment" ${many ? '' : 'disabled'}>\u2212</button></div>
        ${many ? `<div class="asmtabsmall">${A.fragments.map((_, k) => `<button class="${k === i ? 'on' : ''}" data-fsel="${k}">${k + 1}</button>`).join('')}</div>` : ''}
        <select id="asmFdoc">${docs().map(x => `<option value="${x.id}" ${String(x.id) === String(f.docId) ? 'selected' : ''}>${esc(x.name)}</option>`).join('')}</select>
        ${prepRadios(f, 'f')}
        ${prepBody(f, 'f')}
        ${f.prep === 'pcr' && feats.length ? `<select id="asmFfeat"><option value="">use a feature\u2026</option>${feats.map((x, k) => `<option value="${k}">${esc(x.name)} (${fmt(x.start + 1)}..${fmt(x.end)})</option>`).join('')}</select>` : ''}
        <div class="asmrow"><span>Orientation of Fragment:</span>
          <button class="asmarrow${f.orient === 1 ? ' on' : ''}" data-fori="1">\u2192</button>
          <button class="asmarrow${f.orient === -1 ? ' on' : ''}" data-fori="-1">\u2190</button></div>
      </div>`;
    }

    function productPanel(st) {
      const rows = [];
      if (st.pieces) for (const p of st.pieces) rows.push(`<div class="asmsrc"><b>${esc(p.name)}</b> ${fmt(p.seq.length)} bp<br><span class="note">Source: ${esc(p.doc ? p.doc.name : '—')} · ${esc(p.label)}</span></div>`);
      return `<div class="asmside">
        ${meth().joins === 'overlap' ? '<button class="btn" id="asmChoose" style="width:100%">Choose Overlapping PCR Primers…</button>' : ''}
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
      const methodSel = `<select id="asmMethod" class="asmmethod" title="Which cloning method this window runs">${Object.keys(ASM_METHODS).map(k => `<option value="${k}" ${k === A.method ? 'selected' : ''}>${esc(ASM_METHODS[k].label)}</option>`).join('')}</select>`;
      const picked = A.pick && A.pick.tab === A.tab ? A.pick : null;
      const body = `
        <div class="asmtabs">${tabBtn('vector', 'Vector')}${tabBtn('fragment', 'Fragment')}${tabBtn('product', 'Product')}<span class="asmspacer"></span>${methodSel}</div>
        <div class="asmmain">
          <div class="asmleft">
            ${A.tab === 'product'
          ? `<div class="asmmapbox" id="asmMapBox">${productMap(st, 640)}</div>
                 <div class="asmsel">${selLine(picked, st)}</div>
                 <div class="asmprod">${junctionDiagram(st, 640)}</div>`
          : `<div class="asmmapbox" id="asmMapBox">${pieceMap(mapDoc, region, 640)}</div>
                 <div class="asmsel">Selected: <b>${esc(selText)}</b>${picked ? ` \u00b7 ${selLine(picked, st, true)}` : ''}</div>`}
          </div>
          <div class="asmright">${A.tab === 'vector' ? vectorPanel() : (A.tab === 'fragment' ? fragmentPanel() : productPanel(st))}</div>
        </div>
        <div class="asmbottom">
          <div class="asmbl">${A.tab === 'product' ? `<div class="asmprimers">${(st.plan && st.plan.junctions || []).map(j => `<div class="asmjrow"><span class="asmjdot ${j.ok ? 'ok' : 'bad'}"></span>${esc(j.from)} → ${esc(j.to)}: ${j.ok ? esc(meth().joins === 'overlap' ? `${j.len} bp overlap, Tm ${j.tm.toFixed(1)} °C` : (j.overhang || 'blunt')) : `<span style="color:var(--danger)">${esc(j.why)}</span>`}</div>`).join('')}</div>`
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
        bg.innerHTML = `<div class="asmdlg${ctx.standalone ? ' standalone' : ''}" id="asmWin">${ctx.standalone ? '' : `<div class="hd">${esc(meth().title)}</div>`}<div class="asmbody"></div></div>`;
        document.body.appendChild(bg);
        box = $('#asmWin');
      }
      const hd = box.querySelector('.hd'); if (hd) hd.textContent = meth().title;
      box.querySelector('.asmbody').innerHTML = body;
      bind();
      const mdoc = mapDocOf();
      if (mdoc) TOOLUI().bindMapPane(box, mdoc, viewState(mapStateKey()), () => render());
    }
    /* Escape closes the window wherever the focus happens to be, but not while a dialog
     * of its own is on top of it. */
    const onKey = e => {
      if (e.key !== 'Escape' || !$('#asmWin')) return;
      const modal = document.getElementById('modal');
      if (modal && !modal.hidden) return;
      // the strip's panel closes first, so Escape does not take the window with it
      const pop = document.getElementById('tmspop');
      if (pop && !pop.hidden) { e.preventDefault(); pop.hidden = true; pop.innerHTML = ''; return; }
      e.preventDefault(); close();
    };
    document.addEventListener('keydown', onKey);
    const close = () => {
      document.removeEventListener('keydown', onKey);
      const bg = $('#asmBg'); if (bg) bg.remove();
      if (ctx.onClose) ctx.onClose();
    };

    const copyText = (txt, msg) => {
      if (ctx.copy) { ctx.copy(txt, msg); return; }
      navigator.clipboard.writeText(txt).then(() => ctx.toast(msg || 'Copied')).catch(() => ctx.toast('Could not copy'));
    };
    function copySelection() {
      const pick = A.pick;
      if (!pick) { ctx.toast('Nothing is selected on the map.'); return; }
      let seq = '';
      if (pick.tab === 'product') { const st = A.state; seq = st && st.built ? st.built.seq.slice(pick.start, pick.end) : ''; }
      else {
        const doc = pick.tab === 'vector' ? docOf(A.vector.docId) : docOf(A.fragments[A.activeFrag].docId);
        if (doc) seq = DNA.sub(doc.seq, pick.start, pick.end);
      }
      if (!seq) { ctx.toast('Nothing to copy.'); return; }
      copyText(seq, `Copied ${fmt(seq.length)} bp`);
    }
    // a map title reads "name 123..456"; the coordinates are printed separately
    const titleName = el => String((el.querySelector('title') || {}).textContent || '')
      .replace(/\s+[\d,]+\.\.[\d,]+\s*$/, '').trim();

    /* A click on the map is the same choice as a click in the lists beside it. */
    function mapPick(el, shift) {
      const id = el.dataset.pick || '';
      const a = +el.dataset.a, b = +el.dataset.b;
      const onVector = A.tab === 'vector';
      const sel = onVector ? A.vector : A.fragments[A.activeFrag];
      if (A.tab === 'product') {
        A.pick = { tab: 'product', start: Math.min(a, b), end: Math.max(a, b), name: titleName(el) };
        render(); return;
      }
      if (id.startsWith('e:')) {
        const nm = id.split(':')[1];
        if (sel.prep === 'typeIIs') { sel.enzyme = nm; sel.fragment = 0; render(); return; }
        if (!prepsFor().includes('enzymes')) { ctx.toast(`${meth().label} does not cut with restriction enzymes.`); return; }
        sel.prep = 'enzymes';
        const chosen = (sel.cuts || []).filter(c => c && c.enzyme);
        if (shift) { if (!chosen.some(c => c.enzyme === nm)) chosen.push({ enzyme: nm, site: 0 }); }
        else { chosen.length = 0; chosen.push({ enzyme: nm, site: 0 }); }
        sel.cuts = chosen;
        if (chosen.length < 2) sel.keepInner = false;
        render(); return;
      }
      if (id.startsWith('p:')) {
        if (sel.prep !== 'pcr') { ctx.toast('Switch this piece to PCR to choose primers.'); return; }
        const pid = +id.split(':')[1];
        const list = ((A.pick && A.pick.primers) || []).slice();
        if (shift) { if (!list.includes(pid)) list.push(pid); while (list.length > 2) list.shift(); }
        else list.splice(0, list.length, pid);
        if (shift || list.length === 1) {
          if (!shift) sel.from = Math.min(a, b);
          else sel.to = Math.max(a, b);
        }
        A.pick = { tab: A.tab, start: sel.from, end: sel.to, name: 'PCR region', primers: list };
        render(); return;
      }
      if (id.startsWith('f:')) {
        if (sel.prep === 'pcr') { sel.from = Math.min(a, b); sel.to = Math.max(a, b); }
        A.pick = { tab: A.tab, start: Math.min(a, b), end: Math.max(a, b), name: titleName(el) };
        render(); return;
      }
    }
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
        const vd = t.closest('[data-vdel]'); if (vd) { A.vector.cuts = A.vector.cuts.filter((c, k) => k !== +vd.dataset.vdel); if (A.vector.cuts.length < 2) A.vector.keepInner = false; render(); return; }
        const fd2 = t.closest('[data-fdel]'); if (fd2) { const f = A.fragments[A.activeFrag]; f.cuts = f.cuts.filter((c, k) => k !== +fd2.dataset.fdel); if (f.cuts.length < 2) f.keepInner = false; render(); return; }
        const fs = t.closest('[data-fsel]'); if (fs) { A.activeFrag = +fs.dataset.fsel; render(); return; }
        if (t.id === 'asmCopySel' || t.closest('#asmCopySel')) { copySelection(); return; }
        if (t.id === 'asmCopyAll' || t.closest('#asmCopyAll')) { const st2 = A.state; if (st2 && st2.built) copyText(st2.built.seq, `Copied ${fmt(st2.built.seq.length)} bp`); return; }
        const pk = t.closest('#asmMapBox [data-pick]');
        if (pk) { mapPick(pk, e.shiftKey); return; }
        const vo = t.closest('[data-vori]'); if (vo) { A.vector.orient = +vo.dataset.vori; render(); return; }
        const fo = t.closest('[data-fori]'); if (fo) { A.fragments[A.activeFrag].orient = +fo.dataset.fori; render(); return; }
      };
      box.onchange = e => {
        const t = e.target;
        if (t.name === 'vprep') { A.vector.prep = t.value; render(); return; }
        if (t.name === 'fprep') { A.fragments[A.activeFrag].prep = t.value; render(); return; }
        if (t.id === 'asmVdoc') { A.vector.docId = t.value; const d2 = docOf(t.value); A.vector.prep = prepsFor().includes(d2 && d2.circular ? 'enzymes' : 'direct') ? (d2 && d2.circular ? 'enzymes' : 'direct') : prepsFor()[0]; A.vector.cuts = []; A.vector.from = 0; A.vector.to = d2 ? d2.seq.length : 0; render(); return; }
        if (t.id === 'asmFdoc') { const f = A.fragments[A.activeFrag]; f.docId = t.value; const d2 = docOf(t.value); f.cuts = []; const c = d2 && d2.features.find(x => x.type === 'CDS'); f.from = c ? c.start : 0; f.to = c ? c.end : Math.min(d2 ? d2.seq.length : 0, 800); render(); return; }
        if (t.id === 'asmClose') { A.product.close = t.checked; return; }
        if (t.id === 'asmStrain') { A.product.strain = t.value; return; }
        if (t.id === 'asmFfeat') { const d2 = docOf(A.fragments[A.activeFrag].docId); const f2 = d2.features[+t.value]; if (f2) { A.fragments[A.activeFrag].from = f2.start; A.fragments[A.activeFrag].to = f2.end; } render(); return; }
        // the enzyme chain: setting a row's enzyme replaces it, clearing one drops the row
        const setCut = (sel, k, enzyme) => {
          const chosen = (sel.cuts || []).filter(c => c && c.enzyme);
          if (!enzyme) chosen.splice(k, 1);
          else if (k < chosen.length) chosen[k] = { enzyme, site: 0 };
          else chosen.push({ enzyme, site: 0 });
          sel.cuts = chosen;
          if (chosen.length < 2) sel.keepInner = false;
        };
        const ve = t.closest('[data-venz]'); if (ve) { setCut(A.vector, +ve.dataset.venz, t.value); render(); return; }
        const fe = t.closest('[data-fenz]'); if (fe) { setCut(A.fragments[A.activeFrag], +fe.dataset.fenz, t.value); render(); return; }
        const vs = t.closest('[data-vsite]'); if (vs) { const k = +vs.dataset.vsite; const c = A.vector.cuts[k]; if (c) c.site = +t.value; render(); return; }
        const fsi = t.closest('[data-fsite]'); if (fsi) { const k = +fsi.dataset.fsite; const c = A.fragments[A.activeFrag].cuts[k]; if (c) c.site = +t.value; render(); return; }
        if (t.hasAttribute('data-vinner')) { A.vector.keepInner = t.checked; render(); return; }
        if (t.hasAttribute('data-finner')) { A.fragments[A.activeFrag].keepInner = t.checked; render(); return; }
        if (t.hasAttribute('data-vt2s')) { A.vector.enzyme = t.value; A.vector.fragment = 0; render(); return; }
        if (t.hasAttribute('data-ft2s')) { A.fragments[A.activeFrag].enzyme = t.value; A.fragments[A.activeFrag].fragment = 0; render(); return; }
        if (t.hasAttribute('data-vt2f')) { A.vector.fragment = +t.value; render(); return; }
        if (t.hasAttribute('data-ft2f')) { A.fragments[A.activeFrag].fragment = +t.value; render(); return; }
        if (t.id === 'asmMethod') { setMethod(t.value); return; }
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
        histLabel: `${meth().ready.replace('Ready for ', '')} of ${st.pieces.length} pieces`,
        provenance: {
          method: meth().title.split(':')[0],
          detail: st.plan.junctions.map(j => `${j.from} → ${j.to}: ${meth().joins === 'overlap' ? `${j.len} bp overlap` : (j.overhang || 'blunt')}`).join('  ·  '),
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
