/* Dogma Studio — the desktop pass.
 *
 * Opens app/index.html in a desktop-sized window and walks the application: every view,
 * the dialogs, the side strip, the maps, the cloning window and the simulation. It is
 * the half of `npm test` that answers "is the PC build still right".
 */
const path = require('path');
const { chromium } = require('playwright');
const { runner } = require('./lib');
const INDEX = 'file://' + path.join(__dirname, '..', 'app', 'index.html');
(async () => {
  const b = await chromium.launch();
  const p = await b.newPage({ viewport: { width: 1440, height: 900 } });
  const errs = [];
  p.on('pageerror', e => errs.push('PAGEERROR: ' + e.message));
  p.on('console', m => { if (m.type() === 'error') errs.push('CONSOLE: ' + m.text()); });
  await p.goto(INDEX);
  await p.waitForTimeout(900);

  const step = runner(errs);

  // every view
  for (const v of ['map', 'seq', 'enz', 'feat', 'primer', 'align', 'compare', 'hist']) {
    await step('view ' + v, async () => { await p.click('.vt[data-v="' + v + '"]'); await p.waitForTimeout(350); });
  }
  await step('linear map', async () => { await p.click('.vt[data-v="map"]'); await p.click('.tb[data-a="lin"]'); await p.waitForTimeout(300); });
  await step('circular map', async () => { await p.click('.tb[data-a="circ"]'); await p.waitForTimeout(300); });
  // the side strip replaced the old left panel
  await step('ORFs on', async () => { await p.click('[data-ms="orf"]'); await p.waitForTimeout(400); });
  await step('ORFs off', async () => { await p.click('[data-ms="orf"]'); await p.waitForTimeout(300); });
  await step('enzyme panel opens', async () => { await p.click('[data-mschev="enz"]'); await p.waitForTimeout(350); if (!(await p.$('#mspop [data-msv="enzMode"]'))) throw new Error('no enzyme panel'); });
  await step('enzyme mode all-unique', async () => { await p.selectOption('#mspop [data-msv="enzMode"]', 'unique'); await p.waitForTimeout(400); });
  await step('enzyme mode back', async () => { await p.selectOption('#mspop [data-msv="enzMode"]', 'unique6'); await p.waitForTimeout(300); await p.evaluate(() => { const x = document.querySelector('#mspop'); if (x) x.hidden = true; }); });

  // editing + history
  await step('sequence typing', async () => {
    await p.click('.vt[data-v="seq"]'); await p.waitForTimeout(300);
    await p.click('#seqscroll .srow .nt', { position: { x: 40, y: 8 } });
    await p.keyboard.type('GATTACA'); await p.waitForTimeout(300);
  });
  await step('undo', async () => { await p.keyboard.press('Control+z'); await p.waitForTimeout(300); });
  await step('find', async () => { await p.selectOption('#findmode','dna').catch(()=>{}); await p.fill('#findq', 'GAATTC'); await p.waitForTimeout(400); });

  // dialogs that must open and close
  const dialogs = [
    ['Actions', 'Restriction Digest & Gel…'],
    ['Actions', 'Translate Selection…'],
    ['Tools', 'Design PCR Primers…'],
    ['Tools', 'Design Sequencing Primers…'],
    ['Tools', 'Codon Optimization…'],
    ['Analyze', 'Choose Enzymes…'],
    ['Analyze', 'Enzyme Database…'],
    ['Analyze', 'Feature Database…'],
    ['Analyze', 'RNA Secondary Structure…'],
    ['Analyze', 'Protein Properties…'],
    ['File', 'Print…'],
    ['File', 'Open from NCBI…'],
    ['Help', 'Keyboard Shortcuts'],
    ['Help', 'About Dogma Studio'],
    ['Edit', 'Go To Position…'],
    ['Features', 'Add Feature…'],
    ['Primers', 'Add Primer…']
  ];
  await step('a click lands on the base it was aimed at', async () => {
    await p.evaluate(() => window.__gc.setView('seq'));
    await p.waitForTimeout(600);
    const g = await p.evaluate(() => {
      const V = window.__gc.S.seq, nt = document.querySelector('#seqrows .nt');
      const r = nt.closest('.srow').getBoundingClientRect();
      return { G: V.G, left: r.left, top: r.top, CW: nt.getBoundingClientRect().width / nt.textContent.length, bpr: V.bpr };
    });
    const bad = [];
    for (const k of [0, 1, 40, 90, g.bpr - 1]) {
      await p.mouse.click(g.left + g.G + (k + 0.5) * g.CW, g.top + 30);
      await p.waitForTimeout(80);
      const got = await p.evaluate(() => window.__gc.D().sel.start);
      if (got !== k) bad.push(`${k}\u2192${got}`);
    }
    if (bad.length) throw new Error('off by: ' + bad.join(', '));
  });
  await step('the bases keep clear of the floating strip', async () => {
    const r = await p.evaluate(() => {
      const strip = document.querySelector('#mapside').getBoundingClientRect();
      let worst = 0;
      for (const nt of document.querySelectorAll('#seqrows .nt')) worst = Math.max(worst, nt.getBoundingClientRect().right);
      return { strip: strip.left, worst };
    });
    if (r.worst > r.strip) throw new Error(`the bases reach ${Math.round(r.worst)} and the strip starts at ${Math.round(r.strip)}`);
  });
  await step('a primer is drawn as a primer, not as a feature', async () => {
    const n = await p.$$eval('#seqrows .primbox', es => es.length);
    if (!n) throw new Error('no primer boxes');
    const r = await p.evaluate(() => {
      const box = document.querySelector('#seqrows .primbox');
      const txt = box.textContent, d = window.__gc.D();
      const key = box.dataset.k, id = +box.dataset.prim;
      const pr = d.primers.find(x => x.id === id);
      const top = document.querySelector('#seqrows .nt').getBoundingClientRect();
      return { txt, name: pr && pr.name, bases: /^[ACGTN]+$/i.test(txt), above: box.getBoundingClientRect().top < top.top };
    });
    if (!r.bases) throw new Error('the box does not carry bases: ' + r.txt);
    if (!r.name) throw new Error('the box is not tied to a primer');
    // and it is not in the feature lanes
    const asFeature = await p.$$eval('#seqrows .fbar', (es, nm) => es.some(e => e.textContent.trim() === nm), r.name);
    if (asFeature) throw new Error('it is still drawn as a feature bar too');
  });
  await step('a primer box stands exactly over the bases it binds', async () => {
    const r = await p.evaluate(() => {
      const out = [];
      for (const box of document.querySelectorAll('#seqrows .primbox')) {
        const row = box.closest('.srow'), nts = row.querySelectorAll('.nt');
        const top = nts[0], bot = nts[1] || nts[0];
        const rb = box.getBoundingClientRect(), rt = top.getBoundingClientRect();
        const CW = rt.width / top.textContent.length;
        const at = (rb.left + 1 - rt.left) / CW;          // the 1px is the box's own border
        const idx = Math.round(at);
        const strand = rb.top > rt.top ? bot.textContent : top.textContent;
        out.push({
          off: Math.abs(at - idx),
          same: strand.substr(idx, box.textContent.length).toUpperCase() === box.textContent.toUpperCase()
        });
      }
      return out;
    });
    if (!r.length) throw new Error('no primer boxes to check');
    const slip = r.filter(x => x.off > 0.05);
    if (slip.length) throw new Error(`${slip.length} box(es) start between two bases`);
    const wrong = r.filter(x => !x.same);
    if (wrong.length) throw new Error(`${wrong.length} box(es) do not read as the strand under them`);
  });
  await step('a primer keeps the case it was typed in', async () => {
    await p.evaluate(() => {
      const d = window.__gc.D();
      d.primers.push({ id: 987654, name: 'case probe', seq: 'GTAAAACGACGGCCAGT', seqCase: 'gtaaaaCGACGGCCAGT' });
      d.cache = null;                       // the binding sites are worked out again
      window.__gc.renderAll();
    });
    await p.waitForTimeout(700);
    const shown = await p.$$eval('#seqrows .primbox', es => es.map(e => e.textContent));
    if (!shown.some(t => /[a-z]/.test(t) && /[A-Z]/.test(t))) throw new Error('the case was flattened: ' + shown.join(' | '));
    await p.evaluate(() => {
      const d = window.__gc.D();
      d.primers = d.primers.filter(x => x.id !== 987654);
      d.cache = null;
      window.__gc.renderAll();
    });
    await p.waitForTimeout(400);
  });
  await step('a translated feature carries a bar behind its residues', async () => {
    const n = await p.$$eval('#seqrows .aaband', es => es.length);
    if (!n) throw new Error('no band behind any amino acid row');
    const clickable = await p.evaluate(() => {
      const b = document.querySelector('#seqrows .aaband');
      return getComputedStyle(b).pointerEvents === 'none';
    });
    if (!clickable) throw new Error('the band swallows clicks meant for the residues');
    await p.evaluate(() => window.__gc.setView('map'));
    await p.waitForTimeout(400);
  });
  await step('the map draws primers as oligos, not as arrows', async () => {
    const n = await p.$$eval('#mapwrap g.prim', es => es.length);
    if (!n) throw new Error('no primers on the map');
    const parts = await p.evaluate(() => {
      const g2 = document.querySelector('#mapwrap g.prim');
      return { line: !!g2.querySelector('.primline'), hook: !!g2.querySelector('.primhook'), head: !!g2.querySelector('.primhead') };
    });
    if (!parts.line || !parts.hook || !parts.head) throw new Error('missing ' + JSON.stringify(parts));
  });

  // hovering an enzyme lights up every copy of it, with a card beside the name
  await step('enzyme hover', async () => {
    await p.evaluate(() => { const g = window.__gc, d = g.D(); d.chosen = ['EcoRI', 'AluI', 'SmaI']; d.opts.enzMode = 'chosen'; g.setView('map'); g.renderAll(); });
    await p.waitForTimeout(500);
    await p.hover('#mapwrap text.lbl.e[data-enz="AluI"]');
    await p.waitForTimeout(400);
    const all = await p.$$eval('#mapwrap [data-enz="AluI"]', e => e.length);
    const hot = await p.$$eval('#mapwrap [data-enz="AluI"].hov', e => e.length);
    if (!all || hot !== all) throw new Error(`${hot} of ${all} parts lit up`);
    const tip = (await p.textContent('#enztip')).replace(/\s+/g, ' ');
    if (!/AGCT/.test(tip) || !/cuts \d+ site/.test(tip)) throw new Error('card says: ' + tip.slice(0, 60));
    if (!await p.$$eval('#mapwrap text.lbl.e.bl', e => e.length)) throw new Error('no blunt enzyme is green');
    if (!await p.$$eval('#mapwrap text.lbl.e.u', e => e.length)) throw new Error('no single cutter is bold');
    await p.mouse.move(10, 950); await p.waitForTimeout(300);
    if (await p.$$eval('.hov', e => e.length)) throw new Error('the highlight did not clear');
  });
  await step('the cut drawn on the sequence', async () => {
    await p.evaluate(() => window.__gc.setView('seq'));
    await p.waitForTimeout(600);
    await p.hover('#seqrows .enz[data-enz="EcoRI"]');
    await p.waitForTimeout(400);
    const lines = await p.$$eval('#seqrows .cutline', es => es.map(e => Math.round(parseFloat(e.style.left))));
    if (lines.length !== 2 || lines[0] === lines[1]) throw new Error('5′ overhang not drawn staggered: ' + lines.join(','));
    await p.hover('#seqrows .enz[data-enz="SmaI"]');
    await p.waitForTimeout(400);
    const blunt = await p.$$eval('#seqrows .cutline', es => es.map(e => Math.round(parseFloat(e.style.left))));
    if (blunt.length !== 2 || blunt[0] !== blunt[1]) throw new Error('blunt cut not drawn level: ' + blunt.join(','));
    await p.mouse.move(10, 950); await p.waitForTimeout(300);
    await p.evaluate(() => window.__gc.setView('map'));
    await p.waitForTimeout(400);
  });
  await step('primers are purple, forward out and reverse in', async () => {
    const got = await p.$$eval('#mapwrap g.prim .primline', es => es.map(e => e.getAttribute('stroke')));
    if (!got.length || got.some(c => c !== '#7b3fb0' && c !== '#b48ad6')) throw new Error('primer colours: ' + got.join(','));
    // the forward ones ride outside the backbone, the reverse ones inside
    const sides = await p.evaluate(() => {
      const g2 = window.__gc, d = g2.D(), A = g2.analyze ? null : null;
      const box = document.querySelector('#mapwrap svg').getBoundingClientRect();
      return [...document.querySelectorAll('#mapwrap g.prim')].map(el => {
        const t = el.querySelector('title').textContent;
        const r = el.getBoundingClientRect();
        return { rev: /reverse/.test(t), d: Math.hypot(r.left + r.width / 2 - (box.left + box.width / 2), r.top + r.height / 2 - (box.top + box.height / 2)) };
      });
    });
    const fwd = sides.filter(x => !x.rev), rev = sides.filter(x => x.rev);
    if (fwd.length && rev.length && Math.min(...fwd.map(x => x.d)) <= Math.max(...rev.map(x => x.d))) {
      throw new Error('a reverse primer is further out than a forward one');
    }
  });
  await step('the pieces of a digest are listed', async () => {
    await p.evaluate(() => window.__gc.menuAct('Cloning', 'Restriction Cloning…'));
    await p.waitForTimeout(700);
    for (const nm of ['EcoRI', 'HindIII']) {
      await p.evaluate(e => { const ss = [...document.querySelectorAll('[data-venz]')]; const s2 = ss[ss.length - 1]; s2.value = e; s2.dispatchEvent(new Event('change', { bubbles: true })); }, nm);
      await p.waitForTimeout(450);
    }
    const rows = await p.$$eval('.asmfrags .asmlistitem', es => es.map(e => e.textContent.replace(/\s+/g, ' ').trim()));
    if (rows.length !== 2) throw new Error(rows.length + ' pieces listed');
    if (!await p.$('#asmMapBox circle.bb')) throw new Error('the circular vector is not drawn as a circle');
    const first = await p.$eval('#asmMapBox .sellbl', e => e.textContent);
    await p.click('.asmfrags .asmlistitem:nth-child(2) .aflead');
    await p.waitForTimeout(600);
    const second = await p.$eval('#asmMapBox .sellbl', e => e.textContent);
    if (first === second) throw new Error('choosing another piece changed nothing (' + first + ')');
    await p.click('.asmfrags .asmlistitem:nth-child(2) [data-useas^="f:"]');
    await p.waitForTimeout(700);
    if ((await p.$eval('.asmtab.on', e => e.textContent)) !== 'Fragment') throw new Error('"as Fragment" did not move to the Fragment tab');
    await p.keyboard.press('Escape'); await p.waitForTimeout(250);
  });

  // the one cloning window, and its method dropdown
  for (const item of ['Restriction Cloning…', 'Golden Gate Assembly…', 'Gibson Assembly…', 'In-Fusion Cloning…']) {
    await step('cloning window ' + item, async () => {
      await p.evaluate(i => window.__gc.menuAct('Cloning', i), item);
      await p.waitForTimeout(600);
      if (!await p.$('#asmWin')) throw new Error('did not open');
      if (!await p.$('#asmMethod')) throw new Error('no method dropdown');
      for (const tab of ['fragment', 'product', 'vector']) {
        await p.click(`[data-tab="${tab}"]`); await p.waitForTimeout(250);
      }
      await p.keyboard.press('Escape');
      await p.waitForTimeout(250);
      if (await p.$('#asmWin')) throw new Error('did not close');
    });
  }
  await step('cloning map has the side strip', async () => {
    await p.evaluate(() => window.__gc.menuAct('Cloning', 'Restriction Cloning…'));
    await p.waitForTimeout(600);
    const rows = await p.$$eval('.tmapside .msrow', e => e.length);
    if (rows !== 7) throw new Error('strip has ' + rows + ' rows');
    await p.click('.tmapside [data-mschev="enz"]'); await p.waitForTimeout(400);
    if (!await p.$('#tmspop [data-msv="enzMode"]')) throw new Error('no enzyme panel');
    const n0 = await p.$$eval('#asmMapBox text[data-pick^="e:"]', e => e.length);
    await p.selectOption('#tmspop [data-msv="enzMode"]', 'none'); await p.waitForTimeout(500);
    const n1 = await p.$$eval('#asmMapBox text[data-pick^="e:"]', e => e.length);
    if (!(n0 > 0 && n1 === 0)) throw new Error('enzymes did not hide: ' + n0 + ' -> ' + n1);
    await p.selectOption('#tmspop [data-msv="enzMode"]', 'unique6'); await p.waitForTimeout(500);
    // hide one primer through the strip
    await p.click('.tmapside [data-mschev="primer"]'); await p.waitForTimeout(400);
    const pr = await p.$$eval('#tmspop [data-msitem="primer"]', es => es.map(e => e.value));
    const d0 = await p.$$eval('#asmMapBox .prim', e => e.length);
    if (pr.length) {
      await p.uncheck(`#tmspop [data-msitem="primer"][value="${pr[0]}"]`); await p.waitForTimeout(500);
      const d1 = await p.$$eval('#asmMapBox .prim', e => e.length);
      if (d1 !== d0 - 1) throw new Error('primer did not hide: ' + d0 + ' -> ' + d1);
    }
    await p.keyboard.press('Escape'); await p.waitForTimeout(200);   // closes the panel
    await p.keyboard.press('Escape'); await p.waitForTimeout(250);   // closes the window
  });
  await step('simulation map has the side strip', async () => {
    await p.click('#toolbar [data-a="gel"]'); await p.waitForTimeout(800);
    const rows = await p.$$eval('#simBody .tmapside .msrow', e => e.length);
    if (rows !== 7) throw new Error('strip has ' + rows + ' rows');
    await p.click('#simBody .tmapside [data-ms="feat"]'); await p.waitForTimeout(500);
    if (await p.$$eval('#simBody .feat', e => e.length)) throw new Error('features did not hide');
    await p.click('#simBody .tmapside [data-ms="feat"]'); await p.waitForTimeout(500);
    if (!await p.$$eval('#simBody .feat', e => e.length)) throw new Error('features did not come back');
    await p.click('#doctabs .dtab:last-child .x'); await p.waitForTimeout(400);
  });
  await step('cloning window switches method', async () => {
    await p.evaluate(() => window.__gc.menuAct('Cloning', 'Restriction Cloning…'));
    await p.waitForTimeout(600);
    for (const m of ['goldengate', 'gibson', 'gateway', 'restriction']) {
      await p.selectOption('#asmMethod', m); await p.waitForTimeout(400);
      const got = await p.$eval('#asmMethod', e => e.value);
      if (got !== m) throw new Error('method stuck on ' + got);
    }
    await p.keyboard.press('Escape'); await p.waitForTimeout(250);
  });

  await p.evaluate(() => window.__gc.setSel(100, 400, { anchor: 100 }));
  for (const [m, i] of dialogs) {
    await step('dialog ' + m + ' > ' + i, async () => {
      await p.evaluate(([mm, ii]) => window.__gc.menuAct(mm, ii), [m, i]);
      await p.waitForTimeout(450);
      const open = await p.$('#dlg .hd');
      if (!open) throw new Error('did not open');
      await p.keyboard.press('Escape');
      await p.waitForTimeout(250);
      if (await p.$eval('#modal', e => !e.hidden)) throw new Error('did not close');
    });
  }

  await step('detect common features', async () => {
    await p.evaluate(() => window.__gc.menuAct('Features', 'Detect Common Features\u2026'));
    await p.waitForTimeout(400);
  });
  await step('mark att sites', async () => {
    await p.evaluate(() => window.__gc.menuAct('Cloning', 'Show att Sites in This Sequence'));
    await p.waitForTimeout(400);
  });
  await step('new agarose gel', async () => {
    await p.evaluate(() => window.__gc.menuAct('Analyze', 'New Agarose Gel\u2026'));
    await p.waitForTimeout(700);
    if (!(await p.$('#glImg svg'))) throw new Error('no gel drawn');
    await p.click('.dtab');
    await p.waitForTimeout(400);
  });
  await step('save genbank round-trip', async () => {
    const r = await p.evaluate(() => {
      const d = window.__gc.D();
      const gb = window.DNA.toGenBank(d, () => []);
      const back = window.DNA.parseGenBank(gb);
      return { same: back.seq === d.seq, circ: back.circular === d.circular, f: back.features.length, orig: d.features.length };
    });
    if (!r.same || !r.circ) throw new Error('round-trip mismatch ' + JSON.stringify(r));
    console.log('     round-trip features ' + r.f + '/' + r.orig);
  });

  await step('the View menu fits a short window', async () => {
    await p.setViewportSize({ width: 1200, height: 620 });
    await p.waitForTimeout(250);
    await p.click('.menubar .mi[data-i="2"]');
    await p.waitForTimeout(250);
    const m = await p.evaluate(() => {
      const el = document.querySelector('#menupop');
      const r = el.getBoundingClientRect();
      return { top: r.top, bottom: r.bottom, h: r.height, win: window.innerHeight,
               scrolls: el.scrollHeight > el.clientHeight + 1,
               groups: [...el.querySelectorAll('.it.has-sub')].map(x => x.textContent.trim()) };
    });
    if (m.bottom > m.win) throw new Error('it runs off the bottom: ' + JSON.stringify(m));
    for (const g of ['Map', 'Sequence', 'Show', 'Colours', 'Window', 'Map and Sequence Size']) {
      if (!m.groups.includes(g)) throw new Error('no ' + g + ' group: ' + m.groups.join(', '));
    }
    await p.keyboard.press('Escape');
    await p.setViewportSize({ width: 1280, height: 860 });
    await p.waitForTimeout(250);
  });

  await step('a long menu scrolls instead of running off', async () => {
    await p.setViewportSize({ width: 1200, height: 360 });
    await p.waitForTimeout(250);
    await p.click('.menubar .mi[data-i="0"]');            // File, the longest one
    await p.waitForTimeout(250);
    const m = await p.evaluate(() => {
      const el = document.querySelector('#menupop');
      const r = el.getBoundingClientRect();
      return { bottom: r.bottom, win: window.innerHeight, scrolls: el.scrollHeight > el.clientHeight + 1 };
    });
    if (m.bottom > m.win) throw new Error('it still runs off: ' + JSON.stringify(m));
    if (!m.scrolls) throw new Error('it was not given a scrollbar: ' + JSON.stringify(m));
    await p.keyboard.press('Escape');
    await p.setViewportSize({ width: 1280, height: 860 });
    await p.waitForTimeout(250);
  });

  await step('the interface size fits more on the screen, it does not magnify', async () => {
    await p.evaluate(() => window.__gc.setView('seq'));
    await p.waitForTimeout(500);
    const read = async f => {
      await p.evaluate(x => window.__gc.setUiScale(x, true), f);
      await p.waitForTimeout(450);
      return p.evaluate(() => {
        const nt = document.querySelector('.srow .nt');
        const bar = document.querySelector('.srow .fbar');
        return { px: parseFloat(getComputedStyle(nt).fontSize), perRow: nt.textContent.trim().length,
                 barH: bar ? Math.round(bar.getBoundingClientRect().height) : 0 };
      });
    };
    const small = await read(0.8), big = await read(1.25);
    if (!(small.px < big.px)) throw new Error('the text did not change size: ' + JSON.stringify([small, big]));
    if (!(small.perRow > big.perRow)) throw new Error('a smaller size has to fit more bases per row: ' + JSON.stringify([small, big]));
    if (!(small.barH < big.barH)) throw new Error('the feature bars did not follow: ' + JSON.stringify([small, big]));
    console.log(`     ${small.perRow} bases a row at 80%, ${big.perRow} at 125%`);
    await p.evaluate(() => window.__gc.setUiScale(1, true));
    await p.waitForTimeout(400);
  });

  await step('every primer points the way it reads', async () => {
    await p.evaluate(() => window.__gc.setView('seq'));
    await p.waitForTimeout(500);
    const arrows = await p.evaluate(() => [...document.querySelectorAll('.primarrow')].map(a => {
      // a primer split over several rows has one box per row; take the one beside this head
      const box = [...a.closest('.srow').querySelectorAll('.primbox')].find(x => x.dataset.prim === a.dataset.prim);
      const ar = a.getBoundingClientRect(), br = box.getBoundingClientRect();
      const top = box.closest('.srow').querySelector('.nt').getBoundingClientRect();
      return { right: a.className.trim().endsWith('r'), after: ar.left > br.left, bottomStrand: br.top > top.top,
               touching: Math.abs(a.className.trim().endsWith('r') ? ar.left - br.right : br.left - ar.right) < 2 };
    }));
    if (!arrows.length) throw new Error('no primer was given a head');
    for (const a of arrows) {
      if (a.bottomStrand === a.right) throw new Error('a head points the wrong way: ' + JSON.stringify(a));
      if (a.right !== a.after) throw new Error('a head is on the wrong end: ' + JSON.stringify(a));
      if (!a.touching) throw new Error('a head is adrift from its box: ' + JSON.stringify(a));
    }
    console.log('     ' + arrows.length + ' primer heads, each at its 3\u2032 end');
  });

  await step('the enzyme reference opens on both of its tabs', async () => {
    await p.evaluate(() => window.__gc.menuAct('Enzymes', 'Enzyme Database\u2026'));
    await p.waitForTimeout(800);
    const sups = await p.$$eval('#esSup option', o => o.length);
    if (sups < 10) throw new Error('only ' + sups + ' suppliers in the picker');
    // a supplier that publishes buffer figures shows one column per buffer
    await p.selectOption('#esSup', 'Yeasen').catch(() => {});
    await p.waitForTimeout(400);
    const row = await p.$$eval('.esgrid tbody tr:first-child td', td => td.map(x => x.textContent.trim()));
    if (row.length < 6) throw new Error('no activity columns: ' + JSON.stringify(row));
    if (!/^\d+%$/.test(row[4])) throw new Error('the first activity is not a percentage: ' + JSON.stringify(row));
    if (!/\u00b0C$/.test(row[2])) throw new Error('no incubation temperature: ' + JSON.stringify(row));
    // the filter narrows it
    await p.fill('#esQ', 'BamHI');
    await p.waitForTimeout(350);
    const n = await p.$$eval('.esgrid tbody tr', r => r.length);
    if (!n) throw new Error('the filter found no BamHI');
    await p.fill('#esQ', '');
    await p.click('[data-t="byenz"]');
    await p.waitForTimeout(400);
    const rows2 = await p.$$eval('.esgrid tbody tr', r => r.length);
    if (rows2 < 500) throw new Error('the supplier list is short: ' + rows2);
    console.log('     ' + sups + ' suppliers, ' + rows2 + ' enzymes listed');
    await p.keyboard.press('Escape');
    await p.waitForTimeout(300);
  });

  await step('every enzyme cuts where a palindrome says it must', async () => {
    const r = await p.evaluate(() => {
      const E = window.DNA.ENZYMES;
      // a palindromic site has to be cut symmetrically; anything else is a typo in the table
      const lop = E.filter(e => window.DNA.rc(e.site) === e.site && e.bot !== e.site.length - e.top)
        .map(e => `${e.name} ${e.site} (${e.top}/${e.bot})`);
      const odd = E.filter(e => !(e.top >= 0 && e.bot >= 0 && e.top <= e.site.length + 30 && e.bot <= e.site.length + 30))
        .map(e => e.name);
      const dup = E.length - new Set(E.map(e => e.name)).size;
      return { total: E.length, common: E.filter(e => e.common).length, lop, odd, dup };
    });
    if (r.lop.length) throw new Error('cut not mirrored on a palindrome: ' + r.lop.join(', '));
    if (r.odd.length) throw new Error('cut position out of range: ' + r.odd.join(', '));
    if (r.dup) throw new Error(r.dup + ' enzymes share a name');
    if (r.common >= r.total) throw new Error('the catalogue additions are not marked');
    console.log(`     ${r.total} enzymes, ${r.common} of them in common use`);
  });

  await step('the catalogue stays out of the way until it is asked for', async () => {
    const r = await p.evaluate(() => {
      const d = window.__gc.D();
      const sites = window.DNA.findSites(d.seq, d.circular);
      const u6 = wide => window.DNA.ENZYMES.filter(e => {
        const s = sites[e.name];
        if (!s || (!wide && !e.common)) return false;
        return s.length === 1 && e.core >= 6;
      }).length;
      return { narrow: u6(false), wide: u6(true) };
    });
    if (!(r.wide > r.narrow)) throw new Error('the catalogue adds nothing: ' + JSON.stringify(r));
    if (r.narrow > 60) throw new Error('too many labels drawn by default: ' + r.narrow);
    console.log(`     ${r.narrow} unique 6+ cutters drawn, ${r.wide} with the whole catalogue`);
  });

  await step('an enzyme the catalogue lists knows who sells it', async () => {
    const r = await p.evaluate(() => {
      // the supplier tables are what fills e.suppliers, so they are what to check against
      const listed = new Set();
      for (const s of window.ENZYME_DB.suppliers) for (const x of s.enzymes) listed.add(window.DNA.baseEnzymeName(x.name));
      const bad = window.DNA.ENZYMES.filter(e => listed.has(e.name) && !(e.suppliers || []).length).map(e => e.name);
      const unlisted = window.DNA.ENZYMES.filter(e => !listed.has(e.name)).map(e => e.name);
      let wrong = 0;
      for (const e of window.DNA.ENZYMES) for (const a of e.aliases || []) {
        if (window.DNA.baseEnzymeName(a) !== e.name) wrong++;
      }
      const smaI = window.DNA.ENZ_BY_NAME['SmaI'];
      return { bad, unlisted, wrong, neo: (smaI.aliases || []).some(a => /XmaI/.test(a)) };
    });
    if (r.bad.length) throw new Error(r.bad.length + ' listed enzymes have no supplier: ' + r.bad.slice(0, 5).join(', '));
    if (r.wrong) throw new Error(r.wrong + ' branded names hang off the wrong enzyme');
    if (r.neo) throw new Error('SmaI took XmaI as one of its own names, but they cut the same site differently');
    console.log(`     every listed enzyme has its suppliers; ${r.unlisted.length} of ours the catalogue does not list (${r.unlisted.join(', ')})`);
  });

  await step('no label on the map lands on another', async () => {
    // a crowded sequence, so the labels have to compete for room
    await p.evaluate(() => {
      const d = window.__gc.D(), n = d.seq.length;
      d.features = Array.from({ length: 24 }, (_, i) => ({
        id: 'sf' + i, name: 'a feature with a long name ' + (i + 1), type: i % 3 ? 'misc_feature' : 'CDS',
        start: Math.round(n * i / 24), end: Math.round(n * (i + 0.7) / 24) % n,
        strand: i % 2 ? 1 : -1, color: '#6b8ed6', translated: i % 3 === 0
      }));
      d.primers = Array.from({ length: 18 }, (_, i) => {
        const at = Math.round(n * i / 18);
        const part = window.DNA.sub(d.seq, at, at + 20);
        // every third one reads off the bottom strand, so both bands are exercised
        return { id: 'sp' + i, name: 'screening primer number ' + (i + 1),
          seq: i % 3 === 2 ? window.DNA.rc(part) : part };
      });
      d.opts.showLocations = true;
      d.cache = null;
      window.__gc.setView('map');
    });
    for (const shape of ['circular', 'linear']) {
      await p.evaluate(m => { window.__gc.D().opts.mapMode = m; window.__gc.renderAll(true); }, shape);
      await p.waitForTimeout(800);
      const r = await p.evaluate(() => {
        const svg = document.querySelector('#view .mapwrap svg') || document.querySelector('#view svg');
        const t = [...svg.querySelectorAll('text')]
          .map(e => { const b = e.getBoundingClientRect(); return { s: e.textContent.trim(), x: b.x, y: b.y, w: b.width, h: b.height }; })
          .filter(e => e.w > 0.5 && e.h > 0.5);
        const hit = [];
        for (let i = 0; i < t.length; i++) for (let j = i + 1; j < t.length; j++) {
          const a = t[i], c = t[j];
          if (Math.min(a.x + a.w, c.x + c.w) - Math.max(a.x, c.x) > 1 &&
              Math.min(a.y + a.h, c.y + c.h) - Math.max(a.y, c.y) > 1) hit.push(`"${a.s}" over "${c.s}"`);
        }
        return { texts: t.length, hit };
      });
      if (r.hit.length) throw new Error(shape + ': ' + r.hit.length + ' overlapping, e.g. ' + r.hit.slice(0, 3).join('; '));
      console.log(`     ${shape}: ${r.texts} labels, none touching`);
    }
  });

  await step('features lie against the sequence and primers outside them', async () => {
    await p.evaluate(() => { window.__gc.setView('seq'); window.__gc.renderAll(true); });
    await p.waitForTimeout(800);
    const r = await p.evaluate(() => {
      const out = [];
      for (const row of document.querySelectorAll('.srow')) {
        const base = row.getBoundingClientRect().top;
        const low = sel => {
          let v = null;
          for (const e of row.querySelectorAll(sel)) { const t = e.getBoundingClientRect().top - base; if (v === null || t < v) v = t; }
          return v;
        };
        const strand = low('.nt'), feat = low('.fbar');
        let rev = null;
        for (const e of row.querySelectorAll('.primbox')) {
          const t = e.getBoundingClientRect().top - base;
          if (strand !== null && t > strand && (rev === null || t < rev)) rev = t;
        }
        if (feat !== null && rev !== null) out.push({ feat, rev });
      }
      return out;
    });
    if (!r.length) throw new Error('no row held both a feature and a bottom-strand primer');
    const bad = r.filter(x => x.feat >= x.rev);
    if (bad.length) throw new Error(bad.length + ' rows put a primer nearer the sequence than a feature');
    console.log(`     ${r.length} rows, the feature nearer in every one`);
  });

  await step('a name too wide for its feature is written beside it', async () => {
    await p.evaluate(() => {
      const d = window.__gc.D();
      d.features = [{ id: 'tiny', name: 'a name far wider than its feature', type: 'misc_feature',
        start: 320, end: 328, strand: 1, color: '#c45b9b' }];
      d.cache = null;
    });
    const read = async mode => {
      await p.evaluate(m => { window.__gc.D().opts.labelMode = m; window.__gc.renderAll(true); }, mode);
      await p.waitForTimeout(600);
      return p.evaluate(() => {
        const bar = [...document.querySelectorAll('.fbar')].find(e => (e.getAttribute('title') || '').includes('a name far wider'));
        const out = [...document.querySelectorAll('.fname')].find(e => e.textContent.includes('a name far wider'));
        return { inBar: bar ? bar.textContent.trim() : null, beside: !!out };
      });
    };
    const auto = await read('auto'), inside = await read('inside'), outside = await read('outside');
    if (auto.inBar || !auto.beside) throw new Error('auto should put it beside the bar: ' + JSON.stringify(auto));
    if (!inside.inBar || inside.beside) throw new Error('internal only should keep it in the bar: ' + JSON.stringify(inside));
    if (outside.inBar || !outside.beside) throw new Error('external only should put it beside: ' + JSON.stringify(outside));
    await p.evaluate(() => { window.__gc.D().opts.labelMode = 'auto'; });
    console.log('     inside, beside and clipped all follow the setting');
  });

  await step('a primer name keeps clear of its own box', async () => {
    await p.evaluate(() => { const d = window.__gc.D(); d.features = []; d.cache = null; });
    for (const sc of [0.8, 1, 1.25]) {
      await p.evaluate(x => window.__gc.setUiScale(x, true), sc);
      await p.waitForTimeout(650);
      const gaps = await p.evaluate(() => {
        const out = [];
        for (const nm of document.querySelectorAll('.primname')) {
          const row = nm.closest('.srow');
          const box = [...row.querySelectorAll('.primbox')].find(x => x.dataset.prim === nm.dataset.prim);
          if (!box) continue;
          const a = nm.getBoundingClientRect(), c = box.getBoundingClientRect();
          out.push(a.top > c.top ? a.top - c.bottom : c.top - a.bottom);
        }
        return out;
      });
      if (!gaps.length) throw new Error('no primer name was drawn');
      const tight = gaps.filter(g => g < 2);
      if (tight.length) throw new Error(`at ${sc}, ${tight.length} names sit within 2px of their box`);
    }
    await p.evaluate(() => window.__gc.setUiScale(1, true));
    await p.waitForTimeout(500);
    console.log('     clear at 80%, 100% and 125%');
  });

  await step('the longest translated feature stacks nearest the sequence', async () => {
    await p.evaluate(() => {
      const d = window.__gc.D();
      d.features = [
        { id: 'a', name: 'short early', type: 'misc_feature', start: 100, end: 300, strand: 1, color: '#c45b9b' },
        { id: 'b', name: 'long plain', type: 'misc_feature', start: 100, end: 900, strand: 1, color: '#4fa86a' },
        { id: 'c', name: 'long translated', type: 'CDS', start: 100, end: 900, strand: 1, color: '#e8a33d', translated: true },
        { id: 'd', name: 'short late', type: 'misc_feature', start: 100, end: 300, strand: 1, color: '#6b8ed6' }
      ];
      d.cache = null;
      window.__gc.renderAll(true);
    });
    await p.waitForTimeout(700);
    const order = () => p.evaluate(() => {
      const row = [...document.querySelectorAll('.srow')].find(r => r.querySelectorAll('.fbar').length > 2);
      return row ? [...row.querySelectorAll('.fbar')]
        .map(e => ({ t: (e.getAttribute('title') || '').replace(/ \d.*$/, ''), y: e.getBoundingClientRect().top }))
        .sort((a, b) => a.y - b.y).map(x => x.t) : [];
    });
    const want = ['long translated', 'long plain', 'short early', 'short late'];
    const got = await order();
    if (got.join('|') !== want.join('|')) throw new Error('default order is ' + got.join(' -> '));
    await p.evaluate(() => {
      const d = window.__gc.D(), f = d.features.find(x => x.id === 'd');
      for (let i = 0; i < 3; i++) window.__gc.stackMove(f, -1);
    });
    await p.waitForTimeout(700);
    const moved = await order();
    if (moved[0] !== 'short late') throw new Error('moving it nearer did nothing: ' + moved.join(' -> '));
    await p.evaluate(() => window.__gc.stackReset());
    await p.waitForTimeout(700);
    const back = await order();
    if (back.join('|') !== want.join('|')) throw new Error('resetting did not restore the order: ' + back.join(' -> '));
    console.log('     ' + want.join(' -> ') + ', and a move sticks until reset');
  });

  console.log('\nERRORS:', errs.length);
  errs.forEach(e => console.log('  ' + e));
  await b.close();
  process.exit(errs.length ? 1 : 0);
})();
