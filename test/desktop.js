/* GeneCanvas — the desktop pass.
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
    ['Help', 'About GeneCanvas'],
    ['Edit', 'Go To Position…'],
    ['Features', 'Add Feature…'],
    ['Primers', 'Add Primer…']
  ];
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
    const got = await p.$$eval('#mapwrap path.prim', es => es.map(e => e.getAttribute('fill')));
    if (!got.length || got.some(c => c !== '#7b3fb0' && c !== '#b48ad6')) throw new Error('primer colours: ' + got.join(','));
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

  console.log('\nERRORS:', errs.length);
  errs.forEach(e => console.log('  ' + e));
  await b.close();
  process.exit(errs.length ? 1 : 0);
})();
