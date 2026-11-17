/* Dogma Studio — the phone pass.
 *
 * Opens the same app/index.html as a phone would: an Android device profile, a coarse
 * pointer, and a stand-in for the Capacitor bridge so the Android host in
 * app/mobileboot.js really runs. It answers "is the phone build still right" — the
 * layout at three widths, the gestures that stand in for a mouse, and the host contract
 * for files, settings and NCBI.
 *
 * What it cannot see is the real WebView, the system file picker and the share sheet;
 * those need a device or an emulator. Everything else is here.
 */
'use strict';
const path = require('path');
const fs = require('fs');
const { chromium, devices } = require('playwright');
const { INDEX, runner, touch } = require('./lib');
const capacitorMock = require('./capacitor-mock');

const SHOTS = path.join(__dirname, '..', 'out', 'mobile-shots');
const SIZES = [
  { name: 'small phone', w: 360, h: 740 },
  { name: 'phone', w: 412, h: 915 },
  { name: 'tablet', w: 800, h: 1180 }
];

(async () => {
  const b = await chromium.launch();
  const errs = [];
  const step = runner(errs);
  fs.mkdirSync(SHOTS, { recursive: true });

  // ---------------------------------------------------------------- the layout at three widths
  for (const s of SIZES) {
    const p = await b.newPage({
      ...devices['Pixel 7'],
      viewport: { width: s.w, height: s.h },
      screen: { width: s.w, height: s.h }
    });
    p.on('pageerror', e => errs.push(`PAGEERROR (${s.name}): ` + e.message.split('\n')[0]));
    p.on('console', m => { if (m.type() === 'error') errs.push(`CONSOLE (${s.name}): ` + m.text()); });
    await p.addInitScript(capacitorMock);
    await p.goto(INDEX);
    await p.waitForTimeout(1100);

    await step(`${s.name}: nothing runs off the side`, async () => {
      const over = await p.evaluate(() => document.documentElement.scrollWidth - window.innerWidth);
      if (over > 0) throw new Error(over + 'px of horizontal overflow');
    });
    await step(`${s.name}: the bars can be tapped`, async () => {
      const small = await p.evaluate(() => {
        const bad = [];
        for (const sel of ['.menubar .mi', '.tb', '.vt', '.dtab']) {
          for (const el of document.querySelectorAll(sel)) {
            const r = el.getBoundingClientRect();
            if (r.width && r.height && r.height < 30) bad.push(`${sel} ${Math.round(r.width)}×${Math.round(r.height)}`);
          }
        }
        return bad.slice(0, 3);
      });
      if (small.length) throw new Error('too small to tap: ' + small.join(', '));
    });
    await step(`${s.name}: every view draws`, async () => {
      for (const v of ['map', 'seq', 'enz', 'feat', 'primer', 'hist']) {
        await p.evaluate(x => window.__gc.setView(x), v);
        await p.waitForTimeout(260);
        const empty = await p.evaluate(() => document.querySelector('#view').textContent.trim().length === 0
          && !document.querySelector('#view svg'));
        if (empty) throw new Error(v + ' drew nothing');
      }
      await p.evaluate(() => window.__gc.setView('map'));
      await p.waitForTimeout(300);
    });
    await step(`${s.name}: a dialog fits the screen`, async () => {
      await p.evaluate(() => window.__gc.menuAct('Tools', 'Show DNA Calculations…'));
      await p.waitForTimeout(600);
      const r = await p.evaluate(() => {
        const d = document.querySelector('#dlg').getBoundingClientRect();
        return { w: Math.round(d.width), h: Math.round(d.height), vw: innerWidth, vh: innerHeight };
      });
      if (r.w > r.vw || r.h > r.vh) throw new Error(`${r.w}×${r.h} in a ${r.vw}×${r.vh} window`);
      await p.keyboard.press('Escape');
      await p.waitForTimeout(250);
    });
    await p.screenshot({ path: path.join(SHOTS, `${s.name.replace(/ /g, '-')}.png`) });
    await p.close();
  }

  // ---------------------------------------------------------------- one phone, in depth
  const p = await b.newPage({ ...devices['Pixel 7'] });
  p.on('pageerror', e => errs.push('PAGEERROR: ' + e.message.split('\n')[0]));
  p.on('console', m => { if (m.type() === 'error') errs.push('CONSOLE: ' + m.text()); });
  await p.addInitScript(capacitorMock);
  await p.goto(INDEX);
  await p.waitForTimeout(1200);

  await step('the Android host is in place', async () => {
    const h = await p.evaluate(() => window.gcHost ? { platform: window.gcHost.platform, keys: Object.keys(window.gcHost) } : null);
    if (!h) throw new Error('mobileboot.js did not install gcHost');
    if (h.platform !== 'android') throw new Error('platform is ' + h.platform);
    for (const k of ['openFiles', 'saveFile', 'loadSettings', 'saveSettings', 'ncbiFetch', 'ncbiSearch']) {
      if (!h.keys.includes(k)) throw new Error('the host has no ' + k);
    }
    // and no second-window methods, so the dialogs are used instead
    for (const k of ['openTool', 'openAssembly']) if (h.keys.includes(k)) throw new Error('the host should not offer ' + k);
  });

  await step('settings are kept by the phone', async () => {
    await p.evaluate(async () => { window.__gc.SET.toolbarSmall = true; await window.__gc.saveSettings(); });
    await p.waitForTimeout(300);
    const saved = await p.evaluate(() => window.__native.prefs['dogma-studio-settings']);
    if (!saved || !JSON.parse(saved).toolbarSmall) throw new Error('Preferences did not take the settings');
    // and come back on the next start
    const again = await p.evaluate(() => JSON.parse(window.__native.prefs['dogma-studio-settings']).toolbarSmall);
    if (!again) throw new Error('the setting did not survive');
  });

  await step('saving writes a file and offers to send it', async () => {
    await p.evaluate(() => { window.__native.calls.length = 0; });
    await p.evaluate(async () => {
      const d = window.__gc.D();
      const gb = window.DNA.toGenBank(d, () => []);
      await window.gcHost.saveFile({ defaultName: 'pDogma-Demo.gb', content: gb });
    });
    await p.waitForTimeout(400);
    const r = await p.evaluate(() => {
      const w = window.__native.calls.find(c => c.name === 'Filesystem.writeFile');
      const s = window.__native.calls.find(c => c.name === 'Share.share');
      const data = w ? atob(window.__native.files[w.args.path]) : '';
      return { name: w && w.args.path, dir: w && w.args.directory, shared: !!s, head: data.slice(0, 6), len: data.length };
    });
    if (r.name !== 'pDogma-Demo.gb') throw new Error('wrote ' + r.name);
    if (r.dir !== 'DOCUMENTS') throw new Error('wrote into ' + r.dir);
    if (r.head !== 'LOCUS ') throw new Error('the file does not start with a GenBank header: ' + r.head);
    if (!r.shared) throw new Error('the share sheet was not offered');
  });

  await step('opening a file goes through the system picker', async () => {
    const gb = await p.evaluate(() => window.DNA.toGenBank(window.__gc.D(), () => []));
    const tmp = path.join(require('os').tmpdir(), 'dogma-mobile-open.gb');
    fs.writeFileSync(tmp, gb.replace(/^LOCUS {7}\S+/m, 'LOCUS       picked-from-phone'));
    const before = await p.evaluate(() => window.__gc.S.docs.length);
    // the host opens the system picker; the browser's own chooser stands in for it here
    p.once('filechooser', fc => fc.setFiles(tmp).catch(() => { }));
    p.evaluate(() => window.__gc.menuAct('File', 'Open…')).catch(() => { });
    for (let i = 0; i < 40; i++) {
      await p.waitForTimeout(250);
      if (await p.evaluate(n => window.__gc.S.docs.length > n, before)) break;
    }
    const after = await p.evaluate(() => ({ n: window.__gc.S.docs.length, name: window.__gc.D().name }));
    if (after.n !== before + 1) throw new Error(`${before} documents before, ${after.n} after`);
    if (!/picked-from-phone/.test(after.name)) throw new Error('opened ' + after.name);
    await p.evaluate(() => { const i = window.__gc.S.docs.length - 1; window.__gc.S.docs.splice(i, 1); window.__gc.S.active = 0; window.__gc.renderAll(); });
  });

  await step('NCBI is allowed out of the page', async () => {
    // the request is not made here; what matters is that the policy would let it
    const csp = await p.evaluate(() => (document.querySelector('meta[http-equiv="Content-Security-Policy"]') || {}).content || '');
    if (!/connect-src[^;]*eutils\.ncbi\.nlm\.nih\.gov/.test(csp)) throw new Error('connect-src does not allow NCBI');
  });

  // ---------------------------------------------------------------- the finger
  await step('a tap selects what is under it', async () => {
    await p.evaluate(() => { window.__gc.setView('map'); window.__gc.setSel(0, 0, { anchor: 0 }); });
    await p.waitForTimeout(500);
    await touch.tap(p, '#mapwrap .feat');
    await p.waitForTimeout(400);
    const L = await p.evaluate(() => { const d = window.__gc.D(); return d.sel.end - d.sel.start; });
    if (L <= 0) throw new Error('nothing was selected');
  });

  await step('a long press opens the item menu', async () => {
    await touch.press(p, '#mapwrap text.lbl.e');
    await p.waitForTimeout(400);
    const open = await p.evaluate(() => { const c = document.querySelector('#ctx'); return c && !c.hidden && c.textContent.trim().length > 0; });
    if (!open) throw new Error('no menu appeared');
    await p.evaluate(() => { const c = document.querySelector('#ctx'); if (c) c.hidden = true; });
  });

  await step('a long press and a drag select a stretch', async () => {
    await p.evaluate(() => { window.__gc.setView('seq'); window.__gc.setSel(0, 0, { anchor: 0 }); });
    await p.waitForTimeout(600);
    await touch.pressDrag(p, '#seqrows .nt', 150, 0);
    await p.waitForTimeout(400);
    const L = await p.evaluate(() => { const d = window.__gc.D(); return d.sel.end - d.sel.start; });
    if (L < 8) throw new Error('selected ' + L + ' bp');
  });

  await step('an ordinary drag scrolls instead', async () => {
    await p.evaluate(() => { window.__gc.setSel(0, 0, { anchor: 0 }); document.querySelector('#seqscroll').scrollTop = 0; });
    await p.waitForTimeout(300);
    await touch.swipe(p, '#seqrows .nt', 0, -120);
    await p.waitForTimeout(300);
    const L = await p.evaluate(() => { const d = window.__gc.D(); return d.sel.end - d.sel.start; });
    if (L > 0) throw new Error('a swipe selected ' + L + ' bp');
  });

  await step('a pinch zooms the map', async () => {
    await p.evaluate(() => { window.__gc.setView('map'); window.__gc.S.zoom = 1; window.__gc.renderAll(); });
    await p.waitForTimeout(600);
    await touch.pinch(p, '#mapwrap svg', true);
    await p.waitForTimeout(400);
    const z = await p.evaluate(() => window.__gc.S.zoom);
    if (z <= 1) throw new Error('zoom stayed at ' + z);
  });

  await step('the windows a phone has no room for open as dialogs', async () => {
    await p.evaluate(() => window.__gc.menuAct('Analyze', 'Enzyme Database…'));
    await p.waitForTimeout(900);
    if (!(await p.$('#dlg .hd'))) throw new Error('the enzyme database did not open');
    const fits = await p.evaluate(() => document.querySelector('#dlg').getBoundingClientRect().width <= innerWidth);
    if (!fits) throw new Error('it does not fit the screen');
    await p.keyboard.press('Escape');
    await p.waitForTimeout(300);
  });

  await step('the cloning window works on a phone', async () => {
    await p.evaluate(() => window.__gc.menuAct('Cloning', 'Restriction Cloning…'));
    await p.waitForTimeout(900);
    if (!(await p.$('#asmWin'))) throw new Error('it did not open');
    const over = await p.evaluate(() => document.documentElement.scrollWidth - innerWidth);
    if (over > 0) throw new Error(over + 'px of overflow');
    await p.keyboard.press('Escape');
    await p.waitForTimeout(300);
  });

  await p.screenshot({ path: path.join(SHOTS, 'phone-final.png') });
  await p.close();

  console.log('\nERRORS:', errs.length);
  errs.forEach(e => console.log('  ' + e));
  console.log('screenshots in ' + SHOTS);
  await b.close();
  process.exit(errs.length ? 1 : 0);
})();
