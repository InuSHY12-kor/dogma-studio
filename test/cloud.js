/* GeneCanvas — what the optional account does with two machines.
 *
 * The account's job is to merge three libraries between places without losing an edit
 * or resurrecting a deletion. That is logic, not network, so this pass runs it against
 * a stand-in server held in memory: two "machines" sync against the same rows and the
 * result has to be what a person would expect.
 */
'use strict';
const path = require('path');
const { chromium } = require('playwright');
const { INDEX, runner } = require('./lib');

/* The table lives here, in node, so both pages write to the same one. In the page the
 * only stand-in is fetch: it answers the handful of requests app/cloud.js makes and
 * hands the rows to this side, which sets updated_at the way the real table's trigger
 * does — later for every write, so "who edited last" is decided by the server. */
const rows = [];
let clock = Date.now();
function dbGet() { return rows; }
function dbUpsert(list) {
  for (const row of list) {
    clock = Math.max(clock, Date.now()) + 1000;
    const rec = Object.assign({}, row, { updated_at: new Date(clock).toISOString() });
    const at = rows.findIndex(r => r.user_id === rec.user_id && r.kind === rec.kind && r.name === rec.name);
    if (at >= 0) rows[at] = rec; else rows.push(rec);
  }
  return true;
}

function pageStub() {
  window.fetch = async (url, opts) => {
    const u = String(url), o = opts || {};
    const reply = (body, status) => ({
      ok: (status || 200) < 400, status: status || 200, statusText: 'x',
      text: async () => (body === null ? '' : JSON.stringify(body))
    });
    if (u.includes('/auth/v1/token') || u.includes('/auth/v1/signup')) {
      const b = JSON.parse(o.body || '{}');
      return reply({ access_token: 'tok', refresh_token: 'ref', user: { id: 'u1', email: b.email } });
    }
    if (u.includes('/rest/v1/genecanvas_items')) {
      if ((o.method || 'GET') === 'GET') return reply(await window.__dbGet());
      await window.__dbUpsert(JSON.parse(o.body));
      return reply(null, 201);
    }
    return reply({ message: 'unexpected ' + u }, 404);
  };
}

(async () => {
  const b = await chromium.launch();
  const errs = [];
  const step = runner(errs);
  const ctx = await b.newContext({ viewport: { width: 1280, height: 860 } });

  const machine = async name => {
    const p = await ctx.newPage();
    p.on('pageerror', e => errs.push(`PAGEERROR (${name}): ` + e.message.split('\n')[0]));
    p.on('console', m => { if (m.type() === 'error') errs.push(`CONSOLE (${name}): ` + m.text()); });
    await p.exposeFunction('__dbGet', dbGet);
    await p.exposeFunction('__dbUpsert', dbUpsert);
    await p.addInitScript(pageStub);
    await p.goto(INDEX);
    await p.waitForTimeout(900);
    // each machine keeps its own settings, as two computers would
    await p.evaluate(n => { window.__gc.SET.enzymeSets = {}; window.__gc.SET.cloud = { items: {}, tombs: {}, lastSync: 0, email: '' }; }, name);
    return p;
  };

  const a = await machine('A');
  const c = await machine('B');

  const signIn = p => p.evaluate(async () => { await window.CLOUD.signIn('two@example.com', 'secret'); });
  const syncOf = p => p.evaluate(async () => {
    const s = await window.CLOUD.sync(window.__gc.SET);
    await window.__gc.saveSettings();
    return s;
  });
  const sets = p => p.evaluate(() => Object.keys(window.__gc.SET.enzymeSets || {}).sort());
  const addSet = (p, name, list) => p.evaluate(async ([n, l]) => {
    window.__gc.SET.enzymeSets[n] = l;
    await window.__gc.saveSettings();
  }, [name, list]);
  const dropSet = (p, name) => p.evaluate(async n => {
    delete window.__gc.SET.enzymeSets[n];
    await window.__gc.saveSettings();
  }, name);

  await step('signing in does not touch what is already here', async () => {
    await addSet(a, 'My Favourites', ['EcoRI', 'HindIII']);
    await signIn(a);
    const after = await sets(a);
    if (!after.includes('My Favourites')) throw new Error('the local set disappeared: ' + after.join(','));
  });

  await step('a set made here turns up on the other machine', async () => {
    const s1 = await syncOf(a);
    if (!s1.sent) throw new Error('nothing was sent');
    await signIn(c);
    await syncOf(c);
    const there = await sets(c);
    if (!there.includes('My Favourites')) throw new Error('the other machine has ' + there.join(','));
  });

  await step('an edit on the other machine comes back', async () => {
    await addSet(c, 'My Favourites', ['EcoRI', 'HindIII', 'BamHI']);
    await syncOf(c);
    await syncOf(a);
    const list = await a.evaluate(() => window.__gc.SET.enzymeSets['My Favourites']);
    if (list.length !== 3) throw new Error('this machine still has ' + list.join(','));
  });

  await step('each machine keeps what only it has', async () => {
    await addSet(a, 'Cloning kit', ['BsaI']);
    await addSet(c, 'Screening', ['AluI']);
    await syncOf(a); await syncOf(c); await syncOf(a);
    const here = await sets(a), there = await sets(c);
    for (const n of ['My Favourites', 'Cloning kit', 'Screening']) {
      if (!here.includes(n)) throw new Error('this machine lost ' + n);
      if (!there.includes(n)) throw new Error('the other machine lost ' + n);
    }
  });

  await step('a deletion travels instead of coming back', async () => {
    await dropSet(a, 'Cloning kit');
    await syncOf(a);
    await syncOf(c);
    if ((await sets(c)).includes('Cloning kit')) throw new Error('it came back on the other machine');
    await syncOf(a);
    if ((await sets(a)).includes('Cloning kit')) throw new Error('it came back here');
  });

  await step('primers and my own features travel too', async () => {
    await a.evaluate(async () => {
      const S = window.__gc.SET;
      S.primerDb = S.primerDb || { primers: [], labels: [] };
      S.primerDb.primers.push({ id: 1, name: 'T7 fwd', seq: 'TAATACGACTCACTATAGGG', labels: ['stock'] });
      S.primerDb.labels = ['stock'];
      S.myFeatures = S.myFeatures || [];
      S.myFeatures.push({ name: 'my tag', type: 'misc_feature', seq: 'GGATCCGAA', color: '#ff8800' });
      await window.__gc.saveSettings();
    });
    await syncOf(a);
    await syncOf(c);
    const got = await c.evaluate(() => {
      const S = window.__gc.SET;
      return {
        primers: (S.primerDb.primers || []).map(p => p.name),
        labels: (S.primerDb.labels || []).map(l => typeof l === 'string' ? l : l.name),
        feats: (S.myFeatures || []).map(f => f.name)
      };
    });
    if (!got.primers.includes('T7 fwd')) throw new Error('primers: ' + got.primers.join(','));
    if (!got.labels.includes('stock')) throw new Error('labels: ' + got.labels.join(','));
    if (!got.feats.includes('my tag')) throw new Error('features: ' + got.feats.join(','));
  });

  await step('an edit on both machines keeps both copies', async () => {
    await addSet(a, 'Shared', ['EcoRI']);
    await syncOf(a); await syncOf(c);
    await addSet(a, 'Shared', ['EcoRI', 'SalI']);        // here
    await addSet(c, 'Shared', ['EcoRI', 'XbaI']);        // and there, before either syncs
    await syncOf(c);                                      // the other machine gets in first
    const st = await syncOf(a);
    const here = await a.evaluate(() => ({
      shared: window.__gc.SET.enzymeSets['Shared'],
      mine: window.__gc.SET.enzymeSets['Shared (this device)']
    }));
    if (!here.shared || !here.shared.includes('XbaI')) throw new Error('the account copy was not taken: ' + JSON.stringify(here.shared));
    if (!here.mine || !here.mine.includes('SalI')) throw new Error('what this machine had was lost: ' + JSON.stringify(here.mine));
    if (!st.kept) throw new Error('the sync did not report keeping a copy');
  });

  await step('signing out leaves everything on the machine', async () => {
    await a.evaluate(() => window.CLOUD.signOut());
    const here = await sets(a);
    if (!here.includes('My Favourites')) throw new Error('signing out took the sets away');
    if (await a.evaluate(() => !!window.CLOUD.session())) throw new Error('still signed in');
  });

  await step('the app works with the account switched off', async () => {
    const p = await b.newPage({ viewport: { width: 1280, height: 860 } });
    p.on('pageerror', e => errs.push('PAGEERROR (no account): ' + e.message.split('\n')[0]));
    await p.goto(INDEX);
    await p.waitForTimeout(900);
    // a build shipped without an account server: the settings file is simply empty
    await p.evaluate(() => { window.GC_CLOUD = { url: '', anonKey: '' }; });
    if (await p.evaluate(() => window.CLOUD.configured())) throw new Error('it thinks it has a server');
    await p.evaluate(() => window.__gc.menuAct('File', 'Sign In…'));
    await p.waitForTimeout(500);
    const said = await p.textContent('#dlg .bd');
    if (!/kept on this machine/.test(said)) throw new Error('it did not say what happens without a server');
    await p.close();
  });

  console.log('\nERRORS:', errs.length);
  errs.forEach(e => console.log('  ' + e));
  await b.close();
  process.exit(errs.length ? 1 : 0);
})();
