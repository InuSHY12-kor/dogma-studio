/* Dogma Studio — what the optional account does with two machines.
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
let clock = Date.now(), nextId = 1;
function dbGet() { return rows; }
function dbUpsert(list) {
  for (const row of list) {
    clock = Math.max(clock, Date.now()) + 1000;
    const at = rows.findIndex(r => r.user_id === row.user_id && r.kind === row.kind && r.name === row.name);
    const id = at >= 0 ? rows[at].id : 'id-' + (nextId++);
    const rec = Object.assign({ id }, row, { updated_at: new Date(clock).toISOString() });
    if (at >= 0) rows[at] = rec; else rows.push(rec);
  }
  return true;
}
/* A row taken out of the account is really gone, the way the real table deletes it —
 * no tombstone, nothing left behind. */
function dbDelete(ids) {
  let n = 0;
  for (const id of ids) {
    const at = rows.findIndex(r => r.id === id);
    if (at >= 0) { rows.splice(at, 1); n++; }
  }
  return n;
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
    if (u.includes('/rest/v1/dogma_items')) {
      const m = o.method || 'GET';
      if (m === 'GET') return reply(await window.__dbGet());
      if (m === 'DELETE') {
        const ids = (u.match(/id=in\.\(([^)]*)\)/) || [, ''])[1].split(',').filter(Boolean);
        await window.__dbDelete(ids);
        return reply(null, 204);
      }
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
    await p.exposeFunction('__dbDelete', dbDelete);
    await p.addInitScript(pageStub);
    await p.goto(INDEX);
    await p.waitForTimeout(900);
    // each machine keeps its own settings, as two computers would
    await p.evaluate(n => { const S = window.__gc.SET; S.enzymeSets = {}; S.cloudPick = {}; S.cloud = { items: {}, tombs: {}, lastSync: 0, email: '' }; }, name);
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
  // a primer or a feature is only shared where it is ticked, as the dialog does it
  const tick = (p, kind, on) => p.evaluate(async ([k, yes]) => {
    const S = window.__gc.SET;
    S.cloudPick = S.cloudPick || {};
    S.cloudPick[k] = yes ? window.CLOUD.choices(S, k).map(c => c.name) : [];
    await window.__gc.saveSettings();
  }, [kind, on]);
  const untick = (p, kind, has) => p.evaluate(async ([k, h]) => {
    const S = window.__gc.SET;
    S.cloudPick[k] = (S.cloudPick[k] || []).filter(n => !n.includes(h));
    await window.__gc.saveSettings();
  }, [kind, has]);
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

  await step('primers and my own features travel with their sequences', async () => {
    await a.evaluate(async () => {
      const S = window.__gc.SET;
      S.primerDb = S.primerDb || { primers: [], labels: [] };
      S.primerDb.primers.push({ id: 1, name: 'T7 fwd', seq: 'TAATACGACTCACTATAGGG', notes: 'stock tube 4', labels: ['stock'] });
      S.primerDb.labels = ['stock'];
      S.myFeatures = S.myFeatures || [];
      S.myFeatures.push({ name: 'my tag', type: 'misc_feature', kind: 'dna', seq: 'GGATCCGAAGGTGGA', molecule: 'DNA', color: '#ff8800', translated: false, matchProtein: false });
      S.myFeatures.push({ name: 'my tag', type: 'CDS', kind: 'protein', seq: 'MGSSHHHH', molecule: 'DNA', color: '#3366cc', translated: true, matchProtein: true });
      await window.__gc.saveSettings();
    });
    await tick(a, 'primer', true);
    await tick(a, 'feature', true);
    await syncOf(a);
    await syncOf(c);
    const got = await c.evaluate(() => {
      const S = window.__gc.SET;
      return {
        primers: (S.primerDb.primers || []).map(p => [p.name, p.seq, p.notes]),
        labels: (S.primerDb.labels || []).map(l => typeof l === 'string' ? l : l.name),
        feats: (S.myFeatures || []).map(f => [f.name, f.seq, f.kind, f.color, !!f.matchProtein])
      };
    });
    const t7 = got.primers.find(p => p[0] === 'T7 fwd');
    if (!t7) throw new Error('the primer did not arrive: ' + JSON.stringify(got.primers));
    if (t7[1] !== 'TAATACGACTCACTATAGGG') throw new Error('the primer arrived without its sequence: ' + JSON.stringify(t7));
    if (t7[2] !== 'stock tube 4') throw new Error('the note was dropped: ' + JSON.stringify(t7));
    if (!got.labels.includes('stock')) throw new Error('labels: ' + got.labels.join(','));
    // two features share a name and differ only in what they match: both have to arrive
    const dna = got.feats.find(f => f[1] === 'GGATCCGAAGGTGGA');
    const pro = got.feats.find(f => f[1] === 'MGSSHHHH');
    if (!dna || !pro) throw new Error('features: ' + JSON.stringify(got.feats));
    if (dna[2] !== 'dna' || dna[3] !== '#ff8800') throw new Error('the DNA feature lost something: ' + JSON.stringify(dna));
    if (pro[2] !== 'protein' || pro[4] !== true) throw new Error('the protein feature lost something: ' + JSON.stringify(pro));
  });

  await step('a saved primer set travels, primers and all', async () => {
    await a.evaluate(async () => {
      const S = window.__gc.SET;
      S.primerSets = S.primerSets || {};
      S.primerSets['pET screening'] = [
        { name: 'T7 fwd', seq: 'TAATACGACTCACTATAGGG', notes: '' },
        { name: 'T7 term', seq: 'GCTAGTTATTGCTCAGCGG', notes: 'reverse' }
      ];
      await window.__gc.saveSettings();
    });
    await syncOf(a);
    await syncOf(c);
    const set = await c.evaluate(() => (window.__gc.SET.primerSets || {})['pET screening']);
    if (!set) throw new Error('the set did not arrive: ' + JSON.stringify(await c.evaluate(() => Object.keys(window.__gc.SET.primerSets || {}))));
    if (set.length !== 2) throw new Error('the set arrived with ' + set.length + ' primers');
    if (set[0].seq !== 'TAATACGACTCACTATAGGG' || set[1].seq !== 'GCTAGTTATTGCTCAGCGG') throw new Error('a primer in the set lost its sequence: ' + JSON.stringify(set));
    if (set[1].notes !== 'reverse') throw new Error('a note in the set was dropped: ' + JSON.stringify(set));
  });

  await step('a feature type and its colour travel with the features that use it', async () => {
    await a.evaluate(async () => {
      const S = window.__gc.SET;
      S.featureTypes = S.featureTypes || [];
      S.featureTypes.push({ name: 'my linker', color: '#aa44cc', scope: 'dna' });
      S.typeColors = S.typeColors || {};
      S.typeColors['CDS'] = '#118822';
      await window.__gc.saveSettings();
    });
    await syncOf(a);
    await syncOf(c);
    const got = await c.evaluate(() => ({
      custom: (window.__gc.SET.featureTypes || []).map(t => [t.name, t.color]),
      std: window.__gc.SET.typeColors || {}
    }));
    const mine = got.custom.find(t => t[0] === 'my linker');
    if (!mine) throw new Error('the custom type did not arrive: ' + JSON.stringify(got.custom));
    if (mine[1] !== '#aa44cc') throw new Error('the custom type lost its colour: ' + JSON.stringify(mine));
    if (got.std.CDS !== '#118822') throw new Error('the recoloured standard type did not arrive: ' + JSON.stringify(got.std));
  });

  await step('a deleted primer set does not come back', async () => {
    await a.evaluate(async () => { delete window.__gc.SET.primerSets['pET screening']; await window.__gc.saveSettings(); });
    await syncOf(a);
    await syncOf(c);
    if (await c.evaluate(() => !!(window.__gc.SET.primerSets || {})['pET screening'])) throw new Error('it came back on the other machine');
    await syncOf(a);
    if (await a.evaluate(() => !!(window.__gc.SET.primerSets || {})['pET screening'])) throw new Error('it came back here');
  });

  await step('a primer that is not ticked never leaves the machine', async () => {
    await a.evaluate(async () => {
      const S = window.__gc.SET;
      S.primerDb.primers.push({ id: 9, name: 'Private fwd', seq: 'CCCGGGAAATTTCCCGGG', notes: 'not for the account' });
      await window.__gc.saveSettings();
    });
    await syncOf(a);                       // ticked list still holds only what was ticked
    await syncOf(c);
    const there = await c.evaluate(() => (window.__gc.SET.primerDb.primers || []).map(p => p.name));
    if (there.includes('Private fwd')) throw new Error('it went up without being ticked: ' + there.join(','));
    const here = await a.evaluate(() => (window.__gc.SET.primerDb.primers || []).map(p => p.name));
    if (!here.includes('Private fwd')) throw new Error('it was lost on the machine that made it');
  });

  await step('ticking it later sends it', async () => {
    await tick(a, 'primer', true);
    await syncOf(a);
    await syncOf(c);
    const there = await c.evaluate(() => (window.__gc.SET.primerDb.primers || []).map(p => [p.name, p.seq]));
    const got = there.find(p => p[0] === 'Private fwd');
    if (!got) throw new Error('ticking did not send it: ' + JSON.stringify(there));
    if (got[1] !== 'CCCGGGAAATTTCCCGGG') throw new Error('it arrived without its sequence');
  });

  await step('unticking takes it out of the account and off the other machine\u2019s list', async () => {
    await untick(a, 'primer', 'Private fwd');
    const st = await syncOf(a);
    if (!st.unshared) throw new Error('nothing was taken out of the account: ' + JSON.stringify(st));
    // it is gone from the table, not tombstoned
    const rowsNow = rows.filter(r => !r.deleted && r.kind === 'primer' && r.name.startsWith('Private fwd'));
    if (rowsNow.length) throw new Error('the row is still in the table: ' + JSON.stringify(rowsNow));
    // and the machine that made it still has it
    const here = await a.evaluate(() => (window.__gc.SET.primerDb.primers || []).map(p => p.name));
    if (!here.includes('Private fwd')) throw new Error('unticking deleted the primer here');
  });

  await step('the other machine keeps its copy but stops sharing it', async () => {
    const st = await syncOf(c);
    if (!st.unsharedThere) throw new Error('it did not notice: ' + JSON.stringify(st));
    const there = await c.evaluate(() => (window.__gc.SET.primerDb.primers || []).map(p => p.name));
    if (!there.includes('Private fwd')) throw new Error('unticking elsewhere deleted it here');
    const ticked = await c.evaluate(() => (window.__gc.SET.cloudPick.primer || []).some(n => n.startsWith('Private fwd')));
    if (ticked) throw new Error('it is still ticked on the other machine');
  });

  await step('and it does not creep back on the next sync', async () => {
    await syncOf(c); await syncOf(a); await syncOf(c);
    const back = rows.filter(r => !r.deleted && r.kind === 'primer' && r.name.startsWith('Private fwd'));
    if (back.length) throw new Error('it came back into the account');
    for (const [who, p] of [['A', a], ['B', c]]) {
      const has = await p.evaluate(() => (window.__gc.SET.primerDb.primers || []).map(x => x.name).includes('Private fwd'));
      if (!has) throw new Error('machine ' + who + ' lost its copy');
    }
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
