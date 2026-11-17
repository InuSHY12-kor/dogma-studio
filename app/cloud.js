/* Dogma Studio — the optional account.
 *
 * Signing in is never required: everything works exactly as before without it, kept on
 * the machine it was made on. With an account, three libraries follow the user instead
 * — the enzyme sets, the primer database and the common features they added themselves.
 * Sequences are not uploaded.
 *
 * It talks to Supabase over plain HTTP, so the application needs no build step: GoTrue
 * for the account, PostgREST for the rows. Which rows a request may touch is decided by
 * the database, not here: every row carries its owner and the policies only ever let a
 * user see their own.
 */
(function (root) {
  'use strict';
  const CFG = () => root.GC_CLOUD || {};
  const SESSION_KEY = 'dogma-studio-cloud-session';
  const TABLE = 'dogma_items';

  /* What is kept in the cloud, and how each library is taken apart into items. One item
   * is one thing a person made: a set, a primer, a label, a feature. They are matched up
   * between machines by kind and name, so the same set edited in two places is one set. */
  const KINDS = {
    enzyme_set: {
      list: SET => Object.entries(SET.enzymeSets || {}).map(([name, enzymes]) => ({ name, payload: { enzymes } })),
      apply: (SET, name, payload) => { (SET.enzymeSets = SET.enzymeSets || {})[name] = (payload && payload.enzymes) || []; },
      drop: (SET, name) => { if (SET.enzymeSets) delete SET.enzymeSets[name]; }
    },
    primer: {
      list: SET => ((SET.primerDb || {}).primers || []).map(p => ({ name: p.name + '\u0001' + (p.seq || ''), payload: { ...p, id: undefined } })),
      apply: (SET, name, payload) => {
        const db = (SET.primerDb = SET.primerDb || { primers: [], labels: [] });
        const key = (payload.name || '') + '\u0001' + (payload.seq || '');
        const at = db.primers.findIndex(p => (p.name + '\u0001' + (p.seq || '')) === key);
        const id = at >= 0 ? db.primers[at].id : Math.max(0, ...db.primers.map(x => x.id || 0)) + 1;
        const rec = { ...payload, id };
        if (at >= 0) db.primers[at] = rec; else db.primers.push(rec);
      },
      drop: (SET, name) => {
        const db = SET.primerDb; if (!db) return;
        db.primers = (db.primers || []).filter(p => (p.name + '\u0001' + (p.seq || '')) !== name);
      }
    },
    primer_label: {
      list: SET => (((SET.primerDb || {}).labels) || []).map(l => ({ name: typeof l === 'string' ? l : l.name, payload: typeof l === 'string' ? { name: l } : l })),
      apply: (SET, name, payload) => {
        const db = (SET.primerDb = SET.primerDb || { primers: [], labels: [] });
        db.labels = db.labels || [];
        const same = l => (typeof l === 'string' ? l : l.name) === name;
        if (!db.labels.some(same)) db.labels.push(payload && payload.color !== undefined ? payload : name);
      },
      drop: (SET, name) => {
        const db = SET.primerDb; if (!db || !db.labels) return;
        db.labels = db.labels.filter(l => (typeof l === 'string' ? l : l.name) !== name);
      }
    },
    feature: {
      /* A feature is what it matches, so its own sequence is part of what tells two of
       * them apart: the same name over a different stretch of DNA is a different feature
       * and both are kept. The payload is the whole record — sequence, type, colour,
       * whether it is matched as protein — so a machine that takes it can detect with it
       * exactly as the machine that made it. */
      list: SET => (SET.myFeatures || []).map(f => ({ name: featKey(f), payload: f })),
      apply: (SET, name, payload) => {
        SET.myFeatures = SET.myFeatures || [];
        const at = SET.myFeatures.findIndex(f => featKey(f) === name);
        if (at >= 0) SET.myFeatures[at] = payload; else SET.myFeatures.push(payload);
      },
      drop: (SET, name) => { SET.myFeatures = (SET.myFeatures || []).filter(f => featKey(f) !== name); }
    },
    /* A saved primer set is a named list of primers, each with its own sequence. It is
     * what Primers › Saved Sets keeps, and it travels whole. */
    primer_set: {
      list: SET => Object.entries(SET.primerSets || {}).map(([name, primers]) => ({ name, payload: { primers } })),
      apply: (SET, name, payload) => { (SET.primerSets = SET.primerSets || {})[name] = (payload && payload.primers) || []; },
      drop: (SET, name) => { if (SET.primerSets) delete SET.primerSets[name]; }
    },
    /* Feature types the user made, and the colours they gave standard ones. Without
     * these a feature that followed the account would be drawn in the default grey on
     * the other machine, under a type it does not know. The two are told apart by the
     * key's first part, so a custom type never collides with a recoloured standard one. */
    feature_type: {
      list: SET => [
        ...(SET.featureTypes || []).map(t => ({ name: 'custom\u0001' + t.name, payload: t })),
        ...Object.entries(SET.typeColors || {}).map(([name, color]) => ({ name: 'colour\u0001' + name, payload: { name, color } }))
      ],
      apply: (SET, key, payload) => {
        const [what, name] = key.split('\u0001');
        if (what === 'custom') {
          SET.featureTypes = SET.featureTypes || [];
          const at = SET.featureTypes.findIndex(t => t.name === name);
          if (at >= 0) SET.featureTypes[at] = payload; else SET.featureTypes.push(payload);
        } else {
          (SET.typeColors = SET.typeColors || {})[name] = payload && payload.color;
        }
      },
      drop: (SET, key) => {
        const [what, name] = key.split('\u0001');
        if (what === 'custom') SET.featureTypes = (SET.featureTypes || []).filter(t => t.name !== name);
        else if (SET.typeColors) delete SET.typeColors[name];
      }
    }
  };

  // what makes one of the user's own features that one and not another
  const featKey = f => (f.name || '') + '\u0001' + (f.seq || '');

  /* Two kinds are chosen one at a time. A primer and one of the user's own features go
   * to the account only where the user has ticked them, because a primer database is a
   * personal thing and not all of it is meant to leave the machine. Everything else —
   * the named sets, the primer labels, the feature types — follows the account whole.
   *
   * Unticking is not deleting: the row is taken out of the account for good, on every
   * machine, while the primer or the feature stays where it is on this one. */
  const CHOSEN = { primer: true, feature: true };
  const pickOf = SET => (SET.cloudPick = SET.cloudPick || {});
  const isShared = (SET, kind, name) => !CHOSEN[kind] || (pickOf(SET)[kind] || []).indexOf(name) >= 0;
  function setShared(SET, kind, name, on) {
    if (!CHOSEN[kind]) return;
    const p = pickOf(SET);
    const s = new Set(p[kind] || []);
    if (on) s.add(name); else s.delete(name);
    p[kind] = [...s];
  }
  /* Everything of one kind that is on this machine, each said to be ticked or not, so a
   * dialog can show the list without knowing how any of it is stored. */
  const choices = (SET, kind) => (KINDS[kind] ? KINDS[kind].list(SET) : [])
    .map(it => ({ name: it.name, payload: it.payload, shared: isShared(SET, kind, it.name) }));

  // ---------------------------------------------------------------- small helpers
  const hash = o => {
    const s = JSON.stringify(o);
    let h = 5381;
    for (let i = 0; i < s.length; i++) h = ((h * 33) ^ s.charCodeAt(i)) >>> 0;
    return h.toString(36) + ':' + s.length;
  };
  const stateOf = SET => (SET.cloud = SET.cloud || { items: {}, tombs: {}, lastSync: 0, email: '' });
  const keyOf = (kind, name) => kind + '\u0000' + name;
  const localItems = SET => {
    const out = new Map();
    for (const kind of Object.keys(KINDS)) for (const it of KINDS[kind].list(SET)) out.set(keyOf(kind, it.name), { kind, name: it.name, payload: it.payload });
    return out;
  };

  /* What this machine knows about each item: the payload as it stood at the last sync,
   * and the server's own stamp for that version. An item is "changed here" when its
   * payload no longer hashes to what was last synced, and "changed there" when the
   * server's stamp is not the one that was last seen. Two clocks are never compared —
   * a phone and a laptop rarely agree on the time, and the server sets its own stamps. */
  function stamp(SET) {
    const st = stateOf(SET);
    const seen = localItems(SET);
    for (const [k, it] of seen) {
      const h = hash(it.payload);
      const was = st.items[k];
      if (!was) st.items[k] = { h: null, v: null, now: h };        // new here, never synced
      else was.now = h;
      delete st.tombs[k];
    }
    for (const k of Object.keys(st.items)) {
      if (seen.has(k)) continue;
      const was = st.items[k];
      // a tombstone remembers the server stamp it was last seen at, so a deletion here
      // can be told apart from an edit made somewhere else
      if (was && was.v) st.tombs[k] = { v: was.v };
      delete st.items[k];
    }
    return st;
  }

  // ---------------------------------------------------------------- the account
  const configured = () => { const c = CFG(); return !!(c.url && c.anonKey); };
  let session = null;
  // the key is named after the program; one saved under the old name still signs in
  try {
    session = JSON.parse(localStorage.getItem(SESSION_KEY) || localStorage.getItem('genecanvas-cloud-session') || 'null');
  } catch (e) { session = null; }
  const saveSession = s => {
    session = s;
    try { s ? localStorage.setItem(SESSION_KEY, JSON.stringify(s)) : localStorage.removeItem(SESSION_KEY); } catch (e) { /* a private window */ }
  };

  async function api(path, opts) {
    if (!configured()) throw new Error('This build has no account server configured.');
    const o = opts || {};
    const headers = Object.assign({ apikey: CFG().anonKey, 'Content-Type': 'application/json' }, o.headers || {});
    if (o.auth !== false && session && session.access_token) headers.Authorization = 'Bearer ' + session.access_token;
    let r;
    try {
      r = await fetch(CFG().url + path, { method: o.method || 'GET', headers, body: o.body ? JSON.stringify(o.body) : undefined });
    } catch (e) {
      throw new Error('Could not reach the account server. Check the connection.');
    }
    if (r.status === 401 && o.retry !== false && session && session.refresh_token) {
      await refresh();
      return api(path, Object.assign({}, o, { retry: false }));
    }
    const text = await r.text();
    let data = null;
    try { data = text ? JSON.parse(text) : null; } catch (e) { data = text; }
    if (!r.ok) {
      const msg = (data && (data.error_description || data.msg || data.message || data.error)) || `${r.status} ${r.statusText}`;
      throw new Error(String(msg));
    }
    return data;
  }

  const useToken = d => {
    if (!d || !d.access_token) throw new Error('The server sent no session.');
    saveSession({
      access_token: d.access_token, refresh_token: d.refresh_token,
      email: (d.user && d.user.email) || (session && session.email) || '',
      user_id: (d.user && d.user.id) || (session && session.user_id) || ''
    });
    return session;
  };
  async function signIn(email, password) {
    return useToken(await api('/auth/v1/token?grant_type=password', { method: 'POST', auth: false, body: { email, password } }));
  }
  async function signUp(email, password) {
    const d = await api('/auth/v1/signup', { method: 'POST', auth: false, body: { email, password } });
    if (d && d.access_token) return useToken(d);
    // the project asks for the address to be confirmed first
    return { pending: true, email };
  }
  async function refresh() {
    const d = await api('/auth/v1/token?grant_type=refresh_token', { method: 'POST', auth: false, retry: false, body: { refresh_token: session.refresh_token } })
      .catch(() => null);
    if (!d) { saveSession(null); throw new Error('The session has expired. Sign in again.'); }
    return useToken(d);
  }
  function signOut() { saveSession(null); }

  // ---------------------------------------------------------------- syncing
  /* One pass in both directions. Each item is matched by kind and name; whichever side
   * changed it later wins, and a deletion is an item too, so removing a set on one
   * machine removes it on the other instead of coming back on the next sync. */
  async function sync(SET) {
    if (!session) throw new Error('Sign in first.');
    const st = stamp(SET);
    const rows = await api(`/rest/v1/${TABLE}?select=id,kind,name,payload,updated_at,deleted`);
    const server = new Map();
    for (const r of rows || []) server.set(keyOf(r.kind, r.name), r);

    const mine = localItems(SET);
    const push = [];
    const drop = [];                       // rows to take out of the account for good
    const stats = { sent: 0, taken: 0, removedHere: 0, removedThere: 0, kept: 0, unshared: 0, unsharedThere: 0 };
    const send = (kind, name, payload, deleted) => { push.push({ kind, name, payload: payload || {}, deleted: !!deleted }); };

    const keys = new Set([...mine.keys(), ...server.keys(), ...Object.keys(st.tombs)]);
    for (const k of keys) {
      const it = mine.get(k);
      const row = server.get(k);
      const was = st.items[k] || null;
      const [kind, name] = k.split('\u0000');
      if (!KINDS[kind]) continue;

      /* Unticked here: the row leaves the account altogether, which is a plain delete
       * and not a tombstone — the other machines are meant to keep their copy and
       * simply stop sharing it, not to lose it. */
      if (it && !isShared(SET, kind, name)) {
        if (row && !row.deleted) { drop.push(row.id); stats.unshared++; }
        continue;
      }
      const changedHere = it ? (!was || was.h !== hash(it.payload)) : false;
      const changedThere = row ? (!was || was.v !== row.updated_at) : false;

      if (it && !row) {
        /* A row this machine had seen is gone, with no tombstone in its place: it was
         * unticked somewhere else. The copy here stays, and stops being shared. */
        if (was && was.v) { setShared(SET, kind, name, false); stats.unsharedThere++; }
        else { send(kind, name, it.payload, false); stats.sent++; }
        continue;
      }
      if (!it && st.tombs[k]) {
        // gone from this machine since the last sync
        const seenAt = st.tombs[k].v;
        if (row && !row.deleted && row.updated_at === seenAt) { send(kind, name, row.payload, true); stats.removedThere++; }
        else if (row && !row.deleted) { KINDS[kind].apply(SET, name, row.payload); setShared(SET, kind, name, true); stats.taken++; }
        continue;
      }
      if (!it && row) {
        // it is in the account, so it is shared by definition once it lands here
        if (!row.deleted) { KINDS[kind].apply(SET, name, row.payload); setShared(SET, kind, name, true); stats.taken++; }
        continue;
      }
      if (row.deleted) {
        if (changedHere) { send(kind, name, it.payload, false); stats.sent++; }      // edited here after it went
        else { KINDS[kind].drop(SET, name); stats.removedHere++; }
        continue;
      }
      if (changedHere && !changedThere) { send(kind, name, it.payload, false); stats.sent++; continue; }
      if (!changedHere && changedThere) { KINDS[kind].apply(SET, name, row.payload); stats.taken++; continue; }
      if (changedHere && changedThere) {
        /* Both sides moved. For the named lists nothing is thrown away: the account's
         * copy becomes the one under this name, and what this machine had is kept beside
         * it so it can be looked at and merged by hand. The other kinds are each one
         * thing rather than a list — a primer, a feature, a type's colour — and two of
         * them under one name would only read as a duplicate, so there the account's
         * copy is taken. A feature carries its sequence in its key, so an edit that
         * changes what it matches is a new feature and both are kept anyway. */
        if (hash(it.payload) !== hash(row.payload)) {
          if (kind === 'enzyme_set' || kind === 'primer_set') {
            const alt = name + ' (this device)';
            KINDS[kind].apply(SET, alt, it.payload);
            send(kind, alt, it.payload, false);
            stats.kept++;
          }
          KINDS[kind].apply(SET, name, row.payload);
          stats.taken++;
        }
      }
    }

    if (drop.length) {
      await api(`/rest/v1/${TABLE}?id=in.(${drop.join(',')})`, { method: 'DELETE', headers: { Prefer: 'return=minimal' } });
    }
    if (push.length) {
      const body = push.map(p => Object.assign({ user_id: session.user_id }, p));
      await api(`/rest/v1/${TABLE}?on_conflict=user_id,kind,name`, {
        method: 'POST', headers: { Prefer: 'resolution=merge-duplicates,return=minimal' }, body
      });
    }

    /* The merged state is read back, so what this machine remembers about each item is
     * the server's own stamp rather than a guess. */
    const after = await api(`/rest/v1/${TABLE}?select=id,kind,name,payload,updated_at,deleted`);
    st.items = {};
    st.tombs = {};
    const now = localItems(SET);
    for (const r of after || []) {
      const k = keyOf(r.kind, r.name);
      if (r.deleted || !now.has(k)) continue;
      st.items[k] = { h: hash(now.get(k).payload), v: r.updated_at };
    }
    for (const [k, it] of now) if (!st.items[k]) st.items[k] = { h: hash(it.payload), v: null };
    st.lastSync = Date.now();
    st.email = session.email;
    return stats;
  }

  root.CLOUD = {
    configured, stamp, sync, signIn, signUp, signOut,
    session: () => session,
    kinds: Object.keys(KINDS),
    chosenKinds: Object.keys(CHOSEN),
    choices, isShared, setShared,
    countShared: SET => [...localItems(SET)].filter(([k, it]) => isShared(SET, it.kind, it.name)).length,
    countLocal: SET => localItems(SET).size
  };
})(typeof window !== 'undefined' ? window : globalThis);
