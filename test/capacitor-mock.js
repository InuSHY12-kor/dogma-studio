/* A stand-in for the Capacitor bridge, for the phone pass.
 *
 * In a plain browser window.Capacitor does not exist, so app/mobileboot.js steps aside
 * and the application falls back to its browser paths — which means the Android host
 * would never be exercised by a test. This script is injected before the page's own
 * scripts and puts a bridge there: it answers like the real one and writes down every
 * call, so a check can say "saving really did write that file and open the share sheet".
 *
 * It is a test double, not a shim that ships: nothing in app/ loads it.
 */
module.exports = function capacitorMock() {
  const store = { files: {}, prefs: {}, calls: [] };
  const note = (name, args) => { store.calls.push({ name, args, at: Date.now() }); };
  const ok = v => Promise.resolve(v);

  window.__native = store;
  window.Capacitor = {
    isNativePlatform: () => true,
    getPlatform: () => 'android',
    Plugins: {
      Filesystem: {
        writeFile: o => { note('Filesystem.writeFile', o); store.files[o.path] = o.data; return ok({ uri: 'file:///Documents/' + o.path }); },
        readFile: o => {
          note('Filesystem.readFile', o);
          if (!(o.path in store.files)) return Promise.reject(new Error('no such file: ' + o.path));
          return ok({ data: store.files[o.path] });
        },
        getUri: o => ok({ uri: 'file:///Documents/' + o.path }),
        mkdir: o => ok({})
      },
      Preferences: {
        get: o => { note('Preferences.get', o); return ok({ value: store.prefs[o.key] === undefined ? null : store.prefs[o.key] }); },
        set: o => { note('Preferences.set', o); store.prefs[o.key] = o.value; return ok({}); },
        remove: o => { delete store.prefs[o.key]; return ok({}); }
      },
      Share: {
        share: o => { note('Share.share', o); return ok({ activityType: 'test' }); }
      },
      StatusBar: {
        setOverlaysWebView: o => { note('StatusBar.setOverlaysWebView', o); return ok({}); },
        setStyle: o => ok({}),
        setBackgroundColor: o => ok({})
      },
      App: {
        addListener: (ev, cb) => { store['on' + ev] = cb; return ok({ remove: () => { } }); },
        exitApp: () => { note('App.exitApp', {}); return ok({}); }
      }
    }
  };
};
