/* Dogma Studio — the Android host.
 *
 * The application code talks to its surroundings through one object, window.gcHost:
 * on the desktop that object comes from the Electron preload, and here it comes from
 * Capacitor. Nothing else in app/ knows which one it is talking to, so the two builds
 * run the same code.
 *
 * What this file does NOT provide is deliberate. There are no openTool or openAssembly
 * methods, because a phone has no second window: leaving them out makes the application
 * fall back to the dialogs it already uses in a plain browser.
 */
(function () {
  'use strict';
  const C = window.Capacitor;
  if (!C || !C.isNativePlatform || !C.isNativePlatform()) return;   // desktop or browser: nothing to do
  if (window.gcHost) return;                                        // never override the Electron preload
  const P = C.Plugins || {};
  const EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
  const TOOL = 'DogmaStudio';
  const q = o => Object.entries(o).filter(([, v]) => v !== undefined && v !== null && v !== '')
    .map(([k, v]) => `${k}=${encodeURIComponent(v)}`).join('&');

  // ---------------------------------------------------------------- files
  /* Android has no "open a file" API of its own that shows the places people keep
   * files; the system picker the WebView gives an <input type="file"> does, and it
   * reaches Drive, Downloads and anything else with a document provider. */
  function pickFiles(accept, multiple) {
    return new Promise(resolve => {
      const inp = document.createElement('input');
      inp.type = 'file';
      inp.accept = accept || '';
      if (multiple) inp.multiple = true;
      inp.style.cssText = 'position:fixed;left:-9999px;top:0';
      document.body.appendChild(inp);
      let done = false;
      const finish = list => { if (done) return; done = true; inp.remove(); resolve(list); };
      inp.onchange = async () => {
        const out = [];
        for (const f of [...(inp.files || [])]) {
          out.push({ name: f.name, path: f.name, data: await asBase64(f) });
        }
        finish(out);
      };
      // the picker can be dismissed without firing change, so the page coming back
      // into focus with nothing chosen ends the wait
      window.addEventListener('focus', () => setTimeout(() => { if (!inp.files || !inp.files.length) finish([]); }, 1200), { once: true });
      inp.click();
    });
  }
  const asBase64 = file => new Promise((res, rej) => {
    const r = new FileReader();
    r.onerror = () => rej(new Error('Could not read ' + file.name));
    r.onload = () => res(String(r.result).split(',')[1] || '');
    r.readAsDataURL(file);
  });

  const DOCS = 'DOCUMENTS';
  const b64 = text => {
    const bytes = new TextEncoder().encode(text);
    let s = '';
    for (const b of bytes) s += String.fromCharCode(b);
    return btoa(s);
  };
  /* Saving writes into the app's Documents folder and then offers the file to whatever
   * the phone can send it to, which is how a file leaves an Android app. The name that
   * comes back is what the application reports to the user. */
  async function saveFile(opts) {
    const o = opts || {};
    const name = (o.path && String(o.path).split('/').pop()) || o.defaultName || 'DogmaStudio.txt';
    const data = o.encoding === 'base64' ? o.content : b64(o.content || '');
    if (!P.Filesystem) throw new Error('This build cannot write files.');
    await P.Filesystem.writeFile({ path: name, data, directory: DOCS, recursive: true });
    const { uri } = await P.Filesystem.getUri({ path: name, directory: DOCS });
    if (P.Share) {
      try { await P.Share.share({ title: name, url: uri, dialogTitle: 'Send ' + name }); } catch (e) { /* the user closed the sheet; the file is still saved */ }
    }
    return name;
  }
  async function readFile(path) {
    const r = await P.Filesystem.readFile({ path, directory: DOCS });
    return { name: String(path).split('/').pop(), path, data: r.data };
  }

  // ---------------------------------------------------------------- settings
  const KEY = 'dogma-studio-settings';
  async function loadSettings() {
    if (!P.Preferences) return null;
    const { value } = await P.Preferences.get({ key: KEY });
    try { return value ? JSON.parse(value) : null; } catch (e) { return null; }
  }
  async function saveSettings(data) {
    if (!P.Preferences) return false;
    await P.Preferences.set({ key: KEY, value: JSON.stringify(data) });
    return true;
  }

  // ---------------------------------------------------------------- NCBI
  /* The desktop asks the main process so it can keep to NCBI's three-requests-a-second
   * rule; here the requests go straight out of the WebView, one at a time, with the same
   * spacing and the same reading of the answers. */
  let last = 0;
  async function paced(url) {
    const wait = Math.max(0, 340 - (Date.now() - last));
    if (wait) await new Promise(r => setTimeout(r, wait));
    last = Date.now();
    const r = await fetch(url, { headers: { Accept: 'text/plain' } });
    if (!r.ok) throw new Error(`NCBI answered ${r.status}. Try again in a moment.`);
    return r.text();
  }
  async function ncbiSearch({ db, term, retmax }) {
    const s = await paced(`${EUTILS}esearch.fcgi?${q({ db: db || 'nuccore', term, retmax: retmax || 20, retmode: 'json', tool: TOOL })}`);
    let ids = [];
    try { ids = (JSON.parse(s).esearchresult || {}).idlist || []; } catch (e) { throw new Error('NCBI sent a search result that could not be read.'); }
    if (!ids.length) return { ids: [], items: [] };
    const sum = await paced(`${EUTILS}esummary.fcgi?${q({ db: db || 'nuccore', id: ids.join(','), retmode: 'json', tool: TOOL })}`);
    const items = [];
    try {
      const r = JSON.parse(sum).result || {};
      for (const id of r.uids || []) {
        const it = r[id] || {};
        items.push({ uid: id, accession: it.caption || it.accessionversion || id, title: it.title || '', slen: it.slen || 0, topology: it.topology || '', organism: it.organism || '', updated: it.updatedate || '' });
      }
    } catch (e) { /* fall back to bare ids */ }
    return { ids, items };
  }
  async function ncbiFetch({ db, id, rettype }) {
    const text = await paced(`${EUTILS}efetch.fcgi?${q({ db: db || 'nuccore', id, rettype: rettype || 'gbwithparts', retmode: 'text', tool: TOOL })}`);
    const t = text.trim();
    if (!t) throw new Error(`NCBI returned nothing for “${id}”. Check the accession number.`);
    if (/^(Error|<\?xml|<!DOCTYPE)/i.test(t) || /id list is empty|cannot get document summary|Failed to retrieve/i.test(t.slice(0, 400))) {
      throw new Error(`NCBI could not find “${id}”. Check the accession number or try a search instead.`);
    }
    return text;
  }

  // ---------------------------------------------------------------- printing
  /* A phone has no print dialog to hand the page to, so the printable page is written
   * out as a file and offered to whatever can open or print it. */
  async function printHtml(html) {
    const name = `DogmaStudio-${Date.now()}.html`;
    await saveFile({ defaultName: name, content: html });
    return true;
  }

  window.gcHost = {
    platform: 'android',
    version: (C.getPlatform && C.getPlatform()) || 'android',
    openFiles: () => pickFiles('.gb,.gbk,.genbank,.dna,.fa,.fasta,.seq,.txt,.ab1,.json,.csv', true),
    saveFile,
    readFile,
    loadSettings,
    saveSettings,
    ncbiSearch,
    ncbiFetch,
    printHtml,
    onOpenFile: () => { },
    onCloseRequest: () => { },
    confirmClose: () => { }
  };

  // the status bar is part of the window on Android, so the chrome is told to leave room
  document.documentElement.classList.add('on-android');
  if (P.StatusBar && P.StatusBar.setOverlaysWebView) {
    P.StatusBar.setOverlaysWebView({ overlay: false }).catch(() => { });
  }
})();
