// GeneCanvas — Electron main process
const { app, BrowserWindow, ipcMain, dialog, Menu, shell, net } = require('electron');
const path = require('path');
const fs = require('fs');
const os = require('os');

let win = null;
let allowClose = false;
const pendingFiles = [];

const OPENABLE = /\.(gb|gbk|genbank|gbff|ape|fa|fasta|fna|seq|txt|dna|ab1)$/i;

function fileArgs(argv) {
  return argv.slice(1).filter(a => !a.startsWith('-') && OPENABLE.test(a) && fs.existsSync(a));
}

function readForRenderer(p) {
  return { path: p, name: path.basename(p), data: fs.readFileSync(p).toString('base64') };
}

function sendFile(p) {
  if (!win || win.webContents.isLoading()) { pendingFiles.push(p); return; }
  try { win.webContents.send('open-file', readForRenderer(p)); } catch (e) { dialog.showErrorBox('GeneCanvas', `Could not open ${p}\n${e.message}`); }
}

function createWindow() {
  win = new BrowserWindow({
    width: 1440,
    height: 920,
    minWidth: 900,
    minHeight: 600,
    title: 'GeneCanvas',
    backgroundColor: '#ffffff',
    icon: path.join(__dirname, 'build', 'icon.png'),
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      contextIsolation: true,
      nodeIntegration: false,
      spellcheck: false
    }
  });
  Menu.setApplicationMenu(null); // the app draws its own menu bar
  win.loadFile(path.join(__dirname, 'app', 'index.html'));

  win.webContents.on('did-finish-load', () => {
    const files = pendingFiles.splice(0);
    for (const f of files) sendFile(f);
  });
  // external links open in the default browser
  win.webContents.setWindowOpenHandler(({ url }) => { if (/^https?:/.test(url)) shell.openExternal(url); return { action: 'deny' }; });
  win.webContents.on('will-navigate', e => e.preventDefault());

  win.on('close', e => {
    if (allowClose) return;
    e.preventDefault();
    win.webContents.send('close-request');
  });
  win.on('closed', () => { win = null; });
}

// one window; double-clicked files go to the running instance
if (!app.requestSingleInstanceLock()) {
  app.quit();
} else {
  app.on('second-instance', (_e, argv) => {
    if (win) { if (win.isMinimized()) win.restore(); win.focus(); }
    for (const f of fileArgs(argv)) sendFile(f);
  });
  app.whenReady().then(() => {
    pendingFiles.push(...fileArgs(process.argv));
    createWindow();
  });
  app.on('window-all-closed', () => app.quit());
}

ipcMain.handle('open-files', async () => {
  const r = await dialog.showOpenDialog(win, {
    title: 'Open DNA File',
    properties: ['openFile', 'multiSelections'],
    filters: [
      { name: 'DNA files', extensions: ['gb', 'gbk', 'genbank', 'gbff', 'dna', 'fa', 'fasta', 'fna', 'seq', 'txt', 'ape'] },
      { name: 'Sequencing traces', extensions: ['ab1'] },
      { name: 'All files', extensions: ['*'] }
    ]
  });
  if (r.canceled) return [];
  return r.filePaths.map(readForRenderer);
});

ipcMain.handle('save-file', async (_e, { path: p, defaultName, content, filters }) => {
  let target = p;
  if (!target) {
    const r = await dialog.showSaveDialog(win, { title: 'Save', defaultPath: defaultName, filters: [...(filters || []), { name: 'All files', extensions: ['*'] }] });
    if (r.canceled || !r.filePath) return null;
    target = r.filePath;
  }
  fs.writeFileSync(target, content, 'utf8');
  return target;
});

ipcMain.on('confirm-close', () => { allowClose = true; if (win) win.close(); });
ipcMain.on('quit', () => { if (win) win.close(); });
ipcMain.on('version', e => { e.returnValue = app.getVersion(); });

// ---------------------------------------------------------------- NCBI E-utilities
// Downloads sequence records by accession. NCBI asks unauthenticated clients to stay
// under three requests a second, so every call goes through one paced queue.
const EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
const TOOL = 'GeneCanvas';
let lastCall = 0;
let chain = Promise.resolve();

function httpGet(url) {
  return new Promise((resolve, reject) => {
    const req = net.request({ url, method: 'GET' });
    const timer = setTimeout(() => { req.abort(); reject(new Error('NCBI did not respond in time.')); }, 30000);
    req.on('response', res => {
      let body = '';
      res.on('data', c => { body += c.toString('utf8'); });
      res.on('end', () => {
        clearTimeout(timer);
        if (res.statusCode >= 400) reject(new Error(`NCBI returned HTTP ${res.statusCode}.`));
        else resolve(body);
      });
      res.on('error', e => { clearTimeout(timer); reject(e); });
    });
    req.on('error', e => {
      clearTimeout(timer);
      reject(new Error(`Could not reach NCBI: ${e.message}. Check the internet connection.`));
    });
    req.end();
  });
}
function paced(url) {
  chain = chain.then(async () => {
    const wait = Math.max(0, 350 - (Date.now() - lastCall));
    if (wait) await new Promise(r => setTimeout(r, wait));
    lastCall = Date.now();
    return httpGet(url);
  }, async () => {
    lastCall = Date.now();
    return httpGet(url);
  });
  return chain;
}
const q = o => Object.entries(o).map(([k, v]) => `${k}=${encodeURIComponent(v)}`).join('&');

ipcMain.handle('ncbi-search', async (_e, { db, term, retmax }) => {
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
});

ipcMain.handle('ncbi-fetch', async (_e, { db, id, rettype }) => {
  const text = await paced(`${EUTILS}efetch.fcgi?${q({ db: db || 'nuccore', id, rettype: rettype || 'gbwithparts', retmode: 'text', tool: TOOL })}`);
  const t = text.trim();
  if (!t) throw new Error(`NCBI returned nothing for “${id}”. Check the accession number.`);
  if (/^(Error|<\?xml|<!DOCTYPE)/i.test(t) || /id list is empty|cannot get document summary|Failed to retrieve/i.test(t.slice(0, 400))) {
    throw new Error(`NCBI could not find “${id}”. Check the accession number or try a search instead.`);
  }
  return text;
});

// ---------------------------------------------------------------- printing
function withPrintWindow(html, fn) {
  const file = path.join(os.tmpdir(), `genecanvas-print-${Date.now()}.html`);
  fs.writeFileSync(file, html, 'utf8');
  const w = new BrowserWindow({ show: false, webPreferences: { sandbox: true, javascript: false } });
  const cleanup = () => { try { fs.unlinkSync(file); } catch (e) { } if (!w.isDestroyed()) w.destroy(); };
  return w.loadFile(file).then(() => fn(w)).then(
    r => { cleanup(); return r; },
    e => { cleanup(); throw e; }
  );
}

ipcMain.handle('print-html', async (_e, { html }) =>
  withPrintWindow(html, w => new Promise(resolve => {
    w.webContents.print({ silent: false, printBackground: true }, (ok, reason) => resolve({ ok, reason }));
  })));

ipcMain.handle('pdf-html', async (_e, { html, defaultName, landscape }) => {
  const r = await dialog.showSaveDialog(win, { title: 'Export PDF', defaultPath: defaultName, filters: [{ name: 'PDF document', extensions: ['pdf'] }] });
  if (r.canceled || !r.filePath) return null;
  const data = await withPrintWindow(html, w => w.webContents.printToPDF({ printBackground: true, landscape: !!landscape, margins: { marginType: 'default' } }));
  fs.writeFileSync(r.filePath, data);
  return r.filePath;
});
