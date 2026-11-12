// GeneCanvas — Electron main process
const { app, BrowserWindow, ipcMain, dialog, Menu, shell } = require('electron');
const path = require('path');
const fs = require('fs');

let win = null;
let allowClose = false;
const pendingFiles = [];

const OPENABLE = /\.(gb|gbk|genbank|gbff|ape|fa|fasta|fna|seq|txt|dna)$/i;

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
