// Exposes a small, explicit API to the page (window.gcHost)
const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('gcHost', {
  version: ipcRenderer.sendSync('version'),
  openFiles: () => ipcRenderer.invoke('open-files'),
  saveFile: opts => ipcRenderer.invoke('save-file', opts),
  onOpenFile: cb => ipcRenderer.on('open-file', (_e, f) => cb(f)),
  onCloseRequest: cb => ipcRenderer.on('close-request', () => cb()),
  confirmClose: () => ipcRenderer.send('confirm-close'),
  quit: () => ipcRenderer.send('quit'),
  ncbiSearch: opts => ipcRenderer.invoke('ncbi-search', opts),
  ncbiFetch: opts => ipcRenderer.invoke('ncbi-fetch', opts),
  loadSettings: () => ipcRenderer.invoke('load-settings'),
  saveSettings: data => ipcRenderer.invoke('save-settings', data),
  printHtml: html => ipcRenderer.invoke('print-html', { html }),
  pdfHtml: opts => ipcRenderer.invoke('pdf-html', opts)
});
