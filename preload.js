// Exposes a small, explicit API to the page (window.gcHost)
const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('gcHost', {
  version: ipcRenderer.sendSync('version'),
  openFiles: () => ipcRenderer.invoke('open-files'),
  saveFile: opts => ipcRenderer.invoke('save-file', opts),
  readFile: path => ipcRenderer.invoke('read-file', { path }),
  onOpenFile: cb => ipcRenderer.on('open-file', (_e, f) => cb(f)),
  onCloseRequest: cb => ipcRenderer.on('close-request', () => cb()),
  confirmClose: () => ipcRenderer.send('confirm-close'),
  quit: () => ipcRenderer.send('quit'),
  ncbiSearch: opts => ipcRenderer.invoke('ncbi-search', opts),
  ncbiFetch: opts => ipcRenderer.invoke('ncbi-fetch', opts),
  loadSettings: () => ipcRenderer.invoke('load-settings'),
  saveSettings: data => ipcRenderer.invoke('save-settings', data),
  openAssembly: payload => ipcRenderer.invoke('assembly-open', payload),
  assemblyInit: () => ipcRenderer.invoke('assembly-init'),
  assemblySync: docs => ipcRenderer.send('assembly-sync', docs),
  assemblyCreate: src => ipcRenderer.send('assembly-create', src),
  onAssemblyDocs: cb => ipcRenderer.on('assembly-docs', (_e, d) => cb(d)),
  onAssemblyReload: cb => ipcRenderer.on('assembly-reload', (_e, p) => cb(p)),
  onAssemblyProduct: cb => ipcRenderer.on('assembly-product', (_e, s) => cb(s)),
  closeSelf: () => ipcRenderer.send('assembly-close'),
  printHtml: html => ipcRenderer.invoke('print-html', { html }),
  pdfHtml: opts => ipcRenderer.invoke('pdf-html', opts)
});
