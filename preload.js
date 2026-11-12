// Exposes a small, explicit API to the page (window.gcHost)
const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('gcHost', {
  version: ipcRenderer.sendSync('version'),
  openFiles: () => ipcRenderer.invoke('open-files'),
  saveFile: opts => ipcRenderer.invoke('save-file', opts),
  onOpenFile: cb => ipcRenderer.on('open-file', (_e, f) => cb(f)),
  onCloseRequest: cb => ipcRenderer.on('close-request', () => cb()),
  confirmClose: () => ipcRenderer.send('confirm-close'),
  quit: () => ipcRenderer.send('quit')
});
