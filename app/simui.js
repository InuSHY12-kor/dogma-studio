// placeholder so the tool window loads before SIMUI and CLONEUI exist
window.SIMUI = window.SIMUI || { open: (el) => { el.innerHTML = '<div class="hint">Simulator coming next.</div>'; } };
window.CLONEUI = window.CLONEUI || { open: (el) => { el.innerHTML = '<div class="hint">Cloning window coming next.</div>'; } };
