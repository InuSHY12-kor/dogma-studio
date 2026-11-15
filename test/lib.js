/* Shared pieces of the two test passes. */
'use strict';
const path = require('path');

const INDEX = 'file://' + path.join(__dirname, '..', 'app', 'index.html');

/* Each check prints one line and adds to the error list rather than stopping, so one
 * broken thing never hides the state of everything after it. A check also fails when
 * the page logged an error of its own while it ran. */
function runner(errs) {
  return async function step(name, fn) {
    const before = errs.length;
    try { await fn(); } catch (e) { errs.push('STEP ' + name + ': ' + e.message); }
    console.log((errs.length > before ? 'FAIL ' : 'ok   ') + name);
  };
}

/* A finger, for the phone pass. The gestures are the ones app/touch.js reads:
 * a tap, a long press, a long press that then drags, and a pinch. */
const touch = {
  async tap(p, sel) { await gesture(p, sel, { hold: 0, dx: 0, dy: 0 }); },
  async press(p, sel, ms) { await gesture(p, sel, { hold: ms === undefined ? 600 : ms, dx: 0, dy: 0 }); },
  async pressDrag(p, sel, dx, dy) { await gesture(p, sel, { hold: 600, dx, dy: dy || 0 }); },
  async swipe(p, sel, dx, dy) { await gesture(p, sel, { hold: 0, dx, dy: dy || 0, steps: 4 }); },
  async pinch(p, sel, out) {
    await p.evaluate(([s, o]) => {
      const el = document.querySelector(s);
      const r = el.getBoundingClientRect();
      const cx = r.left + r.width / 2, cy = r.top + r.height / 2;
      const mk = (tp, pts) => {
        const ts = pts.map((q, i) => new Touch({ identifier: i, target: el, clientX: q[0], clientY: q[1] }));
        el.dispatchEvent(new TouchEvent(tp, { bubbles: true, cancelable: true, touches: tp === 'touchend' ? [] : ts, changedTouches: ts, targetTouches: ts }));
      };
      const near = [[cx - 50, cy], [cx + 50, cy]], far = [[cx - 120, cy], [cx + 120, cy]];
      mk('touchstart', o ? near : far);
      mk('touchmove', o ? far : near);
      mk('touchend', o ? far : near);
    }, [sel, !!out]);
  }
};

async function gesture(p, sel, o) {
  await p.evaluate(async ([s, hold, dx, dy, steps]) => {
    const el = document.querySelector(s);
    if (!el) throw new Error('no element for ' + s);
    const r = el.getBoundingClientRect();
    const x = r.left + Math.min(40, r.width / 2), y = r.top + r.height / 2;
    const mk = (tp, cx, cy) => {
      const t = new Touch({ identifier: 1, target: el, clientX: cx, clientY: cy });
      el.dispatchEvent(new TouchEvent(tp, {
        bubbles: true, cancelable: true,
        touches: tp === 'touchend' ? [] : [t], changedTouches: [t], targetTouches: tp === 'touchend' ? [] : [t]
      }));
    };
    mk('touchstart', x, y);
    if (hold) await new Promise(r2 => setTimeout(r2, hold));
    const n = steps || (dx || dy ? 6 : 0);
    for (let k = 1; k <= n; k++) { mk('touchmove', x + dx * k / n, y + dy * k / n); await new Promise(r2 => setTimeout(r2, 30)); }
    mk('touchend', x + dx, y + dy);
  }, [sel, o.hold, o.dx, o.dy, o.steps || 0]);
}

module.exports = { INDEX, runner, touch };
