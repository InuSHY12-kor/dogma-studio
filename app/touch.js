/* Dogma Studio — using the map and the sequence with a finger.
 *
 * The application is written against mouse events. Rather than rewrite those handlers,
 * this layer reads touch gestures and hands the application the mouse events it already
 * understands, so one set of code serves both builds:
 *
 *   tap                     click: select a feature, an enzyme, a position
 *   long press on an item   its context menu
 *   long press on the page  then drag: select a stretch of DNA
 *   drag                    scroll, as a phone should
 *   pinch                   zoom the map in and out
 *
 * It only runs where the pointer is coarse, so a desktop never sees it.
 */
(function () {
  'use strict';
  if (!window.matchMedia || !window.matchMedia('(pointer: coarse)').matches) return;

  const LONG_MS = 420;          // how long a press has to be to mean "long press"
  const SLOP = 12;              // how far a finger may wander and still be a tap
  let press = null;             // the gesture in progress
  let pinch = null;

  const send = (type, t, opts) => {
    const ev = new MouseEvent(type, Object.assign({
      bubbles: true, cancelable: true, view: window,
      clientX: t.clientX, clientY: t.clientY,
      button: 0, buttons: type === 'mouseup' ? 0 : 1
    }, opts || {}));
    (document.elementFromPoint(t.clientX, t.clientY) || document.body).dispatchEvent(ev);
  };
  const sendAt = (type, el, t, opts) => {
    const ev = new MouseEvent(type, Object.assign({
      bubbles: true, cancelable: true, view: window,
      clientX: t.clientX, clientY: t.clientY, button: 0, buttons: 1
    }, opts || {}));
    el.dispatchEvent(ev);
  };

  // what the finger came down on: an item has its own menu, open canvas can be dragged
  const itemUnder = el => el && el.closest && (el.closest('[data-k]') || el.closest('[data-pick]') || el.closest('[data-enz]'));
  const surfaceUnder = el => el && el.closest && el.closest('#mapwrap svg, #seqscroll, .imapbox svg');

  function clearPress() {
    if (press && press.timer) clearTimeout(press.timer);
    if (press && press.dragging) document.documentElement.classList.remove('touch-drag');
    press = null;
  }

  document.addEventListener('touchstart', e => {
    if (e.touches.length === 2) {
      const [a, b] = e.touches;
      const el = surfaceUnder(document.elementFromPoint(a.clientX, a.clientY));
      if (el && el.closest('#mapwrap')) {
        pinch = { d: Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY), done: 0 };
        clearPress();
      }
      return;
    }
    if (e.touches.length !== 1) return;
    const t = e.touches[0];
    const el = document.elementFromPoint(t.clientX, t.clientY);
    const surface = surfaceUnder(el);
    if (!surface) return;                                   // ordinary page: leave it alone
    press = {
      x: t.clientX, y: t.clientY, t0: Date.now(), surface,
      item: itemUnder(el), moved: false, dragging: false, long: false
    };
    press.timer = setTimeout(() => {
      if (!press || press.moved) return;
      press.long = true;
      if (press.item) {
        // an item: its context menu, the same one a right-click gives
        const ev = new MouseEvent('contextmenu', { bubbles: true, cancelable: true, view: window, clientX: press.x, clientY: press.y, button: 2 });
        press.item.dispatchEvent(ev);
        clearPress();
      } else {
        // open canvas: a press that is held starts a selection the finger then drags out.
        // The press goes to whatever is under the finger, because the sequence view reads
        // the row from the event's own target.
        press.dragging = true;
        document.documentElement.classList.add('touch-drag');
        send('mousedown', { clientX: press.x, clientY: press.y });
      }
    }, LONG_MS);
  }, { passive: true });

  document.addEventListener('touchmove', e => {
    if (pinch && e.touches.length === 2) {
      const [a, b] = e.touches;
      const d = Math.hypot(a.clientX - b.clientX, a.clientY - b.clientY);
      const r = d / pinch.d;
      // one step per half-step of the pinch, so the map does not race away
      if (r > 1.25 || r < 0.8) {
        const btn = document.querySelector(r > 1 ? '.map-zoom [data-z="in"]' : '.map-zoom [data-z="out"]');
        if (btn) btn.click();
        pinch.d = d;
      }
      e.preventDefault();
      return;
    }
    if (!press || e.touches.length !== 1) return;
    const t = e.touches[0];
    if (!press.moved && Math.hypot(t.clientX - press.x, t.clientY - press.y) > SLOP) press.moved = true;
    if (press.dragging) {
      // a selection in progress: the application wants mousemove on the window
      window.dispatchEvent(new MouseEvent('mousemove', { bubbles: true, cancelable: true, view: window, clientX: t.clientX, clientY: t.clientY, buttons: 1 }));
      e.preventDefault();
    } else if (press.moved && !press.long) {
      clearPress();                                        // an ordinary scroll
    }
  }, { passive: false });

  const end = e => {
    if (pinch && e.touches.length < 2) pinch = null;
    if (!press) return;
    const t = (e.changedTouches && e.changedTouches[0]) || { clientX: press.x, clientY: press.y };
    if (press.dragging) {
      window.dispatchEvent(new MouseEvent('mouseup', { bubbles: true, cancelable: true, view: window, clientX: t.clientX, clientY: t.clientY, buttons: 0 }));
      clearPress();
      return;
    }
    if (!press.moved && !press.long) {
      // a tap: the press and release the application reads as a click
      send('mousedown', t);
      send('mouseup', t, { buttons: 0 });
    }
    clearPress();
  };
  document.addEventListener('touchend', end, { passive: true });
  document.addEventListener('touchcancel', () => { clearPress(); pinch = null; }, { passive: true });

  /* A finger has no hover, so the enzyme card a tap would never show is opened by the
   * tap itself: the application's own hover handler is given a mousemove first. */
  document.addEventListener('touchstart', e => {
    if (e.touches.length !== 1) return;
    const t = e.touches[0];
    const el = document.elementFromPoint(t.clientX, t.clientY);
    const enz = el && el.closest && el.closest('[data-enz]');
    if (enz) sendAt('mousemove', enz, t);
  }, { passive: true });
})();
