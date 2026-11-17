/* Dogma Studio — agarose gel simulation. Exposes window.GEL
 *
 * Migration is modelled, not measured. For each agarose percentage the model holds the
 * two ends of the useful resolving range: the size that barely leaves the well and the
 * size that has run off the end of the gel. Everything between moves with the logarithm
 * of its length, which is how DNA behaves over that range. Band sizes for the ladders
 * follow the suppliers' published values and are worth checking against the datasheet
 * in the box, since suppliers do revise them.
 */
(function (root) {
  'use strict';

  // percent agarose -> [largest resolved, smallest resolved] in bp
  const RANGE = [
    [0.5, 30000, 700], [0.7, 20000, 400], [0.8, 15000, 300], [1.0, 12000, 200],
    [1.2, 8000, 120], [1.5, 5000, 80], [2.0, 3000, 50], [2.5, 1500, 30], [3.0, 1000, 20]
  ];
  // TBE resolves small fragments a little better; TAE lets large ones run further.
  const BUFFERS = {
    TAE: { name: 'TAE', big: 1.15, small: 1.2, speed: 1.0 },
    TBE: { name: 'TBE', big: 0.85, small: 0.8, speed: 0.9 }
  };

  function rangeFor(pct) {
    if (pct <= RANGE[0][0]) return [RANGE[0][1], RANGE[0][2]];
    if (pct >= RANGE[RANGE.length - 1][0]) return [RANGE[RANGE.length - 1][1], RANGE[RANGE.length - 1][2]];
    for (let i = 1; i < RANGE.length; i++) {
      if (pct <= RANGE[i][0]) {
        const [p0, b0, s0] = RANGE[i - 1], [p1, b1, s1] = RANGE[i];
        const t = (pct - p0) / (p1 - p0);
        return [b0 + (b1 - b0) * t, s0 + (s1 - s0) * t];
      }
    }
    return [RANGE[0][1], RANGE[0][2]];
  }

  // Fraction of the run distance a fragment of `size` has travelled, 0 at the well.
  function mobility(size, opts) {
    const buf = BUFFERS[opts.buffer] || BUFFERS.TAE;
    let [big, small] = rangeFor(opts.percent);
    big *= buf.big; small *= buf.small;
    const run = (opts.run === undefined ? 0.75 : opts.run) * buf.speed;
    const hi = Math.log10(big * 2.2), lo = Math.log10(small / 2.6);
    const f = (hi - Math.log10(Math.max(1, size))) / (hi - lo);
    return Math.max(0, Math.min(1.06, f)) * run;
  }

  // name -> { brand, bands: [bp], bright: [bp], note }
  const LADDERS = {
    'neb-1kb': { brand: 'NEB', name: '1 kb DNA Ladder', bands: [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500], bright: [3000] },
    'neb-1kb-plus': { brand: 'NEB', name: '1 kb Plus DNA Ladder', bands: [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], bright: [1000, 3000] },
    'neb-2log': { brand: 'NEB', name: '2-Log DNA Ladder', bands: [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], bright: [3000, 1000] },
    'neb-100bp': { brand: 'NEB', name: '100 bp DNA Ladder', bands: [1517, 1200, 1000, 900, 800, 700, 600, 517, 400, 300, 200, 100], bright: [517] },
    'neb-50bp': { brand: 'NEB', name: '50 bp DNA Ladder', bands: [1350, 916, 766, 700, 650, 500, 350, 300, 250, 200, 150, 100, 50], bright: [350] },
    'thermo-1kb': { brand: 'Thermo', name: 'GeneRuler 1 kb', bands: [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 750, 500, 250], bright: [6000, 3000, 1000] },
    'thermo-1kb-plus': { brand: 'Thermo', name: 'GeneRuler 1 kb Plus', bands: [20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300, 200, 75], bright: [5000, 1500, 500] },
    'thermo-100bp': { brand: 'Thermo', name: 'GeneRuler 100 bp', bands: [1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], bright: [500] },
    'thermo-100bp-plus': { brand: 'Thermo', name: 'GeneRuler 100 bp Plus', bands: [3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100], bright: [500, 1500] },
    'thermo-mix': { brand: 'Thermo', name: 'GeneRuler DNA Ladder Mix', bands: [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1200, 1031, 900, 800, 700, 600, 500, 400, 300, 200, 100], bright: [6000, 3000, 1000] },
    'invitrogen-1kb-plus': { brand: 'Invitrogen', name: '1 Kb Plus DNA Ladder', bands: [12000, 11000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1650, 1000, 850, 650, 500, 400, 300, 200, 100], bright: [1650] },
    'promega-1kb': { brand: 'Promega', name: '1 kb DNA Ladder', bands: [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250], bright: [3000] },
    'lambda-hindiii': { brand: 'Classic', name: 'λ DNA / HindIII', bands: [23130, 9416, 6557, 4361, 2322, 2027, 564, 125], bright: [] },
    'lambda-ecori-hindiii': { brand: 'Classic', name: 'λ DNA / EcoRI + HindIII', bands: [21226, 5148, 4973, 4268, 3530, 2027, 1904, 1584, 1375, 947, 831, 564, 125], bright: [] },
    'phix-haeiii': { brand: 'Classic', name: 'φX174 / HaeIII', bands: [1353, 1078, 872, 603, 310, 281, 271, 234, 194, 118, 72], bright: [] }
  };

  function ladderList(custom) {
    const out = Object.entries(LADDERS).map(([key, l]) => ({ key, ...l, custom: false }));
    for (const [key, l] of Object.entries(custom || {})) out.push({ key: 'custom:' + key, brand: 'My ladders', name: key, bands: l.bands, bright: l.bright || [], custom: true });
    return out;
  }
  function ladderByKey(key, custom) {
    if (key && key.startsWith('custom:')) {
      const l = (custom || {})[key.slice(7)];
      return l ? { name: key.slice(7), brand: 'My ladders', bands: l.bands, bright: l.bright || [] } : null;
    }
    return LADDERS[key] || null;
  }
  // Accepts "10000, 8000, 6000" or one size per line, with an optional * to mark a
  // reference band, and returns the sizes in descending order.
  function parseLadder(text) {
    const bright = [], bands = [];
    for (const tok of String(text).split(/[\s,;]+/)) {
      if (!tok) continue;
      const star = tok.includes('*');
      const m = tok.replace(/\*/g, '').match(/^([\d.]+)\s*(kb|bp)?$/i);
      if (!m) throw new Error(`“${tok}” is not a band size.`);
      let v = parseFloat(m[1]);
      if (/kb/i.test(m[2] || '') || (!m[2] && v < 50 && v > 0)) v *= 1000;
      v = Math.round(v);
      if (!(v > 0)) throw new Error(`“${tok}” is not a band size.`);
      bands.push(v);
      if (star) bright.push(v);
    }
    if (!bands.length) throw new Error('Enter at least one band size.');
    bands.sort((a, b) => b - a);
    return { bands, bright };
  }

  root.GEL = { RANGE, BUFFERS, rangeFor, mobility, LADDERS, ladderList, ladderByKey, parseLadder };
})(typeof window !== 'undefined' ? window : globalThis);
