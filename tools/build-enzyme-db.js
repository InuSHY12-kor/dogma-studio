#!/usr/bin/env node
/* Build app/enzyme-db.js from the three enzyme tables.
 *
 *   node tools/build-enzyme-db.js <folder with master.csv, buffers.csv, suppliers.csv>
 *
 * This is the reference side of the enzyme data — who sells which enzyme, what it
 * recognises, the temperature it works at, the buffer its supplier recommends and how
 * well it works in that supplier's other buffers. It is what the Enzyme Database window
 * shows, and it carries no cut positions: the enzymes the program actually cuts with
 * live in dna.js, where each one says where in its site the cut falls.
 *
 *   master.csv     Supplier, Enzyme, Recognition Sequence, Incubation Temp, Recommended Buffer
 *   buffers.csv    the same, plus Buffer Name and Activity (%)
 *   suppliers.csv  Enzyme, Suppliers (one comma-separated list)
 *
 * The tables are not part of the repository — the generated file is.
 */
'use strict';
const fs = require('fs');
const path = require('path');

const SRC = process.argv[2];
const OUT = process.argv[3] || path.join(__dirname, '..', 'app', 'enzyme-db.js');
if (!SRC) {
  console.error('usage: node tools/build-enzyme-db.js <folder of csv files> [out]');
  process.exit(2);
}

function readCsv(file) {
  const text = fs.readFileSync(file, 'utf8');
  const rows = [];
  let row = [], cur = '', quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"') { if (text[i + 1] === '"') { cur += '"'; i++; } else quoted = false; }
      else cur += c;
    } else if (c === '"') quoted = true;
    else if (c === ',') { row.push(cur); cur = ''; }
    else if (c === '\n') { row.push(cur); rows.push(row); row = []; cur = ''; }
    else if (c !== '\r') cur += c;
  }
  if (cur || row.length) { row.push(cur); rows.push(row); }
  return rows.filter(r => r.some(c => c !== ''));
}

const tidy = s => String(s == null ? '' : s).replace(/ /g, ' ').trim();
// "37 °C" is the same fact as "37°C"; the window adds the unit itself
const temp = s => tidy(s).replace(/\s*°?\s*C$/i, '').trim();

const master = readCsv(path.join(SRC, 'master.csv')).slice(1);
const buffers = readCsv(path.join(SRC, 'buffers.csv')).slice(1);
const suppliers = readCsv(path.join(SRC, 'suppliers.csv')).slice(1);

/* One entry per supplier and enzyme. The activity table is keyed the same way, so the
 * two are joined on that pair rather than on the row order, which differs between
 * the two tables. */
const bySupplier = new Map();
const keyOf = (sup, enz) => sup + '\u0000' + enz;
const entries = new Map();

for (const r of master) {
  const [sup, enz, site, t, buf] = r.map(tidy);
  if (!sup || !enz) continue;
  const k = keyOf(sup, enz);
  if (!entries.has(k)) entries.set(k, { sup, enz, site, temp: temp(t), buffer: buf, act: {} });
}
let orphans = 0;
const bufNames = new Map();          // supplier -> the buffers it sells, in the order seen
for (const r of buffers) {
  const [sup, enz, site, t, buf, bufName, act] = r.map(tidy);
  if (!sup || !enz || !bufName) continue;
  const k = keyOf(sup, enz);
  if (!entries.has(k)) { entries.set(k, { sup, enz, site, temp: temp(t), buffer: buf, act: {} }); orphans++; }
  const e = entries.get(k);
  if (!e.site && site) e.site = site;
  if (!e.temp && t) e.temp = temp(t);
  if (!e.buffer && buf) e.buffer = buf;
  if (act) e.act[bufName] = act;
  if (!bufNames.has(sup)) bufNames.set(sup, []);
  const list = bufNames.get(sup);
  if (!list.includes(bufName)) list.push(bufName);
}

for (const e of entries.values()) {
  if (!bySupplier.has(e.sup)) bySupplier.set(e.sup, []);
  bySupplier.get(e.sup).push(e);
}

const supplierList = [...bySupplier.keys()].sort((a, b) => a.localeCompare(b)).map(name => {
  const bufs = (bufNames.get(name) || []).slice();
  const rows = bySupplier.get(name)
    .sort((a, b) => a.enz.localeCompare(b.enz, undefined, { numeric: true }))
    .map(e => {
      const o = { name: e.enz, site: e.site, temp: e.temp };
      if (e.buffer) o.buffer = e.buffer;
      if (bufs.length) o.act = bufs.map(b => e.act[b] || '');
      return o;
    });
  return { name, buffers: bufs, enzymes: rows };
});

const byEnzyme = suppliers
  .map(r => [tidy(r[0]), tidy(r[1])])
  .filter(([n, s]) => n && s)
  .map(([n, s]) => {
    // the list keeps its order but not its repeats: a name listed twice is one supplier
    const seen = [];
    for (const part of s.split(',').map(tidy)) if (part && !seen.includes(part)) seen.push(part);
    return [n, seen];
  });

const head = `/* Dogma Studio — the enzyme reference.
 *
 * Who sells which restriction enzyme, what it recognises, the temperature it works at,
 * the buffer its supplier recommends, and how well it works in that supplier's other
 * buffers. This is what the Enzyme Database window shows.
 *
 * It holds no cut positions, on purpose. The enzymes the program cuts with live in
 * dna.js, where each one records where in its site the cut falls; a name here that is
 * not there is a name the program can tell you about but cannot digest with.
 *
 * Generated by tools/build-enzyme-db.js; it is not edited by hand.
 *
 * ${supplierList.length} suppliers, ${entries.size} supplier-and-enzyme entries, ${byEnzyme.length} enzymes listed by supplier.
 */
(function (root) {
  'use strict';
  root.ENZYME_DB = {
    suppliers: [
`;
const body = supplierList.map(s =>
  '      ' + JSON.stringify({ name: s.name, buffers: s.buffers, enzymes: s.enzymes })
).join(',\n');
const tail = `
    ],
    byEnzyme: [
${byEnzyme.map(e => '      ' + JSON.stringify(e)).join(',\n')}
    ]
  };
})(typeof window !== 'undefined' ? window : globalThis);
`;
fs.writeFileSync(OUT, head + body + tail);

console.log(`${supplierList.length} suppliers, ${entries.size} entries, ${byEnzyme.length} enzymes by supplier` +
  (orphans ? `, ${orphans} seen only in the activity table` : '') +
  ` -> ${path.relative(process.cwd(), OUT)} (${(fs.statSync(OUT).size / 1024).toFixed(0)} KB)`);
