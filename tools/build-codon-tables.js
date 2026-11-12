#!/usr/bin/env node
/*
 * Builds app/codon-tables.js from real coding sequences at NCBI.
 *
 * Bacteria: every CDS of the reference genome.
 * Everything else: the RefSeq Select transcript of up to `LIMIT` genes, which is
 * one representative transcript per gene, so no gene is counted twice.
 *
 * Run with:  node tools/build-codon-tables.js
 */
'use strict';
const fs = require('fs');
const path = require('path');
const https = require('https');
const http = require('http');
const tls = require('tls');
const { URL } = require('url');

// Honour HTTPS_PROXY so the script also works from behind a proxy.
const PROXY = process.env.HTTPS_PROXY || process.env.https_proxy || '';

function request(url) {
  if (!PROXY) return new Promise((resolve, reject) => bind(https.get(url), resolve, reject));
  const target = new URL(url), prox = new URL(PROXY);
  return new Promise((resolve, reject) => {
    const req = http.request({
      host: prox.hostname, port: prox.port || 80, method: 'CONNECT',
      path: `${target.hostname}:443`,
      headers: prox.username ? { 'Proxy-Authorization': 'Basic ' + Buffer.from(`${prox.username}:${prox.password}`).toString('base64') } : {}
    });
    req.on('connect', (res, socket) => {
      if (res.statusCode !== 200) { socket.destroy(); return reject(new Error('proxy CONNECT returned ' + res.statusCode)); }
      const s = tls.connect({ socket, servername: target.hostname }, () => {
        bind(https.get({ host: target.hostname, path: target.pathname + target.search, createConnection: () => s }), resolve, reject);
      });
      s.on('error', reject);
    });
    req.on('error', reject);
    req.end();
  });
}
function bind(req, resolve, reject) {
  req.on('response', res => {
    if (res.statusCode !== 200) { res.resume(); return reject(new Error('HTTP ' + res.statusCode)); }
    let body = '';
    res.setEncoding('utf8');
    res.on('data', c => { body += c; });
    res.on('end', () => resolve(body));
  });
  req.setTimeout(180000, () => req.destroy(new Error('timeout')));
  req.on('error', reject);
  return req;
}

const EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
const TOOL = 'GeneCanvas';
const LIMIT = 600;          // genes per organism; codon usage converges well before this
const BATCH = 100;

const ORGANISMS = [
  { key: 'ecoli', name: 'Escherichia coli K-12', short: 'E. coli', genome: 'NC_000913.3' },
  { key: 'yeast', name: 'Saccharomyces cerevisiae', short: 'S. cerevisiae', organism: 'Saccharomyces cerevisiae' },
  { key: 'human', name: 'Homo sapiens', short: 'Human', organism: 'Homo sapiens' },
  { key: 'mouse', name: 'Mus musculus', short: 'Mouse', organism: 'Mus musculus' },
  { key: 'cho', name: 'Cricetulus griseus (CHO)', short: 'CHO', organism: 'Cricetulus griseus' },
  { key: 'pichia', name: 'Komagataella phaffii (Pichia pastoris)', short: 'P. pastoris', organism: 'Komagataella phaffii' },
  { key: 'insect', name: 'Drosophila melanogaster', short: 'Insect', organism: 'Drosophila melanogaster' },
  { key: 'plant', name: 'Arabidopsis thaliana', short: 'Plant', organism: 'Arabidopsis thaliana' }
];

let last = 0;
const sleep = ms => new Promise(r => setTimeout(r, ms));
async function get(url) {
  const wait = Math.max(0, 420 - (Date.now() - last));
  if (wait) await sleep(wait);
  last = Date.now();
  for (let attempt = 1; attempt <= 4; attempt++) {
    try {
      return await request(url);
    } catch (e) {
      if (attempt === 4) throw e;
      process.stderr.write(`  retry ${attempt} (${e.message})\n`);
      await sleep(1500 * attempt);
    }
  }
}
const q = o => Object.entries(o).map(([k, v]) => `${k}=${encodeURIComponent(v)}`).join('&');

async function search(term) {
  const s = await get(`${EUTILS}esearch.fcgi?${q({ db: 'nuccore', term, retmax: LIMIT, retmode: 'json', tool: TOOL })}`);
  return (JSON.parse(s).esearchresult || {}).idlist || [];
}
// RefSeq Select gives one transcript per gene but only exists for a few well-curated
// organisms, so fall back to plain RefSeq mRNA records for everything else.
async function searchIds(organism) {
  let ids = await search(`"${organism}"[Organism] AND refseq_select[filter] AND biomol_mrna[prop]`);
  if (ids.length) return { ids, how: 'RefSeq Select' };
  ids = await search(`"${organism}"[Organism] AND srcdb_refseq[prop] AND biomol_mrna[prop]`);
  return { ids, how: 'RefSeq mRNA' };
}
async function fetchCds(ids) {
  const out = [];
  for (let i = 0; i < ids.length; i += BATCH) {
    const chunk = ids.slice(i, i + BATCH);
    const text = await get(`${EUTILS}efetch.fcgi?${q({ db: 'nuccore', id: chunk.join(','), rettype: 'fasta_cds_na', retmode: 'text', tool: TOOL })}`);
    out.push(text);
    process.stderr.write(`  fetched ${Math.min(i + BATCH, ids.length)}/${ids.length}\n`);
  }
  return out.join('\n');
}

// Counts codons over every complete, in-frame CDS in a FASTA stream.
function countCodons(fasta) {
  const counts = Object.create(null);
  let genes = 0, skipped = 0, codons = 0;
  let cur = null;
  const flush = () => {
    if (cur === null) return;
    const s = cur.toUpperCase().replace(/[^ACGTU]/g, '').replace(/U/g, 'T');
    cur = null;
    // only clean, complete reading frames, so partial records cannot shift the counts
    if (s.length < 60 || s.length % 3 !== 0) { skipped++; return; }
    if (!/^(ATG|GTG|TTG)/.test(s)) { skipped++; return; }
    const stop = s.slice(-3);
    if (stop !== 'TAA' && stop !== 'TAG' && stop !== 'TGA') { skipped++; return; }
    let internalStop = false;
    for (let i = 3; i < s.length - 3; i += 3) {
      const c = s.substr(i, 3);
      if (c === 'TAA' || c === 'TAG' || c === 'TGA') { internalStop = true; break; }
    }
    if (internalStop) { skipped++; return; }
    genes++;
    // the initiator and the stop codon are not a free choice, so they are left out
    for (let i = 3; i < s.length - 3; i += 3) {
      const c = s.substr(i, 3);
      counts[c] = (counts[c] || 0) + 1;
      codons++;
    }
  };
  for (const line of fasta.split('\n')) {
    if (line.startsWith('>')) { flush(); cur = ''; }
    else if (cur !== null) cur += line.trim();
  }
  flush();
  return { counts, genes, skipped, codons };
}

(async () => {
  const tables = {};
  for (const o of ORGANISMS) {
    process.stderr.write(`\n${o.name}\n`);
    let fasta;
    if (o.genome) {
      process.stderr.write(`  genome ${o.genome}\n`);
      fasta = await get(`${EUTILS}efetch.fcgi?${q({ db: 'nuccore', id: o.genome, rettype: 'fasta_cds_na', retmode: 'text', tool: TOOL })}`);
    } else {
      const { ids, how } = await searchIds(o.organism);
      process.stderr.write(`  ${ids.length} transcripts (${how})\n`);
      if (!ids.length) { process.stderr.write('  SKIPPED: nothing found\n'); continue; }
      o.how = how;
      fasta = await fetchCds(ids);
    }
    const r = countCodons(fasta);
    if (r.codons < 20000) { process.stderr.write(`  SKIPPED: only ${r.codons} codons\n`); continue; }
    process.stderr.write(`  ${r.genes} genes used, ${r.skipped} skipped, ${r.codons} codons\n`);
    tables[o.key] = { name: o.name, short: o.short, source: o.genome ? `NCBI ${o.genome}, all CDS` : `NCBI ${o.how}, ${r.genes} records`, genes: r.genes, codons: r.codons, counts: r.counts };
  }

  const BASES = 'TCAG';
  const AAS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const order = [];
  for (let i = 0; i < 64; i++) order.push(BASES[i >> 4] + BASES[(i >> 2) & 3] + BASES[i & 3]);

  let out = `/* GeneCanvas — codon usage tables. Exposes window.CODON_TABLES.
 *
 * GENERATED FILE — do not edit by hand. Rebuild with:
 *     node tools/build-codon-tables.js
 *
 * Counts were computed from coding sequences downloaded from NCBI on ${new Date().toISOString().slice(0, 10)}.
 * Only complete, in-frame CDS with a start codon, a stop codon and no internal stop
 * were counted. The initiator and stop codons themselves are excluded, since neither
 * is a free synonymous choice.
 */
(function (root) {
  'use strict';
  var CODONS = ${JSON.stringify(order)};
  var T = {};
`;
  for (const [key, t] of Object.entries(tables)) {
    const counts = order.map(c => t.counts[c] || 0);
    out += `  T.${key} = { name: ${JSON.stringify(t.name)}, short: ${JSON.stringify(t.short)}, source: ${JSON.stringify(t.source)}, genes: ${t.genes}, codons: ${t.codons},\n    counts: ${JSON.stringify(counts)} };\n`;
  }
  out += `  root.CODON_TABLES = { codons: CODONS, tables: T };
})(typeof window !== 'undefined' ? window : globalThis);
`;
  const dest = path.join(__dirname, '..', 'app', 'codon-tables.js');
  fs.writeFileSync(dest, out);
  process.stderr.write(`\nwrote ${dest} with ${Object.keys(tables).length} tables\n`);
})().catch(e => { console.error('FAILED', e); process.exit(1); });
