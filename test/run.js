/* Both passes, one after the other: the PC build and the phone build.
 * This is what to run before calling a change done. */
'use strict';
const { spawnSync } = require('child_process');
const path = require('path');

let bad = 0;
for (const [name, file] of [['desktop', 'desktop.js'], ['phone', 'mobile.js']]) {
  console.log(`\n──────── ${name} ────────`);
  const r = spawnSync(process.execPath, [path.join(__dirname, file)], { stdio: 'inherit' });
  if (r.status !== 0) bad++;
}
console.log(bad ? `\n${bad} pass(es) reported errors.` : '\nBoth passes are clean.');
process.exit(bad ? 1 : 0);
