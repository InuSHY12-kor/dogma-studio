#!/usr/bin/env node
/* Puts the version from package.json into the Android project, so the APK reports the
 * same version as the desktop build instead of Capacitor's placeholder 1.0.
 *
 * versionName is the version as people read it (0.3.0); versionCode is the integer
 * Android orders updates by, built from the same three numbers (0.3.0 -> 300) so it
 * rises whenever the version does.
 */
const fs = require('fs');
const path = require('path');

const root = path.join(__dirname, '..');
const version = require(path.join(root, 'package.json')).version;
const m = /^(\d+)\.(\d+)\.(\d+)/.exec(version);
if (!m) { console.error(`package.json has no three-part version: ${version}`); process.exit(1); }
const code = (+m[1]) * 10000 + (+m[2]) * 100 + (+m[3]);

const file = path.join(root, 'android', 'app', 'build.gradle');
if (!fs.existsSync(file)) { console.error('No android/app/build.gradle — run "npx cap add android" first.'); process.exit(1); }
const before = fs.readFileSync(file, 'utf8');
const after = before
  .replace(/versionCode\s+\d+/, `versionCode ${code}`)
  .replace(/versionName\s+"[^"]*"/, `versionName "${version}"`);
if (after === before) {
  console.log(`android/app/build.gradle already at ${version} (${code})`);
} else {
  fs.writeFileSync(file, after);
  console.log(`android/app/build.gradle set to versionName ${version}, versionCode ${code}`);
}
