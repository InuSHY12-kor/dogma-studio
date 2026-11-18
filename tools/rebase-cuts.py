#!/usr/bin/env python3
"""Where each restriction enzyme cuts, from REBASE.

    python3 tools/rebase-cuts.py check     compare the program's table with REBASE
    python3 tools/rebase-cuts.py lines     print every usable enzyme as "name site top bot"
    python3 tools/rebase-cuts.py trade     the same, for catalogue names REBASE has not heard of

REBASE is the reference list of restriction enzymes and their cut positions. Its own
server sends a broken certificate chain, so this reads the same data from the copy
Biopython ships (pip install biopython), which is generated from REBASE.

The program writes a cut as two numbers counted from the start of the recognition
sequence: `top` is how many bases of the top strand come before the cut, `bot` the same
for the bottom strand. EcoRI is "GAATTC 1 5" — G^AATTC, with a four-base 5' overhang.
Biopython counts the top cut the same way and the bottom one back from the end of the
site, so `bot = len(site) + fst3`.

Two kinds of enzyme are left out, because one cut per site is all the program can
represent: those that cut on both sides of their site (BcgI, BaeI and the like) and
those that nick a single strand.
"""
import sys

try:
    from Bio.Restriction import AllEnzymes
except ImportError:
    sys.exit("biopython is not installed: pip install biopython")


def usable():
    """Every enzyme that makes one double-strand cut, as (name, site, top, bot)."""
    out = {}
    for e in AllEnzymes:
        name, site = str(e), e.site
        try:
            fst5, fst3, scd5, scd3 = e.fst5, e.fst3, e.scd5, e.scd3
        except Exception:
            continue
        if scd5 is not None or scd3 is not None:
            continue                       # cuts on both sides of its site
        if fst5 is None or fst3 is None:
            continue                       # no published cut position
        top, bot = fst5, len(site) + fst3
        if top < 0 or bot < 0:
            continue                       # cuts before its own site; the program counts from its start
        out[name] = (site, top, bot)
    return out


def table():
    """What dna.js currently says, read straight out of its ENZ_RAW block."""
    import os
    import re
    path = os.path.join(os.path.dirname(__file__), '..', 'app', 'dna.js')
    text = open(path, encoding='utf8').read()
    m = re.search(r'const ENZ_RAW = `(.*?)`;', text, re.S)
    if not m:
        sys.exit('could not find the enzyme table in dna.js')
    rows = {}
    for line in m.group(1).strip().split('\n'):
        parts = line.split()
        # "name site top bot", with a trailing + on the wider catalogue's entries
        if len(parts) in (4, 5):
            rows[parts[0]] = (parts[1], int(parts[2]), int(parts[3]))
    return rows


COMP = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')


def rc(site):
    return site.translate(COMP)[::-1]


def catalogue():
    """Every name a supplier lists, with the recognition sequence it publishes."""
    import json
    import os
    import re
    path = os.path.join(os.path.dirname(__file__), '..', 'app', 'enzyme-db.js')
    out = {}
    for line in open(path, encoding='utf8'):
        line = line.strip().rstrip(',')
        if not line.startswith('{"name":') or '"enzymes":' not in line:
            continue
        for e in json.loads(line)['enzymes']:
            if e.get('site'):
                out.setdefault(e['name'], e['site'])
    return out


def trade():
    """Cut positions for the old trade names.

    A catalogue still sells enzymes under names REBASE no longer lists — SstII for
    SacII, AatI for StuI. The catalogue says what each one recognises, and that is
    enough to place its cut, but only where it is not in doubt: every enzyme REBASE
    knows that reads the same sequence has to agree on where it falls. Four enzymes
    read GGCGCC and cut it in four different places, so a name that only says GGCGCC
    could be any of them and is left out. Nicking enzymes are left out too — they
    break one strand, which a cut position cannot describe.
    """
    import re
    reb = usable()
    ours = table()
    bysite = {}
    for name, (site, top, bot) in reb.items():
        bysite.setdefault(site, {}).setdefault((top, bot), []).append(name)

    from_name = re.compile(r'^(FastDigest®?|Anza™|QuickCut™|EZ-CleanCut™|Fast\s*Digest|FuniCut™|SureCut™|WarmStart®?|Speedy)\s+', re.I)

    def base(n):
        n = from_name.sub('', n.strip())
        n = re.sub(r'\s*[-\u2013]?\s*HF®?\b', '', n, flags=re.I)
        n = re.sub(r'\s*-?v\d+$', '', n, flags=re.I)
        return n.replace('®', '').replace('™', '').strip()

    wanted = {}
    for listed, site in catalogue().items():
        n = base(listed)
        if n in ours or n in reb:
            continue
        wanted.setdefault(n, site)

    out, unsure, unknown = [], [], []
    for name in sorted(wanted, key=str.lower):
        site = wanted[name]
        if name.startswith(('Nt.', 'Nb.')):
            unsure.append((name, site, 'nicks one strand'))
            continue
        table_for = bysite.get(site) or bysite.get(rc(site))
        key = site if site in bysite else rc(site)
        if not table_for:
            unknown.append((name, site))
        elif len(table_for) > 1:
            unsure.append((name, site, 'read %d different ways' % len(table_for)))
        else:
            (top, bot), who = next(iter(table_for.items()))
            out.append((name, key, top, bot, who[0]))
    return out, unsure, unknown


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else 'check'
    reb = usable()
    if what == 'trade':
        out, unsure, unknown = trade()
        for name, site, top, bot, like in out:
            print('%s %s %d %d    # cuts as %s does' % (name, site, top, bot, like))
        print('\n# left out, in doubt: %s' % ', '.join('%s (%s)' % (n, why) for n, _, why in unsure), file=sys.stderr)
        print('# left out, REBASE does not read that sequence: %s'
              % ', '.join(n for n, _ in unknown), file=sys.stderr)
        return
    if what == 'lines':
        for name in sorted(reb, key=str.lower):
            site, top, bot = reb[name]
            print('%s %s %d %d' % (name, site, top, bot))
        return

    ours = table()
    same = diff_site = diff_cut = 0
    for name, (site, top, bot) in sorted(ours.items()):
        if name not in reb:
            continue
        rsite, rtop, rbot = reb[name]
        if rsite != site:
            diff_site += 1
            print('site   %-10s ours %-14s rebase %s' % (name, site, rsite))
        elif (rtop, rbot) != (top, bot):
            diff_cut += 1
            print('cut    %-10s %-14s ours (%d/%d)  rebase (%d/%d)' % (name, site, top, bot, rtop, rbot))
        else:
            same += 1
    missing = [n for n in ours if n not in reb]
    fresh = [n for n in reb if n not in ours]
    print('\n%d of the program\'s %d enzymes agree with REBASE' % (same, len(ours)))
    print('%d differ in recognition sequence, %d in cut position' % (diff_site, diff_cut))
    print('%d of ours are not in this REBASE copy: %s' % (len(missing), ', '.join(sorted(missing))))
    print('%d enzymes in REBASE that the program does not have' % len(fresh))

    # a palindrome has to be cut symmetrically; anything else is an error in the table
    comp = str.maketrans('ACGTRYSWKMBDHVN', 'TGCAYRSWMKVHDBN')
    odd = [n for n, (s, t, b) in ours.items()
           if s.translate(comp)[::-1] == s and b != len(s) - t]
    print('%d of ours break the mirror symmetry a palindrome must have%s'
          % (len(odd), (': ' + ', '.join(sorted(odd))) if odd else ''))


main()
