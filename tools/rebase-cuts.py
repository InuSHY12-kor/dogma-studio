#!/usr/bin/env python3
"""Where each restriction enzyme cuts, from REBASE.

    python3 tools/rebase-cuts.py check     compare the program's table with REBASE
    python3 tools/rebase-cuts.py lines     print every usable enzyme as "name site top bot"

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


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else 'check'
    reb = usable()
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
