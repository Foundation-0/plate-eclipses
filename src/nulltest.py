"""Null test: run the blind event search on every comparison star (each checked against the others).
Constant stars should show no events; this measures the false-alarm rate of events.py."""
import json, os
from dasch import DATA
from core import plates, events, year
xm = json.load(open(os.path.join(DATA, 'xmatch.json'))); comps = json.load(open(os.path.join(DATA, 'comps.json')))
nstar = nev = 0; hits = []
for name, cl in comps.items():
    refs = [c['ref_number'] for c in cl]
    for c in cl:
        others = [r for r in refs if r != c['ref_number']]
        try: S = plates(c['ref_number'], others)
        except Exception as e: continue
        if S is None: continue
        nstar += 1; ev = events(S)
        if ev:
            nev += len(ev); hits.append((name, c['ref_text'], c['stdmag'], [(round(year(e[0]['jd']),2), len(e)) for e in ev]))
print('comparison stars tested', nstar, 'stars with >=1 event', len(hits), 'events', nev)
for h in hits: print('  ', h)
