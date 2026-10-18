"""Step 2: blind search for historical dimming events (no period assumed). Prints events of >= 3 nights.
Usage: python events.py [target names...]   (default: all targets)"""
import json, os, sys
from dasch import DATA
from core import plates, events, year
from targets import T
xm = json.load(open(os.path.join(DATA, 'xmatch.json'))); comps = json.load(open(os.path.join(DATA, 'comps.json')))
for name in (sys.argv[1:] or [t[0] for t in T]):
    if name not in xm: continue
    S = plates(xm[name]['target']['ref_number'], [c['ref_number'] for c in comps[name]])
    if S is None: print(f"{name}: too few good detections"); continue
    ev = events(S)
    print(f"{name}: median B {S['med']:.2f}, plate scatter {S['sig']:.2f}, {S['ngood']} good plates, {len(ev)} events")
    for e in ev:
        print(f"   {year(e[0]['jd']):9.3f}-{year(e[-1]['jd']):9.3f} plates={len(e)} nights={len({int(r['jd']+.5) for r in e})} series={','.join(sorted({r['series'] for r in e}))} :: " +
              ' '.join((f"F{r['dm']:+.2f}" if r['kind']=='D' else f"N>{r['lim']:+.2f}") + f"[{r['cdet']}/{r['cn']}]" for r in e[:8]))
