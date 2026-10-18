"""For each dipper: fraction of 1890-1990 in which a dimming like the modern one (depth, duration)
would have been caught, i.e. the dimming window would contain a good plate showing normal brightness
(a 'normal' plate is one within depth/2 of the median)."""
import json, os
import numpy as np
from targets import T
from dasch import DATA
from core import plates, year
xm = json.load(open(os.path.join(DATA, 'xmatch.json'))); comps = json.load(open(os.path.join(DATA, 'comps.json')))
J0, J1 = 2411368.5, 2447892.5   # 1890.0 - 1990.0
res = {}
for name, ra, dec, cls, g, ndip, dur, pmin, last, depth in T:
    S = plates(xm[name]['target']['ref_number'], [c['ref_number'] for c in comps[name]])
    if S is None: continue
    if depth < 3 * S['sig']:
        res[name] = None; continue
    normal = np.array(sorted(r['jd'] for r in S['plates'] if r['kind'] == 'D' and r['ok'] and r['dm'] < depth / 2))
    for D in sorted({dur if dur >= 1 else 1, 30, 100, 365}):
        t = np.arange(J0, J1, 1.0)
        i = np.searchsorted(normal, t)          # first normal plate at or after t
        caught = (i < len(normal)) & (normal[np.minimum(i, len(normal) - 1)] <= t + D)
        res.setdefault(name, {})[str(int(D))] = round(float(caught.mean()), 3)
    res[name]['n_normal'] = int(len(normal)); res[name]['depth'] = depth; res[name]['dur'] = dur; res[name]['class'] = cls
    print(f"{name} {cls:3s} depth {depth:4.2f} dur {dur:7.1f} sig {S['sig']:.2f} normal plates {len(normal):4d} | caught fraction: " + ', '.join(f"{k} d: {v:.0%}" for k, v in res[name].items() if k.isdigit()))
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'coverage.json'), 'w'), indent=1)
