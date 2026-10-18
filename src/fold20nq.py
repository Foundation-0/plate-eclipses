"""Trials-free folded test for ASASSN-20nq (J005437+644347) at its published period.
JoHantgen et al. (arXiv:2507.19594) found P = 752.3 d from four identical modern dips
(~0.6 mag, ~8 d). No period is searched here: every scored plate that falls within
+-4.5 d of a predicted mid-eclipse (anchor: 2024 dip start 01/12/2024 + 4.5 d) is listed.
Usage: python fold20nq.py"""
import json, math, os
import numpy as np
from dasch import DATA
from core import plates

xm = json.load(open(os.path.join(DATA, 'xmatch.json'))); comps = json.load(open(os.path.join(DATA, 'comps.json')))
S = plates(xm['J005437+644347']['target']['ref_number'], [c['ref_number'] for c in comps['J005437+644347']])
P, T0, D = 752.3, 2460650.0, 9.0
for lab, minc in [("repo cuts (core ok)", 0), ("strict (cdet >= 4 comparison residuals)", 4)]:
    det = [r for r in S['plates'] if r['kind'] == 'D' and r['ok'] and r['cdet'] >= minc]
    nond = [r for r in S['plates'] if r['kind'] == 'N' and r['lim'] > 0.5 and r['cn'] >= 3 and r['cdet'] >= 0.4 * r['cn']]
    dm = np.array([r['dm'] for r in det]); sig = 1.4826 * float(np.median(np.abs(dm - np.median(dm))))
    print(f"--- {lab}: {len(det)} detections, {len(nond)} informative non-detections, plate sigma {sig:.3f}")
    hits = []
    for r in det + nond:
        ph = (T0 - r['jd']) / P; k = round(ph); dist = abs((ph - k) * P)
        if dist < D / 2 and k >= 1:
            hits.append((r['plate'], round(r['jd'], 1), k, round(dist, 1), r['kind'], None if r['dm'] is None else round(r['dm'], 2), round(r['lim'], 2), f"{r['cdet']}/{r['cn']}"))
    print("   plates in predicted eclipses (plate, jd, cycle back, dist to mid d, kind, dB, lim, comparisons):")
    for h in sorted(hits, key=lambda x: x[1]): print("   ", h)
    print(f"   expected detections in windows from plate placement alone: {len(det) * D / P:.1f}; observed: {sum(1 for h in hits if h[4] == 'D')}")
    faint = [(r['plate'], round(r['jd'], 1), round(r['dm'], 2), round(abs(((T0 - r['jd']) / P - round((T0 - r['jd']) / P)) * P), 1)) for r in det if r['dm'] > 0.4]
    print(f"   all faint detections (dB > 0.4) and their distance to the nearest predicted mid-eclipse: {faint}")
print("conclusion: the two plates in predicted windows (+0.27, +0.28 mag, same cycle, ~1.7 sigma each against a 0.6 mag modern depth) do not confirm the period on plates; no faint plate lands near any other predicted eclipse.")
