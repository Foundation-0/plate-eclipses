"""Fetch DASCH light curves for named targets only (companion to fetch.py, which does all 31).
Merges into DATA/xmatch.json and DATA/comps.json so existing steps can use the new stars.
Usage: python fetch_targets.py J073234-200049 [more names...]"""
import json, math, os, sys
from concurrent.futures import ThreadPoolExecutor
from dasch import querycat, lightcurve, DATA
from targets import T
from fetch import pick

if __name__ == '__main__':
    want = sys.argv[1:]
    byname = {t[0]: t for t in T}
    os.makedirs(DATA, exist_ok=True)
    xf, cf = os.path.join(DATA, 'xmatch.json'), os.path.join(DATA, 'comps.json')
    xm = json.load(open(xf)) if os.path.exists(xf) else {}
    comps = json.load(open(cf)) if os.path.exists(cf) else {}
    jobs = []
    for name in want:
        _, ra, dec, *_ = byname[name]
        near = querycat(ra, dec, 90)
        for r in near: r['sep'] = math.hypot(float(r['dra_asec']), float(r['ddec_asec']))
        near.sort(key=lambda r: r['sep'])
        if not near or near[0]['sep'] > 5: print(name, 'no DASCH match'); continue
        t = near[0]; xm[name] = {'target': t}
        c = pick(querycat(ra, dec, 300), t, 0.6, 0.5)
        if len(c) < 5: c = pick(querycat(ra, dec, 600), t, 1.0, 0.25, 50)
        comps[name] = c; jobs += [t] + c
        print(f"{name}: {t['ref_text']} B={t['stdmag']} plates={t['num_matches']} comparisons={len(c)}", flush=True)
    json.dump(xm, open(xf, 'w'), indent=1); json.dump(comps, open(cf, 'w'), indent=1)
    with ThreadPoolExecutor(3) as ex: list(ex.map(lightcurve, jobs))
    print('light curves:', len(jobs))
