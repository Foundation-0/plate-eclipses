"""Step 1: match each target to the DASCH APASS catalogue, pick comparison stars, download all light curves.
Comparison stars: from a 300" catalogue query (the API returns stars somewhat beyond the nominal radius; all are kept),
at least 20" away, within 0.6 mag of the target, with at least half as many plate matches; where that gives fewer
than 5, a 600" query with 1.0 mag and a quarter of the matches (at least 50). Up to 12, nearest first."""
import json, math, os
from concurrent.futures import ThreadPoolExecutor
from dasch import querycat, lightcurve, DATA
from targets import T

def pick(cat, t, dmag, frac, nmin=0):
    B = float(t['stdmag']) if t['stdmag'] != '99.9' else 14.5; out = []
    for r in cat:
        if r['ref_number'] == t['ref_number'] or r['stdmag'] == '99.9': continue
        sep = math.hypot(float(r['dra_asec']), float(r['ddec_asec']))
        if sep >= 20 and abs(float(r['stdmag']) - B) <= dmag and int(r['num_matches']) >= max(nmin, frac * int(t['num_matches'])):
            r['sep'] = sep; out.append(r)
    out.sort(key=lambda r: (abs(float(r['stdmag']) - B) > 0.6, r['sep']))
    return out[:12]

if __name__ == '__main__':
    os.makedirs(DATA, exist_ok=True); xm, comps, jobs = {}, {}, []
    for name, ra, dec, *_ in T:
        near = querycat(ra, dec, 90)
        for r in near: r['sep'] = math.hypot(float(r['dra_asec']), float(r['ddec_asec']))
        near.sort(key=lambda r: r['sep'])
        if not near or near[0]['sep'] > 5: print(name, 'no DASCH match'); continue
        t = near[0]; xm[name] = {'target': t}
        c = pick(querycat(ra, dec, 300), t, 0.6, 0.5)
        if len(c) < 5: c = pick(querycat(ra, dec, 600), t, 1.0, 0.25, 50)
        comps[name] = c; jobs += [t] + c
        print(f"{name}: {t['ref_text']} B={t['stdmag']} plates={t['num_matches']} comparisons={len(c)}", flush=True)
    json.dump(xm, open(os.path.join(DATA, 'xmatch.json'), 'w'), indent=1)
    json.dump(comps, open(os.path.join(DATA, 'comps.json'), 'w'), indent=1)
    with ThreadPoolExecutor(3) as ex: list(ex.map(lightcurve, jobs))
    print('light curves:', len(jobs))
