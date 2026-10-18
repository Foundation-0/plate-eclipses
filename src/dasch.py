"""Minimal client for the public DASCH DR7 API (no key needed). Light curves are cached as JSON under DATA/lc."""
import json, os, time, urllib.request
API = 'https://api.starglass.cfa.harvard.edu/public/dasch/dr7/'
DATA = os.environ.get('DASCH_PILOT_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data'))

def call(ep, p, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(API + ep, data=json.dumps(p).encode(),
                headers={'Content-Type': 'application/json', 'User-Agent': 'foundation0-plate-eclipses'})
            return json.load(urllib.request.urlopen(req, timeout=120))
        except Exception:
            if k == tries - 1: raise
            time.sleep(2 ** (k + 1))

def rows(csvlines):
    h = csvlines[0].split(',')
    return [dict(zip(h, x.split(','))) for x in csvlines[1:]]

def querycat(ra, dec, rad):
    return rows(call('querycat', {"refcat": "apass", "ra_deg": ra, "dec_deg": dec, "radius_arcsec": rad}))

def lightcurve(rec):
    f = os.path.join(DATA, 'lc', f"{rec['ref_number']}.json")
    if os.path.exists(f): return json.load(open(f))
    d = call('lightcurve', {"refcat": "apass", "gsc_bin_index": int(rec['gsc_bin_index']), "ref_number": int(rec['ref_number'])})
    os.makedirs(os.path.dirname(f), exist_ok=True); json.dump(d, open(f, 'w'))
    return d
