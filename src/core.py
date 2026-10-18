"""Comparison-checked DASCH plate photometry and dimming-event search.
Quality mask = daschlab's standard rejections + suspected defect + rejected blend + bad plate + multiple-exposure confusion."""
import json, math, os, statistics as st
from collections import defaultdict
from dasch import DATA, rows
JD2000 = 2451545.0
def year(jd): return 2000 + (jd - JD2000) / 365.25
BAD = (1<<6)|(1<<11)|(1<<12)|(1<<13)|(1<<15) | (1<<25)|(1<<27)|(1<<7) | (1<<8)|(1<<10)

def parse(ref):
    d = json.load(open(os.path.join(DATA, 'lc', f'{ref}.json')))
    det, lim = [], []
    for x in rows(d):
        jd = float(x['date_jd']); pid = f"{x['series']}{x['plate_number']}"
        if x['magcal_magdep']:
            af = int(x['aflags'] or 0)
            det.append(dict(jd=jd, m=float(x['magcal_magdep']), lim=float(x['limiting_mag_local'] or 0), series=x['series'], plate=pid,
                            ok=(x['reject_flag'] in ('0', '') and not (af & BAD))))
        else:
            lim.append(dict(jd=jd, lim=float(x['limiting_mag_local'] or 0), series=x['series'], plate=pid))
    return det, lim

def base(det):
    g = [p['m'] for p in det if p['ok']]
    if len(g) < 20: return None, None
    med = st.median(g); return med, st.median([abs(m - med) for m in g]) * 1.4826

def plates(ref, comp_refs):
    det, lim = parse(ref); med, mad = base(det)
    if med is None: return None
    pp = defaultdict(lambda: {'res': [], 'n': 0, 'ndet': 0})
    for c in comp_refs:
        cd, cl = parse(c); cmed, _ = base(cd)
        if cmed is None: continue
        for p in cd:
            q = pp[p['plate']]; q['n'] += 1
            if p['ok']: q['ndet'] += 1; q['res'].append(p['m'] - cmed)
        for p in cl:
            if p['lim'] > cmed + 0.3: pp[p['plate']]['n'] += 1
    out = []
    for p in det:
        q = pp.get(p['plate'], {'res': [], 'n': 0, 'ndet': 0}); zp = st.median(q['res']) if len(q['res']) >= 2 else None
        out.append(dict(jd=p['jd'], kind='D', ok=p['ok'] and zp is not None, dm=(p['m'] - (zp or 0)) - med, lim=p['lim'] - med,
                        series=p['series'], plate=p['plate'], cn=q['n'], cdet=q['ndet']))
    for p in lim:
        q = pp.get(p['plate'], {'res': [], 'n': 0, 'ndet': 0})
        out.append(dict(jd=p['jd'], kind='N', ok=True, dm=None, lim=p['lim'] - med, series=p['series'], plate=p['plate'], cn=q['n'], cdet=q['ndet']))
    out.sort(key=lambda r: r['jd'])
    good = [r['dm'] for r in out if r['kind'] == 'D' and r['ok']]
    c0 = st.median(good); sig = st.median([abs(g - c0) for g in good]) * 1.4826
    for r in out:
        if r['dm'] is not None: r['dm'] -= c0
        r['lim'] -= c0
    return dict(med=med + c0, sig=sig, ngood=len(good), plates=out)

def faint_plates(S, min_dm=0.5, nsig=3.5, lim_margin=0.5, min_cfrac=0.4):
    """Plates where the star looks faint: a good detection > max(min_dm, nsig*sig) below normal, or a
    non-detection on a plate that reached lim_margin below normal and on which >= min_cfrac of comparisons were seen."""
    th = max(min_dm, nsig * S['sig']); out = []
    for r in S['plates']:
        if r['kind'] == 'D' and r['ok'] and r['dm'] > th: out.append((r, r['dm']))
        elif r['kind'] == 'N' and r['lim'] > max(lim_margin, th) and r['cn'] >= 3 and r['cdet'] >= min_cfrac * r['cn']: out.append((r, r['lim']))
    return out

def events(S, gap=30, min_nights=3, min_series=1, **kw):
    """Group faint plates into events; an event is broken by any good normal detection in between."""
    th = max(kw.get('min_dm', 0.5), kw.get('nsig', 3.5) * S['sig'])
    fp = {id(r) for r, _ in faint_plates(S, **kw)}
    ev, cur = [], []
    for r in S['plates']:
        if id(r) in fp:
            if cur and r['jd'] - cur[-1]['jd'] > gap: ev.append(cur); cur = []
            cur.append(r)
        elif r['kind'] == 'D' and r['ok'] and abs(r['dm']) < max(0.25, 2.5 * S['sig']):
            if cur: ev.append(cur); cur = []
    if cur: ev.append(cur)
    out = []
    for e in ev:
        nights = len({int(r['jd'] + 0.5) for r in e}); ser = {r['series'] for r in e}
        if nights >= min_nights and len(ser) >= min_series: out.append(e)
    return out
