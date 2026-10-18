"""WAVE2: real (not proxy) DASCH search on the top Tier-A targets of TARGETS-V2.

Method = the ASASSN-21co template, unchanged thresholds:
  fetch target + up to 12 comparison stars (fetch.py `pick` rule verbatim),
  per-plate comparison zero-point (core.plates), blind event search
  (core.events: faint = detection > max(0.5, 3.5 sigma) or qualifying
  non-detection; event = >=3 nights, broken by a normal plate).

Tier-A periods are unknown, so the blind search is primary (no period scan
is possible without a modern epoch; none is invented). Every candidate
event (>=3 nights) gets the ASASSN-21co sanity gauntlet:
  * per-plate table: plate ID, JD, date, dm or limit, comparisons detected;
  * dlogL of a constant-depth box over the event window vs no event
    (clipped-Gaussian ll as asassn21co.py; non-detections scored as
    period.py: log(1/eps) when the limit is shallower than depth-0.2);
    window = event span extended to the bracketing good normal plates,
    so normal plates inside the window penalise the box;
  * shifted-epoch null: the same-duration box placed at 200 seeded random
    epochs over the plate baseline; fraction reaching the candidate dlogL;
  * comparison-star null: core.events run on each comparison star.
Sensitivity per target: detection threshold, catalogue depth vs threshold,
and caught fractions (coverage.py definition, corrected normal plates)
plus the stricter fraction of windows containing >=3 good-plate nights
(what the blind search actually needs).

Controls (run first, must recover the known results before targets count):
  ASASSN-24fw = J070519+061219 blind events (published 1937 + 1981);
  ASASSN-21co = J181752-580749 blind events (known: none) + period scan
  anchored on the 2021 eclipse (period.py settings of the README command,
  40 shifted-epoch nulls; known best peak 4343.5 d, dlogL 18.0).

Usage (from src/):
  SSL_CERT_FILE=/etc/ssl/cert.pem python3 wave2_search.py [--only NAME ...]
Outputs: ../results/wave2/<slug>.json + .md per star, summary.md/json.
Light curves cache in ../data/lc (git-ignored), shared with fetch.py.
"""
import argparse, csv, datetime as dtm, json, math, os, re, sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dasch import DATA, querycat, lightcurve          # noqa: E402
from core import plates, events, faint_plates, year   # noqa: E402
from fetch import pick                                # noqa: E402
from period import scan, peaks                        # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
OUT = os.path.join(ROOT, 'results', 'wave2')
J0, J1 = 2411368.5, 2447892.5   # 1890.0-1990.0, as coverage.py
CONTROLS = [dict(name='J070519+061219 (ASASSN-24fw control)', ra=106.328990, dec=6.205306,
                 depth=3.84, kind='control_blind'),
            dict(name='J181752-580749 (ASASSN-21co control)', ra=274.4647083, dec=-58.13033333,
                 depth=0.7, kind='control_period')]


def slug(name):
    return re.sub(r'[^A-Za-z0-9]+', '_', name).strip('_').lower()


def jd2date(j):
    return dtm.date.fromordinal(int(j - 1721424.5)).isoformat()


def tier_a_targets(n=8):
    with open(os.path.join(ROOT, 'results', 'targets_v2.csv')) as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith('#')))
    # CSV rank order is authoritative: proxy-scored rows first by score, then
    # PRIOR-only rows (whose nominal score can exceed proxy scores). Do NOT
    # re-sort by score alone (that promoted no-DASCH-match rows; fixed here).
    a = [r for r in rows if r['tier'] == 'A']
    a.sort(key=lambda r: int(r['rank']))
    return [dict(name=r['name'], ra=float(r['ra_deg']), dec=float(r['dec_deg']),
                 depth=float(r['depth_mag']), kind='tierA', rank=r['rank'],
                 score=r['score'], dasch_ref_csv=r['dasch_ref'],
                 type=r['type'], period_status=r['period_status']) for r in a[:n]]


def fetch_star(tg):
    """fetch.py matching verbatim, except the nearest-match radius is 15"
    (VSX coordinates; targets_v2.py precedent) and the separation recorded."""
    near = querycat(tg['ra'], tg['dec'], 90)
    for r in near:
        r['sep'] = math.hypot(float(r['dra_asec']), float(r['ddec_asec']))
    near.sort(key=lambda r: r['sep'])
    if not near or near[0]['sep'] > 15:
        raise RuntimeError(f"no DASCH APASS match within 15 arcsec (nearest "
                           f"{near[0]['sep']:.1f} arcsec)" if near else "querycat empty")
    t = near[0]
    comps = pick(querycat(tg['ra'], tg['dec'], 300), t, 0.6, 0.5)
    if len(comps) < 5:
        comps = pick(querycat(tg['ra'], tg['dec'], 600), t, 1.0, 0.25, 50)
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(lightcurve, [t] + comps))
    return t, comps


def box_dlogl(S, jd_lo, jd_hi, depth, eps):
    """dlogL (box of given depth on [jd_lo, jd_hi]) vs no event, over plates
    in the bracketing window [prev normal, next normal]. Detections:
    clipped Gaussian (asassn21co.py ll). Qualifying non-detections
    (cn>=3, cdet>=40%): period.py rule, log(1/eps) if lim < depth-0.2."""
    sig = math.hypot(S['sig'], 0.05)
    normals = [r['jd'] for r in S['plates'] if r['kind'] == 'D' and r['ok']
               and abs(r['dm']) < max(0.25, 2.5 * S['sig'])]
    lo = max([j for j in normals if j < jd_lo], default=jd_lo - 1)
    hi = min([j for j in normals if j > jd_hi], default=jd_hi + 1)
    ll = lambda r: -0.5 * min((r / sig) ** 2, 16.0)  # noqa: E731
    d = 0.0
    for r in S['plates']:
        if not (lo <= r['jd'] <= hi):
            continue
        inside = jd_lo <= r['jd'] <= jd_hi
        if r['kind'] == 'D' and r['ok']:
            d += ll(r['dm'] - (depth if inside else 0.0)) - ll(r['dm'])
        elif r['kind'] == 'N' and r['cn'] >= 3 and r['cdet'] >= 0.4 * r['cn']:
            if inside and r['lim'] < depth - 0.2:
                d += math.log(1.0 / eps)
    return d


def miss_rate(S):
    ok = [r for r in S['plates'] if r['cn'] >= 3 and r['cdet'] >= 0.4 * r['cn']]
    deep = [r for r in ok if (r['kind'] == 'D' and r['ok']) or (r['kind'] == 'N' and r['lim'] > 0.5)]
    return max(0.02, sum(r['kind'] == 'N' for r in deep) / max(1, len(deep)))


def gauntlet(S, ev, cat_depth, eps, rng):
    jd_lo, jd_hi = ev[0]['jd'], ev[-1]['jd']
    det_dm = [r['dm'] for r in ev if r['kind'] == 'D']
    depth_fit = float(np.median(det_dm)) if det_dm else cat_depth
    d_fit = box_dlogl(S, jd_lo, jd_hi, depth_fit, eps)
    d_cat = box_dlogl(S, jd_lo, jd_hi, cat_depth, eps)
    dur = max(jd_hi - jd_lo, 1.0)
    jds = [r['jd'] for r in S['plates']]
    null = []
    for _ in range(200):
        c = rng.uniform(min(jds) + dur, max(jds) - dur)
        null.append(box_dlogl(S, c - dur / 2, c + dur / 2, depth_fit, eps))
    null = np.array(null)
    return dict(jd_start=round(jd_lo, 2), jd_end=round(jd_hi, 2),
                date_start=jd2date(jd_lo), date_end=jd2date(jd_hi),
                span_d=round(jd_hi - jd_lo, 1), n_plates=len(ev),
                n_nights=len({int(r['jd'] + 0.5) for r in ev}),
                series=sorted({r['series'] for r in ev}),
                depth_fit_mag=round(depth_fit, 2),
                dlogL_box_fitdepth=round(float(d_fit), 1),
                dlogL_box_catdepth=round(float(d_cat), 1),
                shifted_epoch_null=dict(n=200, median=round(float(np.median(null)), 1),
                                        p95=round(float(np.percentile(null, 95)), 1),
                                        max=round(float(null.max()), 1),
                                        frac_ge_candidate=round(float(np.mean(null >= d_fit)), 3)),
                plates=[dict(plate=r['plate'], jd=round(r['jd'], 2), date=jd2date(r['jd']),
                             kind=r['kind'], dm=(round(r['dm'], 2) if r['dm'] is not None else None),
                             lim=round(r['lim'], 2), comparisons_detected=f"{r['cdet']}/{r['cn']}")
                        for r in ev])


def sensitivity(S, cat_depth):
    th = max(0.5, 3.5 * S['sig'])
    normal = np.array(sorted(r['jd'] for r in S['plates']
                             if r['kind'] == 'D' and r['ok'] and r['dm'] < cat_depth / 2))
    good_nights = np.array(sorted({int(r['jd'] + 0.5) for r in S['plates']
                                   if r['kind'] == 'D' and r['ok']}), float)
    out = dict(plate_sigma_mag=round(S['sig'], 3), faint_threshold_mag=round(th, 2),
               catalogue_depth_mag=cat_depth,
               catalogue_depth_detectable=bool(cat_depth > th),
               n_normal_plates=int(len(normal)), caught_frac={}, caught_frac_3nights={})
    t = np.arange(J0, J1, 1.0)
    for D in (10, 30, 100, 365):
        i = np.searchsorted(normal, t)
        caught = (i < len(normal)) & (normal[np.minimum(i, len(normal) - 1)] <= t + D)
        out['caught_frac'][str(D)] = round(float(caught.mean()), 3)
        lo = np.searchsorted(good_nights, t); hi = np.searchsorted(good_nights, t + D)
        out['caught_frac_3nights'][str(D)] = round(float(np.mean(hi - lo >= 3)), 3)
    return out


def comp_null(comp_refs):
    nstar, hits = 0, []
    for c in comp_refs:
        try:
            Sc = plates(c, [r for r in comp_refs if r != c])
        except Exception:  # noqa: BLE001 - recorded as untested below via count gap
            continue
        if Sc is None:
            continue
        nstar += 1
        ev = events(Sc)
        if ev:
            hits.append(dict(ref=c, n_events=len(ev),
                             events=[dict(year_start=round(year(e[0]['jd']), 2), n_plates=len(e))
                                     for e in ev]))
    return dict(comparisons_tested=nstar, with_events=len(hits),
                total_events=sum(h['n_events'] for h in hits), hits=hits)


def period_control(S):
    """period.py prepare/scan with the README ASASSN-21co settings + 40 nulls."""
    ok = [r for r in S['plates'] if r['cn'] >= 3 and r['cdet'] >= 0.4 * r['cn']]
    eps = miss_rate(S)
    obs = [(r['jd'], 'D', r['dm']) for r in S['plates'] if r['kind'] == 'D' and r['ok']]
    obs += [(r['jd'], 'N', r['lim']) for r in ok if r['kind'] == 'N']
    T0, D, depth, Pmin = 2459289.5, 66.0, 0.7, 3000.0
    Ps, L, NN = scan(obs, S['sig'], eps, T0, D, depth, Pmin)
    pk = peaks(Ps, L, NN)
    rng = np.random.default_rng(1); null = []
    for _ in range(40):
        off = rng.uniform(3 * D, 0.5 * Pmin) * rng.choice([-1, 1])
        _, L2, N2 = scan(obs, S['sig'], eps, T0 + off, D, depth, Pmin)
        m = N2 >= 2
        null.append(float(L2[m].max()) if m.any() else 0.0)
    null = np.array(null); best = pk[0][1] if pk else 0.0
    return dict(settings=dict(T0=T0, D=D, depth=depth, Pmin=Pmin, Pmax=45000),
                sigma=round(S['sig'], 3), miss_rate=round(float(eps), 3),
                n_obs=len(obs), n_trial_periods=len(Ps),
                peaks=[list(p) for p in pk],
                null40=dict(median=round(float(np.median(null)), 1),
                            max=round(float(null.max()), 1),
                            frac_ge_best=round(float(np.mean(null >= best)), 2)),
                known='best peak 4343.5 d, dlogL 18.0, 3 nights; null40 max 15.4 '
                      '(results/period_scan_asassn21co.txt)')


def run_star(tg, rng):
    res = dict(name=tg['name'], kind=tg['kind'], ra_deg=tg['ra'], dec_deg=tg['dec'],
               catalogue_depth_mag=tg['depth'])
    for k in ('rank', 'score', 'dasch_ref_csv', 'type', 'period_status'):
        if k in tg:
            res[k] = tg[k]
    try:
        t, comps = fetch_star(tg)
    except Exception as e:  # noqa: BLE001 - task: skip and record why
        res.update(status='fetch_failed', error=f'{type(e).__name__}: {e}')
        return res
    res.update(dasch_ref=t['ref_text'], dasch_ref_number=t['ref_number'],
               dasch_sep_arcsec=round(float(t['sep']), 2), dasch_B=t['stdmag'],
               dasch_num_matches=t['num_matches'], n_comparisons=len(comps),
               comparisons=[c['ref_text'] for c in comps])
    if 'dasch_ref_csv' in res:
        res['ref_matches_targets_v2'] = bool(res['dasch_ref_csv'] == t['ref_text'])
    comp_refs = [c['ref_number'] for c in comps]
    S = plates(t['ref_number'], comp_refs)
    if S is None:
        res.update(status='too_few_good_detections')
        return res
    eps = miss_rate(S)
    ev = events(S)
    sub = events(S, min_nights=1)   # sub-threshold context, NOT candidates
    res.update(status='searched', median_B_corrected=round(S['med'], 2),
               plate_sigma_mag=round(S['sig'], 3), n_good_plates=S['ngood'],
               n_plate_rows=len(S['plates']), miss_rate=round(float(eps), 3),
               jd_range=[round(min(r['jd'] for r in S['plates']), 1),
                         round(max(r['jd'] for r in S['plates']), 1)],
               n_faint_plates=len(faint_plates(S)),
               n_events_ge3nights=len(ev), n_clusters_ge1night=len(sub),
               sensitivity=sensitivity(S, tg['depth']))
    ev_sorted = sorted(ev, key=lambda e: (-len({int(r['jd'] + 0.5) for r in e}), -len(e)))
    res['gauntlet_capped_at_25'] = len(ev) > 25
    res['candidates'] = [gauntlet(S, e, tg['depth'], eps, rng) for e in ev_sorted[:25]]
    res['comparison_null'] = comp_null(comp_refs)
    if tg['kind'] == 'control_period':
        res['period_scan'] = period_control(S)
    return res


def write_md(res):
    L = [f"# {res['name']} — WAVE2 DASCH search", "",
         f"Kind: {res['kind']}; status: {res.get('status')}"]
    if res.get('status') != 'searched':
        L += [f"Error: {res.get('error', 'see JSON')}", ""]
        return "\n".join(L)
    s = res['sensitivity']
    L += [f"DASCH {res['dasch_ref']} (sep {res['dasch_sep_arcsec']} arcsec, B={res['dasch_B']}, "
          f"{res['n_comparisons']} comparisons). Median corrected B {res['median_B_corrected']}, "
          f"plate sigma {res['plate_sigma_mag']} mag, {res['n_good_plates']} good plates "
          f"({res['n_plate_rows']} rows). Faint threshold {s['faint_threshold_mag']} mag; "
          f"catalogue depth {s['catalogue_depth_mag']} mag "
          f"({'detectable' if s['catalogue_depth_detectable'] else 'BELOW threshold'}).",
          f"Caught fraction (1890-1990, >=1 normal plate in window): {s['caught_frac']}; "
          f">=3 good nights in window: {s['caught_frac_3nights']}.",
          f"Blind events >=3 nights: {res['n_events_ge3nights']} "
          f"(>=1-night faint clusters: {res['n_clusters_ge1night']}, context only). "
          f"Comparison null: {res['comparison_null']['with_events']}/"
          f"{res['comparison_null']['comparisons_tested']} comparisons show "
          f"{res['comparison_null']['total_events']} event(s).", ""]
    for c in res['candidates']:
        n = c['shifted_epoch_null']
        L += [f"## Candidate {c['date_start']} to {c['date_end']} "
              f"(JD {c['jd_start']}-{c['jd_end']}, span {c['span_d']} d)",
              f"{c['n_plates']} plates, {c['n_nights']} nights, series {','.join(c['series'])}; "
              f"fit depth {c['depth_fit_mag']} mag; dlogL box vs no-event "
              f"{c['dlogL_box_fitdepth']} (catalogue depth: {c['dlogL_box_catdepth']}); "
              f"shifted-epoch null (n=200): median {n['median']}, p95 {n['p95']}, max {n['max']}, "
              f"fraction >= candidate {n['frac_ge_candidate']}.", "",
              "| plate | JD | date | kind | dm | lim | comps det |",
              "|---|---|---|---|---|---|---|"]
        L += [f"| {p['plate']} | {p['jd']} | {p['date']} | {p['kind']} | {p['dm']} | {p['lim']} "
              f"| {p['comparisons_detected']} |" for p in c['plates']]
        L += [""]
    if 'period_scan' in res:
        p = res['period_scan']
        L += [f"Period scan (control): peaks (P d, dlogL, nights) {p['peaks']}; "
              f"null40 median {p['null40']['median']}, max {p['null40']['max']}, "
              f"frac>=best {p['null40']['frac_ge_best']}. Known: {p['known']}", ""]
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', nargs='*', help='substring filter on star name')
    ap.add_argument('--n-tier-a', type=int, default=8)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    stars = CONTROLS + tier_a_targets(args.n_tier_a)
    if args.only:
        stars = [s for s in stars if any(o.lower() in s['name'].lower() for o in args.only)]
    rng = np.random.default_rng(20261018)
    summary = []
    for tg in stars:
        print(f"=== {tg['name']} ===", flush=True)
        res = run_star(tg, rng)
        base = os.path.join(OUT, slug(tg['name']))
        json.dump(res, open(base + '.json', 'w'), indent=1)
        open(base + '.md', 'w').write(write_md(res))
        summary.append({k: res.get(k) for k in
                        ('name', 'kind', 'status', 'error', 'dasch_ref', 'n_good_plates',
                         'plate_sigma_mag', 'n_events_ge3nights')})
        print(json.dumps(summary[-1]), flush=True)
    spath = os.path.join(OUT, 'summary.json')   # merge so --only runs compose
    if os.path.exists(spath):
        old = {s['name']: s for s in json.load(open(spath))}
        old.update({s['name']: s for s in summary})
        summary = list(old.values())
    json.dump(summary, open(spath, 'w'), indent=1)
    print('done', len(summary))


if __name__ == '__main__':
    main()
