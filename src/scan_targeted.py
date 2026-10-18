"""Targeted period scan (second-generation), anchored on a modern dip.

Differences from period.py (kept for the record): plate zero points must rest on >= 4 comparison
residuals (cdet >= 4); only informative non-detections are scored, and a deep non-detection
inside a predicted eclipse counts AGAINST it (an eclipsed star that deep should have been seen);
detection log-likelihood ratios use a trapezoid matched filter capped at +-9; a peak needs
support from >= 2 nights in >= 2 different cycles. Null calibration shifts the anchor epoch.

Usage: python scan_targeted.py NAME [--pmin D] [--pmax D] [--dur D] [--depth M] [--nulls N]
Anchor, duration, depth and minimum period default to the targets.py row (mid-dip = dip
start + duration/2)."""
import datetime as dtm, json, math, os, sys
import numpy as np
from dasch import DATA
from core import plates
from targets import T

def jd_of(date_str, add=0.0):
    d, m, y = (int(x) for x in date_str.split('/'))
    return dtm.date(y, m, d).toordinal() + 1721424.5 + add

def prepare(name):
    xm = json.load(open(os.path.join(DATA, 'xmatch.json')))
    comps = json.load(open(os.path.join(DATA, 'comps.json')))
    S = plates(xm[name]['target']['ref_number'], [c['ref_number'] for c in comps[name]])
    det = [r for r in S['plates'] if r['kind'] == 'D' and r['ok'] and r['cdet'] >= 4]
    dm = np.array([r['dm'] for r in det])
    sig = math.hypot(1.4826 * float(np.median(np.abs(dm - np.median(dm)))), 0.05)
    nond = [r for r in S['plates'] if r['kind'] == 'N' and r['lim'] > 0.5 and r['cn'] >= 4 and r['cdet'] >= max(3, 0.4 * r['cn'])]
    obs = [(r['jd'], True, r['dm'], r['plate']) for r in det] + [(r['jd'], False, r['lim'], r['plate']) for r in nond]
    return S, sig, obs

def scan(obs, sig, T0, D, depth, Pmin, Pmax):
    jd = np.array([o[0] for o in obs]); isdet = np.array([o[1] for o in obs]); val = np.array([o[2] for o in obs])
    span = T0 - jd.min(); Ps, L, NN, CC = [], [], [], []
    P = float(Pmin)
    while P < Pmax:
        ph = (T0 - jd) / P; k = np.round(ph); dist = np.abs((ph - k) * P)
        m = (dist < D / 2) & (k >= 1)
        w = np.clip((D / 2 - dist[m]) / (D / 4), 0, 1)
        v, d = val[m], isdet[m]
        llr_det = np.clip((v**2 - (v - depth * w)**2) / (2 * sig**2), -9, 9)
        llr_nd = np.where(v <= depth - 0.2, 3.0, np.where(v >= depth + 0.2, -3.0, 0.0))
        llr = np.where(d, llr_det, llr_nd)
        pos = llr > 1
        Ps.append(P); L.append(float(llr.sum()))
        NN.append(len(set(np.floor(jd[m][pos])))); CC.append(len(set(k[m][pos].astype(int))))
        P += (D / 4) / max(1.0, span / P) / 3
    return np.array(Ps), np.array(L), np.array(NN), np.array(CC)

def peaks(Ps, L, NN, CC, n=6):
    out = []
    for i in np.argsort(-L):
        if L[i] < 3 or len(out) >= n: break
        if NN[i] < 2 or CC[i] < 2 or any(abs(Ps[i] - p) < 0.01 * p for p, *_ in out): continue
        out.append((round(float(Ps[i]), 1), round(float(L[i]), 1), int(NN[i]), int(CC[i])))
    return out

def contributors(obs, sig, T0, D, depth, P):
    rows = []
    for jd0, isdet, v, pid in obs:
        ph = (T0 - jd0) / P; k = round(ph); dist = abs((ph - k) * P)
        if dist < D / 2 and k >= 1:
            w = min(1.0, max(0.0, (D / 2 - dist) / (D / 4)))
            llr = float(np.clip((v**2 - (v - depth * w)**2) / (2 * sig**2), -9, 9)) if isdet else (3.0 if v <= depth - 0.2 else (-3.0 if v >= depth + 0.2 else 0.0))
            rows.append((pid, round(jd0, 1), k, round(dist, 1), 'D' if isdet else 'N', round(v, 2), round(llr, 1)))
    return sorted(rows, key=lambda r: -r[6])

if __name__ == '__main__':
    args = sys.argv[1:]; name = args[0]
    opt = dict(zip(args[1::2], args[2::2]))
    row = next(t for t in T if t[0] == name)
    _, _, _, cls, g, ndip, dur, pmin, last, depth = row
    D = float(opt.get('--dur', max(dur, 5.0))); depth = float(opt.get('--depth', depth))
    Pmin = float(opt.get('--pmin', pmin if pmin else 500)); Pmax = float(opt.get('--pmax', 45000))
    nnull = int(opt.get('--nulls', 30))
    T0 = jd_of(last, D / 2)
    S, sig, obs = prepare(name)
    print(f"{name} ({cls}): {len(obs)} scored plates ({sum(o[1] for o in obs)} det / {sum(not o[1] for o in obs)} informative non-det), sigma {sig:.3f}, ngood(core) {S['ngood']}")
    print(f"  anchor mid-dip JD {T0:.1f} (dip start {last} + {D/2:.1f} d), D={D:.1f} d, depth={depth:.2f}, P range {Pmin:.0f}-{Pmax:.0f} d")
    Ps, L, NN, CC = scan(obs, sig, T0, D, depth, Pmin, Pmax)
    pk = peaks(Ps, L, NN, CC)
    print(f"  {len(Ps)} trial periods; excluded (dlogL<-3): {np.mean(L < -3):.0%}; peaks (P, dlogL, nights, cycles): {pk}")
    for off in (-D / 4, D / 4):
        P2, L2, N2, C2 = scan(obs, sig, T0 + off, D, depth, Pmin, Pmax)
        print(f"  anchor offset {off:+.1f} d: peaks {peaks(P2, L2, N2, C2, 3)}")
    if pk:
        print(f"  contributors at best P={pk[0][0]} (plate, jd, cycle k, dist d, kind, value, llr):")
        for r in contributors(obs, sig, T0, D, depth, pk[0][0])[:14]: print("   ", r)
    if nnull:
        rng = np.random.default_rng(3); null = []
        for _ in range(nnull):
            off = rng.uniform(3 * D, 3 * D + 0.5 * max(Pmin, 600)) * rng.choice([-1, 1])
            _, L2, N2, C2 = scan(obs, sig, T0 + off, D, depth, Pmin, Pmax)
            m = (N2 >= 2) & (C2 >= 2); null.append(float(L2[m].max()) if m.any() else 0.0)
        null = np.array(null); best = pk[0][1] if pk else 0.0
        print(f"  nulls ({nnull} shifted anchors): median best {np.median(null):.1f}, max {null.max():.1f}; fraction >= observed best ({best}): {np.mean(null >= best):.2f}")
