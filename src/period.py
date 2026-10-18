"""Blind period scan for a binary with one well-timed modern eclipse (T0, full duration D, depth).
For each trial period P the past eclipses are at T0 - kP; plates inside the eclipse core are scored
against 'in eclipse' vs 'out of eclipse' with comparison-star-corrected magnitudes."""
import json, math, os, sys
import numpy as np
from dasch import DATA
from core import plates

def prepare(name):
    xm = json.load(open(os.path.join(DATA, 'xmatch.json'))); comps = json.load(open(os.path.join(DATA, 'comps.json')))
    S = plates(xm[name]['target']['ref_number'], [c['ref_number'] for c in comps[name]])
    ok = [r for r in S['plates'] if r['cn'] >= 3 and r['cdet'] >= 0.4 * r['cn']]
    deep = [r for r in ok if (r['kind'] == 'D' and r['ok']) or (r['kind'] == 'N' and r['lim'] > 0.5)]
    eps = max(0.02, sum(r['kind'] == 'N' for r in deep) / max(1, len(deep)))   # empirical miss rate
    obs = [(r['jd'], 'D', r['dm']) for r in S['plates'] if r['kind'] == 'D' and r['ok']]
    obs += [(r['jd'], 'N', r['lim']) for r in ok if r['kind'] == 'N']
    return S, eps, obs

def scan(obs, sig, eps, T0, D, depth, Pmin, Pmax=45000, core=0.6):
    jd = np.array([o[0] for o in obs]); det = np.array([o[1] == 'D' for o in obs]); val = np.array([o[2] for o in obs])
    s = math.hypot(sig, 0.05)
    llr = np.where(det, np.clip((val**2 - (val - depth)**2) / (2 * s**2), -6, 6), np.where(val < depth - 0.2, math.log(1 / eps), 0.0))
    h = core * D / 2; span = T0 - jd.min()
    Ps, L, NN = [], [], []
    P = Pmin
    while P < Pmax:
        ph = (T0 - jd) / P; k = np.round(ph); m = (np.abs((ph - k) * P) < h) & (k >= 1)
        Ps.append(P); L.append(llr[m].sum()); NN.append(len(set(np.floor(jd[m & (llr > 1)]))))
        P += h / max(1, span / P) / 3
    return np.array(Ps), np.array(L), np.array(NN)

def peaks(Ps, L, NN, k=6, min_nights=2):
    out = []
    for i in np.argsort(-L):
        if L[i] < 3 or len(out) >= k: break
        if NN[i] < min_nights or any(abs(Ps[i] - p) < 0.01 * p for p, _, _ in out): continue
        out.append((round(float(Ps[i]), 1), round(float(L[i]), 1), int(NN[i])))
    return out

if __name__ == '__main__':
    name, T0, D, depth, Pmin = sys.argv[1], *map(float, sys.argv[2:6])
    ntrial = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    S, eps, obs = prepare(name)
    Ps, L, NN = scan(obs, S['sig'], eps, T0, D, depth, Pmin)
    print(f"{name}: sigma {S['sig']:.2f}, miss rate {eps:.3f}, {len(obs)} plates, {len(Ps)} trial periods {Pmin:.0f}-45000 d")
    print(f"  excluded (dlogL < -3): {np.mean(L < -3):.0%}; no plate in any predicted eclipse core: {np.mean(L == 0):.0%}")
    print('  peaks (P d, dlogL, nights):', peaks(Ps, L, NN))
    if ntrial:   # false-alarm calibration: same scan with the modern eclipse moved off its true phase
        rng = np.random.default_rng(1); null = []
        for _ in range(ntrial):
            off = rng.uniform(3 * D, 0.5 * Pmin) * rng.choice([-1, 1])
            Ps2, L2, N2 = scan(obs, S['sig'], eps, T0 + off, D, depth, Pmin)
            m = N2 >= 2; null.append(float(L2[m].max()) if m.any() else 0.0)
        null = np.array(null); best = peaks(Ps, L, NN)[0][1] if peaks(Ps, L, NN) else 0
        print(f"  null best-peak dlogL over {ntrial} shifted epochs: median {np.median(null):.1f}, max {null.max():.1f}; fraction >= observed best ({best}): {np.mean(null >= best):.2f}")
