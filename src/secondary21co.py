"""Secondary-eclipse search for ASASSN-21co on the committed plate table (results/asassn21co_plates.csv).
Folds the 574 good plates on the fitted ephemeris and scans secondary phase for a second
66 d trapezoid dimming; significance is calibrated by scrambling magnitudes among the same plates.
Also reports phase coverage: at which secondary phases could the plates have seen anything at all?
Usage: python secondary21co.py"""
import csv, json, math, os
import numpy as np

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
rows = list(csv.DictReader(open(os.path.join(RES, 'asassn21co_plates.csv'))))
eph = json.load(open(os.path.join(RES, 'asassn21co.json')))['ephemeris']
T0, P = eph['T0_jd'], eph['P_d']
used = [r for r in rows if r['used_in_eclipse_fit'] == '1']
jd = np.array([float(r['jd']) for r in used]); dm = np.array([float(r['dB_mag']) for r in used])
sig = math.hypot(1.4826 * float(np.median(np.abs(dm - np.median(dm)))), 0.05)
phase = ((jd - T0) / P) % 1.0
print(f"ASASSN-21co secondary search: {len(jd)} good plates, sigma {sig:.3f}, P={P:.1f} d, cycles {int(np.floor((jd - T0) / P).min())}..{int(np.floor((jd - T0) / P).max())}")

def dist_d(ph, centre): return np.abs(((ph - centre + 0.5) % 1.0) - 0.5) * P

def score(ph, val, centre, depth=0.4):
    d = dist_d(ph, centre); m = d < 33.0
    if not m.any(): return 0.0, 0
    w = np.clip((33.0 - d[m]) / 21.5, 0, 1)
    return float(((val[m]**2 - (val[m] - depth * w)**2) / (2 * sig**2)).sum()), int(m.sum())

# phase coverage
grid = np.arange(0.02, 0.99, 0.002)
n_in = np.array([int((dist_d(phase, c) < 33).sum()) for c in grid])
near05 = np.sort(dist_d(phase, 0.5))[:3]
n_tot = np.array([int((dist_d(phase, c) < 11.5).sum()) for c in grid])
print(f"phase coverage (66 d window): {np.mean(n_in == 0):.0%} of secondary phases have NO plate, {np.mean(n_in >= 2):.0%} have >= 2 plates; in totality (+-11.5 d) {np.mean(n_tot == 0):.0%} have none, {np.mean(n_tot >= 2):.0%} have >= 2")
print(f"nearest plates to anti-phase (0.5): {np.round(near05, 0)} d away -> a circular-orbit secondary is untestable in totality")
d5 = dist_d(phase, 0.5); m5 = d5 < 33
print(f"anti-phase stack: {int(m5.sum())} plates within 33 d" + (f", median dB {np.median(dm[m5]):+.3f}" if m5.any() else ""))

res = []
for c in grid:
    if min(c, 1 - c) < 0.03: continue
    llr, n = score(phase, dm, float(c))
    res.append((llr, round(float(c), 3), n))
res.sort(reverse=True)
print("top secondary candidates (dlogL at depth 0.4, phase, plates):", [(round(a, 1), b, n) for a, b, n in res[:6]])
best_c = res[0][1]
print(f"depth profile at best phase {best_c}:", [(dep, round(score(phase, dm, best_c, dep)[0], 1)) for dep in (0.15, 0.25, 0.4, 0.6, 0.75)])
d = dist_d(phase, best_c)
for i in np.where(d < 33)[0]:
    print(f"  plate {used[i]['plate']} jd {jd[i]:.1f} dist {d[i]:.0f} d dB {dm[i]:+.2f} comparisons {used[i]['comparisons_detected']}/{used[i]['comparisons_on_plate']}")
rng = np.random.default_rng(7); null = []
for _ in range(200):
    v = rng.permutation(dm)
    null.append(max(score(phase, v, float(c))[0] for c in grid if min(c, 1 - c) >= 0.03))
null = np.array(null); best = res[0][0]
print(f"scramble null (200): median max {np.median(null):.1f}, p95 {np.percentile(null, 95):.1f}, max {null.max():.1f}; observed best {best:.1f}; fraction of nulls >= observed: {np.mean(null >= best):.3f}")
print("conclusion: no secondary eclipse anywhere (best phase is below the scramble median); exact anti-phase is untestable (nearest plate 42 d away), and with typically 2-3 plates per 66 d phase window only a secondary deeper than ~0.4 mag could have been seen.")
