"""Step 4: ASASSN-21co (J181752-580749). Historical eclipse times from DASCH plates, period from plates alone,
checks against the 2009 (ASAS) and 2021 (ASAS-SN) eclipses, a linear ephemeris over all four, and a figure.
Eclipse shape fixed to Rowan et al. 2021 (RNAAS 5, 147): total duration 66 d, totality 23 d."""
import datetime as dtm, json, math, os
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from dasch import DATA
from core import plates, year
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
xm = json.load(open(os.path.join(DATA, 'xmatch.json'))); comps = json.load(open(os.path.join(DATA, 'comps.json')))
NAME = 'J181752-580749'
S = plates(xm[NAME]['target']['ref_number'], [c['ref_number'] for c in comps[NAME]])
sig = math.hypot(S['sig'], 0.05)
good = [r for r in S['plates'] if r['kind'] == 'D' and r['ok']]
jd = np.array([r['jd'] for r in good]); dm = np.array([r['dm'] for r in good])
tT, tF = 66.0, 23.0

def model(dt, depth):
    a = np.abs(dt); m = np.where(a <= tF / 2, depth, 0.0)
    ing = (a > tF / 2) & (a < tT / 2); m[ing] = depth * (tT / 2 - a[ing]) / ((tT - tF) / 2)
    return m
def ll(r): return -0.5 * np.minimum((r / sig) ** 2, 16)          # Gaussian, clipped at 4 sigma
def jd2date(j): return dtm.date.fromordinal(int(j - 1721424.5)).isoformat()

def fit_epoch(guess, halfwin=120):
    """Profile likelihood of mid-eclipse time (depth free) from plates within +-halfwin d of a guess."""
    sel = np.abs(jd - guess) < halfwin; prof = []
    for Tm in np.arange(guess - 60, guess + 60, 0.25):
        prof.append((Tm, max(ll(dm[sel] - model(jd[sel] - Tm, d)).sum() for d in np.arange(0.3, 1.21, 0.05))))
    prof = np.array(prof); b = prof[np.argmax(prof[:, 1])]; ok = prof[:, 1] > b[1] - 0.5
    return dict(T_best=float(b[0]), T_lo=float(prof[ok, 0].min()), T_hi=float(prof[ok, 0].max()),
                dlogL_vs_no_eclipse=round(float(b[1] - ll(dm[sel]).sum()), 1), plates_in_window=int(sel.sum()),
                plates_in_eclipse=[dict(plate=good[i]['plate'], jd=round(float(jd[i]), 2), date=jd2date(jd[i]), dB=round(float(dm[i]), 2),
                                        comparisons_detected=f"{good[i]['cdet']}/{good[i]['cn']}")
                                   for i in np.where(sel & (np.abs(jd - b[0]) < tT / 2))[0]])

T21 = (2459255.5 + 2459257.5) / 2 + tT / 2    # ASAS-SN: eclipse began 2021 Feb 10-12 (AAVSO Alert Notice 732)
E85 = fit_epoch(T21 - 3 * 4344); E02 = fit_epoch(T21 - 10 * 4344)   # windows located with the published period
mid = lambda E: (E['T_lo'] + E['T_hi']) / 2
err = lambda E: max(2.0, (E['T_hi'] - E['T_lo']) / 2)
# period from the plates alone, and what it predicts for the two modern eclipses
Pp = (mid(E85) - mid(E02)) / 7; sPp = math.hypot(err(E85), err(E02)) / 7
pred = {y: (mid(E85) + c * Pp, math.hypot(err(E85), c * sPp)) for y, c in ((2009, 2), (2021, 3))}
# linear ephemeris over all four eclipses (modern timing uncertainties are assumed, see 'source')
ep = [(-10, mid(E02), err(E02), 'DASCH plate b29521 (this work)'),
      (-3, mid(E85), err(E85), 'DASCH plates dsb1716, dsb1725 (this work)'),
      (-1, 2454944.5, 5.0, 'ASAS, "around UT 2009-04-23" (Rowan et al. 2021); 5 d sigma assumed'),
      (0, T21, 3.0, 'ASAS-SN, began 2021 Feb 10-12 (AAVSO Alert Notice 732) + 33 d; 3 d sigma assumed')]
E = np.array([e[0] for e in ep], float); Tm = np.array([e[1] for e in ep]); s = np.array([e[2] for e in ep])
A = np.vstack([np.ones_like(E), E]).T; W = np.diag(1 / s**2); cov = np.linalg.inv(A.T @ W @ A); T0, P = cov @ A.T @ W @ Tm
chi2 = float(((Tm - A @ [T0, P]) / s) @ ((Tm - A @ [T0, P]) / s)); Tn = T0 + P; sTn = float(np.sqrt([1, 1] @ cov @ [1, 1]))
res = dict(star='ASASSN-21co = ASASSN-V J181752.04-580749.3', dasch_source=xm[NAME]['target']['ref_text'],
           good_plates=len(good), plate_sigma_mag=round(sig, 3), comparison_stars=len(comps[NAME]),
           eclipse_1985=E85, eclipse_1902=E02,
           period_from_plates_only=dict(P_d=round(Pp, 1), sigma_d=round(sPp, 1)),
           plates_only_predicts={str(y): dict(mid_jd=round(t, 1), date=jd2date(t), sigma_d=round(u, 1)) for y, (t, u) in pred.items()},
           observed={'2009': 'around UT 2009-04-23 (ASAS; Rowan et al. 2021)', '2021': f'began 2021 Feb 10-12 (ASAS-SN), mid about {jd2date(T21)}'},
           ephemeris=dict(epochs=[dict(cycle=int(c), mid_jd=round(t, 1), sigma_d=round(u, 1), source=src) for c, t, u, src in ep],
                          T0_jd=round(float(T0), 1), T0_date=jd2date(T0), P_d=round(float(P), 1), sigma_P_d=round(float(np.sqrt(cov[1, 1])), 1),
                          chi2=round(chi2, 2), dof=2,
                          next_primary_eclipse=dict(mid_jd=round(float(Tn), 1), mid_date=jd2date(Tn), sigma_d=round(sTn, 1),
                                                    first_contact=jd2date(Tn - tT / 2), last_contact=jd2date(Tn + tT / 2))))
os.makedirs(OUT, exist_ok=True); json.dump(res, open(os.path.join(OUT, 'asassn21co.json'), 'w'), indent=1)
print(json.dumps({k: res[k] for k in ('period_from_plates_only', 'plates_only_predicts', 'observed')}, indent=1))
print(json.dumps({k: v for k, v in res['ephemeris'].items() if k != 'epochs'}, indent=1))

fig, ax = plt.subplots(2, 1, figsize=(9, 6.4))
yrs = np.array([year(j) for j in jd])
ax[0].plot(yrs, dm, '.', ms=3, color='0.35'); ax[0].invert_yaxis()
for x in (mid(E02), mid(E85)): ax[0].axvline(year(x), color='crimson', lw=0.8, alpha=0.6)
ax[0].set_ylabel('B − median (mag)')
ax[0].set_title(f'ASASSN-21co on Harvard plates (DASCH DR7): {len(good)} good plates, comparison-star corrected', fontsize=10)
dtd = (((jd - T0) / P + 0.5) % 1 - 0.5) * P; s150 = np.abs(dtd) < 150
ax[1].plot(dtd[s150], dm[s150], 'o', ms=4, color='0.25', label='plates, 1889–1990')
xx = np.linspace(-150, 150, 601); ax[1].plot(xx, model(xx, 0.75), color='crimson', lw=1, label='eclipse shape of Rowan et al. 2021, depth 0.75')
for i in np.where(s150 & (dm > 0.4))[0]:
    ax[1].annotate(jd2date(jd[i]), (dtd[i], dm[i]), fontsize=7, xytext=(4, -2), textcoords='offset points')
ax[1].invert_yaxis(); ax[1].set_xlabel(f'days from mid-eclipse, folded on P = {P:.1f} d'); ax[1].set_ylabel('B − median (mag)'); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(OUT, 'asassn21co.png'), dpi=120)
