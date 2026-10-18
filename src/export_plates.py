"""Export the per-plate evidence for ASASSN-21co (J181752-580749): every plate row with the comparison-star corrected
brightness (or the plate's limit when the star was not detected), how many comparison stars were detected, whether
the row was used in the eclipse fit (good detections only) or is an informative non-detection (plate reached 0.5 mag
below normal and >= 40% of comparison stars were seen), and the eclipse model's prediction at the fitted ephemeris. Writes results/asassn21co_plates.csv.
DASCH plate IDs let anyone re-fetch the scans and cutouts."""
import csv, json, os
import numpy as np
from dasch import DATA
from core import plates, year
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
xm = json.load(open(os.path.join(DATA, 'xmatch.json'))); comps = json.load(open(os.path.join(DATA, 'comps.json')))
NAME = 'J181752-580749'
S = plates(xm[NAME]['target']['ref_number'], [c['ref_number'] for c in comps[NAME]])
eph = json.load(open(os.path.join(OUT, 'asassn21co.json')))['ephemeris']
T0, P = eph['T0_jd'], eph['P_d']; tT, tF, depth = 66.0, 23.0, 0.75
def model(dt):
    a = abs(dt)
    return depth if a <= tF / 2 else (depth * (tT / 2 - a) / ((tT - tF) / 2) if a < tT / 2 else 0.0)
with open(os.path.join(OUT, 'asassn21co_plates.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['plate', 'series', 'jd', 'year', 'kind', 'used_in_eclipse_fit', 'informative_nondetection', 'dB_mag', 'limit_minus_normal_mag', 'comparisons_detected', 'comparisons_on_plate',
                'days_from_predicted_mid_eclipse', 'model_dB_mag'])
    for r in S['plates']:
        dt = ((r['jd'] - T0) / P + 0.5) % 1 * P - P / 2
        w.writerow([r['plate'], r['series'], f"{r['jd']:.5f}", f"{year(r['jd']):.4f}", 'detection' if r['kind'] == 'D' else 'non-detection',
                    int(r['kind'] == 'D' and r['ok']), int(r['kind'] == 'N' and r['lim'] > 0.5 and r['cn'] >= 3 and r['cdet'] >= 0.4 * r['cn']),
                    '' if r['dm'] is None else f"{r['dm']:.3f}", f"{r['lim']:.3f}", r['cdet'], r['cn'], f"{dt:.1f}", f"{model(dt):.2f}"])
print('normal level', round(S['med'], 3), 'B; plate sigma', round(S['sig'], 3), '; rows', len(S['plates']))
