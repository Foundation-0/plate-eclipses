"""Independent re-measurement on DASCH cutouts: background-subtracted plate density in a small aperture,
calibrated per plate against APASS stars in the same cutout (target excluded), with a quadratic fit of magnitude against log summed density."""
import json, math, os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from dasch import DATA, querycat
ra, dec, ref_text = float(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
files = sys.argv[4:]
cat = [r for r in querycat(ra, dec, 900) if r['stdmag'] != '99.9']
def phot(img, x, y, rap, rin, rout):
    yy, xx = np.indices(img.shape); d = np.hypot(xx - x, yy - y)
    ann = img[(d > rin) & (d < rout)]; bg = np.median(ann)
    ap = img[d < rap]
    return float(np.sum(ap - bg)), float(np.std(ann))
for f in files:
    h = fits.open(f)[0]; w = WCS(h.header); img = h.data.astype(float)
    scale = abs(w.proj_plane_pixel_scales()[0].value) * 3600
    rap, rin, rout = 20 / scale, 35 / scale, 55 / scale
    ny, nx = img.shape
    pts = []; tgt = None
    for r in cat:
        x, y = w.world_to_pixel_values(float(r['ra_deg']), float(r['dec_deg']))
        x, y = float(x), float(y)
        if not (rout < x < nx - rout and rout < y < ny - rout): continue
        # isolation: no catalogued star brighter than mag+2 within 40"
        flux, _ = phot(img, x, y, rap, rin, rout)
        if r['ref_text'] == ref_text: tgt = (flux, float(r['stdmag'])); continue
        m = float(r['stdmag'])
        near = [s for s in cat if s is not r and math.hypot((float(s['ra_deg']) - float(r['ra_deg'])) * math.cos(math.radians(dec)), float(s['dec_deg']) - float(r['dec_deg'])) * 3600 < 40 and float(s['stdmag']) < m + 2]
        if near or flux <= 0: continue
        pts.append((m, math.log10(flux)))
    pts = np.array(pts)
    sel = (pts[:, 0] > tgt[1] - 1.5) & (pts[:, 0] < tgt[1] + 1.5)
    p = pts[sel]
    c = np.polyfit(p[:, 1], p[:, 0], 2)  # mag as a function of log density sum
    res = p[:, 0] - np.polyval(c, p[:, 1]); s = 1.4826 * np.median(np.abs(res - np.median(res)))
    mt = np.polyval(c, math.log10(max(tgt[0], 1e-3)))
    print(f"{os.path.basename(f)[:20]:20s} calib stars {len(p):3d} scatter {s:.2f} | target B_cat {tgt[1]:.2f} measured {mt:.2f} (dm {mt - tgt[1]:+.2f})")
