"""Fetch DASCH plate cutouts around a position and render a strip of PNG thumbnails, star marked."""
import base64, gzip, io, json, os, sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from dasch import call, DATA, rows
ra, dec = float(sys.argv[1]), float(sys.argv[2]); ref = sys.argv[3]; out = sys.argv[4]; want = sys.argv[5].split(',')
labels = dict(x.split('=') for x in sys.argv[6].split(';')) if len(sys.argv) > 6 else {}
lc = rows(json.load(open(os.path.join(DATA, 'lc', f'{ref}.json'))))
sol = {}
for x in lc:
    pid = f"{x['series']}{x['plate_number']}"
    if pid in want: sol[pid] = (x['series'], int(x['plate_number']), int(x['solution_number'] or 0))
fig, axs = plt.subplots(1, len(want), figsize=(3.2 * len(want), 3.6))
for ax, pid in zip(np.atleast_1d(axs), want):
    s, n, sn = sol[pid]
    f = os.path.join(DATA, 'cutouts', f'{s}{n:05d}_{sn}_{ra:.4f}_{dec:.4f}.fits')
    if not os.path.exists(f):
        r = call('cutout', {"plate_id": f"{s}{n:05d}", "solution_number": sn, "center_ra_deg": ra, "center_dec_deg": dec})
        open(f, 'wb').write(gzip.decompress(base64.b64decode(r)))
    h = fits.open(f)[0]; w = WCS(h.header); img = h.data.astype(float)
    x, y = w.world_to_pixel_values(ra, dec)
    # crop ~ +-90 arcsec
    scale = abs(w.proj_plane_pixel_scales()[0].value) * 3600
    r = int(90 / scale); xi, yi = int(round(float(x))), int(round(float(y)))
    sub = img[max(0, yi - r):yi + r, max(0, xi - r):xi + r]
    lo, hi = np.percentile(sub, [1, 99.5])
    ax.imshow(sub, origin='lower', cmap='gray_r', vmin=lo, vmax=hi)
    cx, cy = xi - max(0, xi - r), yi - max(0, yi - r)
    ax.add_patch(plt.Circle((cx, cy), 12 / scale, fill=False, color='crimson', lw=1.2))
    ax.set_title(f"{pid}\n{labels.get(pid, '')}", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
plt.tight_layout(); plt.savefig(out, dpi=110); print('wrote', out, 'pixel scale', round(scale, 2), 'arcsec')
