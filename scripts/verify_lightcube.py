"""Standalone verification (cloudfield/scripts). The Zucker et al. (2021) Table 1
extents below are used ONLY for the Ophiuchus/Taurus check against Bish et al.;
they are not cloud boundaries for the pipeline.

Standalone verification of the LightCube sampling / units / density
conversion against published numbers in Bish et al. (2026, arXiv:2608.27389).
Changes nothing; imports cloudfield read-only.

Check 1  solar value: mean I_UV within 25 pc of the grid origin   (Sec III.1: 63.3%)
Check 2  vertical profile: exponential scale height                (Sec III.3: ~400 pc)
Check 3  Ophiuchus / Taurus mean ISRF over n > 30 cm^-3 gas       (Sec IV.3: ~3x)
"""
import sys
import numpy as np
from cloudfield.cubes import open_cubes, origin_index
from cloudfield.conversions import NH_CONVERSIONS, DEFAULT_CONVERSION as NH_CONVERSION, nh_per_tau
from cloudfield.paths import data_path

ORIGIN_MODE = "header"
R_CLIP = 500.0               # pc; our files are clipped at r = 500 pc

# Check 1
SUN_R = 25.0
SUN_EXPECTED = 0.633         # Bish+2026 Sec III.1: "63.3% of the Standard ISRF"
SUN_TOL = 0.01               # fractional

# Check 2
CYL_R = 140.0                # pc; fixed horizontal footprint -> |z| up to sqrt(500^2-140^2)
ZBIN = 10.0                  # pc
H_EXPECTED = 400.0           # pc, Bish+2026 Sec III.3

# Check 3: cloud volumes. Bish+2026 give NO boundaries for the Sec IV.3 / Fig 10
# comparison ("molecular clouds within 400 pc ... averaging over the entire
# n>30 cm-3 cloud volume"). The published per-cloud 3D extents available here
# are Zucker et al. 2021 (ApJ 919, 35) Table 1: min/max heliocentric Galactic
# Cartesian x, y, z of each cloud's 3D SKELETON (Leike+2020 map, n=35 cm^-3
# features). These are skeleton extents, not envelope boundaries; no padding
# is added. The proper boundaries are the Zucker+21 3D masks
# (doi:10.7910/DVN/NMRLP1), not available locally.
CLOUDS = {   # name: (x_min, x_max, y_min, y_max, z_min, z_max) pc, Zucker+21 Table 1
    "Ophiuchus": (104, 144, -26, 5, 24, 48),
    "Taurus": (-166, -122, 6, 44, -45, -26),
}
N_CUT = 30.0                 # cm^-3
RATIO_EXPECTED = 3.0
N_CUT_ALT = (20.0, 40.0)


def hdr(t):
    print("=" * 78)
    print(t)


isrf, dust, h = open_cubes(data_path("lightcube"), data_path("edenhofer_tau"))
o = origin_index(isrf.shape, h, ORIGIN_MODE)
print(f"cube {isrf.shape}, origin [z,y,x] = {o} ({ORIGIN_MODE}); conversion '{NH_CONVERSION}' "
      f"= {nh_per_tau():.2f} cm^-3 per tau_1565/pc")
if NH_CONVERSION != "lightcube":
    sys.exit("Check 3 requires the 'lightcube' conversion (Bish+2026 Fig 7 density axis).")
RESULT = {}

# ----------------------------------------------------------------------------
hdr("CHECK 1  solar value")
k = int(np.ceil(SUN_R)) + 1
c = np.floor(o).astype(int)
sl = tuple(slice(ci - k, ci + k + 2) for ci in c)
blk = np.asarray(isrf[sl], dtype=np.float64)
zz, yy, xx = [np.arange(s.start, s.stop) - oi for s, oi in zip(sl, o)]
R = np.sqrt(zz[:, None, None] ** 2 + yy[None, :, None] ** 2 + xx[None, None, :] ** 2)
v = blk[R <= SUN_R]
m1 = v.mean()
# also the integer-centre convention used in the first run
o2 = origin_index(isrf.shape, h, "n//2")
k2 = int(SUN_R)
b2 = np.asarray(isrf[tuple(slice(int(oi) - k2, int(oi) + k2 + 1) for oi in o2)], dtype=np.float64)
g = np.arange(-k2, k2 + 1)
r2 = np.sqrt(g[:, None, None] ** 2 + g[None, :, None] ** 2 + g[None, None, :] ** 2)
m2 = b2[r2 <= SUN_R].mean()
print(f"  mean I_UV, voxels with centre within {SUN_R:g} pc of origin:")
print(f"    header origin (499.5): {m1:.4f}  ({v.size} voxels, median {np.median(v):.4f})")
print(f"    n//2 origin (500):     {m2:.4f}  (first-run convention)")
dev = m1 / SUN_EXPECTED - 1
RESULT[1] = abs(dev) <= SUN_TOL
print(f"  expected {SUN_EXPECTED} -> deviation {dev:+.2%} (tolerance {SUN_TOL:.0%}): "
      f"{'PASS' if RESULT[1] else 'FAIL'}")
print("  (Check 1 does not involve the density conversion, so it is unchanged by it.)")

# ----------------------------------------------------------------------------
hdr("CHECK 2  vertical profile")
cr = int(np.ceil(CYL_R)) + 1
ys = slice(int(c[1]) - cr, int(c[1]) + cr + 2)
xs = slice(int(c[2]) - cr, int(c[2]) + cr + 2)
yv, xv = np.arange(ys.start, ys.stop) - o[1], np.arange(xs.start, xs.stop) - o[2]
Rh = np.sqrt(yv[:, None] ** 2 + xv[None, :] ** 2)
zmax = np.sqrt(R_CLIP ** 2 - CYL_R ** 2) - 2.0
edges = np.arange(0, zmax + ZBIN, ZBIN)
nb = len(edges) - 1
z_c, mean_z, med_z, cnt = 0.5 * (edges[1:] + edges[:-1]), [], [], []
vals = [[] for _ in range(nb)]
zs_all = np.arange(isrf.shape[0]) - o[0]
for z0 in range(0, isrf.shape[0], 50):
    zsl = slice(z0, min(z0 + 50, isrf.shape[0]))
    zv = zs_all[zsl]
    if np.all(np.abs(zv) > zmax):
        continue
    cube = np.asarray(isrf[zsl, ys, xs], dtype=np.float64)
    for iz, zc_ in enumerate(np.abs(zv)):
        j = np.searchsorted(edges, zc_, side="right") - 1
        if 0 <= j < nb and np.sqrt(zc_ ** 2 + CYL_R ** 2) <= R_CLIP:
            vals[j].append(cube[iz][Rh <= CYL_R])
for j in range(nb):
    a = np.concatenate(vals[j])
    mean_z.append(a.mean()); med_z.append(np.median(a)); cnt.append(a.size)
mean_z, med_z = np.array(mean_z), np.array(med_z)
print(f"  footprint: cylinder of horizontal radius {CYL_R:g} pc about the Sun, both hemispheres "
      f"folded; |z| = 0-{edges[-1]:.0f} pc in {ZBIN:g} pc bins")
print(f"  LIMIT: the r{R_CLIP:.0f} clip caps |z| at {edges[-1]:.0f} pc for this footprint; "
      f"the paper's 400 pc scale height (and its 800 pc point) cannot be sampled out to 1-2 h.")


def fit(z, y, lo, hi):
    s = (z >= lo) & (z <= hi)
    A = np.vstack([np.ones(s.sum()), z[s]]).T
    coef, *_ = np.linalg.lstsq(A, np.log(y[s]), rcond=None)
    res = np.log(y[s]) - A @ coef
    r2 = 1 - res.var() / np.log(y[s]).var()
    return -1 / coef[1], np.exp(coef[0]), np.sqrt(np.mean(res ** 2)), r2


print(f"  {'stat':<7} {'|z| range':>12} {'h [pc]':>8} {'I0':>7} {'rms ln-resid':>13} {'R^2':>6}")
fits_ = {}
for name, y in [("mean", mean_z), ("median", med_z)]:
    for lo, hi in [(0, edges[-1]), (0, 150), (150, 300), (300, edges[-1])]:
        hh, i0, rms, r2 = fit(z_c, y, lo, hi)
        fits_[(name, lo, hi)] = hh
        print(f"  {name:<7} {lo:5.0f}-{hi:5.0f} pc {hh:8.0f} {i0:7.3f} {rms:13.3f} {r2:6.3f}")
for zz_ in (0, 100, 200, 300, 400, 470):
    j = min(int(zz_ // ZBIN), nb - 1)
    print(f"    |z|~{z_c[j]:4.0f} pc: mean {mean_z[j]:.3f}, median {med_z[j]:.3f}  ({cnt[j]} voxels)")
h_full = fits_[("median", 0, edges[-1])]
RESULT[2] = 0.75 * H_EXPECTED <= h_full <= 1.25 * H_EXPECTED
print(f"  paper: 0.633 at the Sun and 12% of Standard at 800 pc imply an effective "
      f"h = 800/ln(0.633/0.12) = {800/np.log(0.633/0.12):.0f} pc over 0-800 pc")
print(f"  median-profile fit over full range: h = {h_full:.0f} pc vs {H_EXPECTED:.0f} (+-25%): "
      f"{'PASS' if RESULT[2] else 'FAIL'}")

# ----------------------------------------------------------------------------
hdr("CHECK 3  Ophiuchus vs Taurus, n > 30 cm^-3")
print("  Volumes: Zucker et al. 2021 Table 1 skeleton x/y/z extents (NOT envelope boundaries;")
print("  Bish+2026 state no boundaries; Zucker+21 3D masks doi:10.7910/DVN/NMRLP1 not available).")
stats = {}
for name, (x0, x1, y0, y1, z0, z1) in CLOUDS.items():
    corners = np.array([[x, y, z] for x in (x0, x1) for y in (y0, y1) for z in (z0, z1)], float)
    rmax = np.sqrt((corners ** 2).sum(1)).max()
    if rmax > R_CLIP - 2:
        sys.exit(f"{name} box reaches r={rmax:.0f} pc, outside the r{R_CLIP:.0f} clip. "
                 "Fix: use the full unclipped LightCube/Edenhofer cubes (hannahbish.me/LightCubes).")
    # voxels whose centres lie inside the box
    ix = np.arange(int(np.ceil(o[2] + x0)), int(np.floor(o[2] + x1)) + 1)
    iy = np.arange(int(np.ceil(o[1] + y0)), int(np.floor(o[1] + y1)) + 1)
    iz = np.arange(int(np.ceil(o[0] + z0)), int(np.floor(o[0] + z1)) + 1)
    sl = (slice(iz[0], iz[-1] + 1), slice(iy[0], iy[-1] + 1), slice(ix[0], ix[-1] + 1))
    I = np.nan_to_num(np.asarray(isrf[sl], dtype=np.float64))
    tau = np.nan_to_num(np.asarray(dust[sl], dtype=np.float64))
    n = tau * nh_per_tau()
    sel = n > N_CUT
    s = I[sel]
    stats[name] = dict(n_box=I.size, n_sel=sel.sum(), mean=s.mean() if s.size else np.nan,
                       median=np.median(s) if s.size else np.nan,
                       p5=np.percentile(s, 5) if s.size else np.nan,
                       p95=np.percentile(s, 95) if s.size else np.nan, nmax=n.max(), I=I, n=n,
                       tau=tau)
    print(f"  {name}: box x[{x0},{x1}] y[{y0},{y1}] z[{z0},{z1}] pc, max r = {rmax:.0f} pc (inside clip)")
    print(f"    voxels in box {I.size}; n>{N_CUT:g}: {sel.sum()} ({sel.mean():.2%}); max n {n.max():.1f}")
    print(f"    I_UV over n>{N_CUT:g}: mean {stats[name]['mean']:.4f}, median {stats[name]['median']:.4f}, "
          f"5-95% {stats[name]['p5']:.4f}-{stats[name]['p95']:.4f}")
    print(f"    I_UV over whole box: mean {I.mean():.4f}, median {np.median(I):.4f}")
O, T = stats["Ophiuchus"], stats["Taurus"]
r_mean, r_med = O["mean"] / T["mean"], O["median"] / T["median"]
print(f"  Ophiuchus/Taurus: ratio of means {r_mean:.2f}, ratio of medians {r_med:.2f}  (expected ~{RATIO_EXPECTED:g})")

print("  sensitivity (diagnosing a miss):")
print("    (c) density cut placement, voxel counts and ratio of means vs cut and conversion:")
for conv in NH_CONVERSIONS:
    row = []
    for cut in (N_CUT,) + N_CUT_ALT:
        cnts, means = [], []
        for nm in ("Ophiuchus", "Taurus"):
            d = stats[nm]
            sel = d["tau"] * nh_per_tau(conv) > cut
            cnts.append(sel.sum()); means.append(d["I"][sel].mean() if sel.any() else np.nan)
        row.append(f"n>{cut:g}: O {cnts[0]:5d} T {cnts[1]:5d} ratio {means[0]/means[1]:.2f}")
    print(f"      {conv:<17} " + " | ".join(row))
print("    (a) boundary dependence, whole-box (no density cut) ratio of means: "
      f"{O['I'].mean()/T['I'].mean():.2f}, of medians {np.median(O['I'])/np.median(T['I']):.2f}")
RESULT[3] = 2.0 <= r_mean <= 4.0
print(f"  ratio of means {r_mean:.2f} within 2-4 (\"roughly 3x\"): {'PASS' if RESULT[3] else 'FAIL'}")

hdr("SUMMARY")
for i, t in [(1, "solar value"), (2, "vertical scale height"), (3, "Oph/Tau ratio")]:
    print(f"  Check {i} ({t}): {'PASS' if RESULT[i] else 'FAIL'}")
