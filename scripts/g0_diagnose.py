"""Diagnostic: is the structure in G0_surface caused by the boundary finder
(r0 jumping between density features) or by the ISRF cube itself?

Compares, for the Edenhofer maps (the only map source), G0 sampled at r0 (G0_surface), at r_peak
(G0_peak) and at a fixed distance (field-median r_peak, G0_fixed).
Reads the existing eos_maps outputs; changes nothing.
"""
import numpy as np
from astropy.io import fits
from scipy.stats import spearmanr, pearsonr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cloudfield.cubes import open_cubes, origin_index, los_xyz, sample_los
from cloudfield.regions import healpix_on_grid
from cloudfield.paths import data_path
from cloudfield.config import load_cloud
from cloudfield.maps import Run

OUT_DIR = "output"           # where `cloudfield run eos` wrote its FITS
R = Run(load_cloud("eos"), verbose=False)


class _EM:
    pass


em = _EM()
g = R.grid
em.LL2, em.BB2, em.PIX, em.NB, em.NL = g.LL2, g.BB2, g.pix, g.NB, g.NL
em.lvals, em.bvals, em.PAD, em.ORIGIN_MODE = g.lvals, g.bvals, R.cfg["pad"], R.cfg["origin_mode"]
_reg = R.cfg["region"]
em.FIMS_BOX, em.FIG_MARGIN = _reg["pick_box"], R.cfg["figure"]["margin_deg"]
em.healpix_on_grid = healpix_on_grid
em.H2_MAP, em.H2_COLUMN, em.H2_LEVEL = str(data_path(_reg["data_key"])), _reg["column"], _reg["level"]

SOURCES = ("edenhofer",)   # DUSTRIBUTION is not a map source
JUMP_PC = 5.0
G0_LIM = (0.30, 0.70)        # fixed, linear
SEAMS = dict(l43=(42.5, 43.5), l36=(35.8, 37.4))   # columns to check, deg
LOBES = {"(40,44)": (40.0, 44.0), "(40,45)": (40.0, 45.0)}
OUT_PNG = f"{OUT_DIR}/eos_g0_diagnostic.png"

CLOUD = fits.getdata(f"{OUT_DIR}/eos_maps.fits", "CLOUD_H2_10PCT").astype(bool)
isrf_mm, dust_mm, hdr = open_cubes(data_path("lightcube"), data_path("edenhofer_tau"))
origin = origin_index(isrf_mm.shape, hdr, em.ORIGIN_MODE)
LL, BB = em.LL2, em.BB2
PIX = em.PIX


def grad_mag(a):
    """|grad| per pixel (central differences along l and b; NaN-propagating)."""
    gb, gl = np.gradient(a)
    return np.hypot(gl, gb)


def sample_fixed(r_fix):
    out = np.full(LL.shape, np.nan)
    for j0 in range(0, em.NB, 40):
        j1 = min(j0 + 40, em.NB)
        X, Y, Z = los_xyz(LL[j0:j1].ravel(), BB[j0:j1].ravel(), np.array([r_fix]))
        _, I = sample_los(isrf_mm, None, origin, X, Y, Z, pad=em.PAD)
        out[j0:j1] = I[:, 0].reshape(j1 - j0, em.NL)
    return out


D = {}
for s in SOURCES:
    r0 = fits.getdata(f"{OUT_DIR}/eos_maps.fits", "R0").astype(float)
    rp = fits.getdata(f"{OUT_DIR}/eos_maps.fits", "RPEAK").astype(float)
    g0s = fits.getdata(f"{OUT_DIR}/eos_G0.fits", "G0_SURFACE").astype(float)
    g0p = fits.getdata(f"{OUT_DIR}/eos_G0.fits", "G0_PEAK").astype(float)
    r_fix = float(np.nanmedian(rp[CLOUD]))
    g0f = sample_fixed(r_fix)
    D[s] = dict(r0=r0, rp=rp, g0s=g0s, g0p=g0p, g0f=g0f, r_fix=r_fix)

for s in SOURCES:
    d = D[s]
    print("=" * 78)
    print(f"[{s}]  (statistics inside the 10% contour, {CLOUD.sum()} px; pixel = {PIX} deg)")
    # 1. gradient correlation
    gr0 = grad_mag(d["r0"])
    m = CLOUD & np.isfinite(gr0)
    print(f"  1. |grad r0| vs |grad G0| (Spearman rho / Pearson r):")
    for k, lab in [("g0s", "G0_surface (at r0)"), ("g0p", "G0_peak (at r_peak)"),
                   ("g0f", f"G0_fixed (r={d['r_fix']:.1f} pc)")]:
        g = grad_mag(d[k])
        mm = m & np.isfinite(g)
        print(f"     {lab:<26} rho={spearmanr(gr0[mm], g[mm])[0]:+.3f}  r={pearsonr(gr0[mm], g[mm])[0]:+.3f}"
              f"   median |grad G0| = {np.median(g[mm]):.4f}/px, 95th pct {np.percentile(g[mm],95):.4f}")
    # 2. adjacent jumps
    dl = np.abs(np.diff(d["r0"], axis=1))           # along l
    db = np.abs(np.diff(d["r0"], axis=0))           # along b
    ml = CLOUD[:, 1:] & CLOUD[:, :-1]
    mb = CLOUD[1:, :] & CLOUD[:-1, :]
    jl, jb = dl[ml & np.isfinite(dl)], db[mb & np.isfinite(db)]
    allj = np.concatenate([jl, jb])
    q = np.percentile(allj, [50, 90, 99, 99.9])
    print(f"  2. |r0(i+1)-r0(i)| over {allj.size} adjacent pairs: median {q[0]:.2f}, 90% {q[1]:.2f}, "
          f"99% {q[2]:.2f}, 99.9% {q[3]:.2f}, max {allj.max():.1f} pc")
    print(f"     fraction > {JUMP_PC:g} pc: along l {np.mean(jl>JUMP_PC):.2%}, along b {np.mean(jb>JUMP_PC):.2%}, "
          f"all {np.mean(allj>JUMP_PC):.2%}")
    jump_px = np.zeros_like(CLOUD)
    jump_px[:, 1:] |= np.nan_to_num(dl) > JUMP_PC
    jump_px[:, :-1] |= np.nan_to_num(dl) > JUMP_PC
    jump_px[1:, :] |= np.nan_to_num(db) > JUMP_PC
    jump_px[:-1, :] |= np.nan_to_num(db) > JUMP_PC
    print(f"     pixels adjacent to a >{JUMP_PC:g} pc jump: {np.mean(jump_px[CLOUD]):.2%} of contour")
    # coincidence of strongest G0_surface edges with r0 jumps
    g = grad_mag(d["g0s"])
    top = CLOUD & np.isfinite(g) & (g >= np.nanpercentile(g[CLOUD], 95))
    print(f"     of the top-5% |grad G0_surface| pixels, {np.mean(jump_px[top]):.1%} sit on an r0 jump "
          f"(base rate {np.mean(jump_px[CLOUD]):.1%})")
    # specific features
    for nm, (a, b) in SEAMS.items():
        col = (LL >= a) & (LL <= b) & (BB >= 38) & (BB <= 62)
        print(f"     seam {nm} (l={a}-{b}): r0-jump fraction {np.mean(jump_px[col]):.1%}; "
              f"median |grad G0_surface| {np.nanmedian(g[col]):.4f} vs contour {np.nanmedian(g[CLOUD]):.4f}")
    for nm, (l0, b0) in LOBES.items():
        k = (np.abs(LL - l0) <= 0.75) & (np.abs(BB - b0) <= 0.75)
        print(f"     lobe {nm}: G0_surface {np.nanmedian(d['g0s'][k]):.3f}, G0_peak {np.nanmedian(d['g0p'][k]):.3f}, "
              f"G0_fixed {np.nanmedian(d['g0f'][k]):.3f};  r0 {np.nanmedian(d['r0'][k]):.1f} pc "
              f"(range {np.nanmin(d['r0'][k]):.1f}-{np.nanmax(d['r0'][k]):.1f}), r_peak {np.nanmedian(d['rp'][k]):.1f}")
    # 3. smoothness comparison
    print(f"  3. smoothness (median |grad| inside contour, per 0.1 deg px) and range:")
    for k, lab in [("g0s", "G0_surface"), ("g0p", "G0_peak"), ("g0f", "G0_fixed")]:
        v = d[k][CLOUD]
        print(f"     {lab:<11} median|grad| {np.nanmedian(grad_mag(d[k])[CLOUD]):.4f}   "
              f"5-95% {np.nanpercentile(v,5):.3f}-{np.nanpercentile(v,95):.3f}  "
              f"full {np.nanmin(v):.3f}-{np.nanmax(v):.3f}")

# ----------------------------------------------------------------------------
fl = (em.FIMS_BOX["l"][0] - em.FIG_MARGIN, em.FIMS_BOX["l"][1] + em.FIG_MARGIN)
fb = (em.FIMS_BOX["b"][0] - em.FIG_MARGIN, em.FIMS_BOX["b"][1] + em.FIG_MARGIN)
aspect = 1.0 / np.cos(np.radians(np.mean(fb)))
extent = [em.lvals[0] + PIX / 2, em.lvals[-1] - PIX / 2, em.bvals[0] - PIX / 2, em.bvals[-1] + PIX / 2]
H2 = em.healpix_on_grid(em.H2_MAP, em.H2_COLUMN, em.lvals, em.bvals)
fig, axs = plt.subplots(len(SOURCES), 4, figsize=(18, 6.5 * len(SOURCES)), sharex=True, sharey=True,
                        squeeze=False)
for row, s in enumerate(SOURCES):
    d = D[s]
    panels = [("r0", "r0 [pc]", "cividis", (85, 135)),
              ("g0s", "G0_surface = I_UV(r0)", "magma", G0_LIM),
              ("g0p", "G0_peak = I_UV(r_peak)", "magma", G0_LIM),
              ("g0f", f"G0_fixed = I_UV({d['r_fix']:.1f} pc)", "magma", G0_LIM)]
    for ax, (k, lab, cm, lim) in zip(axs[row], panels):
        cmap = plt.get_cmap(cm).copy(); cmap.set_bad("0.85")
        im = ax.imshow(d[k], origin="lower", extent=extent, aspect=aspect, cmap=cmap,
                       vmin=lim[0], vmax=lim[1], interpolation="nearest")
        ax.contour(em.lvals, em.bvals, H2, levels=[em.H2_LEVEL], colors="c", linewidths=1.0)
        v = d[k][CLOUD]
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.06, fraction=0.05)
        cb.set_label(f"{lab}\ncontour range {np.nanmin(v):.3f}-{np.nanmax(v):.3f}, "
                     f"5-95% {np.nanpercentile(v,5):.3f}-{np.nanpercentile(v,95):.3f}", fontsize=8)
        ax.set_title(f"{s}: {k}", fontsize=9)
    axs[row, 0].set_ylabel("b [deg]")
for ax in axs[-1]:
    ax.set_xlabel("l [deg]")
axs[0, 0].set_xlim(fl[1], fl[0]); axs[0, 0].set_ylim(*fb)
fig.suptitle(f"G0 diagnostic: linear fixed scale {G0_LIM[0]}-{G0_LIM[1]} (Draine) for all G0 panels; "
             "cyan = 10% H2/FUV contour", fontsize=11)
fig.savefig(OUT_PNG, dpi=130, bbox_inches="tight")
print(f"wrote {OUT_PNG}")
