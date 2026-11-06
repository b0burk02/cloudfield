"""Diagnostic: is the dust-column vs GALFA N_HI deficit a comparison artifact?

1. Restrict GALFA N_HI to the Eos velocity component, chosen from the data:
   median spectrum inside vs outside the 10% FIMS/SPEAR H2/FUV contour.
2. Truncation: extend Edenhofer rays from 260 to ~495 pc.
No conversion constants are changed.
"""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cloudfield.constants import PC_CM, R_NO_DUST
from cloudfield.cubes import open_cubes, origin_index, los_xyz, sample_los, rng
from cloudfield.conversions import nh_per_tau
from cloudfield.regions import healpix_on_grid
from cloudfield.paths import data_path
from cloudfield.config import load_cloud
from cloudfield.maps import Run
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import map_coordinates, uniform_filter

OUT_DIR = "output"           # where `cloudfield run eos` wrote eos_maps.fits
R = Run(load_cloud("eos"), verbose=False)
GALFA_FILE = str(data_path("galfa_nhi"))
GALFA_CHECK = str(data_path("galfa_check"))
tau_to_n = lambda t: t * R.nh_per_tau


class _EM:                    # names the original script used from eos_maps
    pass


em = _EM()

CUBE = str(data_path("galfa_cube"))
WIN_FRAC = 0.10              # window: contiguous channels with diff >= WIN_FRAC * peak
WIN_FRAC_ALT = (0.05, 0.20)  # sensitivity
NHI_PER_KKMS = 1.8224e18     # cm^-2 / (K km/s), optically thin
R_EXT = 495.0                # pc; Edenhofer cube here is clipped at r=500 pc
DR = R.dr
OUT_PNG = f"{OUT_DIR}/eos_hi_spectra.png"

def galfa_on_grid(lvals, bvals):
    """GALFA-HI DR2 stray-corrected N_HI (-90..+90 km/s) onto the (l, b) grid.
    The file's CTYPEs are bare 'RA'/'DEC' (linear in RA, Dec: plate carree, J2000),
    so pixel coordinates are computed from the linear header relation explicitly.
    Boxcar-averaged to the output pixel scale (NaN-aware), then bilinear."""
    with fits.open(GALFA_FILE, memmap=True) as hl:
        h = hl[0].header
        if (h["CTYPE1"].strip(), h["CTYPE2"].strip()) != ("RA", "DEC"):
            raise SystemExit(f"GALFA CTYPEs {h['CTYPE1']!r}, {h['CTYPE2']!r}: expected bare RA/DEC")
        LL, BB = np.meshgrid(lvals, bvals)
        c = SkyCoord(l=LL * u.deg, b=BB * u.deg, frame="galactic").fk5
        dra = (c.ra.deg - h["CRVAL1"] + 180.0) % 360.0 - 180.0
        x = dra / h["CDELT1"] + h["CRPIX1"] - 1.0
        y = (c.dec.deg - h["CRVAL2"]) / h["CDELT2"] + h["CRPIX2"] - 1.0
        nbox = max(int(round(R.grid.pix / abs(h["CDELT2"]))), 1)
        x0, x1 = max(int(np.floor(x.min())) - nbox, 0), min(int(np.ceil(x.max())) + nbox + 1, h["NAXIS1"])
        y0, y1 = max(int(np.floor(y.min())) - nbox, 0), min(int(np.ceil(y.max())) + nbox + 1, h["NAXIS2"])
        sub = np.asarray(hl[0].data[y0:y1, x0:x1], dtype=np.float64)
    print(f"  GALFA {GALFA_FILE}: {h['NAXIS1']}x{h['NAXIS2']}, BUNIT={h.get('BUNIT')!r}, "
          f"EQUINOX={h.get('EQUINOX')}; crop [{y0}:{y1}, {x0}:{x1}], boxcar {nbox} px")
    good = np.isfinite(sub) & (sub > 0)
    num = uniform_filter(np.where(good, sub, 0.0), nbox)
    den = uniform_filter(good.astype(float), nbox)
    with np.errstate(invalid="ignore", divide="ignore"):
        sm = np.where(den > 0.5, num / den, np.nan)
    out = map_coordinates(sm, [y.ravel() - y0, x.ravel() - x0], order=1,
                          mode="constant", cval=np.nan).reshape(LL.shape)
    rng("GALFA N_HI", out)
    # frame check against the Galactic-TAN GALFA mosaic (different velocity cut)
    try:
        from reproject import reproject_interp
        with fits.open(GALFA_CHECK) as hl:
            d2 = np.squeeze(hl[0].data).astype(float)
            w2 = WCS(hl[0].header).celestial
        chk, _ = reproject_interp((d2, w2), R.grid.wcs(), shape_out=out.shape)
        m = np.isfinite(chk) & np.isfinite(out) & (chk > 0) & (out > 0)
        print(f"  frame check vs {GALFA_CHECK} (GLON-TAN): corr(log) = "
              f"{np.corrcoef(np.log10(chk[m]), np.log10(out[m]))[0,1]:.3f} over {m.sum()} px; "
              f"median ratio {np.median(out[m]/chk[m]):.3f} (velocity ranges differ)")
    except Exception as e:  # check only; the main product does not depend on it
        print(f"  frame check skipped: {type(e).__name__}: {e}")
    return out



em.map_wcs = R.grid.wcs
em.lvals, em.bvals, em.NL, em.NB = R.grid.lvals, R.grid.bvals, R.grid.NL, R.grid.NB
em.R_MIN, em.R_MAX, em.ORIGIN_MODE, em.PAD = R.r[0], R.r[-1], R.cfg["origin_mode"], R.cfg["pad"]
em.healpix_on_grid = healpix_on_grid
em.galfa_on_grid = galfa_on_grid
reg = R.cfg["region"]
em.H2_MAP, em.H2_COLUMN, em.H2_LEVEL = str(data_path(reg["data_key"])), reg["column"], reg["level"]
em.PIX = R.grid.pix
PIX = R.grid.pix

# ----------------------------------------------------------------------------
# masks on the map grid
# ----------------------------------------------------------------------------
w_map = em.map_wcs()
CLOUD = fits.getdata(f"{OUT_DIR}/eos_maps.fits", "CLOUD_H2_10PCT").astype(bool)
NHF = {"edenhofer": fits.getdata(f"{OUT_DIR}/eos_maps.fits", "NH_FULL").astype(float)}   # Edenhofer only (DUSTRIBUTION columns: eos_column_systematics.py)
print("H2 map for the outside mask:")
H2 = em.healpix_on_grid(em.H2_MAP, em.H2_COLUMN, em.lvals, em.bvals)
OUTSIDE = np.isfinite(H2) & (H2 < em.H2_LEVEL)
print(f"inside 10% contour: {CLOUD.sum()} px; outside (H2/FUV < 10%): {OUTSIDE.sum()} px")

# ----------------------------------------------------------------------------
# 1. GALFA cube: spectra inside / outside
# ----------------------------------------------------------------------------
print("=" * 78)
with fits.open(CUBE, memmap=True) as hl:
    h = hl[0].header
    cube = np.asarray(hl[0].data, dtype=np.float32)
cube[~np.isfinite(cube)] = np.nan
nv = h["NAXIS3"]
vel = (h["CRVAL3"] + (np.arange(nv) + 1 - h["CRPIX3"]) * h["CDELT3"]) / 1e3   # km/s
dv = abs(h["CDELT3"]) / 1e3
w_cube = WCS(h).celestial
print(f"GALFA cube {CUBE}: shape {cube.shape}, BUNIT={h['BUNIT']!r}, {h['CTYPE1']}/{h['CTYPE2']}, "
      f"v = {vel[0]:.2f}..{vel[-1]:.2f} km/s, dv = {dv:.4f} km/s ({h['SPECSYS']})")
rng("cube Tb", cube[::10])

# map each cube pixel to the map grid (nearest) for mask lookup
ny, nx = cube.shape[1:]
yy, xx = np.mgrid[0:ny, 0:nx]
lon, lat = w_cube.wcs_pix2world(xx, yy, 0)
lon = np.mod(lon, 360.0)
px, py = w_map.wcs_world2pix(lon, lat, 0)
ix, iy = np.rint(px).astype(int), np.rint(py).astype(int)
ongrid = (ix >= 0) & (ix < em.NL) & (iy >= 0) & (iy < em.NB)
cin = np.zeros((ny, nx), bool)
cout = np.zeros((ny, nx), bool)
cin[ongrid] = CLOUD[iy[ongrid], ix[ongrid]]
cout[ongrid] = OUTSIDE[iy[ongrid], ix[ongrid]]
print(f"cube pixels: inside contour {cin.sum()}, outside {cout.sum()}, "
      f"map-grid coverage by cube {np.mean(ongrid):.1%} of cube px")
cov = np.zeros_like(CLOUD)
cov[iy[ongrid], ix[ongrid]] = True
print(f"fraction of the contour covered by the cube: {np.mean(cov[CLOUD]):.1%}")

spec_in = np.nanmedian(cube[:, cin], axis=1)
spec_out = np.nanmedian(cube[:, cout], axis=1)
diff = spec_in - spec_out
ipk = int(np.nanargmax(diff))


def window(frac):
    thr = frac * diff[ipk]
    a = ipk
    while a > 0 and diff[a - 1] >= thr:
        a -= 1
    b = ipk
    while b < nv - 1 and diff[b + 1] >= thr:
        b += 1
    return a, b


i0, i1 = window(WIN_FRAC)
print(f"difference spectrum peak: {diff[ipk]:.3f} K at v = {vel[ipk]:.2f} km/s")
print(f"Eos window ({WIN_FRAC:.0%} of peak): v = {vel[i0]:.2f} .. {vel[i1]:.2f} km/s "
      f"({i1-i0+1} ch)")
for f in WIN_FRAC_ALT:
    a, b = window(f)
    print(f"  sensitivity ({f:.0%} of peak): v = {vel[a]:.2f} .. {vel[b]:.2f} km/s")
wtot = lambda s, a, b: np.nansum(s[a:b + 1]) * dv * NHI_PER_KKMS
print(f"median-spectrum columns: inside full {wtot(spec_in,0,nv-1):.3g}, window "
      f"{wtot(spec_in,i0,i1):.3g}; outside full {wtot(spec_out,0,nv-1):.3g}, window "
      f"{wtot(spec_out,i0,i1):.3g} cm^-2")
print(f"  share of inside median spectrum in window: {wtot(spec_in,i0,i1)/wtot(spec_in,0,nv-1):.1%}")

fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                             gridspec_kw=dict(height_ratios=[2, 1]))
a1.plot(vel, spec_in, color="C3", label=f"inside {em.H2_LEVEL:g}% H2/FUV contour (median)")
a1.plot(vel, spec_out, color="C0", label="outside (H2/FUV < 10%, median)")
a2.plot(vel, diff, color="k", label="inside - outside")
for ax in (a1, a2):
    ax.axvspan(vel[i0], vel[i1], color="gold", alpha=0.3, lw=0)
    ax.axhline(0, color="0.6", lw=0.6)
a2.axhline(WIN_FRAC * diff[ipk], color="0.4", ls=":", lw=0.8,
           label=f"{WIN_FRAC:.0%} of peak -> window {vel[i0]:.1f}..{vel[i1]:.1f} km/s")
a1.set_ylabel(r"$T_b$ [K]"); a2.set_ylabel(r"$\Delta T_b$ [K]")
a2.set_xlabel(r"$v_{\rm LSR}$ [km s$^{-1}$]")
a1.legend(fontsize=8); a2.legend(fontsize=8)
a1.set_xlim(-60, 40)
a1.set_title(f"GALFA-HI ({CUBE.split('/')[-1]})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=150)
print(f"wrote {OUT_PNG}")

# integrated maps -> map grid
mom_full = np.nansum(cube, axis=0) * dv * NHI_PER_KKMS
NHI_full_cube, _ = reproject_interp((mom_full, w_cube), w_map, shape_out=CLOUD.shape)
NHI_ALT = {}
for f in (WIN_FRAC,) + WIN_FRAC_ALT:
    a, b = window(f)
    mom = np.nansum(cube[a:b + 1], axis=0) * dv * NHI_PER_KKMS
    NHI_ALT[f], _ = reproject_interp((mom, w_cube), w_map, shape_out=CLOUD.shape)
NHI_win = NHI_ALT[WIN_FRAC]
del cube
print("GALFA N_HI used previously (stray-corrected, -90..+90 km/s, RA/Dec):")
NHI_old = em.galfa_on_grid(em.lvals, em.bvals)
m = CLOUD & np.isfinite(NHI_full_cube) & np.isfinite(NHI_old)
print(f"check: cube full-range ({vel[0]:.0f}..{vel[-1]:.0f} km/s) / stray-corrected -90..+90 map "
      f"inside contour: median {np.median(NHI_full_cube[m]/NHI_old[m]):.3f}")

# ----------------------------------------------------------------------------
# 2. Truncation: Edenhofer rays to R_EXT
# ----------------------------------------------------------------------------
print("=" * 78)
print(f"Truncation test: Edenhofer rays {em.R_MIN:g} -> {R_EXT:g} pc (inside contour only)")
_, dust_mm, hdr = open_cubes(data_path("lightcube"), data_path("edenhofer_tau"))
origin = origin_index(dust_mm.shape, hdr, em.ORIGIN_MODE)
r = np.arange(em.R_MIN, R_EXT + DR / 2, DR)
jj, ii = np.where(CLOUD)                      # row-major -> compact chunks in b
NH_ext = np.full(CLOUD.shape, np.nan)
NH_260 = np.full(CLOUD.shape, np.nan)
NH_350 = np.full(CLOUD.shape, np.nan)
k260, k350 = r <= em.R_MAX, r <= 350.0
CH = 3000
for s in range(0, jj.size, CH):
    j, i = jj[s:s + CH], ii[s:s + CH]
    X, Y, Z = los_xyz(em.lvals[i], em.bvals[j], r)
    # sample_los returns (dust, isrf); passing the dust cube in the isrf slot with
    # dust_mm=None samples only that cube, with identical crop/interp settings.
    _, tau = sample_los(dust_mm, None, origin, X, Y, Z, pad=em.PAD)
    if s == 0:
        rng("tau_pc (Edenhofer)", tau)
    n = tau_to_n(tau)
    n[:, r < R_NO_DUST] = 0.0
    col = lambda k: n[:, k].sum(1) * DR * PC_CM
    NH_ext[j, i] = col(slice(None))
    NH_260[j, i] = col(k260)
    NH_350[j, i] = col(k350)
chk = NH_260[CLOUD] / NHF["edenhofer"][CLOUD]
print(f"consistency: re-integrated 69-260 / NH_full from eos_maps.fits: "
      f"median {np.median(chk):.5f}, max |dev| {np.max(np.abs(chk-1)):.2e}")
g = NH_ext[CLOUD] / NH_260[CLOUD]
print(f"NH(69-{R_EXT:g}) / NH(69-260): median {np.median(g):.3f} (16-84% "
      f"{np.percentile(g,16):.3f}-{np.percentile(g,84):.3f});  NH(69-350)/NH(69-260): "
      f"{np.median(NH_350[CLOUD]/NH_260[CLOUD]):.3f}")
print(f"median NH: 69-260 {np.median(NH_260[CLOUD]):.3g}, 69-350 {np.median(NH_350[CLOUD]):.3g}, "
      f"69-{R_EXT:g} {np.median(NH_ext[CLOUD]):.3g} cm^-2")

# ----------------------------------------------------------------------------
# 3. Comparison table
# ----------------------------------------------------------------------------
print("=" * 78)
print(f"Dust column / GALFA N_HI inside the {em.H2_LEVEL:g}% contour")
HI = {"N_HI -90..+90 (stray-corr., used before)": NHI_old,
      f"N_HI cube {vel[0]:.0f}..{vel[-1]:.0f}": NHI_full_cube,
      f"N_HI Eos window {vel[i0]:.1f}..{vel[i1]:.1f}": NHI_win}
for k, v in HI.items():
    print(f"  {k:<42} median {np.nanmedian(v[CLOUD]):.3g} cm^-2")
DUST = {"edenhofer NH 69-260": NHF["edenhofer"],
        f"edenhofer NH 69-{R_EXT:g}": NH_ext}
base = list(HI)[0]
print(f"\n  {'dust column':<26} {'N_HI':<42} {'med NH/NHI':>10} {'NH<NHI':>7} "
      f"{'deficit closed':>14}")
res = {}
for dk, D in DUST.items():
    D0 = NHF[dk.split()[0]]                    # same source, 69-260 pc: the original baseline
    m0 = CLOUD & np.isfinite(D0) & np.isfinite(HI[base])
    r_old = np.median(D0[m0] / HI[base][m0])
    for hk, Hm in HI.items():
        m = CLOUD & np.isfinite(D) & np.isfinite(Hm)
        q = D[m] / Hm[m]
        rm = np.median(q)
        closed = np.log(rm / r_old) / np.log(1 / r_old)
        res[(dk, hk)] = rm
        print(f"  {dk:<26} {hk:<42} {rm:10.3f} {np.mean(q<1):7.1%} {closed:14.0%}")
    print()
print("  'deficit closed' = log(new/old) / log(1/old); old = same source, 69-260 pc, vs the "
      "-90..+90 map")
print("  window sensitivity, Edenhofer NH 69-495 / N_HI(window):")
for f, Hm in NHI_ALT.items():
    a, b = window(f)
    m = CLOUD & np.isfinite(NH_ext) & np.isfinite(Hm)
    q = NH_ext[m] / Hm[m]
    print(f"    {f:4.0%} of peak ({vel[a]:6.2f}..{vel[b]:5.2f} km/s): median {np.median(q):.3f}, "
          f"NH<NHI {np.mean(q<1):.1%}")
ed_old = res[("edenhofer NH 69-260", base)]
ed_win = res[("edenhofer NH 69-260", list(HI)[2])]
ed_ext_old = res[(f"edenhofer NH 69-{R_EXT:g}", base)]
ed_both = res[(f"edenhofer NH 69-{R_EXT:g}", list(HI)[2])]
L = lambda a, b: np.log(a / b) / np.log(1 / ed_old)
print(f"Edenhofer decomposition of the deficit (log space, total = 1/{ed_old:.3f} = {1/ed_old:.2f}x):")
print(f"  velocity restriction alone : {L(ed_win, ed_old):6.0%}")
print(f"  truncation alone           : {L(ed_ext_old, ed_old):6.0%}")
print(f"  both                       : {L(ed_both, ed_old):6.0%}   -> residual NH/NHI = {ed_both:.3f}")
