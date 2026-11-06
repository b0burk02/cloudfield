"""Column-density systematics between the Edenhofer map (our map source) and
DUSTRIBUTION (Burkhart et al. 2025's map). COMPARISON ONLY: prints numbers,
writes no map, computes no G0 or alphaG. DUSTRIBUTION is never a map source
(LightCube's I_UV was computed through the Edenhofer map).

All statistics inside the 10% FIMS/SPEAR H2/FUV contour (CLOUD_H2_10PCT in
eos_maps.fits). Requires eos_maps.py to have been run.
"""
import numpy as np
from astropy.io import fits
from cloudfield.constants import PC_CM
from cloudfield.columns import (NH_PER_MAGPC, load_dustribution_density, sample_dustribution,
                                dustribution_index, DUSTRIBUTION_SHAPE, DUSTRIBUTION_DMAX)
from cloudfield.paths import data_path
from cloudfield.config import load_cloud
from cloudfield.maps import Run

OUT_DIR = "output"           # where `cloudfield run eos` wrote eos_maps.fits
R = Run(load_cloud("eos"), verbose=False)


class _EM:
    pass


em = _EM()
em.LL2, em.BB2, em.R_MIN, em.R_MAX, em.DR = R.grid.LL2, R.grid.BB2, R.r[0], R.r[-1], R.dr

EXT_CUBE = str(data_path("dustribution_cumulative"))   # DUSTRIBUTION cumulative extinction (mag)
# Burkhart et al. (2025): Eos mass from DUSTRIBUTION vs from Edenhofer (quoted, not recomputed)
MASS_DUSTRIBUTION, MASS_EDENHOFER = 5.5e3, 1.6e3       # Msun

CLOUD = fits.getdata(f"{OUT_DIR}/eos_maps.fits", "CLOUD_H2_10PCT").astype(bool)
NH_ED = fits.getdata(f"{OUT_DIR}/eos_maps.fits", "NH_FULL").astype(float)
L, B = em.LL2[CLOUD], em.BB2[CLOUD]
print(f"inside the 10% contour: {CLOUD.sum()} px; Edenhofer NH_full = integral over "
      f"{em.R_MIN:g}-{em.R_MAX:g} pc under the default conversion")

dens = load_dustribution_density(data_path("dustribution_density"))                            # n_H [cm^-3], (l, b, d)
ext = np.load(EXT_CUBE, mmap_mode="r")                # mag at d_bounds, (l, b, 118)
il, ib, _ = dustribution_index(L, B, 0)
il, ib = np.rint(il).astype(int), np.rint(ib).astype(int)

# (1) Edenhofer NH_full vs the DUSTRIBUTION cumulative extinction product (0-350 pc)
cum_NH = np.asarray(ext[..., -1])[il, ib] * NH_PER_MAGPC * PC_CM
q = NH_ED[CLOUD] / cum_NH
print("\n(1) Edenhofer NH_full(69-260) vs DUSTRIBUTION cumulative column (ext_med_cube 0-350 pc x "
      f"{NH_PER_MAGPC:.0f} cm^-3/(mag/pc) x pc):")
print(f"    medians {np.median(NH_ED[CLOUD]):.3g} vs {np.median(cum_NH):.3g} cm^-2 -> ratio of medians "
      f"{np.median(NH_ED[CLOUD])/np.median(cum_NH):.3f}; per-pixel median {np.median(q):.3f} "
      f"(16-84% {np.percentile(q,16):.3f}-{np.percentile(q,84):.3f})")

# (2) published mass ratio
print(f"\n(2) Burkhart et al. 2025 (quoted): Eos mass DUSTRIBUTION {MASS_DUSTRIBUTION:.2g} vs "
      f"Edenhofer {MASS_EDENHOFER:.2g} Msun -> {MASS_DUSTRIBUTION/MASS_EDENHOFER:.2f}x")

# (3) DUSTRIBUTION median-density integral vs its own cumulative product, native units (mag)
cells = np.unique(np.stack([il, ib], 1), axis=0)
cl, cb = cells.T
dd = DUSTRIBUTION_DMAX / DUSTRIBUTION_SHAPE[2]
integ = (dens[cl, cb] / NH_PER_MAGPC).sum(1) * dd       # back to mag (no gas conversion)
cumm = np.asarray(ext[..., -1])[cl, cb]
r3 = cumm / integ
print(f"\n(3) DUSTRIBUTION cumulative extinction / integral of its median density cube, 0-350 pc, "
      f"native mag, {len(cells)} 1-deg cells:")
print(f"    median {np.median(r3):.3f} (16-84% {np.percentile(r3,16):.2f}-{np.percentile(r3,84):.2f}) "
      "-> use the cumulative product for columns")

# (4) same-geometry check: DUSTRIBUTION density integral over 69-260 pc vs Edenhofer NH_full
r = np.arange(em.R_MIN, em.R_MAX + em.DR / 2, em.DR)
nd = np.concatenate([sample_dustribution(dens, L[s:s+4000], B[s:s+4000], r)
                     for s in range(0, L.size, 4000)])
NH_du = nd.sum(1) * em.DR * PC_CM
q4 = NH_du / NH_ED[CLOUD]
print(f"\n(4) same geometry ({em.R_MIN:g}-{em.R_MAX:g} pc): DUSTRIBUTION median-density integral / "
      f"Edenhofer NH_full: ratio of medians {np.median(NH_du)/np.median(NH_ED[CLOUD]):.3f}, "
      f"per-pixel median {np.median(q4):.3f}")
