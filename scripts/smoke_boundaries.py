"""Smoke test of every boundary method on the Eos field (real data, no mocks).

- half_max / fixed_range / nearest_peak: full run, statistics + flag fractions.
- mask3d: (1) the loader must raise NotImplementedError (no mask file yet);
          (2) the segment logic is exercised with a Mask3D built from the
              Edenhofer dust itself (isosurface tau_1565/pc >= threshold, no
              file format involved). This tests the contract, not a product.
Writes no FITS.
"""
import copy
import numpy as np
from cloudfield.config import load_cloud
from cloudfield.maps import Run
from cloudfield.boundaries import (Mask3D, Mask3DBoundary, load_mask3d, FLAG_NAMES,
                                   METHOD_CODES)
from cloudfield.grid import pct

base = load_cloud("eos")


def run_with(boundary, label, mask_obj=None):
    cfg = copy.deepcopy(base)
    cfg["boundary"] = boundary
    R = Run(cfg, verbose=False).open()
    if mask_obj is not None:
        R.boundary = Mask3DBoundary(mask_obj(R), min_chord_pc=2.0, segment="peak")
    R.build_region()
    out = R.compute()
    M, reg, Qn = out["M"], R.region.mask, out["Qn"]
    print(f"--- {label}: passing {Qn.sum()} px, region {Qn[reg].mean():.1%}")
    for k in ("r0", "rpeak", "r1", "G0_surface", "alphaG_peak"):
        (m, a, b), nf = pct(M[k], reg)
        print(f"    {k:<12} median {m:8.4g}  5-95% {a:8.4g}-{b:8.4g}  N={nf}")
    fl = out["FLAGS"][reg]
    print("    flags (all region px): " + ", ".join(
        f"{n} {np.mean((fl & b) > 0):.1%}" for b, n in sorted(FLAG_NAMES.items()) if np.any(fl & b)))
    mc = out["METHOD"][reg]
    print("    method per pixel: " + ", ".join(f"{n} {np.mean(mc == c):.1%}"
                                              for n, c in METHOD_CODES.items() if np.any(mc == c)))
    d = np.abs(np.diff(M["r0"], axis=1))
    d = d[(reg[:, 1:] & reg[:, :-1]) & np.isfinite(d)]
    print(f"    adjacent |dr0| > 5 pc: {np.mean(d > 5):.3%}, max {d.max():.1f} pc")
    return out


run_with({"method": "half_max", "window": [90.0, 145.0]}, "half_max (Eos default)")
run_with({"method": "fixed_range", "r0_pc": 95.0, "r1_pc": 140.0}, "fixed_range 95-140 pc")
run_with({"method": "nearest_peak", "window": [90.0, 145.0], "target_pc": 110.0, "min_frac": 0.5},
         "nearest_peak target 110 pc")

print("--- mask3d loader")
try:
    load_mask3d({"path": "anything.fits"})
    print("    ERROR: loader did not raise")
except NotImplementedError as e:
    print(f"    raises NotImplementedError as intended: {str(e)[:90]}...")


class DustIsosurface(Mask3D):
    """TEST ONLY: cloud = Edenhofer tau_1565/pc >= thr, defined for 69 <= r <= 200 pc
    (narrower than the 69-260 pc sampling, to exercise coverage truncation)."""
    voxel_pc = 1.0

    def __init__(self, run, thr=2.0 * 9.196e-4):
        self.d, self.o, self.thr = run.dust_mm, run.origin, thr

    def _idx(self, x, y, z):
        return tuple(np.clip(np.rint(np.asarray(v) + o).astype(int), 0, s - 1)
                     for v, o, s in zip((z, y, x), self.o, self.d.shape))

    def contains(self, x, y, z):
        iz, iy, ix = self._idx(x, y, z)
        lo = [int(a.min()) for a in (iz, iy, ix)]
        hi = [int(a.max()) + 1 for a in (iz, iy, ix)]
        sub = np.asarray(self.d[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]])
        return sub[iz - lo[0], iy - lo[1], ix - lo[2]] >= self.thr

    def defined(self, x, y, z):
        r = np.sqrt(np.asarray(x) ** 2 + np.asarray(y) ** 2 + np.asarray(z) ** 2)
        return (r >= 69.0) & (r <= 200.0)


run_with({"method": "half_max", "window": [90.0, 145.0]},
         "mask3d segment logic on a dust isosurface (TEST ONLY)", mask_obj=DustIsosurface)
