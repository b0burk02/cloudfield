"""Sky grid (Galactic plate carree), WCS, radial sampling, and resolution report.

Pixel scale modes (cloud config `grid.pixel`):
  {mode: fixed,    deg: 0.1}                       fixed angular pixel
  {mode: physical, pc: 0.5, ref_distance_pc: 120}  pixel = pc / ref_distance (rad -> deg)

The dust reconstruction has a finite effective resolution. At distance d it is
max(cube voxel size, d x native angular resolution of the dust map). A pixel
smaller than that oversamples a smoothed reconstruction: neighbouring pixels are
not independent. The oversampling factor is reported and written to the FITS
headers, and a warning is printed when it exceeds 1.
"""
import warnings
from dataclasses import dataclass, field
import numpy as np
from astropy.wcs import WCS


class OversamplingWarning(UserWarning):
    pass


def pixel_deg_from_config(pix_cfg):
    mode = pix_cfg.get("mode", "fixed")
    if mode == "fixed":
        return float(pix_cfg["deg"])
    if mode == "physical":
        return float(np.degrees(pix_cfg["pc"] / pix_cfg["ref_distance_pc"]))
    raise ValueError(f"unknown grid.pixel.mode {mode!r} (fixed | physical)")


@dataclass
class SkyGrid:
    l_range: tuple
    b_range: tuple
    pix: float                                   # deg
    NL: int = field(init=False)
    NB: int = field(init=False)

    def __post_init__(self):
        L_RANGE, B_RANGE, PIX = self.l_range, self.b_range, self.pix
        self.NL = int(round((L_RANGE[1] - L_RANGE[0]) / PIX)) + 1
        self.NB = int(round((B_RANGE[1] - B_RANGE[0]) / PIX)) + 1
        self.lvals = L_RANGE[1] - PIX * np.arange(self.NL)     # descending: l increases to the left
        self.bvals = B_RANGE[0] + PIX * np.arange(self.NB)
        self.LL2, self.BB2 = np.meshgrid(self.lvals, self.bvals)

    @property
    def shape(self):
        return (self.NB, self.NL)

    def wcs(self):
        L_RANGE, B_RANGE, PIX = self.l_range, self.b_range, self.pix
        w = WCS(naxis=2)
        w.wcs.ctype = ["GLON-CAR", "GLAT-CAR"]
        w.wcs.cunit = ["deg", "deg"]
        w.wcs.cdelt = [-PIX, PIX]
        w.wcs.crval = [0.5 * (L_RANGE[0] + L_RANGE[1]), 0.0]   # CRVAL2=0 -> pixels uniform in b
        w.wcs.crpix = [1 + (L_RANGE[1] - w.wcs.crval[0]) / PIX, 1 - B_RANGE[0] / PIX]
        return w

    def check_wcs(self):
        w = self.wcs()
        wl, wb = w.wcs_pix2world([[0, 0], [self.NL - 1, self.NB - 1]], 0).T
        if not (np.allclose(wl, [self.lvals[0], self.lvals[-1]])
                and np.allclose(wb, [self.bvals[0], self.bvals[-1]])):
            raise RuntimeError(f"WCS does not reproduce the pixel grid: {wl}, {wb}")


def radial_grid(r_min, r_max, dr):
    return np.arange(r_min, r_max + dr / 2, dr)


@dataclass
class Resolution:
    distance_pc: float
    pixel_deg: float
    pixel_pc: float
    effective_pc: float
    effective_deg: float
    oversampling: float

    def header_cards(self):
        return [("RESDIST", round(self.distance_pc, 2), "pc; distance resolution is quoted at"),
                ("PIXPC", round(self.pixel_pc, 4), "pc per pixel at RESDIST"),
                ("EFFRESPC", round(self.effective_pc, 4), "pc; dust-map effective resolution at RESDIST"),
                ("OVERSAMP", round(self.oversampling, 3), "effective res / pixel (>1: oversampled)")]


def resolution(pixel_deg, distance_pc, voxel_pc, native_ang_deg=None, warn=True):
    pixel_pc = distance_pc * np.radians(pixel_deg)
    eff = voxel_pc if native_ang_deg is None else max(voxel_pc, distance_pc * np.radians(native_ang_deg))
    res = Resolution(distance_pc, pixel_deg, pixel_pc, eff, float(np.degrees(eff / distance_pc)),
                     eff / pixel_pc)
    if warn and res.oversampling > 1.0:
        warnings.warn(f"pixel {pixel_deg:g} deg = {pixel_pc:.3g} pc at {distance_pc:g} pc is finer "
                      f"than the dust map's effective resolution {eff:.3g} pc "
                      f"({res.effective_deg:.3g} deg): oversampled x{res.oversampling:.1f}. "
                      "Neighbouring pixels are not independent.", OversamplingWarning)
    return res


def pct(a, mask=None):
    f = a[np.isfinite(a) & (True if mask is None else mask)]
    return (np.percentile(f, [50, 5, 95]) if f.size else np.full(3, np.nan)), f.size
