"""Cube access, cropping, line-of-sight sampling and the origin convention.

Cubes: LightCube grid, heliocentric Galactic Cartesian (x toward l=0, y toward
l=90, z toward b=90), axis order [z, y, x], 1 pc voxels. Both the ISRF cube
and the Edenhofer tau_1565 cube must be on the same grid.
"""
import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates

from .constants import ISRF_BUNIT, SOLAR_I_UV, SOLAR_RADIUS_PC


class CubeError(RuntimeError):
    pass


def open_cubes(isrf_path, dust_path):
    """Memmapped (isrf, dust, isrf_header). Nothing is read yet."""
    hi = fits.open(isrf_path, memmap=True)
    hd = fits.open(dust_path, memmap=True)
    isrf, dust = hi[0].data, hd[0].data
    if isrf.shape != dust.shape:
        raise CubeError(f"Cube shape mismatch: ISRF {isrf.shape} vs dust {dust.shape}")
    return isrf, dust, hi[0].header


def open_isrf(isrf_path):
    hi = fits.open(isrf_path, memmap=True)
    return hi[0].data, hi[0].header


def origin_index(shape, header, mode):
    """0-based [z, y, x] voxel index of the Sun.
    'n//2'  : shape//2
    'header': from the ISRF WCS (CRVAL=-500 at CRPIX=0.5 -> 499.5)."""
    if mode == "n//2":
        return np.array([s // 2 for s in shape], dtype=float)
    if mode == "header":
        return np.array([(0.0 - header[f"CRVAL{k}"]) / header[f"CDELT{k}"]
                         + header[f"CRPIX{k}"] - 1.0 for k in (3, 2, 1)])
    raise CubeError(f"Unknown origin mode {mode}")


def cube_radius_pc(shape, origin):
    """Largest radius (pc) fully inside the cube box (for clipped cubes this is
    an upper bound; a spherical clip is reported by the caller)."""
    return float(min(min(o, s - 1 - o) for o, s in zip(origin, shape)))


def los_xyz(l_deg, b_deg, r):
    """Heliocentric Galactic Cartesian (pc) for rays: l_deg, b_deg shape [P], r [R]
    -> X, Y, Z each [P, R]."""
    lr = np.radians(np.asarray(l_deg, float))[:, None]
    br = np.radians(np.asarray(b_deg, float))[:, None]
    X = r[None, :] * np.cos(br) * np.cos(lr)
    Y = r[None, :] * np.cos(br) * np.sin(lr)
    Z = r[None, :] * np.sin(br)
    return X, Y, Z


def sample_los(isrf_mm, dust_mm, origin, X, Y, Z, pad=2, verbose=None):
    """Crop both cubes to the padded integer bbox of the sample points, read as
    float32 + nan_to_num, and trilinearly sample (order=1, constant, cval=0).
    Returns (tau_pc, isrf) as float64 arrays shaped like X; tau_pc is None
    if dust_mm is None (ISRF only).
    verbose: None, or a callable(name, array) for range printing."""
    shape = isrf_mm.shape
    g = (Z + origin[0], Y + origin[1], X + origin[2])
    lo = [int(np.floor(a.min())) - pad for a in g]
    hi = [int(np.ceil(a.max())) + pad + 1 for a in g]
    for ax, a, b, n in zip("zyx", lo, hi, shape):
        if a < 0 or b > n:
            raise CubeError(f"Crop on {ax} ({a}:{b}) leaves the cube (0:{n}): the requested "
                            "rays extend beyond the cube volume.")
    sl = tuple(slice(a, b) for a, b in zip(lo, hi))
    isrf_c = np.nan_to_num(np.asarray(isrf_mm[sl], dtype=np.float32))
    dust_c = (None if dust_mm is None else
              np.nan_to_num(np.asarray(dust_mm[sl], dtype=np.float32)))
    if verbose is not None:
        print(f"  crop [z,y,x] = [{lo[0]}:{hi[0]}, {lo[1]}:{hi[1]}, {lo[2]}:{hi[2]}]")
        verbose("ISRF crop", isrf_c)
        if dust_c is not None:
            verbose("dust crop", dust_c)
    coords = np.stack([g[0] - lo[0], g[1] - lo[1], g[2] - lo[2]]).reshape(3, -1)
    samp = lambda c: map_coordinates(c, coords, order=1, mode="constant",
                                     cval=0.0).reshape(X.shape).astype(np.float64)
    return (None if dust_c is None else samp(dust_c)), samp(isrf_c)


def solar_mean(isrf_mm, header, radius=SOLAR_RADIUS_PC, origin_mode="header"):
    """Mean I_UV over voxels whose centres lie within `radius` pc of the origin
    (Bish+2026 Sec III.1 convention). Expected 0.633 for the combined cube."""
    o = origin_index(isrf_mm.shape, header, origin_mode)
    k = int(np.ceil(radius)) + 1
    c = np.floor(o).astype(int)
    sl = tuple(slice(ci - k, ci + k + 2) for ci in c)
    blk = np.asarray(isrf_mm[sl], dtype=np.float64)
    zz, yy, xx = [np.arange(s.start, s.stop) - oi for s, oi in zip(sl, o)]
    R = np.sqrt(zz[:, None, None] ** 2 + yy[None, :, None] ** 2 + xx[None, None, :] ** 2)
    v = blk[R <= radius]
    return float(v.mean()), int(v.size)


def check_isrf_units(header, solar, tol=0.05):
    """Raise unless the cube is the Draine-normalized combined cube."""
    bunit = str(header.get("BUNIT", "")).strip()
    if bunit != ISRF_BUNIT:
        raise CubeError(f"ISRF BUNIT is {bunit!r}, expected {ISRF_BUNIT!r} (the combined "
                        "912-2000 A Draine-normalized LightCube cube). A per-band cube gives "
                        "alphaG wrong by orders of magnitude; see data/README.md.")
    if abs(solar / SOLAR_I_UV - 1) > tol:
        raise CubeError(f"Solar-neighbourhood mean I_UV = {solar:.4g}, expected {SOLAR_I_UV} "
                        f"(+-{tol:.0%}). Wrong cube? See data/README.md.")


def voxel_snr(isrf_mm, origin, center_xyz_pc, half=20):
    """Per-voxel signal-to-noise of the ISRF cube from the robust (MAD) scatter
    of relative first differences in a (2*half)^3 pc block centred on
    center_xyz_pc (heliocentric x, y, z). Real gradients bias this low (conservative)."""
    zc, yc, xc = (np.array([center_xyz_pc[2], center_xyz_pc[1], center_xyz_pc[0]])
                  + origin).round().astype(int)
    hb = int(half)
    blk_ = np.asarray(isrf_mm[zc-hb:zc+hb, yc-hb:yc+hb, xc-hb:xc+hb], dtype=np.float64)
    rel = np.concatenate([(np.diff(blk_, axis=a) /
                           (0.5 * (blk_[(slice(None),)*a + (slice(1, None),)] +
                                   blk_[(slice(None),)*a + (slice(None, -1),)]))).ravel()
                          for a in range(3)])
    return 1.0 / (1.4826 * np.median(np.abs(rel - np.median(rel))) / np.sqrt(2))


def rng(name, a):
    a = np.asarray(a, dtype=float)
    f = a[np.isfinite(a)]
    if f.size == 0:
        print(f"  {name:<16s} shape={a.shape}  all non-finite")
    else:
        print(f"  {name:<16s} shape={a.shape}  min={f.min():.4g}  "
              f"max={f.max():.4g}  median={np.median(f):.4g}  "
              f"nonfinite={a.size - f.size}")
