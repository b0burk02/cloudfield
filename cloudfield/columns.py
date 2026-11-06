"""DUSTRIBUTION (Dharmawardena et al.; the 3D map used by Burkhart et al. 2025).

COLUMN-DENSITY COMPARISON ONLY. Never used for G0, alphaG or any map:
LightCube's I_UV was computed through the Edenhofer map, so n must come from
Edenhofer too. cloudfield.maps raises if asked to use any other dust source.

For columns, prefer the CUMULATIVE extinction product: integrating the median
density cube comes out ~2.4x lower inside clouds (median of integrals vs
integral of medians); see README Limitations.
"""
import zipfile
import numpy as np
from scipy.ndimage import map_coordinates

from .constants import PC_CM

DUSTRIBUTION_SHAPE = (360, 180, 117)    # (n_l, n_b, n_d), row-major
DUSTRIBUTION_DMAX = 350.0               # pc; d bins linear 0-350
# Dust properties at 547 nm (Draine 2003), as in Burkhart et al. (2025)
KAPPA_ABS_547 = 8.551e3      # cm^2 per g of dust
ALBEDO_547 = 0.6735
GAS_TO_DUST = 124.0          # gas/dust mass ratio (Burkhart+25 notebook: 123.6)
MU_H = 1.4                   # gas mass per H nucleus, in m_H
M_H = 1.6735e-24             # g
MAG_PER_TAU = 1.086
# n_H [cm^-3] = rho_ext [mag/pc] / (1.086 kappa_ext pc) * g2d / (mu m_H)
NH_PER_MAGPC = GAS_TO_DUST / (MAG_PER_TAU * KAPPA_ABS_547 / (1 - ALBEDO_547)
                              * PC_CM * MU_H * M_H)


def load_dustribution_density(path, verbose=True):
    """Median density cube -> n_H [cm^-3], shape (n_l, n_b, n_d). Read from the
    torch zip's raw float64 storage without importing torch."""
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith("/data/0")]
        if len(names) != 1:
            raise RuntimeError(f"Unexpected DUSTRIBUTION archive layout: {z.namelist()}")
        raw = np.frombuffer(z.read(names[0]), dtype="<f8")
    if raw.size != np.prod(DUSTRIBUTION_SHAPE):
        raise RuntimeError(f"DUSTRIBUTION size {raw.size} != {DUSTRIBUTION_SHAPE}")
    dens = raw.reshape(DUSTRIBUTION_SHAPE)
    if verbose:
        print(f"  DUSTRIBUTION {path}: shape {dens.shape}, {dens.min():.3g}-{dens.max():.3g} mag/pc; "
              f"n_H = rho x {NH_PER_MAGPC:.4g} cm^-3 per (mag/pc)")
    return dens * NH_PER_MAGPC


def load_dustribution_cumulative(path):
    """Cumulative extinction [mag] at the 118 distance bounds (0..350 pc), (l, b, 118)."""
    return np.load(path, mmap_mode="r")


def dustribution_index(l_deg, b_deg, r):
    """Fractional (il, ib, id) cell-centre indices. Grid: uniform in l (1 deg),
    uniform in sin(b) (180 bins over [-1, 1]), uniform in d (350/117 pc)."""
    nl, nb, nd = DUSTRIBUTION_SHAPE
    il = np.mod(l_deg, 360.0) * nl / 360.0 - 0.5
    ib = (np.sin(np.radians(b_deg)) + 1.0) * nb / 2.0 - 0.5
    idd = r * nd / DUSTRIBUTION_DMAX - 0.5
    return il, ib, idd


def sample_dustribution(cube, l_deg, b_deg, r):
    """Trilinear n_H [cm^-3] along rays: l_deg, b_deg [P], r [R] -> [P, R]."""
    if r.max() > DUSTRIBUTION_DMAX:
        raise ValueError(f"r up to {r.max()} pc exceeds DUSTRIBUTION extent {DUSTRIBUTION_DMAX} pc")
    il, ib, _ = dustribution_index(np.asarray(l_deg)[:, None], np.asarray(b_deg)[:, None], 0)
    _, _, idd = dustribution_index(0, 0, r[None, :])
    il = np.broadcast_to(il, (il.shape[0], r.size))
    ib = np.broadcast_to(ib, (ib.shape[0], r.size))
    idd = np.broadcast_to(idd, il.shape)
    lo = [max(int(np.floor(a.min())) - 2, 0) for a in (il, ib, idd)]
    hi = [min(int(np.ceil(a.max())) + 3, n) for a, n in zip((il, ib, idd), cube.shape)]
    sub = cube[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
    coords = np.stack([il - lo[0], ib - lo[1], idd - lo[2]]).reshape(3, -1)
    return map_coordinates(sub, coords, order=1, mode="nearest").reshape(il.shape)
