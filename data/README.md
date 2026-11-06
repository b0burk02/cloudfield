# Input data (none of it is in this repository)

This repository contains **no data products**. Every input is publicly available
from its authors, and you must download it yourself. The code finds files through
a config file, never through hardcoded paths:

```bash
cp config.example.yaml config.yaml        # then edit the paths; config.yaml is gitignored
export CLOUDFIELD_CONFIG=/path/to/config.yaml   # optional; default ./config.yaml
export CLOUDFIELD_DATA=/path/to/data_dir        # optional; base for relative paths
cloudfield check-data                     # ALWAYS run this first
```

If a required file is missing, every command stops with an error that names the
config key and the path, and points back to this file.

| Key | Needed for | If missing |
|---|---|---|
| `lightcube` | `cloudfield run` (required) | run aborts |
| `edenhofer_tau` | `cloudfield run` (required) | run aborts |
| `fims_h2` | only clouds whose `region` uses it (Eos) | that cloud's run aborts; other clouds unaffected |
| `dustribution_density`, `dustribution_cumulative` | `scripts/column_systematics.py` only | that script aborts |
| `galfa_nhi`, `galfa_cube`, `galfa_check` | `scripts/hi_check.py` only | that script aborts |

---

## LightCube (Bish et al. 2026): I_UV  [required]

* Paper: Bish, Peek, Zucker, Gordon, Murray, Clark, Willems & Hamden 2026,
  arXiv:2608.27389.
* Download: https://hannahbish.me/LightCubes, also via Globus (origin id
  `f2d103b7-ee4a-4f0b-bf1f-dc7e6df06b6f`).
* **You need the COMBINED 912–2000 Å cube, normalized to the Draine (1978)
  Standard ISRF** (header `BUNIT = 'Draine78 ISRF'`). Do **not** use one of the
  single TD1-band cubes (1565, 1965, 2365 or 2740 Å; erg cm^-3 Å^-1) or the
  uncertainty cubes.
  * Pointing the code at a band cube makes alphaG wrong by orders of magnitude
    with nothing else raising an error.
  * `cloudfield check-data` and every run reject it, because they check BUNIT and
    the solar value.
* **Check on download.** The mean over voxels within 25 pc of the cube origin must
  be **0.633** (Bish et al. 2026 §III.1). With the r = 500 pc product used during
  development we measure 0.6307.
* **Resolution.** Released on a 1 pc and a 4 pc grid. The underlying model has
  1 pc voxels only within 200 pc of the Sun: 2 pc at 200–400 pc, 4 pc at
  400–800 pc, 8 pc beyond. So the 1 pc grid adds no information beyond 200 pc.
* **File and size.** Development used `ISRF_StandardUV_1pc-res_clipped-r500.fits`:
  1000³ float32, **4.0 GB**, values clipped to a sphere of r = 500 pc. Set
  `cube_radius_pc: 500` in config.yaml when using a clipped cube; runs that would
  sample beyond it then raise instead of silently reading zeros.
* **Full ±1.25 kpc sizes (estimated from the grid, not measured):**
  * 1 pc grid: 2500³ voxels ≈ **62.5 GB per cube** (float32).
  * 4 pc grid: 625³ ≈ **1.0 GB**.

## Edenhofer et al. (2024) 3D dust: n_H  [required]

* Paper: Edenhofer et al. 2024, A&A 685, A82. Data: Zenodo (see the paper or the
  `dustmaps` documentation).
* **Native release.** E per pc (Zhang, Green & Rix 2023 units) on a spherical grid:
  * HEALPix Nside 256 (~14′) × logarithmic distance bins (0.4 pc at 69 pc to
    ~7 pc at 1.25 kpc);
  * a mean map, a standard-deviation map and posterior samples (12 reported);
  * a lower-fidelity extension out to 2 kpc.
* **What this code expects.** The dust as **tau_1565 per pc on the LightCube
  Cartesian grid**: same shape, origin and voxel size as the LightCube cube. This
  is the form used to build LightCube.
  * Development used `Edenhofer_tau_per_pc_1565_1pc-res_clipped-r500.fits`:
    1000³ float64, **8.0 GB**, r ≤ 500 pc.
  * The full ±1.25 kpc cube on the 1 pc grid would be ≈ 125 GB as float64
    (62.5 GB as float32), estimated.
* **This Cartesian tau_1565 product is not part of the Edenhofer Zenodo release.**
  Obtain it from the LightCube authors, or build it from the native release:
  1. resample E/pc to the LightCube grid;
  2. multiply by 2.8 × A_1565/A_V / 1.086 = 6.9695 (docs/METHOD.md).

  The code does **not** do this resampling.
* Which Edenhofer statistic (mean, or a sample) the LightCube-grid cube was built
  from is not recorded in its header. Confirm with the LightCube authors.
* The map has no dust inside 69 pc (the cells are zero, not missing). The code
  masks r < 69 pc.

## DUSTRIBUTION (Dharmawardena et al.): column comparison only  [optional]

* Source: www.mwdust.com (3D dust density cube, solar neighbourhood to 350 pc).
  Grid (l, b, d) = (360, 180, 117): uniform in l, sin b and d (1° × 1° × 3 pc),
  extinction density in mag/pc.
* `dustribution_density`: the median density cube. Development used
  `gpy_dens_median.out`, a torch zip of 7,581,600 float64 values (61 MB). It is
  read without importing torch.
* `dustribution_cumulative`: the cumulative extinction cube, (360, 180, 118) in
  mag at the distance bounds. Development used `ext_med_cube.pkl.npy` (61 MB).
* **For columns, use the CUMULATIVE product,** not an integral of the median
  density cube. The integral comes out ~2.4× lower inside clouds; see the README,
  Limitations.
* Never used for G0 or alphaG.

## GALFA-HI and FIMS/SPEAR: Eos validation only  [optional]

These are needed only to reproduce the Eos checks (`scripts/hi_check.py`, the Eos
region contour), not to run the pipeline on other clouds.

| Key | File used | Size | Notes |
|---|---|---|---|
| `galfa_nhi` | `GALFA-HI_NHISRCORR_VLSR-90+90kms.fits` | 420 MB | GALFA-HI DR2 (Peek et al. 2018), stray-corrected N_HI, −90..+90 km/s, RA/Dec |
| `galfa_cube` | `GALFAmosaic_newhd.GAL.fits` | 642 MB | PPV mosaic, GLON/GLAT-TAN, −100..+100 km/s |
| `galfa_check` | `GALFAmosaic_newhd.GAL.HI_col.fits` | 1.8 MB | Galactic-frame column map (frame check) |
| `fims_h2` | `f03c.fits` | 1 MB | FIMS/SPEAR H2/FUV ratio (%), HEALPix Nside 64, RING, COORDSYS G. The Eos boundary is its 10% contour (Burkhart et al. 2025) |

## Practical note on cube extent

The development cubes are clipped at r = 500 pc. Rays that stay inside that
sphere work (e.g. Eos, 69–260 pc; Taurus and Ophiuchus, ≲180 pc). Anything
reaching beyond needs the unclipped LightCube and Edenhofer cubes. On disk:
* about 1 GB each for the 4 pc products;
* about 60–125 GB for the 1 pc products.

The code reads them by memory map and crops before reading, so RAM is not the
limit.
