# Method: the dust-to-gas conversion chain, line by line

Every number below comes from `cloudfield/conversions.py` and
`cloudfield/columns.py`. They were evaluated with `dust_extinction` 1.5,
numpy 1.21.5 and astropy 5.0.4. Re-derive any of them with:

```python
from cloudfield import conversions as c
c.A1565_AV, c.nh_per_tau("lightcube"), c.nh_per_E("lightcube"), c.nh_per_av("lightcube")
```

## 1. What the dust cube contains

The Edenhofer et al. (2024) map is released in units of **E per pc**. E is the
extinction unit of Zhang, Green & Rix (2023). The cube used here is
**tau_1565 per pc** on the LightCube grid (1 pc voxels, heliocentric Galactic
Cartesian). LightCube (Bish et al. 2026) built it from E/pc by the chain below.
The chain is invertible, so a gas conversion quoted per (E/pc) or per (A_V/pc)
maps to one per (tau_1565/pc).

| Step | Relation | Value | Source |
|---|---|---|---|
| E → A_V | A_V = 2.8 E | `E_TO_AV = 2.8` | Zhang et al. (2023), adopted by Edenhofer et al. (2024) |
| A_V → A(1565 Å) | A_1565 = (A_1565/A_V) A_V | **A_1565/A_V = 2.7032** (2.703165…) | Gordon et al. (2023) "G23" curve, R_V = 3.1, via `dust_extinction` 1.5: `G23(Rv=3.1)(0.1565 µm)` |
| A_1565 → tau_1565 | tau = A / 1.086 | `MAG_PER_TAU = 1.086` | 2.5 log10 e |
| combined | tau_1565 = E × 2.8 × 2.7032 / 1.086 | **tau_1565 = 6.9695 E** | |

The G23 model returns 0.984 at 0.55 µm, so its V normalization is not exactly at
0.55 µm. The model output is used directly as A(1565)/A(V), with no
renormalization. An earlier version of this pipeline assumed A_1565/A_V ≈ 2.6;
G23 gives 4% more.

## 2. Gas density

n_H is the **total** hydrogen-nucleus density, n_HI + 2 n_H2, in cm^-3:

    n_H = (tau_1565 per pc) × K_tau,      K_tau = K_E × 1.086 / (2.8 × A_1565/A_V)
    K_E = K_AV × 2.8,                     K_AV = n_H per (A_V per pc)
    N_H / A_V = K_AV × (3.086e18 cm/pc)   [cm^-2 mag^-1]

| Option | K_AV [cm^-3 per A_V/pc] | K_E [per E/pc] | **K_tau [per tau_1565/pc]** | implied N_H/A_V [cm^-2/mag] | vs Bohlin 1.87e21 | Source |
|---|---|---|---|---|---|---|
| `lightcube` (default) | 812.00 | 2273.60 | **326.222** | 2.506e21 | ×1.340 | Bish et al. 2026 §III.3 |
| `zucker21` | 590.00 | 1652.00 | 237.033 | 1.821e21 | ×0.974 | Zucker+2021 / O'Neill+2024 / McCallum+2025 route |
| `hensley_draine23` | 964.29 | 2700.00 | 387.403 | 2.976e21 | ×1.591 | Hensley & Draine (2023) route (approximate) |
| `sigma1565_legacy` | 640.14 | 1792.40 | 257.178 | 1.976e21 | ×1.056 | σ_1565 = 1.26e-21 cm²/H: Bohlin N_H/A_V = 1.87e21 divided by an assumed A_1565/A_V = 2.6 |

**The Bohlin caveat on `lightcube`.** Bish et al. (2026) attribute 812 cm^-3 per
(A_V/pc) to Bohlin, Savage & Drake (1978). However, 812 × 3.086e18 =
2.51e21 cm^-2 mag^-1, about 34% above the standard Bohlin value of 1.87e21 for
R_V = 3.1. It is adopted because it is internally consistent with the LightCube
cube being sampled, not because it is the Bohlin value.

**This is a single multiplicative scalar on n,** and therefore on alphaG ∝ 1/n.
Across the four options it spans a factor of 1.63 (237 to 387 per tau/pc). It is
the largest arbitrary factor in the pipeline.

## 3. Quality cut (conversion-independent)

A sightline counts as having a cloud if its peak dust density satisfies

    tau_1565,peak ≥ 0.3 / K_tau(lightcube) = 9.1962e-4 per pc   (= 1.3195e-4 E per pc)

The cut is applied in dust units, upstream of any gas conversion, so the
quality mask is bit-identical for every conversion option
(`scripts/conversion_compare.py`). The equivalent n_peak cut is 0.300
(`lightcube`), 0.218 (`zucker21`), 0.356 (`hensley_draine23`) and 0.237
(`sigma1565_legacy`) cm^-3.

## 4. alphaG

    alphaG = 0.59 I_UV (100 cm^-3 / n) [9.9 / (1 + 8.9 σ̃_g)]^0.37
           = 59 I_UV / n          for σ̃_g = 1 (σ_g = 1.9e-21 cm², solar metallicity)

(Sternberg et al. 2014; Bialy & Sternberg 2016). I_UV is the free-space FUV
field in units of the Draine (1978) standard ISRF. n is the total
hydrogen-nucleus density.

* **I_UV.** The LightCube combined 912–2000 Å cube, used as is. It is already
  normalized to the Draine field (header BUNIT `'Draine78 ISRF'`). Check: the mean
  over voxels within 25 pc of the Sun is 0.6307, against 0.633 in Bish et al.
  (2026) §III.1.
* **No rescaling is applied to I_UV.**

## 5. Geometry

* **Rays.** Rays are cast from the Sun through each (l, b) pixel centre.
  Samples run every `dr` (0.25 pc) from r_min to r_max, at positions
  x = r cos b cos l, y = r cos b sin l, z = r sin b.
* **Interpolation.** Both cubes are trilinearly interpolated (order 1, constant
  0 outside), after a padded crop.
* **Origin.** The Sun is at 0-based voxel index 499.5 on every axis, from the
  LightCube header (CRVAL = −500 at CRPIX = 0.5, CDELT = 1).
* **Inner masking.** r < 69 pc is masked: the Edenhofer map has no dust there.
* **Boundaries.** r0, r_peak and r1 come from the configured boundary method
  (docs/BOUNDARIES.md). G0_surface = I_UV(r0) and G0_peak = I_UV(r_peak).
* **Columns.**
  * NH_cloud = ∫_{r0}^{r1} n dr.
  * NH_full = ∫ n dr over the whole radial grid.

## 6. DUSTRIBUTION (column comparison only)

`cloudfield/columns.py` follows Burkhart et al. (2025) and Draine (2003):

| Quantity | Value |
|---|---|
| κ_abs(547 nm) | 8.551e3 cm²/g |
| albedo(547 nm) | 0.6735 |
| κ_ext = κ_abs/(1−albedo) | 26189.9 cm²/g |
| gas/dust mass ratio | 124 (notebook: 123.6) |
| mean gas mass per H nucleus | 1.4 m_H (m_H = 1.6735e-24 g) |
| n_H per (mag/pc) = g2d / (1.086 κ_ext pc × 1.4 m_H) | **602.99 cm^-3** |

DUSTRIBUTION is never used for G0 or alphaG maps.
