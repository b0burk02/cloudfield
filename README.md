# cloudfield

Maps of the incident far-UV field **G0** and the H I-to-H2 transition parameter
**αG** for a molecular cloud with a defined 3D volume. The inputs are the
LightCube 3D model of the local FUV interstellar radiation field (Bish et al.
2026) and the Edenhofer et al. (2024) 3D dust map. Each cloud is described by a
YAML config file. Eos (Burkhart et al. 2025) is the worked example.

> **No data products are in this repository.** Download them as described in
> [data/README.md](data/README.md), then run `cloudfield check-data` first.

**Contents:**
[1. Quantities](#1-what-is-computed) ·
[2. Inputs](#2-inputs-and-why-both-come-from-edenhofer) ·
[3. Method](#3-method-step-by-step) ·
[4. Dust-to-gas conversion](#4-dust-to-gas-conversion) ·
[5. Eos](#5-worked-example-eos) ·
[6. Limitations](#6-limitations) ·
[7. Verification](#7-verification-against-bish-et-al-2026) ·
[8. Running](#8-installation-and-running) ·
[9. New clouds](#9-adding-a-new-cloud) ·
[Citing](#citing)

---

## 1. What is computed

**G0** is the far-UV (912–2000 Å) radiation field incident on the cloud, written
I_UV and measured in units of the **Draine (1978) Standard ISRF** (I_UV = 1 is
the Draine field).

**αG** sets the structure of the H I-to-H2 transition in the analytic theory of
Sternberg et al. (2014). We use Bialy & Sternberg (2016), **Eq. 22**, for H2
formation on dust grains:

```math
\alpha G \;=\; 0.59\; I_{\rm UV} \left(\frac{100\ {\rm cm^{-3}}}{n}\right)\left(\frac{9.9}{1+8.9\,\tilde\sigma_g}\right)^{0.37}
\;=\; \frac{59\, I_{\rm UV}}{n/{\rm cm^{-3}}}\qquad (\tilde\sigma_g = 1)
```

* σ̃_g = σ_g / 1.9×10⁻²¹ cm² is the Lyman-Werner dust absorption cross-section per
  hydrogen nucleus, relative to the standard Galactic value (Bialy & Sternberg
  2016, Eqs. 5–6). We take σ̃_g = 1: σ_g = 1.9×10⁻²¹ cm², solar metallicity.
* **n is the TOTAL hydrogen-nucleus density, n = n_HI + 2 n_H2**, not the atomic
  density.
* The total H I column of the transition layer is
  N_HI,tot = σ_g⁻¹ ln(αG/2 + 1) (Bialy & Sternberg 2016, Eq. 18).
  * αG ≫ 1: H I self-shielding dominates (a thick H I layer).
  * αG ≲ 1: dust shielding dominates.

For every pixel of a Galactic (l, b) grid the pipeline writes:

| Map (FITS extension) | Meaning |
|---|---|
| `G0_SURFACE` | I_UV at the cloud's near-side boundary r0: the field incident on the surface |
| `G0_PEAK` | I_UV at the density peak r_peak |
| `N_SURFACE`, `N_PEAK` | n at r0 and at r_peak [cm⁻³] |
| `ALPHAG_PEAK` (primary), `ALPHAG_SURFACE` | 59 I_UV/n evaluated at r_peak and at r0 |
| `NH_CLOUD` | ∫ n dr from r0 to r1 [cm⁻²] |
| `NH_FULL` | ∫ n dr over the whole sampled ray [cm⁻²] |
| `R0`, `RPEAK`, `R1`, `THICKNESS` | boundary radii [pc] |
| `ALPHAG_LOSMEAN`, `RATIO_SURF_LOSMEAN` | the old slab estimator, kept for comparison only |
| `QUAL_NPEAK`, `QUAL_RPEAK_INTERIOR`, region | quality masks, and the pixels statistics are quoted over |
| `BOUNDARY_FLAGS`, `BOUNDARY_METHOD` | per-pixel diagnostics of the boundary finder |

## 2. Inputs, and why both come from Edenhofer

| Quantity | Source | Form used |
|---|---|---|
| I_UV | **LightCube** (Bish et al. 2026) | combined 912–2000 Å cube, already Draine-normalized (`BUNIT = 'Draine78 ISRF'`), 1 pc voxels |
| n | **Edenhofer et al. (2024)** 3D dust | τ₁₅₆₅ per pc on the LightCube grid, converted to n (Section 4) |

* **Both I_UV and n come from Edenhofer dust.** LightCube's radiative transfer was
  computed *through the Edenhofer map*, so the field it predicts at a point
  already includes extinction by exactly that dust. Pairing LightCube's I_UV with
  a density from another dust map would put the cloud surface, and the
  attenuation the field has experienced, in inconsistent places. The code raises
  if asked to make maps from any other dust source.
  * DUSTRIBUTION (Dharmawardena et al.) is supported only for column-density
    comparison.
* **No rescaling of I_UV.** The LightCube combined cube is already in Draine units.
  To check a new file: the mean over voxels within 25 pc of the cube origin must
  be **0.633** (Bish et al. 2026 §III.1); we measure 0.6307.
  * Using one of the single-band TD1 cubes instead (erg cm⁻³ Å⁻¹) would make αG
    wrong by orders of magnitude, with nothing else raising an error.
    `cloudfield check-data` and every run check BUNIT and the solar value to
    catch this.

## 3. Method, step by step

1. **Ray casting.** Every pixel centre (l, b) defines a ray from the Sun. The cubes
   are heliocentric Galactic Cartesian, and a point at distance r on the ray is

   ```math
   x = r\cos b\cos l,\qquad y = r\cos b\sin l,\qquad z = r\sin b .
   ```

   The Sun is at voxel index 499.5 on each axis (from the LightCube header).
2. **Radial sampling.** Each ray is sampled every Δr (0.25 pc for Eos) from r_min
   to r_max. Both cubes are interpolated trilinearly. r < 69 pc is masked, because
   the Edenhofer map contains no dust there.
3. **Geometry is decided in dust units.** The peak and both boundaries are found
   on the dust profile τ₁₅₆₅(r), never on n. So is the quality cut: a sightline
   has a cloud if τ₁₅₆₅,peak ≥ 9.196×10⁻⁴ pc⁻¹, i.e. n_peak = 0.3 cm⁻³ under the
   default conversion. The gas conversion therefore cannot move a boundary or
   change which pixels pass.
4. **Boundaries r0, r_peak, r1.** The default method for a single clean cloud is
   `half_max`:
   * r_peak = argmax τ(r) within a distance window (90–145 pc for Eos);
   * r0 = the largest r < r_peak with τ < τ_peak/2 (the near-side half-maximum
     crossing);
   * r1 = the smallest r > r_peak with τ < τ_peak/2 (the far side).

   r0 and r1 may lie outside the search window. Other methods are `fixed_range`,
   `nearest_peak` (multi-component profiles) and `mask3d` (an externally supplied
   3D cloud volume). Every pixel records which method produced its boundary, plus
   flags for multi-peaked profiles and boundaries on a search edge. See
   [docs/BOUNDARIES.md](docs/BOUNDARIES.md).
5. **Read the field AT the surface.** G0_surface = I_UV(r0) and
   αG_surface = 59 I_UV(r0)/n(r0). αG describes the field *arriving at* the cloud.
   A density-weighted average of I_UV along the sightline (the earlier slab
   estimator, 59 ⟨I_UV⟩ₙ/⟨n⟩ₙ) mixes in field values from gas in front of and
   behind the cloud. For Eos the two differ by a factor of 1.24 (median), and the
   difference correlates with cloud thickness (Spearman ρ = 0.70), so it is a
   structured bias, not noise. The slab estimator is kept only as
   `ALPHAG_LOSMEAN`.
6. **Peak values.** G0_peak = I_UV(r_peak) and αG_peak = 59 I_UV(r_peak)/n_peak.
   αG_peak is the primary map.
7. **Attenuation test.** G0_surface/G0_peak is compared with two things: the
   per-voxel noise of LightCube (a robust estimate from relative first
   differences), and the ratio expected from the field's large-scale radial
   gradient alone.
8. **Columns.** NH_cloud = ∫_{r0}^{r1} n dr, and NH_full = ∫ n dr over the whole
   sampled ray.
9. **Resolution.** The sky pixel is set either in degrees or in parsecs at a
   reference distance. At distance d the dust map's effective resolution is
   max(voxel size, d × 14′). A run warns, and writes `OVERSAMP` to every header,
   when pixels are finer than that.

## 4. Dust-to-gas conversion

The cube holds τ₁₅₆₅ per pc. LightCube built it from Edenhofer's native E per pc,
where E is the extinction unit of Zhang, Green & Rix (2023), via:

| Step | Constant | Value | Source |
|---|---|---|---|
| E → A_V | A_V / E | **2.8** | Zhang et al. (2023), adopted by Edenhofer et al. (2024) |
| A_V → A₁₅₆₅ | A₁₅₆₅ / A_V | **2.7032** | Gordon et al. (2023) "G23" curve at R_V = 3.1, via `dust_extinction` 1.5 |
| A → τ | A / τ | **1.086** | 2.5 log₁₀ e |
| E → τ₁₅₆₅ | τ₁₅₆₅ / E | **6.9695** | product of the above |
| A_V → n_H | K_AV = n_H / (A_V/pc) | **812 cm⁻³** (default) | Bish et al. (2026) §III.3 |
| pc → cm | | 3.086×10¹⁸ | |

```math
n_{\rm H} \;=\; K_\tau\,\frac{\tau_{1565}}{\rm pc},\qquad
K_\tau \;=\; K_{A_V}\,\frac{1.086}{A_{1565}/A_V}
\;=\; 812\times\frac{1.086}{2.7032} \;=\; 326.2\ {\rm cm^{-3}}
```

| Option | n_H per (τ₁₅₆₅/pc) | per (E/pc) | per (A_V/pc) | implied N_H/A_V [cm⁻² mag⁻¹] | Source |
|---|---|---|---|---|---|
| **`lightcube`** (default) | **326.2** | 2273.6 | 812 | 2.51×10²¹ (**1.34 × Bohlin**) | Bish et al. 2026 §III.3 |
| `zucker21` | 237.0 | 1652 | 590 | 1.82×10²¹ | Zucker et al. 2021 / O'Neill et al. 2024 / McCallum et al. 2025 route |
| `hensley_draine23` | 387.4 | 2700 | 964 | 2.98×10²¹ | Hensley & Draine (2023) route |
| `sigma1565_legacy` | 257.2 | 1792 | 640 | 1.98×10²¹ | σ₁₅₆₅ = 1.26×10⁻²¹ cm²/H (Bohlin with an assumed A₁₅₆₅/A_V = 2.6) |

* **Bohlin caveat.** Bish et al. (2026) attribute 812 cm⁻³ per (A_V/pc) to Bohlin,
  Savage & Drake (1978). But it implies N_H/A_V = 2.51×10²¹ cm⁻² mag⁻¹, **~34%
  above the standard 1.87×10²¹** for R_V = 3.1. It is the default because it is
  internally consistent with the LightCube cube being sampled, not because it is
  the Bohlin value.
* **The conversion is a pure scalar on n, and therefore on αG ∝ 1/n.** It changes
  no boundary and no quality mask: geometry is decided in dust units, and
  `scripts/conversion_compare.py` shows the masks are bit-identical across all
  four options. Across these options it spans a factor of 1.63. **It is the
  largest arbitrary factor in the pipeline.**

Every constant and the full arithmetic: [docs/METHOD.md](docs/METHOD.md).

## 5. Worked example: Eos

Configuration `cloudfield/clouds/eos.yaml`:
* sky grid l = 20–50°, b = 30–70° at 0.1°;
* rays 69–260 pc at Δr = 0.25 pc;
* `half_max` boundaries with r_peak searched in 90–145 pc;
* statistics quoted inside the FIMS/SPEAR 10% H2/FUV contour, the published Eos
  boundary (Burkhart et al. 2025, Fig. 1): 32,219 pixels.

| Quantity (inside the contour, default conversion) | Median | 5–95% |
|---|---|---|
| G0_surface [Draine] | 0.493 | 0.435–0.554 |
| G0_peak [Draine] | 0.470 | 0.404–0.538 |
| n_peak [cm⁻³] | 2.67 | 1.07–9.53 |
| n_surface [cm⁻³] | 1.32 | 0.53–4.64 |
| **αG_peak** | **10.3** | 2.6–28.1 |
| αG_surface | 22.3 | 5.8–57.7 |
| NH_cloud [cm⁻²] | 1.46×10²⁰ | 0.65–2.24×10²⁰ |
| NH_full [cm⁻²] | 2.26×10²⁰ | 1.08–3.78×10²⁰ |
| r0 / r_peak / r1 [pc] | 102.5 / 113.8 / 129.5 | 94–124 / 106–130 / 118–136 |
| thickness [pc] | 22.0 | 7.8–35.5 |

* **Coverage:** every pixel in the contour passes the quality cut, and in none is
  r_peak pinned to the edge of the search window.
* **Pixels below αG = 1:** 339 in the contour (αG_peak), around MBM 40.
* **Attenuation:** G0_surface/G0_peak has a median of 1.05 (16–84%: 0.95–1.16).
  The voxel-noise floor is 0.12 and the radial gradient alone predicts 1.02, so
  **no attenuation is detected between r0 and the density peak.**
* **Oversampling:** the grid is oversampled ×4.9 (0.1° = 0.21 pc at 117.5 pc,
  against a 1 pc effective dust resolution).

The package reproduces the original Eos analysis bit-for-bit (all 17 output data
arrays).

## 6. Limitations

* **αG values are upper limits in clumped gas.** 3D dust reconstructions smooth
  out dense structure. At MBM 40 the Edenhofer map gives n_peak = 74.7 cm⁻³
  (nearest pixel) and 51.5 cm⁻³ (0.5° aperture), against the published
  n_H = 280 cm⁻³ (Karoly et al.). That is a density dilution of **3.75–5.4×**.
  Since αG ∝ 1/n, αG there is overestimated by about that factor, and the same
  applies to any unresolved condensation.
* **Per-pixel G0 is noise-dominated.**
  * LightCube's per-voxel signal-to-noise is **~11.5**.
  * G0 sampled on two surfaces **1.6 pc apart** along the line of sight
    correlates at only **ρ = 0.44**.
  * Sub-degree structure in G0_surface mostly reflects the boundary surface
    cutting obliquely through that voxel noise (`scripts/g0_diagnose.py`).
  * **The region median is the usable quantity.**
* **Chemical steady state is assumed.** αG and the Sternberg/Bialy theory describe
  H I/H2 equilibrium. That does not hold for every cloud; Eos is net dissociating
  (Burkhart et al. 2025).
* **αG_surface ≈ 2 × αG_peak under `half_max`, by construction.** r0 is the
  half-maximum crossing, so n_surface ≈ n_peak/2 (Eos median 0.988). The two maps
  are not independent measurements.
* **The dust-to-gas conversion** is a single arbitrary scalar (factor 1.63 across
  the options), and the default sits 34% above standard Bohlin (Section 4).
* **Edenhofer vs DUSTRIBUTION columns** (`scripts/column_systematics.py`; a
  systematic, not used in any map). Inside the Eos contour:

  | Comparison | Value |
  |---|---|
  | Edenhofer NH_full (69–260 pc) / DUSTRIBUTION cumulative column (0–350 pc) | 0.27 |
  | published Eos mass, DUSTRIBUTION / Edenhofer (Burkhart et al. 2025: 5.5×10³ vs 1.6×10³ M☉) | ~3.4 |
  | DUSTRIBUTION cumulative extinction / integral of its own median density cube | 2.43 |
  | the two density cubes over the same 69–260 pc | 0.93 |

  Over the same distances the two density cubes nearly agree. The ~3–4× gap comes
  from three things:
  * a median-of-integrals effect in DUSTRIBUTION, which builds up entirely inside
    the cloud;
  * foreground dust inside 69 pc, which only DUSTRIBUTION contains;
  * different conversions.

  Use DUSTRIBUTION's cumulative product for columns.
* **The GALFA H I comparison depends on geometry.** Edenhofer N_H / GALFA N_HI is
  0.60 against the full −90 to +90 km/s integral. It is 1.28 against the Eos
  velocity window (−8.6 to +1.4 km/s) with rays extended to 495 pc, and the window
  choice moves that within 1.19–1.34. This is **not** a measurement of the
  molecular fraction: the 34% Bohlin excess of the default conversion could
  account for all of it.
* **Cube extent.** The development cubes are clipped at r = 500 pc. Runs beyond the
  configured `cube_radius_pc` raise rather than read zeros. Distant clouds need the
  unclipped cubes ([data/README.md](data/README.md)).

## 7. Verification against Bish et al. (2026)

`scripts/verify_lightcube.py`:

| Check | Bish et al. 2026 | This pipeline | Result |
|---|---|---|---|
| Solar I_UV, mean over r < 25 pc | 0.633 (§III.1) | **0.6307** | pass (−0.4%) |
| Vertical scale height of I_UV | ~400 pc (§III.3) | **376 pc** (median profile, \|z\| = 0–480 pc) | pass |
| Ophiuchus / Taurus ISRF over n > 30 cm⁻³ | "roughly 3×" (§IV.3) | **4.0** in medians (5.3 in means) | close; depends on cloud volume |

* **Scale height.** The local scale height grows from ~370 pc to **459 pc** at
  larger |z|. That is consistent with the paper's own two reference points (0.633
  at the Sun, 12% of Standard at 800 pc), which imply an effective h ≈ 481 pc. The
  r = 500 pc cube clip limits the fit to |z| ≤ 480 pc.
* **Ophiuchus/Taurus.** The ratio is **insensitive to both the conversion option
  and the density cut**: 5.0–5.3 in means across all four conversions and cuts at
  20, 30 and 40 cm⁻³. What remains is the cloud volume:
  * Bish et al. do not state theirs.
  * This check uses the Zucker et al. (2021) Table 1 skeleton extents, only here
    and not as pipeline boundaries.

## 8. Installation and running

```bash
git clone https://github.com/b0burk02/cloudfield.git && cd cloudfield
pip install -e .
cp config.example.yaml config.yaml      # fill in the data paths (data/README.md)
cloudfield check-data                   # verifies every input, incl. the 0.633 solar test
cloudfield run eos --out output/        # FITS + PNG in output/
```

The validation scripts in `scripts/` read the same `config.yaml` and the
`output/` of `cloudfield run eos`.

## 9. Adding a new cloud

1. **Write `cloudfield/clouds/<name>.yaml`.** Copy `eos.yaml` and set:
   * `name`, `distance_pc`
   * `grid` (l, b and the pixel scale; `{mode: physical, pc, ref_distance_pc}`
     keeps the physical resolution fixed)
   * `radial` (r_min, r_max, Δr), staying inside the cubes
   * `boundary`
   * `region` (`lb_box`, `healpix_contour` or `all`)
   * optional `annotations` and `validation_points`
2. **For a 3D cloud mask,** [docs/BOUNDARIES.md §4](docs/BOUNDARIES.md) lists
   exactly what the mask file must specify. The `mask3d` interface and its rules
   are fixed; the file loader will be implemented against the first real mask.
3. **Before trusting the output, check:**
   * the `OVERSAMP` factor;
   * the `BOUNDARY_FLAGS` fractions printed by the run;
   * the r_peak histogram (no pile-up on the window edge);
   * adjacent-pixel r0 jumps (`scripts/g0_diagnose.py`).

## Citing

**Cite the data products, not just this code:**

* **LightCube:** Bish, H. V., Peek, J. E. G., Zucker, C., et al. 2026,
  [arXiv:2608.27389](https://arxiv.org/abs/2608.27389).
* **Edenhofer, G., et al. 2024,** A&A 685, A82 (3D dust;
  [arXiv:2308.01295](https://arxiv.org/abs/2308.01295)).
* **Any dust map you use for comparison,** e.g. DUSTRIBUTION (Dharmawardena et
  al.; [mwdust.com](http://www.mwdust.com)).
* **Theory:** Sternberg, A., et al. 2014, ApJ 790, 10; Bialy, S. & Sternberg, A.
  2016, ApJ 822, 83 ([arXiv:1601.02608](https://arxiv.org/abs/1601.02608)).
* **Extinction curve:** Gordon, K. D., et al. 2023, ApJ 950, 86.
* **Eos:** Burkhart, B., et al. 2025.
* **Other references:**
  * Draine 1978, ApJS 36, 595
  * Bohlin, Savage & Drake 1978, ApJ 224, 132
  * Zhang, Green & Rix 2023
  * Zucker et al. 2021, ApJ 919, 35
  * Hensley & Draine 2023
  * Peek et al. 2018 (GALFA-HI DR2)

Software citation: [CITATION.cff](CITATION.cff). License: [BSD-3-Clause](LICENSE).
