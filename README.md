# cloudfield

Maps of the incident far-UV field **G0** and the H I-to-H2 transition parameter
**alphaG** for a molecular cloud with a defined 3D volume. Inputs are the
LightCube 3D model of the local FUV interstellar radiation field (Bish et al.
2026) and the Edenhofer et al. (2024) 3D dust map. Any cloud is described by a
YAML config; Eos (Burkhart et al. 2025) is the worked example.

> **First thing to run:** `cloudfield check-data` (see [Data](#2-inputs) and
> [data/README.md](data/README.md)). No data products are in this repository.

---

## 1. What it produces

For every sightline (pixel) on a Galactic (l, b) grid, the code finds where the
ray enters the cloud (r0), where its dust density peaks (r_peak) and where it
leaves (r1). It then reports:

| Map | Meaning |
|---|---|
| `G0_SURFACE` | I_UV at r0: the field incident on the cloud surface |
| `G0_PEAK` | I_UV at r_peak |
| `N_SURFACE`, `N_PEAK` | n_H at r0 and at r_peak [cm^-3] |
| `ALPHAG_PEAK` (primary), `ALPHAG_SURFACE` | alphaG = 59 I_UV / n at r_peak and at r0 |
| `NH_CLOUD` | ∫ n dr from r0 to r1 [cm^-2] |
| `NH_FULL` | ∫ n dr over the whole sampled ray [cm^-2] |
| `R0`, `RPEAK`, `R1`, `THICKNESS` | boundary radii [pc] |
| `ALPHAG_LOSMEAN`, `RATIO_SURF_LOSMEAN` | old slab estimator 59 ⟨I⟩_n/⟨n⟩_n (comparison only) |
| `QUAL_NPEAK`, `QUAL_RPEAK_INTERIOR` | quality masks |
| `BOUNDARY_FLAGS`, `BOUNDARY_METHOD` | per-pixel diagnostics of the boundary finder |
| region extension (e.g. `CLOUD_H2_10PCT`) | the pixels statistics are quoted over |

**alphaG.** Following Sternberg et al. (2014) and Bialy & Sternberg (2016),

    alphaG = 0.59 I_UV (100 cm^-3 / n) [9.9 / (1 + 8.9 σ̃_g)]^0.37  =  59 I_UV / n

for σ̃_g = 1 (σ_g = 1.9e-21 cm² per H nucleus, solar metallicity).

* **I_UV** is the free-space FUV field in units of the **Draine (1978) Standard
  ISRF**.
* **n** is the **total hydrogen-nucleus density, n_HI + 2 n_H2**, not the atomic
  density.
* **Reading the value.** alphaG ≫ 1 means the H I shielding layer dominates (a
  sharp H I-to-H2 transition, large H I column). alphaG ≲ 1 means dust shielding
  dominates.

Outputs: `<prefix>_G0.fits`, `<prefix>_alphaG.fits`, `<prefix>_maps.fits` (all
with a Galactic CAR WCS), and `<prefix>_maps.png`.

## 2. Inputs

* **I_UV:** LightCube, the **combined 912–2000 Å cube**. It is already normalized
  to the Draine field (header `BUNIT = 'Draine78 ISRF'`) and needs **no
  rescaling**.
  * Diagnostic for any new file: the mean over voxels within 25 pc of the cube
    origin must be **0.633** (Bish et al. 2026 §III.1). We measure 0.6307.
  * A per-band cube instead of the combined one gets alphaG wrong by orders of
    magnitude with no error raised. `cloudfield check-data` and every run check
    for this.
* **n:** the Edenhofer et al. (2024) map, as tau_1565 per pc on the LightCube grid.
* **Both I_UV and n come from Edenhofer dust.** LightCube's radiative transfer was
  computed through the Edenhofer map, so pairing its I_UV with any other dust map
  would be inconsistent. The code raises if asked to.
  * DUSTRIBUTION (Dharmawardena et al.) is supported **only for column-density
    comparison** (`cloudfield/columns.py`).

Where to download each file, and which one: [data/README.md](data/README.md).

## 3. Method

1. **Rays.** Each (l, b) pixel centre defines a ray from the Sun. It is sampled
   every `dr` (0.25 pc for Eos) between `r_min` and `r_max`. Both cubes are
   trilinearly interpolated at x = r cos b cos l, y = r cos b sin l, z = r sin b.
   r < 69 pc is masked (no Edenhofer dust there).
2. **Dust units first.** All geometry (the peak, r0, r1) and the quality cut are
   decided on the dust profile tau_1565(r), never on n. The gas conversion
   (Section 4) only scales values, so changing it cannot move a boundary or change
   which pixels pass. The quality cut is tau_1565,peak ≥ 9.196e-4 pc^-1 (n_peak =
   0.3 cm^-3 under the default conversion).
3. **Boundaries.** A pluggable method returns r0, r_peak and r1 per pixel, plus
   flags; see [docs/BOUNDARIES.md](docs/BOUNDARIES.md).
   * `half_max`: peak in a distance window; r0 and r1 are the near and far
     half-maximum crossings. The Eos default.
   * `fixed_range`, `nearest_peak` (multi-component profiles).
   * `mask3d`: an externally supplied 3D volume. Its interface and segment rules
     are fixed; its file loader is pending a real mask.

   Every pixel records which method produced its boundary and whether the profile
   was multi-peaked, or the boundary hit a search edge.
4. **Read the field at the surface.** G0_surface = I_UV(r0). alphaG describes the
   field that *arrives at* the cloud. A density-weighted average of I_UV along
   the sightline (the old slab estimator, kept as `ALPHAG_LOSMEAN`) mixes in field
   values from in front of and behind the cloud. For Eos it differs from the
   surface value by a factor of 1.24 (median; 16–84%: 1.10–1.38). That difference
   correlates with cloud thickness (Spearman ρ = 0.70), so it is a structured bias,
   not noise.
5. **Attenuation test.** G0_surface/G0_peak is compared with the per-voxel noise
   of the LightCube cube, which is estimated robustly from relative first
   differences. It is also compared with the ratio expected from the large-scale
   radial gradient alone. For Eos: median 1.05 (16–84%: 0.95–1.16) against a
   two-voxel noise of 0.12 and a trend prediction of 1.02. **No attenuation is
   detected between r0 and the density peak.**
6. **Grid and resolution.** The sky pixel is configurable, either fixed in degrees
   or fixed in parsecs at a reference distance. The effective resolution of the
   dust map at distance d is max(voxel size, d × native angular resolution). The
   run warns, and writes `OVERSAMP` to the headers, when pixels are finer than
   that. Eos uses 0.1° = 0.21 pc at 117.5 pc against a 1 pc effective resolution,
   so it is **oversampled ×4.9**: neighbouring pixels are not independent.

## 4. Dust-to-gas conversion

The cube holds tau_1565/pc. A gas conversion quoted per E/pc or per A_V/pc maps to
it through tau_1565 = E × 2.8 × A_1565/A_V / 1.086. Here A_1565/A_V = **2.703**,
from Gordon et al. (2023), R_V = 3.1, via `dust_extinction` 1.5. Full arithmetic:
[docs/METHOD.md](docs/METHOD.md).

| Option | n_H per (tau_1565/pc) | per (E/pc) | implied N_H/A_V [cm^-2/mag] | Source |
|---|---|---|---|---|
| **`lightcube`** (default) | **326.2** | 2273.6 | 2.51e21 (1.34 × Bohlin) | Bish et al. 2026 §III.3: 812 cm^-3 per (A_V/pc) |
| `zucker21` | 237.0 | 1652 | 1.82e21 | Zucker 2021 / O'Neill 2024 / McCallum 2025 route |
| `hensley_draine23` | 387.4 | 2700 | 2.98e21 | Hensley & Draine (2023) route |
| `sigma1565_legacy` | 257.2 | 1792 | 1.98e21 | σ_1565 = 1.26e-21 cm²/H (Bohlin with an assumed A_1565/A_V = 2.6) |

* **The default, `lightcube`, is not the standard Bohlin value.** Bish et al.
  attribute 812 cm^-3 per (A_V/pc) to Bohlin et al. (1978), but it implies
  N_H/A_V = 2.51e21, **~34% above the standard 1.87e21** for R_V = 3.1. It is the
  default because it is internally consistent with the LightCube cube we sample,
  not because it matches Bohlin.
* **The conversion is a single multiplicative scalar on n,** and therefore on
  alphaG ∝ 1/n. Across these options it spans a factor of 1.63. **It is the
  largest arbitrary factor in the pipeline.**
* `scripts/conversion_compare.py` shows that masks and geometry are bit-identical
  across options, and tabulates what changes.

## 5. Limitations

* **alphaG values are upper limits in clumped gas.** 3D dust reconstructions smooth
  out dense structure. At MBM 40 the Edenhofer map recovers only 1/3.75 (nearest
  pixel) to 1/5.4 (0.5° aperture) of the published n_H = 280 cm^-3 (Karoly et
  al.). alphaG ∝ 1/n there is overestimated by about that factor, and the same
  applies to any unresolved condensation.
* **Per-pixel G0 is not a measurement.**
  * LightCube has a per-voxel signal-to-noise of ~11.5 (robust estimate).
  * G0 sampled 1.6 pc apart along the sightline correlates at only ρ = 0.44.
  * Structure in G0_surface at sub-degree scales mostly reflects the boundary
    surface cutting through that voxel noise (`scripts/g0_diagnose.py`).
  * **Use the region median**; Eos: 0.493 (5–95%: 0.435–0.554).
* **alphaG_surface ≈ 2 × alphaG_peak under `half_max`, by construction.** r0 is
  the half-maximum crossing, so n_surface ≈ n_peak/2 (Eos median 0.988). The two
  are not independent measurements.
* **Edenhofer vs DUSTRIBUTION (column systematic, not used in any map).** Inside
  the Eos 10% contour (`scripts/column_systematics.py`):

  | Comparison | Value |
  |---|---|
  | Edenhofer NH_full (69–260 pc) / DUSTRIBUTION cumulative column (0–350 pc) | 0.27 (per-pixel median 0.29) |
  | published Eos mass, DUSTRIBUTION / Edenhofer (Burkhart et al. 2025: 5.5e3 vs 1.6e3 M☉) | ~3.4 |
  | DUSTRIBUTION cumulative extinction / integral of its own median density cube | 2.43 (16–84%: 2.06–2.74) |
  | the two density cubes over the same 69–260 pc | 0.93 |

  Over the same distances the density cubes nearly agree. The ~3–4× gap comes from
  three things:
  * a median-of-integrals effect in DUSTRIBUTION, building up entirely inside the
    cloud (90–140 pc);
  * foreground dust inside 69 pc, which only DUSTRIBUTION contains;
  * different conversions.

  The disagreement is unresolved. **Use DUSTRIBUTION's cumulative product, not an
  integral of its median density cube, for columns.**
* **alphaG assumes H I/H2 chemical steady state.** That does not hold for every
  cloud; Eos is net dissociating (Burkhart et al. 2025).
* **Cube extent.** The development cubes are clipped at r = 500 pc. Runs beyond
  the configured `cube_radius_pc` raise. The unclipped cubes are needed for
  distant clouds (data/README.md).
* **The comparison with GALFA H I is sensitive to geometry.** Edenhofer N_H
  (69–495 pc) / GALFA N_HI is 1.28 against the Eos velocity window
  (−8.6 to +1.4 km/s), and 0.60 against the full −90..+90 km/s integral.
  * The window choice moves it within 1.19–1.34.
  * This is **not** a measurement of the molecular fraction. The 34% excess of
    the `lightcube` conversion over standard Bohlin could account for all of it.

## 6. Verification against Bish et al. (2026)

`scripts/verify_lightcube.py`:

| Check | Bish et al. 2026 | This pipeline | Result |
|---|---|---|---|
| Solar I_UV (mean, r < 25 pc) | 0.633 (§III.1) | **0.6307** | pass (−0.4%) |
| Vertical scale height | ~400 pc (§III.3) | **376 pc** (median profile, \|z\| = 0–480 pc) | pass |
| Ophiuchus / Taurus ISRF over n > 30 cm^-3 | "roughly 3×" (§IV.3) | **4.0** (medians), 5.3 (means) | outside 2–4 in means |

* **Scale height.** The local scale height grows from ~370 pc to 459 pc at larger
  |z|. That is consistent with the paper's own two reference points (0.633 at the
  Sun, 12% of Standard at 800 pc), which imply an effective h ≈ 481 pc. The r500
  clip limits the fit to |z| ≤ 480 pc.
* **Ophiuchus/Taurus.** The ratio is **insensitive to both the conversion option
  and the density cut** (5.0–5.3 in means across all four conversions and cuts at
  20, 30 and 40 cm^-3). The remaining difference lies in the cloud volumes:
  * Bish et al. do not state theirs.
  * The check uses the Zucker et al. (2021) Table 1 skeleton extents. These are
    used there only, not as pipeline boundaries.

## 7. Worked example: Eos

```bash
pip install -e .
cp config.example.yaml config.yaml && $EDITOR config.yaml   # set data paths
cloudfield check-data
cloudfield run eos --out output/
```

`cloudfield/clouds/eos.yaml` (abridged):

```yaml
name: Eos
grid:   {l: [20, 50], b: [30, 70], pixel: {mode: fixed, deg: 0.1}}
radial: {r_min: 69.0, r_max: 260.0, dr: 0.25}
boundary: {method: half_max, window: [90.0, 145.0]}
region:
  type: healpix_contour          # FIMS/SPEAR 10% H2/FUV contour (Burkhart+2025 Fig. 1)
  data_key: fims_h2
  level: 10.0
  pick_box: {l: [25, 45], b: [40, 63]}
  extname: CLOUD_H2_10PCT
validation_points:
  - {name: "MBM 40", l: 37.52, b: 44.57, radius_deg: 0.5, published_n_H: 280.0}
```

Eos results inside the 10% contour (32,219 pixels), default conversion:

| | median | 5–95% |
|---|---|---|
| G0_surface | 0.493 | 0.435–0.554 |
| n_peak [cm^-3] | 2.67 | 1.07–9.53 |
| alphaG_peak | 10.3 | 2.6–28.1 |
| alphaG_surface | 22.3 | 5.8–57.7 |
| NH_full [cm^-2] | 2.26e20 | 1.08–3.78e20 |
| r0 / r_peak / r1 [pc] | 102.5 / 113.8 / 129.5 | |

This refactor reproduces the original single-purpose Eos script's maps
**bit-for-bit** (all 17 data arrays).

## 8. Adding a new cloud

1. **Write `cloudfield/clouds/<name>.yaml`** (or any path). The fields:
   * `name`, `distance_pc` (a range, used for the resolution report)
   * `grid` (`l`, `b`, `pixel`: `{mode: fixed, deg}` or
     `{mode: physical, pc, ref_distance_pc}`)
   * `radial` (`r_min`, `r_max`, `dr`); `r_max` must stay inside the cube
   * `boundary` (method plus parameters; optional `fallback`)
   * `region` (`lb_box`, `healpix_contour` or `all`)
   * `conversion`; optional `annotations`, `validation_points`,
     `noise_block_xyz_pc`
2. **For a 3D mask,** read [docs/BOUNDARIES.md](docs/BOUNDARIES.md), Section 4, for
   exactly what the file must specify. The `mask3d` loader is implemented against
   the first real file supplied.
3. **Run it:** `cloudfield run <name>`.
4. **Check before trusting the output:**
   * the printed `OVERSAMP`;
   * the `BOUNDARY_FLAGS` fractions (MULTIPEAK, RPEAK_EDGE, MASK_*);
   * the r_peak histogram (no pile-up on the window edge);
   * adjacent-pixel r0 jumps (`scripts/g0_diagnose.py`);
   * that `QUAL_NPEAK` covers the region.

## Citing

If you use this code, **cite the data products it is built on**, not just the code:

* **LightCube:** Bish, H. V., et al. 2026, arXiv:2608.27389.
* **Edenhofer et al. 2024,** A&A 685, A82 (3D dust).
* **The dust map you select for any comparison,** e.g. DUSTRIBUTION
  (Dharmawardena et al.; www.mwdust.com).
* **The H I-to-H2 theory:** Sternberg et al. 2014, ApJ 790, 10; Bialy & Sternberg
  2016, ApJ 822, 83.
* **The extinction curve:** Gordon et al. 2023, ApJ 950, 86.
* **Eos:** Burkhart et al. 2025.

Other references: Draine 1978, ApJS 36, 595; Bohlin, Savage & Drake 1978,
ApJ 224, 132; Zhang, Green & Rix 2023; Zucker et al. 2021, ApJ 919, 35; Hensley &
Draine 2023; Peek et al. 2018 (GALFA-HI DR2).

Software citation: [CITATION.cff](CITATION.cff).
