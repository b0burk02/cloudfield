# Cloud boundaries: the plugin contract and the `mask3d` file specification

`cloudfield` reads the UV field at the cloud **surface**, so where a sightline enters
the cloud (r0), where its density peaks (r_peak) and where it exits (r1) set every
G0 and alphaG value. Boundary finders are plugins (`cloudfield/boundaries.py`) with
one shared contract. This document defines that contract. It also defines exactly
what a 3D mask file must provide so it can be used with `mask3d`.

## 1. The plugin contract (every method)

**Input.** The method is called once per chunk of P sightlines, with:

| Argument | Shape | Meaning |
|---|---|---|
| `D` | (P, R) float | Dust quantity along each ray: tau_1565 per pc for Edenhofer. NaN where undefined (r < 69 pc). |
| `ctx.l`, `ctx.b` | (P,) | Galactic longitude/latitude of each sightline [deg]. These are the pixel centres of the (l, b) grid. |
| `ctx.r` | (R,) | The radial sampling array [pc] (uniform; `radial: {r_min, r_max, dr}`). |
| `ctx.X, ctx.Y, ctx.Z` | (P, R) | Heliocentric Galactic Cartesian position of every sample [pc]: x toward l=0, y toward l=90°, z toward b=+90°. |

Boundaries are decided on **D, never on the gas density n**. The geometry, and
therefore G0, is independent of the dust-to-gas conversion. `mask3d`
additionally receives a mask object (Section 2).

**Output.** A `BoundaryResult` with one entry per ray:

| Field | Type | Meaning |
|---|---|---|
| `i0`, `i1` | int | Indices into `r` of the near (r0) and far (r1) boundary samples |
| `ipk` | int | Index of the density peak used (r_peak) |
| `has0`, `has1` | bool | Whether r0 / r1 exist. False makes that boundary NaN in the maps |
| `valid` | bool | Method-level validity (for example, the ray actually crosses the mask). AND-ed with the dust-unit quality cut |
| `rpeak_edge` | bool | The peak sits on the edge of the search range |
| `flags` | uint16 | Bitmask below. Flags **report**; they never silently change a value |

Per pixel, the pipeline writes `BOUNDARY_FLAGS` (the bitmask) and `BOUNDARY_METHOD`
(which method produced the boundary: 0 none, 1 half_max, 2 fixed_range,
3 nearest_peak, 4 mask3d) to `<cloud>_maps.fits`. A cloud config can name a
`fallback` method: it is used only where the primary returns `valid=False`, and
`BOUNDARY_METHOD` records which method was used in each pixel.

| Bit | Flag | Set when |
|---|---|---|
| 1 | MULTIPEAK | the profile has >1 separate run with D ≥ D_peak/2 |
| 2 | RPEAK_EDGE | the peak sits on the edge of the search window / range / mask segment |
| 4 | R0_MISSING | no near-side boundary found |
| 8 | R1_MISSING | no far-side boundary found |
| 16 | R0_RAY_START | r0 is the first radial sample (the true boundary may be closer) |
| 32 | R1_RAY_END | r1 is the last radial sample (the true boundary may be farther) |
| 64 | VALLEY_BOUND | nearest_peak placed a boundary at a valley between components |
| 128 | MASK_GRAZING | mask3d: chosen in-mask chord shorter than `min_chord_pc` |
| 256 | MASK_MULTISEG | mask3d: ray enters and exits the mask more than once |
| 512 | MASK_TRUNCATED | mask3d: segment touches the sampling range or the mask coverage edge |
| 1024 | MASK_MISS | mask3d: ray never enters the mask |

## 2. Methods

| Method | Parameters | r0 / r_peak / r1 | Use |
|---|---|---|---|
| `half_max` | `window` [pc] | peak of D in the window; near/far half-max crossings (may lie outside the window) | Single clean component. Eos default |
| `fixed_range` | `r0_pc`, `r1_pc` | the supplied radii; peak = argmax between them | Known distance range |
| `nearest_peak` | `window`, `target_pc`, `min_frac`=0.5, `smooth_pc`=2, `min_dip`=0.2 | peak of the smoothed profile nearest `target_pc`; half-max crossings, or the valley minimum if a *distinct* neighbouring component (dip below (1−min_dip)×the lower peak) is not separated by a half-max crossing | Several components along the sightline |
| `mask3d` | `mask`, `min_chord_pc`, `segment`, `target_pc` | entry/exit of the ray through the mask; peak = argmax inside the chosen segment | Primary method once masks exist |

How fragile is half_max? On Eos/Edenhofer, r0 never jumps more than 6.5 pc
between adjacent pixels. On DUSTRIBUTION it jumped by up to 28 pc, where the
half-max crossing snapped to a different density feature. `mask3d` removes that
failure mode by taking the volume from an external definition.

## 3. `mask3d` behaviour (implemented; only the file loader is pending)

`segments_to_boundaries()` implements the rules below and has been exercised on
real data. `load_mask3d()` raises `NotImplementedError` until a real file is
supplied. The format will be matched to that file, not guessed.

* **Segments.** A segment is a maximal run of consecutive radial samples that are
  inside the mask and inside its coverage.
* **Multiple entries and exits.** When a ray crosses the mask more than once, one
  segment is chosen and MASK_MULTISEG is set. Segments are never merged. The
  choice is set by `segment`:
  * `peak` (default): the segment containing the maximum of D among in-mask samples.
  * `nearest`: midpoint nearest `target_pc`.
  * `longest`: the segment with the most samples.
  * `first`: the segment nearest the Sun.
* **Boundaries.** r0 = the first in-mask sample of the segment (entry) and
  r1 = the last (exit). Unlike half_max, r0 is inside the cloud.
* **Grazing rays (short path).** If the segment's chord, n_samples × dr, is shorter
  than `min_chord_pc` (default 2 × the mask voxel size), the ray only clips the
  mask edge tangentially. The pixel gets `valid=False` (NaN, or the fallback
  method) and MASK_GRAZING. Rationale: at a tangent point, nearest-neighbour
  resampling of a voxelized surface decides entry and exit by single voxels, so
  r0 and r1 are not measurements.
* **Masks that do not cover the full radial range.** The mask object reports
  coverage separately from membership (`defined()` vs `contains()`). Samples outside
  coverage are treated as unknown, not "outside the cloud". If the chosen segment
  touches an undefined sample, or the first or last radial sample, MASK_TRUNCATED
  is set. The corresponding boundary is then a limit, not a measurement.
* **Rays that miss.** The pixel gets `valid=False` and MASK_MISS.
* **Mask grid differs from the dust/ISRF grid.** The mask is never regridded onto
  the dust grid. It is evaluated directly at every ray sample position
  (`ctx.X, ctx.Y, ctx.Z`, at the radial sampling `dr`, 0.25 pc for Eos) in the
  mask's own coordinates:
  * **Boolean masks:** nearest-neighbour lookup on the mask grid. Never interpolated.
  * **Float masks** (probability, density, label): trilinear interpolation on the
    mask grid, then the stated threshold or label selection.
  * **Spherical grids** (e.g. the native Edenhofer HEALPix × distance grid): the
    sample's (l, b, d) is converted to (HEALPix pixel, distance bin) with the
    file's own Nside, ordering and bin edges.

  Choose `dr` well below the mask voxel size, so entry and exit are located to
  better than one voxel.

## 4. What a mask file must provide (please send this section to the mask supplier)

Any format is fine (FITS preferred). The file, or a short note with it, must state:

1. **Grid type.** Either:
   * Cartesian, or
   * spherical (HEALPix pixel × distance bin, as in the Edenhofer et al. 2024 map).
2. **Frame and origin.**
   * Cartesian: heliocentric or Galactocentric; the axis directions; and the
     position of the Sun in pixel coordinates (0- or 1-based, voxel centre or
     corner).
   * Spherical: Galactic or equatorial; HEALPix Nside, ordering (RING/NESTED) and
     coordinate system.
3. **Voxel size and axis order.**
   * Cartesian: the voxel size in pc per axis, and the array axis order
     (e.g. [z, y, x]).
   * Spherical: the distance bin edges in pc.
4. **Semantics of the values.** Boolean (1 = cloud)? A probability or density
   with a threshold to apply? Integer labels, and which label is the cloud?
5. **Coverage.** The volume over which the mask is defined: the radial range, and
   any blanked region. Use a separate "undefined" value (e.g. NaN or −1),
   distinct from 0 = "not cloud".
6. **Provenance.**
   * which dust map (and which sample or statistic: mean, median, sample k);
   * what threshold or algorithm produced the mask;
   * the reference to cite.

If the mask is on the same grid as the LightCube/Edenhofer cube used here (1 pc,
heliocentric Galactic Cartesian, [z, y, x], Sun at index 499.5 in the
r = 500 pc product), say so; that is the simplest case.

## 5. Checking a boundary result before trusting it

* `BOUNDARY_FLAGS`: the fraction of pixels in the region with MULTIPEAK,
  RPEAK_EDGE, MASK_MULTISEG, MASK_TRUNCATED and MASK_GRAZING. `cloudfield run`
  prints these.
* `BOUNDARY_METHOD`: how many pixels came from the fallback.
* **Adjacent-pixel r0 jumps.** Jumps of more than a few pc between neighbouring
  pixels mean the boundary snapped between features (`scripts/g0_diagnose.py`).
* **r_peak histogram.** A pile-up on the window edge means the window is too
  narrow.
