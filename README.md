# cloudfield

Maps of the far-UV field incident on a molecular cloud (**G0**) and the H I-to-H2
transition parameter **αG**. The field comes from LightCube (Bish et al. 2026);
the gas density comes from the Edenhofer et al. (2024) 3D dust map.

```math
\alpha G = 0.59\, I_{\rm UV}\left(\frac{100\ {\rm cm^{-3}}}{n}\right)\left(\frac{9.9}{1+8.9\,\tilde\sigma_g}\right)^{0.37} = \frac{59\, I_{\rm UV}}{n}
```

This is Bialy & Sternberg (2016), Eq. 22, with σ_g = 1.9×10⁻²¹ cm². I_UV is in
units of the Draine (1978) field, and n is the total hydrogen-nucleus density
(n_HI + 2n_H2).

## Method

1. **Inputs.**
   * **I_UV** is LightCube's combined 912–2000 Å cube, already in Draine units.
   * **n** comes from the Edenhofer dust map.
   * Both must come from Edenhofer dust, because LightCube's radiative transfer
     was computed through that map.
2. **Rays.** A ray is cast from the Sun through each (l, b) pixel. The dust
   density and I_UV are sampled along it at fine radial steps.
3. **Cloud boundaries.** Along each ray, the code finds where the cloud starts
   (r0), where it peaks (r_peak) and where it ends (r1). The default method uses
   the half-maximum points of the dust profile around its peak. A 3D cloud mask
   can be used instead.
4. **Read the field at the surface.** G0 = I_UV(r0) is the field arriving at the
   cloud. It is not averaged along the line of sight, which would mix in gas in
   front of and behind the cloud.
5. **Compute αG.** αG = 59 I_UV / n is evaluated at the density peak (the primary
   map) and at the surface.
6. **Dust to gas.** τ₁₅₆₅ per pc is converted to n with
   n = 812 cm⁻³ per (A_V/pc) × 1.086 / (A₁₅₆₅/A_V), where A₁₅₆₅/A_V = 2.70
   (Gordon et al. 2023). This gives n = 326 cm⁻³ per (τ₁₅₆₅/pc). The conversion is
   a single scalar on n, and therefore on αG.

αG values are upper limits in dense, clumpy gas, because 3D dust maps smooth out
small condensations.

## Usage

```bash
pip install -e .
cp config.example.yaml config.yaml   # set paths to the data (see data/README.md)
cloudfield check-data
cloudfield run eos --out output/
```

Each cloud is a YAML file in `cloudfield/clouds/`; Eos is the example.

## More

* [docs/METHOD.md](docs/METHOD.md): the conversion chain, constant by constant.
* [docs/BOUNDARIES.md](docs/BOUNDARIES.md): the boundary methods and the 3D mask format.
* [docs/DETAILS.md](docs/DETAILS.md): full reference, covering the Eos results,
  limitations and verification.
* [data/README.md](data/README.md): where to download the inputs. No data is
  included in this repo.

## Citing

Please cite the data this builds on:
* LightCube: Bish et al. 2026, [arXiv:2608.27389](https://arxiv.org/abs/2608.27389)
* Edenhofer et al. 2024, A&A 685, A82
* Bialy & Sternberg 2016, ApJ 822, 83

See also [CITATION.cff](CITATION.cff). License: BSD-3-Clause.
