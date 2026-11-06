"""Dust-to-gas conversion for the Edenhofer tau_1565 cube.

The cube is tau_1565 per pc. LightCube built it from Edenhofer+24 E (Zhang,
Green & Rix 2023 units) per pc via
  E -> A_V        : x E_TO_AV (Zhang+23, adopted by Edenhofer+24)
  A_V -> A_1565   : Gordon et al. (2023) extinction curve, R_V = RV
  A_1565 -> tau   : / 1.086
so tau_1565 = E * E_TO_AV * A1565_AV / 1.086, and a gas conversion quoted per
(E/pc) maps to per (tau_1565/pc) as  n/tau = (n/E) * 1.086 / (E_TO_AV * A1565_AV).

The choice is a single multiplicative scalar on n (and therefore on alphaG);
it is the largest arbitrary factor in the pipeline. See docs/METHOD.md.
"""
from .constants import PC_CM

E_TO_AV = 2.8
RV = 3.1
LAMBDA_UV_UM = 0.1565
MAG_PER_TAU = 1.086


def _a1565_over_av():
    try:
        import astropy.units as u
        from dust_extinction.parameter_averages import G23
    except ImportError as e:
        raise ImportError("dust_extinction>=1.5 (G23 model) is required for A(1565)/A(V): "
                          "pip install dust_extinction==1.5") from e
    return float(G23(Rv=RV)(LAMBDA_UV_UM * u.micron))


def dust_extinction_version():
    try:
        import dust_extinction
        return dust_extinction.__version__
    except Exception:
        return "unknown"


A1565_AV = _a1565_over_av()          # 2.7032 for G23, R_V=3.1 (dust_extinction 1.5)
N_H_PER_AV = 812.0                   # cm^-3 per (A_V/pc); Bish et al. 2026 Sec 3.3
NH_PER_AV_BOHLIN_STD = 1.87e21       # cm^-2 mag^-1; standard Bohlin+78 value for R_V = 3.1

# Named options. Each gives n_H in cm^-3 per (E/pc) or per (tau_1565/pc).
NH_CONVERSIONS = {
    # DEFAULT. "lightcube": 812 cm^-3 per (A_V/pc), Bish et al. 2026 Sec 3.3
    # -> 812 * 2.8 = 2273.6 per E/pc. They attribute this to Bohlin et al. (1978),
    # but it implies N_H/A_V = 812 pc = 2.51e21 cm^-2/mag, ~34% above the standard
    # 1.87e21 for R_V = 3.1. Adopted for internal consistency with the LightCube
    # cube we sample, NOT because it is the Bohlin value.
    "lightcube": dict(
        per_E=N_H_PER_AV * E_TO_AV,
        source="Bish+2026 Sec 3.3: 812 cm^-3 per (A_V/pc) x 2.8",
        caveat=["Bish+2026 attribute 812 per A_V/pc to Bohlin+1978, but it",
                f"implies N_H/A_V = {N_H_PER_AV * PC_CM:.3g} cm^-2/mag, "
                f"~{100 * (N_H_PER_AV * PC_CM / NH_PER_AV_BOHLIN_STD - 1):.0f}% above the "
                "standard",
                "1.87e21 for R_V=3.1. Adopted for consistency with the LightCube",
                "cube we sample, NOT because it is the Bohlin value."]),
    # Zucker+21 / O'Neill+24 / McCallum+25 route
    "zucker21": dict(per_E=1652.0,
                     source="Zucker+2021 / O'Neill+2024 / McCallum+2025: 1652 cm^-3 per (E/pc)"),
    # Hensley & Draine 2023 route (approximate)
    "hensley_draine23": dict(per_E=2700.0,
                             source="Hensley & Draine 2023 route: ~2700 cm^-3 per (E/pc)"),
    # Previous default: sigma_1565 = 1.26e-21 cm^2/H, i.e. Bohlin N_H/A_V = 1.87e21
    # divided by an ASSUMED A(1565)/A(V) ~ 2.6 (~1722 per E/pc with that ratio).
    "sigma1565_legacy": dict(per_tau=1.0 / (1.26e-21 * PC_CM),
                             source="legacy: sigma_1565 = 1.26e-21 cm^2/H (Bohlin / assumed "
                                    "A1565/AV=2.6)"),
}
DEFAULT_CONVERSION = "lightcube"


def nh_per_tau(name=None):
    """n_H [cm^-3] per (tau_1565/pc) for a named conversion."""
    c = NH_CONVERSIONS[name or DEFAULT_CONVERSION]
    if "per_tau" in c:
        return c["per_tau"]
    return c["per_E"] * MAG_PER_TAU / (E_TO_AV * A1565_AV)


def nh_per_E(name=None):
    """n_H [cm^-3] per (E/pc) for a named conversion (via the same chain)."""
    c = NH_CONVERSIONS[name or DEFAULT_CONVERSION]
    if "per_E" in c:
        return c["per_E"]
    return c["per_tau"] * E_TO_AV * A1565_AV / MAG_PER_TAU


def nh_per_av(name=None):
    """n_H [cm^-3] per (A_V/pc)."""
    return nh_per_E(name) / E_TO_AV


def conversion_summary(name=None):
    name = name or DEFAULT_CONVERSION
    return (f"n_H = {nh_per_tau(name):.4g} cm^-3 per (tau_1565/pc) "
            f"[{nh_per_E(name):.4g} per (E/pc)]; {name}: {NH_CONVERSIONS[name]['source']}; "
            f"A1565/AV={A1565_AV:.4f} (G23, Rv={RV}), E->A_V x{E_TO_AV}")


def tau_peak_min(n_peak_min=0.3, reference="lightcube"):
    """Quality cut in DUST units (tau_1565/pc): the tau that gives n_peak_min
    under the reference conversion. Applied identically for every conversion."""
    return n_peak_min / nh_per_tau(reference)
