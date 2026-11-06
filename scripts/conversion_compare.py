"""Run the Edenhofer maps under every named gas conversion and verify that the
quality mask and all geometry (r0, r_peak, r1, G0) are bit-identical across
conversions; tabulate the conversion-dependent quantities on that common mask.
Writes no FITS.
"""
import sys
import numpy as np
from cloudfield.config import load_cloud
from cloudfield.maps import Run
from cloudfield.conversions import NH_CONVERSIONS, DEFAULT_CONVERSION, nh_per_tau, nh_per_E

cloud = load_cloud(sys.argv[1] if len(sys.argv) > 1 else "eos")
RES = {}
region = None
for k in NH_CONVERSIONS:
    R = Run(cloud, conversion=k, verbose=False).open()
    if region is None:
        region = R.build_region().mask        # conversion-independent
    RES[k] = R.compute()
print(f"quality cut (dust units): tau_1565,peak >= {R.tau_peak_min:.6g} /pc")
for k in NH_CONVERSIONS:
    print(f"  {k:<18} {nh_per_tau(k):8.2f} cm^-3 per tau/pc  {nh_per_E(k):7.1f} per E/pc  "
          f"-> equivalent n_peak cut {R.tau_peak_min * nh_per_tau(k):.4f} cm^-3")
ref = RES[DEFAULT_CONVERSION]
GEOM = ["r0", "rpeak", "r1", "thickness", "G0_surface", "G0_peak", "G0_trend_ratio"]
print(f"bit-identity vs '{DEFAULT_CONVERSION}':")
for k, Rk in RES.items():
    same = {"Qn": np.array_equal(Rk["Qn"], ref["Qn"]), "Qr": np.array_equal(Rk["Qr"], ref["Qr"]),
            "IPK": np.array_equal(Rk["IPK"], ref["IPK"]),
            "FLAGS": np.array_equal(Rk["FLAGS"], ref["FLAGS"])}
    same.update({g: Rk["M"][g].tobytes() == ref["M"][g].tobytes() for g in GEOM})
    print(f"  {k:<18} " + "  ".join(f"{a}={'Y' if b else 'N'}" for a, b in same.items()))
print(f"common mask: {ref['Qn'].sum()} px pass; {ref['Qn'][region].sum()}/{region.sum()} in region")
print(f"medians inside the region [5-95%]:")
print(f"  {'map':<15}" + "".join(f"{k:>28}" for k in NH_CONVERSIONS))
for m in ["n_peak", "n_surface", "alphaG_peak", "alphaG_surface", "NH_full", "NH_cloud", "alphaG_losmean"]:
    row = ""
    for Rk in RES.values():
        v = Rk["M"][m][region]
        v = v[np.isfinite(v)]
        p50, p5, p95 = np.percentile(v, [50, 5, 95])
        row += f"{p50:10.4g} [{p5:7.3g}-{p95:7.3g}]"
    print(f"  {m:<15}{row}")
print("pixels with alphaG < 1 (full field / region):")
for m in ["alphaG_peak", "alphaG_surface"]:
    print(f"  {m:<15}" + "".join(f"{np.sum(Rk['M'][m] < 1):>10d} / {np.sum(Rk['M'][m][region] < 1):<5d}       "
                                 for Rk in RES.values()))
