"""`cloudfield check-data`: preflight for every expected input.

The single most consequential failure mode is pointing the code at a per-band
LightCube cube instead of the combined, Draine-normalized one: alphaG comes out
wrong by orders of magnitude and nothing else raises. This check catches it via
BUNIT and the 0.633 solar-neighbourhood test.
"""
import numpy as np
from astropy.io import fits

from .constants import ISRF_BUNIT, SOLAR_I_UV
from .paths import DATA_KEYS, load_data_config, config_file, DATA_README
from .cubes import solar_mean, origin_index


def _status(ok, msg):
    return ("PASS" if ok else "FAIL"), msg


def check_lightcube(p):
    with fits.open(p, memmap=True) as hl:
        h, d = hl[0].header, hl[0].data
        bunit = str(h.get("BUNIT", "")).strip()
        msgs = [f"shape {d.shape}, dtype {d.dtype}, BUNIT {bunit!r}"]
        if bunit != ISRF_BUNIT:
            return _status(False, "; ".join(msgs) + f" -- expected {ISRF_BUNIT!r} (combined "
                           "912-2000 A Draine-normalized cube, not a TD1 band cube)")
        m, nv = solar_mean(d, h)
        ok = abs(m / SOLAR_I_UV - 1) < 0.02
        o = origin_index(d.shape, h, "header")
        msgs.append(f"origin {o}; solar mean (r<25 pc, {nv} voxels) = {m:.4f} vs {SOLAR_I_UV}")
        return _status(ok, "; ".join(msgs))


def check_edenhofer(p, lightcube_path):
    with fits.open(p, memmap=True) as hl:
        d = hl[0].data
        msg = f"shape {d.shape}, dtype {d.dtype}"
        shape_ok = True
        if lightcube_path is not None:
            with fits.open(lightcube_path, memmap=True) as hl2:
                shape_ok = hl2[0].data.shape == d.shape
                msg += f"; LightCube shape {hl2[0].data.shape}"
        c = [s // 2 for s in d.shape]
        v = float(np.asarray(d[c[0], c[1], c[2]]))
        sample = np.asarray(d[c[0] + 100, c[1], c[2]:c[2] + 50], dtype=float)
        msg += f"; value at origin {v:.3g} (Edenhofer: ~0 inside 69 pc); tau/pc sample max {sample.max():.3g}"
        return _status(shape_ok and np.isfinite(sample).all(), msg)


def check_healpix(p):
    with fits.open(p) as hl:
        hx = next((h for h in hl if h.header.get("PIXTYPE") == "HEALPIX"), None)
        if hx is None:
            return _status(False, "no HEALPix extension")
        h = hx.header
        return _status(True, f"Nside {h['NSIDE']}, {h['ORDERING']}, COORDSYS {h.get('COORDSYS')!r}")


def check_fits_header(p):
    with fits.open(p, memmap=True) as hl:
        h = hl[0].header
        return _status(True, f"NAXIS {h.get('NAXIS')}, BUNIT {h.get('BUNIT')!r}, "
                             f"CTYPE {h.get('CTYPE1')!r}/{h.get('CTYPE2')!r}")


def check_dustribution_density(p):
    import zipfile
    with zipfile.ZipFile(p) as z:
        names = [n for n in z.namelist() if n.endswith("/data/0")]
        sz = z.getinfo(names[0]).file_size if names else 0
    from .columns import DUSTRIBUTION_SHAPE
    ok = sz == 8 * np.prod(DUSTRIBUTION_SHAPE)
    return _status(ok, f"torch zip, storage {sz} bytes (expect {8*np.prod(DUSTRIBUTION_SHAPE)} for "
                       f"{DUSTRIBUTION_SHAPE} float64)")


def check_dustribution_cumulative(p):
    a = np.load(p, mmap_mode="r")
    return _status(a.shape == (360, 180, 118), f"shape {a.shape} (expect (360, 180, 118))")


CHECKS = {"lightcube": check_lightcube, "fims_h2": check_healpix,
          "galfa_nhi": check_fits_header, "galfa_cube": check_fits_header,
          "galfa_check": check_fits_header,
          "dustribution_density": check_dustribution_density,
          "dustribution_cumulative": check_dustribution_cumulative}


def check_data():
    print(f"data config: {config_file()}  (see {DATA_README})")
    paths = load_data_config()
    worst = 0
    for key, (desc, required) in DATA_KEYS.items():
        p = paths.get(key)
        tag = "required" if required else "optional"
        if p is None or not p.exists():
            st = "MISSING" if required else "absent"
            print(f"  [{st:7}] {key:<24} ({tag}) {desc}" + (f": {p}" if p else ": not configured"))
            worst = max(worst, 2 if required else 0)
            continue
        try:
            if key == "edenhofer_tau":
                st, msg = check_edenhofer(p, paths.get("lightcube"))
            else:
                st, msg = CHECKS.get(key, check_fits_header)(p)
        except Exception as e:      # report, don't hide
            st, msg = "FAIL", f"{type(e).__name__}: {e}"
        print(f"  [{st:7}] {key:<24} ({tag}) {p}\n            {msg}")
        if st == "FAIL":
            worst = max(worst, 2 if required else 1)
    print("result:", {0: "OK", 1: "OK (optional inputs have problems)", 2: "FAILED"}[worst])
    return worst
