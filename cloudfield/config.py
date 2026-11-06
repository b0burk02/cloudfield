"""Cloud configuration (YAML). One file per cloud in cloudfield/clouds/ (or any path).

Required keys: name, grid, radial, boundary, region.
See cloudfield/clouds/eos.yaml for a fully annotated example and README
"Adding a new cloud".
"""
from pathlib import Path
import yaml

CLOUD_DIR = Path(__file__).parent / "clouds"

DEFAULTS = {
    "dust_source": "edenhofer",
    "conversion": "lightcube",
    "origin_mode": "header",
    "quality": {"n_peak_min": 0.3, "reference_conversion": "lightcube"},
    "losmean_window_pc": [60.0, 200.0],
    "chunk_rows": 20,
    "pad": 2,
    "noise_block_half_pc": 20,
    "annotations": [],
    "validation_points": [],
    "published": {},
    "figure": {"g0_lim": [0.30, 0.70], "margin_deg": 3.0},
    "resolution": {"voxel_pc": 1.0, "native_ang_deg": None},
}


def list_clouds():
    return sorted(p.stem for p in CLOUD_DIR.glob("*.yaml"))


def load_cloud(name_or_path):
    p = Path(name_or_path)
    if not p.suffix:
        p = CLOUD_DIR / f"{name_or_path}.yaml"
    if not p.exists():
        raise FileNotFoundError(f"cloud config {p} not found; available: {list_clouds()}")
    with open(p) as f:
        cfg = yaml.safe_load(f)
    out = {**DEFAULTS, **cfg}
    for k in ("quality", "figure", "resolution"):
        out[k] = {**DEFAULTS[k], **(cfg.get(k) or {})}
    for k in ("name", "grid", "radial", "boundary", "region"):
        if k not in out:
            raise ValueError(f"cloud config {p}: missing required key '{k}'")
    out.setdefault("output_prefix", out["name"].lower())
    out["_path"] = str(p)
    return out
