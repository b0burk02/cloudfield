"""Data-path configuration. Nothing is hardcoded: paths come from a YAML file.

Resolution order for the config file:
  1. $CLOUDFIELD_CONFIG
  2. ./config.yaml
Relative paths inside it are resolved against $CLOUDFIELD_DATA if set,
otherwise against the config file's directory.

See config.example.yaml and data/README.md.
"""
import os
from pathlib import Path
import yaml

DATA_README = "data/README.md"

# key -> (description, required for `cloudfield run`)
DATA_KEYS = {
    "lightcube": ("LightCube combined 912-2000 A cube, Draine-normalized (Bish+2026)", True),
    "edenhofer_tau": ("Edenhofer+2024 dust as tau_1565 per pc on the LightCube grid", True),
    "dustribution_density": ("DUSTRIBUTION median density cube (column comparison only)", False),
    "dustribution_cumulative": ("DUSTRIBUTION cumulative extinction cube (column comparison only)", False),
    "fims_h2": ("FIMS/SPEAR H2/FUV ratio HEALPix map (Eos region / validation)", False),
    "galfa_nhi": ("GALFA-HI integrated N_HI map (validation scripts)", False),
    "galfa_cube": ("GALFA-HI PPV cube (validation scripts)", False),
    "galfa_check": ("GALFA-HI Galactic-frame column map (frame check, validation scripts)", False),
}


class DataPathError(FileNotFoundError):
    pass


def config_file():
    p = os.environ.get("CLOUDFIELD_CONFIG")
    if p:
        return Path(p)
    return Path("config.yaml")


def load_data_config():
    cf = config_file()
    if not cf.exists():
        raise DataPathError(
            f"No data config found at {cf}. Copy config.example.yaml to config.yaml "
            f"(or set CLOUDFIELD_CONFIG) and fill in the paths; see {DATA_README}.")
    with open(cf) as f:
        cfg = yaml.safe_load(f) or {}
    base = Path(os.environ.get("CLOUDFIELD_DATA", cf.resolve().parent))
    out = {}
    for k, v in (cfg.get("data") or {}).items():
        if v is None:
            continue
        p = Path(os.path.expanduser(str(v)))
        out[k] = p if p.is_absolute() else (base / p)
    return out


def data_setting(key, default=None):
    """Non-path settings in the data config (e.g. cube_radius_pc for clipped cubes)."""
    cf = config_file()
    if not cf.exists():
        return default
    with open(cf) as f:
        cfg = yaml.safe_load(f) or {}
    return cfg.get(key, default)


def data_path(key, required=True):
    """Absolute path for a data key; raises DataPathError naming the file if missing."""
    paths = load_data_config()
    desc = DATA_KEYS.get(key, (key, required))[0]
    p = paths.get(key)
    if p is None:
        if required:
            raise DataPathError(f"data key '{key}' ({desc}) is not set in {config_file()}; "
                                f"see {DATA_README}.")
        return None
    if not p.exists():
        if required:
            raise DataPathError(f"missing file for '{key}' ({desc}): {p}; see {DATA_README}.")
        return None
    return p
