"""Command line: `cloudfield check-data | list | run <cloud>`."""
import argparse
import sys


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cloudfield", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check-data", help="verify inputs (run this first)")
    sub.add_parser("list", help="list cloud configs")
    r = sub.add_parser("run", help="compute maps for a cloud")
    r.add_argument("cloud", help="config name in cloudfield/clouds/ or a path to a YAML file")
    r.add_argument("--out", default="output", help="output directory (default: output/)")
    r.add_argument("--conversion", default=None, help="override the dust-to-gas conversion")
    r.add_argument("--no-figure", action="store_true")
    a = ap.parse_args(argv)
    if a.cmd == "check-data":
        from .check import check_data
        return 2 if check_data() == 2 else 0
    if a.cmd == "list":
        from .config import list_clouds
        print("\n".join(list_clouds()))
        return 0
    if a.cmd == "run":
        from .config import load_cloud
        from .maps import run_cloud
        run_cloud(load_cloud(a.cloud), a.out, conversion=a.conversion, make_figure=not a.no_figure)
        return 0


if __name__ == "__main__":
    sys.exit(main())
