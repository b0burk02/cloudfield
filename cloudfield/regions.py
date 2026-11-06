"""Sky regions: where statistics are quoted and what is overlaid on figures.

The map is always computed over the full l/b grid; the region only selects the
pixels for quoted statistics (and is written to FITS as a 0/1 extension).

Region types (cloud config `region`):
  {type: lb_box, l: [l0, l1], b: [b0, b1]}
  {type: healpix_contour, data_key: <paths key>, column: SIGNAL, level: 10.0,
   pick_box: {l: [..], b: [..]}, label: "...", context_levels: [...]}
     -> the connected component of (map >= level) with most pixels inside pick_box.
  {type: all}
"""
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import label

from .cubes import rng
from .paths import data_path


def healpix_on_grid(path, column, lvals, bvals, spot=()):
    """Sample a HEALPix FITS map onto the Galactic (l, b) grid.
    Frame, ordering and indexing are read from the header. Bilinear
    interpolation (hp.get_interp_val); bad pixels -> NaN."""
    import healpy as hp
    with fits.open(path) as hl:
        ext = next(i for i, h in enumerate(hl) if h.header.get("PIXTYPE") == "HEALPIX")
        h, d = hl[ext].header, hl[ext].data
        nside = int(h["NSIDE"])
        nest = h["ORDERING"].strip().upper() == "NESTED"
        frame = h["COORDSYS"].strip().upper()[0]
        sig = np.asarray(d[column], dtype=np.float64).ravel()
        if h.get("INDXSCHM", "IMPLICIT").strip().upper() == "EXPLICIT":
            pix = np.asarray(d["PIXEL"]).ravel()
        else:
            pix = np.arange(sig.size)
    full = np.full(hp.nside2npix(nside), np.nan)
    full[pix] = sig
    bad = h.get("BAD_DATA", -1.6375e30)
    full[(full <= 0.5 * bad) | ~np.isfinite(full)] = np.nan
    print(f"  HEALPix {path}: HDU {ext}, Nside={nside}, "
          f"{'NESTED' if nest else 'RING'}, COORDSYS={h['COORDSYS']!r} -> frame {frame}")

    def to_frame(l, b):
        g = SkyCoord(l=np.ravel(l) * u.deg, b=np.ravel(b) * u.deg, frame="galactic")
        if frame == "G":
            return g.l.deg, g.b.deg
        if frame == "C":
            c = g.icrs
            return c.ra.deg, c.dec.deg
        if frame == "E":
            c = g.barycentrictrueecliptic
            return c.lon.deg, c.lat.deg
        raise ValueError(f"Unknown HEALPix COORDSYS {h['COORDSYS']!r}")

    print(f"  frame conversion applied: Galactic (l,b) -> "
          f"{ {'G': 'none (map is Galactic)', 'C': 'ICRS (RA,Dec)', 'E': 'ecliptic'}[frame] }")
    for name, (l, b) in spot:
        lon, lat = to_frame([l], [b])
        ip = hp.ang2pix(nside, lon[0], lat[0], lonlat=True, nest=nest)
        plon, plat = hp.pix2ang(nside, ip, lonlat=True, nest=nest)
        if frame == "G":
            gl, gb = plon, plat
        else:
            src = dict(C="icrs", E="barycentrictrueecliptic")[frame]
            gg = SkyCoord(plon * u.deg, plat * u.deg, frame=src).galactic
            gl, gb = gg.l.deg, gg.b.deg
        sep = SkyCoord(l * u.deg, b * u.deg, frame="galactic").separation(
            SkyCoord(gl * u.deg, gb * u.deg, frame="galactic")).deg
        print(f"    spot {name:<12} (l,b)=({l:.2f},{b:.2f}) -> map ({lon[0]:.2f},{lat[0]:.2f})"
              f" -> pix {ip} centre (l,b)=({gl:.2f},{gb:.2f}), sep {sep:.2f} deg "
              f"(pixel {hp.nside2resol(nside, arcmin=True)/60:.2f} deg), value {full[ip]:.4g}")
    LL, BB = np.meshgrid(lvals, bvals)
    lon, lat = to_frame(LL, BB)
    out = hp.get_interp_val(full, lon, lat, lonlat=True, nest=nest).reshape(LL.shape)
    ipx = hp.ang2pix(nside, lon, lat, lonlat=True, nest=nest)
    print(f"  raw HEALPix pixels in field: max {np.nanmax(full[ipx]):.4g}, "
          f"min {np.nanmin(full[ipx]):.4g}")
    rng("HEALPix on grid", out)
    return out


class Region:
    """mask: bool (NB, NL); overlay: optional dict for figures
    (map, level, context_levels, label); extname: FITS extension name."""

    def __init__(self, mask, extname="REGION", description="", overlay=None, bbox=None):
        self.mask, self.extname, self.description = mask, extname, description
        self.overlay, self.bbox = overlay, bbox


def build_region(cfg, grid, spots=()):
    t = cfg.get("type", "all")
    LL2, BB2 = grid.LL2, grid.BB2
    if t == "all":
        return Region(np.ones(grid.shape, bool), cfg.get("extname", "REGION"), "full grid",
                      bbox=dict(l=list(grid.l_range), b=list(grid.b_range)))
    if t == "lb_box":
        m = ((LL2 >= cfg["l"][0]) & (LL2 <= cfg["l"][1]) &
             (BB2 >= cfg["b"][0]) & (BB2 <= cfg["b"][1]))
        return Region(m, cfg.get("extname", "REGION"),
                      f"l={cfg['l']}, b={cfg['b']} box", bbox=dict(l=cfg["l"], b=cfg["b"]))
    if t == "healpix_contour":
        path = data_path(cfg["data_key"])
        level = float(cfg["level"])
        H = healpix_on_grid(path, cfg.get("column", "SIGNAL"), grid.lvals, grid.bvals, spot=spots)
        lab, nlab = label(H >= level)
        pb = cfg["pick_box"]
        inbox = ((LL2 >= pb["l"][0]) & (LL2 <= pb["l"][1]) &
                 (BB2 >= pb["b"][0]) & (BB2 <= pb["b"][1]))
        counts = np.bincount(lab[inbox].ravel(), minlength=nlab + 1)
        counts[0] = 0
        m = lab == counts.argmax()
        lc, bc_ = LL2[m], BB2[m]
        print(f"  {level:g} region: {nlab} components; region = largest within pick box: "
              f"{m.sum()} px, l={lc.min():.1f}-{lc.max():.1f}, b={bc_.min():.1f}-{bc_.max():.1f}; "
              f"median inside {np.median(H[m]):.2f} vs outside {np.nanmedian(H[~m]):.2f}")
        lbl = cfg.get("label", cfg["data_key"])
        return Region(m, cfg.get("extname", "REGION"),
                      f"largest {lbl} >= {level:g} component within {pb}",
                      overlay=dict(map=H, level=level, label=lbl,
                                   context_levels=list(cfg.get("context_levels", []))),
                      bbox=pb)
    raise ValueError(f"unknown region type {t!r}")
