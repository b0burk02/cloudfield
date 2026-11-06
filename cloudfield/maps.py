"""Per-pixel G0 / alphaG map computation, statistics, FITS output and figure.

Edenhofer only: LightCube's I_UV was computed through the Edenhofer map, so n
must come from it too. Any other dust source raises.
"""
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import spearmanr

from .constants import PC_CM, ALPHA_PREF, ISRF_BUNIT, R_NO_DUST
from .cubes import CubeError, open_cubes, origin_index, los_xyz, sample_los, rng, solar_mean, \
    check_isrf_units, voxel_snr
from .conversions import (NH_CONVERSIONS, A1565_AV, RV, nh_per_tau, nh_per_E,
                          conversion_summary, tau_peak_min, dust_extinction_version)
from .grid import SkyGrid, radial_grid, resolution, pixel_deg_from_config, pct
from .regions import build_region
from .boundaries import (RayContext, BoundaryResult, build_boundary, METHODS, METHOD_CODES,
                         FLAG_NAMES)
from .paths import data_path, data_setting, config_file

MAP_NAMES = ["G0_surface", "G0_peak", "n_surface", "n_peak", "alphaG_surface",
             "alphaG_peak", "NH_cloud", "NH_full", "r0", "rpeak", "r1",
             "thickness", "alphaG_losmean", "G0_trend_ratio"]


def require_edenhofer(src):
    if src != "edenhofer":
        raise ValueError(f"dust source {src!r} is not allowed for maps: LightCube's I_UV was "
                         "computed through the Edenhofer map, so n must come from it too. "
                         "DUSTRIBUTION is for column-density comparison only (cloudfield.columns).")


def _merge(primary, fallback, res_p, res_f):
    """Use the fallback where the primary method is not valid."""
    use = ~res_p.valid
    pick = lambda a, b: np.where(use, b, a)
    merged = BoundaryResult(pick(res_p.ipk, res_f.ipk), pick(res_p.i0, res_f.i0),
                            pick(res_p.i1, res_f.i1), pick(res_p.has0, res_f.has0),
                            pick(res_p.has1, res_f.has1), res_p.valid | res_f.valid,
                            pick(res_p.rpeak_edge, res_f.rpeak_edge),
                            res_p.flags | np.where(use, res_f.flags, 0).astype(np.uint16))
    code = np.where(use, np.where(res_f.valid, fallback.code, 0), primary.code)
    return merged, code


class Run:
    """Everything needed to compute one cloud."""

    def __init__(self, cloud_cfg, conversion=None, verbose=True):
        self.cfg = c = cloud_cfg
        require_edenhofer(c["dust_source"])
        self.conversion = conversion or c["conversion"]
        if self.conversion not in NH_CONVERSIONS:
            raise ValueError(f"unknown conversion {self.conversion!r}; options {list(NH_CONVERSIONS)}")
        self.nh_per_tau = nh_per_tau(self.conversion)
        q = c["quality"]
        self.tau_peak_min = tau_peak_min(q["n_peak_min"], q["reference_conversion"])
        g = c["grid"]
        self.grid = SkyGrid(tuple(g["l"]), tuple(g["b"]), pixel_deg_from_config(g["pixel"]))
        self.grid.check_wcs()
        rd = c["radial"]
        self.r = radial_grid(rd["r_min"], rd["r_max"], rd["dr"])
        self.dr = rd["dr"]
        self.boundary = build_boundary(c["boundary"])
        fb = c["boundary"].get("fallback")
        self.fallback = build_boundary(fb) if fb else None
        rres = c["resolution"]
        dist = c.get("distance_pc", [rd["r_min"], rd["r_max"]])
        self.res = resolution(self.grid.pix, 0.5 * (dist[0] + dist[1]), rres["voxel_pc"],
                              rres.get("native_ang_deg"), warn=verbose)
        self.verbose = verbose

    # ------------------------------------------------------------------
    def open(self):
        isrf_path, dust_path = data_path("lightcube"), data_path("edenhofer_tau")
        clip = data_setting("cube_radius_pc")
        if clip is not None and self.r[-1] > float(clip):
            raise CubeError(f"radial sampling reaches {self.r[-1]:g} pc but the cubes are clipped "
                            f"at r = {clip} pc (cube_radius_pc in {config_file()}); samples beyond "
                            "it would silently read zeros. Use the unclipped cubes (data/README.md).")
        self.isrf_mm, self.dust_mm, self.hdr = open_cubes(isrf_path, dust_path)
        solar, nvox = solar_mean(self.isrf_mm, self.hdr, origin_mode=self.cfg["origin_mode"])
        check_isrf_units(self.hdr, solar)
        self.solar = solar
        self.origin = origin_index(self.isrf_mm.shape, self.hdr, self.cfg["origin_mode"])
        return self

    # ------------------------------------------------------------------
    def compute(self):
        """Per-pixel maps over the full grid. Geometry and the quality cut are
        decided on the dust quantity D (tau_1565/pc); n = D * nh_per_tau is used
        only for values."""
        grid, r, DR = self.grid, self.r, self.dr
        NB, NL, NR = grid.NB, grid.NL, r.size
        c = self.cfg
        lw = c["losmean_window_pc"]
        w_los = (r >= lw[0]) & (r <= lw[1]) & (r >= R_NO_DUST)
        rc = r - r.mean()
        M = {k: np.full((NB, NL), np.nan) for k in MAP_NAMES}
        Qn = np.zeros((NB, NL), bool)
        Qr = np.zeros((NB, NL), bool)
        IPK = np.full((NB, NL), -1, int)
        FLAGS = np.zeros((NB, NL), np.uint16)
        METHOD = np.zeros((NB, NL), np.uint8)
        pick = lambda a, i: a[np.arange(a.shape[0]), i]
        CH = c["chunk_rows"]
        for j0 in range(0, NB, CH):
            j1 = min(j0 + CH, NB)
            BB, LL = np.meshgrid(grid.bvals[j0:j1], grid.lvals, indexing="ij")
            X, Y, Z = los_xyz(LL.ravel(), BB.ravel(), r)
            first = j0 == 0 and self.verbose
            D, I = sample_los(self.isrf_mm, self.dust_mm, self.origin, X, Y, Z, pad=c["pad"],
                              verbose=rng if first else None)            # tau_1565/pc
            scale, d_min = self.nh_per_tau, self.tau_peak_min
            D[:, r < R_NO_DUST] = np.nan
            n = D * scale
            if first:
                rng("n [edenhofer]", n)
            P = n.shape[0]

            ctx = RayContext(LL.ravel(), BB.ravel(), r, X, Y, Z)
            res = self.boundary.find(D, ctx)
            code = np.full(P, self.boundary.code)
            if self.fallback is not None:
                res, code = _merge(self.boundary, self.fallback, res, self.fallback.find(D, ctx))
            ipk, i0, i1, has0, has1 = res.ipk, res.i0, res.i1, res.has0, res.has1
            dpk = pick(D, ipk)
            npk = pick(n, ipk)

            ok = (dpk >= d_min) & res.valid
            s0, s1 = ok & has0, ok & has1
            nan = np.full(P, np.nan)
            G0s = np.where(s0, pick(I, i0), nan)
            G0p = np.where(ok, pick(I, ipk), nan)
            ns = np.where(s0, pick(n, i0), nan)
            npk_m = np.where(ok, npk, nan)
            C = np.cumsum(np.nan_to_num(n), axis=1) * DR * PC_CM
            nw = np.where(w_los[None, :], n, 0.0)
            sw = nw.sum(1)
            with np.errstate(invalid="ignore", divide="ignore"):
                I_n = (nw * I).sum(1) / sw
                n_n = (nw * nw).sum(1) / sw
            slope = (I * rc[None, :]).sum(1) / (rc ** 2).sum()
            icpt = I.mean(1)
            trend = (icpt + slope * (r[i0] - r.mean())) / (icpt + slope * (r[ipk] - r.mean()))
            r0 = np.where(s0, r[i0], nan)
            r1 = np.where(s1, r[i1], nan)
            vals = dict(G0_surface=G0s, G0_peak=G0p, n_surface=ns, n_peak=npk_m,
                        alphaG_surface=ALPHA_PREF * G0s / ns,
                        alphaG_peak=ALPHA_PREF * G0p / npk_m,
                        NH_cloud=np.where(s0 & s1, pick(C, i1) - pick(C, i0), nan),
                        NH_full=np.where(ok, C[:, -1], nan), r0=r0,
                        rpeak=np.where(ok, r[ipk], nan), r1=r1, thickness=r1 - r0,
                        alphaG_losmean=np.where(ok, ALPHA_PREF * I_n / n_n, nan),
                        G0_trend_ratio=np.where(s0, trend, nan))
            sh = (j1 - j0, NL)
            for k in MAP_NAMES:
                M[k][j0:j1] = vals[k].reshape(sh)
            Qn[j0:j1] = ok.reshape(sh)
            Qr[j0:j1] = (ok & ~res.rpeak_edge).reshape(sh)
            IPK[j0:j1] = np.where(ok, ipk, -1).reshape(sh)
            FLAGS[j0:j1] = res.flags.reshape(sh)
            METHOD[j0:j1] = np.where(ok, code, 0).astype(np.uint8).reshape(sh)
        M["ratio_surface_losmean"] = M["alphaG_surface"] / M["alphaG_losmean"]
        if self.verbose:
            print(f"  done: quality cut passed {Qn.sum()}/{NB*NL}, r0 found {np.isfinite(M['r0']).sum()}")
        self.out = dict(M=M, Qn=Qn, Qr=Qr, IPK=IPK, FLAGS=FLAGS, METHOD=METHOD)
        return self.out

    # ------------------------------------------------------------------
    def build_region(self):
        spots = [(p["name"], (p["l"], p["b"])) for p in self.cfg["annotations"]]
        self.region = build_region(self.cfg["region"], self.grid, spots=spots)
        return self.region

    def noise_snr(self):
        c = self.cfg
        cen = c.get("noise_block_xyz_pc")
        if cen is None:     # default: centre of the region at the middle of the distance range
            lb = np.radians([np.mean(self.region.bbox["l"]), np.mean(self.region.bbox["b"])])
            d = 0.5 * sum(c.get("distance_pc", [self.r[0], self.r[-1]]))
            cen = [d * np.cos(lb[1]) * np.cos(lb[0]), d * np.cos(lb[1]) * np.sin(lb[0]),
                   d * np.sin(lb[1])]
        self.snr = voxel_snr(self.isrf_mm, self.origin, cen, c["noise_block_half_pc"])
        return self.snr

    # ------------------------------------------------------------------
    def report(self):
        M, Qn, Qr, IPK = self.out["M"], self.out["Qn"], self.out["Qr"], self.out["IPK"]
        reg = self.region.mask
        g = self.grid
        S = {}
        print("=" * 78)
        print(f"[{self.cfg['name']}] statistics inside region '{self.region.extname}' ({reg.sum()} px)")
        win = getattr(self.boundary, "window", None)
        if win is not None:
            sel = Qn & reg
            rp = M["rpeak"][sel]
            edges = np.arange(win[0], win[1] + 5.0, 5.0)
            h_, _ = np.histogram(rp, bins=edges)
            print("  r_peak histogram (5 pc bins):")
            for a, b_, cc in zip(edges[:-1], edges[1:], h_):
                print(f"    {a:5.0f}-{b_:5.0f} pc  {cc:6d}  {cc/max(rp.size,1):6.1%}  "
                      + "#" * int(50 * cc / max(h_.max(), 1)))
        print(f"  {'map':<22} {'median':>10} {'p5':>10} {'p95':>10}  N")
        for k, v in M.items():
            (p50, p5, p95), nf = pct(v, reg)
            print(f"  {k:<22} {p50:10.4g} {p5:10.4g} {p95:10.4g}  {nf}")
        print(f"  quality inside region: peak cut passed {Qn[reg].mean():.1%} ({Qn.sum()} px full "
              f"field); r_peak not on search edge {Qr[reg].mean():.1%}")
        fl = self.out["FLAGS"][reg & Qn]
        print("  boundary flags inside region (fraction of passing pixels): " +
              ", ".join(f"{n} {np.mean((fl & b) > 0):.1%}" for b, n in sorted(FLAG_NAMES.items())
                        if np.any(fl & b)) or "none")
        mc = self.out["METHOD"][reg & Qn]
        print("  boundary method per pixel: " + ", ".join(
            f"{n} {np.mean(mc == code):.1%}" for n, code in METHOD_CODES.items() if np.any(mc == code)))
        gg = reg & np.isfinite(M["alphaG_surface"]) & np.isfinite(M["alphaG_peak"])
        gr = (M["G0_surface"] / M["G0_peak"])[gg]
        hf = (M["n_surface"] / (M["n_peak"] / 2))[gg]
        S.update(GR=np.median(gr), GR16=np.percentile(gr, 16), GR84=np.percentile(gr, 84),
                 TR=np.nanmedian(M["G0_trend_ratio"][gg]), HF=np.median(hf))
        print(f"  n_surface/(n_peak/2) median {S['HF']:.4f};  G0_surface/G0_peak median "
              f"{S['GR']:.4f} (16-84% {S['GR16']:.3f}-{S['GR84']:.3f}); radial trend predicts "
              f"{S['TR']:.4f}; voxel SNR ~{self.snr:.1f} -> ratio noise {np.sqrt(2)/self.snr:.3f}")
        ratio = M["ratio_surface_losmean"]
        gd = np.isfinite(ratio) & reg
        a, b_, cc = np.percentile(ratio[gd], [16, 50, 84])
        rho, _ = spearmanr(np.log10(ratio[gd]), M["thickness"][gd])
        print(f"  alphaG_surface/alphaG_losmean median {b_:.3f} (16-84% {a:.3f}-{cc:.3f}); "
              f"Spearman vs thickness {rho:+.3f}")
        S["validation"] = []
        for vp in self.cfg["validation_points"]:
            im = int(np.argmin(np.abs(g.lvals - vp["l"])))
            jm = int(np.argmin(np.abs(g.bvals - vp["b"])))
            sep = SkyCoord(g.LL2 * u.deg, g.BB2 * u.deg, frame="galactic").separation(
                SkyCoord(vp["l"] * u.deg, vp["b"] * u.deg, frame="galactic")).deg
            ap = sep <= vp.get("radius_deg", 0.5)
            print(f"  {vp['name']} nearest pixel ({g.lvals[im]:.1f},{g.bvals[jm]:.1f}), "
                  f"aperture {ap.sum()} px:")
            for k in ["n_peak", "n_surface", "G0_surface", "alphaG_peak", "alphaG_surface",
                      "rpeak", "NH_full"]:
                print(f"    {k:<16} nearest {M[k][jm, im]:10.4g}   median(ap) {np.nanmedian(M[k][ap]):10.4g}")
            if "published_n_H" in vp:
                d_px = vp["published_n_H"] / M["n_peak"][jm, im]
                d_ap = vp["published_n_H"] / np.nanmedian(M["n_peak"][ap])
                S["validation"].append(dict(name=vp["name"], dil_px=d_px, dil_ap=d_ap,
                                            pub=vp["published_n_H"]))
                print(f"  * {vp['name']} dilution {vp['published_n_H']:g}/n_peak: nearest "
                      f"{d_px:.2f}, aperture {d_ap:.2f}  ({vp.get('reference', '')})")
        self.stats = S
        return S

    # ------------------------------------------------------------------
    def header_lines(self):
        conv = NH_CONVERSIONS[self.conversion]
        hist = [
            "Edenhofer+2024 dust is used for BOTH I_UV and n: LightCube's RT",
            "(Bish+2026) was computed through the Edenhofer map.",
            f"n_H conversion '{self.conversion}': {self.nh_per_tau:.4g} cm^-3 per tau_1565/pc",
            f"= {nh_per_E(self.conversion):.5g} per E/pc. E->A_V x2.8 (Zhang+2023);",
            f"A1565/AV = {A1565_AV:.4f} from Gordon+2023 (G23, R_V={RV:g}),",
            f"dust_extinction {dust_extinction_version()}.",
        ] + conv.get("caveat", [conv["source"][:72]]) + [
            f"Quality cut in dust units: tau_1565,peak >= {self.tau_peak_min:.5g} /pc",
            f"(= n_peak {self.cfg['quality']['n_peak_min']:g} cm^-3 under "
            f"'{self.cfg['quality']['reference_conversion']}'); same mask for every",
            "conversion."]
        rd = self.cfg["radial"]
        comments = self.boundary.describe() + (
            ["fallback where primary invalid:"] + self.fallback.describe() if self.fallback else []) + [
            f"Pixels with tau_1565,peak < {self.tau_peak_min:.4g}/pc are NaN (conversion-indep.).",
            "ISRF: LightCube (Bish et al. 2026), Draine (1978) normalized 912-2000 A.",
            f"Rays from Sun r={rd['r_min']:g}-{rd['r_max']:g} pc, dr={rd['dr']:g} pc.",
            f"n: Edenhofer+2024 tau_1565/pc x {self.nh_per_tau:.4g} cm^-3, total H nuclei",
            f"({nh_per_E(self.conversion):.5g} cm^-3 per E/pc; conversion '{self.conversion}').",
            "r < 69 pc masked (no dust in map). Same dust as LightCube's RT."]
        return comments, hist

    def write_fits(self, out_dir):
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        pre = self.cfg["output_prefix"]
        names = dict(G0=out_dir / f"{pre}_G0.fits", alpha=out_dir / f"{pre}_alphaG.fits",
                     maps=out_dir / f"{pre}_maps.fits")
        M, Qn, Qr, S = self.out["M"], self.out["Qn"], self.out["Qr"], self.stats
        w = self.grid.wcs()
        comments, hist = self.header_lines()
        UNIT_NOTE = {"cm-3": "H nuclei per cm^3", "cm-2": "H nuclei per cm^2"}
        name = self.cfg["name"]

        def base_header(extname):
            h = w.to_header()
            h["EXTNAME"] = extname
            h["CLOUD"] = (name, "cloudfield cloud config")
            h["DUSTSRC"] = ("edenhofer", "dust for both I_UV (LightCube RT) and n")
            h["NHCONV"] = (self.conversion, "Edenhofer dust-to-gas conversion")
            h["NHPERTAU"] = (round(self.nh_per_tau, 4), "n_H [cm-3] per tau_1565/pc")
            h["BNDMETH"] = (self.boundary.name, "primary boundary method")
            for k, v, cm in self.res.header_cards():
                h[k] = (v, cm)
            return h

        def hdu(data, extname, bunit, extra=(), primary=False):
            h = base_header(extname)
            if bunit in UNIT_NOTE:
                h["BUNIT"] = (bunit, f"{UNIT_NOTE[bunit]}, n_H via {self.conversion}")
            elif extname.startswith("ALPHAG"):
                h["BUNIT"] = (bunit, f"dimensionless; n_H via {self.conversion}")
            else:
                h["BUNIT"] = bunit
            for cc in list(extra) + comments:
                h.add_comment(cc)
            for cc in hist:
                h.add_history(cc)
            return (fits.PrimaryHDU if primary else fits.ImageHDU)(
                data=np.asarray(data, np.float32), header=h)

        G0_NOTE = ["G0 = I_UV in units of the Draine (1978) ISRF, 912-2000 A,",
                   "taken directly from the LightCube cube values. Per-pixel values",
                   "are not independent measurements (LightCube voxel noise); use",
                   "region medians."]
        fits.HDUList([hdu(M["G0_surface"], "G0_SURFACE", ISRF_BUNIT,
                          G0_NOTE + ["G0_SURFACE = I_UV at r0 (incident field)."], primary=True),
                      hdu(M["G0_peak"], "G0_PEAK", ISRF_BUNIT, G0_NOTE + ["G0_PEAK = I_UV at r_peak."])
                      ]).writeto(names["G0"], overwrite=True)
        below = abs(S["GR"] - 1) < np.sqrt(2) / self.snr
        A_NOTE = [
            "alphaG = 59 * I_UV / n (Draine-normalized I_UV, n = total H nuclei;",
            "Sternberg+2014, Bialy & Sternberg 2016).",
            "ALPHAG_PEAK (primary) uses I_UV(r_peak) and n_peak.",
            "ALPHAG_SURFACE uses I_UV(r0), n(r0).",
            f"G0_surface/G0_peak median {S['GR']:.3f} (16-84% {S['GR16']:.3f}-{S['GR84']:.3f})",
            f"in region, {'below' if below else 'ABOVE'} LightCube voxel noise",
            f"(SNR ~{self.snr:.0f}); radial I_UV trend alone predicts {S['TR']:.3f}.",
            ("No attenuation detected r0 -> r_peak." if below else "Attenuation r0->r_peak exceeds noise."),
        ]
        if self.boundary.name == "half_max":
            A_NOTE += [f"half_max: n_surface ~ n_peak/2 by construction (median {S['HF']:.3f}),",
                       "so alphaG_surface ~ 2 alphaG_peak G0_surface/G0_peak."]
        for v in S["validation"]:
            A_NOTE += [f"{v['name']}: published n_H={v['pub']:g} cm^-3 is {v['dil_px']:.1f}x the",
                       "mapped n_peak: the dust map does not resolve the condensation;",
                       "alphaG in clumped gas is an upper limit."]
        fits.HDUList([hdu(M["alphaG_peak"], "ALPHAG_PEAK", "", A_NOTE, primary=True),
                      hdu(M["alphaG_surface"], "ALPHAG_SURFACE", "", A_NOTE)
                      ]).writeto(names["alpha"], overwrite=True)
        rd = self.cfg["radial"]
        ext = [("N_SURFACE", "n_surface", "cm-3", "n at r0"),
               ("N_PEAK", "n_peak", "cm-3", "n at r_peak"),
               ("NH_CLOUD", "NH_cloud", "cm-2", "integral of n from r0 to r1 (not whole slab)"),
               ("NH_FULL", "NH_full", "cm-2", f"integral of n over {rd['r_min']:g}-{rd['r_max']:g} pc"),
               ("R0", "r0", "pc", "near-side boundary radius"),
               ("RPEAK", "rpeak", "pc", "radius of n_peak"),
               ("R1", "r1", "pc", "far-side boundary radius"),
               ("THICKNESS", "thickness", "pc", "r1 - r0"),
               ("ALPHAG_LOSMEAN", "alphaG_losmean", "", "OLD estimator 59<I>_n/<n>_n, comparison only"),
               ("RATIO_SURF_LOSMEAN", "ratio_surface_losmean", "", "alphaG_surface / alphaG_losmean")]
        hl = [fits.PrimaryHDU()]
        hl[0].header["CLOUD"] = name
        hl[0].header["DUSTSRC"] = "edenhofer"
        hl[0].header.add_comment(f"{name} per-pixel boundary maps; see extensions.")
        for e, k, un, cc in ext:
            hl.append(hdu(M[k], e, un, [cc]))
        qcut = f"1 if tau_1565,peak >= {self.tau_peak_min:.5g}/pc (conversion-independent)"
        for e, q, cc in [("QUAL_NPEAK", Qn, qcut),
                         ("QUAL_RPEAK_INTERIOR", Qr, "1 if r_peak not on the search-range edge"),
                         (self.region.extname, self.region.mask, f"1 inside: {self.region.description}"[:72])]:
            h = base_header(e)
            h["BUNIT"] = ""
            h.add_comment(cc)
            for hc in hist:
                h.add_history(hc)
            hl.append(fits.ImageHDU(q.astype(np.uint8), header=h))
        h = base_header("BOUNDARY_FLAGS"); h["BUNIT"] = ""
        for b, n in sorted(FLAG_NAMES.items()):
            h[f"FLAG{int(np.log2(b))}"] = (n, f"bit value {b}")
        hl.append(fits.ImageHDU(self.out["FLAGS"].astype(np.uint16), header=h))
        h = base_header("BOUNDARY_METHOD"); h["BUNIT"] = ""
        for n, code in METHOD_CODES.items():
            h[f"METH{code}"] = (n, "method code")
        hl.append(fits.ImageHDU(self.out["METHOD"], header=h))
        fits.HDUList(hl).writeto(names["maps"], overwrite=True)
        print(f"  wrote {names['G0']}, {names['alpha']}, {names['maps']}")
        return names

    # ------------------------------------------------------------------
    def figure(self, out_dir):
        from .plotting import map_figure
        return map_figure(self, Path(out_dir) / f"{self.cfg['output_prefix']}_maps.png")


def run_cloud(cloud_cfg, out_dir, conversion=None, make_figure=True):
    R = Run(cloud_cfg, conversion=conversion)
    print(f"cloud {cloud_cfg['name']}: {conversion_summary(R.conversion)}")
    print(f"grid [b,l]=({R.grid.NB},{R.grid.NL}) pixel {R.grid.pix:g} deg; resolution at "
          f"{R.res.distance_pc:g} pc: pixel {R.res.pixel_pc:.3g} pc, dust effective "
          f"{R.res.effective_pc:.3g} pc, oversampling x{R.res.oversampling:.2f}")
    R.open()
    print(f"LightCube solar check: mean I_UV within 25 pc = {R.solar:.4f} (expected 0.633)")
    print("=" * 78)
    R.build_region()
    R.noise_snr()
    R.compute()
    R.report()
    names = R.write_fits(out_dir)
    if make_figure:
        R.figure(out_dir)
    return R, names
