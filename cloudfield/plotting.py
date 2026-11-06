"""Three-panel map figure: NH_full, G0_surface (fixed linear scale), alphaG_peak."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

from .conversions import A1565_AV, nh_per_E


def map_figure(run, path):
    M, reg = run.out["M"], run.region
    g, cfg = run.grid, run.cfg
    figc = cfg["figure"]
    bb = reg.bbox or dict(l=list(g.l_range), b=list(g.b_range))
    mg = figc["margin_deg"]
    fl = (bb["l"][0] - mg, bb["l"][1] + mg)
    fb = (bb["b"][0] - mg, bb["b"][1] + mg)
    view = (g.LL2 >= fl[0]) & (g.LL2 <= fl[1]) & (g.BB2 >= fb[0]) & (g.BB2 <= fb[1])
    aspect = 1.0 / np.cos(np.radians(np.mean(fb)))
    PIX = g.pix
    extent = [g.lvals[0] + PIX / 2, g.lvals[-1] - PIX / 2, g.bvals[0] - PIX / 2, g.bvals[-1] + PIX / 2]
    rd = cfg["radial"]
    G0_LIM = tuple(figc["g0_lim"])
    fig, axs = plt.subplots(1, 3, figsize=(15, 8), sharex=True, sharey=True)
    panels = [("NH_full", rf"$N_{{\rm H}}$ [cm$^{{-2}}$] ({rd['r_min']:g}-{rd['r_max']:g} pc)", "viridis"),
              ("G0_surface", r"$G_0 = I_{\rm UV}(r_0)$ [Draine]", "magma"),
              ("alphaG_peak", r"$\alpha G = 59\,I_{\rm UV}(r_{\rm pk})/n_{\rm pk}$", "viridis")]
    ov = reg.overlay
    present = [lv for lv in (ov or {}).get("context_levels", []) if np.nanmax(ov["map"]) >= lv] if ov else []
    for ax, (k, lab, cm) in zip(axs, panels):
        d = M[k]
        f = d[np.isfinite(d) & (d > 0) & view]
        cmap = plt.get_cmap(cm).copy(); cmap.set_bad("0.85")
        if k.startswith("G0"):
            # Fixed LINEAR scale: a percentile stretch turns LightCube voxel noise
            # into full-scale structure. Actual range is annotated on the bar.
            im_ = ax.imshow(d, origin="lower", extent=extent, aspect=aspect, cmap=cmap,
                            vmin=G0_LIM[0], vmax=G0_LIM[1], interpolation="nearest")
            cb = fig.colorbar(im_, ax=ax, orientation="horizontal", pad=0.07, fraction=0.05)
            vc = d[reg.mask & np.isfinite(d)]
            cb.set_label(f"{lab}\ninside region: {vc.min():.3f}-{vc.max():.3f} "
                         f"(5-95% {np.percentile(vc,5):.3f}-{np.percentile(vc,95):.3f}); "
                         f"linear, fixed {G0_LIM[0]}-{G0_LIM[1]}", fontsize=9)
        else:
            lo, hi = np.percentile(f, [2, 98])
            im_ = ax.imshow(d, origin="lower", extent=extent, aspect=aspect, cmap=cmap,
                            norm=LogNorm(lo, hi), interpolation="nearest")
            cb = fig.colorbar(im_, ax=ax, orientation="horizontal", pad=0.07, fraction=0.05)
            cb.set_label(lab)
            ticks = [t for t in (float(f"{v:.2g}") for v in np.geomspace(lo, hi, 4)) if lo <= t <= hi]
            cb.set_ticks(ticks); cb.set_ticklabels([f"{t:.2g}" for t in ticks])
            cb.ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        if ov:
            ax.contour(g.lvals, g.bvals, ov["map"], levels=[ov["level"]], colors="r", linewidths=1.6)
            if present:
                ax.contour(g.lvals, g.bvals, ov["map"], levels=present, colors="r", linewidths=0.7,
                           alpha=0.45, linestyles="--")
        else:
            ax.contour(g.lvals, g.bvals, reg.mask.astype(float), levels=[0.5], colors="r", linewidths=1.2)
        ax.set_xlabel(r"$l$ [deg]")
    a = M["alphaG_peak"]
    hand = []
    if ov:
        hand.append(Line2D([], [], color="r", lw=1.6, label=f"{ov['label']} = {ov['level']:g}"))
        if present:
            hand.append(Line2D([], [], color="r", lw=0.7, ls="--", alpha=0.45,
                               label=f"{ov['label']} = {', '.join(f'{p:g}' for p in present)}"))
    else:
        hand.append(Line2D([], [], color="r", lw=1.2, label=f"region: {reg.description}"[:60]))
    if np.nanmin(a) < 1 < np.nanmax(a):
        axs[2].contour(g.lvals, g.bvals, a, levels=[1.0], colors="w", linewidths=1.0)
        hand.append(Line2D([], [], color="w", lw=1.0, label=r"$\alpha G = 1$"))
    for p in cfg["annotations"]:
        axs[2].plot(p["l"], p["b"], marker="+", ms=14, mew=1.6, color="k")
        axs[2].annotate(p["name"], (p["l"], p["b"]), xytext=tuple(p.get("offset", (-55, 22))),
                        textcoords="offset points", color="k", fontsize=9, fontweight="bold",
                        arrowprops=dict(arrowstyle="-", color="k", lw=0.8),
                        bbox=dict(fc="w", ec="none", alpha=0.7, pad=1))
    axs[0].set_ylabel(r"$b$ [deg]")
    axs[0].set_xlim(fl[1], fl[0])            # l increasing to the left
    axs[0].set_ylim(*fb)
    leg = axs[2].legend(handles=hand, loc="upper right", fontsize=8, framealpha=0.8)
    leg.get_frame().set_facecolor("0.6")
    missing = [lv for lv in (ov or {}).get("context_levels", []) if lv not in present] if ov else []
    note = f"; {', '.join(f'{p:g}' for p in missing)} not reached in field" if missing else ""
    r = run.res
    conv = (f"\nn_H = {run.nh_per_tau:.1f} cm$^{{-3}}$ per $\\tau_{{1565}}$/pc "
            f"({nh_per_E(run.conversion):.0f} per E/pc; '{run.conversion}'; "
            f"$A_{{1565}}/A_V$={A1565_AV:.3f}, Gordon+2023); "
            f"quality cut $\\tau_{{1565,pk}}\\geq${run.tau_peak_min:.3g}/pc; I_UV and n both from Edenhofer dust"
            f"\npixel {g.pix:g} deg = {r.pixel_pc:.2g} pc at {r.distance_pc:g} pc; dust effective "
            f"resolution {r.effective_pc:.2g} pc (oversampled x{r.oversampling:.1f})")
    fig.suptitle(f"{cfg['name']}: boundary '{run.boundary.name}'; grey = NaN{note}; "
                 f"aspect = 1/cos({np.mean(fb):.1f}°){conv}", fontsize=10)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")
    return path
