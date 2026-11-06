"""Pluggable cloud-boundary finders. Full contract: docs/BOUNDARIES.md.

Every method receives, per chunk of sightlines:
  D    (P, R) float  dust quantity along each ray (tau_1565 per pc for Edenhofer),
                     NaN where undefined (r < 69 pc). All geometry is decided on D,
                     never on the gas density, so it is conversion-independent.
  ctx  RayContext    l, b (P,) in deg; r (R,) in pc; X, Y, Z (P, R) heliocentric
                     Galactic Cartesian pc of every sample.
and returns a BoundaryResult with, per ray (all shape (P,)):
  ipk        index of the density peak used (r_peak = r[ipk])
  i0, i1     indices of the near (r0) and far (r1) boundary samples
  has0, has1 whether r0 / r1 exist (False -> NaN in the maps)
  valid      method-level validity (e.g. the ray actually crosses the mask);
             AND-ed with the dust-unit quality cut
  rpeak_edge peak sits on the edge of the search range
  flags      uint16 bitmask (FLAG_* below), never used to silently fix a value
"""
from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np

# ---------------------------------------------------------------------------
# Flags (bitmask) and method codes
# ---------------------------------------------------------------------------
FLAG_MULTIPEAK = 1        # >1 separate run of D >= D_peak/2 along the ray
FLAG_RPEAK_EDGE = 2       # peak on the edge of the search window / range / mask segment
FLAG_R0_MISSING = 4       # no near-side boundary found
FLAG_R1_MISSING = 8       # no far-side boundary found
FLAG_R0_RAY_START = 16    # r0 at the first radial sample (boundary may lie closer)
FLAG_R1_RAY_END = 32      # r1 at the last radial sample (boundary may lie farther)
FLAG_VALLEY_BOUND = 64    # nearest_peak: boundary placed at a valley between components
FLAG_MASK_GRAZING = 128   # mask3d: chosen in-mask chord shorter than min_chord_pc
FLAG_MASK_MULTISEG = 256  # mask3d: ray enters/exits the mask more than once
FLAG_MASK_TRUNCATED = 512 # mask3d: segment touches the sampling range or mask coverage edge
FLAG_MASK_MISS = 1024     # mask3d: ray never enters the mask

FLAG_NAMES = {v: k[5:] for k, v in globals().items() if k.startswith("FLAG_")}

METHOD_CODES = {"none": 0, "half_max": 1, "fixed_range": 2, "nearest_peak": 3, "mask3d": 4}


def describe_flags(flags):
    return [n for b, n in sorted(FLAG_NAMES.items()) if flags & b]


@dataclass
class RayContext:
    l: np.ndarray
    b: np.ndarray
    r: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    Z: np.ndarray


@dataclass
class BoundaryResult:
    ipk: np.ndarray
    i0: np.ndarray
    i1: np.ndarray
    has0: np.ndarray
    has1: np.ndarray
    valid: np.ndarray
    rpeak_edge: np.ndarray
    flags: np.ndarray


def _pick(a, i):
    return a[np.arange(a.shape[0]), i]


def _halfmax_crossings(D, ipk, dpk):
    """Largest index < ipk with D < dpk/2, smallest index > ipk with D < dpk/2."""
    NR = D.shape[1]
    idx = np.arange(NR)
    below = D < (dpk / 2)[:, None]
    near = below & (idx[None, :] < ipk[:, None])
    far = below & (idx[None, :] > ipk[:, None])
    has0, has1 = near.any(1), far.any(1)
    i0 = NR - 1 - near[:, ::-1].argmax(axis=1)   # largest index < ipk with D < half
    i1 = far.argmax(axis=1)                       # smallest index > ipk with D < half
    return i0, i1, has0, has1


def _multipeak(D, dpk):
    above = np.nan_to_num(D, nan=-np.inf) >= (dpk / 2)[:, None]
    starts = above[:, :1] | (above[:, 1:] & ~above[:, :-1])
    return np.concatenate([above[:, :1], above[:, 1:] & ~above[:, :-1]], axis=1).sum(1) > 1


def _common_flags(D, dpk, i0, i1, has0, has1, rpeak_edge):
    NR = D.shape[1]
    f = np.zeros(D.shape[0], np.uint16)
    f |= np.where(_multipeak(D, dpk), FLAG_MULTIPEAK, 0).astype(np.uint16)
    f |= np.where(rpeak_edge, FLAG_RPEAK_EDGE, 0).astype(np.uint16)
    f |= np.where(~has0, FLAG_R0_MISSING, 0).astype(np.uint16)
    f |= np.where(~has1, FLAG_R1_MISSING, 0).astype(np.uint16)
    f |= np.where(has0 & (i0 == 0), FLAG_R0_RAY_START, 0).astype(np.uint16)
    f |= np.where(has1 & (i1 == NR - 1), FLAG_R1_RAY_END, 0).astype(np.uint16)
    return f


class BoundaryMethod(ABC):
    name = "abstract"

    def __init__(self, **params):
        self.params = params

    @property
    def code(self):
        return METHOD_CODES[self.name]

    @abstractmethod
    def find(self, D, ctx: RayContext) -> BoundaryResult:
        ...

    def describe(self):
        """Lines (<= 72 chars) for FITS headers."""
        return [f"boundary method '{self.name}' {self.params}"[:72]]


# ---------------------------------------------------------------------------
class HalfMax(BoundaryMethod):
    """Peak of D within `window` [pc]; r0/r1 = near/far half-maximum crossings,
    which may lie outside the window. Eos default. Assumes one clean bump."""
    name = "half_max"

    def __init__(self, window):
        super().__init__(window=tuple(window))
        self.window = tuple(window)

    def find(self, D, ctx):
        r = ctx.r
        in_cloud = (r >= self.window[0]) & (r <= self.window[1])
        cr_first, cr_last = np.where(in_cloud)[0][[0, -1]]
        ipk = np.where(in_cloud[None, :], D, -np.inf).argmax(axis=1)
        dpk = _pick(D, ipk)
        i0, i1, has0, has1 = _halfmax_crossings(D, ipk, dpk)
        rpeak_edge = ~((ipk > cr_first) & (ipk < cr_last))
        flags = _common_flags(D, dpk, i0, i1, has0, has1, rpeak_edge)
        return BoundaryResult(ipk, i0, i1, has0, has1, np.ones(D.shape[0], bool), rpeak_edge, flags)

    def describe(self):
        w = self.window
        return [f"boundary: half_max; r_peak = argmax within {w[0]:g}-{w[1]:g} pc;",
                "r0/r1 = near/far half-max crossings (may lie outside window)."]


# ---------------------------------------------------------------------------
class FixedRange(BoundaryMethod):
    """User-supplied r0, r1 [pc]; r_peak = argmax of D between them."""
    name = "fixed_range"

    def __init__(self, r0_pc, r1_pc):
        super().__init__(r0_pc=r0_pc, r1_pc=r1_pc)
        self.r0_pc, self.r1_pc = float(r0_pc), float(r1_pc)

    def find(self, D, ctx):
        r = ctx.r
        sel = (r >= self.r0_pc) & (r <= self.r1_pc)
        if not sel.any():
            raise ValueError(f"fixed_range {self.r0_pc}-{self.r1_pc} pc has no radial samples")
        a, b = np.where(sel)[0][[0, -1]]
        P = D.shape[0]
        ipk = np.where(sel[None, :], D, -np.inf).argmax(axis=1)
        dpk = _pick(D, ipk)
        i0, i1 = np.full(P, a), np.full(P, b)
        has0 = has1 = np.ones(P, bool)
        rpeak_edge = (ipk == a) | (ipk == b)
        flags = _common_flags(D, dpk, i0, i1, has0, has1, rpeak_edge)
        return BoundaryResult(ipk, i0, i1, has0, has1, np.ones(P, bool), rpeak_edge, flags)

    def describe(self):
        return [f"boundary: fixed_range r0={self.r0_pc:g}, r1={self.r1_pc:g} pc;",
                "r_peak = argmax between them."]


# ---------------------------------------------------------------------------
class NearestPeak(BoundaryMethod):
    """Multi-component profiles. Peaks are local maxima of D smoothed over
    `smooth_pc` (kills sub-voxel interpolation wiggles) inside `window`, with
    D >= min_frac * (window max). The peak nearest `target_pc` is chosen. r0/r1
    are the half-max crossings around it, except where the neighbouring peak is
    a DISTINCT component -- the profile between them dips below (1 - min_dip) x
    the lower of the two peaks -- but not below half max: then the boundary is
    placed at that valley's minimum and FLAG_VALLEY_BOUND is set. Rays with no
    interior peak fall back to the window argmax (RPEAK_EDGE if on the edge)."""
    name = "nearest_peak"

    def __init__(self, window, target_pc, min_frac=0.5, smooth_pc=2.0, min_dip=0.2):
        super().__init__(window=tuple(window), target_pc=target_pc, min_frac=min_frac,
                         smooth_pc=smooth_pc, min_dip=min_dip)
        self.window, self.target, self.min_frac = tuple(window), float(target_pc), float(min_frac)
        self.smooth_pc, self.min_dip = float(smooth_pc), float(min_dip)

    def find(self, D, ctx):
        from scipy.ndimage import uniform_filter1d
        r = ctx.r
        P, NR = D.shape
        idx = np.arange(NR)
        win = (r >= self.window[0]) & (r <= self.window[1])
        wf, wl = np.where(win)[0][[0, -1]]
        Dz = np.nan_to_num(D, nan=-np.inf)
        nsm = max(1, int(round(self.smooth_pc / float(np.median(np.diff(r))))))
        Ds = uniform_filter1d(np.nan_to_num(D, nan=0.0), nsm, axis=1, mode="nearest")
        Ds = np.where(np.isnan(D), -np.inf, Ds)
        lm = np.zeros((P, NR), bool)
        lm[:, 1:-1] = (Ds[:, 1:-1] > Ds[:, :-2]) & (Ds[:, 1:-1] >= Ds[:, 2:])
        wmax = np.where(win[None, :], Ds, -np.inf).max(axis=1)
        sig = lm & win[None, :] & (Ds >= self.min_frac * wmax[:, None])
        has_lm = sig.any(1)
        dist = np.where(sig, np.abs(r - self.target)[None, :], np.inf)
        ipk = np.where(has_lm, dist.argmin(axis=1),
                       np.where(win[None, :], Dz, -np.inf).argmax(axis=1))
        dpk = _pick(D, ipk)
        i0, i1, has0, has1 = _halfmax_crossings(D, ipk, dpk)
        flags = np.zeros(P, np.uint16)
        # previous / next significant peak
        prev_m = sig & (idx[None, :] < ipk[:, None])
        next_m = sig & (idx[None, :] > ipk[:, None])
        has_prev, has_next = prev_m.any(1), next_m.any(1)
        prev = NR - 1 - prev_m[:, ::-1].argmax(axis=1)
        nxt = next_m.argmax(axis=1)
        between0 = (idx[None, :] > prev[:, None]) & (idx[None, :] < ipk[:, None])
        between1 = (idx[None, :] > ipk[:, None]) & (idx[None, :] < nxt[:, None])
        val0 = np.where(between0, Ds, np.inf).argmin(axis=1)
        val1 = np.where(between1, Ds, np.inf).argmin(axis=1)
        dps = _pick(Ds, ipk)
        distinct0 = _pick(Ds, val0) < (1 - self.min_dip) * np.minimum(_pick(Ds, prev), dps)
        distinct1 = _pick(Ds, val1) < (1 - self.min_dip) * np.minimum(_pick(Ds, nxt), dps)
        # valley bound only if the neighbour is a distinct component that is NOT
        # already separated from the chosen peak by a half-max crossing
        v0 = has_prev & distinct0 & (~has0 | (prev > i0))
        v1 = has_next & distinct1 & (~has1 | (nxt < i1))
        i0 = np.where(v0, val0, i0)
        i1 = np.where(v1, val1, i1)
        has0, has1 = has0 | v0, has1 | v1
        rpeak_edge = (ipk == wf) | (ipk == wl)
        flags |= _common_flags(D, dpk, i0, i1, has0, has1, rpeak_edge)
        flags |= np.where(v0 | v1, FLAG_VALLEY_BOUND, 0).astype(np.uint16)
        return BoundaryResult(ipk, i0, i1, has0, has1, np.ones(P, bool), rpeak_edge, flags)

    def describe(self):
        w = self.window
        return [f"boundary: nearest_peak to {self.target:g} pc within {w[0]:g}-{w[1]:g} pc,",
                f"peaks >= {self.min_frac:g} x window max; half-max or valley bounds."]


# ---------------------------------------------------------------------------
# mask3d: interface + segment logic. The LOADER is deliberately not implemented.
# ---------------------------------------------------------------------------
class Mask3D(ABC):
    """A 3D cloud volume. Contract (docs/BOUNDARIES.md):

    Coordinates passed in are heliocentric Galactic Cartesian, in pc
    (x toward l=0, y toward l=90, z toward b=+90), as float arrays of any shape.

    contains(x, y, z) -> bool array, True inside the cloud volume.
        Boolean masks MUST be resampled by nearest neighbour on the mask's own
        grid (never interpolated). Float masks (probability / density) are
        trilinearly interpolated on their own grid and then thresholded.
    defined(x, y, z)  -> bool array, True where the mask HAS COVERAGE. Outside
        coverage, contains() is meaningless (not "outside the cloud"); segments
        touching undefined samples are flagged MASK_TRUNCATED.
    voxel_pc          -> float, the mask's native resolution (sets the default
        minimum chord for grazing rays).
    """
    voxel_pc = None

    @abstractmethod
    def contains(self, x, y, z):
        ...

    @abstractmethod
    def defined(self, x, y, z):
        ...


def load_mask3d(spec):
    """Build a Mask3D from a config spec. NOT IMPLEMENTED until a real mask file is
    available: the format will be fixed against the supplied file, not guessed.
    See docs/BOUNDARIES.md for what the file must provide."""
    raise NotImplementedError(
        "mask3d loader is not implemented yet: no mask file has been supplied. "
        "docs/BOUNDARIES.md lists what the file must provide (grid, frame, origin, "
        "voxel size, axis order, mask semantics). Implement a Mask3D subclass for it.")


def segments_to_boundaries(inside, defined, D, r, min_chord_pc, choice="peak", target_pc=None):
    """Per-ray boundary from a sampled mask. inside/defined: (P, R) bool.

    - Segments are maximal runs of consecutive in-mask samples.
    - choice: 'peak'    segment containing the maximum of D among in-mask samples
              'nearest' segment whose midpoint is nearest target_pc
              'longest' segment with the most samples
              'first'   nearest segment to the Sun
    - r0 = r[first sample of the segment], r1 = r[last sample] (entry / exit).
    - chord = n_samples * dr; chord < min_chord_pc -> GRAZING, valid=False.
    - >1 segment -> MULTISEG (the choice above is applied, never merged).
    - segment touching r[0] / r[-1] or an undefined sample -> TRUNCATED
      (the corresponding boundary is a limit, not a measurement).
    - no in-mask sample -> MISS, valid=False.
    """
    P, NR = inside.shape
    idx = np.arange(NR)
    dr = float(np.median(np.diff(r)))
    ins = inside & defined
    starts = ins & ~np.concatenate([np.zeros((P, 1), bool), ins[:, :-1]], axis=1)
    segid = np.cumsum(starts, axis=1) * ins            # 1..k inside, 0 outside
    nseg = starts.sum(1)
    flags = np.zeros(P, np.uint16)
    miss = nseg == 0
    flags |= np.where(miss, FLAG_MASK_MISS, 0).astype(np.uint16)
    flags |= np.where(nseg > 1, FLAG_MASK_MULTISEG, 0).astype(np.uint16)

    Dz = np.nan_to_num(D, nan=-np.inf)
    kmax = int(nseg.max()) if P else 0
    if choice == "peak":
        ipk_all = np.where(ins, Dz, -np.inf).argmax(axis=1)
        seg = _pick(segid, ipk_all)
    else:
        score = np.full((P, max(kmax, 1)), np.inf)
        for k in range(1, kmax + 1):
            m = segid == k
            cnt = m.sum(1)
            if choice == "longest":
                s = -cnt.astype(float)
            elif choice == "first":
                s = np.full(P, float(k))
            elif choice == "nearest":
                if target_pc is None:
                    raise ValueError("choice='nearest' needs target_pc")
                mid = np.where(cnt > 0, (m * r[None, :]).sum(1) / np.maximum(cnt, 1), np.inf)
                s = np.abs(mid - target_pc)
            else:
                raise ValueError(f"unknown segment choice {choice!r}")
            score[:, k - 1] = np.where(cnt > 0, s, np.inf)
        seg = score.argmin(axis=1) + 1
    seg = np.where(miss, 0, seg)
    msk = (segid == seg[:, None]) & ~miss[:, None]
    i0 = msk.argmax(axis=1)
    i1 = NR - 1 - msk[:, ::-1].argmax(axis=1)
    has0 = has1 = ~miss
    ipk = np.where(msk, Dz, -np.inf).argmax(axis=1)
    dpk = _pick(D, ipk)
    chord = msk.sum(1) * dr
    graze = ~miss & (chord < min_chord_pc)
    flags |= np.where(graze, FLAG_MASK_GRAZING, 0).astype(np.uint16)
    und0 = np.where(i0 > 0, ~_pick(defined, np.maximum(i0 - 1, 0)), True)
    und1 = np.where(i1 < NR - 1, ~_pick(defined, np.minimum(i1 + 1, NR - 1)), True)
    trunc = ~miss & (und0 | und1)
    flags |= np.where(trunc, FLAG_MASK_TRUNCATED, 0).astype(np.uint16)
    rpeak_edge = ~miss & ((ipk == i0) | (ipk == i1))
    flags |= _common_flags(D, dpk, i0, i1, has0, has1, rpeak_edge)
    valid = ~miss & ~graze
    return BoundaryResult(ipk, i0, i1, has0, has1, valid, rpeak_edge, flags)


class Mask3DBoundary(BoundaryMethod):
    """r0/r1 = entry/exit of the ray through a 3D cloud mask. See docs/BOUNDARIES.md."""
    name = "mask3d"

    def __init__(self, mask, min_chord_pc=None, segment="peak", target_pc=None):
        super().__init__(mask=mask, min_chord_pc=min_chord_pc, segment=segment, target_pc=target_pc)
        self.mask = mask if isinstance(mask, Mask3D) else load_mask3d(mask)
        self.min_chord_pc = (float(min_chord_pc) if min_chord_pc is not None
                             else 2.0 * float(self.mask.voxel_pc))
        self.segment, self.target_pc = segment, target_pc

    def find(self, D, ctx):
        inside = np.asarray(self.mask.contains(ctx.X, ctx.Y, ctx.Z), bool)
        defined = np.asarray(self.mask.defined(ctx.X, ctx.Y, ctx.Z), bool)
        return segments_to_boundaries(inside, defined, D, ctx.r, self.min_chord_pc,
                                      self.segment, self.target_pc)

    def describe(self):
        return [f"boundary: mask3d entry/exit; segment choice '{self.segment}',",
                f"min chord {self.min_chord_pc:g} pc (shorter = grazing, NaN)."]


METHODS = {c.name: c for c in (HalfMax, FixedRange, NearestPeak, Mask3DBoundary)}


def build_boundary(cfg):
    cfg = dict(cfg)
    name = cfg.pop("method")
    cfg.pop("fallback", None)
    if name not in METHODS:
        raise ValueError(f"unknown boundary method {name!r}; options {sorted(METHODS)}")
    return METHODS[name](**cfg)
