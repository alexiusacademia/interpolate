"""Makes tests/cases.json: what SciPy gives for curves through thousands of
random tables, for the Tessel tests to compare against.

Needs numpy and scipy. Run from the library's folder:

    python3 tools/oracle/oracle.py

Each case is a table of points, how the curve is made (kind, ends,
outside), and SciPy's value, slope, second derivative at probe xs, areas
between pairs of xs, and the xs where the curve meets some ys. NaN is
written as 1e308 (JSON has no NaN).

Author: Alexius Academia
Date: 2026-10-25
"""

import json
import random

import numpy as np
from scipy.interpolate import (Akima1DInterpolator, CubicSpline, PchipInterpolator, PPoly,
                               RegularGridInterpolator, interp1d)

random.seed(20261025)
NAN = 1e308


def code(v):
    v = float(v)
    return NAN if np.isnan(v) else v


def table(n):
    """n points with uneven gaps, and ys of several shapes."""
    xs = np.cumsum([random.uniform(0.2, 3.0) for _ in range(n)]) + random.uniform(-10, 10)
    shape = random.choice(["random", "rising", "falling", "flat parts", "spiky"])
    if shape == "random":
        ys = [random.uniform(-5, 5) for _ in range(n)]
    elif shape == "rising":
        ys = list(np.cumsum([random.uniform(0, 4) for _ in range(n)]))
    elif shape == "falling":
        ys = list(-np.cumsum([random.uniform(0, 4) for _ in range(n)]))
    elif shape == "flat parts":
        ys = [float(random.choice([0, 0, 1, 2, 2])) for _ in range(n)]
    else:
        ys = [random.choice([0.0, 0.0, 0.0, random.uniform(5, 10)]) for _ in range(n)]
    return np.array(xs), np.array([round(y, 6) for y in ys])


def probes(xs):
    lo, hi = xs[0], xs[-1]
    span = hi - lo
    ps = [random.uniform(lo - 0.3 * span, hi + 0.3 * span) for _ in range(12)]
    ps += list(xs[:3]) + [xs[-1], (xs[0] + xs[1]) / 2]
    return np.array(ps)


def ends(n):
    def one():
        k = random.choice(["notAKnot", "natural", "slope", "bending"])
        v = round(random.uniform(-3, 3), 3)
        return k, v
    return one(), one()


def scipy_bc(e):
    k, v = e
    return {"notAKnot": "not-a-knot", "natural": "natural", "slope": (1, v), "bending": (2, v)}[k]


def piecewise(kind, xs, ys, start, end, modified):
    """The curve as a SciPy PPoly that extends its end pieces."""
    if kind == "linear":
        slopes = np.diff(ys) / np.diff(xs)
        return PPoly(np.array([slopes, ys[:-1]]), xs, extrapolate=True)
    if kind == "spline":
        s = CubicSpline(xs, ys, bc_type=(scipy_bc(start), scipy_bc(end)), extrapolate=True)
        return PPoly(s.c, s.x, extrapolate=True)
    if kind == "monotone":
        s = PchipInterpolator(xs, ys, extrapolate=True)
        return PPoly(s.c, s.x, extrapolate=True)
    a = Akima1DInterpolator(xs, ys, method="makima" if modified else "akima")
    return PPoly(a.c, a.x, extrapolate=True)


def curve_case(kind):
    n = random.randint(3 if kind == "akima" else 2, 12)
    xs, ys = table(n)
    start, end = ends(n)
    modified = random.random() < 0.5
    outside = random.choice(["extend", "extend", "hold", "nan"])
    p = piecewise(kind, xs, ys, start, end, modified)
    lo, hi = xs[0], xs[-1]
    at, slope, second = [], [], []
    ps = probes(xs)
    for x in ps:
        inside = lo <= x <= hi
        if inside or outside == "extend":
            at.append(p(x)); slope.append(p(x, 1)); second.append(p(x, 2))
        elif outside == "hold":
            at.append(ys[0] if x < lo else ys[-1]); slope.append(0.0); second.append(0.0)
        else:
            at.append(np.nan); slope.append(np.nan); second.append(np.nan)
    areas = []
    for _ in range(6):
        a, b = random.choice(ps), random.choice(ps)
        if outside == "extend":
            v = p.integrate(a, b)
        elif outside == "nan":
            v = np.nan if min(a, b) < lo or max(a, b) > hi else p.integrate(a, b)
        else:
            u, w = min(a, b), max(a, b)
            v = p.integrate(max(u, lo), min(w, hi)) if max(u, lo) < min(w, hi) else 0.0
            v += ys[0] * max(0.0, min(w, lo) - u) + ys[-1] * max(0.0, w - max(u, hi))
            if a > b:
                v = -v
        areas.append([a, b, v])
    roots = []
    for _ in range(3):
        y = random.choice([random.uniform(min(ys), max(ys)), random.choice(list(ys))])
        found = p.solve(y, discontinuity=False, extrapolate=False)
        found = [r for r in found if lo <= r <= hi]
        if any(np.isnan(found)):
            continue   # the curve is level at y: an interval, not points
        roots.append([float(y)] + sorted(float(r) for r in found))
    return {
        "kind": kind, "xs": list(xs), "ys": list(ys),
        "start": start[0], "startValue": start[1], "end": end[0], "endValue": end[1],
        "modified": modified, "outside": outside, "step": "",
        "probes": list(ps), "at": [code(v) for v in at], "slope": [code(v) for v in slope],
        "second": [code(v) for v in second],
        "areas": [[float(a), float(b), code(v)] for a, b, v in areas],
        "roots": roots,
    }


def steps_case():
    n = random.randint(2, 10)
    xs, ys = table(n)
    step = random.choice(["previous", "next", "nearest"])
    outside = random.choice(["hold", "nan"])
    f = interp1d(xs, ys, kind=step, bounds_error=False,
                 fill_value=(ys[0], ys[-1]) if outside == "hold" else np.nan)
    ps = probes(xs)
    return {
        "kind": "steps", "xs": list(xs), "ys": list(ys), "start": "", "startValue": 0.0, "end": "",
        "endValue": 0.0, "modified": False, "outside": outside, "step": step,
        "probes": list(ps), "at": [code(f(x)) for x in ps], "slope": [], "second": [], "areas": [], "roots": [],
    }


def grid_case():
    nx, ny = random.randint(2, 6), random.randint(2, 6)
    xs, _ = table(nx)
    ys, _ = table(ny)
    zs = np.round(np.random.default_rng(random.randint(0, 10**9)).uniform(-10, 10, (nx, ny)), 4)
    outside = random.choice(["extend", "hold", "nan"])
    g = RegularGridInterpolator((xs, ys), zs, method="linear", bounds_error=False,
                                fill_value=None if outside != "nan" else np.nan)
    px = probes(xs)
    py = probes(ys)
    pts = []
    for x, y in zip(px, py):
        if outside == "hold":
            v = g([[min(max(x, xs[0]), xs[-1]), min(max(y, ys[0]), ys[-1])]])[0]
        else:
            v = g([[x, y]])[0]
        pts.append([float(x), float(y), code(v)])
    return {"xs": list(xs), "ys": list(ys), "zs": zs.tolist(), "outside": outside, "points": pts}


def main():
    curves = []
    for _ in range(600):
        curves.append(curve_case(random.choice(["linear", "spline", "spline", "monotone", "akima"])))
    for _ in range(150):
        curves.append(steps_case())
    grids = [grid_case() for _ in range(150)]
    with open("tests/cases.json", "w") as f:
        json.dump({"curves": curves, "grids": grids}, f)
    checks = sum(len(c["at"]) + len(c["slope"]) + len(c["second"]) + len(c["areas"]) + len(c["roots"]) for c in curves)
    checks += sum(len(g["points"]) for g in grids)
    print(f"{len(curves)} curves, {len(grids)} grids, {checks} values -> tests/cases.json")


main()
