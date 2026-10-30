# interpolate: writes tests/scattered.json, the cases for scattered points,
# with the answers from SciPy, NumPy and Matplotlib. Run from the project
# folder:
#
#   .venv/bin/python tools/oracle/scattered.py
#
# It needs numpy, scipy and matplotlib. The cases:
#   triangulations  random points (even, clustered, stretched, far from the
#                   origin like UTM, jittered grids) and SciPy's Delaunay
#                   triangles and convex hull
#   values          griddata "linear" and "nearest", and inverse distance
#                   weighting worked out with NumPy, at points inside and
#                   outside the outline
#   rbf             RBFInterpolator with every kernel, smoothing and degree
#   contours        Matplotlib's contour lines on the same triangles: how
#                   many lines at each level, how many closed, and their
#                   total length
#
# Author: Alexius Academia
# Date: 2026-10-30

import json
import math
import os
import random

import numpy as np
from scipy.interpolate import griddata, RBFInterpolator
from scipy.spatial import Delaunay, ConvexHull
import matplotlib
matplotlib.use("Agg")
from matplotlib.tri import Triangulation as MplTriangulation  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SEED = 20261030
NAN = 1e308


def enc(v):
    return [NAN if (isinstance(x, float) and math.isnan(x)) else float(x) for x in v]


def points(rng, n, kind):
    if kind == "even":
        xs = [rng.uniform(0, 100) for _ in range(n)]
        ys = [rng.uniform(0, 100) for _ in range(n)]
    elif kind == "clustered":
        xs, ys = [], []
        centres = [(rng.uniform(0, 100), rng.uniform(0, 100)) for _ in range(rng.randint(2, 5))]
        for _ in range(n):
            cx, cy = rng.choice(centres)
            xs.append(cx + rng.gauss(0, 4))
            ys.append(cy + rng.gauss(0, 4))
    elif kind == "stretched":
        xs = [rng.uniform(0, 1000) for _ in range(n)]
        ys = [rng.uniform(0, 1) for _ in range(n)]
    elif kind == "utm":
        xs = [500000 + rng.uniform(0, 2000) for _ in range(n)]
        ys = [1600000 + rng.uniform(0, 2000) for _ in range(n)]
    else:  # jittered grid
        side = max(2, int(math.sqrt(n)))
        xs, ys = [], []
        for i in range(side):
            for j in range(side):
                xs.append(i * 10 + rng.uniform(-2, 2))
                ys.append(j * 10 + rng.uniform(-2, 2))
    return [round(x, 6) for x in xs], [round(y, 6) for y in ys]


def values_at(rng, xs, ys):
    """A surface with hills, for the values at the points."""
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    sx = (x1 - x0) or 1.0
    sy = (y1 - y0) or 1.0
    a, b, c = rng.uniform(-5, 5), rng.uniform(-5, 5), rng.uniform(1, 4)
    return [round(100 + a * (x - x0) / sx + b * (y - y0) / sy +
                  c * math.sin(3 * (x - x0) / sx) * math.cos(2 * (y - y0) / sy) + rng.uniform(-0.2, 0.2), 6)
            for x, y in zip(xs, ys)]


def queries(rng, xs, ys, k):
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    dx = (x1 - x0) or 1.0
    dy = (y1 - y0) or 1.0
    qx = [rng.uniform(x0 - 0.1 * dx, x1 + 0.1 * dx) for _ in range(k)]
    qy = [rng.uniform(y0 - 0.1 * dy, y1 + 0.1 * dy) for _ in range(k)]
    # Some at the points themselves.
    for _ in range(3):
        i = rng.randrange(len(xs))
        qx.append(xs[i])
        qy.append(ys[i])
    return qx, qy


def idw(xs, ys, zs, qx, qy, power, radius, neighbours):
    P = np.array([xs, ys]).T
    Z = np.array(zs)
    out = []
    for x, y in zip(qx, qy):
        d = np.hypot(P[:, 0] - x, P[:, 1] - y)
        idx = np.argsort(d, kind="stable")
        if radius > 0:
            idx = idx[d[idx] <= radius]
        if neighbours > 0:
            idx = idx[:neighbours]
        if len(idx) == 0:
            out.append(float("nan"))
            continue
        if d[idx[0]] == 0:
            out.append(float(Z[idx[0]]))
            continue
        w = 1.0 / d[idx] ** power
        out.append(float(np.sum(w * Z[idx]) / np.sum(w)))
    return out


def contour_stats(xs, ys, zs, triangles, levels):
    tri = MplTriangulation(xs, ys, triangles)
    from matplotlib import _tri
    gen = _tri.TriContourGenerator(tri.get_cpp_triangulation(), np.array(zs, dtype=float))
    out = []
    for level in levels:
        segs, kinds = gen.create_contour(level)
        lines = 0
        closed = 0
        length = 0.0
        # One array of points for each line; kinds: 1 starts a line, 2
        # goes on, 79 closes it (back to its start).
        for pts, kd in zip(segs, kinds):
            prev = None
            for p, k in zip(pts, kd):
                if k == 1:
                    lines += 1
                elif k == 79:
                    closed += 1
                if prev is not None:
                    length += math.hypot(p[0] - prev[0], p[1] - prev[1])
                prev = p
        out.append({"level": level, "lines": lines, "closed": closed, "length": length})
    return out


KERNELS = {
    "thinPlate": "thin_plate_spline", "linear": "linear", "cubic": "cubic", "quintic": "quintic",
    "multiquadric": "multiquadric", "inverseMultiquadric": "inverse_multiquadric",
    "inverseQuadratic": "inverse_quadratic", "gaussian": "gaussian",
}
MIN_DEGREE = {"multiquadric": 0, "linear": 0, "thinPlate": 1, "cubic": 1, "quintic": 2}


def phi(kernel, r):
    if kernel == "thinPlate":
        return np.where(r == 0, 0.0, r * r * np.log(np.where(r == 0, 1.0, r)))
    return {"linear": lambda r: -r, "cubic": lambda r: r ** 3, "quintic": lambda r: -r ** 5,
            "multiquadric": lambda r: -np.sqrt(1 + r * r), "inverseMultiquadric": lambda r: 1 / np.sqrt(1 + r * r),
            "inverseQuadratic": lambda r: 1 / (1 + r * r), "gaussian": lambda r: np.exp(-r * r)}[kernel](r)


def second_opinion(px, py, zs, kernel, epsilon, smoothing, deg, qx, qy):
    """The surface's values at the queries from the same system as SciPy's,
    solved again with NumPy: where they and SciPy's disagree, rounding
    decides the answer, and the case can't test anything."""
    P = np.array([px, py]).T
    n = len(px)
    d = np.hypot(P[:, None, 0] - P[None, :, 0], P[:, None, 1] - P[None, :, 1]) * epsilon
    K = phi(kernel, d) + smoothing * np.eye(n)
    lo, hi = P.min(0), P.max(0)
    shift = (lo + hi) / 2
    scale = np.where(hi > lo, (hi - lo) / 2, 1.0)

    def poly(Q):
        u = (Q - shift) / scale
        cols = []
        for total in range(deg + 1):
            for j in range(total + 1):
                cols.append(u[:, 0] ** (total - j) * u[:, 1] ** j)
        return np.array(cols).T if cols else np.zeros((len(Q), 0))

    M = poly(P)
    m = M.shape[1]
    A = np.block([[K, M], [M.T, np.zeros((m, m))]])
    b = np.concatenate([np.array(zs, dtype=float), np.zeros(m)])
    sol = np.linalg.solve(A, b)
    Q = np.array([qx, qy]).T
    dq = np.hypot(Q[:, None, 0] - P[None, :, 0], Q[:, None, 1] - P[None, :, 1]) * epsilon
    return phi(kernel, dq) @ sol[:n] + poly(Q) @ sol[n:]


def main():
    rng = random.Random(SEED)
    tris = []
    kinds = ["even", "clustered", "stretched", "utm", "grid"]
    for k in range(250):
        kind = kinds[k % len(kinds)]
        n = rng.choice([3, 4, 5, 8, 12, 20, 40, 80, 150, 300, 600])
        xs, ys = points(rng, n, kind)
        # No repeated points in these cases.
        seen = set()
        px, py = [], []
        for x, y in zip(xs, ys):
            if (x, y) not in seen:
                seen.add((x, y))
                px.append(x)
                py.append(y)
        if len(px) < 3:
            continue
        P = np.array([px, py]).T
        try:
            d = Delaunay(P)
            hull = ConvexHull(P)
        except Exception:
            continue
        zs = values_at(rng, px, py)
        qx, qy = queries(rng, px, py, 40)
        Q = np.array([qx, qy]).T
        lin = griddata(P, np.array(zs), Q, method="linear").tolist()
        near = griddata(P, np.array(zs), Q, method="nearest").tolist()
        power = rng.choice([1.0, 2.0, 2.0, 3.0, 1.5])
        span = max(max(px) - min(px), max(py) - min(py))
        radius = rng.choice([0.0, 0.0, 0.3 * span])
        neighbours = rng.choice([0, 0, 5, 12])
        idv = idw(px, py, zs, qx, qy, power, radius, neighbours)
        triangles = sorted(sorted(int(v) for v in s) for s in d.simplices)
        lo, hi = min(zs), max(zs)
        levels = [round(lo + (hi - lo) * f, 6) for f in (0.13, 0.37, 0.5, 0.71, 0.93)]
        cont = contour_stats(px, py, zs, np.array(d.simplices), levels)
        tris.append({"kind": kind, "xs": px, "ys": py, "zs": zs, "triangles": triangles,
                     "hull": sorted(int(v) for v in hull.vertices), "qx": qx, "qy": qy,
                     "linear": enc(lin), "nearest": enc(near), "power": power, "radius": radius,
                     "neighbours": neighbours, "idw": enc(idv), "contours": cont})
    rbfs = []
    skipped = 0
    names = list(KERNELS)
    for k in range(320):
        kernel = names[k % len(names)]
        kind = kinds[(k // len(names)) % len(kinds)]
        n = rng.choice([6, 10, 20, 40, 80, 120])
        xs, ys = points(rng, n, kind)
        seen = set()
        px, py = [], []
        for x, y in zip(xs, ys):
            if (x, y) not in seen:
                seen.add((x, y))
                px.append(x)
                py.append(y)
        zs = values_at(rng, px, py)
        span = max(max(px) - min(px), max(py) - min(py)) or 1.0
        # epsilon about the reciprocal of the points' spacing, as these
        # kernels are used.
        spacing = span / math.sqrt(len(px))
        epsilon = 1.0 if kernel in ("thinPlate", "linear", "cubic", "quintic") else round(rng.uniform(0.3, 2.0) / spacing, 9)
        smoothing = rng.choice([0.0, 0.0, 0.0, 0.01, 1.0])
        least = MIN_DEGREE.get(kernel, -1)
        degree = rng.choice([-2, -2, max(least, 0) + 1, -1 if least < 0 else -2])
        deg = max(least, 0) if degree == -2 else degree
        qx, qy = queries(rng, px, py, 25)
        try:
            f = RBFInterpolator(np.array([px, py]).T, np.array(zs), kernel=KERNELS[kernel], epsilon=epsilon,
                                smoothing=smoothing, degree=deg)
            v = f(np.array([qx, qy]).T).tolist()
        except Exception:
            continue
        if not all(math.isfinite(x) for x in v):
            continue
        # Leave out systems so badly conditioned that rounding alone makes
        # two solvers disagree.
        try:
            with np.errstate(all="ignore"):
                again = second_opinion(px, py, zs, kernel, epsilon, smoothing, deg, qx, qy)
        except Exception:
            skipped += 1
            continue
        size = max(abs(z) for z in zs)
        if not np.all(np.isfinite(again)) or max(abs(a - b) for a, b in zip(again, v)) > 1e-8 * size:
            skipped += 1
            continue
        rbfs.append({"kind": kind, "xs": px, "ys": py, "zs": zs, "kernel": kernel, "epsilon": epsilon,
                     "smoothing": smoothing, "degree": degree, "qx": qx, "qy": qy, "values": v})
    path = os.path.join(ROOT, "tests", "scattered.json")
    with open(path, "w") as f:
        json.dump({"seed": SEED, "triangulations": tris, "rbf": rbfs}, f, separators=(",", ":"))
    values = sum(len(t["qx"]) * 3 for t in tris) + sum(len(r["qx"]) for r in rbfs)
    print("tests/scattered.json: {} triangulations ({} triangles), {} rbf surfaces ({} badly conditioned left out), "
          "{} values, {} contour levels"
          .format(len(tris), sum(len(t["triangles"]) for t in tris), len(rbfs), skipped, values,
                  sum(len(t["contours"]) for t in tris)))


if __name__ == "__main__":
    main()
