"""Embedded constant-curvature polyhedra from the current manuscript.

The returned faces are the manuscript's cyclic polygon boundaries, including
concave faces.  Exporters must triangulate those polygons rather than use a
fan at their first vertex.  Coordinates have no rounding or STL dependencies.

Vertex numbering in the functions is zero based; the paper uses one based
numbering.  ``build_models()`` provides the shipped examples, while the window
families can also be constructed for any genus in their stated ranges.
"""

from itertools import product

import numpy as np


def _model(construction_id, vertices, faces, genus, construction, source,
           parameters=None, historical=False):
    vertices = np.asarray(vertices, dtype=float)
    count = len(vertices)
    result = {
        "name": f"embedded_orientable_g{genus}_v{count}_{construction_id}",
        "vertices": vertices,
        "faces": [[int(i) for i in face] for face in faces],
        "orientable": True,
        "genus": genus,
        "expected_vertices": count,
        "expected_defect": 2 * np.pi * (2 - 2 * genus) / count,
        "construction": construction,
        "source": source,
    }
    if parameters is not None:
        result["parameters"] = parameters
    if historical:
        result["historical"] = True
    return result


def _rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def build_tetrahedron():
    """The regular tetrahedron representing the sphere (defect pi)."""
    vertices = [[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]
    faces = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]
    return _model("tetrahedron", vertices, faces, 0, "Regular tetrahedron",
                  "main.tex: Embedded orientable surfaces")


def build_torus():
    """The paper's nine-vertex embedded flat torus, with its exact diagonals."""
    vertices = []
    for angle, height in [(0, 0), (np.pi / 8, 0.5), (np.pi / 4, 1)]:
        for j in range(3):
            theta = angle + 2 * np.pi * j / 3
            vertices.append([np.cos(theta), np.sin(theta), height])
    faces = []
    for p, q in [(0, 2), (2, 1), (1, 0)]:
        for j in range(3):
            k = (j + 1) % 3
            faces.extend([[3 * p + j, 3 * p + k, 3 * q + k],
                          [3 * p + j, 3 * q + j, 3 * q + k]])
    return _model("flat_torus", vertices, faces, 1,
                  "A flat torus with nine vertices",
                  "main.tex: A flat torus with nine vertices; "
                  "sections/torus-embedding.tex",
                  {"theta": float(np.pi / 4), "ring_heights": [0, 0.5, 1]})


def build_windows(g, shared=False):
    """Theorem 'A prism with windows': 8g vertices, or 6g if shared.

    ``shared=False`` requires g >= 2.  The shared family requires g >= 5;
    at g=4 its inner vertices would coincide and no embedding results.
    """
    if not isinstance(g, int) or isinstance(g, bool) or g < (5 if shared else 2):
        raise ValueError("The window family needs integer g >= "
                         + ("5" if shared else "2"))
    a = np.pi * (g - 1) / ((3 if shared else 4) * g)
    prism_sides = g if shared else 2 * g
    radius = 1 / np.tan(np.pi / prism_sides)
    s, t = np.sin(a) / 2, np.tan(a)
    local = np.array([
        [radius, -1, 1], [radius, -1, -1],
        [radius - t, 0, 1], [radius - t, 0, -1],
        [radius - s * t, -1 + s, 0], [radius - s * t, 1 - s, 0],
        [radius, 1, 1], [radius, 1, -1],
    ])
    per_sector = 6 if shared else 8
    vertices = np.vstack([local[:per_sector] @ _rotation(2 * np.pi * k / g).T
                          for k in range(g)])

    def ix(i, k=0):
        if shared and i in [7, 8]:
            return ((k + 1) % g) * per_sector + i - 7
        return (k % g) * per_sector + i - 1

    faces = []
    for level in ([1, 3], [2, 4]) if shared else ([1, 3, 7], [2, 4, 8]):
        faces.append([ix(i, k) for k in range(g) for i in level])
    for k in range(g):
        for ids in [[1, 3, 4, 2, 5], [7, 3, 4, 8, 6],
                    [1, 7, 6, 5], [2, 8, 6, 5], [1, 2, 8, 7]]:
            faces.append([ix(i, k) for i in ids])
        if not shared:
            faces.append([ix(7, k), ix(8, k), ix(2, k + 1), ix(1, k + 1)])
    return _model(
        "shared_windows" if shared else "windows", vertices, faces, g,
        "A prism with windows (6g vertices, shared vertices)" if shared else
        "A prism with windows (8g vertices)",
        "main.tex: Theorem thm:windows (ii)" if shared else
        "main.tex: Theorem thm:windows (i)",
        {"a": float(a), "R": float(radius), "s": float(s), "t": float(t),
         "prism_sides": prism_sides, "shared": shared},
    )


def build_central_handle(g):
    """The central-handle family with 7g-7 vertices for g=4, 5, 6."""
    if not isinstance(g, int) or isinstance(g, bool) or g not in (4, 5, 6):
        raise ValueError("The central-handle construction needs g in {4, 5, 6}")
    n = g - 1
    a = np.pi / n
    rho = 1 / np.tan(a)
    tau = 1 / np.tan(a + np.pi / 14)
    x0 = np.sin(np.pi / 7) / (3 * np.sqrt(np.sin(a)**2 - np.sin(np.pi / 7)**2))
    y0 = 1 - np.sin(a + np.pi / 14) * np.sin(2 * np.pi / 7) / (
        3 * (np.sin(a + np.pi / 14) + np.cos(2 * np.pi / 7)))
    local = np.array([
        [rho, -1, 1 / 3], [rho, -1, -1 / 3],
        [rho - tau, 0, 1 / 3], [rho - tau, 0, -1 / 3],
        [rho - (1 - y0) * tau, -y0, 0], [rho - (1 - y0) * tau, y0, 0],
        [rho - tau - x0, 0, 0],
    ])
    vertices = np.vstack([local @ _rotation(2 * np.pi * k / n).T
                          for k in range(n)])

    def ix(i, k=0):
        return (k % n) * 7 + i - 1

    faces = []
    for k in range(n):
        specs = [
            [(1, k), (2, k), (2, k + 1), (1, k + 1)],
            [(1, k), (3, k), (3, k - 1)],
            [(2, k), (4, k), (4, k - 1)],
            [(1, k), (5, k), (6, k), (1, k + 1)],
            [(2, k), (5, k), (6, k), (2, k + 1)],
            [(1, k), (3, k), (4, k), (2, k), (5, k)],
            [(1, k + 1), (3, k), (4, k), (2, k + 1), (6, k)],
            [(3, k), (7, k), (7, k + 1), (3, k + 1)],
            [(4, k), (7, k), (7, k + 1), (4, k + 1)],
        ]
        faces.extend([[ix(i, j) for i, j in face] for face in specs])
    return _model("central_handle", vertices, faces, g,
                  "A central handle (7g-7 vertices)",
                  "sections/central-handle.tex: A central handle",
                  {"a": float(a), "rho": float(rho), "tau": float(tau),
                   "x0": float(x0), "y0": float(y0)})


def build_drilled_windows(radius=0.2):
    """One central octagonal tunnel: genus 3, 32 vertices, defect -pi/4.

    The eight-gon fits in the central region of the g=2 window model.
    The punctured caps are subdivided into six simple polygons, without
    adding vertices; their radial boundaries are paired exactly.
    """
    clearance = 1 - np.tan(np.pi / 8)
    if not np.isfinite(radius) or not 0 < radius < clearance / 2:
        raise ValueError(f"Choose a positive tunnel radius below {clearance / 2}")
    base = build_windows(2)
    vertices = list(base["vertices"])
    faces = [face[:] for face in base["faces"][2:]]
    for level, height in enumerate([1, -1]):
        hole = []
        for j in range(8):
            theta = np.pi * j / 4
            hole.append(len(vertices))
            vertices.append([radius * np.cos(theta), radius * np.sin(theta), height])
        outer = base["faces"][level]
        for pos, a in enumerate(outer):
            b = outer[(pos + 1) % len(outer)]
            ja = int(round(np.arctan2(vertices[a][1], vertices[a][0]) / (np.pi / 4))) % 8
            jb = int(round(np.arctan2(vertices[b][1], vertices[b][0]) / (np.pi / 4))) % 8
            count = (jb - ja) % 8
            faces.append([a, b] + [hole[(ja + step) % 8]
                                   for step in range(count, -1, -1)])
    for j in range(8):
        k = (j + 1) % 8
        faces.append([16 + j, 16 + k, 24 + k, 24 + j])
    return _model("drilled_windows", vertices, faces, 3,
                  "Drilling an octagonal tunnel through the genus-two window surface",
                  "sections/embedded-details.tex: Drilling",
                  {"base_genus": 2, "tunnels": 1, "tunnel_sides": 8,
                   "tunnel_radius": float(radius)})


def build_historical_cube(b=0.32, c=0.085):
    """Corrected 24-vertex cube model retained from the earlier manuscript.

    The former STL's coordinates did not give equal defects.  Here the
    C vertices have height 1-4*sqrt(3)*c with parameters in the valid range.
    The present paper uses the window construction instead; this model is
    explicitly marked historical in the catalog.
    """
    if not (0 < c < 1 / (4 * np.sqrt(3)) and c < b < 1 / np.sqrt(3)):
        raise ValueError("Require 0 < c < min(b,1/(4*sqrt(3))) and b < 1/sqrt(3)")
    base = np.array([[1, 1, 1], [1 - np.sqrt(3) * b, 1 - b, 1],
                     [1 - np.sqrt(3) * c, 1 - c, 1 - 4 * np.sqrt(3) * c]])
    signs = list(product([1, -1], repeat=3))
    keys = [(k, *s) for k in range(3) for s in signs]
    lookup = {key: i for i, key in enumerate(keys)}
    vertices = np.array([base[k] * s for k, *s in keys])
    representatives = [
        [(0, 1, 1, 1), (0, 1, -1, 1), (0, 1, -1, -1), (0, 1, 1, -1)],
        [(0, 1, 1, 1), (0, -1, 1, 1), (0, -1, 1, -1), (0, 1, 1, -1)],
        [(0, 1, 1, 1), (0, 1, -1, 1), (1, 1, -1, 1), (1, -1, -1, 1),
         (0, -1, -1, 1), (0, -1, 1, 1), (1, -1, 1, 1), (1, 1, 1, 1)],
        [(0, 1, 1, 1), (1, 1, 1, 1), (1, 1, 1, -1), (0, 1, 1, -1),
         (2, 1, 1, -1), (2, 1, 1, 1)],
        [(1, 1, 1, 1), (1, -1, 1, 1), (1, -1, 1, -1), (1, 1, 1, -1)],
        [(2, 1, 1, 1), (2, -1, 1, 1), (2, -1, 1, -1), (2, 1, 1, -1)],
        [(0, 1, 1, 1), (0, -1, 1, 1), (2, -1, 1, 1), (2, 1, 1, 1)],
    ]
    faces, seen = [], set()
    for representative in representatives:
        for sx, sy, sz in signs:
            face = [lookup[k, x * sx, y * sy, z * sz]
                    for k, x, y, z in representative]
            key = frozenset(face)
            if key not in seen:
                seen.add(key)
                faces.append(face)
    return _model("historical_cube", vertices, faces, 2,
                  "Cube construction (earlier manuscript; corrected coordinates)",
                  "Earlier manuscript cube construction; scripts/fig_embedded.py: cube_mesh",
                  {"b": float(b), "c": float(c)}, historical=True)


def build_models():
    """Return the embedded examples included in the generated model catalog."""
    return ([build_tetrahedron(), build_torus()]
            + [build_windows(g) for g in range(2, 9)]
            + [build_windows(g, shared=True) for g in range(5, 9)]
            + [build_central_handle(g) for g in (4, 5, 6)]
            + [build_drilled_windows(), build_historical_cube()])
