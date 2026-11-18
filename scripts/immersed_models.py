"""Non-orientable immersions from the current paper, using its M(d) bands.

The equations and face incidences are transcribed from
``sections/immersion.tex`` (eq:imm-coordinates, eq:imm-faces,
lem:uniform-core and thm:nonorientable).  This module is self-contained:
it does not import the paper's illustration or verification scripts.

Coordinates of M(d) are numerical solutions, not claimed exact algebraic
coordinates.  A positive epsilon separates the formerly singular sheets;
the inequalities in lem:band-immersion are checked for every solution.
The abstract surfaces are closed, but are non-orientable and have global
self-intersections.  A raw STL of their faces is therefore a mathematical
surface model, not the boundary of an ordinary printable solid.
"""

from functools import lru_cache

import numpy as np
from scipy.optimize import least_squares


def _rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def mobius_band(b, z, epsilon):
    """Return M(b,z,epsilon); indices 0:3, 3:6, 6:9 are A, B, C."""
    if b <= 0 or z <= 0 or epsilon < 0:
        raise ValueError("M requires b,z > 0 and epsilon >= 0")
    k = 1.0 / (b + 1.0)
    base = np.array([
        [1.0, 0.0, 0.0],
        [b, -b / np.sqrt(3.0), np.sqrt(z)],
        [b / 3.0, -(b + epsilon) / np.sqrt(3.0),
         (1.0 + k * epsilon) * np.sqrt(z)],
    ])
    vertices = np.array([
        _rotation(2.0 * np.pi * i / 3.0) @ vertex
        for vertex in base for i in range(3)
    ])
    faces = []
    for i in range(3):
        a = lambda j: j % 3
        b_index = lambda j: 3 + j % 3
        c = 6 + i
        faces.extend([
            [b_index(i), a(i + 1), b_index(i + 1)],
            [b_index(i), b_index(i + 1), c],
            [b_index(i), c, a(i)],
            [a(i - 1), b_index(i + 1), c],
            [c, a(i), a(i - 1)],
        ])
    return vertices, faces


def _angle_sums(vertices, faces):
    """Angle sums for this module's triangles and convex core polygons."""
    sums = np.zeros(len(vertices))
    for face in faces:
        for j, vertex_index in enumerate(face):
            u = vertices[face[j - 1]] - vertices[vertex_index]
            v = vertices[face[(j + 1) % len(face)]] - vertices[vertex_index]
            sums[vertex_index] += np.arctan2(np.linalg.norm(np.cross(u, v)), u @ v)
    return sums


@lru_cache(maxsize=None)
def solve_band(defect, epsilon=0.01):
    """Solve delta_B=delta_C=defect and check the paper's immersion domain.

    Epsilon is fixed rather than rounded away after solving.  The default
    0.01 is small enough for the selected representatives while retaining
    a visible separation in double-precision geometry.  Failure to attain
    the requested accuracy raises instead of silently emitting a model.
    """
    defect, epsilon = float(defect), float(epsilon)
    if abs(defect) > 2.0 * np.pi / 9.0 + 1e-12:
        raise ValueError("The proved band family has |defect| <= 2*pi/9")
    if epsilon <= 0:
        raise ValueError("A positive epsilon is required for an immersion")

    def residual(parameters):
        vertices, faces = mobius_band(*parameters, epsilon)
        defects = 2.0 * np.pi - _angle_sums(vertices, faces)
        return defects[[3, 6]] - defect

    best = None
    for guess in ((0.7, 0.4), (1.0, 1.0), (0.4, 0.15), (1.4, 2.0)):
        result = least_squares(
            residual, guess, bounds=([1e-5, 1e-6], [2.0 - 1e-8, 5.0]),
            xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=300,
        )
        error = float(np.max(np.abs(result.fun)))
        if best is None or error < best[0]:
            best = (error, result.x)
        if error < 2e-12:
            break
    error, (b, z) = best
    k = 1.0 / (b + 1.0)
    if error >= 2e-10 or not (0 < b < 2 and z > 0 and 0 < k * epsilon < 2 * b / 3):
        raise RuntimeError(f"Cannot certify numerical band parameters: {best}")
    vertices, faces = mobius_band(b, z, epsilon)
    sums = _angle_sums(vertices, faces)
    # The boundary turning is forced by Gauss--Bonnet, but check it too.
    errors = np.r_[2 * np.pi - sums[3:] - defect,
                   np.pi - sums[:3] + 2 * defect]
    if np.max(np.abs(errors)) > 2e-10:
        raise RuntimeError("Band symmetry or boundary turning check failed")
    return float(b), float(z), epsilon


def convex_core(genus):
    """Punctured S_g from lem:uniform-core, with equilateral holes of side 2."""
    if isinstance(genus, bool) or int(genus) != genus or genus < 3:
        raise ValueError("The convex core requires an integer genus >= 3")
    genus = int(genus)
    x = np.pi / (2.0 * genus)
    alpha = np.pi / 6.0 - np.pi / (3.0 * genus)
    rho = 1.0 / np.tan(x)
    q = rho - np.tan(alpha)
    h = np.sqrt(3.0 - np.tan(alpha) ** 2)
    sector = np.array([[q, 0.0, h], [rho, -1.0, 0.0], [rho, 1.0, 0.0]])
    vertices = np.vstack([
        sector @ _rotation(2.0 * np.pi * j / genus).T for j in range(genus)
    ])
    faces = [[3 * j for j in range(genus)],
             [i for j in range(genus) for i in (3 * j + 1, 3 * j + 2)]]
    holes = []
    for j in range(genus):
        following = (j + 1) % genus
        faces.append([3 * j, 3 * following, 3 * following + 1, 3 * j + 2])
        holes.append([3 * j, 3 * j + 1, 3 * j + 2])
    return vertices, faces, holes


def _attach_bands(vertices, faces, holes, band_vertices, band_faces):
    """Reuse seam indices exactly; put every band on the core's outer side."""
    center = vertices.mean(axis=0)
    result_vertices = list(vertices)
    result_faces = [list(face) for face in faces]
    for hole in holes:
        triangle = vertices[hole]
        midpoint = triangle.mean(axis=0)
        u = triangle[0] - midpoint
        w = (triangle[1] - triangle[2]) / np.sqrt(3.0)
        scale = np.linalg.norm(u)
        normal = np.cross(u, w)
        normal /= np.linalg.norm(normal)
        if normal @ (midpoint - center) < 0:
            normal = -normal
        mapped = (midpoint + band_vertices[:, 0, None] * u
                  + band_vertices[:, 1, None] * w
                  + band_vertices[:, 2, None] * (scale * normal))
        if not np.allclose(mapped[:3], triangle, rtol=0.0, atol=2e-12):
            raise RuntimeError("Band boundary does not match the core hole")
        if not np.all((mapped[3:] - midpoint) @ normal > 0):
            raise RuntimeError("Band is not on the outer side of the core")
        ids = list(hole) + list(range(len(result_vertices), len(result_vertices) + 6))
        result_vertices.extend(mapped[3:])
        result_faces.extend([[ids[i] for i in face] for face in band_faces])
    return np.array(result_vertices), result_faces


def build_immersed(genus, epsilon=0.01):
    """Build N_g: projective plane, Klein bottle, or S_g with g Mobius bands."""
    if isinstance(genus, bool) or int(genus) != genus or genus < 1:
        raise ValueError("Non-orientable genus must be a positive integer")
    genus = int(genus)
    count = 9 if genus == 1 else 15 if genus == 2 else 9 * genus
    defect = 2.0 * np.pi * (2 - genus) / count
    b, z, epsilon = solve_band(defect, epsilon)
    band_vertices, band_faces = mobius_band(b, z, epsilon)
    if genus == 1:
        vertices = band_vertices
        faces = band_faces + [[0, 1, 2]]
        construction = "Projective plane: M(2*pi/9) capped by its boundary triangle"
        label = "projective_plane"
    elif genus == 2:
        vertices = np.vstack([band_vertices, band_vertices[3:] * [1.0, 1.0, -1.0]])
        faces = band_faces + [[i if i < 3 else i + 6 for i in face] for face in band_faces]
        construction = "Flat Klein bottle: two copies of M(0) reflected across their common boundary"
        label = "klein_bottle"
    else:
        core_vertices, core_faces, holes = convex_core(genus)
        vertices, faces = _attach_bands(core_vertices, core_faces, holes, band_vertices, band_faces)
        construction = "Non-orientable surface N_g: convex core S_g with g immersed Mobius bands"
        label = "mobius_bands"
    error = float(np.max(np.abs(2.0 * np.pi - _angle_sums(vertices, faces) - defect)))
    if len(vertices) != count or error > 2e-10:
        raise RuntimeError(f"Closed immersion fails its count or defect check: {error}")
    return {
        "name": f"immersed_nonorientable_g{genus}_v{count}_{label}",
        "vertices": vertices,
        "faces": faces,
        "orientable": False,
        "genus": genus,
        "expected_vertices": count,
        "expected_defect": defect,
        "construction": construction,
        "source": "sections/immersion.tex: eq:imm-coordinates, eq:imm-faces, lem:uniform-core, thm:nonorientable",
        "parameters": {"b": b, "z": z, "epsilon": epsilon, "k": 1.0 / (b + 1.0)},
    }


def build_models():
    """Representative closed immersions; build_immersed supports any genus."""
    return [build_immersed(genus) for genus in range(1, 7)]
