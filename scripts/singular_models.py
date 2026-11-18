"""Singular polygonal surfaces from Appendix B and the supplement.

Coordinates use the manuscript's scale. These are closed *abstract* surfaces;
their maps into R^3 may have self-intersections and singular vertex stars.
In particular, a non-orientable surface is not the boundary of a solid.
The exporter must retain the polygon incidences before triangulating faces.

No files from the manuscript directory are needed to regenerate these data.
"""

from itertools import combinations, permutations, product

import numpy as np


PI = np.pi
SQRT3 = np.sqrt(3.0)
R_FACES = [[0, 1, 2], [0, 4, 5], [1, 3, 5], [2, 3, 4],
           [1, 2, 4, 5], [2, 0, 5, 3], [0, 1, 3, 4]]


def rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def model(name, vertices, faces, genus, construction, source,
          orientable=False, parameters=None):
    vertices = np.asarray(vertices, dtype=float)
    chi = 2 - (2 if orientable else 1) * genus
    return dict(name=name, vertices=vertices,
                faces=[list(map(int, face)) for face in faces],
                orientable=orientable, genus=genus,
                expected_vertices=len(vertices),
                expected_defect=2 * PI * chi / len(vertices),
                construction=construction, source=source,
                parameters=parameters or {})


def projective_block(r, h):
    """R(r,h), including its triangular cap, with the paper's face cycles."""
    vertices = np.array([[1, 0, 0], [-.5, SQRT3 / 2, 0],
                         [-.5, -SQRT3 / 2, 0], [-r, 0, h],
                         [r / 2, -SQRT3 * r / 2, h],
                         [r / 2, SQRT3 * r / 2, h]])
    return vertices, [face[:] for face in R_FACES]


def tetrahemihexahedron():
    vertices, faces = projective_block(1, np.sqrt(2))
    return model("singular_nonorientable_g1_v6_tetrahemihexahedron",
                 vertices, faces, 1, "Tetrahemihexahedron R(1,sqrt(2))",
                 "Appendix B.1, The projective plane and Klein bottle",
                 parameters={"r": 1.0, "h": float(np.sqrt(2))})


def klein_bottle():
    # The caption's height is used to match the paper's spatial illustration.
    h = .5 * np.sqrt(3 * (1 + SQRT3))
    vertices, faces = projective_block(.5, h)
    vertices = np.vstack([vertices, vertices[3:] * [1, 1, -1]])
    faces = faces[1:] + [[i if i < 3 else i + 3 for i in face]
                        for face in faces[1:]]
    return model("singular_nonorientable_g2_v9_klein_bottle",
                 vertices, faces, 2, "Doubled flat band R^circ(1/2,h)",
                 "Appendix B.1, The projective plane and Klein bottle",
                 parameters={"r": .5, "h": float(h)})


def singular_torus():
    vertices = np.array([[0, SQRT3 / 2, 0], [-.5, 0, 0], [.5, 0, 0],
                         [0, SQRT3 / 2, 1], [-.5, 0, 1], [.5, 0, 1]])
    faces = [[int(i) - 1 for i in word]
             for word in ("1254", "2365", "3146", "345", "156", "246",
                          "126", "234", "135")]
    return model("singular_orientable_g1_v6_flat_torus", vertices, faces, 1,
                 "Six-vertex singular flat torus",
                 "Appendix B.2, A flat torus with six vertices",
                 orientable=True)


def uniform_core(g):
    """Punctured S_g, scaled by sqrt(3)/2 to fit R's boundary."""
    if g < 3:
        raise ValueError("The uniform core requires genus >= 3")
    x, alpha = PI / (2 * g), PI / 6 - PI / (3 * g)
    rho = 1 / np.tan(x)
    q, h = rho - np.tan(alpha), np.sqrt(3 - np.tan(alpha) ** 2)
    local = np.array([[q, 0, h], [rho, -1, 0], [rho, 1, 0]])
    vertices = np.vstack([local @ rotation(2 * PI * j / g).T
                          for j in range(g)]) * SQRT3 / 2
    faces = [[3 * j for j in range(g)],
             [i for j in range(g) for i in (3 * j + 1, 3 * j + 2)]]
    cuts = []
    for j in range(g):
        k = (j + 1) % g
        faces.append([3 * j, 3 * k, 3 * k + 1, 3 * j + 2])
        cuts.append([3 * j, 3 * j + 1, 3 * j + 2])
    return vertices, faces, cuts


def attach_bands(vertices, faces, cuts, r, h):
    """Attach outward copies of R^circ(r,h) along equilateral triangles."""
    vertices = np.asarray(vertices)
    result_vertices = list(vertices)
    result_faces = [face[:] for face in faces]
    block, block_faces = projective_block(r, h)
    center = vertices.mean(axis=0)
    for tri in cuts:
        points = vertices[tri]
        midpoint = points.mean(axis=0)
        u = points[0] - midpoint
        w = (points[1] - points[2]) / SQRT3
        normal = np.cross(u, w)
        normal /= np.linalg.norm(normal)
        if normal @ (midpoint - center) < 0:
            normal = -normal
        if not (np.allclose(np.linalg.norm(u), 1)
                and np.allclose(np.linalg.norm(w), 1)):
            raise ValueError("Attachment triangles must have side sqrt(3)")
        new = [midpoint + p[0] * u + p[1] * w + p[2] * normal
               for p in block[3:]]
        indices = tri + list(range(len(result_vertices), len(result_vertices) + 3))
        result_vertices.extend(new)
        result_faces.extend([[indices[i] for i in face] for face in block_faces[1:]])
    return np.array(result_vertices), result_faces


def core_and_bands(g):
    vertices, faces, cuts = uniform_core(g)
    s = np.sin(PI / (3 * g))
    r = 2 * s / (1 + 2 * s)
    h = np.sqrt(2 * (1 - 2 * s) * (1 + s)) / (1 + 2 * s)
    vertices, faces = attach_bands(vertices, faces, cuts, r, h)
    return model(f"singular_nonorientable_g{g}_v{6*g}_uniform_core",
                 vertices, faces, g, "Uniform core S_g and g projective-plane bands",
                 "Appendix B.3, Higher genus and vertex counts (6g construction)",
                 parameters={"g": g, "r": float(r), "h": float(h)})


def antiprism_and_bands(g):
    if g not in (3, 5, 7, 9, 11):
        raise ValueError("The paper's antiprism construction uses g=3,5,7,9,11")
    a, x = 2 * PI * (g - 2) / (5 * g), PI / g
    rho = SQRT3 / (2 * np.tan(x))
    q = SQRT3 * np.sin(a / 4) / np.sin(x)
    h2 = np.sqrt(9 / 4 - (rho - q) ** 2)
    local = np.array([[q, 0, h2], [rho, -SQRT3 / 2, 0]])
    vertices = np.vstack([local @ rotation(2 * PI * j / g).T
                          for j in range(g)])
    faces = [[2 * j for j in range(g)], [2 * j + 1 for j in range(g)]]
    cuts = []
    for j in range(g):
        k = (j + 1) % g
        faces.append([2 * j, 2 * k, 2 * k + 1])
        cuts.append([2 * j, 2 * j + 1, 2 * k + 1])
    r = .5 * (1 - np.sqrt((11 + 14 * np.cos(a / 2)) /
                         (1 + 2 * np.cos(a / 2))) * np.tan(a / 4))
    vertices, faces = attach_bands(vertices, faces, cuts, r, 1)
    return model(f"singular_nonorientable_g{g}_v{5*g}_antiprism",
                 vertices, faces, g, "Antiprism core and g projective-plane bands",
                 "Appendix B.3 and Supplement, Odd-genus singular models",
                 parameters={"g": g, "r": float(r), "h": 1.0,
                             "a": float(a), "core_height": float(h2)})


def face_cycle(vertices, indices, normal):
    """Cyclic order of a convex planar face, including central polygons."""
    indices = np.asarray(indices, dtype=int)
    points = vertices[indices]
    center = points.mean(axis=0)
    normal = np.asarray(normal, dtype=float)
    normal /= np.linalg.norm(normal)
    u = points[0] - center
    u /= np.linalg.norm(u)
    w = np.cross(normal, u)
    angles = np.arctan2((points - center) @ w, (points - center) @ u)
    return list(map(int, indices[np.argsort(angles)]))


def cubohemioctahedron():
    vertices = np.array(sorted(set(permutations((0, 1, 1))) |
                               set(permutations((0, -1, 1))) |
                               set(permutations((0, -1, -1)))), dtype=float)
    faces = []
    for dim in range(3):
        for sign in (-1, 1):
            indices = np.flatnonzero(vertices[:, dim] == sign)
            faces.append(face_cycle(vertices, indices, np.eye(3)[dim] * sign))
    for normal in np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1], [-1, 1, 1]]):
        indices = np.flatnonzero(np.abs(vertices @ normal) < 1e-10)
        faces.append(face_cycle(vertices, indices, normal))
    return model("singular_nonorientable_g4_v12_cubohemioctahedron",
                 vertices, faces, 4, "Cubohemioctahedron",
                 "Appendix B.3, Higher genus and vertex counts")


def historical_cuboctahedron_and_bands():
    """The earlier 24-vertex genus-four assembly, with corrected exact seams.

    This uses the regular cuboctahedron, not the cubohemioctahedron and not
    the current paper's S_4. All three names must remain distinct. The four
    attachment triangles correspond to normals whose sign product is +1;
    their vertex sets partition the twelve vertices of the cuboctahedron.
    """
    vertices = cubohemioctahedron()["vertices"] * np.sqrt(3 / 2)
    faces = []
    for dim in range(3):
        for sign in (-1, 1):
            normal = np.eye(3)[dim] * sign
            ids = np.flatnonzero(np.isclose(vertices @ normal, np.sqrt(3 / 2)))
            faces.append(face_cycle(vertices, ids, normal))
    cuts = []
    for signs in product((-1, 1), repeat=3):
        normal = np.array(signs, dtype=float)
        ids = np.flatnonzero(np.isclose(vertices @ normal, 2 * np.sqrt(3 / 2)))
        triangle = face_cycle(vertices, ids, normal)
        (cuts if np.prod(signs) == 1 else faces).append(triangle)
    if sorted(i for triangle in cuts for i in triangle) != list(range(12)):
        raise RuntimeError("Cuboctahedron attachment triangles must partition vertices")
    s = np.sin(PI / 12)
    r = 2 * s / (1 + 2 * s)
    h = np.sqrt(2 * (1 - 2 * s) * (1 + s)) / (1 + 2 * s)
    vertices, faces = attach_bands(vertices, faces, cuts, r, h)
    result = model("singular_nonorientable_g4_v24_historical_cuboctahedron",
                   vertices, faces, 4,
                   "Historical regular cuboctahedron with four projective-plane bands",
                   "Earlier manuscript genus-four assembly; band R from Appendix B.1",
                   parameters={"r": float(r), "h": float(h),
                               "core": "regular cuboctahedron"})
    result["historical"] = True
    return result


def small_rhombihexahedron():
    s = 1 + np.sqrt(2)
    vertices = np.array(sorted({p for a, b, c in product((-1, 1), repeat=3)
                                for p in permutations((a, b, c * s))}))
    faces = []
    for zero in range(3):
        for a, b in product((-1, 1), repeat=2):
            normal = np.roll(np.array([a, b, 0.]), zero)
            indices = np.flatnonzero(np.isclose(vertices @ normal, s + 1))
            faces.append(face_cycle(vertices, indices, normal))
    for dim in range(3):
        for sign in (-1, 1):
            indices = np.flatnonzero(np.isclose(vertices[:, dim], sign))
            faces.append(face_cycle(vertices, indices, np.eye(3)[dim]))
    return model("singular_nonorientable_g8_v24_small_rhombihexahedron",
                 vertices, faces, 8, "Small rhombihexahedron",
                 "Appendix B.3, Higher genus and vertex counts")


def small_dodecahemidodecahedron():
    phi = (1 + np.sqrt(5)) / 2
    ico = np.array([(0, a, b * phi) for a, b in product((-1, 1), repeat=2)] +
                   [(a, b * phi, 0) for a, b in product((-1, 1), repeat=2)] +
                   [(b * phi, 0, a) for a, b in product((-1, 1), repeat=2)])
    edges = [(i, j) for i, j in combinations(range(12), 2)
             if np.isclose(np.linalg.norm(ico[i] - ico[j]), 2)]
    vertices = np.array([(ico[i] + ico[j]) / 2 for i, j in edges])
    faces = [face_cycle(vertices, np.flatnonzero(np.isclose(vertices @ normal,
                         max(vertices @ normal))), normal) for normal in ico]
    seen = set()
    for normal in ico:
        indices = np.flatnonzero(np.isclose(vertices @ normal, 0))
        key = tuple(indices)
        if key not in seen:
            faces.append(face_cycle(vertices, indices, normal))
            seen.add(key)
    return model("singular_nonorientable_g14_v30_small_dodecahemidodecahedron",
                 vertices, faces, 14, "Small dodecahemidodecahedron",
                 "Appendix B.3, Higher genus and vertex counts")


def drill_opposite_squares(base, sides, radius=.18):
    """Insert a regular prism through opposite square faces, without new
    subdivision vertices. Four simple planar polygons surround each opening.
    Its interior may cross other sheets, as allowed in Appendix B.3.
    """
    vertices = base["vertices"]
    faces = base["faces"]
    first = next(i for i, face in enumerate(faces) if len(face) == 4)
    square = faces[first]
    center = vertices[square].mean(axis=0)
    opposite = next(i for i, face in enumerate(faces)
                    if len(face) == 4 and i != first
                    and np.allclose(vertices[face].mean(axis=0), -center))
    normal = center / np.linalg.norm(center)
    u = vertices[square[0]] - center
    u /= np.linalg.norm(u)
    w = np.cross(normal, u)
    # Reorder both squares with the same in-plane orientation for vertical sides.
    outer = [face_cycle(vertices, faces[i], normal) for i in (first, opposite)]
    # Rotate the cycles to their first angle in this common coordinate system.
    for k, ids in enumerate(outer):
        c = center if k == 0 else -center
        angles = np.mod(np.arctan2((vertices[ids] - c) @ w,
                                  (vertices[ids] - c) @ u) + 1e-12, 2 * PI)
        outer[k] = [ids[i] for i in np.argsort(angles)]
    rings = []
    result_vertices = list(vertices)
    for c in (center, -center):
        rings.append(list(range(len(result_vertices), len(result_vertices) + sides)))
        result_vertices.extend([c + radius * (np.cos(2 * PI * j / sides) * u +
                                               np.sin(2 * PI * j / sides) * w)
                                for j in range(sides)])
    result_faces = [face[:] for i, face in enumerate(faces)
                    if i not in (first, opposite)]
    for ids, ring in zip(outer, rings):
        anchors = [int(np.floor(j * sides / 4 + .5)) for j in range(4)] + [sides]
        for j in range(4):
            result_faces.append([ids[j], ids[(j + 1) % 4]] +
                                [ring[k % sides] for k in range(anchors[j + 1],
                                                               anchors[j] - 1, -1)])
    for j in range(sides):
        k = (j + 1) % sides
        result_faces.append([rings[0][j], rings[0][k], rings[1][k], rings[1][j]])
    genus = base["genus"] + 2
    n = len(result_vertices)
    name = (f"singular_nonorientable_g{genus}_v{n}_" +
            ("cubohemioctahedron_hexagonal_drilling" if sides == 6 else
             "small_rhombihexahedron_square_drilling"))
    return model(name, result_vertices, result_faces, genus,
                 base["construction"] + f" with one {sides}-gonal tunnel",
                 "Appendix B.3 (singular drilling), using Appendix A.2",
                 parameters={"tunnel_sides": sides, "tunnel_radius": radius,
                             "base_model": base["name"]})


def build_models():
    """Finite catalogue covering every explicit singular construction family."""
    cubo = cubohemioctahedron()
    rhombi = small_rhombihexahedron()
    return [tetrahemihexahedron(), klein_bottle(), singular_torus(),
            *[antiprism_and_bands(g) for g in (3, 5, 7, 9, 11)],
            core_and_bands(3), core_and_bands(4), core_and_bands(13),
            historical_cuboctahedron_and_bands(), cubo,
            drill_opposite_squares(cubo, 6), rhombi,
            drill_opposite_squares(rhombi, 4), small_dodecahemidodecahedron()]
