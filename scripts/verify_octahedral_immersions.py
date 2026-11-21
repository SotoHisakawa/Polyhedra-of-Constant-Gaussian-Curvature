"""Verify the two octahedral immersions in the paper's section on vertex bounds.

Run from the repository root with python scripts/verify_octahedral_immersions.py.
The only inputs are the common vertex orbit and the cyclic face representatives.
"""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import prod

from mesh_geometry import check_vertex_links, edge_incidence, orient_faces


FACE_REPRESENTATIVES = {
    18: [(0, 1, 4), (0, 1, 12), (0, 5, 8), (0, 6, 13)],
    24: [(0, 1, 3, 2), (0, 1, 4), (0, 5, 7), (0, 7, 13), (0, 11, 18)],
}


def cube_rotations():
    for order in permutations(range(3)):
        parity = (-1) ** sum(i > j for i, j in combinations(order, 2))
        for signs in product((-1, 1), repeat=3):
            if prod(signs) == parity:
                yield order, signs


def rotate(point, rotation):
    order, signs = rotation
    return tuple(s * point[i] for i, s in zip(order, signs))


def canonical_cycle(cycle):
    """Identify cyclic shifts and reversals, retaining the sides of a square."""
    return min(seq[i:] + seq[:i] for seq in (cycle, cycle[::-1]) for i in range(len(seq)))


def difference(p, q):
    return tuple(x - y for x, y in zip(p, q))


def dot(p, q):
    return sum(x * y for x, y in zip(p, q))


def cross(p, q):
    x, y, z = p
    a, b, c = q
    return y*c - z*b, z*a - x*c, x*b - y*a


def check_face(points):
    p, q, r = points[:3]
    edge = difference(q, p)
    if cross(edge, difference(r, p)) == (0, 0, 0):
        raise ValueError(f"Degenerate face: {points}")
    if len(points) == 4:
        side = difference(points[3], p)
        if difference(r, q) != side or dot(edge, side) != 0 or dot(edge, edge) != dot(side, side):
            raise ValueError(f"Face is not a square in the given cyclic order: {points}")


def check_geometric_link(vertices, faces):
    """Intersect the edge rays at vertex 0 with the plane (-p_0) . x = 1."""
    rays = [difference(vertex, vertices[0]) for vertex in vertices]
    if any(dot(vertices[0], ray) >= 0 for ray in rays[1:]):
        raise ValueError("A ray fails to meet the link plane in the positive direction")

    def det(i, j, k):
        return dot(cross(rays[i], rays[j]), rays[k])

    sides = []
    for face in faces:
        if 0 in face:
            i = face.index(0)
            sides.append((face[i - 1], face[(i + 1) % len(face)]))

    # The positive projection denominators do not affect determinant signs.
    # At adjacent sides check noncollinearity. At nonadjacent sides require
    # both endpoints of one side to lie strictly to one side of the other line.
    turns, separations = [], []
    for (a, b), (c, d) in combinations(sides, 2):
        indices = sorted({a, b, c, d})
        if len(indices) == 3:
            turns.append(abs(det(*indices)))
        elif len(indices) == 4:
            separations.append(max(det(a, b, c) * det(a, b, d),
                                   det(c, d, a) * det(c, d, b)))
        else:
            raise ValueError("Repeated side in the vertex link")
    if min(turns) <= 0 or min(separations) <= 0:
        raise ValueError("The vertex link failed the separation test")
    return min(turns), min(separations)


def verify(genus, representatives):
    rotations = list(cube_rotations())
    vertices = sorted({rotate((1, 3, 7), rotation) for rotation in rotations})
    if len(vertices) != 24:
        raise ValueError("Expected 24 distinct vertices")
    labels = {point: i for i, point in enumerate(vertices)}
    faces = set()
    for face in representatives:
        check_face([vertices[i] for i in face])
        for rotation in rotations:
            faces.add(canonical_cycle(tuple(labels[rotate(vertices[i], rotation)] for i in face)))
    faces = sorted(faces)
    check_vertex_links(len(vertices), faces)
    _, orientable, components = orient_faces(faces)  # Also checks two faces at every edge.
    edges = edge_incidence(faces)
    chi = len(vertices) - len(edges) + len(faces)
    if orientable or components != 1 or chi != 2 - genus:
        raise ValueError(f"Faces do not give the claimed non-orientable surface N_{genus}")

    # The rotations preserve all faces and act transitively on the vertices.
    # One embedded link thus gives an immersion. All defects agree, so the
    # Gauss--Bonnet theorem determines their exact value.
    turn, separation = check_geometric_link(vertices, faces)
    print(f"N_{genus}: V={len(vertices)}, E={len(edges)}, F={len(faces)}, "
          f"chi={chi}, defect/pi={Fraction(2 * chi, len(vertices))}")
    print(f"  minimum absolute link turn determinant={turn}, "
          f"minimum nonadjacent side separation product={separation}")


if __name__ == "__main__":
    for genus, representatives in FACE_REPRESENTATIVES.items():
        verify(genus, representatives)
