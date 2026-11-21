"""Verify local injectivity of the Appendix A Möbius band over its full parameter domain.

Run python scripts/verify_mobius_links.py with SymPy installed. Calculations
use exact rational functions; no point of the parameter domain is sampled.
"""
from functools import cache
from itertools import combinations

from sympy import Matrix, Poly, cancel, symbols


b, eps, u, v = symbols("b eps u v")


def coefficient_sign(polynomial):
    coefficients = Poly(polynomial, u, v).coeffs()
    if all(c > 0 for c in coefficients):
        return 1
    if all(c < 0 for c in coefficients):
        return -1
    return 0


@cache
def certified_sign(expr):
    """Return a proved strict sign, or zero if the coefficient test is inconclusive."""
    # u,v > 0 parametrise precisely 0 < b < 2, 0 < eps/(b+1) < 2b/3.
    expr = expr.subs(eps, (1 + b) * 2 * b * v / (3 * (1 + v)))
    expr = cancel(expr.subs(b, 2 * u / (1 + u)))
    num, den = expr.as_numer_denom()
    return coefficient_sign(num) * coefficient_sign(den)


def band():
    # Apply (x,y,w) -> (x,sqrt(3)*y,w/sqrt(z)); z>0. The conjugate
    # of the 120-degree rotation has rational entries in these coordinates.
    rotation = Matrix([[-1, -1, 0], [3, -1, 0], [0, 0, 2]]) / 2
    generators = {
        "A": Matrix([1, 0, 0]), "B": Matrix([b, -b, 1]),
        "C": Matrix([b / 3, -b - eps, 1 + eps / (1 + b)]),
    }
    vertices = {name + str(i): rotation**i * point
                for name, point in generators.items() for i in range(3)}
    faces = []
    for i in range(3):
        j, m = (i + 1) % 3, (i - 1) % 3
        ai, aj, am = f"A{i}", f"A{j}", f"A{m}"
        bi, bj, ci = f"B{i}", f"B{j}", f"C{i}"
        faces.extend([(bi, aj, bj), (bi, bj, ci), (bi, ci, ai),
                      (am, bj, ci), (ci, ai, am)])
    # Check the stronger statement for the band capped with its boundary
    # triangle. Its link at A0 is now a cycle, as at B0 and C0.
    faces.append(("A0", "A1", "A2"))
    return vertices, faces


def verify_link(vertices, faces, centre, order, normal):
    edges = {frozenset(v for v in face if v != centre)
             for face in faces if centre in face}
    cycle = list(zip(order, order[1:] + order[:1]))
    if len(order) != len(set(order)) or edges != set(map(frozenset, cycle)):
        raise ValueError(f"{centre}: the given order is not the vertex link")
    rays = {p: vertices[p] - vertices[centre] for p in order}
    for p, ray in rays.items():
        if certified_sign(Matrix(normal).dot(ray)) != 1:
            raise ValueError(f"{centre}: chart denominator is not positive at {p}")

    @cache
    def orientation(a, c, d):
        # Projection onto normal.ray = 1 divides this homogeneous determinant
        # by positive factors and therefore preserves its sign.
        return certified_sign(rays[a].dot(rays[c].cross(rays[d])))

    count = 0
    for (a, c), (d, f) in combinations(cycle, 2):
        if len({a, c, d, f}) != 4:
            continue
        separated = (orientation(a, c, d) * orientation(a, c, f) == 1 or
                     orientation(d, f, a) * orientation(d, f, c) == 1)
        if not separated:
            raise ValueError(f"{centre}: edges {(a, c)} and {(d, f)} were not separated")
        count += 1
    # With at least four edges, an overlap of adjacent edges puts an endpoint
    # on a nonincident edge. Thus pairwise separation above also rules it out.
    print(f"{centre}: {len(cycle)} rays in one chart, {count} nonadjacent pairs separated")


def main():
    vertices, faces = band()
    verify_link(vertices, faces, "A0", ["A2", "C0", "B0", "B2", "C1", "A1"], [-1, 0, 5])
    verify_link(vertices, faces, "B0", ["A1", "B1", "C0", "A0", "B2", "C2"], [-1, 0, -1])
    verify_link(vertices, faces, "C0", ["B0", "B1", "A2", "A0"], [0, 0, -1])
    print("All three vertex-link orbits are simple throughout the parameter domain.")


if __name__ == "__main__":
    main()

