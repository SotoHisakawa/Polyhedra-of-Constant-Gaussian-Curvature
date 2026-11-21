"""Verify the small-genus constructions in the section on vertex bounds.

Run python scripts/check_vertex_bound_seams.py with mpmath installed.
The turns and the genus-eight data are rational. The trigonometric and
genus-four band checks use outward-rounded intervals, including the full
perturbation range 0 <= epsilon <= 1/100000.
"""

from fractions import Fraction

from mpmath import iv


def require(condition, description):
    if condition is not True:
        raise ArithmeticError(description)


def enclose(value):
    """Convert a rational number to an outward-rounded interval."""
    q = Fraction(value)
    return iv.mpf(q.numerator) / q.denominator


def intervals(lower, upper, count):
    """Cover the full closed interval by rational subdivisions."""
    left, right = Fraction(lower), Fraction(upper)
    step = (right - left) / count
    for i in range(count):
        yield iv.mpf([enclose(left + i * step).a, enclose(left + (i + 1) * step).b])


def difference(p, q):
    return tuple(x - y for x, y in zip(p, q))


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def cross(u, v):
    x, y, z = u
    a, b, c = v
    return (y*c - z*b, z*a - x*c, x*b - y*a)


def band_defect(b, z, epsilon, orbit):
    """Angle defect on the indicated interior orbit of M(b,z,epsilon).

    The link orders below are obtained from the fifteen triangles in
    Appendix A. Evaluate each face angle as atan2(|u x v|, u.v).
    """
    h, root3, half = iv.sqrt(z), iv.sqrt(3), enclose("1/2")
    height_c = (1 + epsilon / (1 + b)) * h
    a0, a1, a2 = (1, 0, 0), (-half, root3/2, 0), (-half, -root3/2, 0)
    b0, b1, b2 = (b, -b/root3, h), (0, 2*b/root3, h), (-b, -b/root3, h)
    c0 = (b/3, -(b + epsilon)/root3, height_c)
    c2 = (-2*b/3 - epsilon/2, epsilon/(2*root3), height_c)
    if orbit == "B":
        centre, neighbours = b0, (a1, b1, c0, a0, b2, c2)
    else:
        centre, neighbours = c0, (b0, b1, a2, a0)
    rays = [difference(point, centre) for point in neighbours]
    defect = 2 * iv.pi
    for u, v in zip(rays, rays[1:] + rays[:1]):
        area = iv.sqrt(sum(x**2 for x in cross(u, v)))
        require(area > 0, "Degenerate triangle in the band")
        defect -= iv.atan2(area, dot(u, v))
    return defect


def check_genus_four_band(xi):
    b_min, b_max = "0.595", "0.610"
    z_min, z_max = "0.052", "0.064"
    epsilon = iv.mpf([0, "0.00001"])
    sides = (
        ("left", "B", b_min, z_min, z_max, -53, -6),
        ("right", "B", b_max, z_min, z_max, 4, 50),
        ("bottom", "C", z_min, b_min, b_max, -37, -1),
        ("top", "C", z_max, b_min, b_max, 2, 37),
    )
    for side, orbit, fixed, low, high, lower_bound, upper_bound in sides:
        for varying in intervals(low, high, 12):
            b, z = (enclose(fixed), varying) if orbit == "B" else (varying, enclose(fixed))
            residual = 1000 * (band_defect(b, z, epsilon, orbit) + iv.pi/6)
            require(residual in iv.mpf([lower_bound, upper_bound]), f"N4 {side} sign")
        print(f"N4 {side}: 1000 * (delta_{orbit} + pi/6) in [{lower_bound},{upper_bound}]")

    b, z = iv.mpf([b_min, b_max]), iv.mpf([z_min, z_max])
    # For epsilon < 2b/3 the threefold rotations of each B and C vertex
    # have X >= -b and Z >= sqrt(z). The common rim edge has X=-1/2,Z=0.
    require(epsilon < 2*b/3, "N4 perturbation bound")
    require(0 < xi and xi < iv.pi/2, "N4 angle to the bisector")
    distance = (enclose("1/2") - b) * iv.sin(xi) + iv.sqrt(z) * iv.cos(xi)
    require(distance > 0, "N4 bands must lie on opposite sides of the symmetry plane")
    print("N4 positive separation for the two reflected bands")


def check_triangular_prisms():
    for genus, vertex_count, expected in (
        (3, 17, (20, 43, 39, 51)),
        (4, 24, (14, 17, 17, 24)),
    ):
        # Angles and exterior turns are rational multiples of pi.
        d = Fraction(2*(2 - genus), vertex_count)
        window = 1 + 2*d
        pyramid_b = Fraction(1, 6) + 7*d/4
        pyramid_a = pyramid_b + Fraction(2, 3) + 2*d
        shared = window + pyramid_b - d
        turns = ((pyramid_a, shared, window) if genus == 3
                 else (2*pyramid_a - d, shared, shared))
        a, b, c, denominator = expected
        require(turns == tuple(Fraction(x, denominator) for x in (a, b, c)), "Base turns")
        require(sum(turns) == 2 and all(0 < e < 1 for e in turns), "Convex triangle")
        angle_a, angle_b, angle_c = (iv.pi * enclose(1-e) for e in turns)
        theta = iv.pi * enclose(Fraction(1, 6) - 3*d/4)
        w = 2 * iv.sin(theta)
        length_bc = 2*w * iv.sin(angle_a) / iv.sin(angle_c)
        window_angle = iv.pi * enclose(-d/2)
        require(length_bc > 2*w and 2*w > 2, "Window width")
        require(angle_b > window_angle and angle_c > window_angle, "Window footprint")
        alpha = iv.atan2(iv.sqrt(3-w**2), w)
        if genus == 3:
            require(alpha < iv.pi * enclose(turns[0]), "Pyramid supporting plane")
        else:
            require(turns[1] == turns[2], "Mirror symmetry of the base")
            check_genus_four_band(angle_a/2 + alpha)
        print(f"N{genus} base turns (pi/{denominator}): {a}, {b}, {c}; window fits")


def check_octagonal_mouths():
    # The inward-facing wall of the right-hand window has oriented normal
    # (cos v,-sin v,0), the outer square side (1,0,0). Their sum bisects them.
    v = iv.pi / 8
    tangent = iv.tan(v)
    half = enclose("1/2")
    y_on_wall = half
    x_on_wall = 1 - tangent * (1 - y_on_wall)
    waist = 1 - iv.sin(v)/2
    require(0 < y_on_wall and y_on_wall < waist, "N6 first mouth on the window pentagon")
    # Intersect the ray along the sum of normals with x=1.
    step = (1-x_on_wall) / (1+iv.cos(v))
    target_y = y_on_wall - step * iv.sin(v)
    require(step > 0 and -1 < target_y and target_y < 1, "N6 second mouth on outer rectangle")
    print("N6 reflected octagons: both mouth centres in the face interiors")


def check_triangular_mouths():
    def r(v):
        x, y, z = v
        return -y, z, -x

    def s(v):
        x, y, z = v
        return -x, -y, z

    p = (1, 2, 6)
    q = s(p)
    face_a = (p, r(p), r(r(p)))
    face_b = tuple(s(v) for v in (q, r(q), r(r(q))))
    a, b = (0, 2, 5), (1, 1, 5)
    normal_a, normal_b = (1, -1, -1), (1, -1, 1)
    for face, normal, centre, weights in (
        (face_a, normal_a, a, (30, 3, 4)),
        (face_b, normal_b, b, (41, 5, 3)),
    ):
        require(all(dot(normal, difference(v, centre)) == 0 for v in face), "N8 face normal")
        require(cross(difference(face[1], face[0]), difference(face[2], face[0])) != (0, 0, 0),
                "N8 nondegenerate face")
        require(all(w > 0 for w in weights), "N8 positive barycentric weights")
        for column, component in zip(zip(*face), centre):
            require(dot(weights, column) == sum(weights)*component, "N8 mouth centre")
    # The two r-invariant faces have opposite normals across their common
    # axis. The oriented normal on the s-image of the second is s(-normal_a).
    require(s(tuple(-x for x in normal_a)) == normal_b, "N8 compatible normals")
    require((a[1]-1, a[0]+1, a[2]) == b, "N8 reflection of the mouth centres")
    require((normal_a[1], normal_a[0], normal_a[2]) == tuple(-x for x in normal_b),
            "N8 reflection of the normals")
    direction = difference(b, a)
    require(dot(normal_a, direction) > 0 and dot(normal_b, direction) > 0, "N8 transversality")
    print("N8 reflected triangles: exact face interiors, normals and reflection")


if __name__ == "__main__":
    iv.dps = 40
    check_triangular_prisms()
    check_octagonal_mouths()
    check_triangular_mouths()
