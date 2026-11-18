"""Independent numerical triangle-intersection checks for embedded examples.

Run ``python scripts/check_embedded_intersections.py``; optionally supply
``--output PATH`` to save the report. The checks compare every pair whose
axis-aligned boxes overlap. Triangle intersections are allowed only on the
shared indexed edge or vertex. Coplanar overlaps are checked by polygon
clipping, and noncoplanar intersections by edge/plane intersections.

This is a floating-point check at a stated tolerance, supplementing the
manuscript's embedding proofs. It is not an exact-arithmetic certificate.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def cross2(a, b):
    return a[0] * b[1] - a[1] * b[0]


def inside_triangle(point, triangle, normal, tolerance):
    for i in range(3):
        edge = triangle[(i + 1) % 3] - triangle[i]
        if np.dot(np.cross(edge, point - triangle[i]), normal) < -tolerance * np.linalg.norm(edge):
            return False
    return True


def clip_coplanar(first, second, normal, tolerance):
    """Return vertices of a possibly degenerate coplanar intersection."""
    drop = int(np.argmax(np.abs(normal)))
    keep = [i for i in range(3) if i != drop]
    subject = [p.copy() for p in first]
    clip = second
    xy = clip[:, keep]
    orientation = np.sign(cross2(xy[1] - xy[0], xy[2] - xy[0]))
    for i in range(3):
        if not subject:
            break
        a, b = xy[i], xy[(i + 1) % 3]
        edge = b - a
        threshold = tolerance * np.linalg.norm(edge)
        output = []
        previous = subject[-1]
        old = orientation * cross2(edge, previous[keep] - a)
        for current in subject:
            new = orientation * cross2(edge, current[keep] - a)
            old_inside, new_inside = old >= -threshold, new >= -threshold
            if old_inside != new_inside:
                # The exact supporting-line intersection; tolerance affects
                # classification, not the returned geometric coordinates.
                denominator = old - new
                if abs(denominator) > 1e-30:
                    t = np.clip(old / denominator, 0, 1)
                    output.append(previous + t * (current - previous))
            if new_inside:
                output.append(current)
            previous, old = current, new
        subject = output
    return subject


def triangle_intersection(first, second, normal_first, normal_second, tolerance):
    """Return intersection points spanning the pair's convex intersection."""
    da = (first - second[0]) @ normal_second
    db = (second - first[0]) @ normal_first
    if (np.all(da > tolerance) or np.all(da < -tolerance)
            or np.all(db > tolerance) or np.all(db < -tolerance)):
        return []
    if np.linalg.norm(np.cross(normal_first, normal_second)) < 1e-10:
        if np.max(np.abs(da)) > tolerance:
            return []
        return clip_coplanar(first, second, normal_first, tolerance)
    points = []
    for tri, other, distances, normal in (
            (first, second, da, normal_second),
            (second, first, db, normal_first)):
        for i in range(3):
            j = (i + 1) % 3
            if abs(distances[i]) <= tolerance and inside_triangle(tri[i], other, normal, tolerance):
                points.append(tri[i])
            if (distances[i] < 0 < distances[j]) or (distances[j] < 0 < distances[i]):
                t = distances[i] / (distances[i] - distances[j])
                p = tri[i] + t * (tri[j] - tri[i])
                if inside_triangle(p, other, normal, tolerance):
                    points.append(p)
    return points


def on_shared_simplex(point, common_points, tolerance):
    if len(common_points) == 1:
        return np.linalg.norm(point - common_points[0]) <= tolerance
    if len(common_points) == 2:
        a, b = common_points
        edge = b - a
        t = np.dot(point - a, edge) / np.dot(edge, edge)
        t = np.clip(t, 0, 1)
        return np.linalg.norm(point - (a + t * edge)) <= tolerance
    return False


def check_triangles(vertices, triangles, relative_tolerance=2e-10):
    vertices = np.asarray(vertices, dtype=float)
    triangles = np.asarray(triangles, dtype=int)
    points = vertices[triangles]
    normals = np.cross(points[:, 1] - points[:, 0], points[:, 2] - points[:, 0])
    lengths = np.linalg.norm(normals, axis=1)
    if np.any(lengths <= 1e-16):
        raise ValueError("Intersection check requires nondegenerate triangles")
    normals /= lengths[:, None]
    scale = max(float(np.ptp(vertices, axis=0).max()), 1e-12)
    tolerance = relative_tolerance * scale
    boxes_min = points.min(axis=1)
    boxes_max = points.max(axis=1)
    candidate_count = 0
    illegal = []
    for i in range(len(triangles)):
        candidates = np.flatnonzero(
            np.all(boxes_max[i] + tolerance >= boxes_min, axis=1)
            & np.all(boxes_max + tolerance >= boxes_min[i], axis=1))
        for j in candidates[candidates > i]:
            candidate_count += 1
            common = sorted(set(map(int, triangles[i])) & set(map(int, triangles[j])))
            intersection = triangle_intersection(points[i], points[j], normals[i], normals[j], tolerance)
            if intersection and (len(common) == 3 or
                                 any(not on_shared_simplex(p, vertices[common], 10 * tolerance)
                                     for p in intersection)):
                illegal.append([int(i), int(j)])
    return {"triangles": len(triangles), "candidate_pairs": candidate_count,
            "absolute_tolerance": tolerance,
            "intersection_count": len(illegal), "intersections": illegal,
            "passed": not illegal}


def check_model(model):
    from mesh_geometry import prepare_model
    prepared, _ = prepare_model(model)
    return dict(name=model["name"], **check_triangles(prepared["vertices"], prepared["triangles"]))


def self_test():
    """Exercise crossings and legitimate contact cases independently of models."""
    cases = [
        # Disjoint parallel triangles.
        ([[0, 0, 0], [2, 0, 0], [0, 2, 0], [0, 0, 1], [2, 0, 1], [0, 2, 1]],
         [[0, 1, 2], [3, 4, 5]], True),
        # Noncoplanar crossing with no common indexed vertex.
        ([[0, 0, 0], [2, 0, 0], [0, 2, 0], [.5, .5, -1], [.5, .5, 1], [1.5, .5, 0]],
         [[0, 1, 2], [3, 4, 5]], False),
        # Allowed noncoplanar shared edge.
        ([[0, 0, 0], [2, 0, 0], [0, 2, 0], [0, 0, 2]], [[0, 1, 2], [1, 0, 3]], True),
        # Allowed coplanar shared edge, with opposite interiors.
        ([[0, 0, 0], [2, 0, 0], [0, 2, 0], [0, -2, 0]], [[0, 1, 2], [1, 0, 3]], True),
        # Disallowed coplanar overlap, even though an edge is shared.
        ([[0, 0, 0], [2, 0, 0], [0, 2, 0], [.5, .5, 0]], [[0, 1, 2], [0, 1, 3]], False),
        # Allowed sole common vertex.
        ([[0, 0, 0], [2, 0, 0], [0, 2, 0], [-2, 0, 1], [0, -2, 1]],
         [[0, 1, 2], [0, 3, 4]], True),
        # Intersection beyond a common vertex.
        ([[0, 0, 0], [2, 0, 0], [0, 2, 0], [1, 1, 1], [1, 1, -1]],
         [[0, 1, 2], [0, 3, 4]], False),
        # Contact at an interior point is illegal without shared incidence.
        ([[0, 0, 0], [2, 0, 0], [0, 2, 0], [.5, .5, 0], [.5, .5, 1], [1, .5, 1]],
         [[0, 1, 2], [3, 4, 5]], False),
    ]
    for i, (vertices, triangles, passed) in enumerate(cases):
        result = check_triangles(vertices, triangles)
        if result["passed"] != passed:
            raise RuntimeError(f"Intersection checker self-test {i} failed: {result}")


def main():
    from embedded_models import build_models
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional JSON report path")
    args = parser.parse_args()
    self_test()
    results = []
    for model in build_models():
        result = check_model(model)
        results.append(result)
        print(f"{result['name']}: {result['candidate_pairs']} pairs, "
              f"{result['intersection_count']} forbidden intersections", flush=True)
    report = {"method": "Pairwise triangle intersections after AABB filtering",
              "precision": "floating-point, relative tolerance 2e-10",
              "shared_contact_tolerance_multiplier": 10,
              "self_tests": 8,
              "passed": all(result["passed"] for result in results),
              "models": results}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not report["passed"]:
        raise SystemExit("At least one embedded model has a forbidden intersection")


if __name__ == "__main__":
    main()
