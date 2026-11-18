"""Create closed printable thickness models without changing the source surface.

The input is the exact triangulation of a mathematical surface.  The output
approximates the boundary of the union of radius ``thickness_mm / 2`` balls
centred on that surface.  Intersecting sheets are fused by taking the minimum
unsigned point-to-triangle distance, then a closed isosurface is extracted.
This is a sampled print aid, not a constant-curvature polyhedron: its vertex
count, defects and even topology need not equal those of the source surface.

Requires NumPy, SciPy and scikit-image.  All operations use ordinary arrays;
no external executable or native boolean-engine package is needed.
"""

import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from skimage.measure import marching_cubes


def _triangle_distance_squared(points, triangle):
    """Exact unsigned squared distances from N points to one nondegenerate triangle."""
    a, b, c = triangle
    ab, ac = b - a, c - a
    normal = np.cross(ab, ac)
    normal2 = float(normal @ normal)
    if normal2 <= 1e-24:
        raise ValueError("Cannot thicken a degenerate source triangle")
    pa = points - a
    plane_numerator = pa @ normal
    distance2 = plane_numerator * plane_numerator / normal2
    ab2, ac2, abac = float(ab @ ab), float(ac @ ac), float(ab @ ac)
    paab, paac = pa @ ab, pa @ ac
    # Barycentric coordinates of the orthogonal projection onto the plane.
    determinant = normal2
    bary_b = (ac2 * paab - abac * paac) / determinant
    bary_c = (ab2 * paac - abac * paab) / determinant
    outside = (bary_b < 0.0) | (bary_c < 0.0) | (bary_b + bary_c > 1.0)
    if np.any(outside):
        outside_points = points[outside]
        edge_distance2 = np.full(len(outside_points), np.inf)
        for first, last in ((a, b), (b, c), (c, a)):
            edge = last - first
            relative = outside_points - first
            position = np.clip((relative @ edge) / (edge @ edge), 0.0, 1.0)
            offset = relative - position[:, None] * edge
            edge_distance2 = np.minimum(edge_distance2, np.einsum("ij,ij->i", offset, offset))
        distance2[outside] = edge_distance2
    return distance2


def _verify_closed(vertices, faces):
    """Check edge incidence, face winding and connected vertex links."""
    n_faces = len(faces)
    if n_faces == 0 or not np.all(np.isfinite(vertices)):
        raise RuntimeError("The extracted surface is empty or non-finite")
    face_points = vertices[faces]
    twice_areas = np.linalg.norm(np.cross(face_points[:, 1] - face_points[:, 0],
                                         face_points[:, 2] - face_points[:, 0]), axis=1)
    if np.any(twice_areas <= 1e-12):
        raise RuntimeError("The extracted surface contains degenerate triangles")
    directed_edges = faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2)
    sorted_edges = np.sort(directed_edges, axis=1)
    order = np.lexsort((sorted_edges[:, 1], sorted_edges[:, 0]))
    ordered_edges = sorted_edges[order]
    if len(ordered_edges) % 2 or not np.all(ordered_edges[::2] == ordered_edges[1::2]):
        raise RuntimeError("The printable surface has open edges")
    if len(ordered_edges) > 2 and np.any(np.all(ordered_edges[1:-1:2] == ordered_edges[2::2], axis=1)):
        raise RuntimeError("The printable surface has nonmanifold edges")
    directed_pairs = directed_edges[order].reshape(-1, 2, 2)
    if not np.all(directed_pairs[:, 0] == directed_pairs[:, 1, ::-1]):
        raise RuntimeError("The printable surface has inconsistent face winding")

    # A corner represents a face incident to its vertex.  Across each edge,
    # connect the two corners belonging to the same endpoint.  Every vertex
    # link must be a single cycle, not several cycles pinched at one point.
    starts = np.arange(3 * n_faces, dtype=np.int64)
    ends = (starts // 3) * 3 + (starts % 3 + 1) % 3
    corners = np.column_stack((starts, ends))
    reverse = directed_edges[:, 0] > directed_edges[:, 1]
    corners[reverse] = corners[reverse, ::-1]
    corner_pairs = corners[order].reshape(-1, 2, 2)
    left, right = corner_pairs[:, 0].ravel(), corner_pairs[:, 1].ravel()
    graph = coo_matrix((np.ones(len(left), dtype=np.uint8), (left, right)),
                       shape=(3 * n_faces, 3 * n_faces)).tocsr()
    link_count = connected_components(graph, directed=False, return_labels=False)
    if link_count != len(vertices):
        raise RuntimeError("The printable surface has disconnected vertex links")

    face_pairs = order.reshape(-1, 2) // 3
    face_graph = coo_matrix((np.ones(len(face_pairs), dtype=np.uint8),
                             (face_pairs[:, 0], face_pairs[:, 1])),
                            shape=(n_faces, n_faces)).tocsr()
    shell_count, shell_labels = connected_components(face_graph, directed=False, return_labels=True)
    face_volumes = np.einsum("ij,ij->i", face_points[:, 0],
                            np.cross(face_points[:, 1], face_points[:, 2])) / 6.0
    shell_volumes = np.bincount(shell_labels, weights=face_volumes)
    return ({"vertices": len(vertices), "triangles": n_faces,
             "edges": len(ordered_edges) // 2, "closed": True,
             "manifold_vertex_links": True, "consistent_winding": True,
             "boundary_components": int(shell_count),
             "positive_boundary_components": int(np.sum(shell_volumes > 0)),
             "cavity_boundary_components": int(np.sum(shell_volumes < 0)),
             "shell_volumes_mm3": shell_volumes.tolist(),
             "signed_volume_mm3": float(np.sum(face_volumes))}, shell_labels)


def _float32_mesh(vertices, faces):
    """Weld in the coordinate precision actually stored by binary STL.

    A near-collapsed marching-cubes feature can contain two different
    double-precision points which become identical in STL.  Its two sides
    then become oppositely oriented duplicate triangles; cancel both sides
    of that zero-thickness feature instead of retaining a nonmanifold edge.
    Closure and all vertex links are checked separately after this operation.
    """
    vertices32 = np.asarray(vertices, dtype=np.float32)
    welded, inverse = np.unique(vertices32, axis=0, return_inverse=True)
    faces = inverse[faces]
    repeated = (faces[:, 0] == faces[:, 1]) | (faces[:, 1] == faces[:, 2]) | (faces[:, 2] == faces[:, 0])
    faces = faces[~repeated]
    _, group, counts = np.unique(np.sort(faces, axis=1), axis=0, return_inverse=True, return_counts=True)
    keep = np.ones(len(faces), dtype=bool)
    duplicate_count = 0
    for index in np.flatnonzero(counts > 1):
        ids = np.flatnonzero(group == index)
        triangles = faces[ids]
        inversions = ((triangles[:, 0] > triangles[:, 1]).astype(int)
                      + (triangles[:, 0] > triangles[:, 2])
                      + (triangles[:, 1] > triangles[:, 2]))
        signs = 1 - 2 * (inversions % 2)
        balance = int(np.sum(signs))
        keep[ids] = False
        if balance:
            # A remaining oriented copy is retained; only the subsequent
            # manifold checks can establish whether this weld is valid.
            keep[ids[np.flatnonzero(signs == np.sign(balance))[0]]] = True
        duplicate_count += len(ids) - int(balance != 0)
    faces = faces[keep]
    used, compact = np.unique(faces, return_inverse=True)
    result_vertices = welded[used].astype(float)
    result_faces = compact.reshape(-1, 3)
    return result_vertices, result_faces, {
        "float32_welded_vertices": len(vertices32) - len(welded),
        "float32_collapsed_triangles_removed": int(np.sum(repeated)),
        "float32_coincident_triangles_removed": duplicate_count,
    }


def make_printable(vertices, triangles, size_mm=100.0, thickness_mm=1.6, pitch_mm=0.5):
    """Return (vertices_mm, triangular_faces, report) for a closed thickness model.

    ``size_mm`` is the largest extent of the unthickened source; the final
    result is approximately ``thickness_mm`` larger.  ``pitch_mm`` is the
    isotropic sampling spacing.  Use a pitch no larger than one third of
    the thickness; the default is suitable for a 100 mm demonstration.
    The closed isosurface is validated before it is returned.  Disjoint
    enclosed voids can contribute additional boundary components.
    """
    start_time = time.monotonic()
    vertices = np.asarray(vertices, dtype=float)
    triangles = np.asarray(triangles, dtype=np.int64)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or not len(vertices):
        raise ValueError("vertices must be a nonempty N x 3 array")
    if triangles.ndim != 2 or triangles.shape[1] != 3 or not len(triangles):
        raise ValueError("triangles must be a nonempty F x 3 array")
    if not np.all(np.isfinite(vertices)) or np.min(triangles) < 0 or np.max(triangles) >= len(vertices):
        raise ValueError("The source has non-finite coordinates or invalid indices")
    if min(size_mm, thickness_mm, pitch_mm) <= 0 or not np.all(np.isfinite([size_mm, thickness_mm, pitch_mm])):
        raise ValueError("size, thickness and sampling pitch must be positive and finite")
    if pitch_mm > thickness_mm / 3.0 + 1e-12 * max(1.0, thickness_mm):
        raise ValueError("Sampling pitch must not exceed thickness / 3")
    span = np.ptp(vertices, axis=0)
    if np.max(span) <= 0:
        raise ValueError("The source has zero extent")
    source_center = (vertices.min(axis=0) + vertices.max(axis=0)) / 2.0
    scale = size_mm / float(np.max(span))
    scaled_vertices = (vertices - source_center) * scale
    radius = thickness_mm / 2.0
    # Outside each expanded triangle bounding box the exact distance is
    # greater than the desired level.  The extra two voxels make the global
    # grid boundary strictly exterior and retain a safe interpolation band.
    margin = radius + 2.0 * pitch_mm
    origin = scaled_vertices.min(axis=0) - margin
    shape = np.ceil((scaled_vertices.max(axis=0) + margin - origin) / pitch_mm).astype(int) + 1
    if np.prod(shape, dtype=np.int64) > 100_000_000:
        raise ValueError("Requested sampling grid exceeds 100 million voxels; increase pitch or reduce size")
    axes = [origin[i] + np.arange(shape[i]) * pitch_mm for i in range(3)]
    distances = np.full(tuple(shape), margin ** 2, dtype=np.float32)
    for triangle in scaled_vertices[triangles]:
        lower = np.maximum(0, np.floor((triangle.min(axis=0) - margin - origin) / pitch_mm).astype(int))
        upper = np.minimum(shape, np.ceil((triangle.max(axis=0) + margin - origin) / pitch_mm).astype(int) + 1)
        cross_section = max(1, (upper[1] - lower[1]) * (upper[2] - lower[2]))
        block_width = max(1, 200_000 // cross_section)
        for left in range(lower[0], upper[0], block_width):
            right = min(upper[0], left + block_width)
            block_axes = (axes[0][left:right], axes[1][lower[1]:upper[1]], axes[2][lower[2]:upper[2]])
            points = np.stack(np.meshgrid(*block_axes, indexing="ij"), axis=-1).reshape(-1, 3)
            values = _triangle_distance_squared(points, triangle)
            current = distances[left:right, lower[1]:upper[1], lower[2]:upper[2]]
            np.minimum(current, values.reshape(current.shape), out=current)

    level = radius ** 2
    for i in range(3):
        if np.any(np.take(distances, [0, shape[i] - 1], axis=i) <= level):
            raise RuntimeError("Insufficient padding around the printable model")
    failures = []
    for relative_adjustment in (0.0, 2.0 ** -16, -2.0 ** -16, 2.0 ** -12, -2.0 ** -12):
        result_vertices, result_faces, _, _ = marching_cubes(
            distances, level=level * (1.0 + relative_adjustment), spacing=(pitch_mm,) * 3,
            gradient_direction="descent", allow_degenerate=False, method="lewiner",
        )
        result_vertices = result_vertices.astype(float) + origin
        result_faces = result_faces.astype(np.int64)
        result_vertices, result_faces, rounding_report = _float32_mesh(result_vertices, result_faces)
        try:
            report, shell_labels = _verify_closed(result_vertices, result_faces)
        except RuntimeError as error:
            failures.append(str(error))
        else:
            break
    else:
        raise RuntimeError("No manifold float32 isosurface could be extracted: " + "; ".join(failures))
    if report["signed_volume_mm3"] < 0:
        result_faces = result_faces[:, ::-1].copy()
        report, shell_labels = _verify_closed(result_vertices, result_faces)
    if report["signed_volume_mm3"] <= 0:
        raise RuntimeError("The printable model does not enclose positive volume")
    # Linear interpolation of squared distances can isolate microscopic
    # cavity shells.  Fill only negative-volume cavities smaller than one
    # sampling voxel; retain the larger voids that belong to the geometry.
    shell_volumes = np.array(report["shell_volumes_mm3"])
    tiny_cavities = (shell_volumes < 0) & (shell_volumes > -pitch_mm ** 3)
    tiny_count = int(np.sum(tiny_cavities))
    filled_volume = float(-np.sum(shell_volumes[tiny_cavities]))
    if tiny_count:
        result_faces = result_faces[~tiny_cavities[shell_labels]]
        used, compact_indices = np.unique(result_faces, return_inverse=True)
        result_vertices = result_vertices[used]
        result_faces = compact_indices.reshape(-1, 3)
        report, _ = _verify_closed(result_vertices, result_faces)
    report.update({
        "method": "unsigned triangle distance union; Lewiner marching cubes",
        "source_size_mm": float(size_mm), "source_scale": scale,
        "source_center": source_center.tolist(),
        "thickness_mm": float(thickness_mm), "pitch_mm": float(pitch_mm),
        "effective_thickness_mm": float(thickness_mm * np.sqrt(1.0 + relative_adjustment)),
        "isosurface_level_relative_adjustment": relative_adjustment,
        "grid_shape": shape.tolist(), "extent_mm": np.ptp(result_vertices, axis=0).tolist(),
        "source_vertices": len(vertices), "source_triangles": len(triangles),
        "preserves_source_vertex_counts_and_defects": False,
        "subvoxel_cavities_filled": tiny_count,
        "subvoxel_cavities_filled_mm3": filled_volume,
        "elapsed_seconds": round(time.monotonic() - start_time, 3),
    })
    report.update(rounding_report)
    return result_vertices, result_faces, report
