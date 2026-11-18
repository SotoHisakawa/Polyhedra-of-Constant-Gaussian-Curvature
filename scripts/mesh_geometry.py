"""Polygon triangulation, closed-surface checks, and deterministic STL/OBJ IO.

The polygon presentation is authoritative: in a singular realisation, distinct
abstract diagonals can have the same spatial endpoints. STL cannot record that
distinction, or an orientation on a non-orientable surface.
"""
from collections import Counter, defaultdict, deque
from pathlib import Path
import os
import struct
import tempfile
import time

import numpy as np


def write_bytes_atomic(path, data):
    """Avoid partial artifacts and unnecessary writes in synced directories."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_bytes() == data:
        return
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
        for attempt in range(7):
            try:
                os.replace(temporary, path)
                break
            except OSError as error:
                if attempt == 6 or (error.errno not in (13,22,32) and getattr(error,"winerror",None) not in (5,32,33)):
                    raise
                # A sync client may briefly hold the destination open.
                time.sleep(.05 * 2**attempt)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def cross2(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def face_plane(points):
    """Return oriented orthonormal planar coordinates; reject nonplanar faces."""
    points = np.asarray(points, dtype=float)
    p = points - points[0]
    scale = max(float(np.max(np.linalg.norm(p, axis=1))), 1e-15)
    normal = np.sum(np.cross(p, np.roll(p, -1, axis=0)), axis=0)
    length = np.linalg.norm(normal)
    if length <= 1e-13 * scale**2:
        raise ValueError("Face has zero oriented area")
    normal /= length
    error = float(np.max(np.abs(p @ normal)))
    if error > 2e-10 * scale:
        raise ValueError(f"Nonplanar face: {error:g}")
    ex = p[1] / np.linalg.norm(p[1])
    ey = np.cross(normal, ex)
    return np.column_stack((p @ ex, p @ ey)), normal, error


def triangulate_face(vertices, face):
    """Ear clipping, including concave faces, without adding paper vertices."""
    if len(face) < 3 or len(set(face)) != len(face):
        raise ValueError(f"Invalid polygon indices: {face}")
    xy, _, _ = face_plane(vertices[face])
    scale = max(float(np.ptp(xy, axis=0).max()), 1e-15)
    tol = 1e-12 * scale**2
    # Reject crossings, including a vertex in a nonincident edge.
    for i in range(len(face)):
        a, b = xy[i], xy[(i + 1) % len(face)]
        for j in range(i + 1, len(face)):
            if j == i + 1 or (i == 0 and j == len(face) - 1):
                continue
            c, d = xy[j], xy[(j + 1) % len(face)]
            u, v = cross2(b-a, c-a), cross2(b-a, d-a)
            w, z = cross2(d-c, a-c), cross2(d-c, b-c)
            proper = ((u > tol and v < -tol) or (v > tol and u < -tol)) and ((w > tol and z < -tol) or (z > tol and w < -tol))
            touching = any(abs(t) <= tol and np.dot(p-q, p-r) <= tol for t, p, q, r in
                           [(u, c, a, b), (v, d, a, b), (w, a, c, d), (z, b, c, d)])
            if proper or touching:
                raise ValueError(f"Nonsimple polygon: {face}")
    remaining = list(range(len(face)))
    result = []
    while len(remaining) > 3:
        for k, curr in enumerate(remaining):
            prev, nxt = remaining[k-1], remaining[(k+1) % len(remaining)]
            a, b, c = xy[[prev, curr, nxt]]
            if cross2(b-a, c-a) <= tol:
                continue
            other = [i for i in remaining if i not in (prev, curr, nxt)]
            p = xy[other]
            inside = (cross2(b-a, p-a) >= -tol) & (cross2(c-b, p-b) >= -tol) & (cross2(a-c, p-c) >= -tol)
            if np.any(inside):
                continue
            result.append([face[prev], face[curr], face[nxt]])
            remaining.pop(k)
            break
        else:
            raise ValueError(f"Cannot triangulate polygon: {face}")
    result.append([face[i] for i in remaining])
    tri_area = sum(cross2(xy[face.index(b)]-xy[face.index(a)], xy[face.index(c)]-xy[face.index(a)]) / 2 for a,b,c in result)
    polygon_area = float(np.sum(cross2(xy, np.roll(xy, -1, axis=0))) / 2)
    if abs(tri_area - polygon_area) > 1e-10 * scale**2:
        raise ValueError("Triangulation area differs from polygon area")
    return result


def edge_incidence(faces):
    edges = defaultdict(list)
    for fi, face in enumerate(faces):
        for a, b in zip(face, face[1:] + face[:1]):
            if a == b:
                raise ValueError("Zero-length abstract edge")
            edges[tuple(sorted((a, b)))].append((fi, 1 if a < b else -1))
    return edges


def orient_faces(faces):
    """Return a coherent orientation when possible, plus connectedness."""
    faces = [list(map(int, f)) for f in faces]
    edges = edge_incidence(faces)
    adjacency = defaultdict(list)
    for pair in edges.values():
        if len(pair) != 2:
            raise ValueError(f"Polygon edge has {len(pair)} incident faces, expected 2")
        (a, da), (b, db) = pair
        adjacency[a].append((b, -da * db))
        adjacency[b].append((a, -da * db))
    signs, components, orientable = {}, 0, True
    for start in range(len(faces)):
        if start in signs:
            continue
        components += 1
        signs[start] = 1
        queue = deque([start])
        while queue:
            a = queue.popleft()
            for b, relative in adjacency[a]:
                sign = signs[a] * relative
                if b in signs:
                    orientable &= signs[b] == sign
                else:
                    signs[b] = sign
                    queue.append(b)
    if orientable:
        faces = [f if signs[i] == 1 else f[::-1] for i, f in enumerate(faces)]
    return faces, bool(orientable), components


def check_vertex_links(vertex_count, faces):
    """Every closed-surface vertex must have one cycle as its abstract link."""
    links = [defaultdict(list) for _ in range(vertex_count)]
    for f in faces:
        for k, i in enumerate(f):
            a, b = f[k-1], f[(k+1) % len(f)]
            links[i][a].append(b)
            links[i][b].append(a)
    for i, adjacency in enumerate(links):
        if not adjacency or any(len(v) != 2 for v in adjacency.values()):
            raise ValueError(f"Vertex {i} does not have a cycle as its link")
        seen, todo = set(), [next(iter(adjacency))]
        while todo:
            a = todo.pop()
            if a not in seen:
                seen.add(a)
                todo.extend(adjacency[a])
        if len(seen) != len(adjacency):
            raise ValueError(f"Vertex {i} has a disconnected link")


def triangle_defects(vertices, triangles):
    vertices = np.asarray(vertices, dtype=float)
    triangles = np.asarray(triangles, dtype=int)
    angles = np.zeros(len(vertices))
    for k in range(3):
        u = vertices[triangles[:, (k+1) % 3]] - vertices[triangles[:, k]]
        v = vertices[triangles[:, (k+2) % 3]] - vertices[triangles[:, k]]
        theta = np.arctan2(np.linalg.norm(np.cross(u, v), axis=1), np.sum(u*v, axis=1))
        np.add.at(angles, triangles[:, k], theta)
    return 2*np.pi - angles


def signed_volume(vertices, triangles):
    p = np.asarray(vertices)[np.asarray(triangles)]
    return float(np.sum(p[:, 0] * np.cross(p[:, 1], p[:, 2])) / 6)


def triangle_mesh_report(vertices, triangles):
    """Report geometric welding effects separately from abstract topology."""
    triangles = np.asarray(triangles, dtype=int)
    edge = np.sort(np.vstack((triangles[:, [0,1]], triangles[:, [1,2]], triangles[:, [2,0]])), axis=1)
    _, counts = np.unique(edge, axis=0, return_counts=True)
    areas = np.linalg.norm(np.cross(vertices[triangles[:,1]]-vertices[triangles[:,0]], vertices[triangles[:,2]]-vertices[triangles[:,0]]), axis=1) / 2
    return {
        "vertices": len(vertices), "triangles": len(triangles),
        "boundary_edges": int(np.count_nonzero(counts == 1)),
        "nonmanifold_edges": int(np.count_nonzero(counts > 2)),
        "edge_incidence": {str(k): int(v) for k,v in sorted(Counter(counts).items())},
        "degenerate_triangles": int(np.count_nonzero(areas <= 1e-16)),
        "duplicate_triangles": len(triangles) - len(np.unique(np.sort(triangles, axis=1), axis=0)),
        "signed_volume": signed_volume(vertices, triangles),
    }


def prepare_model(model):
    """Validate paper geometry, then triangulate and orient for export."""
    vertices = np.asarray(model["vertices"], dtype=float)
    faces, orientable, components = orient_faces(model["faces"])
    if not np.all(np.isfinite(vertices)):
        raise ValueError("Nonfinite coordinate")
    if len(vertices) != model["expected_vertices"]:
        raise ValueError("Wrong paper vertex count")
    if components != 1 or orientable != model["orientable"]:
        raise ValueError("Wrong connectedness or orientability")
    if len(np.unique(vertices, axis=0)) != len(vertices):
        raise ValueError("Distinct paper vertices have identical coordinates")
    check_vertex_links(len(vertices), faces)
    edges = edge_incidence(faces)
    chi = len(vertices) - len(edges) + len(faces)
    if chi != 2 - (2 if orientable else 1)*model["genus"]:
        raise ValueError(f"Wrong Euler characteristic: {chi}")
    angles = np.zeros(len(vertices))
    planarity = 0.0
    triangles = []
    for face in faces:
        p = vertices[face]
        _, normal, err = face_plane(p)
        planarity = max(planarity, err)
        for k, idx in enumerate(face):
            u, w = p[k]-p[k-1], p[(k+1) % len(p)]-p[k]
            angles[idx] += np.pi - np.arctan2(np.dot(normal, np.cross(u, w)), np.dot(u, w))
        triangles.extend(triangulate_face(vertices, face))
    triangles = np.asarray(triangles, dtype=int)
    defect = 2*np.pi - angles
    target = model["expected_defect"]
    err = float(np.max(np.abs(defect-target)))
    if err > 2e-8:
        raise ValueError(f"Nonconstant angular defect; max error = {err:g}")
    triangle_error = float(np.max(np.abs(triangle_defects(vertices, triangles)-defect)))
    if triangle_error > 1e-9:
        raise ValueError(f"Triangulation changed angular defects by {triangle_error:g}")
    if orientable and signed_volume(vertices, triangles) < 0:
        faces = [f[::-1] for f in faces]
        triangles = triangles[:, ::-1]
    stats = triangle_mesh_report(vertices, triangles)
    if stats["degenerate_triangles"]:
        raise ValueError("Degenerate triangulation")
    if model["name"].startswith("embedded") and (stats["boundary_edges"] or stats["nonmanifold_edges"] or stats["duplicate_triangles"] or stats["signed_volume"] <= 0):
        raise ValueError("Embedded model is not a closed oriented solid boundary")
    result = dict(model, vertices=vertices, faces=faces, triangles=triangles)
    report = {
        "name": model["name"], "paper_vertices": len(vertices),
        "polygon_edges": len(edges), "polygon_faces": len(faces),
        "euler_characteristic": chi, "genus": model["genus"],
        "orientable": orientable, "connected_components": components,
        "closed_abstract_surface": True, "vertex_links_are_cycles": True,
        "expected_defect_radians": target, "max_defect_error": err,
        "max_planarity_error": planarity, "triangulation_defect_error": triangle_error,
        "triangulation": stats,
    }
    return result, report


def write_stl(path, vertices, triangles, name):
    """Binary float32 STL, sharing the same rounded coordinates on every face."""
    vectors = np.asarray(vertices, dtype=np.float32)[np.asarray(triangles)]
    normals = np.cross(vectors[:,1]-vectors[:,0], vectors[:,2]-vectors[:,0])
    lengths = np.linalg.norm(normals, axis=1)
    if np.any(lengths == 0):
        raise ValueError("STL float32 rounding collapsed a triangle")
    normals /= lengths[:,None]
    dtype = np.dtype([("normal", "<f4", (3,)), ("vertices", "<f4", (3,3)), ("attribute", "<u2")])
    records = np.zeros(len(triangles), dtype=dtype)
    records["normal"] = normals
    records["vertices"] = vectors
    data = name.encode("ascii")[:80].ljust(80, b" ") + struct.pack("<I", len(records)) + records.tobytes()
    write_bytes_atomic(path, data)


def stl_surface(model):
    """Encode parallel abstract diagonals with distinct subdivision vertices.

    Only needed for a few singular polygon presentations. Inserted vertices
    are flat; they are not counted as vertices of the paper construction.
    The original polygons, corner coordinates, and their mapped images stay
    unchanged. Intersections between unrelated sheets are still present.
    """
    vertices, triangles = model["vertices"], model["triangles"]
    if not triangle_mesh_report(vertices, triangles)["nonmanifold_edges"]:
        return vertices, triangles
    result_vertices = list(vertices)
    result_triangles = []
    used = {tuple(np.float32(p)) for p in vertices}

    def unique_point(points):
        for attempt in range(1000):
            t = 0.371 + 0.000413 * attempt
            weights = np.array([1.0, t] if len(points) == 2 else [1.0, t, t*t])
            point = weights @ np.asarray(points) / weights.sum()
            key = tuple(np.float32(point))
            if key not in used:
                used.add(key)
                result_vertices.append(point)
                return len(result_vertices)-1
        raise ValueError("Cannot separate coincident subdivision points")

    for face in model["faces"]:
        boundary = {tuple(sorted((a,b))) for a,b in zip(face, face[1:]+face[:1])}
        midpoints = {}
        for triangle in triangulate_face(vertices, face):
            ring = []
            for a,b in zip(triangle, triangle[1:]+triangle[:1]):
                ring.append(a)
                edge = tuple(sorted((a,b)))
                if edge not in boundary:
                    if edge not in midpoints:
                        midpoints[edge] = unique_point(vertices[list(edge)])
                    ring.append(midpoints[edge])
            center = unique_point(vertices[triangle])
            result_triangles.extend([[a,b,center] for a,b in zip(ring, ring[1:]+ring[:1])])
    result_vertices = np.asarray(result_vertices)
    result_triangles = np.asarray(result_triangles, dtype=int)
    report = triangle_mesh_report(result_vertices, result_triangles)
    if report["boundary_edges"] or report["nonmanifold_edges"] or report["duplicate_triangles"]:
        raise ValueError("Refined STL surface is not closed")
    defects = triangle_defects(result_vertices, result_triangles)
    if np.max(np.abs(defects[:len(vertices)]-model["expected_defect"])) > 2e-8 or np.max(np.abs(defects[len(vertices):])) > 2e-8:
        raise ValueError("Export refinement changed intrinsic curvature")
    check_vertex_links(len(result_vertices), result_triangles.tolist())
    return result_vertices, result_triangles


def read_stl(path):
    data = Path(path).read_bytes()
    count = struct.unpack_from("<I", data, 80)[0] if len(data) >= 84 else -1
    if len(data) == 84 + 50*count:
        dtype = np.dtype([("normal", "<f4", (3,)), ("vertices", "<f4", (3,3)), ("attribute", "<u2")])
        points = np.frombuffer(data, dtype=dtype, count=count, offset=84)["vertices"].reshape(-1,3).astype(float)
    else:
        points = np.array([[float(n) for n in line.split()[1:]] for line in data.decode("ascii").splitlines() if line.strip().startswith("vertex ")])
    vertices, inverse = np.unique(points, axis=0, return_inverse=True)
    return vertices, inverse.reshape(-1,3)


def write_obj(path, model):
    """Preserve polygon incidences and double-precision paper coordinates."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"# {model['name']}", f"# {model['construction']}", f"# {model['source']}"]
    lines += ["v " + " ".join(format(x, ".17g") for x in p) for p in model["vertices"]]
    lines += ["f " + " ".join(str(i+1) for i in f) for f in model["faces"]]
    write_bytes_atomic(path, ("\n".join(lines) + "\n").encode("utf-8"))
