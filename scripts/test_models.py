"""Geometry regressions; run with ``python scripts/test_models.py``.

The checks exercise the source polygons, deliberately damaged inputs, and the
actual float32 STL bytes. Temporary exports never touch the model collection.
"""

from collections import Counter
from contextlib import redirect_stdout
from copy import deepcopy
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np

from embedded_models import (build_historical_cube, build_models as embedded_models,
                             build_tetrahedron)
from immersed_models import build_models as immersed_models
from singular_models import build_models as singular_models
from mesh_geometry import (check_vertex_links, edge_incidence, orient_faces,
                           prepare_model, read_stl, stl_surface, triangle_defects,
                           triangle_mesh_report, triangulate_face, write_stl)
from generate_models import main as generate, verify_export


class ModelGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = embedded_models() + immersed_models() + singular_models()
        cls.prepared = {model["name"]: prepare_model(model) for model in cls.raw}

    def assert_closed(self, vertices, triangles):
        report = triangle_mesh_report(vertices, triangles)
        for key in ("boundary_edges", "nonmanifold_edges", "degenerate_triangles",
                    "duplicate_triangles"):
            self.assertEqual(report[key], 0, (key, report))
        check_vertex_links(len(vertices), triangles.tolist())
        return report

    def test_all_source_surfaces_have_the_claimed_topology_and_curvature(self):
        self.assertEqual(len(self.prepared), len(self.raw), "Model names must be unique")
        self.assertTrue(any(m["name"].startswith("embedded_") for m in self.raw))
        self.assertTrue(any(m["name"].startswith("immersed_") for m in self.raw))
        self.assertTrue(any(m["name"].startswith("singular_") for m in self.raw))
        for raw in self.raw:
            with self.subTest(model=raw["name"]):
                model, report = self.prepared[raw["name"]]
                vertices, faces = model["vertices"], model["faces"]
                edges = edge_incidence(faces)
                self.assertTrue(all(len(incident) == 2 for incident in edges.values()))
                _, orientable, components = orient_faces(faces)
                self.assertEqual(components, 1)
                self.assertEqual(orientable, raw["orientable"])
                chi = len(vertices) - len(edges) + len(faces)
                self.assertEqual(chi, 2 - (2 if orientable else 1) * raw["genus"])
                # Gauss--Bonnet supplies the target independently of the model's
                # expected_defect field; the triangles independently sum angles.
                target = 2 * np.pi * chi / len(vertices)
                np.testing.assert_allclose(
                    triangle_defects(vertices, model["triangles"]), target,
                    atol=2e-8, rtol=0)
                self.assertAlmostEqual(raw["expected_defect"], target, places=12)
                self.assertEqual(len(vertices), raw["expected_vertices"])
                self.assertTrue(report["closed_abstract_surface"])
                self.assertLess(report["max_planarity_error"], 1e-9)

    def test_every_stl_round_trip_is_closed_and_preserves_its_surface(self):
        with tempfile.TemporaryDirectory(prefix="ccp_stl_tests_") as directory:
            for name, (model, _) in self.prepared.items():
                with self.subTest(model=name):
                    vertices, triangles = stl_surface(model)
                    path = Path(directory) / (name + ".stl")
                    write_stl(path, vertices, triangles, name)
                    loaded_vertices, loaded_triangles = read_stl(path)
                    np.testing.assert_array_equal(
                        loaded_vertices[loaded_triangles],
                        vertices[triangles].astype(np.float32).astype(float))
                    self.assert_closed(loaded_vertices, loaded_triangles)
                    _, orientable, components = orient_faces(loaded_triangles.tolist())
                    self.assertEqual(orientable, model["orientable"])
                    self.assertEqual(components, 1)
                    edges = edge_incidence(loaded_triangles.tolist())
                    chi = len(loaded_vertices) - len(edges) + len(loaded_triangles)
                    self.assertEqual(chi, 2 - (2 if orientable else 1) * model["genus"])
                    verify_export(path, vertices, triangles)

    def test_parallel_abstract_diagonals_get_flat_distinct_subdivision_points(self):
        # A naive triangle STL welds these different abstract diagonals. The
        # result must remain the original surface, not be repaired by filling it.
        examples = ["singular_orientable_g1_v6_flat_torus",
                    "singular_nonorientable_g4_v12_cubohemioctahedron"]
        for name in examples:
            with self.subTest(model=name):
                model, report = self.prepared[name]
                self.assertGreater(report["triangulation"]["nonmanifold_edges"], 0)
                vertices, triangles = stl_surface(model)
                source_count = len(model["vertices"])
                self.assertGreater(len(vertices), source_count)
                np.testing.assert_array_equal(vertices[:source_count], model["vertices"])
                self.assertEqual(len(np.unique(vertices.astype(np.float32), axis=0)),
                                 len(vertices))
                self.assert_closed(vertices, triangles)
                defects = triangle_defects(vertices, triangles)
                np.testing.assert_allclose(defects[:source_count],
                                           model["expected_defect"], atol=2e-8, rtol=0)
                np.testing.assert_allclose(defects[source_count:], 0, atol=2e-8, rtol=0)

    def test_concave_polygon_preserves_reflex_angle_and_area(self):
        # The L shape starts outside its kernel, so a fan at vertex zero is
        # invalid. Its area is 3*1 + 1*3 - 1*1 = 5; corner 2 is reflex.
        vertices = np.array([[3, 0, 0], [3, 1, 0], [1, 1, 0],
                             [1, 3, 0], [0, 3, 0], [0, 0, 0]], dtype=float)
        triangles = np.asarray(triangulate_face(vertices, list(range(6))))
        self.assertEqual(len(triangles), 4)
        points = vertices[triangles]
        oriented_areas = np.cross(points[:, 1] - points[:, 0],
                                  points[:, 2] - points[:, 0])[:, 2] / 2
        self.assertTrue(np.all(oriented_areas > 0))
        self.assertAlmostEqual(oriented_areas.sum(), 5.0)
        expected_angles = np.full(6, np.pi / 2)
        expected_angles[2] = 3 * np.pi / 2
        np.testing.assert_allclose(2 * np.pi - triangle_defects(vertices, triangles),
                                   expected_angles, atol=1e-12, rtol=0)
        # Every original boundary survives once; every inserted diagonal twice.
        counts = Counter(tuple(sorted(edge)) for triangle in triangles
                         for edge in zip(triangle, np.roll(triangle, -1)))
        boundary = {tuple(sorted((i, (i + 1) % 6))) for i in range(6)}
        self.assertEqual({edge for edge, count in counts.items() if count == 1}, boundary)
        self.assertTrue(all(count in (1, 2) for count in counts.values()))

    def test_crossed_polygon_with_nonzero_signed_area_is_rejected(self):
        vertices = np.array([[0, 0, 0], [3, 2, 0], [0, 2, 0], [2, 0, 0]], dtype=float)
        with self.assertRaisesRegex(ValueError, "Nonsimple polygon"):
            triangulate_face(vertices, [0, 1, 2, 3])

    def test_missing_face_is_rejected(self):
        broken = deepcopy(build_tetrahedron())
        broken["faces"].pop()
        with self.assertRaisesRegex(ValueError, "incident faces"):
            prepare_model(broken)

    def test_coordinate_perturbation_cannot_pass_constant_defect_validation(self):
        broken = deepcopy(build_tetrahedron())
        broken["vertices"][0, 0] += 0.05
        # All its faces remain planar triangles and its topology is unchanged.
        with self.assertRaisesRegex(ValueError, "Nonconstant angular defect"):
            prepare_model(broken)

    def test_incorrect_genus_or_orientability_is_rejected(self):
        for key, value, diagnostic in [("genus", 1, "Euler characteristic"),
                                       ("orientable", False, "orientability")]:
            with self.subTest(metadata=key):
                broken = deepcopy(build_tetrahedron())
                broken[key] = value
                with self.assertRaisesRegex(ValueError, diagnostic):
                    prepare_model(broken)

    def test_historical_cube_has_corrected_pi_over_six_defects(self):
        model, _ = prepare_model(build_historical_cube())
        self.assertEqual(len(model["vertices"]), 24)
        self.assertEqual(model["genus"], 2)
        self.assertTrue(model["historical"])
        np.testing.assert_allclose(
            triangle_defects(model["vertices"], model["triangles"]),
            -np.pi / 6, atol=1e-12, rtol=0)
        # The old STL used c=1/6, beyond the valid uncollapsed range.
        with self.assertRaises(ValueError):
            build_historical_cube(b=0.5, c=1 / 6)

    def test_check_mode_does_not_rewrite_files_and_detects_changed_geometry(self):
        name = "embedded_orientable_g0_v4_tetrahedron"
        with tempfile.TemporaryDirectory(prefix="ccp_check_tests_") as directory:
            root = Path(directory)
            args = ["--output", str(root), "--model", name]
            with redirect_stdout(io.StringIO()):
                generate(args)

            def snapshot():
                return {str(path.relative_to(root)): (path.read_bytes(), path.stat().st_mtime_ns)
                        for path in root.rglob("*") if path.is_file()}

            before = snapshot()
            with redirect_stdout(io.StringIO()):
                generate(args + ["--check"])
            self.assertEqual(snapshot(), before)
            path = root / "stl" / (name + ".stl")
            vertices, triangles = read_stl(path)
            vertices[0, 0] += 0.05
            write_stl(path, vertices, triangles, name)
            with redirect_stdout(io.StringIO()), self.assertRaisesRegex(
                    ValueError, "STL differs from regenerated geometry"):
                generate(args + ["--check"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
