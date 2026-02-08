import numpy as np
from stl.mesh import Mesh

def rotate_z(vertex, angle):
    """Rotate a vertex around the z-axis by given angle (in radians)"""
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    rotation_matrix = np.array([
        [cos_a, -sin_a, 0],
        [sin_a, cos_a, 0],
        [0, 0, 1]
    ])
    return rotation_matrix @ vertex

# Define the three base vertices
v1 = np.array([1, 0, 0])
v2 = np.array([np.cos(np.pi/8), np.sin(np.pi/8), 0.5])
v3 = np.array([np.cos(np.pi/4), np.sin(np.pi/4), 1])

# Generate all 9 vertices through rotational symmetry
vertices = []

# First set (0 rotation)
vertices.append(v1)      # v_1 (index 0)
vertices.append(v2)      # v_2 (index 1)
vertices.append(v3)      # v_3 (index 2)

# Second set (2π/3 rotation)
vertices.append(rotate_z(v1, 2*np.pi/3))  # v_{1,1} (index 3)
vertices.append(rotate_z(v2, 2*np.pi/3))  # v_{2,1} (index 4)
vertices.append(rotate_z(v3, 2*np.pi/3))  # v_{3,1} (index 5)

# Third set (4π/3 rotation)
vertices.append(rotate_z(v1, 4*np.pi/3))  # v_{1,2} (index 6)
vertices.append(rotate_z(v2, 4*np.pi/3))  # v_{2,2} (index 7)
vertices.append(rotate_z(v3, 4*np.pi/3))  # v_{3,2} (index 8)

vertices = np.array(vertices)

# Print vertices for verification
print("Vertices:")
for i, v in enumerate(vertices):
    print(f"v_{i}: ({v[0]:.6f}, {v[1]:.6f}, {v[2]:.6f})")

# Define faces (triangles) of the polyhedron
# Based on the flat torus construction, we need to connect vertices appropriately
# The polyhedron should have faces that respect the three-fold symmetry

faces = []

# Indexing: v_1=0, v_2=1, v_3=2, v_{1,1}=3, v_{2,1}=4, v_{3,1}=5, v_{1,2}=6, v_{2,2}=7, v_{3,2}=8

# For each of the 3 columns in the diagram
for i in range(3):
    # Current column indices
    v1_curr = i          # v_1, v_{1,1}, v_{1,2}
    v2_curr = i + 3      # v_2, v_{2,1}, v_{2,2}
    v3_curr = i + 6      # v_3, v_{3,1}, v_{3,2}

    # Next column indices (with wrapping)
    v1_next = (i + 1) % 3
    v2_next = ((i + 1) % 3) + 3
    v3_next = ((i + 1) % 3) + 6

    # Bottom row (white) - two triangles in each square
    # Lower triangle: v_1(curr), v_{1,1}(next), v_2(curr)
    faces.append([v1_curr, v1_next, v2_curr])
    # Upper triangle: v_2(curr), v_{1,1}(next), v_{2,1}(next)
    faces.append([v2_curr, v1_next, v2_next])

    # Middle row (cyan) - two triangles in each square
    # Lower triangle: v_2(curr), v_{2,1}(next), v_3(curr)
    faces.append([v2_curr, v2_next, v3_curr])
    # Upper triangle: v_3(curr), v_{2,1}(next), v_{3,1}(next)
    faces.append([v3_curr, v2_next, v3_next])

    # Top row (green) - two triangles in each square
    # Lower triangle: v_3(curr), v_{3,1}(next), v_1(curr) [wraps to bottom]
    faces.append([v3_curr, v3_next, v1_curr])
    # Upper triangle: v_1(curr), v_{3,1}(next), v_{1,1}(next) [wraps to bottom]
    faces.append([v1_curr, v3_next, v1_next])

# Convert faces to numpy array
faces = np.array(faces)

print(f"\nNumber of faces: {len(faces)}")

# Create the mesh
polyhedron_mesh = Mesh(np.zeros(faces.shape[0], dtype=Mesh.dtype))

for i, face in enumerate(faces):
    for j in range(3):
        polyhedron_mesh.vectors[i][j] = vertices[face[j]]

# Save to STL file
output_file = 'flat_torus.stl'
polyhedron_mesh.save(output_file)

print(f"\nSTL file saved to: {output_file}")

# Print edge lengths for verification
print("\nEdge length verification:")
print(f"|v_1 v_3| = {np.linalg.norm(vertices[0] - vertices[2]):.6f}, expected: {np.sqrt(3 - np.sqrt(2)):.6f}")
print(f"|v_1 v_{{3,1}}| = {np.linalg.norm(vertices[0] - vertices[5]):.6f}, expected: {np.sqrt(3 + 2*np.cos(np.pi/12)):.6f}")
print(f"|v_3 v_{{3,1}}| = {np.linalg.norm(vertices[2] - vertices[5]):.6f}, expected: {np.sqrt(3):.6f}")
print(f"|v_2 v_1| = {np.linalg.norm(vertices[1] - vertices[0]):.6f}, expected: {np.sqrt(9/4 - 2*np.cos(np.pi/8)):.6f}")
print(f"|v_2 v_{{1,1}}| = {np.linalg.norm(vertices[1] - vertices[3]):.6f}, expected: {np.sqrt(9/4 + 2*np.cos(11*np.pi/24)):.6f}")
print(f"|v_1 v_{{1,1}}| = {np.linalg.norm(vertices[0] - vertices[3]):.6f}, expected: {np.sqrt(3):.6f}")
