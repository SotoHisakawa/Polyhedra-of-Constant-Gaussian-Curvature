// Step meshes are built from the same oriented polygon cycles as the model catalogue.
const Constructions = (() => {
  const models = {
    windows: "embedded_orientable_g2_v16_windows",
    drilled: "embedded_orientable_g3_v32_drilled_windows",
    torus: "embedded_orientable_g1_v9_flat_torus"
  };

  function intermediate(name, title, vertices, faces, genus, expectedDefect) {
    const model = {name, construction: title, vertices, faces, genus, orientable: true};
    if (expectedDefect !== undefined) model.expected_defect = expectedDefect;
    return model;
  }

  function step(title, description, model, changedFaces = [], changedVertices = []) {
    return {title, description, model, changedFaces, changedVertices};
  }

  function squarePrismStage(completed, windowModel) {
    // Each window sector lists two end pairs and four interior rim/waist vertices.
    // Undo an insertion by omitting its annulus and filling the aperture in each cap.
    const interiors = [[2, 3, 4, 5], [10, 11, 12, 13]];
    const omitted = new Set(interiors.slice(completed).flat());
    const vertices = [];
    const index = new Map();
    windowModel.vertices.forEach((point, id) => {
      if (omitted.has(id)) return;
      index.set(id, vertices.length);
      vertices.push(point);
    });
    const faces = [];
    const changedFaces = [];
    const added = new Set(completed ? interiors[completed - 1] : []);
    windowModel.faces.forEach((original, id) => {
      const face = id < 2 ? original.filter(v => !omitted.has(v)) : original;
      if (face.some(v => omitted.has(v))) return;
      if (id >= 2 && face.some(v => added.has(v))) changedFaces.push(faces.length);
      faces.push(face.map(v => index.get(v)));
    });
    const title = ["Square prism", "Single window", "Two-window surface"][completed];
    const defect = [Math.PI / 2, undefined, -Math.PI / 4][completed];
    const model = completed === 2 ? windowModel : intermediate(
      "square_prism_" + completed, title, vertices, faces, completed, defect
    );
    return {model, changedFaces, changedVertices: [...added].map(v => index.get(v))};
  }

  function windowPath(windowModel) {
    const descriptions = [
      "Start with a square of side 2 and take its prism of half-height 1. Every corner has exterior turn π/2. Mark opposite sides for two windows with a = π/8.",
      "Cut matching triangular apertures from the upper and lower bases. Insert W(π/8), a pinched annulus with two aperture apices and two waist vertices. One handle adds four vertices. The untouched corners still have defect π/2.",
      "Repeat on the opposite side. Every original prism corner now meets one window, giving δ = π/2 − (π − 2a) = −π/4. The apices and waist vertices have the same defect."
    ];
    const titles = ["Start with the square prism", "Insert the first window", "Insert the opposite window"];
    return {
      id: "windows", title: "Prism with windows", reference: "§3 · A prism with windows",
      steps: titles.map((title, i) => {
        const stage = squarePrismStage(i, windowModel);
        return step(title, descriptions[i], stage.model, stage.changedFaces, stage.changedVertices);
      })
    };
  }

  function drillingPath(windowModel, drilledModel) {
    // The generator appends the octagonal rims (top, bottom) and eight lateral faces.
    const n = drilledModel.parameters.tunnel_sides;
    const firstRim = windowModel.vertices.length;
    const upper = Array.from({length: n}, (_, i) => firstRim + i);
    const lower = upper.map(i => i + n);
    const rimVertices = [...upper, ...lower];
    const tubeStart = drilledModel.faces.length - n;
    const tubeFaces = Array.from({length: n}, (_, i) => tubeStart + i);
    // Replace the tube by two planar discs, keeping the surrounding cap partition.
    // Outward normals of the two discs are +z and −z respectively.
    const closedFaces = [...drilledModel.faces.slice(0, tubeStart), upper, lower.slice().reverse()];
    const capModel = intermediate("octagonal_mouths", "Marked octagonal discs",
      drilledModel.vertices, closedFaces, windowModel.genus);
    return {
      id: "drilling", title: "Octagonal drilling", reference: "§3 · Adding handles at fixed defect",
      steps: [
        step("Choose the drilling axis",
          "Begin with the genus-two, sixteen-vertex window surface. The vertical segment through the centres of its horizontal caps misses the rest of the surface. Its two endpoints lie between the windows.",
          windowModel, [0, 1]),
        step("Mark the two octagonal discs",
          "Subdivide each cap around a regular octagon of circumradius 0.2, perpendicular to the chosen axis. The sixteen new rim vertices still have defect zero in these planar faces. The genus remains two.",
          capModel, [tubeStart, tubeStart + 1], rimVertices),
        step("Replace the discs by an annulus",
          "Remove both octagonal interiors and insert the eight rectangles of the prism wall. The operation lowers χ by two. At each rim vertex the removed disc had turning 2π/8, and the wall contributes zero, giving δ = −π/4 throughout.",
          drilledModel, tubeFaces, rimVertices)
      ]
    };
  }

  function subdivideTriangle(base, faceIndex, title) {
    const triangle = base.faces[faceIndex];
    const centre = [0, 1, 2].map(axis =>
      triangle.reduce((sum, vertex) => sum + base.vertices[vertex][axis], 0) / 3
    );
    const vertex = base.vertices.length;
    const faces = base.faces.slice();
    const subfaces = triangle.map((v, i) => [v, triangle[(i + 1) % 3], vertex]);
    faces.splice(faceIndex, 1, ...subfaces);
    return {
      model: intermediate("subdivided_flat_torus_" + (vertex + 1), title,
        [...base.vertices, centre], faces, 1, 0),
      changedFaces: [faceIndex, faceIndex + 1, faceIndex + 2],
      changedVertices: [vertex]
    };
  }

  function vertexPath(torus) {
    // Pick faces on two outer strips to keep the added flat vertices visible.
    const firstFace = 12;
    const first = subdivideTriangle(torus, firstFace, "Ten-vertex flat torus");
    const second = subdivideTriangle(first.model, 7, "Eleven-vertex flat torus");
    return {
      id: "vertices", title: "Adding flat vertices",
      reference: "§3 · Varying the vertex count at fixed genus; Appendix B",
      steps: [
        step("Choose a triangular face",
          "Use the nine-vertex embedded flat torus of Appendix B with twist θ = π/4. Its faces are triangles and every vertex has angle sum 2π. Select a face on the lower short strip.",
          torus, [firstFace]),
        step("Insert a vertex in the face",
          "Place a vertex at the centroid and join it to the three corners in the face plane. The new angle sum is 2π, and each old corner angle is partitioned. Thus flatness is preserved, with ΔV = 1, ΔE = 3 and ΔF = 2.",
          first.model, first.changedFaces, first.changedVertices),
        step("Repeat at fixed genus",
          "Subdivide a face on the upper short strip in the same way. The torus now has eleven vertices and χ = 0. Repeating a triangular subdivision realises every larger vertex count on this same embedded surface.",
          second.model, second.changedFaces, second.changedVertices)
      ]
    };
  }


  // W*(a) from §3. H(u) first prescribes the defect at the extra waist vertex.
  // Bisect the remaining apex-angle equation, using the endpoint limits in the proof.
  function oddWindowParameters(a) {
    const t = Math.tan(2 * a);
    const halfTurn = Math.sin(a / 2);
    function height(u) {
      const leg = t - u;
      const factor = leg / (Math.hypot(1, leg) * halfTurn);
      return u * Math.sqrt(Math.max(0, factor * factor - 1));
    }
    function apexAngles(u) {
      const h = height(u);
      const angle1 = Math.atan2(h, Math.hypot(t, 1));
      // Angle between (t,1,−h) and (u,0,−h) at the upper apex.
      const angle2 = Math.atan2(Math.hypot(h, h * (t - u), u), t * u + h * h);
      return angle1 + angle2;
    }
    let low = 0;
    let high = t - Math.tan(a / 2);
    for (let iteration = 0; iteration < 55; iteration++) {
      const u = (low + high) / 2;
      if (apexAngles(u) > Math.PI / 2 - a) low = u;
      else high = u;
    }
    const u = (low + high) / 2;
    const H = height(u);
    const b = 1 + H * Math.cos(2 * a) * halfTurn / (Math.cos(2 * a) + Math.cos(a / 2));
    return {u, H, b, t};
  }

  function oddWindowPath(windowModel) {
    // Genus 2, K = 0, epsilon = 1. D = 2 keeps the two openings disjoint.
    const a = 2 * Math.PI / 17;
    const {u, H, b, t} = oddWindowParameters(a);
    const h = H / b;
    const leftHalf = 1 + 2 * Math.tan(a / 4);
    const leftPinch = h * Math.sin(a) / 2;
    const leftApex = -1 + leftHalf * Math.tan(a);
    const leftWaist = -1 + leftPinch * Math.tan(a);
    const rightApex = 1 - t;
    const rightWaist = rightApex + t / b;
    // Preserve vertex indices from the two-window catalogue model for all 16 old points.
    const vertices = [
      [1, -1, h], [1, -1, -h], [rightApex, 0, h], [rightApex, 0, -h],
      [rightWaist, -1/b, 0], [rightWaist, 1/b, 0], [1, 1, h], [1, 1, -h],
      [-1, leftHalf, h], [-1, leftHalf, -h], [leftApex, 0, h], [leftApex, 0, -h],
      [leftWaist, leftHalf-leftPinch, 0], [leftWaist, -leftHalf+leftPinch, 0],
      [-1, -leftHalf, h], [-1, -leftHalf, -h]
    ];
    const retuned = intermediate("retuned_even_window", "Retuned genus-two windows",
      vertices, windowModel.faces, 2);
    // Replace the two pentagons meeting at the right-hand inner edge by two fans.
    const waist = vertices.length;
    const sides = [
      [2, 0, 4], [2, 4, waist], [2, waist, 5], [2, 5, 6],
      [3, 4, 1], [3, waist, 4], [3, 5, waist], [3, 7, 5]
    ];
    const retained = windowModel.faces.filter((_, i) => i !== 2 && i !== 3);
    const modified = intermediate("odd_genus_two_window_17", "Seventeen-vertex genus-two surface",
      [...vertices, [rightApex + u/b, 0, 0]], [...retained, ...sides], 2, -2 * a);
    return {
      id: "odd-window", title: "Adding a waist vertex",
      reference: "§3 · Varying the vertex count at fixed genus; Closing the prism base",
      steps: [
        step("Begin with two ordinary windows",
          "The genus-two square-prism model has sixteen vertices of defect −π/4. To keep genus two with seventeen vertices, the Gauss–Bonnet theorem prescribes −4π/17 at every vertex.",
          windowModel, [2, 3, 4, 5]),
        step("Retune the openings and height",
          "Set a = 2π/17. The left aperture has base angle a and the right 2a. With base midpoints two units apart, use half-widths 1 + 2 tan(a/4) and 1. The waist widths and half-height are chosen for the odd window. At this intermediate stage the defects vary.",
          retuned, [2, 3, 4, 5]),
        step("Pinch the inner edge",
          "Replace the two pentagons meeting at the right inner edge by eight triangles through the two old waist corners and one new waist vertex. The waist- and apex-angle equations for W*(a) give common defect −2a = −4π/17. The number of handles is unchanged.",
          modified, sides.map((_, i) => retained.length + i), [waist])
      ]
    };
  }

  function list(modelMap) {
    const get = name => modelMap instanceof Map ? modelMap.get(name) : modelMap[name];
    const windowModel = get(models.windows);
    return [windowPath(windowModel), drillingPath(windowModel, get(models.drilled)),
      vertexPath(get(models.torus)), oddWindowPath(windowModel)];
  }

  return {list};
})();
