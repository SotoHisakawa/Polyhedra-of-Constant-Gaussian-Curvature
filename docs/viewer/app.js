(() => {
  "use strict";
  const byId = id => document.getElementById(id);
  const collection = new Map(MODEL_LIBRARY.map(model => [model.name, model]));
  const routes = Constructions.list(collection);
  const routeById = new Map(routes.map(route => [route.id, route]));
  const titles = {
    tetrahedron: "Regular tetrahedron", flat_torus: "Flat torus",
    windows: "Prism with windows", shared_windows: "Prism with shared windows",
    drilled_windows: "Drilled window surface", central_handle: "Central handle",
    historical_cube: "Cube construction (historical)", mobius_bands: "Möbius bands on a convex core",
    projective_plane: "Projective plane", klein_bottle: "Klein bottle",
    antiprism: "Antiprism with bands", uniform_core: "Uniform core with bands",
    tetrahemihexahedron: "Tetrahemihexahedron", cubohemioctahedron: "Cubohemioctahedron",
    historical_cuboctahedron: "Cuboctahedron with bands (historical)",
    cubohemioctahedron_hexagonal_drilling: "Hexagonally drilled cubohemioctahedron",
    small_rhombihexahedron: "Small rhombihexahedron",
    small_rhombihexahedron_square_drilling: "Square-drilled small rhombihexahedron",
    small_dodecahemidodecahedron: "Small dodecahemidodecahedron"
  };
  const modelTitle = model => titles[model.name.replace(/^(embedded|immersed|singular)_(nonorientable|orientable)_g\d+_v\d+_/, "")] || model.construction;
  const category = model => model.name.split("_")[0];
  const capital = text => text[0].toUpperCase() + text.slice(1);
  const minus = number => String(number).replace("-", "−");
  const text = (id, value) => { byId(id).textContent = value; };
  function options(select, entries) {
    select.replaceChildren(...entries.map(([value, title]) => new Option(title, value)));
  }

  for (const type of ["embedded", "immersed", "singular"]) {
    const group = document.createElement("optgroup");
    group.label = capital(type);
    const models = MODEL_LIBRARY.filter(model => category(model) === type);
    models.sort((a, b) => Number(b.orientable) - Number(a.orientable) || a.genus - b.genus || a.vertices.length - b.vertices.length);
    for (const model of models) {
      group.append(new Option(modelTitle(model) + " · g = " + model.genus + " · " + model.vertices.length + " vertices", model.name));
    }
    byId("model").append(group);
  }
  options(byId("route"), routes.map(route => [route.id, route.title]));
  byId("model").value = "embedded_orientable_g3_v24_windows";

  let mode = "general", index = 0, current, edges = [], viewer;
  function showRenderError() {
    text("render-error", "The 3D view is unavailable in this browser. Measurements and surface data are available below or alongside the view.");
    byId("render-error").hidden = false;
    for (const id of ["wire", "zoom-in", "zoom-out", "reset"]) byId(id).disabled = true;
  }
  try {
    viewer = new MeshViewer(byId("mesh"), {onPick: picked});
  } catch (error) {
    console.error(error);
    showRenderError();
  }

  function picked(selection) {
    if (!selection) return;
    if (selection.kind === "edge") {
      const id = [selection.a, selection.b].sort((a, b) => a - b).join(",");
      byId("edge").value = id;
    } else if (selection.kind === "angle") {
      byId("face").value = String(selection.face);
      setCornerOptions();
      byId("corner").value = String(selection.corner);
    }
    readMeasurement();
  }
  function setCornerOptions() {
    const face = current.faces[Number(byId("face").value)];
    options(byId("corner"), face.map((v, i) => [i, "v" + (v + 1)]));
  }
  function readMeasurement() {
    const kind = byId("measure-kind").value;
    const angle = kind === "angle";
    byId("edge-box").hidden = angle;
    byId("angle-box").hidden = !angle;
    byId("length-units").hidden = angle;
    byId("angle-unit").hidden = !angle;
    const face = Number(byId("face").value), corner = Number(byId("corner").value);
    const [a, b] = byId("edge").value.split(",").map(Number);
    viewer?.setPickMode(kind);
    viewer?.setSelection(angle ? {kind, face, corner} : {kind, a, b});
    let value;
    if (angle) {
      value = MeshViewer.cornerAngle(current, face, corner);
      const radians = byId("angle-unit").value === "radians";
      text("value", radians ? (value * Math.PI / 180).toFixed(5) : value.toFixed(2) + "°");
    } else {
      value = MeshViewer.edgeLength(current, a, b);
      text("value", value.toFixed(4));
    }
    const name = angle ? "F" + (face + 1) + " · v" + (current.faces[face][corner] + 1) : "v" + (a + 1) + "–v" + (b + 1);
    const unit = angle && byId("angle-unit").value === "radians" ? " rad" : "";
    viewer?.setMeasurementLabel(byId("annotate").checked ? name + "   " + byId("value").textContent + unit : null);
    text("measure-help", angle ? "Pick a corner inside a visible face, or choose its face and vertex above." : "Select an edge on the surface or from the list above.");
  }
  function showDefects(model, chi) {
    const values = Array(model.vertices.length).fill(360);
    model.faces.forEach((face, i) => face.forEach((vertex, j) => {
      values[vertex] -= MeshViewer.cornerAngle(model, i, j);
    }));
    const lo = Math.min(...values), hi = Math.max(...values);
    if (hi - lo > 1e-6) {
      text("defect-label", "Vertex defect range");
      text("defect", minus(Number(lo.toFixed(2))) + "° to " + minus(Number(hi.toFixed(2))) + "°");
      return;
    }
    text("defect-label", "Common vertex defect");
    text("defect", fractionOfPi(2 * chi, values.length));
  }

  let targetSequence = null;
  const planner = createTargetFinder({
    models: MODEL_LIBRARY, routes, modelTitle,
    openModel(model) {
      byId("model").value = model.name;
      visitMode("general");
    },
    openRoute(route, stages) {
      byId("route").value = route.id;
      index = 0;
      targetSequence = {...route, steps:stages};
      visitMode("construction");
    }
  });

  function render(preserveView = false) {
    const planning = mode === "target", following = mode === "construction";
    byId("find-target").setAttribute("aria-pressed", String(planning));
    byId("workbench").hidden = planning;
    byId("planner").hidden = !planning;
    byId("general").setAttribute("aria-pressed", String(mode === "general"));
    byId("construction").setAttribute("aria-pressed", String(following));
    if (planning) {
      document.title = "Find a CCP · Polyhedra of Constant Gaussian Curvature";
      planner.restore(location.hash);
      return;
    }
    byId("construction-card").hidden = !following;
    byId("model").hidden = following;
    byId("route").hidden = !following;
    const label = byId("choose-label");
    label.htmlFor = following ? "route" : "model";
    label.textContent = following ? "Construction sequence" : "Model";
    const route = targetSequence ?? routeById.get(byId("route").value);
    const lastStep = route.steps.length - 1;
    const step = following ? route.steps[index] : null;
    current = step ? step.model : collection.get(byId("model").value);
    text("title", following ? route.title : modelTitle(current));
    text("category", following ? "Construction mode" : "Model collection");
    document.title = byId("title").textContent + " · Polyhedra of Constant Gaussian Curvature";
    if (following) {
      text("step-count", "Step " + (index + 1) + " of " + (lastStep + 1));
      text("step-title", step.title);
      text("step-text", step.description);
      text("reference", step.reference ?? route.reference);
      const connections = byId("connections");
      connections.replaceChildren();
      for (const other of routes) {
        if (other.id === route.id) continue;
        const leadsHere = index === 0 && other.steps.at(-1).model.name === step.model.name;
        const continues = index === route.steps.length - 1 && other.steps[0].model.name === step.model.name;
        if (!leadsHere && !continues) continue;
        const button = document.createElement("button");
        button.className = "route-link";
        button.textContent = (leadsHere ? "Build the starting surface · " : "Continue to ") + other.title;
        button.onclick = () => { byId("route").value = other.id; index = 0; targetSequence = null; render(); byId("route").focus(); };
        connections.append(button);
      }
      byId("step-prev").disabled = index === 0;
      byId("step-next").disabled = index === lastStep;
      byId("steps").replaceChildren(...route.steps.map((item, i) => {
        const button = document.createElement("button");
        button.textContent = String(i + 1);
        button.setAttribute("aria-label", "Step " + (i + 1) + ". " + item.title);
        if (i === index) button.setAttribute("aria-current", "step");
        button.addEventListener("click", () => { index = i; render(true); byId("steps").children[i].focus(); });
        return button;
      }));
    }
    edges = MeshViewer.edges(current);
    const genus = current.genus;
    const chi = current.vertices.length - edges.length + current.faces.length;
    const type = following && !collection.has(current.name) ? "Construction stage" : capital(category(current));
    text("surface-type", type + " · " + (current.orientable ? "Orientable" : "Non-orientable (crosscap genus)"));
    text("genus", genus); text("chi", minus(chi));
    text("vertices", current.vertices.length); text("edges", edges.length); text("faces", current.faces.length);
    showDefects(current, chi);
    const highlights = {faces: step?.changedFaces || [], vertices: step?.changedVertices || []};
    byId("legend").hidden = highlights.faces.length + highlights.vertices.length === 0;
    viewer?.setModel(current, highlights, {preserveView});
    viewer?.setRevealChanges(following && !byId("legend").hidden && byId("reveal").checked);
    options(byId("edge"), edges.map(([a, b]) => [[a,b].join(","), "v" + (a + 1) + " — v" + (b + 1)]));
    options(byId("face"), current.faces.map((vs, i) => [i, "F" + (i + 1) + " · " + vs.length + " sides"]));
    if (highlights.faces.length) {
      const faceId = highlights.faces[0], face = current.faces[faceId];
      byId("face").value = String(faceId);
      byId("edge").value = [face[0], face[1]].sort((a,b) => a-b).join(",");
    }
    setCornerOptions();
    readMeasurement();
  }

  function setMode(next, refresh = false) {
    if (next !== mode || refresh) { mode = next; render(); }
  }
  function hashMode() {
    const name = location.hash.slice(1).split("?")[0];
    return ["target", "construction"].includes(name) ? name : "general";
  }
  function visitMode(next) {
    const fragment = next === "target" ? planner.fragment() : "#" + next;
    if (location.hash !== fragment) location.hash = fragment;
    setMode(next);
  }
  byId("general").onclick = () => visitMode("general");
  byId("construction").onclick = () => {
    if (targetSequence) { targetSequence = null; index = 0; if (mode === "construction") render(); }
    visitMode("construction");
  };
  byId("find-target").onclick = () => visitMode("target");
  window.addEventListener("hashchange", () => setMode(hashMode(), hashMode() === "target"));
  byId("route").onchange = () => { targetSequence = null; index = 0; render(); };
  byId("model").onchange = () => render();
  byId("step-prev").onclick = () => { if (index > 0) { index--; render(true); } };
  byId("step-next").onclick = () => {
    const route = targetSequence ?? routeById.get(byId("route").value);
    if (index < route.steps.length - 1) { index++; render(true); }
  };
  for (const id of ["measure-kind", "edge", "corner", "angle-unit"]) byId(id).onchange = readMeasurement;
  byId("face").onchange = () => { setCornerOptions(); readMeasurement(); };
  byId("zoom-in").onclick = () => viewer?.zoom(1.2);
  byId("zoom-out").onclick = () => viewer?.zoom(1/1.2);
  byId("reset").onclick = () => viewer?.reset();
  byId("reveal").onchange = () => viewer?.setRevealChanges(byId("reveal").checked);
  byId("annotate").onchange = readMeasurement;
  byId("wire").onchange = () => viewer?.setWireframe(byId("wire").checked);
  mode = hashMode();
  render();
})();
