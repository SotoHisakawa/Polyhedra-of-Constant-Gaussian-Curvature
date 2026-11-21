const publishedOrientableCounts = new Map([[7,24], [19,60]]);
const nonorientableCoverBases = [[3,17], [4,24], [6,32], [18,24], [24,24]];

// Include the smallest count attained by each family that can improve the minimum.
// For g ≥ 3, the general chain constructions start at 8g (orientable) or 11g
// (non-orientable). Covers of uniform cores have at least 9g vertices.
function recipeVertexCounts(orientable, g) {
  if (orientable) {
    if (g <= 2) return [[4,6], [8,9], [16,17]][g];
    const counts = [6*(g-1), 8*g];
    if (publishedOrientableCounts.has(g)) counts.push(publishedOrientableCounts.get(g));
    return counts;
  }
  if (g === 1) return [9,10];
  if (g === 2) return [15,16];
  const counts = [9*g];
  if (g >= 5) counts.push(9*(g-2));
  if (g >= 8 && g % 2 === 0) counts.push(3*(g-2));
  for (const [baseGenus, baseVertices] of nonorientableCoverBases) {
    const sheets = (g-2)/(baseGenus-2);
    if (Number.isInteger(sheets)) counts.push(sheets*baseVertices);
  }
  return counts;
}

// Constructive cases for embedded orientable surfaces in Section 3.
function orientableRecipe(g, v) {
  if (!Number.isSafeInteger(g) || !Number.isSafeInteger(v) || g < 0 || v < 4) return null;
  const sections = "§3 · Varying the vertex count at fixed genus and Closing the prism base";
  const prismSection = "§3 · A prism with windows";
  const drillingSection = "§3 · Adding handles at fixed defect";
  const recipe = (title, steps, reference) => ({title, steps, reference});

  if (g === 0) {
    const n = v - 1;
    return recipe("A pyramid over a regular polygon", [
      `Place a regular ${n}-gon of circumradius one in a horizontal plane. Place the apex directly above its centre and join it to every base edge.`,
      `Set n = ${n} and v = ${v}. Choose the positive height h satisfying sin(π/n) / √(1+h²) = sin((π−2π/v)/n). Each of the n equal lateral angles at the apex is then (2π−4π/v)/n.`,
      `The apex has defect 4π/${v}. Rotation permutes all ${n} base vertices, and the Gauss–Bonnet theorem makes their common defect equal to the apex defect.`
    ], sections);
  }

  if (g === 1) {
    if (v < 8) return null;
    if (v === 8) return {
      title:"Schwartz’s eight-vertex flat torus",
      steps:[
        "Index eight vertices cyclically modulo 8. For each index j take the triangles (j, j+1, j+3) and (j, j+2, j+3), reducing labels modulo 8.",
        "Use the planar projections and paired heights in equation (1) of Schwartz’s construction. The half-turn about the vertical axis pairs vertex j with vertex 7−j.",
        "Choose the three heights near the stated approximations to satisfy the three angle-sum equations (2). Symmetry and the Gauss–Bonnet theorem give zero defect at all eight vertices. Schwartz proves that the resulting realisation is embedded."
      ],
      reference:"§3 · The torus vertex spectrum",
      source:{label:"Schwartz, Vertex-Minimal Paper Tori, §1.3",url:"https://arxiv.org/pdf/2507.14998#page=4"}
    };
    const steps = [
      "Place equilateral rings at heights 0, 1/2 and 1 on the unit cylinder. Rotate the second and third rings by π/8 and π/4 from the first.",
      "Join rings in the cyclic order 1, 3, 2, 1 by the six triangles per annulus specified in Appendix B. The resulting nine-vertex embedded torus is flat."
    ];
    if (v > 9) steps.push(`Perform ${v - 9} triangular face subdivision${v === 10 ? "" : "s"}. Each time insert a point in the relative interior of one triangle and join it to its corners. The new vertex is flat and all old angle sums are preserved. The final count is ${v}.`);
    return recipe("A flat torus with prescribed vertex count", steps, "§3 · Torus vertex counts · Appendix B · A flat torus with nine vertices");
  }

  if (publishedOrientableCounts.get(g) === v) {
    const group = g === 7 ? "octahedral" : "icosahedral";
    return {
      title:"Grünbaum and Shephard’s vertex-transitive polyhedron",
      steps:[
        `Use the genus-${g} embedded model with ${group} rotational symmetry in Grünbaum and Shephard’s Figure 4. Take the surface with its published face incidences.`,
        `Its vertices form a single orbit of size ${v}. Their defects are equal by symmetry, and the Gauss–Bonnet theorem gives the common value ${g === 7 ? "−π" : "−6π/5"}.`
      ],
      reference:"§5 · Published starting model",
      source:{label:"Grünbaum and Shephard, Figure 4",url:"https://mathreports.ca/article/polyhedra-with-transitivity-properties/"}
    };
  }

  if (g >= 3 && v === 6 * (g - 1)) {
    const steps = [
      "Start with Leopold’s twelve-vertex embedded genus-three polyhedron (Table 4 of Vertex-Transitive Polyhedra of Higher Genus, I). Its common defect is −2π/3."
    ];
    if (g > 3) steps.push(
      "Choose a rotation of order three. The segment joining the centroids of its two invariant triangular faces is perpendicular to those faces and meets the surface only at its ends.",
      `Drill ${g - 3} disjoint, sufficiently narrow triangular tunnels parallel to that segment. Replace each pair of triangular discs by the three prism walls. Each tunnel adds a handle and six vertices of defect −2π/3, giving ${v} vertices.`
    );
    else steps.push("Use all 32 triangular faces of its tetrahedrally symmetric realisation. Each of the twelve vertices has degree eight, and the genus already equals three.");
    return recipe(g === 3 ? "Leopold’s genus-three polyhedron" : "Triangular drilling at fixed defect", steps, drillingSection);
  }

  if (v === 8 * g || (g >= 5 && v === 6 * g)) {
    const shared = v === 6 * g;
    const n = shared ? g : 2 * g;
    return recipe(shared ? "Windows on all prism sides" : "Windows on alternate prism sides", [
      `Build the right prism of half-height one over a regular ${n}-gon of side two. Select ${shared ? "all" : "alternate"} sides, giving ${g} window positions.`,
      `Set a = π(${g}−1)/(${n}+2·${g}) and s = sin(a)/2. On both bases cut an isosceles aperture of base angle a at each selected side. Insert the window W(a), narrowed by the fraction s at mid-height.`,
      `Each window replaces two triangular discs by an annulus and adds four vertices. With ${g} handles the final count is 2·${n} + 4·${g} = ${v}, and all defects are −2a.`
    ], prismSection);
  }

  if (g >= 4 && g <= 7 && v === 7 * (g - 1)) {
    const n = g - 1;
    return recipe("Windows and a central handle", [
      `Take a sufficiently thin right prism over a regular ${n}-gon of side two. Set a = π/${n}, b = a + π/14 and ρ = cot(a).`,
      `Insert a window on every side with aperture depth cot(b) and base-corner turning 2π/7. The ${n} apices form matching regular polygons of radius q = ρ − cot(b) on the two bases.`,
      "Join the central polygons by a tunnel with an inward middle ring. For prism half-height h, move the middle ring inward by h sin(π/7)/√(sin²(a)−sin²(π/7)). Take h small enough that this displacement is less than q.",
      `The extra handle adds ${n} middle-ring vertices to the ${6 * n} window vertices and prism corners, giving ${v} vertices with common defect −4π/7.`
    ], "Companion collection · Central-handle face data");
  }

  if (g === 2 && v >= 16) {
    const odd = v % 2, k = (v - 16 - odd) / 2;
    return recipe("Two windows joined by prescribed chains", [
      `Set a = 2π/${v}. Choose two window sides in a planar base. Join their ends by an upper chain with ${Math.floor(k / 2)} interior vertices and a lower chain with ${Math.ceil(k / 2)}. Each chain vertex turns clockwise by 2a.`,
      `Use endpoint turns π−4a for an ordinary window${odd ? " and π−9a/2 for the modified right window" : " on either side"}. Use the two-window closure in §3, choosing the horizontal separation of the side midpoints slightly greater than D₀ so the aperture triangles are disjoint.`,
      odd ? "Raise the right prism to the half-height prescribed for W*(a). Insert an ordinary W(a) on the left and the three-waist-vertex W*(a) on the right." : "Raise a sufficiently thin right prism over this base and insert the ordinary window W(a) on each of the two sides.",
      `The base has ${4+k} corners. The prism contributes ${8+2*k} vertices and the two windows contribute ${8+odd}, giving ${v} vertices of defect −2a.`
    ], sections);
  }

  const tunnelSides = v / (2 * (g - 1));
  if (Number.isInteger(tunnelSides) && tunnelSides >= 8) {
    return recipe("Two windows followed by polygonal drilling", [
      `Use the genus-two, ${2*tunnelSides}-vertex two-window construction with n = ${tunnelSides}. The window angle is π/n and the total number of interior chain vertices is n−8 = ${tunnelSides-8}.`,
      `Choose ${g-2} disjoint narrow regular ${tunnelSides}-gonal prisms perpendicular to the base faces in the open strip between the two windows.`,
      `Remove both end discs of each prism and insert its lateral wall. Each drilling adds one handle and ${2*tunnelSides} vertices of defect −2π/${tunnelSides}. This gives genus ${g} and ${v} vertices.`
    ], "§3 · Closing the prism base and the orientable vertex-count progression");
  }

  if (v < 8 * g) return null;
  const odd = v % 2, total = (v - 8*g - odd) / 2;
  const k = Math.floor(total/g), extra = total % g;
  const a = 2*Math.PI*((g-1)/v), ordinary = Math.PI-4*a;
  const leastTurn = ordinary-odd*a/2;
  const phiMin = 2*(leastTurn - (k + Number(extra > 0))*a);
  const phiMax = 2*(ordinary - k*a);
  // At each corner, phi=E_i+E_(i+1)-2*k_i*a, and gamma=E_i-k_i*a.
  // These inequalities give 0<gamma<phi<pi and put each aperture strictly
  // inside the fan sector over its side of a tangential convex g-gon.
  // The guard deliberately leaves numerically marginal cases unclassified.
  const margin = 1e-9;
  if (phiMin <= margin || phiMax + 2*(1+odd)*a >= Math.PI - margin) return null;
  if (odd) {
    // Before trimming, a side length is tan(phi_prev/2)+tan(phi_next/2).
    // Rescale the special side to length 2. The lower bound for every other
    // side must exceed h*sin(a), using the odd-window lemma's height estimate.
    if (a >= Math.PI/4 - margin) return null;
    const minLength = 2*Math.tan(phiMin/2)/Math.tan(phiMax/2);
    const required = Math.tan(2*a)**2 * Math.cos(a/2)/2;
    if (minLength <= required + margin) return null;
  }
  const split = extra ? `${extra} chains with ${k+1} turns and ${g-extra} with ${k}` : `${g} chains with ${k} turns each`;
  return recipe("Windows on a base with short corner chains", [
    `Set a = 2π(${g}−1)/${v}. Distribute ${total} clockwise turns of size 2a as ${split}. At each ordinary window side prescribe the endpoint turn Eᵢ = π−4a${odd ? ", and use π−9a/2 at the first side for W*(a)" : ""}.`,
    "With kᵢ turns on chain i, set φᵢ = Eᵢ + Eᵢ₊₁ − 2kᵢa and γᵢ = Eᵢ − kᵢa (cyclic indices). Draw consecutive supporting tangents to a unit circle with turns φᵢ. They bound a convex polygon since the turns sum to 2π and each lies between 0 and π.",
    "Trim each corner by a sufficiently short chord making angle γᵢ with the preceding side. Replace the chord by the circular-tangent chain of §3 with kᵢ clockwise turns. Take the chords small enough to remain in disjoint corner neighbourhoods. The aperture base angles are smaller than both adjacent half-interior angles of the untrimmed polygon, giving clearance for the windows.",
    odd ? `Scale the first trimmed side to length two. Use the half-height of the odd window W*(a) there and fit W(a) on the other ${g-1} sides. The other side lengths exceed the required waist shrink h sin(a).` : `Choose a sufficiently small common prism half-height and insert an ordinary window W(a) on each of the ${g} trimmed sides.`,
    `The prism has ${4*g+2*total} vertices. The ${g} windows add ${4*g+odd}, giving the requested ${v} vertices of common defect −2a.`
  ], "§3 · Window and planar-chain lemmas (tangential base placement)");
}
// Constructive specialisations of the immersed non-orientable families in §§4–5.
// Only a closed locally injective surface counts as a CCP here.
function nonorientableRecipe(g, v) {
  if (!Number.isSafeInteger(g) || !Number.isSafeInteger(v) || g < 1 || v < 9) return null;
  const make = (title, reference, ...steps) => ({title, reference, steps});
  const band = "§4 · Attaching a Möbius band";
  const spectrum = "§4 · Prescribing the vertex count";
  const bounds = "§5 · Vertex bounds";
  if (g === 1) {
    const m = v - 9;
    const steps = [
      "Set d = 2π/" + v + ". Choose the nine-vertex one-sided Möbius band M(d) of Appendix A. Its six interior vertices have defect d and its equilateral rim has turning −2d at each corner."
    ];
    if (m === 0) {
      steps.push("Cap the rim with its flat equilateral triangle. The cap contributes turning 2π/3 = 3d at every corner, giving defect d at all nine vertices.");
    } else if (m % 2 === 1) {
      steps.push("Place a right equilateral pyramid on the other side of the rim. Adjust its height until its apex defect is " + m + "d. Its lateral surface then contributes turning 3d at each rim corner.");
      if (m > 1) steps.push("Make " + ((m - 1) / 2) + (((m - 1) / 2) === 1 ? " shallow apex truncation." : " shallow apex truncations.") + " At each step choose the cutting normal on the median of the remaining isosceles normal triangle, splitting off two spherical triangles of area d.");
    } else {
      const k = m / 2;
      steps.push("Construct the ridge cap of §4. From rim vertex A take two rays symmetric about the angle bisector and mark P and Q at equal distances. Adjust their elevation and distance to give angle sum π−3d at A, B and C. Both ridge endpoints then have defect " + k + "d.");
      if (k > 1) steps.push("Bevel the remaining ridge " + (k - 1) + (k === 2 ? " time." : " times.") + " Choose the cutting normal on the common edge of the two normal triangles so that each cut splits off one defect d at each end.");
    }
    if (m) steps.push("Attach the cap along its triangular rim. Its " + m + " interior vertices and the nine band vertices give " + v + " vertices. At each seam vertex 3d−2d=d.");
    return make("Cap a Möbius band", spectrum + " · The projective plane", ...steps);
  }
  if (g === 2) {
    if (v < 15) return null;
    const steps = ["Choose the flat nine-vertex band M(0) with equilateral rim, six interior vertices and zero boundary turning.",
      "Reflect M(0) through its rim plane and glue along the three rim edges. The two interiors lie on opposite sides. The resulting flat Klein bottle has 6+3+6=15 vertices."];
    if (v > 15) steps.push("Insert one vertex inside a triangular face and join it to all three corners. Repeat for a total of " + (v - 15) + (v === 16 ? " insertion." : " insertions.") + " Each new vertex has defect zero and the old defects are preserved.");
    return make("Double a flat Möbius band", band + " · " + spectrum, ...steps);
  }

  function uniformCore(h) {
    return make("Attach Möbius bands to a convex core", "§4 · A sphere with g triangular holes",
      "Build the convex core S_" + h + " between a regular " + h + "-gon and a regular " + (2*h) + "-gon. Use the upper trapezoid angle (π+π/" + h + ")/3.",
      "Remove the interiors of its " + h + " distinguished vertex-disjoint equilateral triangles. The " + (3*h) + " boundary vertices all have turning 3d, where d = 2π(2−" + h + ")/" + (9*h) + ".",
      "Scale a copy of M(d) to each triangular rim and place it on the exterior side. Each seam has defect 3d−2d=d. The " + h + " bands add " + (6*h) + " interior vertices and one crosscap apiece, giving " + (9*h) + " vertices.");
  }
  if (v === 9*g) return uniformCore(g);
  if (g >= 5 && v === 9*(g-2)) {
    return make("Attach crosscaps to an 18-vertex surface", bounds + " · " + band,
      "Construct the 18-vertex embedded genus-two CCP of §3 with two ordinary windows and one interior chain vertex on its base. Its common defect is d = −2π/9.",
      "Cut " + (g-4) + " pairwise disjoint small equilateral holes in face interiors. At each new rim vertex the punctured face contributes turning −2π/3=3d.",
      "Fit a scaled M(d) to each rim on the opposite side of its face. Each attachment adds three rim vertices and six band vertices, with defect d throughout. The vertex count is 18+9·" + (g-4) + "=" + v + ".");
  }
  function smallPrism(h) {
    const three = h === 3;
    const n = three ? 17 : 24;
    const angles = three ? "(20,43,39)π/51" : "(14,17,17)π/24";
    const d = three ? "−2π/17" : "−π/6";
    return make("Insert a window and attach " + (three ? "one crosscap" : "two crosscaps"), bounds + " · Technical note, triangular prisms",
      "Take a triangular prism of height 2 with base exterior turns (E_A,E_B,E_C)=" + angles + ". Set d=" + d + ", θ=π/6−3d/4 and AB=4 sin θ.",
      "Insert the window W(−d/2) on BC. Its annulus adds one handle and four vertices.",
      three ? "Erect the rectangular pyramid on AB and replace its vertical equilateral face by M(d), on the exterior side of the supporting plane." : "Attach reflected rectangular pyramids to AB and AC. Replace one vertical equilateral face of each by M(d) and its mirror. Use the parameter rectangle and separation check in the technical note to keep the bands disjoint away from their common edge.",
      "Each pyramid and band contributes seven vertices and one crosscap, giving 6+4+" + (three ? "7" : "2·7") + "=" + n + " vertices.");
  }
  function obliqueSix() {
    return make("Drill an oblique octagonal tunnel", bounds + " · Technical note, oblique tunnels",
      "Start with the 16-vertex embedded genus-two model with two windows W(π/8).",
      "Choose the interior points on a window pentagon and the outer wall given in the technical note. Reflect a sufficiently small polar octagonal opening across the perpendicular bisector plane of their joining segment.",
      "Join corresponding rim points by parallel segments. The reflected rims have the same induced orientation, giving N₆. The sixteen new rim vertices have defect −π/4, equal to the original defects, for 32 vertices.");
  }
  function octahedral(h) {
    const reps = h === 18 ? "(0,1,4), (0,1,12), (0,5,8), (0,6,13)" : "(0,1,3,2), (0,1,4), (0,5,7), (0,7,13), (0,11,18)";
    return make("Generate the octahedral face orbits", bounds + " · Technical note, octahedral immersions",
      "Apply the 24 orientation-preserving signed coordinate permutations to (1,3,7). Number the distinct orbit points from 0 to 23 in lexicographic order.",
      "Take all rotated images of the cyclic faces " + reps + ", identifying cyclic shifts and reversals.",
      "Use the face-incidence and vertex-link checks in the technical note. The connected non-orientable immersion has 24 vertices, " + (h === 18 ? "120 edges and 80 faces." : "132 edges and 86 faces."));
  }
  if (g >= 8 && g % 2 === 0 && v === 3*(g-2)) {
    const k = (g-8)/2;
    const steps = ["Start with Leopold’s 12-vertex genus-three polyhedron from Table 4, with defect −2π/3. Reserve the interior points (0,2,5) and (1,1,5) of the two faces specified in §5."];
    if (k) steps.push("Drill " + k + (k === 1 ? " ordinary triangular tunnel" : " ordinary triangular tunnels") + " parallel to the threefold axis, away from the reserved openings. Each adds one handle and six vertices of defect −2π/3.");
    steps.push("Cut small corresponding polar triangles around the reserved points, using the reflection (x,y,z) ↦ (y−1,x+1,z). Join their rims by parallel segments. This oblique tunnel produces N_" + g + " with 12+6·" + k + "+6=" + v + " vertices.");
    return make("Triangular drilling and an oblique tunnel", bounds + " · Technical note, oblique tunnels", ...steps);
  }

  function lift(base, h, k) {
    return make("Take a cyclic cover of N_" + h, bounds + " · " + base.reference,
      ...base.steps,
      "Write N_" + h + " as N_" + (h-2) + " # T². Cut along a meridian of the torus summand and join " + k + " copies of the cut surface cyclically.",
      "Lift its polygonal triangulation and compose the covering map with the immersion. The connected non-orientable cover has χ=" + (2-g) + " and exactly " + v + " vertices. Each lift retains the original defect.");
  }
  for (const [h,n] of nonorientableCoverBases) {
    const k = v/n;
    if (!Number.isInteger(k) || g !== 2+k*(h-2)) continue;
    const base = h <= 4 ? smallPrism(h) : h === 6 ? obliqueSix() : octahedral(h);
    return k === 1 ? base : lift(base,h,k);
  }
  const sheets = (v/9 - (g-2))/2;
  const baseGenus = 2 + (g-2)/sheets;
  if (v % 9 === 0 && Number.isInteger(sheets) && sheets > 1 && Number.isInteger(baseGenus) && baseGenus >= 3)
    return lift(uniformCore(baseGenus),baseGenus,sheets);
  return nonorientableChainRecipe(g,v);
}

function nonorientableChainRecipe(g,v) {
  // Equal-sector targets close by rotation. For the other counts test an explicit
  // instance of the short-chord construction without enumerating its K vertices.
  const pi = Math.PI, tau = 2*pi*(g-2)/v;
  if (!(tau > 0 && tau < 2*pi/9) || v < 11*g) return null;
  const collar = (v-11*g) % 2;
  const K = (v-11*g-3*collar)/2;
  if (K < 0 || !Number.isInteger(K)) return null;
  const uniform = collar === 0 && K % g === 0;
  function sector(r,k) {
    const theta = pi/6+(3+r)*tau/4;
    const B = pi/6-(7+r)*tau/4, A = B+2*pi/3-(2+r)*tau;
    const beta = k*tau/2;
    return {A,B,beta,gamma:B-beta,length:4*Math.sin(theta),theta};
  }
  const margin = 1e-8;
  const good = s => s.theta < pi/3-margin && s.A > pi/2+margin && s.A < pi-margin &&
    s.B > margin && s.B < pi-margin && s.beta < pi/2-margin;
  const quotient = Math.floor(K/g), rem = K % g;
  let closure;
  if (uniform) {
    const s = sector(0,quotient);
    if (!good(s) || 2*pi/g-s.gamma <= margin || 2*pi/g-s.gamma >= pi-margin) return null;
    closure = "Use the same short chord length at every connection. Rotating this sector through successive multiples of 2π/" + g + " closes the base. Take the chords sufficiently short that their chains occupy disjoint neighbourhoods of the regular polygon’s corners.";
  } else if (g > 96) {
    // A uniform sufficient bound for the large-g branch. Compare direction i with
    // 2πi/g, and chord direction i with 2πi/g+π/g−π/3. Balanced flooring has
    // prefix discrepancy <1. The single collar perturbs the cumulative turns by
    // at most 1.5τ, so side headings differ by ≤2.5τ, chord headings by <4.5τ
    // and lengths L by ≤4τ from 2. A common chord length 1/8 leaves residual
    // <10gτ. The two chosen chords have determinant of magnitude >1/2, making
    // both corrections <20gτ. Prefix displacements differ by <50gτ from the
    // regular comparison polygon with sides 2tᵢ⁰+uᵢ⁰/8 (length ≥1.875).
    // For gτ<.001 every chain lies within .36 of its ideal corner by the circle
    // bound below. Each straight side lies within .19 of its ideal side.
    // Nonincident features are separated: nonadjacent edges of the regular
    // polygon are at distance ≥1.875. At consecutive sides, 0<φ,φ−γ<π,
    // since φ differs from 2π/g by ≤2.5τ and γ from π/g−π/3 by <2τ.
    if (g*tau >= 0.001) return null;
    const other = Math.round(g/4)+1;
    closure = "Start with every chord of length 1/8. Adjust chords 1 and " + other +
      " by solving the two planar closure equations ∑(Lᵢtᵢ+ℓᵢuᵢ)=0. Here tᵢ is the side direction and uᵢ is rotated from it by E_B−kᵢτ/2. The inequality gτ<0.001 bounds each correction by 20gτ and keeps the chains in separate corner neighbourhoods.";
  } else {
    const sectors = [];
    for (let i=0; i<g; i++) {
      const extra = Math.floor((i+1)*rem/g)-Math.floor(i*rem/g);
      sectors.push(sector(i === 0 ? collar : 0,quotient+extra));
    }
    const result = nonorientableChordClosure(sectors,tau,good);
    if (!result) return null;
    closure = "Give each chord length " + result.length + " in prism coordinates. Adjust chords 1 and " + (result.other+1) + " by the unique solution of the two planar closure equations ∑(Lᵢtᵢ+ℓᵢuᵢ)=0, where tᵢ and uᵢ are the prescribed side and chord directions. The positive lengths give separate chain neighbourhoods.";
  }
  return {
    title:"Close a prism with polygonal chains",
    reference:"§4 · Prescribing the vertex count · §3 · Closing the prism base",
    steps:[
      "Set d=2π(2−"+g+")/"+v+" and τ=−d. Use "+g+" selected prism sides. " +
        (collar ? "Use M₁(d) on side 1 and M₀(d)=M(d) elsewhere. Form M₁ by attaching a triangular frustum to M(d), choosing turning 3d at the joined rim." : "Use the ordinary band M₀(d)=M(d) on every side."),
      "Distribute K="+K+" interior vertices among the connecting chains, with kᵢ=⌊iK/g⌋−⌊(i−1)K/g⌋ for i=1,…,"+g+". Give each interior vertex a clockwise turn τ. For a band of type r, set θ=π/6−(3+r)d/4, side length L=4 sin θ, E_B=d+π/3−θ and E_A=E_B+2π/3+(2+r)d.",
      closure+" Realise each chain with first and last directions at ±kᵢτ/2 to its chord using the tangent-circle construction of §3.",
      "Erect the prism of height 2, add the specified rectangular pyramid to each selected side, and replace its vertical equilateral face by its band on the exterior side. The result has "+g+" crosscaps and 11·"+g+"+2·"+K+"+3·"+collar+"="+v+" vertices of defect d."
    ]
  };
}

function nonorientableChordClosure(sectors,tau,good) {
  const n=sectors.length, pi=Math.PI;
  const vec=a=>[Math.cos(a),Math.sin(a)];
  const sum=(a,b)=>[a[0]+b[0],a[1]+b[1]];
  const scale=(a,c)=>[a[0]*c,a[1]*c];
  const minus=(a,b)=>[a[0]-b[0],a[1]-b[1]];
  const det=(a,b)=>a[0]*b[1]-a[1]*b[0];
  const dot=(a,b)=>a[0]*b[0]+a[1]*b[1];
  const distance=(a,b)=>Math.hypot(...minus(a,b));
  const margin=1e-7;
  function distanceToSide(p,a,b) {
    const d=minus(b,a), t=Math.max(0,Math.min(1,dot(minus(p,a),d)/dot(d,d)));
    return distance(p,sum(a,scale(d,t)));
  }
  function separated(a,b,c,d) {
    const direction=minus(b,a);
    const x=det(direction,minus(c,a)), y=det(direction,minus(d,a));
    return (x>margin && y>margin) || (x < -margin && y < -margin);
  }
  const sides=[], chords=[];
  let heading=0, sideSum=[0,0], chordSum=[0,0];
  for (let i=0;i<n;i++) {
    const s=sectors[i], next=sectors[(i+1)%n];
    const phi=s.B-2*s.beta+next.A;
    if (!good(s) || phi<=margin || phi>=pi-margin ||
        phi-s.gamma<=margin || phi-s.gamma>=pi-margin) return null;
    const side=scale(vec(heading),s.length), chord=vec(heading+s.gamma);
    sides.push(side); chords.push(chord);
    sideSum=sum(sideSum,side); chordSum=sum(chordSum,chord);
    heading+=phi;
  }
  // Pick the most transverse second chord to avoid an ill-conditioned solve.
  let other=1;
  for (let i=2;i<n;i++)
    if (Math.abs(det(chords[0],chords[i]))>Math.abs(det(chords[0],chords[other]))) other=i;
  const area=det(chords[0],chords[other]);
  if (Math.abs(area)<0.25) return null;
  for (const length of [0.125,0.25,0.5,0.75,1,1.5,2,3]) {
    const residual=scale(sum(sideSum,scale(chordSum,length)),-1);
    const lengths=Array(n).fill(length);
    lengths[0]+=det(residual,chords[other])/area;
    lengths[other]+=det(chords[0],residual)/area;
    if (lengths.some(x=>x<=margin)) continue;
    let p=[0,0];
    const A=[],B=[],radius=[];
    for (let i=0;i<n;i++) {
      A.push(p);p=sum(p,sides[i]);B.push(p);
      p=sum(p,scale(chords[i],lengths[i]));
      const beta=sectors[i].beta;
      // The tangent-circle formula bounds the whole chain by the disc about B_i.
      const factor=(1/Math.cos(beta/2)+Math.tan(beta/2)/Math.cos(tau/2))/2;
      radius.push(lengths[i]*Math.max(1,factor));
    }
    if (Math.hypot(...p)>margin) continue;
    let clear=true;
    for (let i=0;i<n && clear;i++) {
      for (let j=0;j<n && clear;j++) {
        if (j!==i && j!==(i+1)%n && distanceToSide(B[i],A[j],B[j])<=radius[i]+margin) clear=false;
        if (j<=i) continue;
        if (distance(B[i],B[j])<=radius[i]+radius[j]+margin) clear=false;
        // At consecutive sides the signs of phi and phi−gamma give separation;
        // the terminal tangent half-planes separate each chain from its incident sides.
        if (j!==i+1 && !(i===0 && j===n-1) &&
            !separated(A[i],B[i],A[j],B[j]) && !separated(A[j],B[j],A[i],B[i])) clear=false;
      }
    }
    if (clear) return {length,other};
  }
  return null;
}
