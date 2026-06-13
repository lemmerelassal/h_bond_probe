package main

import (
	"fmt"
	"math"
	"runtime"
	"sort"
	"strings"
	"sync"
)

// methylCap appends a methyl C bonded to anchorIdx, placed at the correct
// hybridization angle from the existing bond (awayFrom→anchor).
// For sp2 anchors (amide N, imine N): 120° angle in the plane defined by
// inPlane (a third atom in the same plane, e.g. O or adjacent N).
// For sp3 anchors: 109.5° tetrahedral angle, perpendicular to the plane.
// If inPlane is zero, falls back to collinear (180°) placement.
func methylCap(placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond,
	anchorIdx int, anchor, awayFrom Vec3, inPlane ...Vec3) {
	incoming := anchor.Sub(awayFrom) // direction of existing bond arriving at anchor
	if incoming.Norm() < 1e-6 {
		incoming = Vec3{1, 0, 0}
	}
	incoming = incoming.unit()

	var capDir Vec3
	if len(inPlane) > 0 && inPlane[0].Norm() > 1e-6 {
		// sp2: place cap at 120° from the existing bond (awayFrom→anchor).
		// At anchor, the existing bond arrives from awayFrom (direction: -incoming).
		// Rotate -incoming by 120° in the molecular plane:
		//   capDir = (-incoming)*cos(120°) + planeVec*sin(120°) = incoming*0.5 + planeVec*0.866
		planeVec := inPlane[0].Sub(anchor)
		planeVec = planeVec.Sub(incoming.Scale(planeVec.Dot(incoming)))
		if planeVec.Norm() < 1e-6 {
			planeVec = perpendicular(incoming)
		}
		planeVec = planeVec.unit()
		capDir = incoming.Scale(0.5).Add(planeVec.Scale(0.866))
	} else {
		// Default: collinear extension (180°).
		capDir = incoming
	}
	capPos := anchor.Add(capDir.unit().Scale(1.47))
	capIdx := len(*placedPos) + 1
	*placedPos = append(*placedPos, capPos)
	*placedLabel = append(*placedLabel, "C")
	*bonds = append(*bonds, Bond{anchorIdx, capIdx, 1})
}

func placeArgAcids(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type guanGroup struct {
		NE, CZ, NH1, NH2 Vec3
		found             [4]bool
	}
	groups := map[residueKey]*guanGroup{}
	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "ARG" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &guanGroup{}
		}
		g := groups[key]
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "NE":
			g.NE, g.found[0] = a.Pos, true
		case "CZ":
			g.CZ, g.found[1] = a.Pos, true
		case "NH1":
			g.NH1, g.found[2] = a.Pos, true
		case "NH2":
			g.NH2, g.found[3] = a.Pos, true
		}
	}

	// Tagged heavy atoms supporting per-ARG exclusion during clash checks.
	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] || !g.found[2] || !g.found[3] {
			continue
		}

		// Local frame: Z along NE→CZ (approach direction), X in-plane toward NH1.
		zHat := g.CZ.Sub(g.NE).unit()
		nh1Perp := g.NH1.Sub(g.CZ)
		nh1Perp = nh1Perp.Sub(zHat.Scale(nh1Perp.Dot(zHat)))
		baseX := func() Vec3 {
			if nh1Perp.Norm() < 0.1 {
				return perpendicular(zHat)
			}
			yH := cross3(zHat, nh1Perp).unit()
			return cross3(yH, zHat).unit()
		}()

		buildAcid := func(xHat Vec3, theta float64) (Vec3, Vec3, Vec3, Vec3) {
			cosT, sinT := math.Cos(theta), math.Sin(theta)
			yHat := cross3(zHat, xHat)
			rx := xHat.Scale(cosT).Add(yHat.Scale(sinT))
			o1 := g.CZ.Add(rx.Scale(1.103)).Add(zHat.Scale(3.565))
			o2 := g.CZ.Sub(rx.Scale(1.103)).Add(zHat.Scale(3.565))
			c2 := g.CZ.Add(zHat.Scale(4.151))
			c1 := c2.Add(zHat.Scale(1.52))
			return c1, c2, o1, o2
		}

		found := false
		var bc1, bc2, bo1, bo2 Vec3
		for step := 0; step < 12 && !found; step++ {
			theta := float64(step) * math.Pi / 6
			c1, c2, o1, o2 := buildAcid(baseX, theta)
			if clash.clears(c1, "C", key) && clash.clears(c2, "C", key) && clash.clears(o1, "O", key) && clash.clears(o2, "O", key) {
				bc1, bc2, bo1, bo2 = c1, c2, o1, o2
				found = true
			}
		}
		if !found {
			continue
		}

		base := len(*placedPos) + 1
		*placedPos = append(*placedPos, bc1, bc2, bo1, bo2)
		*placedLabel = append(*placedLabel, "C", "C", "O", "O")
		*bonds = append(*bonds,
			Bond{base, base + 1, 1},     // C1–C2  single
			Bond{base + 1, base + 2, 2}, // C2=O1  double
			Bond{base + 1, base + 3, 1}, // C2–O2  single
		)
		placed++
	}
	return placed
}


// tetrazoleRing builds a planar tetrazolate ring (a carboxylate bioisostere)
// approaching a cationic residue atom.  The ring carbon C5 carries the linker
// (a stem carbon), so both point away from the cation (+zHat); the four ring
// nitrogens face the cation (−zHat) to make the salt-bridge / bidentate H-bond,
// exactly mirroring how the carboxylate oxygens engage.
//
//   anchor — residue atom the ring approaches (ARG CZ or LYS NZ)
//   zHat   — approach axis pointing from the residue outward toward the probe
//   xHat   — in-plane spread axis (rotated during the clash search)
//
// Returns the six atom positions in the order
//   [C5, N1, N2, N3, N4, stemC]
// with ring connectivity C5–N1–N2–N3–N4–C5 and stemC bonded to C5.
func tetrazoleRing(anchor, zHat, xHat Vec3) [6]Vec3 {
	const (
		ringR   = 1.13 // pentagon circumradius (side ≈ 1.33 Å)
		centerD = 4.30 // ring-centre distance from anchor along zHat
		stemLen = 1.47 // C5–stem(C) bond length
	)
	center := anchor.Add(zHat.Scale(centerD))
	// Regular pentagon, C5 apex at +zHat (away from cation); going round by 72°
	// gives N1..N4 so that C5 is bonded to both N1 (72°) and N4 (288°).
	vert := func(angDeg float64) Vec3 {
		a := angDeg * math.Pi / 180
		return center.Add(zHat.Scale(ringR * math.Cos(a))).Add(xHat.Scale(ringR * math.Sin(a)))
	}
	c5 := vert(0)
	stem := c5.Add(zHat.Scale(stemLen)) // exocyclic, continues outward from apex
	return [6]Vec3{c5, vert(72), vert(144), vert(216), vert(288), stem}
}

// appendTetrazole writes a tetrazolate ring (positions from tetrazoleRing) into
// the probe arrays using a Kekulé bond pattern that keeps every atom within its
// valence (so sanitizeValence never strips a ring bond).
func appendTetrazole(placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond, ring [6]Vec3) {
	base := len(*placedPos) + 1
	*placedPos = append(*placedPos, ring[0], ring[1], ring[2], ring[3], ring[4], ring[5])
	*placedLabel = append(*placedLabel, "C", "N", "N", "N", "N", "C")
	*bonds = append(*bonds,
		Bond{base, base + 1, 1},     // C5–N1 single
		Bond{base + 1, base + 2, 1}, // N1–N2 single
		Bond{base + 2, base + 3, 2}, // N2=N3 double
		Bond{base + 3, base + 4, 1}, // N3–N4 single
		Bond{base + 4, base, 2},     // N4=C5 double
		Bond{base, base + 5, 1},     // C5–stem single (linker handle)
	)
}

// placeArgTetrazoles places a tetrazolate near each ARG guanidinium as a
// carboxylate bioisostere — same anchor atoms and approach frame as
// placeArgAcids, but a CN4 ring instead of COO⁻.  Generated in addition to the
// acid so both pharmacophores are explored.  Returns the number placed.
func placeArgTetrazoles(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type guanGroup struct {
		NE, CZ, NH1, NH2 Vec3
		found            [4]bool
	}
	groups := map[residueKey]*guanGroup{}
	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "ARG" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &guanGroup{}
		}
		g := groups[key]
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "NE":
			g.NE, g.found[0] = a.Pos, true
		case "CZ":
			g.CZ, g.found[1] = a.Pos, true
		case "NH1":
			g.NH1, g.found[2] = a.Pos, true
		case "NH2":
			g.NH2, g.found[3] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] || !g.found[2] || !g.found[3] {
			continue
		}

		// Local frame: Z along NE→CZ (approach direction), X in-plane toward NH1.
		zHat := g.CZ.Sub(g.NE).unit()
		nh1Perp := g.NH1.Sub(g.CZ)
		nh1Perp = nh1Perp.Sub(zHat.Scale(nh1Perp.Dot(zHat)))
		baseX := func() Vec3 {
			if nh1Perp.Norm() < 0.1 {
				return perpendicular(zHat)
			}
			yH := cross3(zHat, nh1Perp).unit()
			return cross3(yH, zHat).unit()
		}()

		var ring [6]Vec3
		found := false
		for step := 0; step < 12 && !found; step++ {
			theta := float64(step) * math.Pi / 6
			yHat := cross3(zHat, baseX)
			rx := baseX.Scale(math.Cos(theta)).Add(yHat.Scale(math.Sin(theta)))
			r := tetrazoleRing(g.CZ, zHat, rx)
			ok := true
			for i, p := range r {
				elem := "N"
				if i == 0 || i == 5 {
					elem = "C"
				}
				if !clash.clears(p, elem, key) {
					ok = false
					break
				}
			}
			if ok {
				ring, found = r, true
			}
		}
		if !found {
			continue
		}
		appendTetrazole(placedPos, placedLabel, bonds, ring)
		placed++
	}
	return placed
}


// placeCarboxylateGuanidines places a guanidine molecule near each ASP or GLU
// carboxylate using bidentate H-bond geometry (N···O ≈ 2.9 Å to both oxygens).
//
// The guanidine C is placed 4.15 Å from the carboxylate C; N1/N2 point back
// toward O1/O2 at 120° angles; N3 (imine) points away.  Because the carboxylate
// plane can be approached from either face, and because a buried residue may only
// have a narrow clear window, we search:
//   • both ±zHat approach directions (either face of the carboxylate),
//   • 12 rotations (30° steps) of the whole guanidine around zHat, and
//   • the xHat basis is derived from the actual O1 position but falls back to an
//     arbitrary perpendicular if O1 is nearly collinear with zHat.
// The first candidate that clears all clash checks is accepted.
// Returns the number placed.
func placeCarboxylateGuanidines(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type carboxylGroup struct {
		prevC, carboxylC, O1, O2 Vec3
		found                    [4]bool
	}
	groups := map[residueKey]*carboxylGroup{}

	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		if res != "ASP" && res != "GLU" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &carboxylGroup{}
		}
		g := groups[key]
		switch res + ":" + strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "ASP:CB":
			g.prevC, g.found[0] = a.Pos, true
		case "ASP:CG":
			g.carboxylC, g.found[1] = a.Pos, true
		case "ASP:OD1":
			g.O1, g.found[2] = a.Pos, true
		case "ASP:OD2":
			g.O2, g.found[3] = a.Pos, true
		case "GLU:CG":
			g.prevC, g.found[0] = a.Pos, true
		case "GLU:CD":
			g.carboxylC, g.found[1] = a.Pos, true
		case "GLU:OE1":
			g.O1, g.found[2] = a.Pos, true
		case "GLU:OE2":
			g.O2, g.found[3] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] || !g.found[2] || !g.found[3] {
			continue
		}

		// Base approach axis: prevC→carboxylC.
		baseZ := g.carboxylC.Sub(g.prevC).unit()

		// xHat from O1 position, with fallback if O1 is nearly collinear with zHat.
		buildXHat := func(zHat Vec3) Vec3 {
			o1Vec := g.O1.Sub(g.carboxylC)
			o1Perp := o1Vec.Sub(zHat.Scale(o1Vec.Dot(zHat)))
			if o1Perp.Norm() < 0.1 {
				return perpendicular(zHat)
			}
			yH := cross3(zHat, o1Perp).unit()
			return cross3(yH, zHat).unit()
		}

		// buildGuanidine returns the five atom positions for a given zHat and xHat,
		// rotated by angle θ around zHat.
		// N3 is the imine N pointing away from the carboxylate; stemC is a carbon
		// bonded to N3 that serves as the linker attachment point, preventing the
		// linker from connecting to the H-bond-donor nitrogens N1/N2.
		buildGuanidine := func(zHat, xHat Vec3, theta float64) (Vec3, Vec3, Vec3, Vec3, Vec3) {
			cosT, sinT := math.Cos(theta), math.Sin(theta)
			// Rotate xHat around zHat by theta.
			yHat := cross3(zHat, xHat)
			rx := xHat.Scale(cosT).Add(yHat.Scale(sinT))
			guanC := g.carboxylC.Add(zHat.Scale(4.15))
			n1 := guanC.Add(rx.Scale(1.152)).Sub(zHat.Scale(0.665))
			n2 := guanC.Sub(rx.Scale(1.152)).Sub(zHat.Scale(0.665))
			n3 := guanC.Add(zHat.Scale(1.33))
			stemC := n3.Add(zHat.Scale(1.47)) // C–N single bond, linker handle
			return guanC, n1, n2, n3, stemC
		}

		found := false
		var bestC, bestN1, bestN2, bestN3, bestStemC Vec3
		// Try both faces (±zHat) × 12 rotation steps (30°).
		for _, sign := range []float64{1, -1} {
			if found {
				break
			}
			zHat := baseZ.Scale(sign)
			xHat := buildXHat(zHat)
			for step := 0; step < 12; step++ {
				theta := float64(step) * math.Pi / 6
				guanC, n1, n2, n3, stemC := buildGuanidine(zHat, xHat, theta)
				if clash.clears(guanC, "C", key) && clash.clears(n1, "N", key) && clash.clears(n2, "N", key) && clash.clears(n3, "N", key) && clash.clears(stemC, "C", key) {
					bestC, bestN1, bestN2, bestN3, bestStemC = guanC, n1, n2, n3, stemC
					found = true
					break
				}
			}
		}
		if !found {
			continue
		}

		base := len(*placedPos) + 1
		*placedPos = append(*placedPos, bestC, bestN1, bestN2, bestN3, bestStemC)
		*placedLabel = append(*placedLabel, "C", "N", "N", "N", "C")
		*bonds = append(*bonds,
			Bond{base, base + 1, 1},     // C–N1 single
			Bond{base, base + 2, 1},     // C–N2 single
			Bond{base, base + 3, 2},     // C=N3 double (imine)
			Bond{base + 3, base + 4, 1}, // N3–stemC single (linker handle)
		)
		placed++
	}
	return placed
}


// placeTyrAcetamides places an acetamide (CH3–CO–NH2) near each TYR phenolic OH.
//
// Primary interaction: carbonyl O accepts TYR-OH's H-bond donation at 2.8 Å
// along the CZ→OH axis (linear O–H···O=C geometry).
//
// Molecule orientation: acetamide lies in the TYR ring plane.  C2→N is at 122°
// from C2→O (= −zHat), choosing the ±xHat half-space that avoids clashes.
// C2→C1 (methyl) is on the opposite side to N.
//
// Atom order appended: C1 (methyl), C2 (carbonyl C), O (carbonyl O), N (amide N).
// Bonds: C1–C2 single, C2=O double, C2–N single.
func placeTyrAcetamides(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type tyrGroup struct {
		CZ, OH, CE1 Vec3
		found       [3]bool
	}
	groups := map[residueKey]*tyrGroup{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "TYR" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &tyrGroup{}
		}
		g := groups[key]
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "CZ":
			g.CZ, g.found[0] = a.Pos, true
		case "OH":
			g.OH, g.found[1] = a.Pos, true
		case "CE1":
			g.CE1, g.found[2] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] || !g.found[2] {
			continue
		}

		// Local frame: Z along CZ→OH (H points roughly this way), X in ring plane toward CE1.
		zHat := g.OH.Sub(g.CZ).unit()
		ce1Vec := g.CE1.Sub(g.CZ)
		ce1Perp := ce1Vec.Sub(zHat.Scale(ce1Vec.Dot(zHat)))
		yHat := cross3(zHat, ce1Perp).unit()
		xHat := cross3(yHat, zHat).unit()

		// Phenol C–O–H angle ≈ 109°; O→C = −zHat, so O→H direction in the ring plane:
		//   dH = cos(71°)·zHat ± sin(71°)·xHat  =  0.326·zHat ± 0.945·xHat
		// (cos/sin of 71° = 180°−109°).  The ± chooses which side of the ring the H is on.
		// Carbonyl O placed at OH + dH·2.8 (O···O H-bond); C2 a further 1.24 Å along dH.
		// N and C1 placed in the ring plane at 122° from C2→O; try all 4 orientations.

		type acPlacement struct{ O, C2, N, C1 Vec3 }
		var candidates []acPlacement
		for _, hSign := range []float64{1, -1} {
			// dH: estimated O–H bond direction in the ring plane.
			dH := zHat.Scale(0.326).Add(xHat.Scale(hSign * 0.945))   // unit (0.326²+0.945²≈1)
			dP := zHat.Scale(-0.945).Add(xHat.Scale(hSign * 0.326))  // perpendicular in plane
			aO := g.OH.Add(dH.Scale(2.8))
			aC2 := g.OH.Add(dH.Scale(4.04)) // C2 = O + 1.24 Å along dH
			for _, nSign := range []float64{1, -1} {
				// C2→N at 122° from C2→O (= −dH) in the molecular plane.
				aN := aC2.Add(dH.Scale(1.33 * 0.530)).Add(dP.Scale(nSign * 1.33 * 0.848))
				aC1 := aC2.Add(dH.Scale(1.52 * 0.530)).Sub(dP.Scale(nSign * 1.52 * 0.848))
				candidates = append(candidates, acPlacement{aO, aC2, aN, aC1})
			}
		}

		var chosen acPlacement
		placed2 := false
		for _, c := range candidates {
			if clash.clears(c.O, "O", key) && clash.clears(c.C2, "C", key) && clash.clears(c.N, "N", key) && clash.clears(c.C1, "C", key) {
				chosen = c
				placed2 = true
				break
			}
		}
		if !placed2 {
			continue
		}
		aO, aC2, aN, aC1 := chosen.O, chosen.C2, chosen.N, chosen.C1

		base := len(*placedPos) + 1
		*placedPos = append(*placedPos, aC1, aC2, aO, aN)
		*placedLabel = append(*placedLabel, "C", "C", "O", "N")
		*bonds = append(*bonds,
			Bond{base, base + 1, 1},     // C1–C2 single
			Bond{base + 1, base + 2, 2}, // C2=O  double
			Bond{base + 1, base + 3, 1}, // C2–N  single
		)
		// C1 is the methyl carbon (sp3, bondCount=1) — it is already the natural
		// linker tail without any additional cap.  No N cap needed.
		placed++
	}
	return placed
}


// placeTyrPhenols places a phenol probe near each TYR, capturing the combined
// π-stacking + OH···HO-TYR pharmacophore.
//
// The probe ring is parallel to TYR's aromatic ring at 3.5 Å (standard
// π-stacking distance), with the para carbon (C0) and OH aligned along the
// CG→CZ→OH axis.  The resulting O···O distance is ≈ 3.5 Å — at the outer
// edge of an H-bond but appropriate for a bifunctional probe where stacking is
// the primary constraint.
//
// Atom order: C0–C5 (ring, C0 = para-C bearing OH), O (para-OH).
// Bonds: C0–C1…C5–C0 aromatic (order 4), C0–O single.
func placeTyrPhenols(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type tyrGroup struct {
		ring  [6]Vec3 // CG, CD1, CD2, CE1, CE2, CZ
		CZ    Vec3
		found [7]bool // [6 ring atoms] + CZ duplicate for dPara
	}
	// Atom-name → ring slot index (same 6 as the π-stacking code collects).
	ringSlot := map[string]int{
		"CG": 0, "CD1": 1, "CD2": 2, "CE1": 3, "CE2": 4, "CZ": 5,
	}
	groups := map[residueKey]*tyrGroup{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "TYR" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &tyrGroup{}
		}
		g := groups[key]
		name := strings.TrimSpace(strings.ToUpper(a.Name))
		if idx, ok := ringSlot[name]; ok {
			g.ring[idx] = a.Pos
			g.found[idx] = true
			if name == "CZ" {
				g.CZ = a.Pos
			}
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		// Require all 6 ring atoms.
		allFound := true
		for i := 0; i < 6; i++ {
			if !g.found[i] {
				allFound = false
				break
			}
		}
		if !allFound {
			continue
		}

		// Centroid and ring normal using the same functions as the π-stacking code.
		centroid := Vec3{}
		for _, p := range g.ring {
			centroid = centroid.Add(p)
		}
		centroid = centroid.Scale(1.0 / 6)
		ringNormal := ringPlaneNormal(g.ring[:])

		// Para direction: centroid → CZ (projected into ring plane for exactness).
		dParaRaw := g.CZ.Sub(centroid)
		dPara := dParaRaw.Sub(ringNormal.Scale(dParaRaw.Dot(ringNormal))).unit()
		dPerp := cross3(ringNormal, dPara).unit()

		for _, sign := range []float64{1, -1} {
			phenolCenter := centroid.Add(ringNormal.Scale(sign * 3.5))

			// 6 ring carbons: C0 at angle 0° (dPara direction, the para-C bearing OH).
			var ringPts [6]Vec3
			for i := 0; i < 6; i++ {
				angle := float64(i) * math.Pi / 3
				ringPts[i] = phenolCenter.
					Add(dPara.Scale(1.40 * math.Cos(angle))).
					Add(dPerp.Scale(1.40 * math.Sin(angle)))
			}
			// OH extends from C0 in dPara: ring C at 1.40 Å, O a further 1.36 Å.
			phenolOH := phenolCenter.Add(dPara.Scale(1.40 + 1.36))

			// Centre + OH check only (same policy as the π-stacking probe code).
			if !clash.clears(phenolCenter, "C", key) || !clash.clears(phenolOH, "O", key) {
				continue
			}

			base := len(*placedPos) + 1
			for _, p := range ringPts {
				*placedPos = append(*placedPos, p)
				*placedLabel = append(*placedLabel, "C")
			}
			*placedPos = append(*placedPos, phenolOH)
			*placedLabel = append(*placedLabel, "O")
			for i := 0; i < 6; i++ {
				*bonds = append(*bonds, Bond{base + i, base + (i+1)%6, 4})
			}
			*bonds = append(*bonds, Bond{base, base + 6, 1}) // C0–O

			// Methyl cap on ring C3 (para to OH): radial outward from probe ring centre.
			// Ring C is sp2 with 2 ring bonds already; the cap goes straight radially.
			capDir3 := ringPts[3].Sub(phenolCenter).unit()
			capPos3 := ringPts[3].Add(capDir3.Scale(1.54))
			capIdx3 := len(*placedPos) + 1
			*placedPos = append(*placedPos, capPos3)
			*placedLabel = append(*placedLabel, "C")
			*bonds = append(*bonds, Bond{base + 3, capIdx3, 1})

			placed++
			break
		}
	}
	return placed
}


// placeAsnGlnAcetamides places a formamide-like probe (N–CO–N) near each ASN
// or GLN amide group, capturing the bidentate donor+acceptor pharmacophore.
//
// Geometry: the probe carbonyl O is placed ~2.8 Å from ND2/NE2 (accepts the
// NH donation) and a probe amide N is placed ~3.0 Å from OD1/OE1 (donates
// to the carbonyl acceptor). The probe molecule is H2N–CO–NH2 (urea-like),
// placed so that:
//   - probeC sits 4.1 Å from the residue amide C along the Cprev→CamideC axis
//   - probeO points back toward ND2/NE2 at 2.8 Å
//   - probeN (imine) points back toward OD1/OE1 at 3.0 Å
//
// Atom order: probeC, probeO, probeN1 (toward ND2), probeN2 (imine, toward OD1).
// Bonds: probeC=probeO double, probeC–probeN1 single, probeC–probeN2 single.
func placeAsnGlnAcetamides(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type amideGroup struct {
		prevC, amideC, O, N Vec3
		found               [4]bool
	}
	groups := map[residueKey]*amideGroup{}

	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		if res != "ASN" && res != "GLN" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &amideGroup{}
		}
		g := groups[key]
		tok := res + ":" + strings.TrimSpace(strings.ToUpper(a.Name))
		switch tok {
		case "ASN:CB":
			g.prevC, g.found[0] = a.Pos, true
		case "ASN:CG":
			g.amideC, g.found[1] = a.Pos, true
		case "ASN:OD1":
			g.O, g.found[2] = a.Pos, true
		case "ASN:ND2":
			g.N, g.found[3] = a.Pos, true
		case "GLN:CG":
			g.prevC, g.found[0] = a.Pos, true
		case "GLN:CD":
			g.amideC, g.found[1] = a.Pos, true
		case "GLN:OE1":
			g.O, g.found[2] = a.Pos, true
		case "GLN:NE2":
			g.N, g.found[3] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] || !g.found[2] || !g.found[3] {
			continue
		}

		// Local frame: Z along prevC→amideC, X toward O (in the amide plane).
		// Try both ±zHat faces × 12 rotation steps to handle buried residues.
		baseZ := g.amideC.Sub(g.prevC).unit()
		oVec := g.O.Sub(g.amideC)
		buildBaseX := func(zHat Vec3) Vec3 {
			oPerp := oVec.Sub(zHat.Scale(oVec.Dot(zHat)))
			if oPerp.Norm() < 0.1 {
				return perpendicular(zHat)
			}
			yH := cross3(zHat, oPerp).unit()
			return cross3(yH, zHat).unit()
		}

		buildAmide := func(zHat, xHat Vec3, theta float64) (Vec3, Vec3, Vec3, Vec3) {
			cosT, sinT := math.Cos(theta), math.Sin(theta)
			yHat := cross3(zHat, xHat)
			rx := xHat.Scale(cosT).Add(yHat.Scale(sinT))
			probeC := g.amideC.Add(zHat.Scale(4.10))
			probeO := probeC.Add(rx.Scale(1.103)).Sub(zHat.Scale(0.665))
			probeN1 := probeC.Sub(rx.Scale(1.103)).Sub(zHat.Scale(0.665))
			probeN2 := probeC.Add(zHat.Scale(1.24))
			return probeC, probeO, probeN1, probeN2
		}

		found := false
		var bpC, bpO, bpN1, bpN2 Vec3
		for _, sign := range []float64{1, -1} {
			if found {
				break
			}
			zHat := baseZ.Scale(sign)
			xHat := buildBaseX(zHat)
			for step := 0; step < 12 && !found; step++ {
				theta := float64(step) * math.Pi / 6
				pC, pO, pN1, pN2 := buildAmide(zHat, xHat, theta)
				if clash.clears(pC, "C", key) && clash.clears(pO, "O", key) && clash.clears(pN1, "N", key) && clash.clears(pN2, "N", key) {
					bpC, bpO, bpN1, bpN2 = pC, pO, pN1, pN2
					found = true
				}
			}
		}
		if !found {
			continue
		}

		base := len(*placedPos) + 1
		*placedPos = append(*placedPos, bpC, bpO, bpN1, bpN2)
		*placedLabel = append(*placedLabel, "C", "O", "N", "N")
		*bonds = append(*bonds,
			Bond{base, base + 1, 2}, // C=O  double (carbonyl)
			Bond{base, base + 2, 1}, // C–N1 single (amide N toward residue N)
			Bond{base, base + 3, 1}, // C–N2 single (amino N, away)
		)
		// Methyl caps on both N1 and N2 so linkProbeGroups picks the shortest path.
		methylCap(placedPos, placedLabel, bonds, base+2, bpN1, bpC, bpN2) // N1–CH3
		methylCap(placedPos, placedLabel, bonds, base+3, bpN2, bpC, bpN1) // N2–CH3 (was sp2 120°)
		placed++
	}
	return placed
}


// placeSerThrMethanols places a methanol probe (C–OH) near each SER or THR
// hydroxyl oxygen, capturing the dual donor/acceptor pharmacophore of the
// hydroxyl group.
//
// Geometry: the probe O is placed 2.8 Å from the residue OG/OG1 along the
// CB→OG axis (linear H-bond geometry). The probe methyl C is 1.43 Å further
// along the same axis (C–O bond length). Returns the number placed.
func placeSerThrMethanols(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type hydroxyl struct {
		CB, OG Vec3
		found  [2]bool
	}
	groups := map[residueKey]*hydroxyl{}

	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		if res != "SER" && res != "THR" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &hydroxyl{}
		}
		g := groups[key]
		name := strings.TrimSpace(strings.ToUpper(a.Name))
		switch name {
		case "CB":
			g.CB, g.found[0] = a.Pos, true
		case "OG", "OG1": // SER uses OG; THR uses OG1
			g.OG, g.found[1] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] {
			continue
		}

		// Direction: CB → OG extended (probe is the donor/acceptor partner of the hydroxyl).
		// Try both ±zHat: +zHat captures H-bond donors to the OH, −zHat captures acceptors.
		baseZ := g.OG.Sub(g.CB).unit()

		placed2 := false
		for _, sign := range []float64{1, -1} {
			zHat := baseZ.Scale(sign)
			probeO := g.OG.Add(zHat.Scale(2.80))
			probeC := g.OG.Add(zHat.Scale(2.80 + 1.43))
			if clash.clears(probeO, "O", key) && clash.clears(probeC, "C", key) {
				base := len(*placedPos) + 1
				*placedPos = append(*placedPos, probeO, probeC)
				*placedLabel = append(*placedLabel, "O", "C")
				*bonds = append(*bonds, Bond{base, base + 1, 1})
				placed2 = true
				break
			}
		}
		if !placed2 {
			continue
		}
		placed++
	}
	return placed
}


// placeLysAcetates places a bidentate acetate (acetic acid) near each LYS NZ
// ammonium, mirroring the geometry used for ARG guanidinium.
//
// Geometry: carboxylate C2 is 3.0 Å from NZ along the CE→NZ axis;
// O1/O2 sit at ±1.10 Å perpendicular and 0.59 Å further, giving N···O ≈ 2.9 Å.
// Methyl C1 is 1.52 Å beyond C2. Returns the number placed.
func placeLysAcetates(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type lys struct {
		CE, NZ Vec3
		found  [2]bool
	}
	groups := map[residueKey]*lys{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "LYS" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &lys{}
		}
		g := groups[key]
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "CE":
			g.CE, g.found[0] = a.Pos, true
		case "NZ":
			g.NZ, g.found[1] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] {
			continue
		}

		zHat := g.NZ.Sub(g.CE).unit()
		baseX := perpendicular(zHat)

		buildAcetate := func(xHat Vec3, theta float64) (Vec3, Vec3, Vec3, Vec3) {
			cosT, sinT := math.Cos(theta), math.Sin(theta)
			yHat := cross3(zHat, xHat)
			rx := xHat.Scale(cosT).Add(yHat.Scale(sinT))
			o1 := g.NZ.Add(rx.Scale(1.103)).Add(zHat.Scale(3.565))
			o2 := g.NZ.Sub(rx.Scale(1.103)).Add(zHat.Scale(3.565))
			c2 := g.NZ.Add(zHat.Scale(4.151))
			c1 := c2.Add(zHat.Scale(1.52))
			return c1, c2, o1, o2
		}

		found := false
		var bc1, bc2, bo1, bo2 Vec3
		for step := 0; step < 12 && !found; step++ {
			theta := float64(step) * math.Pi / 6
			c1, c2, o1, o2 := buildAcetate(baseX, theta)
			if clash.clears(c1, "C", key) && clash.clears(c2, "C", key) && clash.clears(o1, "O", key) && clash.clears(o2, "O", key) {
				bc1, bc2, bo1, bo2 = c1, c2, o1, o2
				found = true
			}
		}
		if !found {
			continue
		}

		base := len(*placedPos) + 1
		*placedPos = append(*placedPos, bc1, bc2, bo1, bo2)
		*placedLabel = append(*placedLabel, "C", "C", "O", "O")
		*bonds = append(*bonds,
			Bond{base, base + 1, 1},     // C1–C2 single
			Bond{base + 1, base + 2, 2}, // C2=O1 double
			Bond{base + 1, base + 3, 1}, // C2–O2 single
		)
		placed++
	}
	return placed
}


// placeLysTetrazoles places a tetrazolate near each LYS ammonium as a
// carboxylate bioisostere — same anchor (NZ) and approach axis as
// placeLysAcetates, but a CN4 ring instead of COO⁻.  Generated in addition to
// the acetate.  Returns the number placed.
func placeLysTetrazoles(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type lys struct {
		CE, NZ Vec3
		found  [2]bool
	}
	groups := map[residueKey]*lys{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "LYS" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &lys{}
		}
		g := groups[key]
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "CE":
			g.CE, g.found[0] = a.Pos, true
		case "NZ":
			g.NZ, g.found[1] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] {
			continue
		}

		zHat := g.NZ.Sub(g.CE).unit()
		baseX := perpendicular(zHat)

		var ring [6]Vec3
		found := false
		for step := 0; step < 12 && !found; step++ {
			theta := float64(step) * math.Pi / 6
			yHat := cross3(zHat, baseX)
			rx := baseX.Scale(math.Cos(theta)).Add(yHat.Scale(math.Sin(theta)))
			r := tetrazoleRing(g.NZ, zHat, rx)
			ok := true
			for i, p := range r {
				elem := "N"
				if i == 0 || i == 5 {
					elem = "C"
				}
				if !clash.clears(p, elem, key) {
					ok = false
					break
				}
			}
			if ok {
				ring, found = r, true
			}
		}
		if !found {
			continue
		}
		appendTetrazole(placedPos, placedLabel, bonds, ring)
		placed++
	}
	return placed
}


// placeCysMethanethiols places a methanethiol probe (C–SH) near each CYS SG,
// capturing the thiol donor/acceptor pharmacophore.
//
// Geometry: the probe S is placed 3.6 Å from residue SG along the CB→SG axis
// (S···S or S···N distance at outer H-bond range). Probe C is 1.82 Å further.
// Returns the number placed.
func placeCysMethanethiols(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type thiol struct {
		CB, SG Vec3
		found  [2]bool
	}
	groups := map[residueKey]*thiol{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "CYS" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &thiol{}
		}
		g := groups[key]
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "CB":
			g.CB, g.found[0] = a.Pos, true
		case "SG":
			g.SG, g.found[1] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] {
			continue
		}

		baseZ := g.SG.Sub(g.CB).unit()

		placed2 := false
		for _, sign := range []float64{1, -1} {
			zHat := baseZ.Scale(sign)
			probeS := g.SG.Add(zHat.Scale(3.60))
			probeC := g.SG.Add(zHat.Scale(3.60 + 1.82))
			if clash.clears(probeS, "S", key) && clash.clears(probeC, "C", key) {
				base := len(*placedPos) + 1
				*placedPos = append(*placedPos, probeS, probeC)
				*placedLabel = append(*placedLabel, "S", "C")
				*bonds = append(*bonds, Bond{base, base + 1, 1})
				placed2 = true
				break
			}
		}
		if !placed2 {
			continue
		}
		placed++
	}
	return placed
}


// placeHisImidazoles places an imidazole probe near each HIS residue, capturing
// the combined π-stacking + tautomeric H-bond pharmacophore of the imidazole ring.
//
// Strategy: a five-membered imidazole ring (C3N2, bond lengths N–C ≈ 1.38 Å,
// C–C ≈ 1.37 Å, C=N ≈ 1.32 Å) is placed parallel to the HIS ring at 3.5 Å
// (π-stacking distance), analogous to placeTyrPhenols. Both faces are tried;
// the one with fewer clashes is accepted. Returns the number placed.
func placeHisImidazoles(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type hisGroup struct {
		ring  [5]Vec3 // CG, ND1, CD2, CE1, NE2
		found [5]bool
	}
	ringSlot := map[string]int{
		"CG": 0, "ND1": 1, "CD2": 2, "CE1": 3, "NE2": 4,
	}
	groups := map[residueKey]*hisGroup{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "HIS" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &hisGroup{}
		}
		g := groups[key]
		name := strings.TrimSpace(strings.ToUpper(a.Name))
		if idx, ok := ringSlot[name]; ok {
			g.ring[idx] = a.Pos
			g.found[idx] = true
		}
	}

	clash := newResClashSet(atoms)

	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		allFound := true
		for i := 0; i < 5; i++ {
			if !g.found[i] {
				allFound = false
				break
			}
		}
		if !allFound {
			continue
		}

		// Ring centroid and normal from the 5 heavy atoms.
		centroid := Vec3{}
		for _, p := range g.ring {
			centroid = centroid.Add(p)
		}
		centroid = centroid.Scale(1.0 / 5)
		ringNormal := ringPlaneNormal(g.ring[:])

		// Imidazole probe ring geometry (5-membered, irregular):
		// Approximate as a regular pentagon with inner radius 1.15 Å (C–C ≈ 1.37 Å),
		// labelled C, N, C, C, N in the canonical imidazole order.
		// Atom indices: 0=C2(between two N), 1=N3, 2=C4, 3=C5, 4=N1.
		probeLabels := [5]string{"C", "N", "C", "C", "N"}
		const innerR = 1.15
		u := perpendicular(ringNormal)
		v := cross3(ringNormal, u)

		for _, sign := range []float64{1, -1} {
			centre := centroid.Add(ringNormal.Scale(sign * 3.5))
			if !clash.clears(centre, "C", key) {
				continue
			}

			ringPts := [5]Vec3{}
			for i := 0; i < 5; i++ {
				angle := float64(i) * 2 * math.Pi / 5
				ringPts[i] = centre.
					Add(u.Scale(innerR * math.Cos(angle))).
					Add(v.Scale(innerR * math.Sin(angle)))
			}

			base := len(*placedPos) + 1
			for i, p := range ringPts {
				*placedPos = append(*placedPos, p)
				*placedLabel = append(*placedLabel, probeLabels[i])
			}
			// Bonds: ring closure (aromatic order 4).
			for i := 0; i < 5; i++ {
				*bonds = append(*bonds, Bond{base + i, base + (i+1)%5, 4})
			}
			// Methyl cap on C4: radial outward from ring centre (correct for sp2 ring C).
			capDir4 := ringPts[3].Sub(centre).unit()
			capPos4 := ringPts[3].Add(capDir4.Scale(1.54))
			capIdx4 := len(*placedPos) + 1
			*placedPos = append(*placedPos, capPos4)
			*placedLabel = append(*placedLabel, "C")
			*bonds = append(*bonds, Bond{base + 3, capIdx4, 1})
			placed++
			break
		}
	}
	return placed
}


// placeWaterMethanols places a methanol probe (O–C) at every crystallographic
// water oxygen in the binding site, representing the pharmacophore of a
// displaceable water.
//
// Unlike SER/THR, there is no CB axis to define the approach direction, so
// the methyl C is placed in the direction that maximises clearance from all
// non-water protein atoms.  The probe O sits exactly at the water oxygen
// position (displacement = 0); the methyl C is 1.43 Å away in the best
// direction, chosen from 26 uniformly distributed directions.
//
// Waters that are completely buried (no direction clears for the methyl C)
// are still placed as a bare O–C pair in the least-clashing orientation,
// because any ligand that reaches there will displace the water regardless.
// Returns the number of probes placed.
func placeWaterMethanols(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	// Collect water oxygens.
	type waterO struct {
		pos     Vec3
		chainID string
		resSeq  int
	}
	var waters []waterO
	for _, a := range atoms {
		if !isWater(a) {
			continue
		}
		if strings.ToUpper(a.Element) != "O" {
			continue
		}
		waters = append(waters, waterO{a.Pos, a.ChainID, a.ResSeq})
	}

	// Build non-water heavy atom list for clash checks.
	type heavyAtomW struct {
		pos  Vec3
		vdwR float64
	}
	var heavy []heavyAtomW
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, heavyAtomW{a.Pos, vdw(a.Element)})
	}

	// 26 directions: 6 faces + 12 edges + 8 corners of a unit cube, normalised.
	dirs := make([]Vec3, 0, 26)
	for dx := -1; dx <= 1; dx++ {
		for dy := -1; dy <= 1; dy++ {
			for dz := -1; dz <= 1; dz++ {
				if dx == 0 && dy == 0 && dz == 0 {
					continue
				}
				d := Vec3{float64(dx), float64(dy), float64(dz)}
				dirs = append(dirs, d.Scale(1.0/d.Norm()))
			}
		}
	}

	placed := 0
	for _, w := range waters {
		// Find the direction for the methyl C that maximises clearance.
		bestDir := dirs[0]
		bestClear := -999.0
		for _, d := range dirs {
			cand := w.pos.Add(d.Scale(1.43))
			minClear := 999.0
			for _, h := range heavy {
				gap := cand.Sub(h.pos).Norm() - h.vdwR - vdw("C")
				if gap < minClear {
					minClear = gap
				}
			}
			if minClear > bestClear {
				bestClear = minClear
				bestDir = d
			}
		}

		probeO := w.pos
		probeC := w.pos.Add(bestDir.Scale(1.43))

		base := len(*placedPos) + 1
		*placedPos = append(*placedPos, probeO, probeC)
		*placedLabel = append(*placedLabel, "O", "C")
		*bonds = append(*bonds, Bond{base, base + 1, 1})
		placed++
	}
	return placed
}


// linkProbeGroups connects pairs of nearby probe components with sp3 alkane
// linkers.  Each pair becomes a completely independent ProbeSet containing
// copies of both probe atom sets plus the linker chain — so probes can appear
// in multiple output molecules without sharing atom indices.
// Returns one ProbeSet per linked pair.
// buildChainPath finds a clash-free sp3 carbon chain from start to end using
// beam search. startIncoming is the unit vector of the bond arriving AT start
// (direction from start's predecessor toward start). Returns intermediate atom
// positions between start and end (exclusive), or nil if no path is found.
// buildLinkerChain builds an all-sp3 carbon linker joining two probe atoms at p1
// and p2 whose open-valence (exit) directions are dir1 and dir2. It returns the
// full list of chain atom positions INCLUDING the two end-cap atoms
//
//	c0   = p1 + dir1*1.54     (bonds p1)
//	cEnd = p2 + dir2*1.54     (bonds p2)
//
// such that the bond angle is tetrahedral (~109.5°) at p1, at p2, and at every
// chain atom — including the two end junctions, which the previous planar-zigzag
// builder left uncontrolled (the cause of the 30–70° junctions that the geometry
// gate then discarded). Returns nil if no clash-free chain of ≤ maxInter interior
// atoms can close the gap.
//
// Construction: a beam search grows a tetrahedral chain from c0 with the correct
// outgoing tangent (forward·forward = +1/3 ⇒ 109.5° bond angle), and a candidate
// is accepted only when its closing bond to cEnd is tetrahedral at BOTH ends.
func buildLinkerChain(p1, dir1, p2, dir2 Vec3, grid *clashGrid, maxInter int) []Vec3 {
	const (
		bondLen = 1.54
		cosFwd  = 1.0 / 3.0 // forward·forward dot giving a 109.5° bond angle
		sinFwd  = 0.9428
		minBond = 1.30
		maxBond = 1.65
		beamW   = 40
		nRot    = 12
		angLo   = 95.0  // accept junction angles in [95°,130°]: clear of the 82°
		angHi   = 130.0 // geometry floor and centred on the tetrahedral 109.5°
	)
	deg := func(a, b Vec3) float64 {
		na, nb := a.Norm(), b.Norm()
		if na < 1e-9 || nb < 1e-9 {
			return 180
		}
		c := a.Dot(b) / (na * nb)
		if c > 1 {
			c = 1
		} else if c < -1 {
			c = -1
		}
		return 180 / math.Pi * math.Acos(c)
	}
	c0 := p1.Add(dir1.Scale(bondLen))
	cEnd := p2.Add(dir2.Scale(bondLen))
	if !grid.clashFree(c0, vdw("C"), hardTol) || !grid.clashFree(cEnd, vdw("C"), hardTol) {
		return nil
	}
	angOK := func(a float64) bool { return a >= angLo && a <= angHi }

	// Direct cap-to-cap bond (no interior atoms).
	if L := c0.Sub(cEnd).Norm(); L >= minBond && L <= maxBond {
		if angOK(deg(p1.Sub(c0), cEnd.Sub(c0))) && angOK(deg(p2.Sub(cEnd), c0.Sub(cEnd))) {
			return []Vec3{c0, cEnd}
		}
	}

	// closeOK: penultimate atom `pen` (predecessor penPrev) may bond cEnd if that
	// closing bond is tetrahedral at both pen and cEnd.
	closeOK := func(penPrev, pen Vec3) bool {
		if L := pen.Sub(cEnd).Norm(); L < minBond || L > maxBond {
			return false
		}
		return angOK(deg(penPrev.Sub(pen), cEnd.Sub(pen))) &&
			angOK(deg(pen.Sub(cEnd), p2.Sub(cEnd)))
	}

	// Aim the search at the ideal penultimate-atom position, not at cEnd itself: the
	// penultimate atom should sit at a tetrahedral angle off the cEnd→p2 bond so the
	// chain approaches cEnd from the side (a head-on approach gives a ~180° exit
	// junction, which a pure distance-to-cEnd metric would otherwise drive toward —
	// the failure on near-collinear linkers). Two mirror targets cover both sides.
	toP2 := p2.Sub(cEnd)
	if toP2.Norm() < 1e-9 {
		toP2 = dir2.Scale(-1)
	}
	toP2 = toP2.unit()
	nrm := cross3(toP2, cEnd.Sub(c0))
	if nrm.Norm() < 1e-6 {
		nrm = perpendicular(toP2)
	}
	nrm = nrm.unit()
	tang := cross3(nrm, toP2)
	preA := cEnd.Add(toP2.Scale(-1.0 / 3.0).Add(tang.Scale(sinFwd)).Scale(bondLen))
	preB := cEnd.Add(toP2.Scale(-1.0 / 3.0).Sub(tang.Scale(sinFwd)).Scale(bondLen))
	target := func(v Vec3) float64 {
		da, db := v.Sub(preA).Norm(), v.Sub(preB).Norm()
		if da < db {
			return da
		}
		return db
	}

	type node struct {
		chain []Vec3 // c0 .. last placed atom
		dir   Vec3   // forward bond direction into the last atom
		dist  float64
	}
	beam := []node{{[]Vec3{c0}, dir1, target(c0)}}
	for depth := 0; depth < maxInter; depth++ {
		var next []node
		for _, cur := range beam {
			pos := cur.chain[len(cur.chain)-1]
			bp := perpendicular(cur.dir)
			q := cross3(cur.dir, bp).unit()
			for oi := 0; oi < nRot; oi++ {
				ang := float64(oi) * 2 * math.Pi / float64(nRot)
				p := bp.Scale(math.Cos(ang)).Add(q.Scale(math.Sin(ang)))
				d := cur.dir.Scale(cosFwd).Add(p.Scale(sinFwd)).unit()
				np := pos.Add(d.Scale(bondLen))
				if !grid.clashFree(np, vdw("C"), hardTol) {
					continue
				}
				selfClash := false
				for i := 0; i < len(cur.chain)-1; i++ { // skip immediate predecessor
					if np.Sub(cur.chain[i]).Norm() < 2*vdw("C")-1.0 {
						selfClash = true
						break
					}
				}
				if selfClash {
					continue
				}
				if closeOK(pos, np) {
					return append(append([]Vec3{}, cur.chain...), np, cEnd)
				}
				dd := target(np)
				if remaining := maxInter - depth - 1; dd > float64(remaining)*bondLen*1.3+maxBond {
					continue
				}
				nc := append(append([]Vec3{}, cur.chain...), np)
				next = append(next, node{nc, d, dd})
			}
		}
		if len(next) == 0 {
			break
		}
		sort.Slice(next, func(i, j int) bool { return next[i].dist < next[j].dist })
		if len(next) > beamW {
			next = next[:beamW]
		}
		beam = next
	}
	return nil
}

func linkProbeGroups(placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond,
	proteinHeavy []heavyAtom, maxDist float64, linkableN int) []ProbeSet {
	grid := newClashGrid(proteinHeavy)

	positions := *placedPos
	labels := *placedLabel
	n := len(positions)

	// ── Connected components (BFS) ───────────────────────────────────────────
	adj := make([][]int, n)
	for _, b := range *bonds {
		i, j := b.I-1, b.J-1
		adj[i] = append(adj[i], j)
		adj[j] = append(adj[j], i)
	}
	comp := make([]int, n)
	for i := range comp {
		comp[i] = -1
	}
	nComp := 0
	for start := 0; start < n; start++ {
		if comp[start] >= 0 {
			continue
		}
		comp[start] = nComp
		q := []int{start}
		for len(q) > 0 {
			cur := q[0]
			q = q[1:]
			for _, nb := range adj[cur] {
				if comp[nb] < 0 {
					comp[nb] = nComp
					q = append(q, nb)
				}
			}
		}
		nComp++
	}

	// ── Pre-qualify components: must have ≥3 atoms and contain a heteroatom or ring.
	// This excludes alkane chain segments and isolated single atoms.
	compOK := make([]bool, nComp)
	{
		compSize := make([]int, nComp)
		compHet := make([]bool, nComp)
		for i := 0; i < linkableN; i++ {
			c := comp[i]
			compSize[c]++
			switch strings.ToUpper(labels[i]) {
			case "N", "O", "S", "SE", "F", "CL", "BR":
				compHet[c] = true
			}
		}
		compRing := make([]bool, nComp)
		for _, b := range *bonds {
			if b.I-1 < linkableN {
				compRing[comp[b.I-1]] = compRing[comp[b.I-1]] || (b.Order == 4)
			}
		}
		// Also detect non-aromatic rings via DFS per component.
		cadj := make([][]int, n)
		for _, b := range *bonds {
			i, j := b.I-1, b.J-1
			if i < linkableN && j < linkableN {
				cadj[i] = append(cadj[i], j)
				cadj[j] = append(cadj[j], i)
			}
		}
		visited := make([]bool, n)
		var dfs func(node, parent, c int) bool
		dfs = func(node, parent, c int) bool {
			visited[node] = true
			for _, nb := range cadj[node] {
				if comp[nb] != c {
					continue
				}
				if nb == parent {
					continue
				}
				if visited[nb] {
					return true
				}
				if dfs(nb, node, c) {
					return true
				}
			}
			return false
		}
		for i := 0; i < linkableN; i++ {
			if !visited[i] && dfs(i, -1, comp[i]) {
				compRing[comp[i]] = true
			}
		}
		for c := 0; c < nComp; c++ {
			compOK[c] = compSize[c] >= 3 && (compHet[c] || compRing[c])
		}
	}

	// ── Attachment-atom eligibility & preference ─────────────────────────────
	// A linker may only attach to an atom that still has spare valence; otherwise
	// sanitizeValence would later strip a bond to rebalance it — e.g. removing a
	// tetrazole/aromatic ring bond and breaking the ring.  Among eligible atoms
	// we prefer non-ring carbons (the stem/methyl/chain "handles") so ring
	// systems stay intact and H-bonding heteroatoms stay unblocked.
	valence := atomValences(n, *bonds)
	inCore := ringCoreAtoms(adj)
	// attachPenalty returns a preference penalty (lower = better) and whether the
	// atom is eligible at all (has spare valence for one more single bond).
	attachPenalty := func(idx int) (float64, bool) {
		if valence[idx]+1 > maxValence(labels[idx]) {
			return 0, false
		}
		pen := 0.0
		if strings.ToUpper(labels[idx]) != "C" {
			pen += 1.0 // prefer carbon over heteroatom (keep donors/acceptors free)
		}
		if inCore[idx] {
			pen += 1.5 // prefer non-ring atoms (keep rings intact)
		}
		return pen, true
	}

	// ── Build candidate pairs: best eligible atom between every qualifying pair ──
	// For each (c1, c2) pick the atom pair minimising distance + attachment
	// penalty, considering only atoms that have spare valence.
	type cpair struct {
		c1, c2 int
		t1, t2 int     // chosen attachment atoms (one per component)
		d      float64
	}
	var pairs []cpair
	for c1 := 0; c1 < nComp; c1++ {
		if !compOK[c1] {
			continue
		}
		for c2 := c1 + 1; c2 < nComp; c2++ {
			if !compOK[c2] {
				continue
			}
			bestCost := math.MaxFloat64
			bestD := math.MaxFloat64
			bestI, bestJ := -1, -1
			for i := 0; i < linkableN; i++ {
				if comp[i] != c1 {
					continue
				}
				pi, oki := attachPenalty(i)
				if !oki {
					continue
				}
				for j := 0; j < linkableN; j++ {
					if comp[j] != c2 {
						continue
					}
					pj, okj := attachPenalty(j)
					if !okj {
						continue
					}
					d := positions[i].Sub(positions[j]).Norm()
					if d < 2.0 || d > maxDist {
						continue
					}
					if cost := d + pi + pj; cost < bestCost {
						bestCost = cost
						bestD = d
						bestI, bestJ = i, j
					}
				}
			}
			if bestI >= 0 {
				pairs = append(pairs, cpair{c1, c2, bestI, bestJ, bestD})
			}
		}
	}
	sort.Slice(pairs, func(i, j int) bool { return pairs[i].d < pairs[j].d })

	// ── Connect nearest-neighbour pairs independently (parallel) ────────────
	// Each pair of probe components becomes its own independent ProbeSet
	// containing copies of both probe atom sets plus the linker chain.
	// Probes can appear in multiple output molecules (reuse) as long as
	// their tail atom has remaining valence.
	//
	// The per-pair body is embarrassingly parallel: all inputs (positions,
	// labels, comp, adj, bonds, grid) are read-only after this point.
	// Results are collected via a channel; bondCount is checked and updated
	// serially after all goroutines finish.
	type pairOut struct {
		ps         ProbeSet
		a1, a2     int
		c1, c2     int
		d          float64
	}
	outCh := make(chan pairOut, len(pairs))
	nCPU := runtime.NumCPU()
	fmt.Printf("  linkProbeGroups: %d pairs, %d CPUs\n", len(pairs), nCPU)
	sem := make(chan struct{}, nCPU)
	var wg sync.WaitGroup

	for _, p := range pairs {
		p := p
		sem <- struct{}{}
		wg.Add(1)
		go func() {
			defer func() { <-sem; wg.Done() }()

			a1, a2 := p.t1, p.t2

			pos1, pos2 := positions[a1], positions[a2]
		chainDist := pos1.Sub(pos2).Norm()

		// Build an all-trans sp3 zigzag chain between pos1 and pos2.
		axis := pos2.Sub(pos1).unit()

		// Compute the exit direction from pos1 toward the first chain atom.
		// Geometry depends on hybridization of the connecting atom:
		//   sp2 (2 heavy neighbors, any aromatic/double bond): exocyclic direction
		//       = -normalize(sum of unit vectors to neighbors) — coplanar, 120°.
		//   sp3 (1 heavy neighbor): tetrahedral cone at 109.5°, pick closest to axis.
		//   no neighbors: straight along axis.
		firstAtomDir := axis
		nbs := adj[a1]
		isSp2 := false
		for _, b := range *bonds {
			if (b.I-1 == a1 || b.J-1 == a1) && (b.Order == 2 || b.Order == 4) {
				isSp2 = true
				break
			}
		}
		// Choose the open-valence direction (respecting existing bonds for every
		// neighbour count) that best aligns with the chain axis toward a2.
		{
			nbPos := make([]Vec3, len(nbs))
			for k, nb := range nbs {
				nbPos[k] = positions[nb]
			}
			bestDot := -2.0
			for _, dir := range openValenceDirs(pos1, nbPos, isSp2) {
				if dir.Dot(axis) > bestDot {
					bestDot = dir.Dot(axis)
					firstAtomDir = dir
				}
			}
		}

		// Place the virtual first chain atom at the correct tetrahedral position.
		// The chain is then built from this adjusted start toward pos2.
		firstAtomPos := pos1.Add(firstAtomDir.Scale(1.54))

		// Compute the symmetric exit direction into a2 (same sp2/sp3 logic as a1).
		// lastAtomPos is a new chain atom placed one bond before pos2, so the
		// chain ends at lastAtomPos and bonds pa2 at the correct sp2 angle.
		lastAtomDir := axis.Scale(-1) // default: toward pos1
		nbsA2 := adj[a2]
		isSp2A2 := false
		for _, b := range *bonds {
			if (b.I-1 == a2 || b.J-1 == a2) && (b.Order == 2 || b.Order == 4) {
				isSp2A2 = true
				break
			}
		}
		{
			target := axis.Scale(-1) // a2's linker bond points back toward pos1
			nbPos := make([]Vec3, len(nbsA2))
			for k, nb := range nbsA2 {
				nbPos[k] = positions[nb]
			}
			bestDot := -2.0
			for _, dir := range openValenceDirs(pos2, nbPos, isSp2A2) {
				if dir.Dot(target) > bestDot {
					bestDot = dir.Dot(target)
					lastAtomDir = dir
				}
			}
		}
		lastAtomPos := pos2.Add(lastAtomDir.Scale(1.54))

		// Build a saturated linker as a fully tetrahedral chain that meets both
		// probe atoms at ~109.5°, respecting the entry/exit tangents firstAtomDir /
		// lastAtomDir. buildLinkerChain returns the whole chain INCLUDING the two
		// end-cap atoms (== firstAtomPos / lastAtomPos), with every junction angle
		// controlled — fixing the old zigzag's strained end junctions.
		fullChain := buildLinkerChain(pos1, firstAtomDir, pos2, lastAtomDir, grid, 9)
		if fullChain == nil {
			return
		}
		// Region B below re-adds the two end-cap atoms and bonds the interior; the
		// builder guarantees clash-free geometry, so chainClash is all-false.
		chainPts := fullChain[1 : len(fullChain)-1]
		chainClash := make([]bool, len(chainPts))
		nInter := len(chainPts)

		// ── Try unsaturated linkers first (shorter, straighter) ─────────────
		// Only attempt if the chain axis aligns reasonably with both tail atoms'
		// outgoing bond directions (within 60° of collinear) — otherwise the
		// sp3–sp junction angle would be severely distorted.
		axisOK := func() bool {
			if len(adj[a1]) > 0 {
				nb := adj[a1][0]
				outDir := pos1.Sub(positions[nb]).unit() // direction tail is already facing
				if math.Abs(outDir.Dot(axis)) < 0.5 {    // > 60° off-axis
					return false
				}
			}
			if len(adj[a2]) > 0 {
				nb := adj[a2][0]
				outDir := pos2.Sub(positions[nb]).unit()
				if math.Abs(outDir.Dot(axis)) < 0.5 {
					return false
				}
			}
			return true
		}()
		// atoms linearly along axis and check for clashes.  If valid, emit the
		// pair ProbeSet immediately and skip the sp3 zigzag entirely.
		//
		// Motifs tried (from shortest to longest):
		//   alkyne:  pos1 –C≡C– pos2   (linear, 2 atoms, ≈1.20 Å each bond)
		//   alkene:  pos1 –C=C– pos2   (planar, 2 atoms, ≈1.34 Å each bond)
		//   enyne:   pos1 –C=C–C≡C– pos2 (4 atoms)
		//   diyne:   pos1 –C≡C–C≡C– pos2 (4 atoms)
		type unsatMotif struct {
			atoms []string  // element labels of intermediate atoms
			bonds []int     // bond orders for each bond (len = len(atoms)+1)
			blen  []float64 // ideal bond length for each bond
		}
		motifs := []unsatMotif{
			{ // alkyne bridge: sp3–C≡C–sp3 (linear, 4.14 Å)
				atoms: []string{"C", "C"},
				bonds: []int{1, 3, 1},
				blen:  []float64{1.47, 1.20, 1.47},
			},
			{ // diacetylene: sp3–C≡C–C≡C–sp3 (linear, 7.21 Å)
				atoms: []string{"C", "C", "C", "C"},
				bonds: []int{1, 3, 1, 3, 1},
				blen:  []float64{1.47, 1.20, 1.37, 1.20, 1.47},
			},
		}
		tryUnsaturated := func() bool {
			for _, m := range motifs {
				// Compute ideal end-to-end length for this motif.
				idealTotal := 0.0
				for _, bl := range m.blen {
					idealTotal += bl
				}
				// Accept if distance matches within ±15%.
				if math.Abs(chainDist-idealTotal)/idealTotal > 0.15 {
					continue
				}
				// Place intermediate atoms linearly along axis.
				var pts []Vec3
				cur := pos1
				for i, sym := range m.atoms {
					bl := m.blen[i]
					pt := cur.Add(axis.Scale(bl))
					if !grid.clashFree(pt, vdw(sym), hardTol) {
						pts = nil
						break
					}
					pts = append(pts, pt)
					cur = pt
				}
				if len(m.atoms) > 0 && len(pts) < len(m.atoms) {
					continue // some atom clashed
				}
				// Validate final bond to pos2.
				finalBl := m.blen[len(m.blen)-1]
				if math.Abs(cur.Sub(pos2).Norm()-finalBl) > finalBl*0.15 {
					continue
				}
				// Build independent ProbeSet for this pair.
				var ps ProbeSet
				remapC1 := make([]int, n)
				for i, c := range comp {
					if c == p.c1 {
						remapC1[i] = ps.Add(positions[i], labels[i])
					}
				}
				for _, b := range *bonds {
					if comp[b.I-1] == p.c1 && comp[b.J-1] == p.c1 {
						ps.Bond(remapC1[b.I-1], remapC1[b.J-1], b.Order)
					}
				}
				remapC2 := make([]int, n)
				for i, c := range comp {
					if c == p.c2 {
						remapC2[i] = ps.Add(positions[i], labels[i])
					}
				}
				for _, b := range *bonds {
					if comp[b.I-1] == p.c2 && comp[b.J-1] == p.c2 {
						ps.Bond(remapC2[b.I-1], remapC2[b.J-1], b.Order)
					}
				}
				pa1 := remapC1[a1]
				pa2 := remapC2[a2]
				linkerBase := ps.Len() + 1
				for _, pt := range pts {
					ps.Add(pt, "C")
				}
				prev := pa1
				for i := range pts {
					ps.Bond(prev, linkerBase+i, m.bonds[i])
					prev = linkerBase + i
				}
				ps.Bond(prev, pa2, m.bonds[len(m.bonds)-1])
				outCh <- pairOut{ps: ps, a1: a1, a2: a2, c1: p.c1, c2: p.c2, d: p.d}
				return true
			}
			return false
		}
		if axisOK && tryUnsaturated() {
				return
			}
		if nInter == 0 {
			// No intermediate zigzag atoms: connect via the two end-cap atoms only.
			// Chain: pa1 → firstAtomPos → lastAtomPos → pa2.
			if d := firstAtomPos.Sub(lastAtomPos).Norm(); d < 1.3 || d > 1.65 {
				return
			}
			var ps ProbeSet
			remapC1 := make([]int, n)
			for i, c := range comp {
				if c == p.c1 {
					remapC1[i] = ps.Add(positions[i], labels[i])
				}
			}
			for _, b := range *bonds {
				if comp[b.I-1] == p.c1 && comp[b.J-1] == p.c1 {
					ps.Bond(remapC1[b.I-1], remapC1[b.J-1], b.Order)
				}
			}
			remapC2 := make([]int, n)
			for i, c := range comp {
				if c == p.c2 {
					remapC2[i] = ps.Add(positions[i], labels[i])
				}
			}
			for _, b := range *bonds {
				if comp[b.I-1] == p.c2 && comp[b.J-1] == p.c2 {
					ps.Bond(remapC2[b.I-1], remapC2[b.J-1], b.Order)
				}
			}
			pa1 := remapC1[a1]
			pa2 := remapC2[a2]
			fi := ps.Add(firstAtomPos, "C")
			li := ps.Add(lastAtomPos, "C")
			ps.Bond(pa1, fi, 1)
			ps.Bond(fi, li, 1)
			ps.Bond(li, pa2, 1)
			outCh <- pairOut{ps: ps, a1: a1, a2: a2, c1: p.c1, c2: p.c2, d: p.d}
			return
		} else {
			var kept []int
			for k, cl := range chainClash {
				if !cl {
					kept = append(kept, k)
				}
			}

			bondOK := func(p1, p2 Vec3) bool {
				d := p1.Sub(p2).Norm()
				return d >= 1.3 && d <= 1.65
			}
			bad := false
			prevPos := pos1
			for _, k := range kept {
				if !bondOK(prevPos, chainPts[k]) {
					bad = true
					break
				}
				prevPos = chainPts[k]
			}
			if !bad {
				if prevPos.Sub(lastAtomPos).Norm() < 1.3 && len(kept) > 0 {
					kept = kept[:len(kept)-1]
					prevPos = pos1
					for _, k := range kept {
						prevPos = chainPts[k]
					}
				}
				if !bondOK(prevPos, lastAtomPos) {
					bad = true
				}
			}
			if bad {
				kept = make([]int, nInter)
				for k := range kept { kept[k] = k }
			}

			globalIdx := make([]int, len(chainPts))
			placed2 := 0
			// Build into a local ProbeSet: copy probe atoms from both components,
			// then append the linker chain atoms and bond everything together.
			var ps ProbeSet

			// Remap: copy all atoms in component c1.
			remapC1 := make([]int, n)
			for i, c := range comp {
				if c == p.c1 {
					idx := ps.Add(positions[i], labels[i])
					remapC1[i] = idx
				}
			}
			// Copy bonds within component c1.
			for _, b := range *bonds {
				if comp[b.I-1] == p.c1 && comp[b.J-1] == p.c1 {
					ps.Bond(remapC1[b.I-1], remapC1[b.J-1], b.Order)
				}
			}
			// Remap: copy all atoms in component c2.
			remapC2 := make([]int, n)
			for i, c := range comp {
				if c == p.c2 {
					idx := ps.Add(positions[i], labels[i])
					remapC2[i] = idx
				}
			}
			for _, b := range *bonds {
				if comp[b.I-1] == p.c2 && comp[b.J-1] == p.c2 {
					ps.Bond(remapC2[b.I-1], remapC2[b.J-1], b.Order)
				}
			}
			// Now add linker chain atoms and connect a1→firstAtom→chain→a2.
			pa1 := remapC1[a1]
			pa2 := remapC2[a2]
			chainBase := ps.Len() + 1
			// First: add the tetrahedral entry atom at firstAtomPos.
			firstAtomGlobalIdx := ps.Add(firstAtomPos, "C")
			ps.Bond(pa1, firstAtomGlobalIdx, 1)
			// Then: add remaining zigzag chain atoms.
			placed2 = 0
			for _, k := range kept {
				ps.Add(chainPts[k], "C")
				globalIdx[k] = chainBase + 1 + placed2
				placed2++
			}
			prev := firstAtomGlobalIdx
			for _, k := range kept {
				ps.Bond(prev, globalIdx[k], 1)
				prev = globalIdx[k]
			}
			// Add the sp2-aware exit atom at the a2 end, then bond to pa2.
			lastAtomGlobalIdx := ps.Add(lastAtomPos, "C")
			ps.Bond(prev, lastAtomGlobalIdx, 1)
			ps.Bond(lastAtomGlobalIdx, pa2, 1)
			ps.C1 = p.c1
			ps.C2 = p.c2
			outCh <- pairOut{ps: ps, a1: a1, a2: a2, c1: p.c1, c2: p.c2, d: p.d}
		}
		}()
	}

	go func() { wg.Wait(); close(outCh) }()

	// Collect all results, sort by pair distance, then greedily assign each
	// probe component to at most one output molecule (nearest-first).
	var outs []pairOut
	for out := range outCh {
		outs = append(outs, out)
	}
	sort.Slice(outs, func(i, j int) bool { return outs[i].d < outs[j].d })

	usedComp := make(map[int]bool)
	var result []ProbeSet
	for _, out := range outs {
		if usedComp[out.c1] || usedComp[out.c2] {
			continue
		}
		result = append(result, out.ps)
		usedComp[out.c1] = true
		usedComp[out.c2] = true
	}

	return result
}


// amideVariants generates ProbeSet variants where each adjacent pair of sp3
// chain carbons (degree 2, all single bonds) is replaced with an amide group:
// one C becomes C=O (sp2, 1.22 Å) and the adjacent C becomes N.  Both
// orientations of the amide are tried.  Only variants where the carbonyl O
// clears the protein are kept.
func amideVariants(mol ProbeSet, grid *clashGrid) []ProbeSet {
	n := len(mol.Pos)
	degree := make([]int, n)
	adj := make([][]int, n)
	for _, b := range mol.Bonds {
		i, j := b.I-1, b.J-1
		degree[i]++
		degree[j]++
		adj[i] = append(adj[i], j)
		adj[j] = append(adj[j], i)
	}

	// isChainC: sp3 carbon with exactly 2 single bonds.
	isChainC := func(idx int) bool {
		if mol.Labels[idx] != "C" || degree[idx] != 2 {
			return false
		}
		for _, b := range mol.Bonds {
			if (b.I-1 == idx || b.J-1 == idx) && b.Order != 1 {
				return false
			}
		}
		return true
	}

	var variants []ProbeSet
	seen := make(map[[2]int]bool)

	for _, b := range mol.Bonds {
		if b.Order != 1 {
			continue
		}
		i, j := b.I-1, b.J-1
		if !isChainC(i) || !isChainC(j) {
			continue
		}
		// Try both orientations: (i=C=O, j=N) and (j=C=O, i=N).
		for _, pair := range [][2]int{{i, j}, {j, i}} {
			carbonIdx, nitroIdx := pair[0], pair[1]
			key := [2]int{carbonIdx, nitroIdx}
			if seen[key] {
				continue
			}
			seen[key] = true

			cPos := mol.Pos[carbonIdx]
			nPos := mol.Pos[nitroIdx]

			// Predecessor of carbonIdx (its other neighbour).
			predIdx := -1
			for _, nb := range adj[carbonIdx] {
				if nb != nitroIdx {
					predIdx = nb
					break
				}
			}

			// Skip if carbonIdx is bonded to a probe heteroatom — that junction
			// carbon is sp3 and must not become a C=O.
			if predIdx >= 0 && mol.Labels[predIdx] != "C" {
				continue
			}

			// Carbonyl O direction: sp2 bisector at C (120° from both bonds).
			var v1 Vec3
			if predIdx >= 0 {
				v1 = mol.Pos[predIdx].Sub(cPos).unit()
			} else {
				v1 = nPos.Sub(cPos).Scale(-1).unit()
			}
			v2 := nPos.Sub(cPos).unit()
			sumVec := v1.Add(v2)
			if sumVec.Norm() < 1e-6 {
				continue
			}
			oDir := sumVec.Scale(-1).unit()
			oPos := cPos.Add(oDir.Scale(1.22))

			// Reposition N to the sp2 amide plane: third 120° direction at C.
			// v1 (to pred) + oDir (to O) + vN (to N) = 0 for ideal sp2.
			vN := v1.Add(oDir).Scale(-1)
			if vN.Norm() < 1e-6 {
				continue
			}
			newNPos := cPos.Add(vN.unit().Scale(1.33))

			if !grid.clashFree(oPos, vdw("O"), hardTol) {
				continue
			}
			if !grid.clashFree(newNPos, vdw("N"), hardTol) {
				continue
			}
			// Require the carbonyl O or the amide N to be near a protein atom —
			// otherwise the amide points into solvent and adds no H-bond value.
			if grid.countNearby(oPos, 0.0, 3.5) == 0 && grid.countNearby(newNPos, 0.0, 3.5) == 0 {
				continue
			}

			// Build variant: copy mol, place N at sp2 amide position, add O.
			var v ProbeSet
			for k, p := range mol.Pos {
				lbl := mol.Labels[k]
				pos := p
				if k == nitroIdx {
					lbl = "N"
					pos = newNPos
				}
				v.Add(pos, lbl)
			}
			for _, bond := range mol.Bonds {
				v.Bond(bond.I, bond.J, bond.Order)
			}
			oIdx := v.Add(oPos, "O")
			v.Bond(carbonIdx+1, oIdx, 2) // 1-based, double bond
			variants = append(variants, v)
		}
	}
	return variants
}

// linkMoleculeGroups connects each assembled 2-probe molecule to its nearest
// neighbour molecule by their nearest atoms, producing larger merged molecules.
// Each input molecule is used at most once (greedy nearest-first matching).
// maxDist is the maximum inter-molecule atom distance to attempt linking.
func linkMoleculeGroups(mols []ProbeSet, proteinHeavy []heavyAtom, maxDist float64) []ProbeSet {
	grid := newClashGrid(proteinHeavy)
	if len(mols) < 2 {
		return nil
	}

	// Build local adjacency list and sp2 flag for one ProbeSet.
	molAdj := func(ps ProbeSet) [][]int {
		adj := make([][]int, len(ps.Pos))
		for _, b := range ps.Bonds {
			i, j := b.I-1, b.J-1
			adj[i] = append(adj[i], j)
			adj[j] = append(adj[j], i)
		}
		return adj
	}
	molSp2 := func(ps ProbeSet, idx int) bool {
		for _, b := range ps.Bonds {
			if (b.I-1 == idx || b.J-1 == idx) && (b.Order == 2 || b.Order == 4) {
				return true
			}
		}
		return false
	}

	// exitDir computes the sp2-aware exit direction from atom idx in ps.
	// preferDir is a hint for the preferred half-space (e.g. toward the other probe).
	exitDir := func(ps ProbeSet, idx int, adj [][]int, preferDir Vec3) Vec3 {
		pos := ps.Pos[idx]
		nbs := adj[idx]
		if len(nbs) >= 2 && molSp2(ps, idx) {
			sum := Vec3{}
			for _, nb := range nbs {
				sum = sum.Add(ps.Pos[nb].Sub(pos).unit())
			}
			if sum.Norm() > 1e-6 {
				return sum.Scale(-1).unit()
			}
		} else if len(nbs) == 1 {
			incoming := pos.Sub(ps.Pos[nbs[0]]).unit()
			baseP := perpendicular(incoming)
			p2 := cross3(incoming, baseP)
			const cosT, sinT = -0.333, 0.9428
			best := -2.0
			dir := preferDir
			for oi := 0; oi < 12; oi++ {
				ang := float64(oi) * math.Pi / 6
				perp := baseP.Scale(math.Cos(ang)).Add(p2.Scale(math.Sin(ang)))
				d := incoming.Scale(cosT).Add(perp.Scale(sinT)).unit()
				if d.Dot(preferDir) > best {
					best = d.Dot(preferDir)
					dir = d
				}
			}
			return dir
		}
		return preferDir
	}

	// Attachment eligibility: only atoms with spare valence may be linked (so a
	// later bond never over-fills an atom and forces sanitizeValence to strip a
	// ring bond), preferring non-ring carbons. Same policy as linkProbeGroups.
	valences := make([][]int, len(mols))
	cores := make([][]bool, len(mols))
	for i := range mols {
		valences[i] = atomValences(len(mols[i].Pos), mols[i].Bonds)
		cores[i] = ringCoreAtoms(molAdj(mols[i]))
	}
	penalty := func(mi, idx int) (float64, bool) {
		if valences[mi][idx]+1 > maxValence(mols[mi].Labels[idx]) {
			return 0, false
		}
		pen := 0.0
		if strings.ToUpper(mols[mi].Labels[idx]) != "C" {
			pen += 1.0
		}
		if cores[mi][idx] {
			pen += 1.5
		}
		return pen, true
	}

	type candidate struct {
		i, j   int
		a1, a2 int
		d      float64
	}
	var cands []candidate
	for i := 0; i < len(mols); i++ {
		for j := i + 1; j < len(mols); j++ {
			bestCost := math.MaxFloat64
			bestD := math.MaxFloat64
			bestA1, bestA2 := -1, -1
			for a1 := range mols[i].Pos {
				pi, oki := penalty(i, a1)
				if !oki {
					continue
				}
				for a2 := range mols[j].Pos {
					pj, okj := penalty(j, a2)
					if !okj {
						continue
					}
					d := mols[i].Pos[a1].Sub(mols[j].Pos[a2]).Norm()
					if d < 2.0 || d > maxDist {
						continue
					}
					if cost := d + pi + pj; cost < bestCost {
						bestCost = cost
						bestD = d
						bestA1, bestA2 = a1, a2
					}
				}
			}
			if bestA1 >= 0 {
				cands = append(cands, candidate{i, j, bestA1, bestA2, bestD})
			}
		}
	}
	sort.Slice(cands, func(a, b int) bool { return cands[a].d < cands[b].d })

	type pairOut struct {
		ps   ProbeSet
		i, j int
		d    float64
	}
	outCh := make(chan pairOut, len(cands))
	nCPU := runtime.NumCPU()
	sem := make(chan struct{}, nCPU)
	var wg sync.WaitGroup

	for _, c := range cands {
		c := c
		sem <- struct{}{}
		wg.Add(1)
		go func() {
			defer func() { <-sem; wg.Done() }()

			mol1, mol2 := mols[c.i], mols[c.j]
			adj1 := molAdj(mol1)
			adj2 := molAdj(mol2)
			pos1 := mol1.Pos[c.a1]
			pos2 := mol2.Pos[c.a2]
			axis := pos2.Sub(pos1).unit()

			firstAtomDir := exitDir(mol1, c.a1, adj1, axis)
			lastAtomDir := exitDir(mol2, c.a2, adj2, axis.Scale(-1))

			firstAtomPos := pos1.Add(firstAtomDir.Scale(1.54))
			lastAtomPos := pos2.Add(lastAtomDir.Scale(1.54))

			// Tetrahedral linker respecting both probe atoms' exit tangents; the
			// merge-bonding below re-adds the two end-cap atoms and the interior.
			fullChain := buildLinkerChain(pos1, firstAtomDir, pos2, lastAtomDir, grid, 9)
			if fullChain == nil {
				return
			}
			bestPts := fullChain[1 : len(fullChain)-1]
			nInter := len(bestPts)

			// Validate final bond length.
			var actualFinal Vec3
			if len(bestPts) == 0 {
				actualFinal = firstAtomPos
			} else {
				actualFinal = bestPts[len(bestPts)-1]
			}
			if d := actualFinal.Sub(lastAtomPos).Norm(); d < 1.3 || d > 1.65 {
				return
			}

			// Handle nInter == 0: direct bond between end-cap atoms.
			if nInter == 0 {
				if d := firstAtomPos.Sub(lastAtomPos).Norm(); d < 1.3 || d > 1.65 {
					return
				}
			}

			// Merge mol1 + mol2 + linker into one ProbeSet.
			var ps ProbeSet
			n1 := len(mol1.Pos)
			for i, p := range mol1.Pos {
				ps.Add(p, mol1.Labels[i])
			}
			for _, b := range mol1.Bonds {
				ps.Bond(b.I, b.J, b.Order)
			}
			offset := n1
			for i, p := range mol2.Pos {
				ps.Add(p, mol2.Labels[i])
				_ = offset
			}
			for _, b := range mol2.Bonds {
				ps.Bond(b.I+n1, b.J+n1, b.Order)
			}
			pa1 := c.a1 + 1             // 1-based in merged ps
			pa2 := c.a2 + n1 + 1        // 1-based in merged ps

			fiIdx := ps.Add(firstAtomPos, "C")
			ps.Bond(pa1, fiIdx, 1)

			prev := fiIdx
			for _, pt := range bestPts {
				idx := ps.Add(pt, "C")
				ps.Bond(prev, idx, 1)
				prev = idx
			}

			liIdx := ps.Add(lastAtomPos, "C")
			ps.Bond(prev, liIdx, 1)
			ps.Bond(liIdx, pa2, 1)

			outCh <- pairOut{ps: ps, i: c.i, j: c.j, d: c.d}
		}()
	}
	go func() { wg.Wait(); close(outCh) }()

	// Collect, sort, greedy match.
	var outs []pairOut
	for out := range outCh {
		outs = append(outs, out)
	}
	sort.Slice(outs, func(a, b int) bool { return outs[a].d < outs[b].d })

	usedMol := make(map[int]bool)
	var result []ProbeSet
	for _, out := range outs {
		if usedMol[out.i] || usedMol[out.j] {
			continue
		}
		result = append(result, out.ps)
		usedMol[out.i] = true
		usedMol[out.j] = true
	}
	return result
}

// linkTripleGroups builds tri-pharmacophoric molecules by finding pairs of
// pairwise-linked ProbeSets that share exactly one probe component (the hub),
// then merging them into a single A–hub–C molecule.
// Uses the C1/C2 component indices set by linkProbeGroups.
func linkTripleGroups(pairs []ProbeSet, maxDist float64) []ProbeSet {
	var triples []ProbeSet

	// Index pairs by their component sets for quick lookup.
	type compKey [2]int
	key := func(p ProbeSet) compKey {
		if p.C1 <= p.C2 { return compKey{p.C1, p.C2} }
		return compKey{p.C2, p.C1}
	}

	// For each pair of pairs: check if they share exactly one component.
	for i := 0; i < len(pairs); i++ {
		pi := pairs[i]
		if pi.C1 < 0 || pi.C2 < 0 { continue }
		for j := i + 1; j < len(pairs); j++ {
			pj := pairs[j]
			if pj.C1 < 0 || pj.C2 < 0 { continue }
			// Skip if they connect the same component pair.
			if key(pi) == key(pj) { continue }

			// Find the shared component (hub) and the two outer components.
			var hub, cA, cB int
			switch {
			case pi.C1 == pj.C1: hub, cA, cB = pi.C1, pi.C2, pj.C2
			case pi.C1 == pj.C2: hub, cA, cB = pi.C1, pi.C2, pj.C1
			case pi.C2 == pj.C1: hub, cA, cB = pi.C2, pi.C1, pj.C2
			case pi.C2 == pj.C2: hub, cA, cB = pi.C2, pi.C1, pj.C1
			default: continue // no shared component
			}
			_ = hub

			// Ensure A and C are not directly connected to each other
			// (avoid redundant connections).
			alreadyLinked := false
			for k := 0; k < len(pairs); k++ {
				pk := pairs[k]
				if (pk.C1 == cA && pk.C2 == cB) || (pk.C1 == cB && pk.C2 == cA) {
					alreadyLinked = true
					break
				}
			}
			if alreadyLinked { continue }

			merged, ok := mergeTriple(pi, pj)
			if !ok { continue }

			// Reject if total span is unreasonably large.
			maxD := 0.0
			for a := 0; a < len(merged.Pos); a++ {
				for b := a + 1; b < len(merged.Pos); b++ {
					if d := merged.Pos[a].Sub(merged.Pos[b]).Norm(); d > maxD {
						maxD = d
					}
				}
			}
			if maxD > maxDist*1.5 { continue }

			triples = append(triples, merged)
		}
	}
	return triples
}

// mergeTriple attempts to merge two pairwise-linked ProbeSets that share a
// common probe component.  Returns the merged ProbeSet and true on success.
// The two molecules must have a "hub" region — a set of atoms present in both
// ProbeSets at the same positions — with at least 3 distinct atoms so we can
// confirm it's a real probe, not an accidental overlap.
func mergeTriple(a, b ProbeSet) (ProbeSet, bool) {
	const posTol = 0.05 // Å — atoms at same position are considered identical

	// Find which atoms in b correspond to atoms in a (the shared hub).
	// remapB[i] = index in a of atom b.Pos[i], or -1 if not shared.
	remapB := make([]int, len(b.Pos))
	for i := range remapB {
		remapB[i] = -1
	}
	sharedCount := 0
	for bi, bp := range b.Pos {
		for ai, ap := range a.Pos {
			if bp.Sub(ap).Norm() < posTol && b.Labels[bi] == a.Labels[ai] {
				remapB[bi] = ai
				sharedCount++
				break
			}
		}
	}

	// Need enough shared atoms to be a real probe overlap (not coincidence).
	if sharedCount < 3 {
		return ProbeSet{}, false
	}

	// The shared atoms form the "hub" probe.  Atoms in b that are NOT shared
	// are b's unique probe + its linker chain — append them to a.
	merged := ProbeSet{}
	// Copy all of a.
	for k, p := range a.Pos {
		merged.Add(p, a.Labels[k])
	}
	for _, bond := range a.Bonds {
		merged.Bond(bond.I, bond.J, bond.Order)
	}

	// Map from b's atom index to merged atom index.
	bToMerged := make([]int, len(b.Pos))
	for bi := range b.Pos {
		if remapB[bi] >= 0 {
			bToMerged[bi] = remapB[bi] + 1 // 1-based, same position as in a
		} else {
			// New atom — append.
			idx := merged.Add(b.Pos[bi], b.Labels[bi])
			bToMerged[bi] = idx
		}
	}
	// Copy b's bonds, remapped.
	for _, bond := range b.Bonds {
		newI := bToMerged[bond.I-1]
		newJ := bToMerged[bond.J-1]
		// Skip if this bond already exists in merged (shared hub bonds).
		exists := false
		for _, existing := range merged.Bonds {
			if (existing.I == newI && existing.J == newJ) ||
				(existing.I == newJ && existing.J == newI) {
				exists = true
				break
			}
		}
		if !exists {
			merged.Bond(newI, newJ, bond.Order)
		}
	}
	merged.Name = fmt.Sprintf("triple-%da-%db", len(a.Pos), len(b.Pos))
	// Reject pure-carbon triples (no pharmacophoric value).
	hasHet := false
	for _, lbl := range merged.Labels {
		switch strings.ToUpper(lbl) {
		case "N", "O", "S", "F", "CL":
			hasHet = true
		}
	}
	if !hasHet {
		return ProbeSet{}, false
	}
	return merged, true
}
// placed atom lists.  Intra-residue bonds (N–CA single, CA–C single, C=O
// double) and inter-residue peptide bonds (distance-gated at < 2 Å) are added
// so the backbone appears as one connected SDF entry.
// Returns the number of complete residues added.
func addBackbone(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type resKey struct {
		chain string
		seq   int
	}
	type bbRes struct {
		N, CA, C, O Vec3
		found        [4]bool
	}

	res := map[resKey]*bbRes{}
	var order []resKey

	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		name := strings.TrimSpace(strings.ToUpper(a.Name))
		switch name {
		case "N", "CA", "C", "O":
		default:
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
		if res[key] == nil {
			res[key] = &bbRes{}
			order = append(order, key)
		}
		r := res[key]
		switch name {
		case "N":
			r.N, r.found[0] = a.Pos, true
		case "CA":
			r.CA, r.found[1] = a.Pos, true
		case "C":
			r.C, r.found[2] = a.Pos, true
		case "O":
			r.O, r.found[3] = a.Pos, true
		}
	}

	// Sort by chain then residue number (insertion sort).
	for i := 1; i < len(order); i++ {
		for j := i; j > 0; j-- {
			a, b := order[j-1], order[j]
			if a.chain > b.chain || (a.chain == b.chain && a.seq > b.seq) {
				order[j-1], order[j] = order[j], order[j-1]
			} else {
				break
			}
		}
	}

	// Append atoms + intra-residue bonds; record 1-based atom indices.
	type bbIdx struct{ N, CA, C, O int }
	idx := map[resKey]bbIdx{}
	added := 0
	for _, key := range order {
		r := res[key]
		if !r.found[0] || !r.found[1] || !r.found[2] || !r.found[3] {
			continue
		}
		base := len(*placedPos) + 1
		*placedPos = append(*placedPos, r.N, r.CA, r.C, r.O)
		*placedLabel = append(*placedLabel, "N", "C", "C", "O")
		*bonds = append(*bonds,
			Bond{base, base + 1, 1},     // N–CA single
			Bond{base + 1, base + 2, 1}, // CA–C single
			Bond{base + 2, base + 3, 2}, // C=O double
		)
		idx[key] = bbIdx{base, base + 1, base + 2, base + 3}
		added++
	}

	// Peptide bonds: C(i)–N(i+1), detected by distance (< 2 Å).
	for i := 0; i < len(order)-1; i++ {
		cur, nxt := order[i], order[i+1]
		if cur.chain != nxt.chain {
			continue
		}
		ci, ok1 := idx[cur]
		ni, ok2 := idx[nxt]
		if !ok1 || !ok2 {
			continue
		}
		if (*placedPos)[ci.C-1].Sub((*placedPos)[ni.N-1]).Norm() < 2.0 {
			*bonds = append(*bonds, Bond{ci.C, ni.N, 1})
		}
	}

	return added
}
