package main

import (
	"math"
	"strings"
)

// methylCap appends a methyl C (sp3, bondCount=1) bonded to the atom at
// anchorIdx (1-based), placed 1.47 Å away from anchor in the direction
// opposite to awayFrom.  This gives linkProbeGroups a valid tail atom.
func methylCap(placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond,
	anchorIdx int, anchor, awayFrom Vec3) {
	dir := anchor.Sub(awayFrom)
	if dir.Norm() < 1e-6 {
		dir = Vec3{1, 0, 0}
	}
	capPos := anchor.Add(dir.unit().Scale(1.47))
	capIdx := len(*placedPos) + 1
	*placedPos = append(*placedPos, capPos)
	*placedLabel = append(*placedLabel, "C")
	*bonds = append(*bonds, Bond{anchorIdx, capIdx, 1})
}

func placeArgAcids(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type resKey struct {
		chain string
		seq   int
	}
	type guanGroup struct {
		NE, CZ, NH1, NH2 Vec3
		found             [4]bool
	}
	groups := map[resKey]*guanGroup{}
	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "ARG" {
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
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

	// Collect all non-H atoms tagged with their residue for per-ARG exclusion.
	type taggedHeavy struct {
		pos          Vec3
		vdwR         float64
		chain        string
		seq          int
	}
	var heavy []taggedHeavy
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, taggedHeavy{a.Pos, vdw(a.Element), a.ChainID, a.ResSeq})
	}

	placed := 0
	for key, g := range groups {
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

		clears := func(pos Vec3, elem string) bool {
			r := vdw(elem)
			for _, h := range heavy {
				if h.chain == key.chain && h.seq == key.seq {
					continue
				}
				if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
					return false
				}
			}
			return true
		}

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
			if clears(c1, "C") && clears(c2, "C") && clears(o1, "O") && clears(o2, "O") {
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
	type resKey struct {
		chain string
		seq   int
	}
	type carboxylGroup struct {
		prevC, carboxylC, O1, O2 Vec3
		found                    [4]bool
	}
	groups := map[resKey]*carboxylGroup{}

	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		if res != "ASP" && res != "GLU" {
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
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

	type taggedHeavy struct {
		pos   Vec3
		vdwR  float64
		chain string
		seq   int
	}
	var heavy []taggedHeavy
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, taggedHeavy{a.Pos, vdw(a.Element), a.ChainID, a.ResSeq})
	}

	placed := 0
	for key, g := range groups {
		if !g.found[0] || !g.found[1] || !g.found[2] || !g.found[3] {
			continue
		}

		clears := func(pos Vec3, elem string) bool {
			r := vdw(elem)
			for _, h := range heavy {
				if h.chain == key.chain && h.seq == key.seq {
					continue
				}
				if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
					return false
				}
			}
			return true
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

		// buildGuanidine returns the four atom positions for a given zHat and xHat,
		// rotated by angle θ around zHat.
		buildGuanidine := func(zHat, xHat Vec3, theta float64) (Vec3, Vec3, Vec3, Vec3) {
			cosT, sinT := math.Cos(theta), math.Sin(theta)
			// Rotate xHat around zHat by theta.
			yHat := cross3(zHat, xHat)
			rx := xHat.Scale(cosT).Add(yHat.Scale(sinT))
			guanC := g.carboxylC.Add(zHat.Scale(4.15))
			n1 := guanC.Add(rx.Scale(1.152)).Sub(zHat.Scale(0.665))
			n2 := guanC.Sub(rx.Scale(1.152)).Sub(zHat.Scale(0.665))
			n3 := guanC.Add(zHat.Scale(1.33))
			return guanC, n1, n2, n3
		}

		found := false
		var bestC, bestN1, bestN2, bestN3 Vec3
		// Try both faces (±zHat) × 12 rotation steps (30°).
		for _, sign := range []float64{1, -1} {
			if found {
				break
			}
			zHat := baseZ.Scale(sign)
			xHat := buildXHat(zHat)
			for step := 0; step < 12; step++ {
				theta := float64(step) * math.Pi / 6
				guanC, n1, n2, n3 := buildGuanidine(zHat, xHat, theta)
				if clears(guanC, "C") && clears(n1, "N") && clears(n2, "N") && clears(n3, "N") {
					bestC, bestN1, bestN2, bestN3 = guanC, n1, n2, n3
					found = true
					break
				}
			}
		}
		if !found {
			continue
		}

		base := len(*placedPos) + 1
		*placedPos = append(*placedPos, bestC, bestN1, bestN2, bestN3)
		*placedLabel = append(*placedLabel, "C", "N", "N", "N")
		*bonds = append(*bonds,
			Bond{base, base + 1, 1},     // C–N1 single
			Bond{base, base + 2, 1},     // C–N2 single
			Bond{base, base + 3, 2},     // C=N3 double (imine)
		)
		// Methyl cap on N1: gives linkProbeGroups a valid sp3 C tail.
		// N1 has bondCount=1 after the C–N1 bond, leaving 2 slots free;
		// the methyl C (bondCount=1) becomes the tail.
		methyl := bestN1.Add(bestC.Sub(bestN1).Scale(-1).unit().Scale(1.47))
		methylIdx := len(*placedPos) + 1
		*placedPos = append(*placedPos, methyl)
		*placedLabel = append(*placedLabel, "C")
		*bonds = append(*bonds, Bond{base + 1, methylIdx, 1}) // N1–CH3
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
	type resKey struct {
		chain string
		seq   int
	}
	type tyrGroup struct {
		CZ, OH, CE1 Vec3
		found       [3]bool
	}
	groups := map[resKey]*tyrGroup{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "TYR" {
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
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

	type taggedHeavy struct {
		pos   Vec3
		vdwR  float64
		chain string
		seq   int
	}
	var heavy []taggedHeavy
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, taggedHeavy{a.Pos, vdw(a.Element), a.ChainID, a.ResSeq})
	}

	placed := 0
	for key, g := range groups {
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

		clears := func(pos Vec3, elem string) bool {
			r := vdw(elem)
			for _, h := range heavy {
				if h.chain == key.chain && h.seq == key.seq {
					continue
				}
				if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
					return false
				}
			}
			return true
		}

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
			if clears(c.O, "O") && clears(c.C2, "C") && clears(c.N, "N") && clears(c.C1, "C") {
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
	type resKey struct {
		chain string
		seq   int
	}
	type tyrGroup struct {
		ring  [6]Vec3 // CG, CD1, CD2, CE1, CE2, CZ
		CZ    Vec3
		found [7]bool // [6 ring atoms] + CZ duplicate for dPara
	}
	// Atom-name → ring slot index (same 6 as the π-stacking code collects).
	ringSlot := map[string]int{
		"CG": 0, "CD1": 1, "CD2": 2, "CE1": 3, "CE2": 4, "CZ": 5,
	}
	groups := map[resKey]*tyrGroup{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "TYR" {
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
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

	type taggedHeavy struct {
		pos   Vec3
		vdwR  float64
		chain string
		seq   int
	}
	var heavy []taggedHeavy
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, taggedHeavy{a.Pos, vdw(a.Element), a.ChainID, a.ResSeq})
	}

	clears := func(key resKey, pos Vec3, elem string) bool {
		r := vdw(elem)
		for _, h := range heavy {
			if h.chain == key.chain && h.seq == key.seq {
				continue
			}
			if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
				return false
			}
		}
		return true
	}

	placed := 0
	for key, g := range groups {
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
			if !clears(key, phenolCenter, "C") || !clears(key, phenolOH, "O") {
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

			// Methyl cap on the ring C opposite the OH (C3, index 3 from base).
			// Direction: away from OH, so the cap points outward into solvent.
			methylCap(placedPos, placedLabel, bonds,
				base+3, ringPts[3], phenolOH)

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
	type resKey struct {
		chain string
		seq   int
	}
	type amideGroup struct {
		prevC, amideC, O, N Vec3
		found               [4]bool
	}
	groups := map[resKey]*amideGroup{}

	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		if res != "ASN" && res != "GLN" {
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
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

	type taggedHeavy struct {
		pos   Vec3
		vdwR  float64
		chain string
		seq   int
	}
	var heavy []taggedHeavy
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, taggedHeavy{a.Pos, vdw(a.Element), a.ChainID, a.ResSeq})
	}

	placed := 0
	for key, g := range groups {
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

		clears := func(pos Vec3, elem string) bool {
			r := vdw(elem)
			for _, h := range heavy {
				if h.chain == key.chain && h.seq == key.seq {
					continue
				}
				if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
					return false
				}
			}
			return true
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
				if clears(pC, "C") && clears(pO, "O") && clears(pN1, "N") && clears(pN2, "N") {
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
		// Methyl cap on N2 (the outward-facing amino N), pointing away from C.
		methylCap(placedPos, placedLabel, bonds, base+3, bpN2, bpC)
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
	type resKey struct {
		chain string
		seq   int
	}
	type hydroxyl struct {
		CB, OG Vec3
		found  [2]bool
	}
	groups := map[resKey]*hydroxyl{}

	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		if res != "SER" && res != "THR" {
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
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

	type taggedHeavy struct {
		pos   Vec3
		vdwR  float64
		chain string
		seq   int
	}
	var heavy []taggedHeavy
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, taggedHeavy{a.Pos, vdw(a.Element), a.ChainID, a.ResSeq})
	}

	placed := 0
	for key, g := range groups {
		if !g.found[0] || !g.found[1] {
			continue
		}

		// Direction: CB → OG extended (probe is the donor/acceptor partner of the hydroxyl).
		// Try both ±zHat: +zHat captures H-bond donors to the OH, −zHat captures acceptors.
		baseZ := g.OG.Sub(g.CB).unit()

		clears := func(pos Vec3, elem string) bool {
			r := vdw(elem)
			for _, h := range heavy {
				if h.chain == key.chain && h.seq == key.seq {
					continue
				}
				if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
					return false
				}
			}
			return true
		}

		placed2 := false
		for _, sign := range []float64{1, -1} {
			zHat := baseZ.Scale(sign)
			probeO := g.OG.Add(zHat.Scale(2.80))
			probeC := g.OG.Add(zHat.Scale(2.80 + 1.43))
			if clears(probeO, "O") && clears(probeC, "C") {
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
	type resKey struct {
		chain string
		seq   int
	}
	type lys struct {
		CE, NZ Vec3
		found  [2]bool
	}
	groups := map[resKey]*lys{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "LYS" {
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
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

	type taggedHeavy struct {
		pos   Vec3
		vdwR  float64
		chain string
		seq   int
	}
	var heavy []taggedHeavy
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, taggedHeavy{a.Pos, vdw(a.Element), a.ChainID, a.ResSeq})
	}

	placed := 0
	for key, g := range groups {
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

		clears := func(pos Vec3, elem string) bool {
			r := vdw(elem)
			for _, h := range heavy {
				if h.chain == key.chain && h.seq == key.seq {
					continue
				}
				if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
					return false
				}
			}
			return true
		}

		found := false
		var bc1, bc2, bo1, bo2 Vec3
		for step := 0; step < 12 && !found; step++ {
			theta := float64(step) * math.Pi / 6
			c1, c2, o1, o2 := buildAcetate(baseX, theta)
			if clears(c1, "C") && clears(c2, "C") && clears(o1, "O") && clears(o2, "O") {
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


// placeCysMethanethiols places a methanethiol probe (C–SH) near each CYS SG,
// capturing the thiol donor/acceptor pharmacophore.
//
// Geometry: the probe S is placed 3.6 Å from residue SG along the CB→SG axis
// (S···S or S···N distance at outer H-bond range). Probe C is 1.82 Å further.
// Returns the number placed.
func placeCysMethanethiols(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type resKey struct {
		chain string
		seq   int
	}
	type thiol struct {
		CB, SG Vec3
		found  [2]bool
	}
	groups := map[resKey]*thiol{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "CYS" {
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
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

	type taggedHeavy struct {
		pos   Vec3
		vdwR  float64
		chain string
		seq   int
	}
	var heavy []taggedHeavy
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, taggedHeavy{a.Pos, vdw(a.Element), a.ChainID, a.ResSeq})
	}

	placed := 0
	for key, g := range groups {
		if !g.found[0] || !g.found[1] {
			continue
		}

		baseZ := g.SG.Sub(g.CB).unit()

		clears := func(pos Vec3, elem string) bool {
			r := vdw(elem)
			for _, h := range heavy {
				if h.chain == key.chain && h.seq == key.seq {
					continue
				}
				if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
					return false
				}
			}
			return true
		}

		placed2 := false
		for _, sign := range []float64{1, -1} {
			zHat := baseZ.Scale(sign)
			probeS := g.SG.Add(zHat.Scale(3.60))
			probeC := g.SG.Add(zHat.Scale(3.60 + 1.82))
			if clears(probeS, "S") && clears(probeC, "C") {
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
	type resKey struct {
		chain string
		seq   int
	}
	type hisGroup struct {
		ring  [5]Vec3 // CG, ND1, CD2, CE1, NE2
		found [5]bool
	}
	ringSlot := map[string]int{
		"CG": 0, "ND1": 1, "CD2": 2, "CE1": 3, "NE2": 4,
	}
	groups := map[resKey]*hisGroup{}

	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "HIS" {
			continue
		}
		key := resKey{a.ChainID, a.ResSeq}
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

	type taggedHeavy struct {
		pos   Vec3
		vdwR  float64
		chain string
		seq   int
	}
	var heavy []taggedHeavy
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, taggedHeavy{a.Pos, vdw(a.Element), a.ChainID, a.ResSeq})
	}

	clears := func(key resKey, pos Vec3, elem string) bool {
		r := vdw(elem)
		for _, h := range heavy {
			if h.chain == key.chain && h.seq == key.seq {
				continue
			}
			if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
				return false
			}
		}
		return true
	}

	placed := 0
	for key, g := range groups {
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
			if !clears(key, centre, "C") {
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
			// Methyl cap on C4 (atom index 3, a ring C), pointing away from ring centre.
			methylCap(placedPos, placedLabel, bonds, base+3, ringPts[3], centre)
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


// linkProbeGroups connects distinct probe components with straight carbon-chain
// linkers.  For each component the best tail atom is chosen — preferring C,
// then fewest bonds, then furthest from any polar atom in the same component.
// Only atoms with 0-based index < linkableN are eligible as tail candidates;
// pharmacophoric probes placed after the hydrophobic/aromatic pass are passed
// in the full slice so the BFS sees them, but they are never picked as linker
// attachment points and never merged into a larger fragment.
// Pairs of tail atoms within maxDist are linked in
// minimum-spanning-tree order (closest first, Union-Find to avoid cycles).
// Chain atoms are spaced at 1.54 Å intervals along the straight line; any
// candidate chain atom that clashes with a protein heavy atom is rejected.
// Returns the number of chains added.
func linkProbeGroups(placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond,
	proteinHeavy []heavyAtom, maxDist float64, linkableN int) int {

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

	// ── Per-atom valence count (to enforce element valence limits) ───────────
	// Use bond order as the valence weight; aromatic bonds (order 4) count as 1.
	valenceOf := func(order int) int {
		if order == 4 {
			return 1
		}
		return order
	}
	bondCount := make([]int, len(positions))
	for _, b := range *bonds {
		w := valenceOf(b.Order)
		bondCount[b.I-1] += w
		bondCount[b.J-1] += w
	}

	// ── Tail atom per component ──────────────────────────────────────────────
	// Tail = atom with available valence; prefer C over heteroatom, then
	// fewest bonds, then maximum clearance from the nearest protein heavy atom
	// (so the ghost methyl beats the probe C: the chain grows into solvent).
	maxValence := func(lbl string) int {
		switch strings.ToUpper(lbl) {
		case "N":
			return 3
		case "O", "S":
			return 2
		default:
			return 4
		}
	}
	// Minimum distance from atom idx to any protein heavy atom.
	proteinClearance := func(idx int) float64 {
		best := 999.0
		for _, h := range proteinHeavy {
			if d := positions[idx].Sub(h.pos).Norm(); d < best {
				best = d
			}
		}
		return best
	}
	tailIdx := make([]int, nComp)
	for i := range tailIdx {
		tailIdx[i] = -1
	}
	for i := 0; i < n; i++ {
		if i >= linkableN {
			continue
		}
		// Only sp3 carbon atoms may serve as linker tails.
		// Heteroatoms (N, O, S), aromatic carbons (bond order 4), and
		// sp2 carbons (any double bond, order 2) all produce wrong bond geometry.
		if labels[i] != "C" {
			continue
		}
		isSp3 := true
		for _, b := range *bonds {
			if (b.I-1 == i || b.J-1 == i) && (b.Order == 4 || b.Order == 2) {
				isSp3 = false
				break
			}
		}
		if !isSp3 {
			continue
		}
		if bondCount[i] >= maxValence(labels[i]) {
			continue
		}
		c := comp[i]
		t := tailIdx[c]
		if t < 0 {
			tailIdx[c] = i
			continue
		}
		iIsC := labels[i] == "C"
		tIsC := labels[t] == "C"
		iClear := proteinClearance(i)
		tClear := proteinClearance(t)
		switch {
		case iIsC && !tIsC:
			tailIdx[c] = i
		case !iIsC && tIsC:
			// keep t
		case bondCount[i] < bondCount[t]:
			tailIdx[c] = i
		case bondCount[i] == bondCount[t] && iClear > tClear:
			// prefer the atom further from protein — ghost methyl beats probe C
			tailIdx[c] = i
		}
	}

	// ── Build candidate pairs sorted by probe-C to probe-C distance ─────────
	// Use the probe C (neighbour of the ghost methyl tail) as the reference
	// position — this is the actual pocket surface atom and gives the true
	// intra-pocket distance.
	type cpair struct{ c1, c2 int; d float64 }
	var pairs []cpair
	// Build candidate pairs sorted by tail-atom-to-tail-atom distance.
	for c1 := 0; c1 < nComp; c1++ {
		if tailIdx[c1] < 0 {
			continue
		}
		for c2 := c1 + 1; c2 < nComp; c2++ {
			if tailIdx[c2] < 0 {
				continue
			}
			d := positions[tailIdx[c1]].Sub(positions[tailIdx[c2]]).Norm()
			if d >= 2.0 && d <= maxDist {
				pairs = append(pairs, cpair{c1, c2, d})
			}
		}
	}
	for i := 1; i < len(pairs); i++ {
		for j := i; j > 0 && pairs[j].d < pairs[j-1].d; j-- {
			pairs[j], pairs[j-1] = pairs[j-1], pairs[j]
		}
	}

	// ── Union-Find (MST: closest pairs first, no cycles) ────────────────────
	// compSize tracks how many original probe components have been merged into
	// each super-component.  A merge is refused when it would push either
	// super-component past maxGroupsPerMol original probes.
	const maxGroupsPerMol = 3
	parent := make([]int, nComp)
	compSize := make([]int, nComp) // number of original components in this super-comp
	tailUsed := make([]bool, nComp) // true once a component's tail has accepted a linker
	for i := range parent {
		parent[i] = i
		compSize[i] = 1
	}
	var find func(int) int
	find = func(x int) int {
		if parent[x] != x {
			parent[x] = find(parent[x])
		}
		return parent[x]
	}

	linked := 0
	for _, p := range pairs {
		r1, r2 := find(p.c1), find(p.c2)
		if r1 == r2 {
			continue // already in the same component
		}
		// Refuse if merging would exceed the per-molecule group limit.
		if compSize[r1]+compSize[r2] > maxGroupsPerMol {
			continue
		}
		// Each original component's tail may accept exactly one linker bond.
		if tailUsed[p.c1] || tailUsed[p.c2] {
			continue
		}
		a1, a2 := tailIdx[p.c1], tailIdx[p.c2]

		// Refuse if either tail atom is already at valence limit.
		if bondCount[a1]+1 > maxValence(labels[a1]) || bondCount[a2]+1 > maxValence(labels[a2]) {
			continue
		}

		// Build the chain directly between the two tail atoms (methyl caps).
		// Each tail is an sp3 C pointing outward from its probe; connecting
		// tail-to-tail keeps the linker in the space between probes.
		pos1, pos2 := positions[a1], positions[a2]
		chainDist := pos1.Sub(pos2).Norm()

		// Build an all-trans sp3 zigzag chain between pos1 and pos2.
		// Bond angle: 109.5° → forward pitch cosA=cos(35.25°)=0.8165, side sinA=0.578.
		// We sample 12 perpendicular orientations (30° steps) and pick the one
		// with the fewest clashes.  Clashing intermediate atoms are dropped from
		// the placed set (bridged over with a longer bond) so the chain is always
		// connected; we only reject the whole link if every orientation leaves
		// every intermediate atom inside a protein atom.
		axis := pos2.Sub(pos1).unit()

		const (
			bondLen = 1.54   // Å, standard C–C
			cosA    = 0.8165 // cos(35.25°) — forward component per bond
			sinA    = 0.578  // sin(35.25°) — perpendicular component per bond
		)
		nBase := int(math.Round(chainDist/1.258)) - 1
		if nBase < 0 {
			nBase = 0
		}
		nInter := nBase
		if nInter == 0 && chainDist > 1.60 {
			nInter = 1
		}

		basePerp := func() Vec3 {
			if len(adj[a1]) > 0 {
				nb := adj[a1][0]
				// Vector from the tail's predecessor to the tail (= existing bond direction).
				// Project out the chain-axis component to get the in-plane offset.
				// Negate so chain[0] is placed ANTI to this bond (180° junction torsion).
				away := positions[a1].Sub(positions[nb])
				away = away.Sub(axis.Scale(away.Dot(axis)))
				if away.Norm() > 0.1 {
					return away.Scale(-1).unit() // anti position
				}
			}
			return perpendicular(axis)
		}()
		perp2 := cross3(axis, basePerp).unit()

		// Build a scaled 2D all-trans sp3 chain in the (axis, perp) plane.
		// The 2D frame is rotated so the ideal chain endpoint maps exactly to
		// pos2, giving perfect 180° torsions and ~109.5° bond angles as long as
		// the scale factor s = d/idealLen is close to 1.
		buildZigzag := func(n int, perpV Vec3) ([]Vec3, float64) {
			if n == 0 {
				return nil, pos1.Sub(pos2).Norm()
			}
			type pt2 struct{ u, v float64 }
			local := make([]pt2, n+2)
			for k := 1; k <= n+1; k++ {
				sign := 1.0
				if (k-1)%2 == 1 {
					sign = -1
				}
				local[k] = pt2{local[k-1].u + cosA*bondLen, local[k-1].v + sign*sinA*bondLen}
			}
			eu, ev := local[n+1].u, local[n+1].v
			idealLen := math.Sqrt(eu*eu + ev*ev)
			d := chainDist
			s := d / idealLen
			cosθ := eu / idealLen
			sinθ := ev / idealLen
			pts := make([]Vec3, n)
			for k := 1; k <= n; k++ {
				u, v := local[k].u, local[k].v
				ru := (u*cosθ + v*sinθ) * s
				rv := (-u*sinθ + v*cosθ) * s
				pts[k-1] = pos1.Add(axis.Scale(ru)).Add(perpV.Scale(rv))
			}
			last := pts[n-1]
			return pts, last.Sub(pos2).Norm()
		}

		// estNInter returns the distance from the last chain atom to pos2.
		estNInter := func(n int, perpV Vec3) float64 {
			_, d := buildZigzag(n, perpV)
			return d
		}
		// sRatio returns the scale factor s=d/idealLen for n intermediate atoms.
		// Outside [0.85, 1.15] the bond lengths and angles become unacceptable.
		sRatio := func(n int) float64 {
			type pt2 struct{ u, v float64 }
			local := make([]pt2, n+2)
			for k := 1; k <= n+1; k++ {
				sign := 1.0
				if (k-1)%2 == 1 {
					sign = -1
				}
				local[k] = pt2{local[k-1].u + cosA*bondLen, local[k-1].v + sign*sinA*bondLen}
			}
			eu, ev := local[n+1].u, local[n+1].v
			idealLen := math.Sqrt(eu*eu + ev*ev)
			return chainDist / idealLen
		}

		// Pre-check and build: search nBase±3 × 12 orientations.
		anyValid := false
		for tryN := nInter - 3; tryN <= nInter+3 && !anyValid; tryN++ {
			if tryN < 0 {
				continue
			}
			s := sRatio(tryN)
			if s < 0.85 || s > 1.15 {
				continue
			}
			for oi := 0; oi < 36 && !anyValid; oi++ {
				ang := float64(oi) * math.Pi / 18
				tp := basePerp.Scale(math.Cos(ang)).Add(perp2.Scale(math.Sin(ang)))
				d := estNInter(tryN, tp)
				if d >= 1.3 && d <= 1.65 {
					anyValid = true
					nInter = tryN
				}
			}
		}
		if !anyValid {
			continue
		}

		clashPts := func(pts []Vec3) int {
			n := 0
			for _, pt := range pts {
				for _, h := range proteinHeavy {
					if pt.Sub(h.pos).Norm() < h.vdwR+vdw("C")-hardTol {
						n++
						break
					}
				}
			}
			return n
		}

		// Compute junction torsion deviations from 180° at both chain endpoints.
		junctionTorsionDev := func(pts []Vec3) float64 {
			torsionDeg := func(p1, p2, p3, p4 Vec3) float64 {
				b1 := p2.Sub(p1)
				b2 := p3.Sub(p2)
				b3 := p4.Sub(p3)
				n1 := cross3(b1, b2)
				n2 := cross3(b2, b3)
				if n1.Norm() < 1e-6 || n2.Norm() < 1e-6 {
					return 0
				}
				c := n1.Dot(n2) / (n1.Norm() * n2.Norm())
				if c > 1 {
					c = 1
				}
				if c < -1 {
					c = -1
				}
				t := 180 / math.Pi * math.Acos(c)
				if cross3(n1, n2).Dot(b2) < 0 {
					t = -t
				}
				return t
			}
			maxDev := 0.0
			// a1 junction: predecessor(a1) → a1 → pts[0] → pts[1]
			if len(pts) >= 2 && len(adj[a1]) > 0 {
				pred := positions[adj[a1][0]]
				t := torsionDeg(pred, pos1, pts[0], pts[1])
				dev := math.Abs(math.Abs(t) - 180)
				if dev > maxDev {
					maxDev = dev
				}
			}
			// a2 junction: pts[-2] → pts[-1] → a2 → successor(a2)
			n := len(pts)
			if n >= 2 && len(adj[a2]) > 0 {
				t := torsionDeg(pts[n-2], pts[n-1], pos2, positions[adj[a2][0]])
				dev := math.Abs(math.Abs(t) - 180)
				if dev > maxDev {
					maxDev = dev
				}
			}
			return maxDev
		}

		bestPts := []Vec3(nil)
		bestClash := math.MaxInt32
		bestFinalDev := math.MaxFloat64
		bestTorDev := math.MaxFloat64
		bestNInter := nInter

		for _, tryN := range []int{nInter - 3, nInter - 2, nInter - 1, nInter, nInter + 1, nInter + 2, nInter + 3} {
			if tryN < 0 {
				continue
			}
			// Skip if the scale factor would distort bond geometry too much.
			s := sRatio(tryN)
			if s < 0.85 || s > 1.15 {
				continue
			}
			for oi := 0; oi < 36; oi++ {
				ang := float64(oi) * math.Pi / 18 // 10° steps
				perp := basePerp.Scale(math.Cos(ang)).Add(perp2.Scale(math.Sin(ang)))
				pts, finalDist := buildZigzag(tryN, perp)
				if finalDist < 1.3 || finalDist > 1.65 {
					continue
				}
				nc := clashPts(pts)
				dev := math.Abs(finalDist - bondLen)
				torDev := junctionTorsionDev(pts)
				// Priority: 1) fewest clashes, 2) best junction torsions, 3) best bond length
				if nc < bestClash ||
					(nc == bestClash && torDev < bestTorDev) ||
					(nc == bestClash && torDev == bestTorDev && dev < bestFinalDev) {
					bestClash = nc
					bestFinalDev = dev
					bestTorDev = torDev
					bestPts = pts
					bestNInter = tryN
				}
			}
		}
		if bestPts == nil {
			continue // no valid geometry found (pre-check should have caught this)
		}
		nInter = bestNInter

		// Final validation: confirm the actual endpoint bond is valid.
		var actualFinalPos Vec3
		if len(bestPts) == 0 {
			actualFinalPos = pos1
		} else {
			actualFinalPos = bestPts[len(bestPts)-1]
		}
		finalBondLen := actualFinalPos.Sub(pos2).Norm()
		if finalBondLen < 1.3 || finalBondLen > 1.65 {
			continue
		}

		// Build clash mask for best orientation.
		chainPts := bestPts
		chainClash := make([]bool, len(chainPts))
		for k, pt := range chainPts {
			for _, h := range proteinHeavy {
				if pt.Sub(h.pos).Norm() < h.vdwR+vdw("C")-hardTol {
					chainClash[k] = true
					break
				}
			}
		}

		// Count clashing intermediate atoms.
		nClash := 0
		for _, cl := range chainClash {
			if cl {
				nClash++
			}
		}
		// Reject only if every intermediate atom clashes (completely buried path).
		if nInter > 0 && nClash == nInter {
			continue
		}

		if nInter == 0 {
			*bonds = append(*bonds, Bond{a1 + 1, a2 + 1, 1})
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
				if prevPos.Sub(pos2).Norm() < 1.3 && len(kept) > 0 {
					kept = kept[:len(kept)-1]
					prevPos = pos1
					for _, k := range kept {
						prevPos = chainPts[k]
					}
				}
				if !bondOK(prevPos, pos2) {
					bad = true
				}
			}
			if bad {
				kept = make([]int, nInter)
				for k := range kept {
					kept[k] = k
				}
			}

			globalIdx := make([]int, len(chainPts))
			chainBase := len(positions) + 1
			placed2 := 0
			for _, k := range kept {
				positions = append(positions, chainPts[k])
				labels = append(labels, "C")
				globalIdx[k] = chainBase + placed2
				placed2++
			}
			prev := a1 + 1
			for _, k := range kept {
				*bonds = append(*bonds, Bond{prev, globalIdx[k], 1})
				prev = globalIdx[k]
			}
			*bonds = append(*bonds, Bond{prev, a2 + 1, 1})
		}

		bondCount[a1]++
		bondCount[a2]++
		// Pin both tails to prevent a second linker bond.
		if bondCount[a1] < maxValence(labels[a1]) {
			bondCount[a1] = maxValence(labels[a1])
		}
		if bondCount[a2] < maxValence(labels[a2]) {
			bondCount[a2] = maxValence(labels[a2])
		}
		// Mark both tails as used — each original component's tail accepts
		// exactly one linker bond.
		tailUsed[p.c1] = true
		tailUsed[p.c2] = true
		// Union: merge smaller into larger (by compSize), update size.
		if compSize[r1] < compSize[r2] {
			r1, r2 = r2, r1
		}
		parent[r2] = r1
		compSize[r1] += compSize[r2]
		linked++
	}

	*placedPos = positions
	*placedLabel = labels
	return linked
}


// addBackbone appends every residue's backbone atoms (N, CA, C, O) to the
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
