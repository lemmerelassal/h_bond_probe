package main

import (
	"math"
	"sort"
	"strings"
)

// This file holds the second wave of pharmacophoric placers that fill coverage
// gaps in the original set:
//   • TRP indole-NH and HIS ring-N H-bond acceptor partners
//   • carboxylate salt-bridge toward a protonated HIS imidazolium
//   • amine partner for ASP/GLU carboxylates (symmetry with the guanidine)
//   • MET S–arene, aromatic→cation and cation→aromatic cation–π pairs
//   • hydroxyl / methoxy probes at H-bond hotspots and at backbone carbonyls
//
// All of them reuse the shared scaffolding from placer_scaffold.go
// (residueKey, newResClashSet, sortedResidueKeys).

// ── Shared local helpers ──────────────────────────────────────────────────────

// addCarbonylAcceptor places a ketone-like C=O whose oxygen accepts an H-bond
// from a donor nitrogen (donorN), whose in-plane N–H points outward along the
// external bisector of its two ring neighbours (nbr1, nbr2).  normal is the ring
// plane normal, used to place the methyl in-plane at 120° from the carbonyl.
// Atoms added: O (carbonyl), C (carbonyl C), C (methyl).
func addCarbonylAcceptor(placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond,
	clash *resClashSet, key residueKey, donorN, nbr1, nbr2, normal Vec3) bool {
	dH := donorN.Sub(nbr1).unit().Add(donorN.Sub(nbr2).unit())
	if dH.Norm() < 1e-6 {
		return false
	}
	dH = dH.unit()
	dP := cross3(normal, dH)
	if dP.Norm() < 1e-6 {
		dP = perpendicular(dH)
	}
	dP = dP.unit()

	o := donorN.Add(dH.Scale(2.9))
	c := donorN.Add(dH.Scale(2.9 + 1.22))
	for _, s := range []float64{1, -1} {
		// Methyl at 120° from the C=O bond (C→O = −dH), in the ring plane.
		me := c.Add(dH.Scale(0.5 * 1.50)).Add(dP.Scale(s * 0.866 * 1.50))
		if clash.clears(o, "O", key) && clash.clears(c, "C", key) && clash.clears(me, "C", key) {
			base := len(*placedPos) + 1
			*placedPos = append(*placedPos, o, c, me)
			*placedLabel = append(*placedLabel, "O", "C", "C")
			*bonds = append(*bonds, Bond{base, base + 1, 2}, Bond{base + 1, base + 2, 1})
			return true
		}
	}
	return false
}

// addBenzeneProbe places a benzene ring (centre, normal) plus an in-plane methyl
// cap (linker handle), clash-checking every atom.  Returns false if any clashes.
func addBenzeneProbe(placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond,
	clash *resClashSet, key residueKey, centre, normal Vec3) bool {
	pts, edges := benzeneRing(centre, normal)
	for _, p := range pts {
		if !clash.clears(p, "C", key) {
			return false
		}
	}
	capDir := pts[3].Sub(centre).unit()
	capPos := pts[3].Add(capDir.Scale(1.54))
	if !clash.clears(capPos, "C", key) {
		return false
	}
	base := len(*placedPos) + 1
	for _, p := range pts {
		*placedPos = append(*placedPos, p)
		*placedLabel = append(*placedLabel, "C")
	}
	for _, e := range edges {
		*bonds = append(*bonds, Bond{base + e[0], base + e[1], 4})
	}
	*placedPos = append(*placedPos, capPos)
	*placedLabel = append(*placedLabel, "C")
	*bonds = append(*bonds, Bond{base + 3, len(*placedPos), 1})
	return true
}

// addAmmoniumProbe places a methylammonium (N + capping C) with N at nPos and
// the carbon at cPos.  Returns false on clash.
func addAmmoniumProbe(placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond,
	clash *resClashSet, key residueKey, nPos, cPos Vec3) bool {
	if !clash.clears(nPos, "N", key) || !clash.clears(cPos, "C", key) {
		return false
	}
	base := len(*placedPos) + 1
	*placedPos = append(*placedPos, nPos, cPos)
	*placedLabel = append(*placedLabel, "N", "C")
	*bonds = append(*bonds, Bond{base, base + 1, 1})
	return true
}

// ── TRP indole-NH acceptor ─────────────────────────────────────────────────────

// placeTrpCarbonyls places a carbonyl acceptor pointed at each TRP indole NH
// (NE1), the strongest directional donor the original set ignored.
func placeTrpCarbonyls(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type trp struct {
		CD1, NE1, CE2 Vec3
		found         [3]bool
	}
	groups := map[residueKey]*trp{}
	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "TRP" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &trp{}
		}
		g := groups[key]
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "CD1":
			g.CD1, g.found[0] = a.Pos, true
		case "NE1":
			g.NE1, g.found[1] = a.Pos, true
		case "CE2":
			g.CE2, g.found[2] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)
	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] || !g.found[2] {
			continue
		}
		normal := ringPlaneNormal([]Vec3{g.CD1, g.NE1, g.CE2})
		if addCarbonylAcceptor(placedPos, placedLabel, bonds, clash, key, g.NE1, g.CD1, g.CE2, normal) {
			placed++
		}
	}
	return placed
}

// ── HIS H-bond partner + protonated salt bridge ────────────────────────────────

// hisRing collects the five imidazole atoms of every HIS.
type hisRing struct {
	CG, ND1, CD2, CE1, NE2 Vec3
	found                  [5]bool
}

func collectHisRings(atoms []Atom) map[residueKey]*hisRing {
	groups := map[residueKey]*hisRing{}
	for _, a := range atoms {
		if strings.TrimSpace(strings.ToUpper(a.ResName)) != "HIS" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &hisRing{}
		}
		g := groups[key]
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "CG":
			g.CG, g.found[0] = a.Pos, true
		case "ND1":
			g.ND1, g.found[1] = a.Pos, true
		case "CD2":
			g.CD2, g.found[2] = a.Pos, true
		case "CE1":
			g.CE1, g.found[3] = a.Pos, true
		case "NE2":
			g.NE2, g.found[4] = a.Pos, true
		}
	}
	return groups
}

func (g *hisRing) complete() bool {
	return g.found[0] && g.found[1] && g.found[2] && g.found[3] && g.found[4]
}

func (g *hisRing) centroidNormal() (Vec3, Vec3) {
	pts := []Vec3{g.CG, g.ND1, g.CD2, g.CE1, g.NE2}
	c := Vec3{}
	for _, p := range pts {
		c = c.Add(p)
	}
	return c.Scale(1.0 / 5), ringPlaneNormal(pts)
}

// placeHisCarbonyls places a carbonyl acceptor toward each ring nitrogen (ND1
// and NE2) — whichever bears the NH in the actual tautomer will accept.
func placeHisCarbonyls(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	groups := collectHisRings(atoms)
	clash := newResClashSet(atoms)
	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.complete() {
			continue
		}
		_, normal := g.centroidNormal()
		// ND1 neighbours: CG, CE1.  NE2 neighbours: CD2, CE1.
		if addCarbonylAcceptor(placedPos, placedLabel, bonds, clash, key, g.ND1, g.CG, g.CE1, normal) {
			placed++
		}
		if addCarbonylAcceptor(placedPos, placedLabel, bonds, clash, key, g.NE2, g.CD2, g.CE1, normal) {
			placed++
		}
	}
	return placed
}

// placeHisAcids places a carboxylate that makes a bidentate salt bridge to a
// protonated (cationic) HIS imidazolium, with O1/O2 reaching toward ND1-H/NE2-H.
func placeHisAcids(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	groups := collectHisRings(atoms)
	clash := newResClashSet(atoms)
	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.complete() {
			continue
		}
		centroid, normal := g.centroidNormal()
		nMid := g.ND1.Add(g.NE2).Scale(0.5)
		zHat := nMid.Sub(centroid)
		zHat = zHat.Sub(normal.Scale(zHat.Dot(normal))) // project into ring plane
		if zHat.Norm() < 1e-6 {
			continue
		}
		zHat = zHat.unit()
		perp := g.NE2.Sub(g.ND1)
		perp = perp.Sub(normal.Scale(perp.Dot(normal)))
		if perp.Norm() < 1e-6 {
			continue
		}
		perp = perp.unit()

		o1 := nMid.Add(perp.Scale(1.10)).Add(zHat.Scale(2.80)) // toward NE2-H
		o2 := nMid.Sub(perp.Scale(1.10)).Add(zHat.Scale(2.80)) // toward ND1-H
		c2 := nMid.Add(zHat.Scale(3.40))                       // carboxyl C
		c1 := c2.Add(zHat.Scale(1.52))                         // methyl
		if clash.clears(o1, "O", key) && clash.clears(o2, "O", key) &&
			clash.clears(c2, "C", key) && clash.clears(c1, "C", key) {
			base := len(*placedPos) + 1
			*placedPos = append(*placedPos, c1, c2, o1, o2)
			*placedLabel = append(*placedLabel, "C", "C", "O", "O")
			*bonds = append(*bonds,
				Bond{base, base + 1, 1},     // C1–C2 single
				Bond{base + 1, base + 2, 2}, // C2=O1 double
				Bond{base + 1, base + 3, 1}, // C2–O2 single
			)
			placed++
		}
	}
	return placed
}

// ── ASP/GLU amine partner ──────────────────────────────────────────────────────

// placeCarboxylateAmines places a methylammonium donor bridging each ASP/GLU
// carboxylate, the amine counterpart to placeCarboxylateGuanidines.
func placeCarboxylateAmines(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type carboxyl struct {
		prevC, carboxylC, O1, O2 Vec3
		found                    [4]bool
	}
	groups := map[residueKey]*carboxyl{}
	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		if res != "ASP" && res != "GLU" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &carboxyl{}
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
		baseZ := g.carboxylC.Sub(g.prevC).unit()
		for _, sign := range []float64{1, -1} {
			zHat := baseZ.Scale(sign)
			nPos := g.carboxylC.Add(zHat.Scale(2.90)) // ammonium N bridges both O's
			cPos := nPos.Add(zHat.Scale(1.47))        // capping carbon
			if addAmmoniumProbe(placedPos, placedLabel, bonds, clash, key, nPos, cPos) {
				placed++
				break
			}
		}
	}
	return placed
}

// ── MET S–arene ────────────────────────────────────────────────────────────────

// placeMetArenes places a benzene ring face-on to each MET thioether sulfur
// (SD), capturing the S–arene interaction the hydrophobic chain alone misses.
func placeMetArenes(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type met struct {
		CG, SD, CE Vec3
		found      [3]bool
	}
	groups := map[residueKey]*met{}
	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		if res != "MET" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &met{}
		}
		g := groups[key]
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "CG":
			g.CG, g.found[0] = a.Pos, true
		case "SD":
			g.SD, g.found[1] = a.Pos, true
		case "CE":
			g.CE, g.found[2] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)
	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] || !g.found[2] {
			continue
		}
		// Lone-pair direction: external bisector of CG–SD–CE.
		bis := g.SD.Sub(g.CG).unit().Add(g.SD.Sub(g.CE).unit())
		if bis.Norm() < 1e-6 {
			continue
		}
		bis = bis.unit()
		centre := g.SD.Add(bis.Scale(4.3))
		if addBenzeneProbe(placedPos, placedLabel, bonds, clash, key, centre, bis) {
			placed++
		}
	}
	return placed
}

// ── Cation–π ───────────────────────────────────────────────────────────────────

// placeAromaticCations places a methylammonium probe centred over both faces of
// each PHE/TYR/TRP aromatic ring (cation–π from the ligand side).
func placeAromaticCations(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	wanted := map[string]map[string]bool{
		"PHE": {"CG": true, "CD1": true, "CD2": true, "CE1": true, "CE2": true, "CZ": true},
		"TYR": {"CG": true, "CD1": true, "CD2": true, "CE1": true, "CE2": true, "CZ": true},
		"TRP": {"CD2": true, "CE2": true, "CE3": true, "CZ2": true, "CZ3": true, "CH2": true},
	}
	groups := map[residueKey]map[string]Vec3{}
	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		w, ok := wanted[res]
		if !ok {
			continue
		}
		name := strings.TrimSpace(strings.ToUpper(a.Name))
		if !w[name] {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = map[string]Vec3{}
		}
		groups[key][name] = a.Pos
	}

	clash := newResClashSet(atoms)
	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		m := groups[key]
		if len(m) < 6 {
			continue
		}
		pts := make([]Vec3, 0, 6)
		centroid := Vec3{}
		for _, p := range m {
			pts = append(pts, p)
			centroid = centroid.Add(p)
		}
		centroid = centroid.Scale(1.0 / float64(len(pts)))
		normal := ringPlaneNormal(pts)
		for _, sign := range []float64{1, -1} {
			nPos := centroid.Add(normal.Scale(sign * 3.4))
			cPos := centroid.Add(normal.Scale(sign * (3.4 + 1.47)))
			if addAmmoniumProbe(placedPos, placedLabel, bonds, clash, key, nPos, cPos) {
				placed++
				break
			}
		}
	}
	return placed
}

// placeCationArenes places a benzene ring against each ARG guanidinium plane and
// each LYS ammonium (cation–π from the protein side).
func placeCationArenes(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type arg struct {
		NE, CZ, NH1, NH2 Vec3
		found            [4]bool
	}
	type lys struct {
		CE, NZ Vec3
		found  [2]bool
	}
	args := map[residueKey]*arg{}
	lyss := map[residueKey]*lys{}
	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		key := residueKey{a.ChainID, a.ResSeq}
		name := strings.TrimSpace(strings.ToUpper(a.Name))
		switch res {
		case "ARG":
			if args[key] == nil {
				args[key] = &arg{}
			}
			switch name {
			case "NE":
				args[key].NE, args[key].found[0] = a.Pos, true
			case "CZ":
				args[key].CZ, args[key].found[1] = a.Pos, true
			case "NH1":
				args[key].NH1, args[key].found[2] = a.Pos, true
			case "NH2":
				args[key].NH2, args[key].found[3] = a.Pos, true
			}
		case "LYS":
			if lyss[key] == nil {
				lyss[key] = &lys{}
			}
			switch name {
			case "CE":
				lyss[key].CE, lyss[key].found[0] = a.Pos, true
			case "NZ":
				lyss[key].NZ, lyss[key].found[1] = a.Pos, true
			}
		}
	}

	clash := newResClashSet(atoms)
	placed := 0
	for _, key := range sortedResidueKeys(args) {
		g := args[key]
		if !g.found[0] || !g.found[1] || !g.found[2] || !g.found[3] {
			continue
		}
		pts := []Vec3{g.NE, g.CZ, g.NH1, g.NH2}
		centroid := Vec3{}
		for _, p := range pts {
			centroid = centroid.Add(p)
		}
		centroid = centroid.Scale(0.25)
		normal := ringPlaneNormal(pts)
		for _, sign := range []float64{1, -1} {
			centre := centroid.Add(normal.Scale(sign * 3.7))
			if addBenzeneProbe(placedPos, placedLabel, bonds, clash, key, centre, normal) {
				placed++
				break
			}
		}
	}
	for _, key := range sortedResidueKeys(lyss) {
		g := lyss[key]
		if !g.found[0] || !g.found[1] {
			continue
		}
		zHat := g.NZ.Sub(g.CE).unit()
		centre := g.NZ.Add(zHat.Scale(3.8))
		if addBenzeneProbe(placedPos, placedLabel, bonds, clash, key, centre, zHat) {
			placed++
		}
	}
	return placed
}

// ── Hydroxyl / methoxy probes ──────────────────────────────────────────────────

// clearanceDirs returns up to n of the 26 cube directions from pos, ordered by
// decreasing clearance from the protein, keeping only well-separated ones.
func clearanceDirs(pos Vec3, heavy []heavyAtom, n int) []Vec3 {
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
	clearance := func(d Vec3) float64 {
		cand := pos.Add(d.Scale(1.43))
		minC := 999.0
		for _, h := range heavy {
			if gap := cand.Sub(h.pos).Norm() - h.vdwR - vdw("C"); gap < minC {
				minC = gap
			}
		}
		return minC
	}
	sort.Slice(dirs, func(i, j int) bool { return clearance(dirs[i]) > clearance(dirs[j]) })

	out := []Vec3{}
	for _, d := range dirs {
		if clearance(d) < -0.5 {
			break // remainder are all buried
		}
		ok := true
		for _, e := range out {
			if d.Dot(e) > 0.5 { // keep directions > 60° apart
				ok = false
				break
			}
		}
		if ok {
			out = append(out, d)
		}
		if len(out) == n {
			break
		}
	}
	return out
}

// placeHydroxylVotes places hydroxyl (O–C, donates) probes at H-bond hotspots
// near protein acceptors and methoxy (C–O–C, accepts) probes near protein
// donors, using the same vote-cell machinery as buildHbondProbes but emitting
// full functional groups instead of bare atoms.
func placeHydroxylVotes(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	var heavy []heavyAtom
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, heavyAtom{a.Pos, vdw(a.Element)})
	}
	var donors, acceptors []Atom
	for _, a := range atoms {
		switch hbondRole(a) {
		case "donor":
			donors = append(donors, a)
		case "acceptor":
			acceptors = append(acceptors, a)
		case "dual":
			donors = append(donors, a)
			acceptors = append(acceptors, a)
		}
	}

	var taken []Vec3
	placed := 0
	place := func(group []Atom, nCarbons int) {
		votes := castVotes(group, hbondSphereR, hbondPairDist, nSamples, heavy, "O")
		type cand struct {
			pos Vec3
			pc  int
		}
		var cands []cand
		for _, v := range votes {
			if v.pairCount >= hbondMinPairs {
				cands = append(cands, cand{v.pos, v.pairCount})
			}
		}
		sort.Slice(cands, func(i, j int) bool { return cands[i].pc > cands[j].pc })
		for _, c := range cands {
			if !noHardClash(c.pos, vdw("O"), hardTol, heavy) {
				continue
			}
			tooClose := false
			for _, t := range taken {
				if c.pos.Sub(t).Norm() < hbondProbeSpacing {
					tooClose = true
					break
				}
			}
			if tooClose {
				continue
			}
			dirs := clearanceDirs(c.pos, heavy, 1)
			if len(dirs) < 1 {
				continue
			}
			carbonDirs := []Vec3{dirs[0]}
			if nCarbons == 2 {
				// Second carbon must sit at the ether angle (~111°) from the first,
				// not at an arbitrary clearance direction — a C–O–C ether is bent.
				// Sweep the 111° cone around the first bond and keep the most open one.
				d1 := dirs[0]
				perp := perpendicular(d1)
				q := cross3(d1, perp)
				cosT, sinT := math.Cos(111*math.Pi/180), math.Sin(111*math.Pi/180)
				best := Vec3{}
				bestClear := -1e9
				for k := 0; k < 24; k++ {
					a := float64(k) * 2 * math.Pi / 24
					d2 := d1.Scale(cosT).Add(perp.Scale(math.Cos(a) * sinT)).Add(q.Scale(math.Sin(a) * sinT)).unit()
					cand := c.pos.Add(d2.Scale(1.43))
					minGap := 1e9
					for _, h := range heavy {
						if g := cand.Sub(h.pos).Norm() - h.vdwR - vdw("C"); g < minGap {
							minGap = g
						}
					}
					if minGap > bestClear {
						bestClear = minGap
						best = d2
					}
				}
				if bestClear < -0.5 { // even the best second carbon is buried
					continue
				}
				carbonDirs = append(carbonDirs, best)
			}
			oIdx := len(*placedPos) + 1
			*placedPos = append(*placedPos, c.pos)
			*placedLabel = append(*placedLabel, "O")
			for _, d := range carbonDirs {
				*placedPos = append(*placedPos, c.pos.Add(d.Scale(1.43)))
				*placedLabel = append(*placedLabel, "C")
				*bonds = append(*bonds, Bond{oIdx, len(*placedPos), 1})
			}
			taken = append(taken, c.pos)
			placed++
		}
	}
	place(acceptors, 1) // hydroxyl: O donates to protein acceptor, one carbon
	place(donors, 2)    // methoxy: O accepts from protein donor, two carbons
	return placed
}

// placeBackboneHydroxyls places a hydroxyl donor at every backbone carbonyl
// oxygen — backbone acceptors are excluded from the vote system (hbondRole skips
// backbone), so this fills a major real-world H-bond gap.
func placeBackboneHydroxyls(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	type bb struct {
		C, O  Vec3
		found [2]bool
	}
	groups := map[residueKey]*bb{}
	for _, a := range atoms {
		if a.IsHet || isWater(a) || strings.ToUpper(a.Element) == "H" {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		switch strings.TrimSpace(strings.ToUpper(a.Name)) {
		case "C":
			if groups[key] == nil {
				groups[key] = &bb{}
			}
			groups[key].C, groups[key].found[0] = a.Pos, true
		case "O":
			if groups[key] == nil {
				groups[key] = &bb{}
			}
			groups[key].O, groups[key].found[1] = a.Pos, true
		}
	}

	clash := newResClashSet(atoms)
	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		if !g.found[0] || !g.found[1] {
			continue
		}
		// Outward C=O direction; the probe O–H donates back to the carbonyl O.
		dir := g.O.Sub(g.C).unit()
		probeO := g.O.Add(dir.Scale(2.80))
		probeC := g.O.Add(dir.Scale(2.80 + 1.43))
		if clash.clears(probeO, "O", key) && clash.clears(probeC, "C", key) {
			base := len(*placedPos) + 1
			*placedPos = append(*placedPos, probeO, probeC)
			*placedLabel = append(*placedLabel, "O", "C")
			*bonds = append(*bonds, Bond{base, base + 1, 1})
			placed++
		}
	}
	return placed
}
