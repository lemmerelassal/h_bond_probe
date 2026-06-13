package main

import (
	"math"
	"sort"
	"strings"
)

// This file adds probes for the functional groups that dominate recent (last ~20
// years) FDA-approved small molecules but were missing from the earlier sets:
//
//   • saturated N-heterocycles  : piperidine, pyrrolidine, piperazine, morpholine
//   • aromatic N-heterocycles   : pyridine, pyrimidine, pyrazole, 1,3,4-oxadiazole
//   • halogen bond donors       : C–Cl, C–Br, C–F
//   • trifluoromethyl           : hydrophobic CF3 caps on aliphatic side chains
//   • sulfonamide               : methanesulfonamide (donor NH2 + acceptor SO2)
//   • nitrile                   : acetonitrile (linear C≡N acceptor)
//   • ester / carbamate         : methyl acetate / methyl carbamate (acceptors)
//
// All of the hotspot-driven placers reuse the H-bond vote machinery (castVotes)
// the same way placeHydroxylVotes does: acceptor-emitting groups are placed at
// "donor hotspots" (positions clustered near protein donors, so the probe's
// acceptor faces a donor) and donor-emitting groups at "acceptor hotspots".

// ── Shared hotspot helpers ─────────────────────────────────────────────────────

// hbondSets splits the protein into a clash list plus donor and acceptor atom
// lists, mirroring the classification placeHydroxylVotes uses.
func hbondSets(atoms []Atom) (heavy []heavyAtom, donors, acceptors []Atom) {
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, heavyAtom{a.Pos, vdw(a.Element)})
	}
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
	return
}

// hbondHotspots returns clear vote-cell positions for the given protein group,
// ordered by vote strength (with a positional tiebreak so the order — and hence
// the placement — is deterministic run to run).
func hbondHotspots(group []Atom, heavy []heavyAtom) []Vec3 {
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
	sort.Slice(cands, func(i, j int) bool {
		if cands[i].pc != cands[j].pc {
			return cands[i].pc > cands[j].pc
		}
		if cands[i].pos.X != cands[j].pos.X {
			return cands[i].pos.X < cands[j].pos.X
		}
		if cands[i].pos.Y != cands[j].pos.Y {
			return cands[i].pos.Y < cands[j].pos.Y
		}
		return cands[i].pos.Z < cands[j].pos.Z
	})
	out := make([]Vec3, 0, len(cands))
	for _, c := range cands {
		if noHardClash(c.pos, vdw("O"), hardTol, heavy) {
			out = append(out, c.pos)
		}
	}
	return out
}

// normalsAbout returns a few candidate unit vectors perpendicular to axis,
// obtained by rotating an arbitrary perpendicular about axis in equal steps.
func normalsAbout(axis Vec3, count int) []Vec3 {
	n0 := perpendicular(axis)
	t := cross3(axis, n0).unit()
	out := make([]Vec3, 0, count)
	for k := 0; k < count; k++ {
		ang := float64(k) * math.Pi / float64(count)
		out = append(out, n0.Scale(math.Cos(ang)).Add(t.Scale(math.Sin(ang))).unit())
	}
	return out
}

// ── Heterocyclic rings ─────────────────────────────────────────────────────────

// ringSpec describes one heterocyclic probe ring as a regular polygon.
//
//   elements   — ring atom element symbols, in connectivity order
//   bondOrders — order of ring edge k → (k, (k+1) mod n); chosen as an explicit
//                Kekulé pattern so no ring atom is ever over-valent (which would
//                make sanitizeValence strip a ring bond and break the ring — the
//                exact failure already seen with tetrazoles)
//   nhDonor    — true if the directional atom is an N–H donor (ring is placed at
//                an acceptor hotspot); false if it is a lone-pair acceptor
//                (placed at a donor hotspot)
//   dirIdx     — ring atom that sits at the hotspot and makes the H-bond
//   capIdx     — ring atom that carries the exocyclic methyl linker handle
type ringSpec struct {
	name       string
	elements   []string
	bondOrders []int
	aromatic   bool
	nhDonor    bool
	dirIdx     int
	capIdx     int
}

var heteroRingSpecs = []ringSpec{
	// Aromatic, lone-pair acceptor → donor hotspots.
	{"pyridine", []string{"N", "C", "C", "C", "C", "C"},
		[]int{2, 1, 2, 1, 2, 1}, true, false, 0, 3},
	{"pyrimidine", []string{"N", "C", "N", "C", "C", "C"},
		[]int{2, 1, 2, 1, 2, 1}, true, false, 0, 4},
	{"oxadiazole", []string{"O", "C", "N", "N", "C"},
		[]int{1, 2, 1, 2, 1}, true, false, 2, 4},
	// Aromatic, N–H donor → acceptor hotspots.
	{"pyrazole", []string{"N", "N", "C", "C", "C"},
		[]int{1, 2, 1, 2, 1}, true, true, 0, 3},
	// Saturated, N–H donor → acceptor hotspots.
	{"piperidine", []string{"N", "C", "C", "C", "C", "C"},
		[]int{1, 1, 1, 1, 1, 1}, false, true, 0, 3},
	{"pyrrolidine", []string{"N", "C", "C", "C", "C"},
		[]int{1, 1, 1, 1, 1}, false, true, 0, 3},
	{"piperazine", []string{"N", "C", "C", "N", "C", "C"},
		[]int{1, 1, 1, 1, 1, 1}, false, true, 0, 3},
	// Saturated, ether-O acceptor → donor hotspots.
	{"morpholine", []string{"O", "C", "C", "N", "C", "C"},
		[]int{1, 1, 1, 1, 1, 1}, false, false, 0, 3},
}

// emitHeteroRing places one ring whose directional atom sits at hotspot and whose
// body extends outward along outDir, lying in the plane that contains outDir and
// has the given normal. Returns false on any clash.
func emitHeteroRing(spec ringSpec, hotspot, outDir, normal Vec3,
	placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond, heavy []heavyAtom) bool {
	n := len(spec.elements)
	side := 1.52
	if spec.aromatic {
		side = 1.39
	}
	R := side / (2 * math.Sin(math.Pi/float64(n)))
	center := hotspot.Add(outDir.Scale(R)) // directional vertex sits at hotspot
	a := outDir.Scale(-1)                  // center → directional vertex
	b := cross3(normal, a)
	if b.Norm() < 1e-6 {
		return false
	}
	b = b.unit()
	step := 2 * math.Pi / float64(n)

	pos := make([]Vec3, n)
	for k := 0; k < n; k++ {
		ang := float64(k-spec.dirIdx) * step
		pos[k] = center.Add(a.Scale(R * math.Cos(ang))).Add(b.Scale(R * math.Sin(ang)))
	}
	capDir := pos[spec.capIdx].Sub(center).unit()
	cap := pos[spec.capIdx].Add(capDir.Scale(1.50))

	for k, p := range pos {
		if !noHardClash(p, vdw(spec.elements[k]), hardTol, heavy) {
			return false
		}
	}
	if !noHardClash(cap, vdw("C"), hardTol, heavy) {
		return false
	}

	base := len(*placedPos) + 1
	for k, p := range pos {
		*placedPos = append(*placedPos, p)
		*placedLabel = append(*placedLabel, spec.elements[k])
	}
	for k := 0; k < n; k++ {
		*bonds = append(*bonds, Bond{base + k, base + (k+1)%n, spec.bondOrders[k]})
	}
	*placedPos = append(*placedPos, cap)
	*placedLabel = append(*placedLabel, "C")
	*bonds = append(*bonds, Bond{base + spec.capIdx, len(*placedPos), 1})
	return true
}

// placeHeteroRings places every heterocyclic probe ring at its matching H-bond
// hotspots (donor hotspots for acceptor rings, acceptor hotspots for N–H donor
// rings). Each ring type dedupes against its own placements but is allowed to
// share a hotspot with other ring types, so several scaffold alternatives are
// explored at the same site (as carboxylate and tetrazole already are).
func placeHeteroRings(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	heavy, donors, acceptors := hbondSets(atoms)
	donorHots := hbondHotspots(donors, heavy)
	accHots := hbondHotspots(acceptors, heavy)

	placed := 0
	for _, spec := range heteroRingSpecs {
		hots := donorHots
		if spec.nhDonor {
			hots = accHots
		}
		var taken []Vec3
		for _, hp := range hots {
			tooClose := false
			for _, t := range taken {
				if hp.Sub(t).Norm() < 2.0 {
					tooClose = true
					break
				}
			}
			if tooClose {
				continue
			}
			done := false
			for _, outDir := range clearanceDirs(hp, heavy, 3) {
				for _, normal := range normalsAbout(outDir, 3) {
					if emitHeteroRing(spec, hp, outDir, normal, placedPos, placedLabel, bonds, heavy) {
						done = true
						break
					}
				}
				if done {
					break
				}
			}
			if done {
				taken = append(taken, hp)
				placed++
			}
		}
	}
	return placed
}

// ── Halogen bond donors ─────────────────────────────────────────────────────────

// placeHalogenBonds places C–X probes (X = Cl, Br, F) whose halogen sits at an
// acceptor hotspot with its σ-hole pointing back toward the protein acceptor and
// the carbon extending into open space. Cl and Br are genuine halogen-bond
// donors; F is included as a weak σ-hole / dipole probe.
func placeHalogenBonds(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	heavy, _, acceptors := hbondSets(atoms)
	hots := hbondHotspots(acceptors, heavy)

	type hal struct {
		elem string
		cx   float64 // C–X bond length
	}
	halogens := []hal{{"CL", 1.77}, {"BR", 1.94}, {"F", 1.36}}

	placed := 0
	for _, h := range halogens {
		var taken []Vec3
		for _, hp := range hots {
			tooClose := false
			for _, t := range taken {
				if hp.Sub(t).Norm() < 1.5 {
					tooClose = true
					break
				}
			}
			if tooClose {
				continue
			}
			for _, dir := range clearanceDirs(hp, heavy, 3) {
				c := hp.Add(dir.Scale(h.cx))
				if noHardClash(hp, vdw(h.elem), hardTol, heavy) && noHardClash(c, vdw("C"), hardTol, heavy) {
					base := len(*placedPos) + 1
					*placedPos = append(*placedPos, hp, c)
					*placedLabel = append(*placedLabel, h.elem, "C")
					*bonds = append(*bonds, Bond{base, base + 1, 1})
					taken = append(taken, hp)
					placed++
					break
				}
			}
		}
	}
	return placed
}

// ── Trifluoromethyl ─────────────────────────────────────────────────────────────

// placeTrifluoromethyls extends each aliphatic side-chain methyl/methylene tip
// with a CF3 group pointing outward — the metabolically robust hydrophobic motif
// ubiquitous in modern drugs. The CF3 carbon is the linker handle.
func placeTrifluoromethyls(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	// Each entry: side-chain tip atom and the atom it points away from.
	tipParents := map[string][][2]string{
		"LEU": {{"CD1", "CG"}, {"CD2", "CG"}},
		"VAL": {{"CG1", "CB"}, {"CG2", "CB"}},
		"ILE": {{"CD1", "CG1"}, {"CG2", "CB"}},
		"ALA": {{"CB", "CA"}},
		"THR": {{"CG2", "CB"}},
		"MET": {{"CE", "SD"}},
	}
	type resAtoms struct {
		res string
		pos map[string]Vec3
	}
	groups := map[residueKey]*resAtoms{}
	for _, a := range atoms {
		res := strings.TrimSpace(strings.ToUpper(a.ResName))
		if _, ok := tipParents[res]; !ok {
			continue
		}
		key := residueKey{a.ChainID, a.ResSeq}
		if groups[key] == nil {
			groups[key] = &resAtoms{res: res, pos: map[string]Vec3{}}
		}
		groups[key].pos[strings.TrimSpace(strings.ToUpper(a.Name))] = a.Pos
	}

	clash := newResClashSet(atoms)
	placed := 0
	for _, key := range sortedResidueKeys(groups) {
		g := groups[key]
		for _, tp := range tipParents[g.res] {
			tip, okT := g.pos[tp[0]]
			parent, okP := g.pos[tp[1]]
			if !okT || !okP {
				continue
			}
			dir := tip.Sub(parent).unit()
			cC := tip.Add(dir.Scale(1.52))
			// Three F on a cone ~70.5° off dir (≈109.5° from the C–tip bond).
			p := perpendicular(dir)
			q := cross3(dir, p)
			cosA, sinA := math.Cos(70.5*math.Pi/180), math.Sin(70.5*math.Pi/180)
			fs := make([]Vec3, 3)
			ok := clash.clears(cC, "C", key)
			for i := 0; i < 3 && ok; i++ {
				phi := float64(i) * 2 * math.Pi / 3
				fd := dir.Scale(cosA).Add(p.Scale(math.Cos(phi) * sinA)).Add(q.Scale(math.Sin(phi) * sinA))
				fs[i] = cC.Add(fd.Scale(1.33))
				if !clash.clears(fs[i], "F", key) {
					ok = false
				}
			}
			if !ok {
				continue
			}
			base := len(*placedPos) + 1
			*placedPos = append(*placedPos, cC, fs[0], fs[1], fs[2])
			*placedLabel = append(*placedLabel, "C", "F", "F", "F")
			*bonds = append(*bonds,
				Bond{base, base + 1, 1},
				Bond{base, base + 2, 1},
				Bond{base, base + 3, 1})
			placed++
		}
	}
	return placed
}

// ── Sulfonamide ─────────────────────────────────────────────────────────────────

// placeSulfonamides places a methanesulfonamide (CH3–SO2–NH2) at acceptor
// hotspots, with the amide N donating to the protein acceptor and the sulfonyl
// oxygens accepting from any nearby donors. The sulfur is S(VI) (valence 6).
func placeSulfonamides(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	heavy, _, acceptors := hbondSets(atoms)
	hots := hbondHotspots(acceptors, heavy)

	placed := 0
	var taken []Vec3
	for _, hp := range hots {
		tooClose := false
		for _, t := range taken {
			if hp.Sub(t).Norm() < 2.0 {
				tooClose = true
				break
			}
		}
		if tooClose {
			continue
		}
		done := false
		for _, dir := range clearanceDirs(hp, heavy, 3) {
			s := hp.Add(dir.Scale(1.63)) // S–N ≈ 1.63 Å
			p := perpendicular(dir)
			q := cross3(dir, p)
			// O, O, C on a tetrahedral cone opposite the S→N bond.
			cosA, sinA := math.Cos(70.5*math.Pi/180), math.Sin(70.5*math.Pi/180)
			leg := func(phi float64) Vec3 {
				return dir.Scale(cosA).Add(p.Scale(math.Cos(phi) * sinA)).Add(q.Scale(math.Sin(phi) * sinA))
			}
			o1 := s.Add(leg(0).Scale(1.44))
			o2 := s.Add(leg(2 * math.Pi / 3).Scale(1.44))
			cMe := s.Add(leg(4 * math.Pi / 3).Scale(1.50))
			if noHardClash(hp, vdw("N"), hardTol, heavy) && noHardClash(s, vdw("S"), hardTol, heavy) &&
				noHardClash(o1, vdw("O"), hardTol, heavy) && noHardClash(o2, vdw("O"), hardTol, heavy) &&
				noHardClash(cMe, vdw("C"), hardTol, heavy) {
				base := len(*placedPos) + 1
				*placedPos = append(*placedPos, hp, s, o1, o2, cMe)
				*placedLabel = append(*placedLabel, "N", "S", "O", "O", "C")
				*bonds = append(*bonds,
					Bond{base, base + 1, 1},     // N–S
					Bond{base + 1, base + 2, 2}, // S=O
					Bond{base + 1, base + 3, 2}, // S=O
					Bond{base + 1, base + 4, 1}, // S–CH3
				)
				taken = append(taken, hp)
				placed++
				done = true
				break
			}
		}
		_ = done
	}
	return placed
}

// ── Nitrile ─────────────────────────────────────────────────────────────────────

// placeNitriles places an acetonitrile (CH3–C≡N) at donor hotspots: the nitrile
// nitrogen accepts from the protein donor, the linear body pointing outward.
func placeNitriles(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	heavy, donors, _ := hbondSets(atoms)
	hots := hbondHotspots(donors, heavy)

	placed := 0
	var taken []Vec3
	for _, hp := range hots {
		tooClose := false
		for _, t := range taken {
			if hp.Sub(t).Norm() < 1.5 {
				tooClose = true
				break
			}
		}
		if tooClose {
			continue
		}
		for _, dir := range clearanceDirs(hp, heavy, 3) {
			c := hp.Add(dir.Scale(1.16)) // C≡N
			me := c.Add(dir.Scale(1.46)) // C–CH3
			if noHardClash(hp, vdw("N"), hardTol, heavy) && noHardClash(c, vdw("C"), hardTol, heavy) &&
				noHardClash(me, vdw("C"), hardTol, heavy) {
				base := len(*placedPos) + 1
				*placedPos = append(*placedPos, hp, c, me)
				*placedLabel = append(*placedLabel, "N", "C", "C")
				*bonds = append(*bonds,
					Bond{base, base + 1, 3},     // N≡C
					Bond{base + 1, base + 2, 1}, // C–CH3
				)
				taken = append(taken, hp)
				placed++
				break
			}
		}
	}
	return placed
}

// ── Ester / carbamate ───────────────────────────────────────────────────────────

// placeEsters places a methyl ester (CH3–C(=O)–O–CH3) and a methyl carbamate
// (H2N–C(=O)–O–CH3) at donor hotspots, with the carbonyl oxygen accepting from
// the protein donor. The two differ only in the acyl substituent (C vs N).
func placeEsters(atoms []Atom, placedPos *[]Vec3, placedLabel *[]string, bonds *[]Bond) int {
	heavy, donors, _ := hbondSets(atoms)
	hots := hbondHotspots(donors, heavy)

	type espec struct {
		name     string
		acylElem string
		acylLen  float64
	}
	specs := []espec{{"ester", "C", 1.50}, {"carbamate", "N", 1.36}}

	placed := 0
	for _, sp := range specs {
		var taken []Vec3
		for _, hp := range hots {
			tooClose := false
			for _, t := range taken {
				if hp.Sub(t).Norm() < 2.0 {
					tooClose = true
					break
				}
			}
			if tooClose {
				continue
			}
			done := false
			for _, dir := range clearanceDirs(hp, heavy, 3) {
				for _, normal := range normalsAbout(dir, 3) {
					c := hp.Add(dir.Scale(1.22)) // carbonyl C
					dCO := hp.Sub(c).unit()      // C→O (= −dir)
					p := cross3(normal, dCO)
					if p.Norm() < 1e-6 {
						continue
					}
					p = p.unit()
					cos120, sin120 := math.Cos(2*math.Pi/3), math.Sin(2*math.Pi/3)
					acylDir := dCO.Scale(cos120).Add(p.Scale(sin120))
					estDir := dCO.Scale(cos120).Sub(p.Scale(sin120))
					acyl := c.Add(acylDir.Scale(sp.acylLen))
					estO := c.Add(estDir.Scale(1.34))
					// O–CH3 must be BENT (~120°) off the C–O bond, not collinear: a
					// divalent ether/ester oxygen is sp3-like, never linear. Rotate the
					// O→C direction 120° about the ester-plane normal.
					oToC := c.Sub(estO).unit()
					estMeDir := oToC.Scale(math.Cos(2 * math.Pi / 3)).Add(cross3(normal, oToC).Scale(math.Sin(2 * math.Pi / 3)))
					estMe := estO.Add(estMeDir.unit().Scale(1.43))
					if noHardClash(hp, vdw("O"), hardTol, heavy) && noHardClash(c, vdw("C"), hardTol, heavy) &&
						noHardClash(acyl, vdw(sp.acylElem), hardTol, heavy) &&
						noHardClash(estO, vdw("O"), hardTol, heavy) && noHardClash(estMe, vdw("C"), hardTol, heavy) {
						base := len(*placedPos) + 1
						*placedPos = append(*placedPos, hp, c, acyl, estO, estMe)
						*placedLabel = append(*placedLabel, "O", "C", sp.acylElem, "O", "C")
						*bonds = append(*bonds,
							Bond{base, base + 1, 2},     // O=C
							Bond{base + 1, base + 2, 1}, // C–acyl
							Bond{base + 1, base + 3, 1}, // C–O(ester)
							Bond{base + 3, base + 4, 1}, // O–CH3
						)
						taken = append(taken, hp)
						placed++
						done = true
						break
					}
				}
				if done {
					break
				}
			}
		}
	}
	return placed
}
