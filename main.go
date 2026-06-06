package main

import (
	"fmt"
	"math"
	"os"
	"sort"
	"strings"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Fprintln(os.Stderr, "Usage: pharmacophore <PDB-ID|file.pdb|file.cif>")
		os.Exit(1)
	}

	atoms, pdbID, err := LoadAtoms(os.Args[1])
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
	if len(atoms) == 0 {
		fmt.Fprintln(os.Stderr, "No atoms parsed.")
		os.Exit(1)
	}

	// ── Partition atoms ───────────────────────────────────────────────────────
	// heavy: all non-H, non-water atoms (used for clash checking, includes HET)
	// hydro: protein hydrophobic C atoms only (source for vote system)
	// proteinAtoms: protein-only atom slice fed to all Placers
	var hydro []Atom
	var heavy []heavyAtom
	var proteinAtoms []Atom
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, heavyAtom{a.Pos, vdw(a.Element)})
		if a.IsHet {
			continue
		}
		proteinAtoms = append(proteinAtoms, a)
		if isHydrophobic(a) {
			hydro = append(hydro, a)
		}
	}
	fmt.Printf("Hydrophobic atoms: %d  Heavy atoms total: %d\n", len(hydro), len(heavy))

	// ── Initialise the probe accumulator ─────────────────────────────────────
	ps := &ProbeSet{}

	// ── Alkane chain through hydrophobic intersection points ──────────────────
	buildAlkaneChain(hydro, heavy, ps)

	// ── Benzene / imidazole ring probes for aromatic residues ─────────────────
	buildAromaticProbes(proteinAtoms, heavy, ps)

	// ── H-bond acceptor / donor vote probes ──────────────────────────────────
	buildHbondProbes(proteinAtoms, heavy, ps)

	// ── Pharmacophoric probes (Open/Closed: iterate AllPlacers) ───────────────
	for _, p := range AllPlacers {
		n := p.Place(proteinAtoms, heavy, ps)
		fmt.Printf("%s probes: %d\n", p.Name(), n)
	}

	// ── Link probe groups into multi-fragment molecules ───────────────────────
	linkableN := ps.Len()
	nLinkers := linkProbeGroups(&ps.Pos, &ps.Labels, &ps.Bonds, heavy, 20.0, linkableN)
	fmt.Printf("Linker chains: %d\n", nLinkers)

	// ── Backbone ──────────────────────────────────────────────────────────────
	nBB := addBackbone(proteinAtoms, &ps.Pos, &ps.Labels, &ps.Bonds)
	fmt.Printf("Backbone residues: %d\n", nBB)

	// ── Write output as SDF (not .mol) ────────────────────────────────────────
	outPath := pdbID + "_probes.sdf"
	if err := WriteSDF(outPath, ps); err != nil {
		fmt.Fprintf(os.Stderr, "Error writing SDF: %v\n", err)
		os.Exit(1)
	}
	fmt.Printf("Wrote %d probe atoms -> %s\n", ps.Len(), outPath)
}

// ── Hydrophobic alkane chain ──────────────────────────────────────────────────

func buildAlkaneChain(hydro []Atom, heavy []heavyAtom, ps *ProbeSet) {
	votes := map[[3]int]*voteCell{}
	nPairs := 0
	for i := 0; i < len(hydro); i++ {
		for j := i + 1; j < len(hydro); j++ {
			if hydro[i].Pos.Sub(hydro[j].Pos).Norm() > maxPairDist {
				continue
			}
			pts := sphereCircle(hydro[i].Pos, hydro[j].Pos, sphereR, nSamples)
			if len(pts) == 0 {
				continue
			}
			nPairs++
			pk := [2]int{i, j}
			for _, p := range pts {
				if !noHardClash(p, vdw("C"), hardTol, heavy) {
					continue
				}
				key := gridKey(p)
				if e, ok := votes[key]; ok {
					if !e.pairs[pk] {
						e.pairs[pk] = true
						e.pairCount++
					}
					e.nPts++
					n := float64(e.nPts)
					e.pos = Vec3{
						(e.pos.X*(n-1) + p.X) / n,
						(e.pos.Y*(n-1) + p.Y) / n,
						(e.pos.Z*(n-1) + p.Z) / n,
					}
				} else {
					votes[key] = &voteCell{
						pairCount: 1, pos: p, nPts: 1,
						pairs: map[[2]int]bool{pk: true},
					}
				}
			}
		}
	}
	fmt.Printf("Pairs within %.1f Å: %d  Vote cells: %d\n", maxPairDist, nPairs, len(votes))

	// burialCount returns the number of protein heavy atoms within shellR Å
	// of pos that are also beyond minDist Å (excluding the immediate vdW shell).
	// High burial = deep in a pocket; low burial = exposed on surface.
	const (
		burialShellR  = 8.0 // Å — count atoms within this radius
		burialMinDist = 2.5 // Å — ignore atoms too close (vdW contacts)
		minBurial     = 6   // minimum surrounding atoms to accept as pocket
	)
	burialCount := func(pos Vec3) int {
		n := 0
		for _, h := range heavy {
			d := pos.Sub(h.pos).Norm()
			if d >= burialMinDist && d <= burialShellR {
				n++
			}
		}
		return n
	}

	type candidate struct {
		pos       Vec3
		pairCount int
	}
	var cands []candidate
	for _, v := range votes {
		if v.pairCount >= minPairs {
			cands = append(cands, candidate{v.pos, v.pairCount})
		}
	}
	sort.Slice(cands, func(i, j int) bool { return cands[i].pairCount > cands[j].pairCount })

	var waypoints []Vec3
	nClash, nClose, nSurface := 0, 0, 0
	for _, c := range cands {
		if !noHardClash(c.pos, vdw("C"), hardTol, heavy) {
			nClash++
			continue
		}
		if burialCount(c.pos) < minBurial {
			nSurface++
			continue
		}
		tooClose := false
		for _, wp := range waypoints {
			if c.pos.Sub(wp).Norm() < probeSpacing {
				tooClose = true
				break
			}
		}
		if tooClose {
			nClose++
			continue
		}
		waypoints = append(waypoints, c.pos)
	}
	fmt.Printf("Placed: %d  Rejected (clash): %d  Rejected (too close): %d  Rejected (surface): %d\n",
		len(waypoints), nClash, nClose, nSurface)

	if len(waypoints) == 0 {
		fmt.Println("No hydrophobic vote cells — no alkane chain built.")
		return
	}
	chainAtoms := alkaneChain(waypoints, heavy)
	if len(chainAtoms) == 0 {
		return
	}
	total := 0
	for _, seg := range chainAtoms {
		base := ps.Len() + 1
		for _, p := range seg {
			ps.Add(p, "C")
		}
		for i := 0; i+1 < len(seg); i++ {
			ps.Bond(base+i, base+i+1, 1)
		}
		total += len(seg)
	}
	fmt.Printf("Alkane chain atoms: %d in %d segment(s)\n", total, len(chainAtoms))
}

// alkaneChain orders waypoints by nearest-neighbour and builds all-trans sp3
// chain segments threading through them. Returns a slice of segments; atoms
// within each segment are bonded consecutively, but NO bond is written between
// segments (skipped waypoints would otherwise create impossible long bonds).
func alkaneChain(waypoints []Vec3, heavy []heavyAtom) [][]Vec3 {
	const (
		alkBond = 1.54
		alkCosA = 0.8165
		alkSinA = 0.578
	)

	// Nearest-neighbour path.
	visited := make([]bool, len(waypoints))
	path := []int{0}
	visited[0] = true
	for len(path) < len(waypoints) {
		last := waypoints[path[len(path)-1]]
		bestIdx, bestD := -1, 999.0
		for i, p := range waypoints {
			if visited[i] {
				continue
			}
			if d := last.Sub(p).Norm(); d < bestD {
				bestD, bestIdx = d, i
			}
		}
		path = append(path, bestIdx)
		visited[bestIdx] = true
	}

	// Build one segment at a time, carrying prevDir for torsion continuity.
	buildSeg := func(tip, target Vec3, n int, prevDir Vec3) ([]Vec3, Vec3, float64) {
		axis := target.Sub(tip).unit()
		segLen := tip.Sub(target).Norm()

		basePerp := func() Vec3 {
			if prevDir.Norm() > 0.1 {
				comp := prevDir.Sub(axis.Scale(prevDir.Dot(axis)))
				if comp.Norm() > 0.3 {
					return comp.Scale(-1).unit()
				}
			}
			return perpendicular(axis)
		}()
		perp2 := cross3(axis, basePerp).unit()

		type pt2 struct{ u, v float64 }
		local := make([]pt2, n+2)
		for k := 1; k <= n+1; k++ {
			sign := 1.0
			if (k-1)%2 == 1 {
				sign = -1
			}
			local[k] = pt2{local[k-1].u + alkCosA*alkBond, local[k-1].v + sign*alkSinA*alkBond}
		}
		eu, ev := local[n+1].u, local[n+1].v
		idealLen := math.Sqrt(eu*eu + ev*ev)
		s := segLen / idealLen

		buildWith := func(perpV Vec3) []Vec3 {
			if s < 0.75 || s > 1.25 {
				return nil
			}
			cosθ, sinθ := eu/idealLen, ev/idealLen
			pts := make([]Vec3, n)
			for k := 1; k <= n; k++ {
				u, v := local[k].u, local[k].v
				ru := (u*cosθ + v*sinθ) * s
				rv := (-u*sinθ + v*cosθ) * s
				pts[k-1] = tip.Add(axis.Scale(ru)).Add(perpV.Scale(rv))
			}
			return pts
		}
		clashCt := func(pts []Vec3) int {
			c := 0
			for _, p := range pts {
				for _, h := range heavy {
					if p.Sub(h.pos).Norm() < h.vdwR+vdw("C")-hardTol {
						c++
						break
					}
				}
			}
			return c
		}
		juncDev := func(pts []Vec3) float64 {
			if len(pts) < 2 {
				return 0
			}
			v1 := pts[len(pts)-2].Sub(pts[len(pts)-1])
			v2 := target.Sub(pts[len(pts)-1])
			n1, n2 := v1.Norm(), v2.Norm()
			if n1 < 1e-6 || n2 < 1e-6 {
				return 0
			}
			c := v1.Dot(v2) / (n1 * n2)
			if c > 1 {
				c = 1
			}
			if c < -1 {
				c = -1
			}
			return math.Abs(180/math.Pi*math.Acos(c) - 109.5)
		}

		best := buildWith(basePerp)
		if best == nil {
			best = []Vec3{}
		}
		bestC := clashCt(best)
		bestAD := juncDev(best)
		for oi := 1; oi < 36; oi++ {
			ang := float64(oi) * math.Pi / 18
			pv := basePerp.Scale(math.Cos(ang)).Add(perp2.Scale(math.Sin(ang)))
			pts := buildWith(pv)
			if pts == nil {
				continue
			}
			c := clashCt(pts)
			ad := juncDev(pts)
			if c < bestC || (c == bestC && ad < bestAD) {
				bestC, bestAD, best = c, ad, pts
			}
		}
		var newDir Vec3
		if len(best) >= 2 {
			newDir = best[len(best)-1].Sub(best[len(best)-2]).unit()
		} else if len(best) == 1 {
			newDir = best[0].Sub(tip).unit()
		} else {
			newDir = axis
		}
		return best, newDir, bestAD
	}

	var segments [][]Vec3
	var currentSeg []Vec3
	currentSeg = append(currentSeg, waypoints[path[0]])
	prevDir := Vec3{}

	flushSeg := func() {
		if len(currentSeg) > 1 {
			segments = append(segments, currentSeg)
		}
		currentSeg = nil
	}

	for wi := 1; wi < len(path); wi++ {
		tip := currentSeg[len(currentSeg)-1]
		target := waypoints[path[wi]]
		segLen := tip.Sub(target).Norm()

		n := int(math.Round(segLen/1.258)) - 1
		if n < 2 {
			n = 2
		}
		bestN, bestDev := n, 999.0
		for _, tryN := range []int{n - 1, n, n + 1, n + 2} {
			if tryN < 2 {
				continue
			}
			type pt2 struct{ u, v float64 }
			local := make([]pt2, tryN+2)
			for k := 1; k <= tryN+1; k++ {
				sign := 1.0
				if (k-1)%2 == 1 {
					sign = -1
				}
				local[k] = pt2{local[k-1].u + alkCosA*alkBond, local[k-1].v + sign*alkSinA*alkBond}
			}
			eu, ev := local[tryN+1].u, local[tryN+1].v
			idealLen := math.Sqrt(eu*eu + ev*ev)
			if dev := math.Abs(segLen/idealLen - 1.0); dev < bestDev {
				bestDev, bestN = dev, tryN
			}
		}

		pts, newDir, ad := buildSeg(tip, target, bestN, prevDir)

		if (ad > 45 && len(pts) > 0) || len(pts) == 0 {
			// Bad junction or degenerate segment: flush the current segment
			// and start a fresh one from the target waypoint.  This prevents
			// any bond being written across the discontinuity.
			flushSeg()
			currentSeg = append(currentSeg, target)
			prevDir = Vec3{}
			continue
		}

		for _, p := range pts {
			currentSeg = append(currentSeg, p)
		}
		prevDir = newDir
	}
	flushSeg()
	return segments
}

// ── Aromatic ring probes ──────────────────────────────────────────────────────

func buildAromaticProbes(proteinAtoms []Atom, heavy []heavyAtom, ps *ProbeSet) {
	resAtoms := map[string][]Vec3{}
	resOrder := []string{}
	for _, a := range proteinAtoms {
		if isAromatic(a) {
			k := resKey(a)
			if _, ok := resAtoms[k]; !ok {
				resOrder = append(resOrder, k)
			}
			resAtoms[k] = append(resAtoms[k], a.Pos)
		}
	}
	fmt.Printf("Aromatic residues (PHE/TYR/TRP): %d\n", len(resOrder))

	buildHeavyExcluding := func(rk string) []heavyAtom {
		var out []heavyAtom
		for _, a := range proteinAtoms {
			if strings.ToUpper(a.Element) == "H" || isWater(a) {
				continue
			}
			if resKey(a) == rk {
				continue
			}
			out = append(out, heavyAtom{a.Pos, vdw(a.Element)})
		}
		return out
	}

	nParallel, nTshaped := 0, 0
	ringCentres := []Vec3{}

	for _, rk := range resOrder {
		pts := resAtoms[rk]
		if len(pts) < 3 {
			continue
		}
		centroid := Vec3{}
		for _, p := range pts {
			centroid = centroid.Add(p)
		}
		centroid = centroid.Scale(1.0 / float64(len(pts)))
		normal := ringPlaneNormal(pts)
		excl := buildHeavyExcluding(rk)

		tooClose := func(pos Vec3) bool {
			for _, c := range ringCentres {
				if pos.Sub(c).Norm() < 4.5 {
					return true
				}
			}
			return false
		}

		addRing := func(centre, norm Vec3) bool {
			ringPts, edges := benzeneRing(centre, norm)
			for _, p := range ringPts {
				if !noHardClash(p, vdw("C"), hardTol, excl) {
					return false
				}
			}
			base := ps.Len() + 1
			for _, p := range ringPts {
				ps.Add(p, "C")
			}
			for _, e := range edges {
				ps.Bond(base+e[0], base+e[1], 4)
			}
			// Methyl cap on ring atom 3 (opposite atom 0), pointing away from
			// the source ring centroid so the cap faces into solvent.
			capDir := ringPts[3].Sub(centroid)
			if capDir.Norm() < 1e-6 {
				capDir = norm
			}
			capPos := ringPts[3].Add(capDir.unit().Scale(1.54))
			capIdx := ps.Len() + 1
			ps.Add(capPos, "C")
			ps.Bond(base+3, capIdx, 1) // ring C – methyl cap
			return true
		}

		// Parallel stacking.
		for _, sign := range []float64{1, -1} {
			centre := centroid.Add(normal.Scale(sign * 3.5))
			if tooClose(centre) {
				continue
			}
			if addRing(centre, normal) {
				ringCentres = append(ringCentres, centre)
				nParallel++
				break
			}
		}

		// T-shaped stacking.
		for _, sign := range []float64{1, -1} {
			centre := centroid.Add(normal.Scale(sign * 5.0))
			if tooClose(centre) {
				continue
			}
			tNorm := perpendicular(normal)
			if addRing(centre, tNorm) {
				ringCentres = append(ringCentres, centre)
				nTshaped++
				break
			}
		}
	}
	fmt.Printf("Benzene rings: %d parallel, %d T-shaped\n", nParallel, nTshaped)
}

// ── H-bond vote probes ────────────────────────────────────────────────────────

func buildHbondProbes(proteinAtoms []Atom, heavy []heavyAtom, ps *ProbeSet) {
	var donors, acceptors []Atom
	for _, a := range proteinAtoms {
		role := hbondRole(a)
		if role == "donor" || role == "dual" {
			donors = append(donors, a)
		}
		if role == "acceptor" || role == "dual" {
			acceptors = append(acceptors, a)
		}
	}
	fmt.Printf("H-bond donors: %d  Acceptors: %d\n", len(donors), len(acceptors))

	placeHbondProbes := func(group []Atom, probeElem string) int {
		probeR := vdw(probeElem)
		votes := castVotes(group, hbondSphereR, hbondPairDist, nSamples, heavy, probeElem)
		type candidate struct {
			pos       Vec3
			pairCount int
		}
		var cands []candidate
		for _, v := range votes {
			if v.pairCount >= hbondMinPairs {
				cands = append(cands, candidate{v.pos, v.pairCount})
			}
		}
		sort.Slice(cands, func(i, j int) bool { return cands[i].pairCount > cands[j].pairCount })
		placed := 0
		var placed_pos []Vec3
		for _, c := range cands {
			if !noHardClash(c.pos, probeR, hardTol, heavy) {
				continue
			}
			tooClose := false
			for _, pp := range placed_pos {
				if c.pos.Sub(pp).Norm() < hbondProbeSpacing {
					tooClose = true
					break
				}
			}
			if tooClose {
				continue
			}
			placed_pos = append(placed_pos, c.pos)
			ps.Add(c.pos, probeElem)
			placed++
		}
		return placed
	}

	nAcc := placeHbondProbes(donors, "O")
	fmt.Printf("H-bond acceptor probes (O): %d\n", nAcc)
	nDon := placeHbondProbes(acceptors, "N")
	fmt.Printf("H-bond donor probes (N): %d\n", nDon)
}
