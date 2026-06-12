package main

import (
	"fmt"
	"math"
	"sort"
	"strings"
)

// ScanChain generates two classes of variants for every sp3 linker atom in
// each linked molecule:
//
//  1. Single-atom substitution (chain C only): replace C with N, O, or S.
//     One best variant per position (highest heteroatom burial), emitted only
//     when it scores strictly better than the original.
//
//  2. Double-bond functional group (C, N, S where valence allows): add a new
//     atom via a double bond (C=O carbonyl, C=NH imine, C=S thioketone, etc.).
//     One best variant per position per functional group type.
//     Prioritised for positions that are already heteroatoms (N, O, S) since
//     those chains benefit most from H-bond-capable functional groups.
//
// All results are sorted by heteroatom burial score descending.
func ScanChain(mols []ProbeSet, heavy []heavyAtom) []ProbeSet {
	const clashTol = 0.3

	clashFreeProtein := func(pos Vec3, elem string) bool {
		r := vdw(elem)
		for _, h := range heavy {
			if pos.Sub(h.pos).Norm() < h.vdwR+r-clashTol {
				return false
			}
		}
		return true
	}

	burialScore := func(pos Vec3) int {
		n := 0
		for _, h := range heavy {
			d := pos.Sub(h.pos).Norm()
			if d >= 1.5 && d <= 5.0 {
				n++
			}
		}
		return n
	}

	hetBurialScore := func(mol ProbeSet) int {
		s := 0
		for i, p := range mol.Pos {
			switch strings.ToUpper(mol.Labels[i]) {
			case "N", "O", "S":
				s += burialScore(p)
			}
		}
		return s
	}

	maxValence := map[string]int{"C": 4, "N": 3, "O": 2, "S": 2}

	// sp2Dirs returns candidate directions for placing a double-bond substituent
	// on the atom at pos, given its existing single-bond neighbour indices.
	sp2Dirs := func(pos Vec3, nbIdxs []int, allPos []Vec3) []Vec3 {
		switch len(nbIdxs) {
		case 0:
			return evenDirs(4)
		case 1:
			// Terminal sp2: try a cone of 4 directions at 120° from the bond axis.
			axis := allPos[nbIdxs[0]].Sub(pos).unit()
			perp := perpendicular(axis)
			perp2 := cross3(axis, perp).unit()
			dirs := make([]Vec3, 4)
			for k := range dirs {
				ang := float64(k) * math.Pi / 2
				p := perp.Scale(math.Cos(ang)).Add(perp2.Scale(math.Sin(ang)))
				// 120° from axis: component along axis = -0.5, perp = 0.866
				dirs[k] = axis.Scale(-0.5).Add(p.Scale(0.866)).unit()
			}
			return dirs
		default:
			// Interior sp2: opposite the bisector of all existing bonds.
			sum := Vec3{}
			for _, nb := range nbIdxs {
				sum = sum.Add(allPos[nb].Sub(pos).unit())
			}
			if sum.Norm() < 1e-6 {
				// Near-linear chain: try perpendicular to first bond.
				axis := allPos[nbIdxs[0]].Sub(pos).unit()
				n := perpendicular(axis)
				n2 := cross3(axis, n).unit()
				return []Vec3{n, n.Scale(-1), n2, n2.Scale(-1)}
			}
			return []Vec3{sum.Scale(-1).unit()}
		}
	}

	// dblBondSub describes a double-bond functional group to try.
	type dblBondSub struct {
		newElem string
		blen    float64
		name    string
	}
	dblSubsFor := func(anchorElem string) []dblBondSub {
		switch strings.ToUpper(anchorElem) {
		case "C":
			return []dblBondSub{
				{"O", 1.21, "ketone"},
				{"N", 1.27, "imine"},
			}
		case "N":
			return []dblBondSub{
				{"O", 1.21, "nitroso"},
				{"C", 1.27, "imine-N"},
			}
		case "S":
			return []dblBondSub{
				{"O", 1.44, "sulfoxide"},
			}
		}
		return nil
	}

	var results []ProbeSet

	for molIdx, mol := range mols {
		n := len(mol.Pos)

		bondOrderSum := make([]int, n)
		hasAromBond := make([]bool, n)
		hasDblBond := make([]bool, n)
		adj := make([][]int, n)

		for _, b := range mol.Bonds {
			i, j := b.I-1, b.J-1
			if i < 0 || j < 0 || i >= n || j >= n {
				continue
			}
			adj[i] = append(adj[i], j)
			adj[j] = append(adj[j], i)
			ord := b.Order
			if ord == 4 {
				hasAromBond[i] = true
				hasAromBond[j] = true
				ord = 1
			}
			if b.Order == 2 {
				hasDblBond[i] = true
				hasDblBond[j] = true
			}
			bondOrderSum[i] += ord
			bondOrderSum[j] += ord
		}

		baseScore := hetBurialScore(mol)

		for i := 0; i < n; i++ {
			elem := strings.ToUpper(mol.Labels[i])

			// Only consider sp3 chain atoms: no aromatic or double bonds on this atom.
			if hasAromBond[i] || hasDblBond[i] {
				continue
			}
			// Linker atoms: degree ≤ 2 (interior chain or terminal).
			if bondOrderSum[i] > 2 {
				continue
			}
			// Skip if adjacent to an aromatic ring.
			adjArom := false
			for _, nb := range adj[i] {
				if hasAromBond[nb] {
					adjArom = true
					break
				}
			}
			if adjArom {
				continue
			}

			pos := mol.Pos[i]

			// ── Pass 1: single-atom substitution (C → N / O / S) ─────────────────
			if elem == "C" {
				bestScore := baseScore
				var bestVariant *ProbeSet

				for _, sub := range []string{"N", "O", "S"} {
					mv := maxValence[sub]
					if bondOrderSum[i] > mv {
						continue
					}
					if !clashFreeProtein(pos, sub) {
						continue
					}
					variant := ProbeSet{
						Name:      fmt.Sprintf("chain-pos%d-%s", i+1, sub),
						ParentIdx: molIdx,
					}
					for k, p := range mol.Pos {
						lbl := mol.Labels[k]
						if k == i {
							lbl = sub
						}
						variant.Add(p, lbl)
					}
					for _, b := range mol.Bonds {
						variant.Bond(b.I, b.J, b.Order)
					}
					score := hetBurialScore(variant)
					if score > bestScore {
						bestScore = score
						v := variant
						bestVariant = &v
					}
				}
				if bestVariant != nil {
					results = append(results, *bestVariant)
				}
			}

			// ── Pass 2: double-bond functional group addition ─────────────────────
			// Applies to C, N, S chain atoms where remaining valence ≥ 2.
			if noCarbonyl {
				continue
			}
			mv, ok := maxValence[elem]
			if !ok || mv-bondOrderSum[i] < 2 {
				continue
			}

			// Give a burial bonus to positions that are already heteroatoms,
			// matching the user preference for N/O/S chain atoms.
			isHetAtom := elem == "N" || elem == "O" || elem == "S"
			hetBonus := 0
			if isHetAtom {
				hetBonus = 3 // favour heteroatom positions in tie-breaking
			}

			for _, dsub := range dblSubsFor(elem) {
				dirs := sp2Dirs(pos, adj[i], mol.Pos)
				bestScore := baseScore + hetBonus - 1 // must beat base (adjusted for bonus)
				var bestVariant *ProbeSet

				for _, dir := range dirs {
					newPos := pos.Add(dir.Scale(dsub.blen))
					if !clashFreeProtein(newPos, dsub.newElem) {
						continue
					}
					variant := ProbeSet{
						Name:      fmt.Sprintf("chain-pos%d-%s", i+1, dsub.name),
						ParentIdx: molIdx,
					}
					for k, p := range mol.Pos {
						variant.Add(p, mol.Labels[k])
					}
					for _, b := range mol.Bonds {
						variant.Bond(b.I, b.J, b.Order)
					}
					newIdx := variant.Add(newPos, dsub.newElem)
					variant.Bond(i+1, newIdx, 2) // double bond to new atom

					score := hetBurialScore(variant) + hetBonus
					if score > bestScore {
						bestScore = score
						v := variant
						bestVariant = &v
					}
				}
				if bestVariant != nil {
					results = append(results, *bestVariant)
				}
			}
		}
	}

	sort.Slice(results, func(i, j int) bool {
		return hetBurialScore(results[i]) > hetBurialScore(results[j])
	})

	return results
}

// PolarFG describes a Y-shaped polar functional group that can be appended to a
// chain atom: anchor → central(=DblEl, -SngEl).  SngEl=="" means no single-bond arm.
type PolarFG struct {
	Name      string
	AnchorEl  string  // required anchor element ("C","N","" = any)
	CentralEl string
	CentralBl float64
	DblEl     string
	DblBl     float64
	SngEl     string
	SngBl     float64
}

// polarFGList is the palette of polar functional groups added to linker chains.
var polarFGList = []PolarFG{
	{"ketone",    "C", "C", 1.52, "O", 1.22, "",  0},
	{"amide",     "C", "C", 1.52, "O", 1.22, "N", 1.33},
	{"urea",      "N", "C", 1.36, "O", 1.22, "N", 1.33},
	{"guanidine", "N", "C", 1.36, "N", 1.29, "N", 1.34},
	{"thiourea",  "N", "C", 1.36, "S", 1.71, "N", 1.36},
}

// placeFG attempts to place a PolarFG at anchor+hDir, trying 12 rotations in the
// sp2 plane. Returns the three new atom positions and true on success.
func placeFG(fg PolarFG, anchor, hDir Vec3, clashFree func(Vec3, string) bool) (cPos, dblPos, sngPos Vec3, ok bool) {
	cPos = anchor.Add(hDir.Scale(fg.CentralBl))
	if !clashFree(cPos, fg.CentralEl) {
		return
	}
	axis := hDir
	perp := perpendicular(axis)
	perp2 := cross3(axis, perp).unit()
	for oi := 0; oi < 12; oi++ {
		ang := float64(oi) * math.Pi / 6
		inPlane := perp.Scale(math.Cos(ang)).Add(perp2.Scale(math.Sin(ang)))
		d1 := axis.Scale(-0.5).Add(inPlane.Scale(0.866)).unit()
		dblPos = cPos.Add(d1.Scale(fg.DblBl))
		if !clashFree(dblPos, fg.DblEl) {
			continue
		}
		if fg.SngEl == "" {
			return cPos, dblPos, Vec3{}, true
		}
		d2 := axis.Scale(-0.5).Add(inPlane.Scale(-0.866)).unit()
		sngPos = cPos.Add(d2.Scale(fg.SngBl))
		if !clashFree(sngPos, fg.SngEl) {
			continue
		}
		return cPos, dblPos, sngPos, true
	}
	return
}

// CountPolarFG counts ketone/amide/urea/guanidine/thiourea occurrences in mol.
// Each sp2 carbon with =O, =N, or =S is examined and scored by its neighbours.
func CountPolarFG(mol ProbeSet) int {
	n := len(mol.Pos)
	type edge struct{ order int }
	adjMap := make([]map[int]int, n) // adjMap[i][j] = bond order
	for i := range adjMap {
		adjMap[i] = map[int]int{}
	}
	for _, b := range mol.Bonds {
		i, j := b.I-1, b.J-1
		if i < 0 || j < 0 || i >= n || j >= n {
			continue
		}
		adjMap[i][j] = b.Order
		adjMap[j][i] = b.Order
	}
	el := func(i int) string { return strings.ToUpper(mol.Labels[i]) }
	count := 0
	counted := make([]bool, n)
	for i := 0; i < n; i++ {
		if counted[i] || el(i) != "C" {
			continue
		}
		// Find double-bond heteroatom neighbours.
		dblO, dblN, dblS := false, false, false
		nN := 0
		for j, ord := range adjMap[i] {
			switch {
			case ord == 2 && el(j) == "O":
				dblO = true
			case ord == 2 && el(j) == "N":
				dblN = true
			case ord == 2 && el(j) == "S":
				dblS = true
			case el(j) == "N":
				nN++
			}
		}
		if dblO {
			if nN >= 2 {
				count++ // urea
			} else if nN >= 1 {
				count++ // amide
			} else {
				count++ // ketone
			}
			counted[i] = true
		} else if dblN {
			if nN >= 2 {
				count++ // guanidine
			} else {
				count++ // imine
			}
			counted[i] = true
		} else if dblS {
			if nN >= 2 {
				count++ // thiourea
			} else {
				count++ // thioketone
			}
			counted[i] = true
		}
	}
	return count
}

// tryAddPolarFG finds the best single polar FG addition to mol (highest burial
// of new atoms, clash-free). Returns the updated mol and true, or false.
// Returns immediately when noCarbonyl is set.
func tryAddPolarFG(mol ProbeSet, heavy []heavyAtom) (ProbeSet, bool) {
	if noCarbonyl {
		return mol, false
	}
	n := len(mol.Pos)
	bondSum := make([]int, n)
	adj := make([][]int, n)
	for _, b := range mol.Bonds {
		i, j := b.I-1, b.J-1
		if i < 0 || j < 0 || i >= n || j >= n {
			continue
		}
		adj[i] = append(adj[i], j)
		adj[j] = append(adj[j], i)
		ord := b.Order
		if ord == 4 {
			ord = 1
		}
		bondSum[i] += ord
		bondSum[j] += ord
	}

	clashFree := func(pos Vec3, elem string) bool {
		r := vdw(elem)
		for _, h := range heavy {
			if pos.Sub(h.pos).Norm() < h.vdwR+r-0.3 {
				return false
			}
		}
		for _, ep := range mol.Pos {
			if pos.Sub(ep).Norm() < r+vdw("C")-1.0 {
				return false
			}
		}
		return true
	}

	bur := func(pos Vec3) int {
		n := 0
		for _, h := range heavy {
			d := pos.Sub(h.pos).Norm()
			if d >= 1.5 && d <= 5.0 {
				n++
			}
		}
		return n
	}

	maxVal := map[string]int{"C": 4, "N": 3, "O": 2, "S": 2}

	bestScore := -1
	bestMol := ProbeSet{}

	for i := 0; i < n; i++ {
		elem := strings.ToUpper(mol.Labels[i])
		mv, ok := maxVal[elem]
		if !ok {
			continue
		}
		freeVal := mv - bondSum[i]
		if freeVal < 1 {
			continue
		}
		hDirs := implicitHDirs(mol.Pos[i], adj[i], mol.Pos, false, freeVal)
		for _, fg := range polarFGList {
			if fg.AnchorEl != "" && fg.AnchorEl != elem {
				continue
			}
			for _, hDir := range hDirs {
				cPos, dblPos, sngPos, placed := placeFG(fg, mol.Pos[i], hDir, clashFree)
				if !placed {
					continue
				}
				score := bur(cPos) + bur(dblPos)
				if fg.SngEl != "" {
					score += bur(sngPos)
				}
				if score > bestScore {
					bestScore = score
					bestMol = ProbeSet{Name: mol.Name, ParentIdx: mol.ParentIdx}
					for k, p := range mol.Pos {
						bestMol.Add(p, mol.Labels[k])
					}
					for _, b := range mol.Bonds {
						bestMol.Bond(b.I, b.J, b.Order)
					}
					cIdx := bestMol.Add(cPos, fg.CentralEl)
					dIdx := bestMol.Add(dblPos, fg.DblEl)
					bestMol.Bond(i+1, cIdx, 1)
					bestMol.Bond(cIdx, dIdx, 2)
					if fg.SngEl != "" {
						sIdx := bestMol.Add(sngPos, fg.SngEl)
						bestMol.Bond(cIdx, sIdx, 1)
					}
				}
				break // first valid direction per FG per position
			}
		}
	}

	if bestScore < 0 {
		return mol, false
	}
	return bestMol, true
}

// MaxPolarFGMol greedily decorates a linker molecule with as many polar
// functional groups (ketone/amide/urea/guanidine/thiourea) as possible,
// returning the maximally decorated variant.
func MaxPolarFGMol(mol ProbeSet, heavy []heavyAtom) ProbeSet {
	result := mol
	for iter := 0; iter < 20; iter++ {
		next, added := tryAddPolarFG(result, heavy)
		if !added {
			break
		}
		result = next
	}
	result.Name = fmt.Sprintf("triple-maxfg(%d)", CountPolarFG(result))
	return result
}
