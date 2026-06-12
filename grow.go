package main

import (
	"fmt"
	"math"
	"strings"
)

// substituent describes one fragment to place at an implicit H position.
type substituent struct {
	name  string
	atoms []string  // element labels (first bonded to anchor)
	blen  []float64 // bond lengths: anchor→a0, a0→a1, ...
	bords []int     // bond orders within fragment
	geom  []string  // "linear"|"sp2"|"sp3" for each bond
}

// substituentsFor returns the substituent palette for anchor elem/hybridization.
func substituentsFor(elem string, isSp2 bool) []substituent {
	switch strings.ToUpper(elem) {
	case "C":
		if isSp2 {
			return []substituent{
				{name: "methyl",   atoms: []string{"C"},       blen: []float64{1.51},      geom: []string{"sp2"}},
				{name: "OH",       atoms: []string{"O"},       blen: []float64{1.36},      geom: []string{"sp2"}},
				{name: "NH2",      atoms: []string{"N"},       blen: []float64{1.40},      geom: []string{"sp2"}},
				{name: "SH",       atoms: []string{"S"},       blen: []float64{1.77},      geom: []string{"sp2"}},
				{name: "F",        atoms: []string{"F"},       blen: []float64{1.35},      geom: []string{"sp2"}},
				{name: "Cl",       atoms: []string{"Cl"},      blen: []float64{1.74},      geom: []string{"sp2"}},
				{name: "OMe",      atoms: []string{"O","C"},   blen: []float64{1.36,1.43}, geom: []string{"sp2","sp3"}},
				{name: "NHMe",     atoms: []string{"N","C"},   blen: []float64{1.40,1.47}, geom: []string{"sp2","sp3"}},
				{name: "SMe",      atoms: []string{"S","C"},   blen: []float64{1.77,1.82}, geom: []string{"sp2","sp3"}},
				{name: "vinyl",    atoms: []string{"C","C"},   blen: []float64{1.34,1.34}, bords:[]int{2}, geom: []string{"sp2","sp2"}},
				{name: "ethynyl",  atoms: []string{"C","C"},   blen: []float64{1.47,1.20}, bords:[]int{3}, geom: []string{"sp2","linear"}},
				{name: "CN",       atoms: []string{"C","N"},   blen: []float64{1.47,1.16}, bords:[]int{3}, geom: []string{"sp2","linear"}},
				{name: "CHO",      atoms: []string{"C","O"},   blen: []float64{1.50,1.21}, bords:[]int{2}, geom: []string{"sp2","sp2"}},
				{name: "naphthyl", atoms: []string{"C","C"},   blen: []float64{1.40,1.40}, bords:[]int{4}, geom: []string{"sp2","sp2"}},
			}
		}
		return []substituent{
			{name: "methyl",    atoms: []string{"C"},         blen: []float64{1.54},           geom: []string{"sp3"}},
			{name: "ethyl",     atoms: []string{"C","C"},     blen: []float64{1.54,1.54},      geom: []string{"sp3","sp3"}},
			{name: "isobutyl",  atoms: []string{"C","C","C","C"}, blen: []float64{1.54,1.54,1.54,1.54}, geom: []string{"sp3","sp3","sp3","sp3"}},
			{name: "OH",        atoms: []string{"O"},         blen: []float64{1.43},           geom: []string{"sp3"}},
			{name: "NH2",       atoms: []string{"N"},         blen: []float64{1.47},           geom: []string{"sp3"}},
			{name: "SH",        atoms: []string{"S"},         blen: []float64{1.82},           geom: []string{"sp3"}},
			{name: "F",         atoms: []string{"F"},         blen: []float64{1.39},           geom: []string{"sp3"}},
			{name: "Cl",        atoms: []string{"Cl"},        blen: []float64{1.79},           geom: []string{"sp3"}},
			{name: "OMe",       atoms: []string{"O","C"},     blen: []float64{1.43,1.43},      geom: []string{"sp3","sp3"}},
			{name: "NHMe",      atoms: []string{"N","C"},     blen: []float64{1.47,1.47},      geom: []string{"sp3","sp3"}},
			{name: "SMe",       atoms: []string{"S","C"},     blen: []float64{1.82,1.82},      geom: []string{"sp3","sp3"}},
			{name: "NCH2CH2N",  atoms: []string{"N","C","C","N"}, blen: []float64{1.47,1.54,1.54,1.47}, geom: []string{"sp3","sp3","sp3","sp3"}},
		}
	case "N":
		return []substituent{
			{name: "methyl", atoms: []string{"C"},         blen: []float64{1.47},           geom: []string{"sp3"}},
			{name: "acetyl", atoms: []string{"C","C","O"}, blen: []float64{1.36,1.52,1.21}, bords:[]int{0,2}, geom: []string{"sp2","sp3","sp2"}},
		}
	case "O":
		return []substituent{
			{name: "methyl", atoms: []string{"C"},     blen: []float64{1.43},      geom: []string{"sp3"}},
			{name: "ethyl",  atoms: []string{"C","C"}, blen: []float64{1.43,1.54}, geom: []string{"sp3","sp3"}},
		}
	case "S":
		return []substituent{
			{name: "methyl", atoms: []string{"C"}, blen: []float64{1.82}, geom: []string{"sp3"}},
		}
	}
	return nil
}

// closureBondOrder returns the bond order for a closure bond based on distance
// and the elements involved. O and S are capped at 2; N is capped at 2 as well
// since triple-bond N in ring closures is chemically unrealistic here.
func closureBondOrder(d float64, elemI, elemJ string) int {
	ord := 1
	if d <= 1.28 {
		ord = 3
	} else if d <= 1.44 {
		ord = 2
	}
	// Cap bond order by the lower of the two elements' realistic maxima.
	maxOrd := func(el string) int {
		switch strings.ToUpper(el) {
		case "O", "S", "N":
			return 2
		}
		return 3
	}
	if cap := maxOrd(elemI); ord > cap {
		ord = cap
	}
	if cap := maxOrd(elemJ); ord > cap {
		ord = cap
	}
	return ord
}

// GrowLinked performs depth-first growth from every implicit H on every atom
// of each linked pair. At each H position it tries placing C, O, N, S, F, Cl.
// The atom with the highest burial score (protein contacts within 5 Å) is kept.
// If it is C, the DFS recurses onto that carbon's tetrahedral H positions.
// Growth stops when: (a) no atom scores better than solvent, (b) max depth, or
// (c) the new atom comes within bonding distance of an existing pair atom
// (ring closure) — in which case a closing bond is added.
func GrowLinked(pairs []ProbeSet, heavy []heavyAtom) []ProbeSet {
	const (
		clashTol = 0.8  // protein clash tolerance
		maxDepth = 4    // maximum recursive growth depth
	)

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

	selfClashFree := func(pos Vec3, existingPos []Vec3, bondedIdx int) bool {
		r := vdw("C") // approximate
		for i, ep := range existingPos {
			if i == bondedIdx { continue }
			if pos.Sub(ep).Norm() < r+r-1.0 {
				return false
			}
		}
		return true
	}

	// closureCheck: is pos within bonding distance of any existing atom that is
	// NOT the immediate anchor? Returns (atomIdx, bondOrder) or (-1, 0).
	closureCheck := func(pos Vec3, newElem string, existingPos []Vec3, existingLabels []string, anchorIdx int) (int, int) {
		for i, ep := range existingPos {
			if i == anchorIdx { continue }
			d := pos.Sub(ep).Norm()
			if d >= 1.10 && d <= 1.65 {
				el := ""
				if i < len(existingLabels) { el = existingLabels[i] }
				return i, closureBondOrder(d, newElem, el)
			}
		}
		return -1, 0
	}

	maxVal := map[string]int{"C": 4, "N": 3, "O": 2, "S": 2, "F": 1, "CL": 1}

	// candidates to try at each H position
	candidates := []struct{ elem string; blen float64 }{
		{"C",  1.54}, {"O", 1.43}, {"N", 1.47},
		{"S",  1.82}, {"F", 1.39}, {"Cl", 1.74},
	}

	var results []ProbeSet

	// dfsGrow grows the molecule by placing one atom at a time via DFS.
	// grown: the current ProbeSet being built (copy of original pair + grown atoms)
	// anchorIdx: 0-based index in grown.Pos of the atom we're growing from
	// hDir: direction of the implicit H being filled
	// depth: current recursion depth
	// pairIdx: parent pair index for ParentIdx tagging
	var dfsGrow func(grown ProbeSet, anchorIdx int, hDir Vec3, depth int, pairIdx int)
	dfsGrow = func(grown ProbeSet, anchorIdx int, hDir Vec3, depth int, pairIdx int) {
		if depth > maxDepth { return }

		anchor := grown.Pos[anchorIdx]

		// Try each candidate element.
		bestScore := -1
		bestElem := ""
		bestPos := Vec3{}
		bestClosure := -1
		bestCloseOrder := 0

		for _, cand := range candidates {
			pos := anchor.Add(hDir.Scale(cand.blen))

			if !clashFreeProtein(pos, cand.elem) { continue }
			if !selfClashFree(pos, grown.Pos, anchorIdx) { continue }

			// Check for ring closure.
			closeIdx, closeOrd := closureCheck(pos, cand.elem, grown.Pos, grown.Labels, anchorIdx)

			score := burialScore(pos)
			if closeIdx >= 0 { score += 10 } // bonus for closure

			if score > bestScore {
				bestScore = score
				bestElem = cand.elem
				bestPos = pos
				bestClosure = closeIdx
				bestCloseOrder = closeOrd
			}
		}

		if bestElem == "" { return }

		// Add the best atom.
		newIdx := grown.Add(bestPos, bestElem) - 1 // 0-based index
		grown.Bond(anchorIdx+1, newIdx+1, 1)

		if bestClosure >= 0 {
			grown.Bond(newIdx+1, bestClosure+1, bestCloseOrder)
			grown.Name = fmt.Sprintf("dfs-ring%d", bestCloseOrder)
			results = append(results, grown)
			return
		}

		// If it's C, recurse greedily on its best H position (no branching).
		if strings.ToUpper(bestElem) == "C" {
			// Compute tetrahedral H directions from the new C.
			hDirs2 := implicitHDirs(bestPos, []int{anchorIdx}, grown.Pos[:newIdx], false, 3)
			// Pick the direction with best burial — recurse on just that one.
			bestH := Vec3{}
			bestHScore := -1
			for _, hd2 := range hDirs2 {
				testPos := bestPos.Add(hd2.Scale(1.54))
				if !clashFreeProtein(testPos, "C") { continue }
				s := burialScore(testPos)
				if s > bestHScore { bestHScore = s; bestH = hd2 }
			}
			if bestHScore >= 0 {
				dfsGrow(grown, newIdx, bestH, depth+1, pairIdx)
				return // already emitted in recursive call
			}
		}

		// Emit this grown state as a variant (even if not a ring closure).
		grown.Name = fmt.Sprintf("dfs-grow-%s%d", bestElem, depth)
		grown.ParentIdx = pairIdx
		results = append(results, grown)
	}

	for pairIdx, pair := range pairs {
		n := len(pair.Pos)
		adj := make([][]int, n)
		bondSum := make([]float64, n)
		for _, b := range pair.Bonds {
			i, j := b.I-1, b.J-1
			if i < 0 || j < 0 || i >= n || j >= n { continue }
			adj[i] = append(adj[i], j)
			adj[j] = append(adj[j], i)
			v := float64(b.Order)
			if b.Order == 4 { v = 1.5 }
			bondSum[i] += v; bondSum[j] += v
		}

		for i := 0; i < n; i++ {
			elem := strings.ToUpper(pair.Labels[i])
			mv, ok := maxVal[elem]
			if !ok { continue }
			implicitH := mv - int(math.Round(bondSum[i]))
			if implicitH <= 0 { continue }

			// Skip pure linker carbons (>2 bonds from any heteroatom).
			if elem == "C" {
				hasDbl, hasArom := false, false
				for _, b := range pair.Bonds {
					if (b.I-1 == i || b.J-1 == i) {
						if b.Order == 2 { hasDbl = true }
						if b.Order == 4 { hasArom = true }
					}
				}
				if !hasDbl && !hasArom {
					nearHet := false
					for _, nb := range adj[i] {
						el := strings.ToUpper(pair.Labels[nb])
						if el == "N" || el == "O" || el == "S" { nearHet = true; break }
						for _, nb2 := range adj[nb] {
							el2 := strings.ToUpper(pair.Labels[nb2])
							if el2 == "N" || el2 == "O" || el2 == "S" { nearHet = true; break }
						}
						if nearHet { break }
					}
					if !nearHet { continue }
				}
			}

			isSp2 := false
			for _, b := range pair.Bonds {
				if (b.I-1 == i || b.J-1 == i) && (b.Order == 2 || b.Order == 4) {
					isSp2 = true; break
				}
			}
			// Also treat N as sp2 if any neighbour has a double/aromatic bond
			// (e.g. guanidine N: N–C(=N) → the N is sp2 even though its own bond is single).
			if !isSp2 && elem == "N" {
				for _, nb := range adj[i] {
					for _, b := range pair.Bonds {
						if (b.I-1 == nb || b.J-1 == nb) && (b.Order == 2 || b.Order == 4) {
							isSp2 = true; break
						}
					}
					if isSp2 { break }
				}
			}

			hDirs := implicitHDirs(pair.Pos[i], adj[i], pair.Pos, isSp2, implicitH)
			for _, hd := range hDirs {
				// Start a fresh grown copy for each H position.
				// Preserve the original pharmPair index through ring-scan layers.
				parentIdx := pairIdx
				if pair.ParentIdx > 0 { parentIdx = pair.ParentIdx }
				grown := ProbeSet{ParentIdx: parentIdx}
				for k, p := range pair.Pos { grown.Add(p, pair.Labels[k]) }
				for _, b := range pair.Bonds { grown.Bond(b.I, b.J, b.Order) }
				dfsGrow(grown, i, hd, 0, pairIdx)
			}
		}
	}

	return results
}

// RingClose attempts to form new rings within each linked pair.
func RingClose(pairs []ProbeSet, heavy []heavyAtom) []ProbeSet {
	_ = heavy
	maxVal := map[string]int{"C": 4, "N": 3, "O": 2, "S": 2}
	var results []ProbeSet

	for pairIdx, pair := range pairs {
		n := len(pair.Pos)
		adj := make([][]int, n)
		bondSum := make([]float64, n)
		bonded := make([][]bool, n)
		for i := range bonded { bonded[i] = make([]bool, n) }
		for _, b := range pair.Bonds {
			i, j := b.I-1, b.J-1
			if i < 0 || j < 0 || i >= n || j >= n { continue }
			adj[i] = append(adj[i], j); adj[j] = append(adj[j], i)
			bonded[i][j] = true; bonded[j][i] = true
			v := float64(b.Order)
			if b.Order == 4 { v = 1.5 }
			bondSum[i] += v; bondSum[j] += v
		}
		shortestPath := func(src, dst int) int {
			dist := make([]int, n)
			for i := range dist { dist[i] = -1 }
			dist[src] = 0
			q := []int{src}
			for len(q) > 0 {
				cur := q[0]; q = q[1:]
				if cur == dst { return dist[dst] }
				for _, nb := range adj[cur] {
					if dist[nb] < 0 { dist[nb] = dist[cur] + 1; q = append(q, nb) }
				}
			}
			return dist[dst]
		}

		type candidate struct{ i, j, ringSize, order int }
		var cands []candidate
		for i := 0; i < n; i++ {
			elI := strings.ToUpper(pair.Labels[i])
			mvI, ok := maxVal[elI]; if !ok { continue }
			if int(math.Round(bondSum[i])) >= mvI { continue }
			for j := i + 1; j < n; j++ {
				if bonded[i][j] { continue }
				elJ := strings.ToUpper(pair.Labels[j])
				mvJ, ok := maxVal[elJ]; if !ok { continue }
				if int(math.Round(bondSum[j])) >= mvJ { continue }
				d := pair.Pos[i].Sub(pair.Pos[j]).Norm()
				if d < 1.10 || d > 1.65 { continue }
				path := shortestPath(i, j)
				if path < 3 || path > 5 { continue } // 4, 5, or 6-membered rings
				cands = append(cands, candidate{i, j, path + 1, closureBondOrder(d, elI, elJ)})
			}
		}
		if len(cands) > 6 { cands = cands[:6] }
		for mask := 1; mask < (1 << len(cands)); mask++ {
			uses := map[int]int{}
			valid := true
			var sub []candidate
			for k, c := range cands {
				if mask&(1<<k) == 0 { continue }
				uses[c.i]++; uses[c.j]++
				if uses[c.i] > 1 || uses[c.j] > 1 { valid = false; break }
				sub = append(sub, c)
			}
			if !valid { continue }
			closed := ProbeSet{ParentIdx: pairIdx}
			for k, p := range pair.Pos { closed.Add(p, pair.Labels[k]) }
			for _, b := range pair.Bonds { closed.Bond(b.I, b.J, b.Order) }
			name := ""
			for _, c := range sub {
				closed.Bond(c.i+1, c.j+1, c.order)
				name += fmt.Sprintf("+ring%d-%s%d-%s%d", c.ringSize,
					pair.Labels[c.i], c.i+1, pair.Labels[c.j], c.j+1)
			}
			closed.Name = name[1:]
			results = append(results, closed)
		}
	}
	return results
}

// placeFragment places substituent atoms from anchor in hDir direction.
func placeFragment(anchor, hDir Vec3, sub substituent, anchorSp2 bool, clashFree func(Vec3, string) bool) ([]Vec3, bool) {
	fragPos := make([]Vec3, len(sub.atoms))
	_ = anchorSp2

	pos0 := anchor.Add(hDir.Scale(sub.blen[0]))
	if !clashFree(pos0, sub.atoms[0]) { return nil, false }
	fragPos[0] = pos0
	if len(sub.atoms) == 1 { return fragPos, true }

	for k := 1; k < len(sub.atoms); k++ {
		bl := sub.blen[k]
		geom := "sp3"
		if k < len(sub.geom) { geom = sub.geom[k] }
		curDir := fragPos[k-1].Sub(anchor).unit()
		if k > 1 { curDir = fragPos[k-1].Sub(fragPos[k-2]).unit() }
		basePerp := perpendicular(curDir)
		perp2 := cross3(curDir, basePerp)
		placed := false
		var nextPos Vec3
		switch geom {
		case "linear":
			nextPos = fragPos[k-1].Add(curDir.Scale(bl))
			if clashFree(nextPos, sub.atoms[k]) { placed = true }
		case "sp2":
			for oi := 0; oi < 12 && !placed; oi++ {
				ang := float64(oi) * math.Pi / 6
				inP := basePerp.Scale(math.Cos(ang)).Add(perp2.Scale(math.Sin(ang)))
				d := curDir.Scale(0.5).Add(inP.Scale(0.866))
				nextPos = fragPos[k-1].Add(d.unit().Scale(bl))
				if clashFree(nextPos, sub.atoms[k]) { placed = true }
			}
		default:
			for oi := 0; oi < 12 && !placed; oi++ {
				ang := float64(oi) * math.Pi / 6
				p := basePerp.Scale(math.Cos(ang)).Add(perp2.Scale(math.Sin(ang)))
				d := curDir.Scale(-0.333).Add(p.Scale(0.9428))
				nextPos = fragPos[k-1].Add(d.unit().Scale(bl))
				if clashFree(nextPos, sub.atoms[k]) { placed = true }
			}
		}
		if !placed { return nil, false }
		fragPos[k] = nextPos
	}
	return fragPos, true
}

// implicitHDirs computes implicit H directions for an atom.
func implicitHDirs(pos Vec3, nbIdxs []int, allPos []Vec3, isSp2 bool, nH int) []Vec3 {
	nbs := make([]Vec3, len(nbIdxs))
	for i, idx := range nbIdxs { nbs[i] = allPos[idx] }
	switch len(nbs) {
	case 0:
		return evenDirs(nH)
	case 1:
		axis := nbs[0].Sub(pos).unit()
		perp := perpendicular(axis); q := cross3(axis, perp)
		if isSp2 {
			var dirs []Vec3
			for k := 0; k < nH; k++ {
				ang := float64(k) * 2 * math.Pi / float64(nH+1)
				inP := perp.Scale(math.Cos(ang)).Add(q.Scale(math.Sin(ang)))
				dirs = append(dirs, axis.Scale(0.5).Add(inP.Scale(0.866)).unit())
			}
			return dirs
		}
		const cosA, sinA = 0.333, 0.9428
		var dirs []Vec3
		for k := 0; k < nH; k++ {
			ang := float64(k) * 2 * math.Pi / float64(nH)
			p := perp.Scale(math.Cos(ang)).Add(q.Scale(math.Sin(ang)))
			dirs = append(dirs, axis.Scale(-cosA).Add(p.Scale(sinA)).unit())
		}
		return dirs
	case 2:
		b1 := nbs[0].Sub(pos).unit(); b2 := nbs[1].Sub(pos).unit()
		hDir := b1.Add(b2).Scale(-1).unit()
		if nH == 1 { return []Vec3{hDir} }
		perp := cross3(b1, b2)
		if perp.Norm() < 1e-6 { perp = perpendicular(hDir) }
		perp = perp.unit()
		var dirs []Vec3
		for k := 0; k < nH; k++ {
			ang := float64(k) * 2 * math.Pi / float64(nH)
			dirs = append(dirs, hDir.Add(perp.Scale(math.Sin(ang)*0.3)).unit())
		}
		return dirs
	case 3:
		sum := Vec3{}
		for _, nb := range nbs { sum = sum.Add(nb.Sub(pos).unit()) }
		return []Vec3{sum.Scale(-1).unit()}
	}
	return nil
}

func evenDirs(n int) []Vec3 {
	if n == 1 { return []Vec3{{1, 0, 0}} }
	var dirs []Vec3
	for i := 0; i < n; i++ {
		θ := float64(i) * 2 * math.Pi / float64(n)
		dirs = append(dirs, Vec3{math.Cos(θ), math.Sin(θ), 0})
	}
	return dirs
}

