package main

import (
	"fmt"
	"math"
	"strings"
)

// BidirGrow builds an all-trans sp3 zigzag chain between the two farthest-apart
// tail atoms of each linked pair, growing from both ends simultaneously toward
// the middle.  The chain is split at the midpoint: one half is attributed to
// growing from tail A, the other from tail B.
func BidirGrow(pharmPairs []ProbeSet, heavy []heavyAtom) []ProbeSet {
	const (
		bondLen  = 1.54
		clashTol = 0.8
	)

	clashFree := func(p Vec3) bool {
		for _, h := range heavy {
			if p.Sub(h.pos).Norm() < h.vdwR+vdw("C")-clashTol {
				return false
			}
		}
		return true
	}

	var results []ProbeSet

	for pairIdx, pair := range pharmPairs {
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
			bondSum[i] += v
			bondSum[j] += v
		}

		var growable []int
		for i := 0; i < n; i++ {
			if strings.ToUpper(pair.Labels[i]) != "C" { continue }
			if len(adj[i]) > 2 { continue } // branching point — skip
			sp3 := true
			for _, b := range pair.Bonds {
				if (b.I-1 == i || b.J-1 == i) && (b.Order == 2 || b.Order == 4) {
					sp3 = false; break
				}
			}
			if !sp3 { continue }
			growable = append(growable, i)
		}
		if len(growable) < 2 { continue }

		bfs := func(src int) []int {
			d := make([]int, n)
			for i := range d { d[i] = -1 }
			d[src] = 0
			q := []int{src}
			for len(q) > 0 {
				cur := q[0]; q = q[1:]
				for _, nb := range adj[cur] {
					if d[nb] < 0 { d[nb] = d[cur] + 1; q = append(q, nb) }
				}
			}
			return d
		}

		bestG, gA, gB := -1, 0, 0
		for ti := 0; ti < len(growable); ti++ {
			d := bfs(growable[ti])
			for tj := ti + 1; tj < len(growable); tj++ {
				if d[growable[tj]] > bestG { bestG, gA, gB = d[growable[tj]], growable[ti], growable[tj] }
			}
		}

		fmt.Printf("  pair%d: growable=%d bestG=%d\n", pairIdx, len(growable), bestG)
		if bestG < 3 { continue }

		posA := pair.Pos[gA]
		posB := pair.Pos[gB]
		gap := posA.Sub(posB).Norm()
		axis := posB.Sub(posA).unit()
		basePerp := perpendicular(axis)
		perp2 := cross3(axis, basePerp)

		// 2D all-trans zigzag chain geometry.
		// Build the chain using exact sp3 all-trans geometry.
		// Each atom is placed at: posA + k*axialStep*axis ± lateralStep*perp
		// where axialStep = cos(70.5°)*bondLen = 0.5133 Å per atom along axis
		// and lateralStep = sin(70.5°)*bondLen = 1.4519 Å lateral.
		// We need nAtoms such that nAtoms * axialStep ≈ gap.
		const axialStep = 0.5133 // Å per atom along axis (cos(70.5°) * 1.54)
		const latStep   = 1.4519 // Å lateral (sin(70.5°) * 1.54)

		nNeeded := int(math.Round(gap / axialStep))
		if nNeeded < 2 { nNeeded = 2 }
		if nNeeded > 10 { continue } // gap too large for a sensible bridge

		bestPts := []Vec3(nil)
		bestClash := math.MaxInt32
		bestEndDist := math.MaxFloat64

		for nI := nNeeded - 3; nI <= nNeeded + 3; nI++ {
			if nI < 1 { continue }
			for oi := 0; oi < 36; oi++ {
				ang := float64(oi) * math.Pi / 18
				perp := basePerp.Scale(math.Cos(ang)).Add(perp2.Scale(math.Sin(ang)))

				pts := make([]Vec3, nI)
				for k := 0; k < nI; k++ {
					// All-trans sp3: atom k is at axial offset + lateral offset
					// Lateral alternates: +latStep for k=0,2,4,... and 0 for k=1,3,5,...
					// This gives correct 1.54 Å bond lengths between consecutive atoms.
					axial := float64(k+1) * axialStep
					lat := 0.0
					if k%2 == 0 { lat = latStep }
					pts[k] = posA.Add(axis.Scale(axial)).Add(perp.Scale(lat))
				}

				// Check endpoint distance to posB.
				endDist := pts[nI-1].Sub(posB).Norm()
				if oi == 0 && nI == 19 { fmt.Printf("      nI=%d oi=%d endDist=%.3f pts[-1]=%v posB=%v\n", nI, oi, endDist, pts[nI-1], posB) }
				if endDist > bondLen*1.15 || endDist < bondLen*0.85 { continue }

				nc := 0
				for _, p := range pts {
					if !clashFree(p) { nc++ }
				}
				maxAllowed := nI / 5
				if oi==0 && nI==19 { fmt.Printf("      nc=%d maxAllowed=%d\n", nc, maxAllowed) }
					if nc > maxAllowed { continue }

				if nc < bestClash || (nc == bestClash && endDist < bestEndDist) {
					bestClash = nc
					bestEndDist = endDist
					bestPts = append([]Vec3{}, pts...)
				}
			}
		}

		if bestPts == nil {
			continue
		}

		// Split at midpoint: first half = A's contribution, second = B's.
		mid := len(bestPts) / 2
		chainA := bestPts[:mid]
		chainB := bestPts[mid:]

		grown := ProbeSet{ParentIdx: pairIdx}
		for k, p := range pair.Pos { grown.Add(p, pair.Labels[k]) }
		for _, b := range pair.Bonds { grown.Bond(b.I, b.J, b.Order) }

		bA := grown.Len() + 1
		for _, p := range chainA { grown.Add(p, "C") }
		bB := grown.Len() + 1
		for _, p := range chainB { grown.Add(p, "C") }

		// Bond gA → chainA[0] → chainA[1] → ... → chainA[mid-1]
		prev := gA + 1
		for k := range chainA { grown.Bond(prev, bA+k, 1); prev = bA + k }
		lastA := prev

		// Bond gB → chainB[last] → ... → chainB[0]
		prev = gB + 1
		for k := len(chainB) - 1; k >= 0; k-- { grown.Bond(prev, bB+k, 1); prev = bB + k }

		// Close the two tips. Pick bond order based on distance.
		midDist := bestPts[mid-1].Sub(bestPts[mid]).Norm()
		midOrder := 1
		if midDist < 1.30 { midOrder = 3 } else if midDist < 1.42 { midOrder = 2 }
		grown.Bond(lastA, prev, midOrder)

		grown.Name = fmt.Sprintf("bidir-bridge%d+%d", len(chainA), len(chainB))
		results = append(results, grown)
	}
	return results
}

// outwardDir returns the direction pointing away from a tail atom's neighbours.
func outwardDir(tailIdx int, adj [][]int, pos []Vec3) Vec3 {
	if len(adj[tailIdx]) == 0 { return Vec3{1, 0, 0} }
	sum := Vec3{}
	for _, nb := range adj[tailIdx] { sum = sum.Add(pos[nb].Sub(pos[tailIdx]).unit()) }
	if sum.Norm() < 1e-6 { return perpendicular(pos[adj[tailIdx][0]].Sub(pos[tailIdx])) }
	return sum.Scale(-1).unit()
}
