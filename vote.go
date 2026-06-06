package main

import "math"

const (
	sphereR      = 3.5 // Å — hydrophobic probe sphere radius
	nSamples     = 16  // samples per sphere-intersection circle
	maxPairDist  = 7.0 // Å — max distance between hydrophobic atom pairs
	minPairs     = 2   // minimum pair count to accept a vote cell
	probeSpacing = 1.5 // Å — minimum distance between placed probes

	hbondSphereR    = 2.9 // Å — H-bond donor/acceptor probe sphere radius
	hbondPairDist   = 6.0 // Å — max distance between H-bond atom pairs
	hbondMinPairs   = 2
	hbondProbeSpacing = 1.2 // Å
)

type voteCell struct {
	pos       Vec3
	pairCount int
	nPts      int
	pairs     map[[2]int]bool
}

func gridKey(p Vec3) [3]int {
	return [3]int{
		int(math.Floor(p.X / probeSpacing)),
		int(math.Floor(p.Y / probeSpacing)),
		int(math.Floor(p.Z / probeSpacing)),
	}
}

// castVotes runs the sphere-intersection vote system over a group of atoms.
// Each point that clears the protein clash check casts a vote for its grid cell.
func castVotes(group []Atom, r, pairDist float64, n int, heavy []heavyAtom, probeElem string) map[[3]int]*voteCell {
	probeR := vdw(probeElem)
	v := map[[3]int]*voteCell{}
	for i := 0; i < len(group); i++ {
		for j := i + 1; j < len(group); j++ {
			if group[i].Pos.Sub(group[j].Pos).Norm() > pairDist {
				continue
			}
			pts := sphereCircle(group[i].Pos, group[j].Pos, r, n)
			pk := [2]int{i, j}
			for _, p := range pts {
				if !noHardClash(p, probeR, hardTol, heavy) {
					continue
				}
				key := gridKey(p)
				if e, ok := v[key]; ok {
					if !e.pairs[pk] {
						e.pairs[pk] = true
						e.pairCount++
					}
					e.nPts++
					nn := float64(e.nPts)
					e.pos = Vec3{
						(e.pos.X*(nn-1) + p.X) / nn,
						(e.pos.Y*(nn-1) + p.Y) / nn,
						(e.pos.Z*(nn-1) + p.Z) / nn,
					}
				} else {
					v[key] = &voteCell{
						pairCount: 1,
						pos:       p,
						nPts:      1,
						pairs:     map[[2]int]bool{pk: true},
					}
				}
			}
		}
	}
	return v
}
