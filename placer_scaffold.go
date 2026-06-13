package main

import (
	"sort"
	"strings"
)

// residueKey identifies a residue by chain ID and sequence number.
// (The package already has a resKey *function* in classify.go, so the residue
// grouping key used by the placers gets the longer, unambiguous name.)
type residueKey struct {
	chain string
	seq   int
}

// resClashSet holds every protein heavy atom tagged with the residue it belongs
// to. It answers clash queries that exclude a given residue's own atoms — the
// shared need of every residue-anchored placer, which must let a probe overlap
// the space of the residue it is built against while clearing all others.
type resClashSet struct {
	pos  []Vec3
	vdwR []float64
	key  []residueKey
}

// newResClashSet builds the clash set from all non-hydrogen, non-water atoms.
func newResClashSet(atoms []Atom) *resClashSet {
	s := &resClashSet{}
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		s.pos = append(s.pos, a.Pos)
		s.vdwR = append(s.vdwR, vdw(a.Element))
		s.key = append(s.key, residueKey{a.ChainID, a.ResSeq})
	}
	return s
}

// clears reports whether an atom of element elem placed at pos avoids a hard
// clash with every heavy atom that does not belong to residue exclude.
func (s *resClashSet) clears(pos Vec3, elem string, exclude residueKey) bool {
	r := vdw(elem)
	for i := range s.pos {
		if s.key[i] == exclude {
			continue
		}
		if pos.Sub(s.pos[i]).Norm() < s.vdwR[i]+r-hardTol {
			return false
		}
	}
	return true
}

// maxValence is the bonding budget per element used to decide whether an atom
// can accept one more bond (e.g. a linker) without becoming over-valent.
func maxValence(elem string) int {
	switch strings.ToUpper(elem) {
	case "C":
		return 4
	case "N":
		return 3
	case "O":
		return 2
	case "S":
		// Sulfur is hypervalent in sulfonamides/sulfones (S(VI), valence 6);
		// thioether/thiol probes are terminal so the looser bound is harmless.
		return 6
	case "F", "CL", "BR":
		return 1
	}
	return 4
}

// atomValences returns the summed bond order incident on each of n atoms
// (1-based bonds); aromatic bonds (order 4) count as 2, matching sanitizeValence.
func atomValences(n int, bonds []Bond) []int {
	v := make([]int, n)
	for _, b := range bonds {
		o := b.Order
		if o == 4 {
			o = 2
		}
		v[b.I-1] += o
		v[b.J-1] += o
	}
	return v
}

// ringCoreAtoms marks the 2-core of the bond graph — the atoms that survive
// iterative removal of degree-1 nodes, i.e. the atoms lying on a ring (plus any
// atoms on paths bridging two rings). Used to steer linkers toward non-ring
// "handle" atoms.
func ringCoreAtoms(adj [][]int) []bool {
	n := len(adj)
	inCore := make([]bool, n)
	deg := make([]int, n)
	var q []int
	for i := 0; i < n; i++ {
		deg[i] = len(adj[i])
		inCore[i] = true
	}
	for i := 0; i < n; i++ {
		if deg[i] <= 1 {
			q = append(q, i)
		}
	}
	for len(q) > 0 {
		u := q[0]
		q = q[1:]
		if !inCore[u] {
			continue
		}
		inCore[u] = false
		for _, w := range adj[u] {
			if inCore[w] {
				deg[w]--
				if deg[w] <= 1 {
					q = append(q, w)
				}
			}
		}
	}
	return inCore
}

// sortedResidueKeys returns a residue-group map's keys in deterministic
// (chain, then seq) order. Placers iterate with this instead of ranging the map
// directly so that placement and output ordering are reproducible run-to-run
// (Go randomizes map iteration order).
func sortedResidueKeys[V any](m map[residueKey]V) []residueKey {
	keys := make([]residueKey, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool {
		if keys[i].chain != keys[j].chain {
			return keys[i].chain < keys[j].chain
		}
		return keys[i].seq < keys[j].seq
	})
	return keys
}
