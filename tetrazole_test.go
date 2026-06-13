package main

import "testing"

// buildAdj builds a 0-based adjacency list from 1-based bonds over n atoms.
func buildAdj(n int, bonds []Bond) [][]int {
	adj := make([][]int, n)
	for _, b := range bonds {
		adj[b.I-1] = append(adj[b.I-1], b.J-1)
		adj[b.J-1] = append(adj[b.J-1], b.I-1)
	}
	return adj
}

// TestTetrazoleAttachmentEligibility verifies the invariant behind the
// broken-ring fix: a linker may only attach to an atom with spare valence, and
// every double-bonded ring nitrogen of a tetrazole has none — so a linker can
// never over-fill one and force sanitizeValence to strip a ring bond. The stem
// carbon and the singly-bonded ring N remain eligible.
func TestTetrazoleAttachmentEligibility(t *testing.T) {
	var pos []Vec3
	var lab []string
	var bonds []Bond
	ring := tetrazoleRing(Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0})
	appendTetrazole(&pos, &lab, &bonds, ring)
	// Index map (0-based): 0=C5, 1=N1, 2=N2, 3=N3, 4=N4, 5=stem C.

	val := atomValences(len(pos), bonds)
	core := ringCoreAtoms(buildAdj(len(pos), bonds))

	eligible := func(i int) bool { return val[i]+1 <= maxValence(lab[i]) }

	// Ring atoms 0..4 must be flagged as ring (2-core); stem (5) must not.
	for i := 0; i <= 4; i++ {
		if !core[i] {
			t.Errorf("atom %d (%s) should be a ring atom", i, lab[i])
		}
	}
	if core[5] {
		t.Error("stem carbon (5) should not be a ring atom")
	}

	// Stem carbon is the intended attachment handle: must be eligible.
	if !eligible(5) {
		t.Errorf("stem carbon must be eligible for linking (valence %d)", val[5])
	}

	// The double-bonded ring nitrogens (N2,N3,N4 = indices 2,3,4) are full and
	// must be ineligible — attaching there is exactly what broke the ring.
	for _, i := range []int{2, 3, 4} {
		if eligible(i) {
			t.Errorf("ring nitrogen %d is full (valence %d) and must NOT be linkable", i, val[i])
		}
	}
	// C5 itself is full (two ring bonds + stem) and must be ineligible too.
	if eligible(0) {
		t.Errorf("C5 is full (valence %d) and must NOT be linkable", val[0])
	}
}
