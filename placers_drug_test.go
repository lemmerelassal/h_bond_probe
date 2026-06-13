package main

import (
	"math"
	"testing"
)

// sanitizeKeepsAll builds a ProbeSet from the given atoms/bonds, runs it through
// sanitizeValence, and reports whether every bond survived (i.e. no atom was
// over-valent) — the invariant that keeps rings and S(VI) groups from breaking.
func sanitizeKeepsAll(t *testing.T, name string, pos []Vec3, lab []string, bonds []Bond) {
	t.Helper()
	mol := ProbeSet{Pos: pos, Labels: lab, Bonds: bonds, ParentIdx: -1}
	out := sanitizeValence(mol)
	if len(out.Bonds) != len(bonds) {
		t.Errorf("%s: sanitizeValence stripped %d/%d bonds (over-valent atom)",
			name, len(bonds)-len(out.Bonds), len(bonds))
	}
}

// TestHeteroRingValences verifies that every heterocyclic ring probe is emitted
// within valence, so sanitizeValence never strips a ring bond and breaks it.
func TestHeteroRingValences(t *testing.T) {
	for _, spec := range heteroRingSpecs {
		var pos []Vec3
		var lab []string
		var bonds []Bond
		ok := emitHeteroRing(spec, Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0},
			&pos, &lab, &bonds, nil) // nil heavy → no clashes
		if !ok {
			t.Fatalf("%s: emitHeteroRing failed in empty space", spec.name)
		}
		sanitizeKeepsAll(t, spec.name, pos, lab, bonds)
	}
}

// TestSanitizeProtectsRingBonds reproduces the downstream failure mode: an extra
// (acyclic) bond from linking/fusion over-fills a ring atom. sanitizeValence must
// drop the acyclic bond, never a ring bond — so the ring stays closed.
func TestSanitizeProtectsRingBonds(t *testing.T) {
	for _, spec := range heteroRingSpecs {
		var pos []Vec3
		var lab []string
		var bonds []Bond
		if !emitHeteroRing(spec, Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0},
			&pos, &lab, &bonds, nil) {
			t.Fatalf("%s: emitHeteroRing failed", spec.name)
		}
		nRing := len(spec.elements) // first n bonds are the ring edges
		// Attach two extra carbons to the directional ring atom (1-based dirIdx+1),
		// forcing it over its valence maximum.
		dir := spec.dirIdx + 1
		c1 := len(pos) + 1
		pos = append(pos, Vec3{0.5, 0.5, 0.5})
		lab = append(lab, "C")
		bonds = append(bonds, Bond{dir, c1, 1})
		c2 := len(pos) + 1
		pos = append(pos, Vec3{-0.5, 0.5, 0.5})
		lab = append(lab, "C")
		bonds = append(bonds, Bond{dir, c2, 1})

		out := sanitizeValence(ProbeSet{Pos: pos, Labels: lab, Bonds: bonds, ParentIdx: -1})

		kept := map[[2]int]bool{}
		for _, b := range out.Bonds {
			kept[[2]int{b.I, b.J}] = true
		}
		for k := 0; k < nRing; k++ {
			rb := bonds[k]
			if !kept[[2]int{rb.I, rb.J}] {
				t.Errorf("%s: ring bond %d-%d was stripped (ring broken)", spec.name, rb.I, rb.J)
			}
		}
	}
}

// TestEsterOxygenBent verifies the ester/carbamate ester oxygen is bent (not the
// 180° linear oxygen that an earlier collinear construction produced).
func TestEsterOxygenBent(t *testing.T) {
	var pos []Vec3
	var lab []string
	var bonds []Bond
	n := placeEsters(
		[]Atom{
			{Name: "NZ", ResName: "LYS", Element: "N", Pos: Vec3{0, 0, 0}},
			{Name: "NH1", ResName: "ARG", Element: "N", Pos: Vec3{3.0, 0, 0}},
			{Name: "ND2", ResName: "ASN", Element: "N", Pos: Vec3{1.5, 2.6, 0}},
		},
		&pos, &lab, &bonds)
	if n == 0 {
		t.Skip("no ester placed for synthetic donors")
	}
	adj := buildAdj(len(pos), bonds)
	for i, el := range lab {
		if el != "O" || len(adj[i]) != 2 {
			continue
		}
		u := pos[adj[i][0]].Sub(pos[i])
		v := pos[adj[i][1]].Sub(pos[i])
		ang := 180 / math.Pi * math.Acos(u.unit().Dot(v.unit()))
		if ang > 150 {
			t.Errorf("ester oxygen %d is near-linear (%.1f°); must be bent", i+1, ang)
		}
	}
}

// TestSulfonamideValence verifies the S(VI) sulfonamide keeps both S=O bonds
// (max valence for sulfur must be 6, not 2).
func TestSulfonamideValence(t *testing.T) {
	var pos []Vec3
	var lab []string
	var bonds []Bond
	n := placeSulfonamides(
		[]Atom{
			{Name: "OD1", ResName: "ASP", Element: "O", Pos: Vec3{0, 0, 0}},
			{Name: "OD2", ResName: "ASP", Element: "O", Pos: Vec3{2.2, 0, 0}},
			{Name: "OE1", ResName: "GLU", Element: "O", Pos: Vec3{1.1, 2.0, 0}},
		},
		&pos, &lab, &bonds)
	if n == 0 {
		t.Skip("no sulfonamide placed for synthetic acceptors")
	}
	sanitizeKeepsAll(t, "sulfonamide", pos, lab, bonds)
}
