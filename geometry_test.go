package main

import (
	"math"
	"testing"
)

// angBetween returns the angle (degrees) between vectors a and b.
func angBetween(a, b Vec3) float64 {
	return 180 / math.Pi * math.Acos(a.unit().Dot(b.unit()))
}

// TestLinkerChainTetrahedral verifies buildLinkerChain produces a chain whose
// every junction — including the two end junctions at the probe atoms, which the
// old zigzag left strained — is tetrahedral (no acute angles), across a range of
// gaps and exit-tangent orientations.
func TestLinkerChainTetrahedral(t *testing.T) {
	grid := newClashGrid([]heavyAtom{}) // no protein → never clashes
	cases := []struct {
		p1, d1, p2, d2 Vec3
	}{
		{Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{6, 0, 0}, Vec3{-1, 0, 0}},
		{Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{5, 3, 0}, Vec3{0, -1, 0}},
		{Vec3{0, 0, 0}, Vec3{0, 1, 0}, Vec3{4, 4, 2}, Vec3{-1, 0, 0}},
		{Vec3{0, 0, 0}, Vec3{1, 1, 0}, Vec3{8, 1, 1}, Vec3{-1, 0, 0}},
	}
	for ci, c := range cases {
		chain := buildLinkerChain(c.p1, c.d1.unit(), c.p2, c.d2.unit(), grid, 9)
		if chain == nil {
			t.Errorf("case %d: no chain found", ci)
			continue
		}
		seq := append([]Vec3{c.p1}, chain...)
		seq = append(seq, c.p2)
		for i := 1; i < len(seq)-1; i++ {
			if a := angBetween(seq[i-1].Sub(seq[i]), seq[i+1].Sub(seq[i])); a < 95 {
				t.Errorf("case %d: junction %d angle %.1f° < 95°", ci, i, a)
			}
		}
		for i := 0; i < len(seq)-1; i++ {
			if d := seq[i].Sub(seq[i+1]).Norm(); d < 1.2 || d > 1.8 {
				t.Errorf("case %d: bond %d length %.2f Å out of range", ci, i, d)
			}
		}
	}
}

// TestOpenValenceDirsAngles verifies the linker-attachment helper returns bond
// directions at the correct angle to every existing bond — the fix that stopped
// the linker from attaching at arbitrary (30–70°) angles to atoms with 2 sp3 or
// 3 neighbours.
func TestOpenValenceDirsAngles(t *testing.T) {
	o := Vec3{0, 0, 0}

	// sp3, two existing bonds → both open dirs at ~109.5° to each.
	e0 := Vec3{1, 0, 0}
	e1 := Vec3{-0.33, 0.94, 0}
	for _, d := range openValenceDirs(o, []Vec3{e0, e1}, false) {
		for _, e := range []Vec3{e0, e1} {
			if a := angBetween(d, e); math.Abs(a-109.5) > 3 {
				t.Errorf("sp3/2nbr: open dir %v makes %.1f° with %v (want ~109.5°)", d, a, e)
			}
		}
	}

	// sp3, three existing bonds → single open dir at ~109.5° to each.
	n0 := Vec3{1, 1, 1}.unit()
	n1 := Vec3{1, -1, -1}.unit()
	n2 := Vec3{-1, 1, -1}.unit()
	dirs := openValenceDirs(o, []Vec3{n0, n1, n2}, false)
	if len(dirs) != 1 {
		t.Fatalf("sp3/3nbr: expected 1 open dir, got %d", len(dirs))
	}
	for _, e := range []Vec3{n0, n1, n2} {
		if a := angBetween(dirs[0], e); math.Abs(a-109.5) > 5 {
			t.Errorf("sp3/3nbr: open dir makes %.1f° with %v (want ~109.5°)", a, e)
		}
	}

	// sp2, two existing bonds → exocyclic bisector at ~120° to each.
	s0 := Vec3{1, 0, 0}
	s1 := Vec3{-0.5, 0.866, 0}
	ex := openValenceDirs(o, []Vec3{s0, s1}, true)
	if len(ex) != 1 {
		t.Fatalf("sp2/2nbr: expected 1 open dir, got %d", len(ex))
	}
	for _, e := range []Vec3{s0, s1} {
		if a := angBetween(ex[0], e); math.Abs(a-120) > 3 {
			t.Errorf("sp2/2nbr: open dir makes %.1f° with %v (want ~120°)", a, e)
		}
	}
}
