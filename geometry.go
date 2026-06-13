package main

import "math"

// sphereCircle samples points on the intersection circle of two spheres of
// radius r centred at c1 and c2.  Returns nil if the spheres don't intersect.
func sphereCircle(c1, c2 Vec3, r float64, n int) []Vec3 {
	d := c1.Sub(c2).Norm()
	if d < 1e-6 || d > 2*r {
		return nil
	}
	a := d / 2
	h := math.Sqrt(r*r - a*a)
	mid := c1.Add(c2).Scale(0.5)
	axis := c2.Sub(c1).unit()
	p := perpendicular(axis)
	q := cross3(axis, p)
	pts := make([]Vec3, n)
	for i := range pts {
		θ := 2 * math.Pi * float64(i) / float64(n)
		pts[i] = mid.Add(p.Scale(h * math.Cos(θ))).Add(q.Scale(h * math.Sin(θ)))
	}
	return pts
}

// openValenceDirs returns candidate directions for a new bond off an atom at
// pos, given the positions of its existing bonded neighbours and whether it is
// sp2. The directions respect the existing bond geometry so a linker (or any new
// substituent) joins at a realistic angle instead of an arbitrary one:
//
//   0 neighbours → nil (caller falls back to a free axis)
//   1 neighbour  → a cone of candidates at 120° (sp2) or 109.5° (sp3)
//   2 neighbours → the exocyclic bisector (sp2) or the two open tetrahedral
//                  directions (sp3)
//   3 neighbours → the single remaining tetrahedral direction
//
// The caller picks whichever returned direction best aligns with its target.
func openValenceDirs(pos Vec3, nbPos []Vec3, isSp2 bool) []Vec3 {
	e := make([]Vec3, len(nbPos))
	for i, p := range nbPos {
		d := p.Sub(pos)
		if d.Norm() < 1e-6 {
			return nil
		}
		e[i] = d.unit()
	}
	switch len(e) {
	case 1:
		incoming := e[0]
		baseP := perpendicular(incoming)
		p2 := cross3(incoming, baseP)
		cosT := -0.333 // sp3 tetrahedral
		if isSp2 {
			cosT = -0.5 // sp2 trigonal (120°)
		}
		sinT := math.Sqrt(1 - cosT*cosT)
		dirs := make([]Vec3, 0, 12)
		for k := 0; k < 12; k++ {
			a := float64(k) * math.Pi / 6
			perp := baseP.Scale(math.Cos(a)).Add(p2.Scale(math.Sin(a)))
			dirs = append(dirs, incoming.Scale(cosT).Add(perp.Scale(sinT)).unit())
		}
		return dirs
	case 2:
		sum := e[0].Add(e[1])
		if isSp2 {
			if sum.Norm() < 1e-6 {
				return []Vec3{perpendicular(e[0])}
			}
			return []Vec3{sum.Scale(-1).unit()}
		}
		// sp3 with two bonds: two open directions at 109.5° to each existing bond,
		// in the plane perpendicular to the (e0,e1) plane through their bisector.
		var b Vec3
		if sum.Norm() < 1e-6 {
			b = perpendicular(e[0])
		} else {
			b = sum.unit()
		}
		nrm := cross3(e[0], e[1])
		if nrm.Norm() < 1e-6 {
			nrm = perpendicular(b)
		}
		nrm = nrm.unit()
		const c, s = 0.5774, 0.8165 // cosθ, sinθ for θ ≈ 54.7°
		return []Vec3{
			b.Scale(-c).Add(nrm.Scale(s)).unit(),
			b.Scale(-c).Sub(nrm.Scale(s)).unit(),
		}
	case 3:
		sum := e[0].Add(e[1]).Add(e[2])
		if sum.Norm() < 1e-6 {
			return nil
		}
		return []Vec3{sum.Scale(-1).unit()}
	}
	return nil
}

// benzeneRing returns 6 ring atom positions and their aromatic bond pairs for
// a probe ring centred at centre with the given normal direction.
func benzeneRing(centre, normal Vec3) ([]Vec3, [][2]int) {
	u := perpendicular(normal)
	v := cross3(normal, u)
	const r = 1.40
	pts := make([]Vec3, 6)
	for i := range pts {
		θ := float64(i) * math.Pi / 3
		pts[i] = centre.Add(u.Scale(r * math.Cos(θ))).Add(v.Scale(r * math.Sin(θ)))
	}
	edges := [][2]int{{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}}
	return pts, edges
}
