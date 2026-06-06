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
