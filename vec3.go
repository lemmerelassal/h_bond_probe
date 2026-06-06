package main

import "math"

// Vec3 is an immutable 3-D vector.  All methods return new values.
type Vec3 struct{ X, Y, Z float64 }

func (a Vec3) Sub(b Vec3) Vec3      { return Vec3{a.X - b.X, a.Y - b.Y, a.Z - b.Z} }
func (a Vec3) Add(b Vec3) Vec3      { return Vec3{a.X + b.X, a.Y + b.Y, a.Z + b.Z} }
func (a Vec3) Scale(s float64) Vec3 { return Vec3{a.X * s, a.Y * s, a.Z * s} }
func (a Vec3) Dot(b Vec3) float64   { return a.X*b.X + a.Y*b.Y + a.Z*b.Z }
func (a Vec3) Norm() float64        { return math.Sqrt(a.Dot(a)) }
func (a Vec3) unit() Vec3           { return a.Scale(1.0 / a.Norm()) }

func cross3(a, b Vec3) Vec3 {
	return Vec3{
		a.Y*b.Z - a.Z*b.Y,
		a.Z*b.X - a.X*b.Z,
		a.X*b.Y - a.Y*b.X,
	}
}

// perpendicular returns an arbitrary unit vector perpendicular to v.
func perpendicular(v Vec3) Vec3 {
	if math.Abs(v.X) <= math.Abs(v.Y) && math.Abs(v.X) <= math.Abs(v.Z) {
		return cross3(v, Vec3{1, 0, 0}).unit()
	}
	if math.Abs(v.Y) <= math.Abs(v.Z) {
		return cross3(v, Vec3{0, 1, 0}).unit()
	}
	return cross3(v, Vec3{0, 0, 1}).unit()
}

// ringPlaneNormal returns the least-squares normal of a set of ring atom positions.
func ringPlaneNormal(pts []Vec3) Vec3 {
	n := Vec3{}
	for i := 0; i < len(pts); i++ {
		j := (i + 1) % len(pts)
		n = n.Add(cross3(pts[i], pts[j]))
	}
	if n.Norm() < 1e-9 {
		return Vec3{0, 0, 1}
	}
	return n.unit()
}
