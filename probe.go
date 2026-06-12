package main

// Bond connects two atoms (1-based indices) with a given bond order.
// Order: 1=single, 2=double, 4=aromatic.
type Bond struct{ I, J, Order int }

// ProbeSet accumulates the atoms and bonds of all placed probes.
// It is the single mutable output passed through the placement pipeline.
type ProbeSet struct {
	Pos       []Vec3
	Labels    []string
	Bonds     []Bond
	Name      string // optional label (used for ring scan variants)
	ParentIdx int    // index of the original linked pair this was grown from (-1 = original)
	C1, C2    int    // probe component indices this pair connects (-1 if unknown)
}

// Add appends a new atom and returns its 1-based index.
func (ps *ProbeSet) Add(pos Vec3, label string) int {
	ps.Pos = append(ps.Pos, pos)
	ps.Labels = append(ps.Labels, label)
	return len(ps.Pos)
}

// Bond records a bond between two previously added atoms (1-based).
func (ps *ProbeSet) Bond(i, j, order int) {
	ps.Bonds = append(ps.Bonds, Bond{i, j, order})
}

// Len returns the current atom count.
func (ps *ProbeSet) Len() int { return len(ps.Pos) }

// Placer is the interface every pharmacophoric probe source implements.
// Implementing a new probe type requires only adding a new Placer — the
// runner (main) never needs to change (Open/Closed Principle).
type Placer interface {
	// Name returns a short label used in progress output.
	Name() string
	// Place adds probes to ps and returns the count placed.
	Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int
}
