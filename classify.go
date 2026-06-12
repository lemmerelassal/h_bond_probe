package main

import (
	"math"
	"strings"
)

const hardTol = 0.3 // Å — vdW overlap tolerance for clash detection

func isBackbone(name string) bool {
	switch strings.TrimSpace(strings.ToUpper(name)) {
	case "N", "CA", "C", "O", "OXT":
		return true
	}
	return false
}

// isWater returns true for common water residue names.
func isWater(a Atom) bool {
	switch strings.TrimSpace(strings.ToUpper(a.ResName)) {
	case "HOH", "WAT", "H2O", "DOD", "TIP", "SOL":
		return true
	}
	return false
}

func isHydrophobic(a Atom) bool {
	if strings.ToUpper(a.Element) != "C" {
		return false
	}
	if isBackbone(a.Name) {
		return false
	}
	name := strings.TrimSpace(strings.ToUpper(a.Name))
	res := strings.TrimSpace(strings.ToUpper(a.ResName))
	switch res {
	case "ALA":
		return name == "CB"
	case "VAL":
		switch name {
		case "CB", "CG1", "CG2":
			return true
		}
	case "LEU":
		switch name {
		case "CB", "CG", "CD1", "CD2":
			return true
		}
	case "ILE":
		switch name {
		case "CB", "CG1", "CG2", "CD1":
			return true
		}
	case "PRO":
		switch name {
		case "CB", "CG", "CD":
			return true
		}
	case "MET":
		switch name {
		case "CB", "CG", "CE":
			return true
		}
	case "MSE":
		switch name {
		case "CB", "CG", "CE":
			return true
		}
	case "PHE":
		switch name {
		case "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ":
			return true
		}
	case "TYR":
		switch name {
		case "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ":
			return true
		}
	case "TRP":
		switch name {
		case "CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2":
			return true
		}
	case "SER":
		return name == "CB"
	case "THR":
		switch name {
		case "CB", "CG2":
			return true
		}
	case "LYS":
		switch name {
		case "CB", "CG", "CD", "CE":
			return true
		}
	case "CYS":
		return name == "CB"
	}
	return false
}

// isAromatic returns true for the ring atoms of PHE, TYR, TRP, and HIS.
func isAromatic(a Atom) bool {
	name := strings.TrimSpace(strings.ToUpper(a.Name))
	res := strings.TrimSpace(strings.ToUpper(a.ResName))
	elem := strings.ToUpper(a.Element)
	switch res {
	case "PHE", "TYR":
		if elem != "C" {
			return false
		}
		switch name {
		case "CG", "CD1", "CD2", "CE1", "CE2", "CZ":
			return true
		}
	case "TRP":
		if elem != "C" {
			return false
		}
		switch name {
		case "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2":
			return true
		}
	case "HIS":
		if elem != "C" && elem != "N" {
			return false
		}
		switch name {
		case "CG", "CD2", "CE1", "ND1", "NE2":
			return true
		}
	}
	return false
}

func resKey(a Atom) string {
	return a.ChainID + ":" + strings.TrimSpace(strings.ToUpper(a.ResName)) + ":" + itoa(a.ResSeq)
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	b := make([]byte, 0, 10)
	for n > 0 {
		b = append([]byte{byte('0' + n%10)}, b...)
		n /= 10
	}
	if neg {
		b = append([]byte{'-'}, b...)
	}
	return string(b)
}

// hbondRole classifies an atom's hydrogen-bond capability.
// Returns "donor", "acceptor", "dual", or "" (no role / skip).
func hbondRole(a Atom) string {
	if a.IsHet || isWater(a) || isBackbone(a.Name) {
		return ""
	}
	elem := strings.ToUpper(a.Element)
	if elem == "H" {
		return ""
	}
	name := strings.TrimSpace(strings.ToUpper(a.Name))
	switch elem {
	case "N":
		return "dual"
	case "O":
		if name == "OXT" {
			return "acceptor"
		}
		return "dual"
	case "S":
		if name == "SG" {
			return "dual"
		}
		return "acceptor"
	case "SE":
		if name == "SE" || name == "SEG" {
			return "dual"
		}
		return "acceptor"
	}
	return ""
}

// noHardClash returns true when pos does not overlap any atom in heavy.
func noHardClash(pos Vec3, probeR, tol float64, heavy []heavyAtom) bool {
	for _, h := range heavy {
		if pos.Sub(h.pos).Norm() < h.vdwR+probeR-tol {
			return false
		}
	}
	return true
}

// clashGrid is a spatial hash that makes clash checks O(1) instead of
// O(protein_atoms). Cell size must be ≥ max expected clash distance (~3.5 Å).
type clashGrid struct {
	cells    map[[3]int][]int
	atoms    []heavyAtom
	cellSize float64
	origin   Vec3
}

func newClashGrid(atoms []heavyAtom) *clashGrid {
	const cellSize = 4.0 // Å — larger than max vdW sum (~3.8 Å), so ±1 cell suffices
	g := &clashGrid{cells: make(map[[3]int][]int, len(atoms)), atoms: atoms, cellSize: cellSize}
	if len(atoms) == 0 {
		return g
	}
	g.origin = atoms[0].pos
	for _, a := range atoms {
		if a.pos.X < g.origin.X { g.origin.X = a.pos.X }
		if a.pos.Y < g.origin.Y { g.origin.Y = a.pos.Y }
		if a.pos.Z < g.origin.Z { g.origin.Z = a.pos.Z }
	}
	// Pad by one cell so indices are always non-negative.
	g.origin.X -= cellSize
	g.origin.Y -= cellSize
	g.origin.Z -= cellSize
	for i, a := range atoms {
		k := g.gridKey(a.pos)
		g.cells[k] = append(g.cells[k], i)
	}
	return g
}

func (g *clashGrid) gridKey(pos Vec3) [3]int {
	return [3]int{
		int(math.Floor((pos.X - g.origin.X) / g.cellSize)),
		int(math.Floor((pos.Y - g.origin.Y) / g.cellSize)),
		int(math.Floor((pos.Z - g.origin.Z) / g.cellSize)),
	}
}

// clashFree returns true when pos does not overlap any protein atom.
func (g *clashGrid) clashFree(pos Vec3, probeR, tol float64) bool {
	k := g.gridKey(pos)
	for dx := -1; dx <= 1; dx++ {
		for dy := -1; dy <= 1; dy++ {
			for dz := -1; dz <= 1; dz++ {
				nk := [3]int{k[0] + dx, k[1] + dy, k[2] + dz}
				for _, i := range g.cells[nk] {
					a := g.atoms[i]
					if pos.Sub(a.pos).Norm() < a.vdwR+probeR-tol {
						return false
					}
				}
			}
		}
	}
	return true
}

// countNearby returns the number of protein atoms at distance [minD, maxD] from pos.
// Uses ±2-cell neighbourhood (covers up to 8 Å with default 4 Å cells).
func (g *clashGrid) countNearby(pos Vec3, minD, maxD float64) int {
	k := g.gridKey(pos)
	n := 0
	for dx := -2; dx <= 2; dx++ {
		for dy := -2; dy <= 2; dy++ {
			for dz := -2; dz <= 2; dz++ {
				nk := [3]int{k[0] + dx, k[1] + dy, k[2] + dz}
				for _, i := range g.cells[nk] {
					d := pos.Sub(g.atoms[i].pos).Norm()
					if d >= minD && d <= maxD {
						n++
					}
				}
			}
		}
	}
	return n
}
