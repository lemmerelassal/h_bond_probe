package main

import "strings"

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
