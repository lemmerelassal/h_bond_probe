package main

import (
	"strings"
)

// Atom holds the parsed data for one PDB/mmCIF heavy atom.
type Atom struct {
	Name    string
	ResName string
	ChainID string
	ResSeq  int
	Pos     Vec3
	Element string
	IsHet   bool // true for HETATM records (ligands, cofactors, ions)
}

// heavyAtom is the minimal representation used for clash checking.
type heavyAtom struct {
	pos  Vec3
	vdwR float64
}

var vdwRadius = map[string]float64{
	"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
	"F": 1.47, "CL": 1.75, "BR": 1.85, "P": 1.80,
	"SE": 1.90,
}

func vdw(elem string) float64 {
	if r, ok := vdwRadius[strings.ToUpper(elem)]; ok {
		return r
	}
	return 1.70
}
