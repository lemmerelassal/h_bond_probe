package main

import "strings"

// ── Shared placer helper ──────────────────────────────────────────────────────

// clashFn returns a closure that checks whether a position clashes with any
// protein atom, excluding the source residue identified by (chain, seq).
func clashFn(heavy []heavyAtom) func(pos Vec3, elem string) bool {
	return func(pos Vec3, elem string) bool {
		r := vdw(elem)
		for _, h := range heavy {
			if pos.Sub(h.pos).Norm() < h.vdwR+r-hardTol {
				return false
			}
		}
		return true
	}
}

// clashFnExcluding is like clashFn but skips atoms from a specific residue.
func clashFnExcluding(atoms []Atom, chain string, seq int) func(pos Vec3, elem string) bool {
	var excl []heavyAtom
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) || a.IsHet {
			continue
		}
		if a.ChainID == chain && a.ResSeq == seq {
			continue
		}
		excl = append(excl, heavyAtom{a.Pos, vdw(a.Element)})
	}
	return clashFn(excl)
}


// ── Placer adapters ───────────────────────────────────────────────────────────
// Each concrete Placer wraps one internal placement function.
// Adding a new residue type only requires adding a new adapter; main never changes.

type argAcidPlacer struct{}
func (argAcidPlacer) Name() string { return "ARG acetic acid" }
func (argAcidPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeArgAcids(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type carboxylateGuanidinePlacer struct{}
func (carboxylateGuanidinePlacer) Name() string { return "ASP/GLU guanidine" }
func (carboxylateGuanidinePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeCarboxylateGuanidines(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type tyrAcetamidePlacer struct{}
func (tyrAcetamidePlacer) Name() string { return "TYR acetamide" }
func (tyrAcetamidePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeTyrAcetamides(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type tyrPhenolPlacer struct{}
func (tyrPhenolPlacer) Name() string { return "TYR phenol" }
func (tyrPhenolPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeTyrPhenols(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type asnGlnAcetamidePlacer struct{}
func (asnGlnAcetamidePlacer) Name() string { return "ASN/GLN acetamide" }
func (asnGlnAcetamidePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeAsnGlnAcetamides(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type serThrMethanolPlacer struct{}
func (serThrMethanolPlacer) Name() string { return "SER/THR methanol" }
func (serThrMethanolPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeSerThrMethanols(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type lysAcetatePlacer struct{}
func (lysAcetatePlacer) Name() string { return "LYS acetate" }
func (lysAcetatePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeLysAcetates(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type cysThiolPlacer struct{}
func (cysThiolPlacer) Name() string { return "CYS methanethiol" }
func (cysThiolPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeCysMethanethiols(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type hisImidazolePlacer struct{}
func (hisImidazolePlacer) Name() string { return "HIS imidazole" }
func (hisImidazolePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeHisImidazoles(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type waterMethanolPlacer struct{}
func (waterMethanolPlacer) Name() string { return "water displacement" }
func (waterMethanolPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeWaterMethanols(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

// AllPlacers is the canonical ordered list of pharmacophoric placers.
// To add a new probe type: implement Placer, append here.  main never changes.
var AllPlacers = []Placer{
	argAcidPlacer{},
	carboxylateGuanidinePlacer{},
	tyrAcetamidePlacer{},
	tyrPhenolPlacer{},
	asnGlnAcetamidePlacer{},
	serThrMethanolPlacer{},
	lysAcetatePlacer{},
	cysThiolPlacer{},
	hisImidazolePlacer{},
	waterMethanolPlacer{},
}
