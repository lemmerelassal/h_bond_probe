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

type argTetrazolePlacer struct{}
func (argTetrazolePlacer) Name() string { return "ARG tetrazole" }
func (argTetrazolePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeArgTetrazoles(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type lysTetrazolePlacer struct{}
func (lysTetrazolePlacer) Name() string { return "LYS tetrazole" }
func (lysTetrazolePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeLysTetrazoles(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
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

type trpCarbonylPlacer struct{}
func (trpCarbonylPlacer) Name() string { return "TRP indole-NH carbonyl" }
func (trpCarbonylPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeTrpCarbonyls(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type hisCarbonylPlacer struct{}
func (hisCarbonylPlacer) Name() string { return "HIS carbonyl partner" }
func (hisCarbonylPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeHisCarbonyls(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type hisAcidPlacer struct{}
func (hisAcidPlacer) Name() string { return "HIS acid (protonated)" }
func (hisAcidPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeHisAcids(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type carboxylateAminePlacer struct{}
func (carboxylateAminePlacer) Name() string { return "ASP/GLU amine" }
func (carboxylateAminePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeCarboxylateAmines(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type metArenePlacer struct{}
func (metArenePlacer) Name() string { return "MET S-arene" }
func (metArenePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeMetArenes(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type aromaticCationPlacer struct{}
func (aromaticCationPlacer) Name() string { return "aromatic cation-pi" }
func (aromaticCationPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeAromaticCations(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type cationArenePlacer struct{}
func (cationArenePlacer) Name() string { return "ARG/LYS cation-pi arene" }
func (cationArenePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeCationArenes(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type hydroxylVotePlacer struct{}
func (hydroxylVotePlacer) Name() string { return "hydroxyl/methoxy hotspots" }
func (hydroxylVotePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeHydroxylVotes(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type backboneHydroxylPlacer struct{}
func (backboneHydroxylPlacer) Name() string { return "backbone C=O hydroxyl" }
func (backboneHydroxylPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeBackboneHydroxyls(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type heteroRingPlacer struct{}
func (heteroRingPlacer) Name() string { return "N-heterocycle rings" }
func (heteroRingPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeHeteroRings(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type halogenBondPlacer struct{}
func (halogenBondPlacer) Name() string { return "halogen bond (Cl/Br/F)" }
func (halogenBondPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeHalogenBonds(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type trifluoromethylPlacer struct{}
func (trifluoromethylPlacer) Name() string { return "trifluoromethyl" }
func (trifluoromethylPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeTrifluoromethyls(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type sulfonamidePlacer struct{}
func (sulfonamidePlacer) Name() string { return "sulfonamide" }
func (sulfonamidePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeSulfonamides(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type nitrilePlacer struct{}
func (nitrilePlacer) Name() string { return "nitrile" }
func (nitrilePlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeNitriles(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

type esterPlacer struct{}
func (esterPlacer) Name() string { return "ester/carbamate" }
func (esterPlacer) Place(atoms []Atom, heavy []heavyAtom, ps *ProbeSet) int {
	return placeEsters(atoms, &ps.Pos, &ps.Labels, &ps.Bonds)
}

// AllPlacers is the canonical ordered list of pharmacophoric placers.
// To add a new probe type: implement Placer, append here.  main never changes.
var AllPlacers = []Placer{
	argAcidPlacer{},
	argTetrazolePlacer{},
	carboxylateGuanidinePlacer{},
	tyrAcetamidePlacer{},
	tyrPhenolPlacer{},
	asnGlnAcetamidePlacer{},
	serThrMethanolPlacer{},
	lysAcetatePlacer{},
	lysTetrazolePlacer{},
	cysThiolPlacer{},
	hisImidazolePlacer{},
	waterMethanolPlacer{},
	trpCarbonylPlacer{},
	hisCarbonylPlacer{},
	hisAcidPlacer{},
	carboxylateAminePlacer{},
	metArenePlacer{},
	aromaticCationPlacer{},
	cationArenePlacer{},
	hydroxylVotePlacer{},
	backboneHydroxylPlacer{},
	heteroRingPlacer{},
	halogenBondPlacer{},
	trifluoromethylPlacer{},
	sulfonamidePlacer{},
	nitrilePlacer{},
	esterPlacer{},
}
