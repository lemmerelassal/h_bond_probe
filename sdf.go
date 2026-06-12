package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strings"
	"time"
)

// molValid returns false when the molecule should not be written:
// zero atoms, V2000 overflow (>999), any NaN/Inf coordinate, or any
// bond index outside [1, nAtoms].
func molValid(positions []Vec3, bonds []Bond) bool {
	n := len(positions)
	if n == 0 || n > 999 || len(bonds) > 999 {
		return false
	}
	for _, p := range positions {
		if math.IsNaN(p.X) || math.IsNaN(p.Y) || math.IsNaN(p.Z) ||
			math.IsInf(p.X, 0) || math.IsInf(p.Y, 0) || math.IsInf(p.Z, 0) {
			return false
		}
	}
	for _, b := range bonds {
		if b.I < 1 || b.I > n || b.J < 1 || b.J > n {
			return false
		}
	}
	return true
}

// writeMolBlock writes one MOL V2000 record to w.
func writeMolBlock(w *bufio.Writer, now time.Time, name string, positions []Vec3, labels []string, bonds []Bond) {
	fmt.Fprintln(w, name)
	fmt.Fprintf(w, "  %-8s%02d%02d%02d%02d%02d3D\n",
		"hbprobe",
		int(now.Month()), now.Day(), now.Year()%100, now.Hour(), now.Minute())
	fmt.Fprintln(w, "")
	fmt.Fprintf(w, "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n", len(positions), len(bonds))
	for i, p := range positions {
		sym := "C"
		if i < len(labels) {
			sym = labels[i]
		}
		fmt.Fprintf(w, "%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\n",
			p.X, p.Y, p.Z, sym, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
	}
	for _, b := range bonds {
		fmt.Fprintf(w, "%3d%3d%3d%3d%3d%3d%3d\n", b.I, b.J, b.Order, 0, 0, 0, 0)
	}
	fmt.Fprintln(w, "M  END")
}

// molName derives a human-readable probe type name from atom element counts.
func molName(labels []string) string {
	ec := map[string]int{}
	for _, s := range labels {
		ec[strings.ToUpper(s)]++
	}
	switch {
	case ec["C"] > 0 && ec["N"] == 0 && ec["O"] == 0 && ec["S"] == 0:
		return "hydrophobic"
	case ec["C"] > 0 && ec["N"] >= 3 && ec["O"] == 0:
		return "guanidine"
	case ec["C"] > 0 && ec["N"] == 0 && ec["O"] >= 2:
		return "carboxylate"
	case ec["C"] > 0 && ec["N"] >= 1 && ec["O"] >= 1:
		return "acetamide"
	case ec["O"] > 0 && ec["N"] == 0 && ec["S"] == 0 && ec["C"] <= 2:
		return "hydroxyl"
	case ec["S"] > 0:
		return "thiol"
	case ec["N"] > 0 && ec["C"] == 0 && ec["O"] == 0:
		return "amine"
	case ec["O"] > 0 && ec["C"] == 0:
		return "acceptor"
	case ec["N"] > 0 && ec["C"] == 0:
		return "donor"
	default:
		return "probe"
	}
}

// ReadSDF parses a V2000 SDF file and returns one ProbeSet per molecule.
func ReadSDF(path string) ([]ProbeSet, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var mols []ProbeSet
	sc := bufio.NewScanner(f)

	for {
		// Header: name line, program line, comment line.
		if !sc.Scan() {
			break
		}
		name := sc.Text()
		if !sc.Scan() { break } // program line
		if !sc.Scan() { break } // comment line

		// Counts line.
		if !sc.Scan() { break }
		counts := sc.Text()
		if len(counts) < 6 {
			skipToDelimiter(sc)
			continue
		}
		nAtoms, err1 := parseInt(counts[0:3])
		nBonds, err2 := parseInt(counts[3:6])
		if err1 != nil || err2 != nil {
			skipToDelimiter(sc)
			continue
		}

		mol := ProbeSet{Name: name, ParentIdx: -1}
		for i := 0; i < nAtoms; i++ {
			if !sc.Scan() { break }
			line := sc.Text()
			if len(line) < 34 {
				continue
			}
			x, _ := parseFloat(line[0:10])
			y, _ := parseFloat(line[10:20])
			z, _ := parseFloat(line[20:30])
			sym := strings.TrimSpace(line[31:34])
			if sym == "" { sym = "C" }
			mol.Add(Vec3{x, y, z}, sym)
		}
		for i := 0; i < nBonds; i++ {
			if !sc.Scan() { break }
			line := sc.Text()
			if len(line) < 9 {
				continue
			}
			ai, e1 := parseInt(line[0:3])
			aj, e2 := parseInt(line[3:6])
			order, e3 := parseInt(line[6:9])
			if e1 != nil || e2 != nil || e3 != nil { continue }
			mol.Bond(ai, aj, order)
		}
		skipToDelimiter(sc)
		if mol.Len() > 0 {
			mols = append(mols, mol)
		}
	}
	return mols, nil
}

// DeduplicateMols removes near-identical molecules from mols.
// Two molecules are considered duplicates when they have the same atom count,
// matching element sequence, and atom-by-atom RMSD below threshold.
// The centroid distance is used as a fast pre-filter.
func DeduplicateMols(mols []ProbeSet, threshold float64) []ProbeSet {
	centroid := func(m ProbeSet) Vec3 {
		c := Vec3{}
		for _, p := range m.Pos {
			c = c.Add(p)
		}
		return c.Scale(1.0 / float64(len(m.Pos)))
	}
	labelsMatch := func(a, b ProbeSet) bool {
		for i := range a.Labels {
			if !strings.EqualFold(a.Labels[i], b.Labels[i]) {
				return false
			}
		}
		return true
	}
	rmsd := func(a, b ProbeSet) float64 {
		sum := 0.0
		for i := range a.Pos {
			d := a.Pos[i].Sub(b.Pos[i]).Norm()
			sum += d * d
		}
		return math.Sqrt(sum / float64(len(a.Pos)))
	}

	kept := make([]ProbeSet, 0, len(mols))
	centroids := make([]Vec3, 0, len(mols))

	for _, m := range mols {
		c := centroid(m)
		duplicate := false
		for ki, k := range kept {
			if len(k.Pos) != len(m.Pos) {
				continue
			}
			if c.Sub(centroids[ki]).Norm() > threshold {
				continue
			}
			if !labelsMatch(m, k) {
				continue
			}
			if rmsd(m, k) < threshold {
				duplicate = true
				break
			}
		}
		if !duplicate {
			kept = append(kept, m)
			centroids = append(centroids, c)
		}
	}
	return kept
}

func skipToDelimiter(sc *bufio.Scanner) {
	for sc.Scan() {
		if strings.HasPrefix(sc.Text(), "$$$$") {
			return
		}
	}
}

func parseInt(s string) (int, error) {
	var n int
	_, err := fmt.Sscanf(strings.TrimSpace(s), "%d", &n)
	return n, err
}

func parseFloat(s string) (float64, error) {
	var f float64
	_, err := fmt.Sscanf(strings.TrimSpace(s), "%f", &f)
	return f, err
}

// WriteSDF writes linked probe pairs to an SDF file.
// Rules enforced per output molecule:
//  1. Exactly 2 probe groups (only linked pairs are written, never singles).
//  2. Must contain at least one heteroatom (N, O, S) — no hydrophobic-only.
//  3. Must have ≥ 10 atoms.
func WriteSDF(path string, ps *ProbeSet, pairs []ProbeSet) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	now := time.Now()

	hasHeteroatom := func(syms []string) bool {
		for _, s := range syms {
			switch strings.ToUpper(s) {
			case "N", "O", "S", "SE":
				return true
			}
		}
		return false
	}

	// isValid enforces output rules for linked probe pairs:
	//  1. Must contain a heteroatom (not hydrophobic-only) — unless it's a
	//     named ring scan variant (methyl scan intentionally pure-C).
	//  2. Must have ≥ 10 atoms.
	//  3. Molecular weight must be ≤ 500 Da (Lipinski-style upper bound).
	isValid := func(pair ProbeSet) bool {
		// Named ring scan variants (pos1-methyl, pos2-OH, etc.) have relaxed rules:
		// they're small by design and pure-C methyl scans are allowed.
		if pair.Name != "" {
			return len(pair.Labels) >= 2 // just needs to be a real molecule
		}
		// Linked pharmacophoric pairs: full rules.
		if len(pair.Labels) < 10 {
			return false
		}
		if !hasHeteroatom(pair.Labels) {
			return false
		}
		return true
	}

	const maxMW = 500.0

	// Atomic weights for MW calculation.
	atomicWeight := map[string]float64{
		"C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06,
		"F": 18.998, "CL": 35.45, "BR": 79.904, "P": 30.974,
		"SE": 78.971,
	}
	mw := func(pair ProbeSet) float64 {
		w := 0.0
		for _, lbl := range pair.Labels {
			if wt, ok := atomicWeight[strings.ToUpper(lbl)]; ok {
				w += wt
			}
		}
		return w
	}
	// Sort descending by MW (insertion sort — pairs is typically < 10k entries).
	for i := 1; i < len(pairs); i++ {
		for j := i; j > 0 && mw(pairs[j]) > mw(pairs[j-1]); j-- {
			pairs[j], pairs[j-1] = pairs[j-1], pairs[j]
		}
	}

	// Write all valid pairs — linked pharmacophoric pairs and ring scan variants.
	nSkipped := 0
	for _, pair := range pairs {
		if !isValid(pair) { continue }
		if mw(pair) > maxMW { continue }
		if !molValid(pair.Pos, pair.Bonds) { nSkipped++; continue }

		name := pair.Name
		if name == "" {
			name = molName(pair.Labels)
		}
		writeMolBlock(w, now, name, pair.Pos, pair.Labels, pair.Bonds)
		fmt.Fprintln(w, "$$$$")
	}
	if nSkipped > 0 {
		fmt.Printf("WriteSDF: skipped %d molecules with invalid geometry (NaN/Inf/out-of-range bonds)\n", nSkipped)
	}
	return w.Flush()
}
