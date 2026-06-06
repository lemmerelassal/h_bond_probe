package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"time"
)

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

// WriteSDF writes the probe set to an SDF file (one entry per connected component).
// The file extension is .sdf; the format is MOL V2000 with $$$$ record separators.
func WriteSDF(path string, ps *ProbeSet) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	now := time.Now()

	n := len(ps.Pos)
	adj := make([][]int, n)
	for _, b := range ps.Bonds {
		i, j := b.I-1, b.J-1
		if i >= 0 && i < n && j >= 0 && j < n {
			adj[i] = append(adj[i], j)
			adj[j] = append(adj[j], i)
		}
	}

	// BFS connected-component labelling.
	comp := make([]int, n)
	for i := range comp {
		comp[i] = -1
	}
	numComp := 0
	for start := 0; start < n; start++ {
		if comp[start] >= 0 {
			continue
		}
		comp[start] = numComp
		queue := []int{start}
		for len(queue) > 0 {
			cur := queue[0]
			queue = queue[1:]
			for _, nb := range adj[cur] {
				if comp[nb] < 0 {
					comp[nb] = numComp
					queue = append(queue, nb)
				}
			}
		}
		numComp++
	}

	// One SDF entry per component.
	for c := 0; c < numComp; c++ {
		remap := map[int]int{}
		var cPos []Vec3
		var cSym []string
		for i, p := range ps.Pos {
			if comp[i] != c {
				continue
			}
			remap[i] = len(cPos) + 1
			sym := "C"
			if i < len(ps.Labels) {
				sym = ps.Labels[i]
			}
			cPos = append(cPos, p)
			cSym = append(cSym, sym)
		}
		var cBonds []Bond
		for _, b := range ps.Bonds {
			if comp[b.I-1] != c || comp[b.J-1] != c {
				continue
			}
			cBonds = append(cBonds, Bond{remap[b.I-1], remap[b.J-1], b.Order})
		}
		writeMolBlock(w, now, molName(cSym), cPos, cSym, cBonds)
		fmt.Fprintln(w, "$$$$")
	}
	return w.Flush()
}
