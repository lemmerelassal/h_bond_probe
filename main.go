package main

import (
	"flag"
	"fmt"
	"math"
	"os"
	"sort"
	"strings"
	"time"
)

// noCarbonyl disables placement of carbonyl-containing groups (ketone, amide,
// urea, guanidine, thiourea, imine, etc.) during generation.
// Molecules imported via -input are unaffected.
var noCarbonyl bool

// multiFlag is a repeatable string flag (e.g. -input a.sdf -input b.sdf).
type multiFlag []string

func (m *multiFlag) String() string        { return strings.Join(*m, ",") }
func (m *multiFlag) Set(v string) error    { *m = append(*m, v); return nil }

func main() {
	alphafold := flag.Bool("alphafold", false, "Fetch structure from AlphaFold DB instead of RCSB PDB")
	startRes := flag.Int("start", 0, "First residue number to include (0 = all)")
	endRes := flag.Int("end", 0, "Last residue number to include (0 = all)")
	flag.BoolVar(&noCarbonyl, "nocarbonyl", false, "Skip placing carbonyl-containing groups (urea, amide, guanidine, etc.); -input molecules are unaffected")
	var inputSDFs multiFlag
	flag.Var(&inputSDFs, "input", "SDF file of existing probes to import (may be repeated)")
	flag.Parse()

	args := flag.Args()
	if len(args) < 1 {
		fmt.Fprintln(os.Stderr, "Usage: pharmacophore [flags] <PDB-ID|UniProt-ID|file.pdb|file.cif>")
		fmt.Fprintln(os.Stderr, "Flags:")
		flag.PrintDefaults()
		os.Exit(1)
	}

	atoms, pdbID, err := LoadAtoms(args[0], *alphafold)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
	if len(atoms) == 0 {
		fmt.Fprintln(os.Stderr, "No atoms parsed.")
		os.Exit(1)
	}

	// Apply residue range filter if requested.
	if *startRes > 0 || *endRes > 0 {
		var filtered []Atom
		for _, a := range atoms {
			if *startRes > 0 && a.ResSeq < *startRes {
				continue
			}
			if *endRes > 0 && a.ResSeq > *endRes {
				continue
			}
			filtered = append(filtered, a)
		}
		fmt.Printf("Residue range %d–%d: %d → %d atoms\n", *startRes, *endRes, len(atoms), len(filtered))
		atoms = filtered
	}
	if len(atoms) == 0 {
		fmt.Fprintln(os.Stderr, "No atoms in specified residue range.")
		os.Exit(1)
	}

	// ── Partition atoms ───────────────────────────────────────────────────────
	// heavy: all non-H, non-water atoms (used for clash checking, includes HET)
	// hydro: protein hydrophobic C atoms only (source for vote system)
	// proteinAtoms: protein-only atom slice fed to all Placers
	var hydro []Atom
	var heavy []heavyAtom
	var proteinAtoms []Atom
	for _, a := range atoms {
		if strings.ToUpper(a.Element) == "H" || isWater(a) {
			continue
		}
		heavy = append(heavy, heavyAtom{a.Pos, vdw(a.Element)})
		if a.IsHet {
			continue
		}
		proteinAtoms = append(proteinAtoms, a)
		if isHydrophobic(a) {
			hydro = append(hydro, a)
		}
	}
	fmt.Printf("Hydrophobic atoms: %d  Heavy atoms total: %d\n", len(hydro), len(heavy))

	// ── Initialise the probe accumulator ─────────────────────────────────────
	ps := &ProbeSet{}

	for _, path := range inputSDFs {
		mols, err := ReadSDF(path)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error reading input SDF: %v\n", err)
			os.Exit(1)
		}
		before := len(mols)
		mols = DeduplicateMols(mols, 1.0)
		fmt.Printf("Imported %d molecules from %s (%d duplicates removed)\n", len(mols), path, before-len(mols))
		for _, mol := range mols {
			offset := ps.Len()
			for i, p := range mol.Pos {
				ps.Add(p, mol.Labels[i])
			}
			for _, b := range mol.Bonds {
				ps.Bond(offset+b.I, offset+b.J, b.Order)
			}
		}
	}

	t0 := time.Now()
	tick := func(label string) {
		fmt.Printf("  %-40s %v\n", label, time.Since(t0).Round(time.Millisecond))
		t0 = time.Now()
	}

	// ── Alkane chain through hydrophobic intersection points ──────────────────
	buildAlkaneChain(hydro, heavy, ps)
	tick("alkane chain")

	// ── Benzene / imidazole ring probes for aromatic residues ─────────────────
	buildAromaticProbes(proteinAtoms, heavy, ps)
	tick("aromatic probes")

	// ── H-bond acceptor / donor vote probes ──────────────────────────────────
	buildHbondProbes(proteinAtoms, heavy, ps)
	tick("h-bond probes")

	// ── Pharmacophoric probes (Open/Closed: iterate AllPlacers) ───────────────
	for _, p := range AllPlacers {
		n := p.Place(proteinAtoms, heavy, ps)
		fmt.Printf("%s probes: %d\n", p.Name(), n)
	}
	tick("pharmacophoric probes")

	// ── Ring SAR scan ─────────────────────────────────────────────────────────
	fmt.Println("Running ring SAR scan...")
	ringScanPairs := ScanRings(ps, heavy)
	fmt.Printf("Ring scan variants: %d\n", len(ringScanPairs))
	tick("ring SAR scan")

	// ── Link probe groups into independent pair molecules ─────────────────────
	// Each probe is used at most once (greedy nearest-neighbor matching).
	// The linker chain is capped at 6 carbons (≈10 Å end-to-end).
	linkableN := ps.Len()
	pairs := linkProbeGroups(&ps.Pos, &ps.Labels, &ps.Bonds, heavy, 14.0, linkableN)
	fmt.Printf("Linker chains: %d\n", len(pairs))
	tick("link probe groups")

	merged := linkMoleculeGroups(pairs, heavy, 14.0)
	fmt.Printf("Merged molecules: %d\n", len(merged))
	pairs = append(pairs, merged...)
	tick("link molecule groups")

	// Amide variants: replace adjacent sp3 chain-C pairs with C=O / N.
	amideGrid := newClashGrid(heavy)
	var amidePairs []ProbeSet
	for _, mol := range pairs {
		amidePairs = append(amidePairs, amideVariants(mol, amideGrid)...)
	}
	fmt.Printf("Amide variants: %d\n", len(amidePairs))
	pairs = append(pairs, amidePairs...)
	tick("amide variants")

	pharmPairs := pairs

	// ── Ring SAR scan on linked pairs ─────────────────────────────────────────
	var linkedRingScan []ProbeSet
	for pi, pp := range pharmPairs {
		variants := ScanRings(&pp, heavy)
		for i := range variants {
			variants[i].ParentIdx = pi
		}
		linkedRingScan = append(linkedRingScan, variants...)
	}
	fmt.Printf("Ring-substituted linked pairs: %d\n", len(linkedRingScan))
	tick("ring scan on linked pairs")

	// ── Chain heteroatom scan ─────────────────────────────────────────────────
	chainScanned := ScanChain(pharmPairs, heavy)
	fmt.Printf("Chain heteroatom variants: %d\n", len(chainScanned))
	tick("chain heteroatom scan")

	// ── Fragment growing ────────────────────────────────────────────────────────
	grownPairs := GrowLinked(pharmPairs, heavy)
	fmt.Printf("Grown variants: %d\n", len(grownPairs))
	grownRingScan := GrowLinked(linkedRingScan, heavy)
	fmt.Printf("Grown ring-substituted: %d\n", len(grownRingScan))
	grownPairs = append(grownPairs, grownRingScan...)
	tick("fragment growing")

	// ── Ring closures ────────────────────────────────────────────────────────────
	closedPairs := RingClose(pharmPairs, heavy)
	closedPairs = append(closedPairs, RingClose(grownPairs, heavy)...)
	fmt.Printf("Ring closures: %d\n", len(closedPairs))
	tick("ring closures")

	// ── Bidirectional growing ─────────────────────────────────────────────────
	// NOTE: BidirGrow builds a second bridge between probe endpoints, forming a
	// macrocycle on top of the existing linker. This only makes sense for short
	// gaps (< 8 Å). Disabled until pre-linking architecture is implemented.
	// bidirPairs := BidirGrow(pharmPairs, heavy)
	bidirPairs := []ProbeSet{}
	fmt.Printf("Bidirectional grown: %d\n", len(bidirPairs))

	// selfClashFree: grown atoms must not overlap each other or scaffold.
	selfClashFree := func(v ProbeSet, nOrigAtoms int) bool {
		for i := nOrigAtoms; i < len(v.Pos); i++ {
			ri := vdw(v.Labels[i])
			for j := 0; j < i; j++ {
				bonded := false
				for _, b := range v.Bonds {
					if (b.I-1 == i && b.J-1 == j) || (b.I-1 == j && b.J-1 == i) {
						bonded = true
						break
					}
				}
				if bonded {
					continue
				}
				rj := vdw(v.Labels[j])
				if v.Pos[i].Sub(v.Pos[j]).Norm() < ri+rj-1.0 {
					return false
				}
			}
		}
		return true
	}

	// Keep all grown variants, ring closures, ring-substituted linked pairs.
	allVariants := append(append(append(grownPairs, closedPairs...), bidirPairs...), linkedRingScan...)
	validVariants := make([]ProbeSet, 0, len(pharmPairs)+len(allVariants))
	// Always include the original pharmacophoric pairs.
	validVariants = append(validVariants, pharmPairs...)
	// Add all grown/closed variants that pass self-clash check.
	for _, v := range allVariants {
		pi := v.ParentIdx
		if pi < 0 || pi >= len(pharmPairs) {
			continue
		}
		nOrig := len(pharmPairs[pi].Pos)
		if strings.HasPrefix(v.Name, "bidir-") || selfClashFree(v, nOrig) {
			validVariants = append(validVariants, v)
		}
	}
	fmt.Printf("Total variants (all): %d\n", len(validVariants))
	tick("self-clash filter")

	// ── Deduplicate by linker RMSD ────────────────────────────────────────────
	// Two variants with the same ParentIdx and same atom count are duplicates
	// if the RMSD of their linker atoms (beyond nOrig scaffold atoms) is < 1.0 Å.
	// Keep the one with fewer protein heavy-atom clashes.
	clashCount := func(ps ProbeSet) int {
		n := 0
		for _, p := range ps.Pos {
			for _, h := range heavy {
				if p.Sub(h.pos).Norm() < h.vdwR+vdw("C")-0.8 {
					n++
				}
			}
		}
		return n
	}
	linkerRMSD := func(a, b ProbeSet, start int) float64 {
		if len(a.Pos) != len(b.Pos) || start >= len(a.Pos) {
			return math.MaxFloat64
		}
		sum := 0.0
		n := len(a.Pos) - start
		for i := start; i < len(a.Pos); i++ {
			d := a.Pos[i].Sub(b.Pos[i]).Norm()
			sum += d * d
		}
		return math.Sqrt(sum / float64(n))
	}
	const rmsdThreshold = 1.0 // Å
	keep := make([]bool, len(validVariants))
	for i := range keep { keep[i] = true }
	for i := 0; i < len(validVariants); i++ {
		if !keep[i] { continue }
		vi := validVariants[i]
		pi := vi.ParentIdx
		nOrig := 0
		if pi >= 0 && pi < len(pharmPairs) { nOrig = len(pharmPairs[pi].Pos) }
		for j := i + 1; j < len(validVariants); j++ {
			if !keep[j] { continue }
			vj := validVariants[j]
			if vj.ParentIdx != pi { continue }
			if len(vj.Pos) != len(vi.Pos) { continue }
			if linkerRMSD(vi, vj, nOrig) < rmsdThreshold {
				// Keep the one with fewer clashes.
				if clashCount(vj) < clashCount(vi) {
					keep[i] = false
					break
				}
				keep[j] = false
			}
		}
	}
	deduped := validVariants[:0]
	for i, v := range validVariants {
		if keep[i] { deduped = append(deduped, v) }
	}
	nDupes := len(validVariants) - len(deduped)
	if nDupes > 0 {
		fmt.Printf("Deduplicated: removed %d near-identical variants (RMSD < %.1f Å)\n", nDupes, rmsdThreshold)
	}
	validVariants = deduped
	tick("RMSD deduplication")

	// Combine: all valid variants + standalone ring scan variants.
	pairs = append(validVariants, ringScanPairs...)

	// ── Backbone ──────────────────────────────────────────────────────────────
	nBB := addBackbone(proteinAtoms, &ps.Pos, &ps.Labels, &ps.Bonds)
	fmt.Printf("Backbone residues: %d\n", nBB)

	// ── Write output ──────────────────────────────────────────────────────────
	outPath := pdbID + "_probes.sdf"
	if len(pairs) == 0 && len(ringScanPairs) == 0 {
		fmt.Fprintln(os.Stderr, "No probe pairs generated — nothing to write.")
		os.Exit(1)
	}
	// Combine: pharmacophoric pairs + standalone ring scan variants (those not
	// picked up by linking get written as-is for visual inspection).
	allPairs := append(append(pairs, ringScanPairs...), chainScanned...)

	// Fuse atoms that are within 0.7 Å of each other (overlapping ring/carbonyl
	// atoms), then sanitize any valence overflows introduced by fusion. Fusion is
	// opportunistic: it keeps one of the two overlapping atoms' positions, so it
	// can leave the inherited bonds at strained angles. We accept a fusion only
	// when the result is geometrically sane; otherwise we keep the unfused
	// molecule (which the geometry gate below will still vet).
	nFused := 0
	for i := range allPairs {
		base := sanitizeValence(allPairs[i])
		fused := sanitizeValence(fuseOverlapping(allPairs[i], 0.7))
		if len(fused.Pos) < len(allPairs[i].Pos) && geometryOK(fused) {
			allPairs[i] = fused
			nFused++
		} else {
			allPairs[i] = base
		}
	}
	if nFused > 0 {
		fmt.Printf("Fused overlapping atoms in %d molecules\n", nFused)
	}
	tick("fuse + sanitize valence")

	// Keep only fully connected molecules.
	// Both probe endpoints are already required to contain a ring (noRing filter
	// in linkProbeGroups), so the ring count requirement is satisfied by
	// construction. Grown/closed/chain-scanned variants are kept as long as they
	// are connected.
	var connPairs []ProbeSet
	nBadGeom := 0
	for _, p := range allPairs {
		if !isConnected(p) {
			continue
		}
		if !geometryOK(p) {
			nBadGeom++
			continue
		}
		connPairs = append(connPairs, p)
	}
	fmt.Printf("Connected molecules: %d / %d (dropped %d for strained geometry)\n",
		len(connPairs), len(allPairs), nBadGeom)
	allPairs = connPairs
	tick("connectivity + geometry filter")

	// Surface the individual probes as standalone molecules — otherwise small
	// probes (hydroxyl, methoxy) that are never linked never appear at all.
	standaloneProbes := probeComponents(ps.Pos, ps.Labels, ps.Bonds, linkableN)
	// Hold standalone probes to the same geometry bar as everything else; this
	// drops the hydrophobic alkane-chain probes, which connect hydrophobic hotspots
	// rather than forming a tetrahedral chain and so contain strained angles.
	keptProbes := standaloneProbes[:0]
	for _, sp := range standaloneProbes {
		if geometryOK(sp) {
			keptProbes = append(keptProbes, sp)
		}
	}
	standaloneProbes = keptProbes
	fmt.Printf("Standalone probe molecules: %d\n", len(standaloneProbes))
	allPairs = append(allPairs, standaloneProbes...)

	if err := WriteSDF(outPath, ps, allPairs); err != nil {
		fmt.Fprintf(os.Stderr, "Error writing SDF: %v\n", err)
		os.Exit(1)
	}
	fmt.Printf("Wrote %d probe atoms + %d linked pairs + %d ring variants + %d standalone probes -> %s\n",
		ps.Len(), len(pairs), len(ringScanPairs), len(standaloneProbes), outPath)
	tick("write SDF")
}

// sanitizeValence removes excess bonds from over-valent atoms.
// Cyclic (ring) bonds are processed first so their valence budget is reserved —
// otherwise an extra bond introduced by linking/fusion/growing could over-fill a
// ring atom and cause a ring bond to be dropped, breaking the ring. Within that,
// higher-order bonds are processed first so double bonds (C=O, C=N) are preserved
// over single bonds when an atom is over its element maximum.
func sanitizeValence(mol ProbeSet) ProbeSet {
	maxVal := map[string]int{"C": 4, "N": 3, "O": 2, "S": 6, "F": 1, "CL": 1, "BR": 1}
	n := len(mol.Pos)

	// A bond is cyclic if both endpoints survive in the 2-core of the bond graph
	// (degree-1 atoms iteratively pruned). This protects every true ring bond.
	adj := make([][]int, n)
	for _, b := range mol.Bonds {
		i, j := b.I-1, b.J-1
		if i < 0 || j < 0 || i >= n || j >= n {
			continue
		}
		adj[i] = append(adj[i], j)
		adj[j] = append(adj[j], i)
	}
	inCore := ringCoreAtoms(adj)
	cyclic := func(b Bond) bool {
		i, j := b.I-1, b.J-1
		return i >= 0 && j >= 0 && i < n && j < n && inCore[i] && inCore[j]
	}

	type idxBond struct {
		b   Bond
		idx int
	}
	sorted := make([]idxBond, len(mol.Bonds))
	for k, b := range mol.Bonds {
		sorted[k] = idxBond{b, k}
	}
	// Ring bonds first (reserve their budget), then highest-order first.
	sort.Slice(sorted, func(i, j int) bool {
		ci, cj := cyclic(sorted[i].b), cyclic(sorted[j].b)
		if ci != cj {
			return ci
		}
		return sorted[i].b.Order > sorted[j].b.Order
	})

	used := make([]int, n) // accumulated bond-order sum per atom
	keep := make([]bool, len(mol.Bonds))
	for _, ib := range sorted {
		b := ib.b
		i, j := b.I-1, b.J-1
		if i < 0 || j < 0 || i >= n || j >= n {
			continue
		}
		ord := b.Order
		if ord == 4 {
			ord = 2 // aromatic bond consumes 2 of each atom's valence budget
		}
		elI := strings.ToUpper(mol.Labels[i])
		elJ := strings.ToUpper(mol.Labels[j])
		mvI, okI := maxVal[elI]
		if !okI {
			mvI = 4
		}
		mvJ, okJ := maxVal[elJ]
		if !okJ {
			mvJ = 4
		}
		if used[i]+ord <= mvI && used[j]+ord <= mvJ {
			used[i] += ord
			used[j] += ord
			keep[ib.idx] = true
		}
	}

	result := ProbeSet{Name: mol.Name, ParentIdx: mol.ParentIdx}
	for k, p := range mol.Pos {
		result.Add(p, mol.Labels[k])
	}
	for k, b := range mol.Bonds {
		if keep[k] {
			result.Bond(b.I, b.J, b.Order)
		}
	}
	return result
}

// fuseOverlapping merges atom pairs within thresh Å into a single atom, redirecting
// all bonds to the surviving atom. This produces fused ring systems and
// ring–carbonyl fusions from molecules whose atoms were placed at overlapping positions.
func fuseOverlapping(mol ProbeSet, thresh float64) ProbeSet {
	n := len(mol.Pos)
	parent := make([]int, n)
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		if parent[x] != x {
			parent[x] = find(parent[x])
		}
		return parent[x]
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if mol.Pos[i].Sub(mol.Pos[j]).Norm() < thresh {
				pi, pj := find(i), find(j)
				if pi != pj {
					parent[pi] = pj
				}
			}
		}
	}

	// Map each atom to its canonical new index.
	newIdx := make([]int, n)
	repr := map[int]int{}
	count := 0
	for i := 0; i < n; i++ {
		r := find(i)
		if _, ok := repr[r]; !ok {
			repr[r] = count
			count++
		}
		newIdx[i] = repr[r]
	}

	// Build new atom list using the first atom seen for each group.
	result := ProbeSet{Name: mol.Name, ParentIdx: mol.ParentIdx}
	added := make([]bool, count)
	for i := 0; i < n; i++ {
		ni := newIdx[i]
		if !added[ni] {
			result.Add(mol.Pos[i], mol.Labels[i])
			added[ni] = true
		}
	}

	// Remap bonds; drop self-loops from fusion; keep max order for duplicates.
	type bkey struct{ i, j int }
	best := map[bkey]int{}
	for _, b := range mol.Bonds {
		ni := newIdx[b.I-1] + 1
		nj := newIdx[b.J-1] + 1
		if ni == nj {
			continue
		}
		k := bkey{ni, nj}
		if ni > nj {
			k = bkey{nj, ni}
		}
		if b.Order > best[k] {
			best[k] = b.Order
		}
	}
	for k, ord := range best {
		result.Bond(k.i, k.j, ord)
	}
	return result
}

// isConnected returns true if all atoms in mol are reachable from atom 0 via bonds.
func isConnected(mol ProbeSet) bool {
	n := len(mol.Pos)
	if n <= 1 {
		return true
	}
	adj := make([][]int, n)
	for _, b := range mol.Bonds {
		i, j := b.I-1, b.J-1
		if i < 0 || j < 0 || i >= n || j >= n {
			continue
		}
		adj[i] = append(adj[i], j)
		adj[j] = append(adj[j], i)
	}
	visited := make([]bool, n)
	queue := []int{0}
	visited[0] = true
	count := 1
	for len(queue) > 0 {
		cur := queue[0]
		queue = queue[1:]
		for _, nb := range adj[cur] {
			if !visited[nb] {
				visited[nb] = true
				count++
				queue = append(queue, nb)
			}
		}
	}
	return count == n
}

// probeComponents splits the first `limit` atoms of a probe set (the
// pharmacophore probes, before the backbone is appended) into their connected
// components and returns each as a standalone molecule named by molName. This is
// what surfaces the individual probes — especially the small ones (hydroxyl,
// methoxy) that are too few atoms to ever be linked into a larger molecule and so
// would otherwise never appear in the output at all.
func probeComponents(pos []Vec3, labels []string, bonds []Bond, limit int) []ProbeSet {
	if limit > len(pos) {
		limit = len(pos)
	}
	parent := make([]int, limit)
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	for _, b := range bonds {
		i, j := b.I-1, b.J-1
		if i < 0 || j < 0 || i >= limit || j >= limit {
			continue
		}
		parent[find(i)] = find(j)
	}
	groups := map[int][]int{}
	for i := 0; i < limit; i++ {
		r := find(i)
		groups[r] = append(groups[r], i)
	}
	// Deterministic order: by smallest member index.
	roots := make([]int, 0, len(groups))
	for r := range groups {
		roots = append(roots, r)
	}
	sort.Ints(roots)

	var out []ProbeSet
	for _, r := range roots {
		members := groups[r]
		idxMap := make(map[int]int, len(members))
		var ps ProbeSet
		for _, m := range members {
			idxMap[m] = ps.Add(pos[m], labels[m])
		}
		for _, b := range bonds {
			i, j := b.I-1, b.J-1
			if ni, ok := idxMap[i]; ok {
				if nj, ok2 := idxMap[j]; ok2 {
					ps.Bond(ni, nj, b.Order)
				}
			}
		}
		name := molName(ps.Labels)
		// Distinguish a bridging-O ether (methoxy, C–O–C) from a terminal-O
		// alcohol (hydroxyl, C–OH); molName lumps both as "hydroxyl".
		if name == "hydroxyl" {
			deg := make([]int, ps.Len())
			for _, b := range ps.Bonds {
				deg[b.I-1]++
				deg[b.J-1]++
			}
			for i, l := range ps.Labels {
				if strings.ToUpper(l) == "O" && deg[i] == 2 {
					name = "methoxy"
					break
				}
			}
		}
		ps.Name = "probe-" + name
		ps.ParentIdx = -1
		out = append(out, ps)
	}
	return out
}

// covRadius holds covalent radii (Å) used to validate bond lengths.
var covRadius = map[string]float64{
	"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05,
	"F": 0.57, "CL": 1.02, "BR": 1.20, "P": 1.07, "SE": 1.20,
}

func covR(elem string) float64 {
	if r, ok := covRadius[strings.ToUpper(elem)]; ok {
		return r
	}
	return 0.77
}

// geometryOK rejects molecules whose downstream bond formation (linking, ring
// closure, fusion, growth) produced unphysical geometry: a bond far longer than
// the covalent ideal, or a bond angle so acute it cannot be a real conformation.
// These bonds are added purely on inter-atomic distance with no angular control,
// so this gate is what keeps strained artifacts out of the output.
func geometryOK(mol ProbeSet) bool {
	const (
		minAngleDeg = 82.0 // below this, the angle is unphysical (4-rings ≈ 90° pass)
		lenTol      = 0.40 // Å slack above the covalent-radii sum for a single bond
	)
	n := len(mol.Pos)
	adj := make([][]int, n)
	for _, b := range mol.Bonds {
		i, j := b.I-1, b.J-1
		if i < 0 || j < 0 || i >= n || j >= n {
			continue
		}
		d := mol.Pos[i].Sub(mol.Pos[j]).Norm()
		// Order 2/3 bonds are shorter, so the single-bond ideal is a safe upper bound.
		if d > covR(mol.Labels[i])+covR(mol.Labels[j])+lenTol {
			return false
		}
		adj[i] = append(adj[i], j)
		adj[j] = append(adj[j], i)
	}
	cosMax := math.Cos(minAngleDeg * math.Pi / 180)
	// A divalent oxygen (ether/ester) is sp3-bent; ≥150° is unphysical (linear O).
	cosOLinear := math.Cos(150.0 * math.Pi / 180)
	for c := 0; c < n; c++ {
		nb := adj[c]
		isDivalentO := strings.ToUpper(mol.Labels[c]) == "O" && len(nb) == 2
		for a := 0; a < len(nb); a++ {
			for b := a + 1; b < len(nb); b++ {
				u := mol.Pos[nb[a]].Sub(mol.Pos[c])
				v := mol.Pos[nb[b]].Sub(mol.Pos[c])
				du, dv := u.Norm(), v.Norm()
				if du < 1e-6 || dv < 1e-6 {
					return false
				}
				cosA := u.Dot(v) / (du * dv)
				if cosA > cosMax { // angle < minAngleDeg
					return false
				}
				if isDivalentO && cosA < cosOLinear { // angle > 150° on an ether O
					return false
				}
			}
		}
	}
	return true
}

// ── Hydrophobic alkane chain ──────────────────────────────────────────────────

func buildAlkaneChain(hydro []Atom, heavy []heavyAtom, ps *ProbeSet) {
	votes := map[[3]int]*voteCell{}
	nPairs := 0
	for i := 0; i < len(hydro); i++ {
		for j := i + 1; j < len(hydro); j++ {
			if hydro[i].Pos.Sub(hydro[j].Pos).Norm() > maxPairDist {
				continue
			}
			pts := sphereCircle(hydro[i].Pos, hydro[j].Pos, sphereR, nSamples)
			if len(pts) == 0 {
				continue
			}
			nPairs++
			pk := [2]int{i, j}
			for _, p := range pts {
				if !noHardClash(p, vdw("C"), hardTol, heavy) {
					continue
				}
				key := gridKey(p)
				if e, ok := votes[key]; ok {
					if !e.pairs[pk] {
						e.pairs[pk] = true
						e.pairCount++
					}
					e.nPts++
					n := float64(e.nPts)
					e.pos = Vec3{
						(e.pos.X*(n-1) + p.X) / n,
						(e.pos.Y*(n-1) + p.Y) / n,
						(e.pos.Z*(n-1) + p.Z) / n,
					}
				} else {
					votes[key] = &voteCell{
						pairCount: 1, pos: p, nPts: 1,
						pairs: map[[2]int]bool{pk: true},
					}
				}
			}
		}
	}
	fmt.Printf("Pairs within %.1f Å: %d  Vote cells: %d\n", maxPairDist, nPairs, len(votes))

	// burialCount returns the number of protein heavy atoms within shellR Å
	// of pos that are also beyond minDist Å (excluding the immediate vdW shell).
	// High burial = deep in a pocket; low burial = exposed on surface.
	const (
		burialShellR  = 8.0 // Å — count atoms within this radius
		burialMinDist = 2.5 // Å — ignore atoms too close (vdW contacts)
		minBurial     = 6   // minimum surrounding atoms to accept as pocket
	)
	burialCount := func(pos Vec3) int {
		n := 0
		for _, h := range heavy {
			d := pos.Sub(h.pos).Norm()
			if d >= burialMinDist && d <= burialShellR {
				n++
			}
		}
		return n
	}

	type candidate struct {
		pos       Vec3
		pairCount int
	}
	var cands []candidate
	for _, v := range votes {
		if v.pairCount >= minPairs {
			cands = append(cands, candidate{v.pos, v.pairCount})
		}
	}
	sort.Slice(cands, func(i, j int) bool { return cands[i].pairCount > cands[j].pairCount })

	var waypoints []Vec3
	nClash, nClose, nSurface := 0, 0, 0
	for _, c := range cands {
		if !noHardClash(c.pos, vdw("C"), hardTol, heavy) {
			nClash++
			continue
		}
		if burialCount(c.pos) < minBurial {
			nSurface++
			continue
		}
		tooClose := false
		for _, wp := range waypoints {
			if c.pos.Sub(wp).Norm() < probeSpacing {
				tooClose = true
				break
			}
		}
		if tooClose {
			nClose++
			continue
		}
		waypoints = append(waypoints, c.pos)
	}
	fmt.Printf("Placed: %d  Rejected (clash): %d  Rejected (too close): %d  Rejected (surface): %d\n",
		len(waypoints), nClash, nClose, nSurface)

	if len(waypoints) == 0 {
		fmt.Println("No hydrophobic vote cells — no alkane chain built.")
		return
	}
	chainAtoms := alkaneChain(waypoints, heavy)
	if len(chainAtoms) == 0 {
		return
	}
	total := 0
	for _, seg := range chainAtoms {
		base := ps.Len() + 1
		for _, p := range seg {
			ps.Add(p, "C")
		}
		for i := 0; i+1 < len(seg); i++ {
			ps.Bond(base+i, base+i+1, 1)
		}
		total += len(seg)
	}
	fmt.Printf("Alkane chain atoms: %d in %d segment(s)\n", total, len(chainAtoms))
}

// alkaneChain orders waypoints by nearest-neighbour and builds all-trans sp3
// chain segments threading through them. Returns a slice of segments; atoms
// within each segment are bonded consecutively, but NO bond is written between
// segments (skipped waypoints would otherwise create impossible long bonds).
func alkaneChain(waypoints []Vec3, heavy []heavyAtom) [][]Vec3 {
	const (
		alkBond = 1.54
		alkCosA = 0.8165
		alkSinA = 0.578
	)

	// Nearest-neighbour path.
	visited := make([]bool, len(waypoints))
	path := []int{0}
	visited[0] = true
	for len(path) < len(waypoints) {
		last := waypoints[path[len(path)-1]]
		bestIdx, bestD := -1, 999.0
		for i, p := range waypoints {
			if visited[i] {
				continue
			}
			if d := last.Sub(p).Norm(); d < bestD {
				bestD, bestIdx = d, i
			}
		}
		path = append(path, bestIdx)
		visited[bestIdx] = true
	}

	// Build one segment at a time, carrying prevDir for torsion continuity.
	buildSeg := func(tip, target Vec3, n int, prevDir Vec3) ([]Vec3, Vec3, float64) {
		axis := target.Sub(tip).unit()
		segLen := tip.Sub(target).Norm()

		basePerp := func() Vec3 {
			if prevDir.Norm() > 0.1 {
				comp := prevDir.Sub(axis.Scale(prevDir.Dot(axis)))
				if comp.Norm() > 0.3 {
					return comp.Scale(-1).unit()
				}
			}
			return perpendicular(axis)
		}()
		perp2 := cross3(axis, basePerp).unit()

		type pt2 struct{ u, v float64 }
		local := make([]pt2, n+2)
		for k := 1; k <= n+1; k++ {
			sign := 1.0
			if (k-1)%2 == 1 {
				sign = -1
			}
			local[k] = pt2{local[k-1].u + alkCosA*alkBond, local[k-1].v + sign*alkSinA*alkBond}
		}
		eu, ev := local[n+1].u, local[n+1].v
		idealLen := math.Sqrt(eu*eu + ev*ev)
		s := segLen / idealLen

		buildWith := func(perpV Vec3) []Vec3 {
			if s < 0.75 || s > 1.25 {
				return nil
			}
			cosθ, sinθ := eu/idealLen, ev/idealLen
			pts := make([]Vec3, n)
			for k := 1; k <= n; k++ {
				u, v := local[k].u, local[k].v
				ru := (u*cosθ + v*sinθ) * s
				rv := (-u*sinθ + v*cosθ) * s
				pts[k-1] = tip.Add(axis.Scale(ru)).Add(perpV.Scale(rv))
			}
			return pts
		}
		clashCt := func(pts []Vec3) int {
			c := 0
			for _, p := range pts {
				for _, h := range heavy {
					if p.Sub(h.pos).Norm() < h.vdwR+vdw("C")-hardTol {
						c++
						break
					}
				}
			}
			return c
		}
		juncDev := func(pts []Vec3) float64 {
			if len(pts) < 2 {
				return 0
			}
			v1 := pts[len(pts)-2].Sub(pts[len(pts)-1])
			v2 := target.Sub(pts[len(pts)-1])
			n1, n2 := v1.Norm(), v2.Norm()
			if n1 < 1e-6 || n2 < 1e-6 {
				return 0
			}
			c := v1.Dot(v2) / (n1 * n2)
			if c > 1 {
				c = 1
			}
			if c < -1 {
				c = -1
			}
			return math.Abs(180/math.Pi*math.Acos(c) - 109.5)
		}

		best := buildWith(basePerp)
		if best == nil {
			best = []Vec3{}
		}
		bestC := clashCt(best)
		bestAD := juncDev(best)
		for oi := 1; oi < 36; oi++ {
			ang := float64(oi) * math.Pi / 18
			pv := basePerp.Scale(math.Cos(ang)).Add(perp2.Scale(math.Sin(ang)))
			pts := buildWith(pv)
			if pts == nil {
				continue
			}
			c := clashCt(pts)
			ad := juncDev(pts)
			if c < bestC || (c == bestC && ad < bestAD) {
				bestC, bestAD, best = c, ad, pts
			}
		}
		var newDir Vec3
		if len(best) >= 2 {
			newDir = best[len(best)-1].Sub(best[len(best)-2]).unit()
		} else if len(best) == 1 {
			newDir = best[0].Sub(tip).unit()
		} else {
			newDir = axis
		}
		return best, newDir, bestAD
	}

	var segments [][]Vec3
	var currentSeg []Vec3
	currentSeg = append(currentSeg, waypoints[path[0]])
	prevDir := Vec3{}

	flushSeg := func() {
		if len(currentSeg) > 1 {
			segments = append(segments, currentSeg)
		}
		currentSeg = nil
	}

	for wi := 1; wi < len(path); wi++ {
		tip := currentSeg[len(currentSeg)-1]
		target := waypoints[path[wi]]
		segLen := tip.Sub(target).Norm()

		n := int(math.Round(segLen/1.258)) - 1
		if n < 2 {
			n = 2
		}
		bestN, bestDev := n, 999.0
		for _, tryN := range []int{n - 1, n, n + 1, n + 2} {
			if tryN < 2 {
				continue
			}
			type pt2 struct{ u, v float64 }
			local := make([]pt2, tryN+2)
			for k := 1; k <= tryN+1; k++ {
				sign := 1.0
				if (k-1)%2 == 1 {
					sign = -1
				}
				local[k] = pt2{local[k-1].u + alkCosA*alkBond, local[k-1].v + sign*alkSinA*alkBond}
			}
			eu, ev := local[tryN+1].u, local[tryN+1].v
			idealLen := math.Sqrt(eu*eu + ev*ev)
			if dev := math.Abs(segLen/idealLen - 1.0); dev < bestDev {
				bestDev, bestN = dev, tryN
			}
		}

		pts, newDir, ad := buildSeg(tip, target, bestN, prevDir)

		if (ad > 45 && len(pts) > 0) || len(pts) == 0 {
			// Bad junction or degenerate segment: flush the current segment
			// and start a fresh one from the target waypoint.  This prevents
			// any bond being written across the discontinuity.
			flushSeg()
			currentSeg = append(currentSeg, target)
			prevDir = Vec3{}
			continue
		}

		for _, p := range pts {
			currentSeg = append(currentSeg, p)
		}
		prevDir = newDir
	}
	flushSeg()
	return segments
}

// ── Aromatic ring probes ──────────────────────────────────────────────────────

func buildAromaticProbes(proteinAtoms []Atom, heavy []heavyAtom, ps *ProbeSet) {
	resAtoms := map[string][]Vec3{}
	resOrder := []string{}
	for _, a := range proteinAtoms {
		if isAromatic(a) {
			k := resKey(a)
			if _, ok := resAtoms[k]; !ok {
				resOrder = append(resOrder, k)
			}
			resAtoms[k] = append(resAtoms[k], a.Pos)
		}
	}
	fmt.Printf("Aromatic residues (PHE/TYR/TRP): %d\n", len(resOrder))

	buildHeavyExcluding := func(rk string) []heavyAtom {
		var out []heavyAtom
		for _, a := range proteinAtoms {
			if strings.ToUpper(a.Element) == "H" || isWater(a) {
				continue
			}
			if resKey(a) == rk {
				continue
			}
			out = append(out, heavyAtom{a.Pos, vdw(a.Element)})
		}
		return out
	}

	nParallel, nTshaped := 0, 0
	ringCentres := []Vec3{}

	for _, rk := range resOrder {
		pts := resAtoms[rk]
		if len(pts) < 3 {
			continue
		}
		centroid := Vec3{}
		for _, p := range pts {
			centroid = centroid.Add(p)
		}
		centroid = centroid.Scale(1.0 / float64(len(pts)))
		normal := ringPlaneNormal(pts)
		excl := buildHeavyExcluding(rk)

		tooClose := func(pos Vec3) bool {
			for _, c := range ringCentres {
				if pos.Sub(c).Norm() < 4.5 {
					return true
				}
			}
			return false
		}

		addRing := func(centre, norm Vec3) bool {
			ringPts, edges := benzeneRing(centre, norm)
			for _, p := range ringPts {
				if !noHardClash(p, vdw("C"), hardTol, excl) {
					return false
				}
			}
			base := ps.Len() + 1
			for _, p := range ringPts {
				ps.Add(p, "C")
			}
			for _, e := range edges {
				ps.Bond(base+e[0], base+e[1], 4)
			}
			// Methyl cap on ring atom 3 (para position), placed IN the ring
			// plane by using the radial direction from the probe ring centre.
			// This ensures the cap C is coplanar with the benzene ring.
			capDir := ringPts[3].Sub(centre).unit()
			capPos := ringPts[3].Add(capDir.Scale(1.54))
			capIdx := ps.Len() + 1
			ps.Add(capPos, "C")
			ps.Bond(base+3, capIdx, 1) // ring C – methyl cap (in-plane)
			return true
		}

		// Parallel stacking.
		for _, sign := range []float64{1, -1} {
			centre := centroid.Add(normal.Scale(sign * 3.5))
			if tooClose(centre) {
				continue
			}
			if addRing(centre, normal) {
				ringCentres = append(ringCentres, centre)
				nParallel++
				break
			}
		}

		// T-shaped stacking.
		for _, sign := range []float64{1, -1} {
			centre := centroid.Add(normal.Scale(sign * 5.0))
			if tooClose(centre) {
				continue
			}
			tNorm := perpendicular(normal)
			if addRing(centre, tNorm) {
				ringCentres = append(ringCentres, centre)
				nTshaped++
				break
			}
		}
	}
	fmt.Printf("Benzene rings: %d parallel, %d T-shaped\n", nParallel, nTshaped)
}

// ── H-bond vote probes ────────────────────────────────────────────────────────

func buildHbondProbes(proteinAtoms []Atom, heavy []heavyAtom, ps *ProbeSet) {
	var donors, acceptors []Atom
	for _, a := range proteinAtoms {
		role := hbondRole(a)
		if role == "donor" || role == "dual" {
			donors = append(donors, a)
		}
		if role == "acceptor" || role == "dual" {
			acceptors = append(acceptors, a)
		}
	}
	fmt.Printf("H-bond donors: %d  Acceptors: %d\n", len(donors), len(acceptors))

	placeHbondProbes := func(group []Atom, probeElem string) int {
		probeR := vdw(probeElem)
		votes := castVotes(group, hbondSphereR, hbondPairDist, nSamples, heavy, probeElem)
		type candidate struct {
			pos       Vec3
			pairCount int
		}
		var cands []candidate
		for _, v := range votes {
			if v.pairCount >= hbondMinPairs {
				cands = append(cands, candidate{v.pos, v.pairCount})
			}
		}
		sort.Slice(cands, func(i, j int) bool { return cands[i].pairCount > cands[j].pairCount })
		placed := 0
		var placed_pos []Vec3
		for _, c := range cands {
			if !noHardClash(c.pos, probeR, hardTol, heavy) {
				continue
			}
			tooClose := false
			for _, pp := range placed_pos {
				if c.pos.Sub(pp).Norm() < hbondProbeSpacing {
					tooClose = true
					break
				}
			}
			if tooClose {
				continue
			}
			placed_pos = append(placed_pos, c.pos)
			ps.Add(c.pos, probeElem)
			placed++
		}
		return placed
	}

	nAcc := placeHbondProbes(donors, "O")
	fmt.Printf("H-bond acceptor probes (O): %d\n", nAcc)
	nDon := placeHbondProbes(acceptors, "N")
	fmt.Printf("H-bond donor probes (N): %d\n", nDon)
}
