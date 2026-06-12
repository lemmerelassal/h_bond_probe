package main

import (
	"bufio"
	"fmt"
	"io"
	"net/http"
	"os"
	"strconv"
	"strings"
)

// LoadAtoms reads atoms from a local PDB/mmCIF file, downloads from RCSB,
// or downloads from AlphaFold DB (when alphafold=true).
// Returns the atom slice, the uppercase ID (for output naming), and any error.
func LoadAtoms(arg string, alphafold bool) ([]Atom, string, error) {
	if _, err := os.Stat(arg); err == nil {
		data, e := os.ReadFile(arg)
		if e != nil {
			return nil, "", e
		}
		var atoms []Atom
		var parseErr error
		if strings.HasSuffix(strings.ToLower(arg), ".cif") {
			atoms, parseErr = parseCIF(string(data))
		} else {
			atoms, parseErr = parsePDB(string(data))
		}
		if parseErr != nil {
			return nil, "", parseErr
		}
		base := arg
		for _, suf := range []string{".pdb", ".cif", ".ent"} {
			base = strings.TrimSuffix(base, suf)
		}
		if i := strings.LastIndexAny(base, "/\\"); i >= 0 {
			base = base[i+1:]
		}
		fmt.Printf("Loaded %d atoms from %s\n", len(atoms), arg)
		return atoms, strings.ToUpper(base), nil
	}

	id := strings.ToUpper(strings.TrimSpace(arg))

	if alphafold {
		// AlphaFold DB: try recent model versions, because some entries only
		// expose newer versions such as model_v6.
		for _, version := range []int{6, 5, 4, 3} {
			url := fmt.Sprintf("https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v%d.cif", id, version)
			fmt.Printf("Fetching AlphaFold structure: %s\n", url)
			resp, err := http.Get(url) //nolint:noctx
			if err != nil {
				if resp != nil {
					resp.Body.Close()
				}
				continue
			}
			if resp.StatusCode != 200 {
				resp.Body.Close()
				continue
			}
			data, _ := io.ReadAll(resp.Body)
			resp.Body.Close()
			atoms, _ := parseCIF(string(data))
			if len(atoms) > 0 {
				fmt.Printf("Parsed %d atoms from AlphaFold.\n", len(atoms))
				return atoms, "AF_" + id, nil
			}
		}
		return nil, "", fmt.Errorf("could not fetch AlphaFold structure for: %s", id)
	}

	for _, ext := range []string{"pdb", "cif"} {
		url := fmt.Sprintf("https://files.rcsb.org/download/%s.%s", id, ext)
		fmt.Printf("Trying %s\n", url)
		resp, err := http.Get(url) //nolint:noctx
		if err != nil || resp.StatusCode != 200 {
			if resp != nil {
				resp.Body.Close()
			}
			continue
		}
		data, _ := io.ReadAll(resp.Body)
		resp.Body.Close()
		var atoms []Atom
		if ext == "cif" {
			atoms, _ = parseCIF(string(data))
		} else {
			atoms, _ = parsePDB(string(data))
		}
		if len(atoms) > 0 {
			fmt.Printf("Parsed %d atoms.\n", len(atoms))
			return atoms, id, nil
		}
	}
	return nil, "", fmt.Errorf("could not load structure: %s", arg)
}

func parsePDB(text string) ([]Atom, error) {
	var atoms []Atom
	sc := bufio.NewScanner(strings.NewReader(text))
	for sc.Scan() {
		line := sc.Text()
		if len(line) < 54 {
			continue
		}
		rec := strings.TrimSpace(line[:6])
		if rec != "ATOM" && rec != "HETATM" {
			continue
		}
		x, _ := strconv.ParseFloat(strings.TrimSpace(line[30:38]), 64)
		y, _ := strconv.ParseFloat(strings.TrimSpace(line[38:46]), 64)
		z, _ := strconv.ParseFloat(strings.TrimSpace(line[46:54]), 64)
		elem := ""
		if len(line) >= 78 {
			elem = strings.TrimSpace(line[76:78])
		}
		name := strings.TrimSpace(line[12:16])
		if elem == "" {
			elem = guessElem(name)
		}
		resSeq, _ := strconv.Atoi(strings.TrimSpace(line[22:26]))
		atoms = append(atoms, Atom{
			Name:    name,
			ResName: strings.TrimSpace(line[17:20]),
			ChainID: string(line[21]),
			ResSeq:  resSeq,
			Pos:     Vec3{x, y, z},
			Element: elem,
			IsHet:   rec == "HETATM",
		})
	}
	return atoms, nil
}

func parseCIF(text string) ([]Atom, error) {
	var atoms []Atom
	lines := strings.Split(text, "\n")
	inLoop := false
	colIndex := map[string]int{}
	colCount := 0
	dataStart := -1

	for i, rawLine := range lines {
		line := strings.TrimSpace(rawLine)
		if line == "loop_" {
			if i+1 < len(lines) && strings.HasPrefix(strings.TrimSpace(lines[i+1]), "_atom_site.") {
				inLoop = true
				colIndex = map[string]int{}
				colCount = 0
			}
			continue
		}
		if !inLoop {
			continue
		}
		if strings.HasPrefix(line, "_atom_site.") {
			colIndex[strings.TrimPrefix(line, "_atom_site.")] = colCount
			colCount++
			dataStart = i + 1
			continue
		}
		if strings.HasPrefix(line, "_") || line == "" {
			if dataStart > 0 {
				inLoop = false
			}
			continue
		}
		if dataStart < 0 {
			continue
		}
		fields := splitCIF(line)
		if len(fields) < colCount {
			continue
		}
		get := func(tag string) string {
			if idx, ok := colIndex[tag]; ok && idx < len(fields) {
				return fields[idx]
			}
			return ""
		}
		rec := get("group_PDB")
		if rec != "ATOM" && rec != "HETATM" {
			continue
		}
		x, _ := strconv.ParseFloat(get("Cartn_x"), 64)
		y, _ := strconv.ParseFloat(get("Cartn_y"), 64)
		z, _ := strconv.ParseFloat(get("Cartn_z"), 64)
		resSeq, _ := strconv.Atoi(get("label_seq_id"))
		atoms = append(atoms, Atom{
			Name:    get("label_atom_id"),
			ResName: get("label_comp_id"),
			ChainID: get("label_asym_id"),
			ResSeq:  resSeq,
			Pos:     Vec3{x, y, z},
			Element: strings.TrimSpace(get("type_symbol")),
			IsHet:   rec == "HETATM",
		})
	}
	return atoms, nil
}

func splitCIF(line string) []string {
	var fields []string
	inQ := false
	var cur strings.Builder
	for _, c := range line {
		switch {
		case c == '\'' && !inQ:
			inQ = true
		case c == '\'' && inQ:
			inQ = false
		case c == ' ' && !inQ:
			if cur.Len() > 0 {
				fields = append(fields, cur.String())
				cur.Reset()
			}
		default:
			cur.WriteRune(c)
		}
	}
	if cur.Len() > 0 {
		fields = append(fields, cur.String())
	}
	return fields
}

func guessElem(name string) string {
	n := strings.TrimSpace(name)
	for len(n) > 0 && n[0] >= '0' && n[0] <= '9' {
		n = n[1:]
	}
	if len(n) == 0 {
		return "C"
	}
	if len(n) >= 2 {
		switch strings.ToUpper(n[:2]) {
		case "CL", "BR", "FE", "ZN", "MG", "CA", "MN", "CU", "CO", "NI":
			return strings.ToUpper(n[:2])
		}
	}
	return strings.ToUpper(string(n[0]))
}
