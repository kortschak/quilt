// Copyright ©2013 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// brahma performs annotation of GFF intervals produced by PALS/PILER, taking
// annotation information from a GFF file generated from RepeatMasker output.
package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/biogo/biogo/seq"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
)

var (
	inFile = flag.String("in", "", "Filename for source annotation.")
	help   = flag.Bool("help", false, "Print this usage message.")
)

func main() {
	flag.Parse()
	if *help || *inFile == "" {
		flag.Usage()
		os.Exit(0)
	}

	f, err := os.Open(*inFile)
	if err != nil {
		log.Fatalf("could not open %q: %v", *inFile, err)
	}
	fmt.Fprintf(os.Stderr, "reading repeat features from %q.\n", *inFile)
	defer f.Close()
	source := gff.NewReader(f)

	classes := make(map[partition][]*record)
	for {
		f, err := source.Read()
		if err != nil {
			if err != io.EOF {
				log.Fatalf("failed to read source feature: %v", err)
			}
			break
		}

		gf := f.(*gff.Feature)
		repData := &record{
			genomic: repeat{
				left:   gf.FeatStart,
				right:  gf.FeatEnd,
				loc:    contig(gf.SeqName),
				strand: gf.FeatStrand,
			},
		}
		if gf.FeatScore != nil {
			repData.score = *gf.FeatScore
		}

		ra := gf.FeatAttributes.Get("Repeat")
		if ra == "" {
			log.Fatal("missing repeat tag: file probably not an RM gff.")
		}
		err = repData.parse(ra)
		if err != nil {
			log.Fatalf("failed to parse repeat tag: %v\n", err)
		}

		p := partition{
			chr:    gf.SeqName,
			strand: gf.FeatStrand,
			class:  repData.class,
		}
		classes[p] = append(classes[p], repData)
	}
	for p, c := range classes {
		fmt.Printf("%+v %d\n", p, len(c))
		if len(c) < 2 || c[0].left == none {
			fmt.Println("\tskip")
			continue
		}
		sort.Sort(records(c))
		var splits int
		for i, r := range c[1:] {
			if r.genomic.right-c[i].genomic.right > 5e4 {
				splits++
			}
		}
		fmt.Println("\tpotential splits:", splits)

		a := make([][]element, len(c))
		for i := 0; i < 2; i++ {
			a[i] = make([]element, len(c))
		}
		for i, r := range c {
			a[0][i] = element{link: i, score: r.score}
		}
		a[1][0].score = a[0][0].score

		fn := func(left, right *record) (score float64, ok bool) {
			ok = right.genomic.right-left.genomic.right <= 5e4
			if !ok || right.genomic.strand == seq.None {
				return math.Inf(-1), ok
			}

			score = left.score

			var overlap int

			// Genomic coordinate term.
			overlap = left.genomic.right - right.genomic.left
			if overlap > 0 {
				// Case where there is a back step.
				score -= float64(overlap) * 10
			} else {
				// Case where there is a separation.
				score -= float64(overlap) * -1
			}

			// Repeat coordinate term.
			if right.genomic.strand == seq.Plus {
				overlap = left.right - right.left
			} else {
				overlap = right.right - left.left
			}
			if overlap > 0 {
				// Case where there is a back step.
				score -= float64(overlap) * 100
			} else {
				// Case where there is a separation.
				score -= float64(overlap) * -1
			}

			return score, true
		}

		for i := 1; i < len(a); i++ {
			for j := i; j < len(c); j++ {
				a[i][j] = max(c, a, i, j, fn)
			}
			if i < len(c)-1 {
				a[i+1], a[i-1] = a[i-1], nil
			}
		}

		final := a[len(a)-1]
		wasUsed := make([]bool, len(final))
		for i := len(final) - 1; i >= 0; i-- {
			if wasUsed[i] {
				continue
			}
			var (
				parts []part
				p     int
			)
			for p = i; p != final[p].link; p = final[p].link {
				wasUsed[p] = true
				parts = append(parts, part{name: c[p].name, left: c[p].left, right: c[p].right, genomic: c[p].genomic})

			}
			if p == i {
				continue
			}
			parts = append(parts, part{name: c[p].name, left: c[p].left, right: c[p].right, genomic: c[p].genomic})

			reverse(parts)
			fmt.Printf("\t%v\n", parts)
		}
	}
}

type element struct {
	link  int
	score float64
}

func max(c []*record, a [][]element, i, j int, fn func(left, right *record) (score float64, ok bool)) element {
	e := a[i-1][j]
	right := c[j]
	for k := j - 1; j >= 0; j-- {
		left := c[k]
		s, ok := fn(left, right)
		if !ok {
			break
		}
		s += a[i-1][k].score
		if s > e.score {
			e.score = s
			e.link = k
		}
	}
	return e
}

type partition struct {
	chr    string
	strand seq.Strand
	class  string
}

type part struct {
	name        string
	left, right int
	genomic     repeat
}

func reverse(s []part) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

// contig is a sequence contig with repeats mapped to it.
type contig string

func (c contig) Start() int             { return 0 }
func (c contig) End() int               { return 0 }
func (c contig) Len() int               { return 0 }
func (c contig) Name() string           { return string(c) }
func (c contig) Description() string    { return "contig" }
func (c contig) Location() feat.Feature { return nil }

// repeat is a repeat-matching interval.
type repeat struct {
	left, right int

	loc    feat.Feature
	strand seq.Strand
}

func (r repeat) Start() int             { return r.left }
func (r repeat) End() int               { return r.right }
func (r repeat) Len() int               { return r.right - r.left }
func (r repeat) Name() string           { return fmt.Sprintf("%s:[%d,%d)", r.loc.Name(), r.left, r.right) }
func (r repeat) Description() string    { return "repeat" }
func (r repeat) Location() feat.Feature { return r.loc }

// record is a masked repeat record.
type record struct {
	// genomic is the genomic region matched
	// to the the repeat identified below.
	genomic repeat

	// name and class that the repeat type
	// and class defined by the masker.
	name, class string

	// score is the feature score.
	score float64

	// left and right are the left and right
	// position of the record alignment in
	// consensus-relative coordinates.
	left, right int
}

const none = -1

func (r *record) parse(a string) error {
	fields := strings.Split(a, " ")

	r.name = fields[0]
	r.class = fields[1]
	var err error
	if fields[2] != "." {
		r.left, err = strconv.Atoi(fields[2])
		if err != nil {
			return err
		}
		r.left-- // Convert to 0-based indexing.
	} else {
		r.left = none
	}
	if fields[3] != "." {
		r.right, err = strconv.Atoi(fields[3])
		if err != nil {
			return err
		}
	} else {
		r.right = none
	}

	return nil
}

type records []*record

func (r records) Len() int { return len(r) }
func (r records) Less(i, j int) bool {
	iName := r[i].genomic.loc.Name()
	jName := r[j].genomic.loc.Name()
	return iName < jName || (iName == jName && r[i].genomic.right < r[j].genomic.right)
}
func (r records) Swap(i, j int) { r[i], r[j] = r[j], r[i] }
