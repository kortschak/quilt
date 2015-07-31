// Copyright ©2015 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
)

var (
	inFile  = flag.String("in", "", "Filename for source annotation.")
	workers = flag.Int("workers", 1, "Number of parallel workers.")
	help    = flag.Bool("help", false, "Print this usage message.")
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

	// maxSeparation is the maximum distance between
	// successive element sorted end points that allow the
	// elements to be included in the same analysis block.
	const maxSeparation = 5e4

	done := make(chan []composite)
	limit := make(chan struct{}, *workers)
	var wg sync.WaitGroup
	go func() {
		for p, c := range classes {
			p, c := p, c
			wg.Add(1)
			limit <- struct{}{}
			go func() {
				defer func() {
					<-limit
					wg.Done()
				}()

				if len(c) < 2 || c[0].left == none {
					log.Printf("%v records=%d - skip", p, len(c))
					return
				}

				sort.Sort(records(c))

				var splits int
				for i, r := range c[1:] {
					if r.genomic.right-c[i].genomic.right > maxSeparation {
						splits++
					}
				}
				log.Printf("%v records=%d splits=%d", p, len(c), splits)

				i := 0
				for j, r := range c[1:] {
					if r.genomic.right-c[j].genomic.right > maxSeparation || j == len(c)-2 {
						n := (j + 2) - i
						if n < 2 {
							continue
						}
						if *workers == 1 {
							fmt.Fprintf(os.Stderr, "split size: %d\n", n)
						}
						part := stitch(c[i:j+2], costFunction)
						i = j + 2

						if *workers == 1 {
							sort.Sort(byGenomeLocation(part))
							for _, f := range part {
								fmt.Fprintf(os.Stderr, "\t%+v\n", f)
							}
						}

						done <- part
					}
				}
			}()
		}
		wg.Wait()
		close(done)
	}()

	var all []composite
	for composites := range done {
		all = append(all, composites...)
	}
	log.Println("chaining complete.")

	sort.Sort(byGenomeLocation(all))
	w := gff.NewWriter(os.Stdout, 60, true)
	gf := &gff.Feature{
		Source:    "stitch",
		Feature:   "composite",
		FeatFrame: gff.NoFrame,
	}
	for _, c := range all {
		gf.SeqName = c.parts[0].genomic.loc.Name()
		gf.FeatStart = c.parts[0].genomic.left
		gf.FeatEnd = c.parts[len(c.parts)-1].genomic.right
		score := c.score
		gf.FeatScore = &score
		gf.FeatStrand = c.parts[0].genomic.strand

		gf.FeatAttributes = gff.Attributes{
			{Tag: "Class", Value: `"` + c.class + `"`},
			{Tag: "Parts", Value: c.parts.String()},
		}

		w.Write(gf)
	}
}

func abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

type composite struct {
	class string
	score float64
	parts parts
}

type parts []part

func (p parts) String() string {
	var buf bytes.Buffer
	for i, e := range p {
		if i == 0 {
			buf.WriteByte('"')
		} else {
			buf.WriteByte('|')
		}
		fmt.Fprintf(&buf, `%s %d %d %d %d`,
			e.name,
			feat.ZeroToOne(e.left), e.right,
			feat.ZeroToOne(e.genomic.left), e.genomic.right,
		)
	}
	buf.WriteByte('"')
	return buf.String()
}

type byGenomeLocation []composite

func (c byGenomeLocation) Len() int { return len(c) }
func (c byGenomeLocation) Less(i, j int) bool {
	iName := c[i].parts[0].genomic.loc.Name()
	jName := c[j].parts[0].genomic.loc.Name()
	return iName < jName || (iName == jName && c[i].parts[0].genomic.left < c[j].parts[0].genomic.left)
}
func (c byGenomeLocation) Swap(i, j int) { c[i], c[j] = c[j], c[i] }

func stitch(c []*record, fn func(left, right *record) (score float64, ok bool)) []composite {
	if len(c) < 2 {
		return nil
	}

	// Dynamic programming to maximise the score of chained features.
	a := make([][]element, len(c))
	for i := 0; i < 2; i++ {
		a[i] = make([]element, len(c))
	}
	for i, r := range c {
		a[0][i] = element{link: i, score: r.score}
	}
	a[1][0].score = a[0][0].score

	for i := 1; i < len(a); i++ {
		for j := i; j < len(c); j++ {
			a[i][j] = max(c, a, i, j, fn)
		}
		if i < len(c)-1 {
			a[i+1], a[i-1] = a[i-1], nil
		}
	}

	var composites []composite

	final := a[len(a)-1]
	wasUsed := make([]bool, len(final))
	for i := len(final) - 1; i >= 0; i-- {
		if wasUsed[i] {
			continue
		}
		var (
			cmp composite
			p   int
		)
		cmp.score = final[i].score
		cmp.class = c[0].class
		for p = i; p != final[p].link; p = final[p].link {
			wasUsed[p] = true
			cmp.parts = append(cmp.parts, part{name: c[p].name, left: c[p].left, right: c[p].right, genomic: c[p].genomic})

		}
		if p == i {
			continue
		}
		cmp.parts = append(cmp.parts, part{name: c[p].name, left: c[p].left, right: c[p].right, genomic: c[p].genomic})

		reverse(cmp.parts)
		composites = append(composites, cmp)
	}

	return composites
}

type element struct {
	link  int
	score float64
}

func max(c []*record, a [][]element, i, j int, fn func(left, right *record) (score float64, ok bool)) element {
	e := a[i-1][j]
	right := c[j]
	for k := j - 1; k >= 0; k-- {
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

func (p partition) String() string {
	return fmt.Sprintf("chr:%s strand:(%v) class:%s", p.chr, p.strand, p.class)
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
	loc feat.Feature

	left, right int

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
