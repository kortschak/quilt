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

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/store/interval"
)

var (
	inFile = flag.String("in", "", "Filename for stitch TinT analysis.")
	help   = flag.Bool("help", false, "Print this usage message.")
)

func main() {
	flag.Parse()
	if *help || *inFile == "" || len(flag.Args()) == 0 {
		flag.Usage()
		os.Exit(0)
	}

	f, err := os.Open(*inFile)
	if err != nil {
		log.Fatalf("could not open %q: %v", *inFile, err)
	}
	fmt.Fprintf(os.Stderr, "reading repeat features from %q.\n", *inFile)
	in := gff.NewReader(f)

	trees := make(map[string]*interval.IntTree)
	for {
		f, err := in.Read()
		if err != nil {
			if err != io.EOF {
				log.Fatalf("failed to read source feature: %v", err)
			}
			break
		}

		c := composite{Feature: f.(*gff.Feature)}
		if c.Source != "stitch" {
			continue
		}

		parts := c.FeatAttributes.Get("Parts")
		if parts == "" {
			log.Fatal("missing parts tag: file not a stitch gff.")
		}
		parts, err = strconv.Unquote(parts)
		if err != nil {
			log.Fatalf("failed to unquote parts: %v", err)
		}

		for _, p := range strings.Split(parts, "|") {
			fields := strings.Fields(p)
			if len(fields) != 5 {
				log.Fatalf("unexpected number of fields in %q: %v", p, c.Feature)
			}
			left, err := strconv.Atoi(fields[3])
			if err != nil {
				log.Fatalf("failed to parse left coordinate: %v", fields[3])
			}
			right, err := strconv.Atoi(fields[3])
			if err != nil {
				log.Fatalf("failed to parse right coordinate: %v", fields[4])
			}
			c.parts = append(c.parts, part{left: left, right: right})
		}

		t, ok := trees[c.SeqName]
		if !ok {
			t = &interval.IntTree{}
			trees[c.Feature.SeqName] = t
		}
		c.id = uintptr(t.Len())
		t.Insert(c, true)
	}
	f.Close()

	var chroms []string
	for chr, t := range trees {
		t.AdjustRanges()
		chroms = append(chroms, chr)
	}
	sort.Strings(chroms)

	haveInsertion := make(map[*gff.Feature][]*gff.Feature)
	for _, q := range flag.Args() {
		if q == *inFile {
			// We already have this in memory, so don't read it again.
			for _, t := range trees {
				t.Do(func(i interval.IntInterface) (done bool) {
					noteComposite(i.(composite).Feature, trees, haveInsertion)
					return
				})
			}
			continue
		}

		f, err := os.Open(q)
		if err != nil {
			log.Fatalf("could not open %q: %v", q, err)
		}
		fmt.Fprintf(os.Stderr, "reading repeat features from %q.\n", q)
		in := gff.NewReader(f)

		for {
			f, err := in.Read()
			if err != nil {
				if err != io.EOF {
					log.Fatalf("failed to read source feature: %v", err)
				}
				break
			}

			noteComposite(f.(*gff.Feature), trees, haveInsertion)
		}
		f.Close()
	}

	w := gff.NewWriter(os.Stdout, 60, true)
	for _, chr := range chroms {
		trees[chr].Do(func(i interval.IntInterface) (done bool) {
			in := i.(composite).Feature
			if f, ok := haveInsertion[in]; ok {
				sort.Sort(byGenomeLocation(f))
				in.Source = "patch"
				in.FeatAttributes = append(in.FeatAttributes, gff.Attribute{
					Tag:   "TinT",
					Value: formatInsertions(f),
				})
			}
			w.Write(in)
			return
		})
	}
}

func noteComposite(f *gff.Feature, trees map[string]*interval.IntTree, hasInsertion map[*gff.Feature][]*gff.Feature) {
	t, ok := trees[f.SeqName]
	if !ok {
		return
	}
	for _, m := range t.Get((*query)(f)) {
		inGap := true
		in := m.(composite)
		for _, p := range in.parts {
			if f.FeatStart < p.right && f.FeatEnd > p.left {
				inGap = false
				break
			}
		}
		if inGap {
			hasInsertion[in.Feature] = append(hasInsertion[in.Feature], f)
		}
	}
}

func formatInsertions(p []*gff.Feature) string {
	var buf bytes.Buffer
	fmt.Fprintf(&buf, "%d ", len(p))
	for i, f := range p {
		if i == 0 {
			buf.WriteByte('"')
		} else {
			buf.WriteByte('|')
		}
		name := f.FeatAttributes.Get("Repeat")
		if name != "" {
			np := strings.Fields(name)[:2]
			name = fmt.Sprintf("%s/%s", np[1], np[0])
		} else if name = f.FeatAttributes.Get("Class"); name != "" {
			name, _ = strconv.Unquote(name)
		}
		fmt.Fprintf(&buf, `%s %d %d %s`,
			name, feat.ZeroToOne(f.FeatStart), f.FeatEnd, f.FeatStrand,
		)
	}
	buf.WriteByte('"')
	return buf.String()
}

type composite struct {
	*gff.Feature
	id uintptr

	parts []part
}

type part struct {
	left, right int
}

// Overlap returns whether q overlaps b.
func (c composite) Overlap(b interval.IntRange) bool {
	return c.FeatEnd > b.Start && c.FeatStart < b.End
}
func (c composite) ID() uintptr              { return c.id }
func (c composite) Range() interval.IntRange { return interval.IntRange{c.FeatStart, c.FeatEnd} }

type query gff.Feature

// Overlap returns whether q is entirely within b.
func (q *query) Overlap(b interval.IntRange) bool {
	return q.FeatStart >= b.Start && q.FeatEnd <= b.End
}
func (q *query) ID() uintptr              { return 0 }
func (q *query) Range() interval.IntRange { return interval.IntRange{q.FeatStart, q.FeatEnd} }

type byGenomeLocation []*gff.Feature

func (f byGenomeLocation) Len() int { return len(f) }
func (f byGenomeLocation) Less(i, j int) bool {
	iName := f[i].SeqName
	jName := f[j].SeqName
	return iName < jName ||
		(iName == jName && f[i].FeatStart < f[j].FeatStart) ||
		(iName == jName && f[i].FeatStart == f[j].FeatStart && f[i].FeatEnd > f[j].FeatEnd)
}
func (f byGenomeLocation) Swap(i, j int) { f[i], f[j] = f[j], f[i] }
