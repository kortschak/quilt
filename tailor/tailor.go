// Copyright ©2015 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"

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
	fmt.Fprintf(os.Stderr, "reading stitched repeat features from %q.\n", *inFile)
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
			log.Fatal("missing parts tag: file not a valid stitch gff.")
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
			trees[c.SeqName] = t
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

	w := gff.NewWriter(os.Stdout, 60, true)
	for _, q := range flag.Args() {
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

			if !hitsComposite(f.(*gff.Feature), trees) {
				w.Write(f)
			}
		}
		f.Close()
	}
}

func hitsComposite(f *gff.Feature, trees map[string]*interval.IntTree) bool {
	t, ok := trees[f.SeqName]
	if !ok {
		return false
	}
	for _, m := range t.Get((*query)(f)) {
		in := m.(composite)
		for _, p := range in.parts {
			if f.FeatStart < p.right && f.FeatEnd > p.left {
				return true
			}
		}
	}
	return false
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
