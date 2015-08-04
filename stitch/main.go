// Copyright ©2015 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// stitch is a program that joins repeat annotation features based on end points.
//
// The stitch program takes a collection of repeat features in GFF format with a
// feature attribute field "Repeat". The Repeat attribute holds two string fields,
// repeat type and repeat class, and three integer fields, start and end of alignment
// relative to the repeat consensus and the number of bases the consensus extends
// beyond the alignment end. For example:
//  Repeat AluJr SINE/Alu 3 295 17
// indicates the repeat is an AluJr, a SINE/Alu, and that the repeat annotation
// includes bases 3 to 295 of the AluJr consensus, stopping 17 bases before the end
// of the consensus.
//
// Stitch assembles repeat features by first grouping alignments by chromosome, strand
// and repeat class. It then performs a dynamic programming extension of a repeat
// feature chain by adding score improving segments while considering a cost function
// based on genomic and repeats consensus alignment end points. To extend a chain,
// the cost function must be be outweighed by the score gained by including the
// chain prefix.
//
package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"sort"

	"github.com/biogo/biogo/io/featio/gff"
)

var (
	inFile  = flag.String("in", "", "filename of a GFF file containing repeat annotations")
	workers = flag.Int("workers", 0, "number of parallel workers to use for stitching repeats (if 0 use GOMAXPROCS)")
)

func main() {
	flag.Parse()
	if *inFile == "" {
		flag.Usage()
		os.Exit(0)
	}
	if *workers == 0 {
		*workers = runtime.GOMAXPROCS(0)
	}

	f, err := os.Open(*inFile)
	if err != nil {
		log.Fatalf("could not open %q: %v", *inFile, err)
	}
	fmt.Fprintf(os.Stderr, "reading repeat features from %q\n", *inFile)
	defer f.Close()
	in := gff.NewReader(f)

	classes := make(map[partition][]*simple)
	for {
		f, err := in.Read()
		if err != nil {
			if err != io.EOF {
				log.Fatalf("failed to read source feature: %v", err)
			}
			break
		}

		gf := f.(*gff.Feature)
		r, err := newSimpleRepeat(gf)
		if err != nil {
			log.Fatal(err)
		}
		p := partition{
			chrom:  gf.SeqName,
			strand: gf.FeatStrand,
			class:  r.class,
		}
		classes[p] = append(classes[p], r)
	}

	// maxSeparation is the maximum distance between
	// successive element sorted end points that allow the
	// elements to be included in the same analysis block.
	const maxSeparation = 5e4

	var all []composite
	for c := range compositesFrom(classes, maxSeparation, *workers) {
		all = append(all, c...)
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
		right := c.parts[0].genomic.right
		for _, p := range c.parts[1:] {
			if p.right > right {
				right = p.right
			}
		}
		gf.SeqName = c.parts[0].genomic.chrom
		gf.FeatStart = c.parts[0].genomic.left
		gf.FeatEnd = right
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
