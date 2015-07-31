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

	"github.com/biogo/biogo/io/featio/gff"
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
		gf.SeqName = c.parts[0].genomic.chrom
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
