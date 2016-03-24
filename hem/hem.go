// Copyright ©2015 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// hem performs a basic analysis of the result of running stitch on a set
// of repeat annotations.
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
	"text/tabwriter"

	"github.com/biogo/biogo/io/featio/gff"
)

var (
	inFile   = flag.String("in", "", "Filename for stitch checking.")
	discords = flag.Bool("discords", false, "Output discordant features to stderr.")
	help     = flag.Bool("help", false, "Print this usage message.")
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
	in := gff.NewReader(f)

	var w *gff.Writer
	if *discords {
		w = gff.NewWriter(os.Stderr, 60, false)
	}

	var (
		n      int
		merged int

		partsIn   = make(map[string]int)
		discordIn = make(map[string]int)
		allIn     = make(map[string]int)
	)
	for {
		f, err := in.Read()
		if err != nil {
			if err != io.EOF {
				log.Fatalf("failed to read source feature: %v", err)
			}
			break
		}

		gf := f.(*gff.Feature)
		if gf.Source != "stitch" {
			continue
		}
		n++

		class := gf.FeatAttributes.Get("Class")
		if class == "" {
			log.Fatal("missing class tag: file not a valid stitch gff.")
		}
		class, err = strconv.Unquote(class)
		if err != nil {
			log.Fatalf("failed to unquote class: %v", err)
		}
		p := gf.FeatAttributes.Get("Parts")
		if p == "" {
			log.Fatal("missing parts tag: file not a valid stitch gff.")
		}
		p, err = strconv.Unquote(p)
		if err != nil {
			log.Fatalf("failed to unquote parts: %v", err)
		}
		allIn[class]++

		parts := strings.Split(p, "|")
		partsIn[class] += len(parts)
		merged += len(parts)
		first := strings.Fields(parts[0])[0]
		for _, p := range parts[1:] {
			if strings.Fields(p)[0] != first {
				discordIn[class]++
				if w != nil {
					w.Write(gf)
				}
				break
			}
		}
	}

	names := make([]string, 0, len(allIn))
	for c := range allIn {
		names = append(names, c)
	}
	sort.Strings(names)
	tw := tabwriter.NewWriter(os.Stdout, 0, 0, 1, ' ', tabwriter.Debug)
	fmt.Fprintln(tw, "class\t parts merged\t composites\t condensation\t discord count\t discord freq")
	for _, c := range names {
		fmt.Fprintf(tw, "%s\t%12d\t%10d\t    % .2f\t%13d\t   % .3f\n",
			c, partsIn[c], allIn[c], float64(partsIn[c])/float64(allIn[c]),
			discordIn[c], float64(discordIn[c])/float64(allIn[c]),
		)
	}
	tw.Flush()
	fmt.Printf("\ntotal composited: %d comprising: %d\n", n, merged)
}
