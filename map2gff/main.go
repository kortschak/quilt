// Copyright ©2015 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// rm2gff converts Censor map files to GFF including the stitch-required Repeat attribute.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq"
	"github.com/biogo/biogo/seq/linear"
)

const (
	queryNameField = iota
	queryStartField
	queryEndField
	repeatTypeField
	repeatStartField
	repeatEndField
	strandField
	_ // alignment similarity
	_ // alignment positive fraction
	scoreField
	_ // query coverage fraction - not used because we need class name anyway.
	_ // repeat coverage fraction - not used because we need class name anyway.

	numberOfFields
)

var (
	libFile   = flag.String("lib", "", "fasta file to use to define repeat family lengths and classes")
	defFile   = flag.String("defs", "", "tab delimted file to use to define repeat family lengths and classes")
	defHeader = flag.Bool("defs-header", true, "defs file has header")
)

func main() {
	flag.Parse()
	if *libFile == "" && *defFile == "" {
		flag.Usage()
		os.Exit(1)
	}

	var (
		classes map[string]record
		err     error
	)
	switch {
	case *libFile != "":
		classes, err = readClassesFromFasta(*libFile)
	case *defFile != "":
		classes, err = readClassesFromDefs(*defFile, *defHeader)
	}
	if err != nil {
		log.Fatal(err)
	}

	f := &gff.Feature{
		Source:         "Censor",
		Feature:        "repeat",
		FeatFrame:      gff.NoFrame,
		FeatAttributes: gff.Attributes{{Tag: "Repeat"}},
	}
	w := gff.NewWriter(os.Stdout, 60, true)
	sc := bufio.NewScanner(os.Stdin)
	for n := 1; sc.Scan(); n++ {
		err := fill(f, strings.Fields(sc.Text()), classes)
		if err != nil {
			log.Fatalf("parse error on line %d: %v", n, err)
		}
		w.Write(f)
	}
	err = sc.Err()
	if err != nil {
		log.Fatal(err)
	}
}

type record struct {
	name   string
	length int
}

func readClassesFromFasta(file string) (map[string]record, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, fmt.Errorf("failed to open library %q: %v", file, err)
	}
	defer f.Close()

	table := make(map[string]record)
	r := fasta.NewReader(f, linear.NewSeq("", nil, alphabet.DNA))
	sc := seqio.NewScanner(r)
	for sc.Next() {
		var class string
		desc := sc.Seq().Description()
		if len(desc) == 0 {
			class = sc.Seq().Name()
		} else {
			class = strings.Replace(strings.Split(desc, "\t")[0], " ", "_", -1)
		}
		table[sc.Seq().Name()] = record{
			name:   class,
			length: sc.Seq().Len(),
		}
	}
	if sc.Error() != nil {
		return nil, fmt.Errorf("failed during library read: %v", sc.Error())
	}

	return table, nil
}

func readClassesFromDefs(file string, header bool) (map[string]record, error) {
	f, err := os.Open(file)
	if err != nil {
		return nil, fmt.Errorf("failed to open library %q: %v", file, err)
	}
	defer f.Close()

	table := make(map[string]record)
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		if header {
			header = false
			continue
		}
		f := strings.Split(sc.Text(), "\t")
		if len(f) != 3 {
			log.Fatalf("line does not have 3 fields: %q", sc.Text())
		}
		length, err := strconv.Atoi(f[2])
		if err != nil {
			log.Fatalf("failed to parse length for %q: %v", sc.Text(), err)
		}

		table[f[0]] = record{
			name:   strings.Replace(f[1], " ", "_", -1),
			length: length,
		}
	}

	return table, nil
}

func handlePanic(err *error) {
	r := recover()
	if r != nil {
		switch r := r.(type) {
		case error:
			*err = r
		default:
			panic(r)
		}
	}
}

func mustAtoi(s string) int {
	i, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return i
}

func mustAtofp(s string) *float64 {
	f, err := strconv.ParseFloat(s, 64)
	if err != nil {
		panic(err)
	}
	return &f
}

func mustOutToSane(s string) seq.Strand {
	switch s {
	case "c":
		return seq.Minus
	case "d":
		return seq.Plus
	default:
		panic(fmt.Errorf("illegal strand: %q", s))
	}
}

func fill(f *gff.Feature, data []string, classes map[string]record) (err error) {
	defer handlePanic(&err)
	f.FeatScore = mustAtofp(data[scoreField])
	f.SeqName = data[queryNameField]
	f.FeatStart = feat.OneToZero(mustAtoi(data[queryStartField]))
	f.FeatEnd = mustAtoi(data[queryEndField])
	f.FeatStrand = mustOutToSane(data[strandField])
	f.FeatAttributes[0].Value = repeatAttribute(data, classes)
	return
}

func repeatAttribute(data []string, classes map[string]record) string {
	name := data[repeatTypeField]
	left := mustAtoi(data[repeatStartField])
	right := mustAtoi(data[repeatEndField])
	class, ok := classes[name]
	if !ok {
		panic(fmt.Errorf("no record for %q", name))
	}
	remains := class.length - right
	if remains < 0 {
		panic(fmt.Errorf("remaining sequence less than zero: %d", remains))
	}
	return fmt.Sprintf("%s %s %d %d %d", name, class.name, left, right, remains)
}
