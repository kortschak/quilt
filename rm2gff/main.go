// Copyright ©2015 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// rm2gff converts RM out files to GFF including the stitch-required Repeat attribute.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/biogo/biogo/seq"
)

const firstDataLine = 4

const (
	swScoreField = iota
	fracDivergeField
	fracDelField
	fracInsField
	queryNameField
	queryStartField
	queryEndField
	queryRemainingField
	strandField
	repeatTypeField
	repeatClassField
	repeatPosField1
	repeatPosField2
	repeatPosField3
	idField
	otherMatchField

	numberOfFields
)

var otherMatchAttribute = flag.Bool("mark-other", false, "mark features where RM indicates another higher score match overlaps")

func main() {
	flag.Parse()
	f := &gff.Feature{
		Source:         "RepeatMasker",
		Feature:        "repeat",
		FeatFrame:      gff.NoFrame,
		FeatAttributes: gff.Attributes{{Tag: "Repeat"}},
	}
	if *otherMatchAttribute {
		f.FeatAttributes = append(f.FeatAttributes, gff.Attribute{Tag: "OtherMatch", Value: "yes"})
	}
	w := gff.NewWriter(os.Stdout, 60, true)
	sc := bufio.NewScanner(os.Stdin)
	for n := 1; sc.Scan(); n++ {
		if n < firstDataLine {
			continue
		}
		err := fill(f, strings.Fields(sc.Text()), *otherMatchAttribute)
		if err != nil {
			log.Fatalf("parse error on line %d: %v", n, err)
		}
		w.Write(f)
	}
	err := sc.Err()
	if err != nil {
		log.Fatal(err)
	}
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

func mustRMtoSane(s string) seq.Strand {
	switch s {
	case "C":
		return seq.Minus
	case "+":
		return seq.Plus
	default:
		panic(fmt.Errorf("illegal strand: %q", s))
	}
}

func fill(f *gff.Feature, data []string, markOther bool) (err error) {
	defer handlePanic(&err)
	f.FeatScore = mustAtofp(data[swScoreField])
	f.SeqName = data[queryNameField]
	f.FeatStart = feat.OneToZero(mustAtoi(data[queryStartField]))
	f.FeatEnd = mustAtoi(data[queryEndField])
	f.FeatStrand = mustRMtoSane(data[strandField])
	f.FeatAttributes[0].Value = repeatAttribute(data)
	if markOther && len(data) == numberOfFields && data[otherMatchField] == "*" {
		f.FeatAttributes = f.FeatAttributes[:2]
	} else {
		f.FeatAttributes = f.FeatAttributes[:1]
	}
	return
}

func repeatAttribute(data []string) string {
	var left, right, remains int
	switch {
	case data[repeatPosField1][0] == '(':
		left = mustAtoi(data[repeatPosField3])
		right = mustAtoi(data[repeatPosField2])
		remains = mustAtoi(data[repeatPosField1][1 : len(data[repeatPosField1])-1])
	case data[repeatPosField3][0] == '(':
		left = mustAtoi(data[repeatPosField1])
		right = mustAtoi(data[repeatPosField2])
		remains = mustAtoi(data[repeatPosField3][1 : len(data[repeatPosField3])-1])
	default:
		panic(fmt.Errorf("illegal repeat coordinates: %q", data[repeatPosField1:repeatPosField3+1]))
	}
	return fmt.Sprintf("%s %s %d %d %d", data[repeatTypeField], data[repeatClassField], left, right, remains)
}
