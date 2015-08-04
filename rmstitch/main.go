// Copyright ©2015 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

// rmstitch performs repeat annotation chaining on RM out data using the repeat
// id field to identify chainable features.
//
// The output of rmstitch is the same format as the output of stitch and is
// intended to be used as a comparison between stitch and the RM repeat chains.
package main

import (
	"bufio"
	"bytes"
	"fmt"
	"log"
	"os"
	"sort"
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

func main() {
	groups := [][]*gff.Feature{nil} // Gain one-based index for this case.
	sc := bufio.NewScanner(os.Stdin)
	for n := 1; sc.Scan(); n++ {
		if n < firstDataLine {
			continue
		}
		f := &gff.Feature{
			Source:         "RepeatMasker",
			Feature:        "repeat",
			FeatFrame:      gff.NoFrame,
			FeatAttributes: gff.Attributes{{Tag: "Repeat"}, {Tag: "ID"}},
		}
		data := strings.Fields(sc.Text())
		err := fill(f, data)
		if err != nil {
			log.Fatalf("parse error on line %d: %v", n, err)
		}
		id, err := strconv.Atoi(data[idField])
		if err != nil {
			log.Fatalf("id parse error on line %d: %v", n, err)
		}
		f.FeatAttributes[1].Value = fmt.Sprint(id)
		switch {
		case id == len(groups):
			groups = append(groups, []*gff.Feature{f})
		case id < len(groups):
			groups[id] = append(groups[id], f)
		default:
			log.Fatalf("id out of order: expect <= %d got %d", len(groups), id)
		}
	}
	err := sc.Err()
	if err != nil {
		log.Fatal(err)
	}

	w := gff.NewWriter(os.Stdout, 60, true)
	for _, g := range groups {
		if len(g) < 2 {
			continue
		}
		sort.Sort(byGenomeLocation(g))
		right := g[0].FeatEnd
		for _, p := range g[1:] {
			if p.FeatEnd > right {
				right = p.FeatEnd
			}
		}
		w.Write(&gff.Feature{
			Source:         "stitch",
			Feature:        "composite",
			SeqName:        g[0].SeqName,
			FeatStart:      g[0].FeatStart,
			FeatEnd:        right,
			FeatFrame:      gff.NoFrame,
			FeatAttributes: attributes(g),
		})
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

func fill(f *gff.Feature, data []string) (err error) {
	defer handlePanic(&err)
	f.FeatScore = mustAtofp(data[swScoreField])
	f.SeqName = data[queryNameField]
	f.FeatStart = feat.OneToZero(mustAtoi(data[queryStartField]))
	f.FeatEnd = mustAtoi(data[queryEndField])
	f.FeatStrand = mustRMtoSane(data[strandField])
	f.FeatAttributes[0].Value = repeatAttribute(data)
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

type byGenomeLocation []*gff.Feature

func (g byGenomeLocation) Len() int { return len(g) }
func (g byGenomeLocation) Less(i, j int) bool {
	iName := g[i].SeqName
	jName := g[j].SeqName
	return iName < jName || (iName == jName && g[i].FeatStart < g[j].FeatStart)
}
func (g byGenomeLocation) Swap(i, j int) { g[i], g[j] = g[j], g[i] }

func attributes(p []*gff.Feature) gff.Attributes {
	var (
		buf   bytes.Buffer
		class string
	)
	for i, e := range p {
		ra := e.FeatAttributes.Get("Repeat")
		f := strings.Fields(ra)
		if i == 0 {
			class = f[1]
			buf.WriteByte('"')
		} else {
			buf.WriteByte('|')
		}
		fmt.Fprintf(&buf, `%s %d %d %d %d`,
			f[0],
			mustAtoi(f[2]), mustAtoi(f[3]),
			feat.ZeroToOne(e.FeatStart), e.FeatEnd,
		)
	}
	buf.WriteByte('"')
	return gff.Attributes{
		{Tag: "Class", Value: `"` + class + `"`},
		{Tag: "Parts", Value: buf.String()},
		{Tag: "ID", Value: p[0].FeatAttributes.Get("ID")},
	}
}
