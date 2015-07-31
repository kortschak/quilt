// Copyright ©2015 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package main

import (
	"bytes"
	"fmt"
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

func compositesFrom(parts map[partition][]*simple, maxSeparation, workers int) <-chan []composite {
	done := make(chan []composite)
	limit := make(chan struct{}, workers)
	var wg sync.WaitGroup
	go func() {
		for p, recs := range parts {
			p, recs := p, recs
			wg.Add(1)
			limit <- struct{}{}
			go func() {
				defer func() {
					<-limit
					wg.Done()
				}()

				if len(recs) < 2 || recs[0].left == none {
					log.Printf("%v records=%d - skip", p, len(recs))
					return
				}

				sort.Sort(byRightEnd(recs))

				var splits int
				for i, r := range recs[1:] {
					if r.genomic.right-recs[i].genomic.right > maxSeparation {
						splits++
					}
				}
				log.Printf("%v records=%d splits=%d", p, len(recs), splits)

				i := 0
				for j, r := range recs[1:] {
					if r.genomic.right-recs[j].genomic.right > maxSeparation || j == len(recs)-2 {
						n := (j + 2) - i
						if n < 2 {
							continue
						}
						if workers == 1 {
							fmt.Fprintf(os.Stderr, "split size:%d from(right end):%d to:%d\n",
								n, recs[i].genomic.right, recs[j+1].genomic.right)
						}
						c := stitch(recs[i:j+2], costFunction)
						i = j + 2

						if workers == 1 {
							sort.Sort(byGenomeLocation(c))
							for _, f := range c {
								fmt.Fprintf(os.Stderr, "\t%+v\n", f)
							}
						}

						done <- c
					}
				}
			}()
		}
		wg.Wait()
		close(done)
	}()

	return done
}

func stitch(repeats []*simple, cost func(left, right *simple) (score float64, ok bool)) []composite {
	if len(repeats) < 2 {
		return nil
	}

	// Dynamic programming to maximise the score of chained features.
	a := make([][]element, len(repeats))
	for i := 0; i < 2; i++ {
		a[i] = make([]element, len(repeats))
	}
	for i, r := range repeats {
		a[0][i] = element{link: i, score: r.score}
	}
	a[1][0].score = a[0][0].score

	for i := 1; i < len(a); i++ {
		for j := i; j < len(repeats); j++ {
			a[i][j] = max(repeats, a, i, j, cost)
		}
		if i < len(repeats)-1 {
			a[i+1], a[i-1] = a[i-1], nil
		}
	}

	var cmp []composite

	final := a[len(a)-1]
	wasUsed := make([]bool, len(final))
	for i := len(final) - 1; i >= 0; i-- {
		if wasUsed[i] {
			continue
		}
		var (
			c composite
			p int
		)
		c.score = final[i].score
		c.class = repeats[0].class
		for p = i; p != final[p].link; p = final[p].link {
			wasUsed[p] = true
			c.parts = append(c.parts, part{
				name:    repeats[p].name,
				left:    repeats[p].left,
				right:   repeats[p].right,
				genomic: repeats[p].genomic,
			})

		}
		if p == i {
			continue
		}
		c.parts = append(c.parts, part{
			name:    repeats[p].name,
			left:    repeats[p].left,
			right:   repeats[p].right,
			genomic: repeats[p].genomic,
		})

		reverse(c.parts)
		cmp = append(cmp, c)
	}

	return cmp
}

func max(repeats []*simple, a [][]element, i, j int, fn func(left, right *simple) (score float64, ok bool)) element {
	e := a[i-1][j]
	right := repeats[j]
	for k := j - 1; k >= 0; k-- {
		left := repeats[k]
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

func reverse(s []part) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

type element struct {
	link  int
	score float64
}

type composite struct {
	class string
	score float64
	parts parts
}

type partition struct {
	chrom  string
	strand seq.Strand
	class  string
}

func (p partition) String() string {
	return fmt.Sprintf("chr:%s strand:(%v) class:%s", p.chrom, p.strand, p.class)
}

type part struct {
	name        string
	left, right int
	genomic     location
}

// location is a repeat-matching interval.
type location struct {
	chrom  string
	left   int
	right  int
	strand seq.Strand
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

// simple is a masked repeat record.
type simple struct {
	// genomic is the genomic region matched
	// to the the repeat identified below.
	genomic location

	// name and class that the repeat type
	// and class defined by the masker.
	name, class string

	// score is the feature score.
	score float64

	// left and right are the left and right
	// position of the simple alignment in
	// consensus-relative coordinates.
	left, right int
}

const none = -1

func newSimpleRepeat(f *gff.Feature) (*simple, error) {
	repeat := &simple{
		genomic: location{
			left:   f.FeatStart,
			right:  f.FeatEnd,
			chrom:  f.SeqName,
			strand: f.FeatStrand,
		},
	}
	if f.FeatScore != nil {
		repeat.score = *f.FeatScore
	}

	ra := f.FeatAttributes.Get("Repeat")
	if ra == "" {
		return nil, fmt.Errorf("missing repeat tag: file probably not an RM gff.")
	}
	err := repeat.parse(ra)
	if err != nil {
		return nil, fmt.Errorf("failed to parse repeat tag: %v\n", err)
	}
	return repeat, nil
}

func (r *simple) parse(a string) error {
	fields := strings.Split(a, " ")

	r.name = fields[0]
	r.class = fields[1]
	var err error
	if fields[2] != "." {
		r.left, err = strconv.Atoi(fields[2])
		if err != nil {
			return err
		}
		r.left = feat.OneToZero(r.left)
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

type byRightEnd []*simple

func (r byRightEnd) Len() int { return len(r) }
func (r byRightEnd) Less(i, j int) bool {
	iName := r[i].genomic.chrom
	jName := r[j].genomic.chrom
	return iName < jName || (iName == jName && r[i].genomic.right < r[j].genomic.right)
}
func (r byRightEnd) Swap(i, j int) { r[i], r[j] = r[j], r[i] }

type byGenomeLocation []composite

func (c byGenomeLocation) Len() int { return len(c) }
func (c byGenomeLocation) Less(i, j int) bool {
	iName := c[i].parts[0].genomic.chrom
	jName := c[j].parts[0].genomic.chrom
	return iName < jName || (iName == jName && c[i].parts[0].genomic.left < c[j].parts[0].genomic.left)
}
func (c byGenomeLocation) Swap(i, j int) { c[i], c[j] = c[j], c[i] }
