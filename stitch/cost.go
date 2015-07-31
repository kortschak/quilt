package main

import (
	"math"

	"github.com/biogo/biogo/seq"
)

// maxSpan is the maximum distance we will examine left of our current right element.
const maxSpan = 1e5

// Tolerance values specify the width of troughs in the cost function.
// Values for tolerance are greater than or equal to zero.
const (
	gOverlapTolerance = 2
	rOverlapTolerance = 1
	concordTolerance  = 0.5
)

func costFunction(left, right *record) (score float64, ok bool) {
	// Short circuit if we got here without a strand or
	// if the distance between the sorted ends is greater
	// than our maximum span.
	if right.genomic.strand == seq.None || right.genomic.right-left.genomic.right > maxSpan {
		return math.Inf(-1), false
	}

	gOverlap := left.genomic.right - right.genomic.left
	var rOverlap int
	if right.genomic.strand == seq.Plus {
		rOverlap = left.right - right.left
	} else {
		rOverlap = right.right - left.left
	}

	cost := math.Pow(float64(abs(gOverlap)), gOverlapTolerance) *
		math.Pow(float64(abs(rOverlap)), rOverlapTolerance) *
		math.Pow(float64(abs(gOverlap-rOverlap)), concordTolerance)
	switch {
	// Special-case immediately adjacent intervals.
	case rOverlap == 0:
		if gOverlap < 0 {
			cost = float64(gOverlap) * 2
		} else {
			cost = float64(gOverlap) * 100
		}
	case gOverlap == 0:
		if rOverlap < 0 {
			cost = float64(rOverlap) * 40
		} else {
			cost = float64(rOverlap) * 100
		}
	case rOverlap < 0 && gOverlap < 0:
		// Separated parts.
		cost *= 10
	case rOverlap < 0 && gOverlap > 0:
		// Overlapping genomic segments from different non-overlapping element parts.
		cost *= 100
	case rOverlap > 0 && gOverlap < 0:
		// Overlapping element parts from different non-overlapping genome segments.
		cost *= 100
	default:
		// Co-overlaps.
		cost *= 10
	}

	return left.score - math.Abs(cost), true
}
