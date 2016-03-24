# Quilt pipeline

## Run repeat annotator

Example with CENSOR

Note that that this stage CENSOR is not completely supported.
CENSOR does not keep repeat coordinates defined by the consensus coordinate system which is required for stitch operation.
This is a work in progress.

```bash
censor -lib ${LIB}-${VER}.fa ${GENOME}.fa
```

Example with RepeatMasker

```bash
RepeatMasker -e wublast -lib ${LIB}-${VER}.fa ${GENOME}.fa
```

## Convert annotator output to GFF

CENSOR

```bash
TYPE=censor
map2gff -lib ${LIB}.fa <${GENOME}-${VER}.map >${GENOME}-${VER}.rep.${TYPE}.gff
```

RepeatMasker

```bash
TYPE=rm
rm2gff <${GENOME}-${VER}.map >${GENOME}-${VER}.rep.${TYPE}.gff
```

## Perform dynamic programming segment collation

```bash
stitch -in ${GENOME}-${VER}.rep.${TYPE}.gff >${GENOME}-${VER}.rep.${TYPE}.stitch.gff
```

## Remove collated segments from input

```bash
tailor -in ${GENOME}-${VER}.rep.${TYPE}.stitch.gff ${GENOME}-${VER}.rep.${TYPE}.gff >${GENOME}-${VER}.rep.${TYPE}.stitch.tailor.gff
```

## Overlay simple and collated segments with collated repeat breaks

```bash
patchwork -in ${GENOME}-${VER}.rep.${TYPE}.stitch.gff ${GENOME}-${VER}.rep.${TYPE}.stitch.gff ${GENOME}-${VER}.rep.${TYPE}.stitch.tailor.gff >${GENOME}-${VER}.rep.${TYPE}.stitch.patchwork.gff
```