#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate bigsnprenv_1911
outdir="/path/to/outdir"
#assume the output folder has the same tree structure from bgen_parallel.sbatch

head -1 $outdir/*/*.prs.tsv | grep -v "^=" | tail -1 | cut -f 1-3 > $outdir/all.prs.tsv
tail -n+2 $outdir/*/*.prs.tsv | grep -v "^=" | sort -V | sed '1{/^ *$/d}' | datamash groupby 1,2 sum 3 >> $outdir/all.prs.tsv