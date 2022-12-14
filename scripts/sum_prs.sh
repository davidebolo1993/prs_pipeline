#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate bigsnprenv_1911
outdir="/path/to/outdir"
matchid="/path/to/bgen.sample"

#assume the output folder has the same tree structure from bgen_parallel.sbatch

#run bgen_parallel
sentence=$(sbatch bgen_parallel.sbatch)

#wait until all the job ids are consumed
stringarray=($sentence)
jobid=(${stringarray[3]})
VAR=0
while [ "$VAR" -eq 0 ]; do
    sleep 1
    nrow=$(squeue -j $jobid | wc -l)
    if [ $nrow -eq 1 ] ; then #only the header of the jobid table
        VAR=1
    fi
done

#sum
head -1 $outdir/*/*.prs.tsv | grep -v "^=" | tail -1 | cut -f 1-3 > $outdir/all.prs.tsv
tail -n+2 $outdir/*/*.prs.tsv | grep -v "^=" | sed '/^$/d' | sort | datamash groupby 1,2 sum 3 | sort -V | cut -f 3 > $outdir/vals.tsv
tail -n+3 $matchid | awk '{print $1,$2}' OFS="\t" > $outdir/ids.tsv
paste -d "\t" $outdir/ids.tsv $outdir/vals.tsv >> $outdir/all.prs.tsv && rm $outdir/vals.tsv $outdir/ids.tsv

