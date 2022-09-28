#!/bin/bash
#SBATCH --job-name=bgen_parallel
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --array=1-?
#SBATCH --mem=450000M
#SBATCH --output=bgen_parallel.log

echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo "RUNNING ON: $(hostname)"
echo "STARTED: $(date)"
echo '----------------------------'

eval "$(conda shell.bash hook)"
conda activate bigsnprenv_1911



summary="/path/to/summary.tsv"
json="/path/to/columns.json"

### assuming .bgen files are in a folder named "test"
### one can simply do ls test/*.bgen > bgen.fofn
### this will generate the proper .fofn file to use as input

bgen="/path/to/bgen.fofn"
outdir="/path/to/output.dir"

line=$(head -n $SLURM_ARRAY_TASK_ID $bgen | tail -1)
ext="${line%%.*}" #get the path without the extension 
var=$(basename $ext) #get the basename

Rscript prs.r --input $line --summary $summary --summarycols $json --threads 32 --output $outdir/$var

echo '----------------------------'
echo "COMPLETED: $(date)"
echo '----------------------------'
