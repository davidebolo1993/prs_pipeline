#!/bin/bash
#SBATCH --job-name=prs_pipeline
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=450000M
#SBATCH --output=prs_pipeline.log

echo '----------------------------'
echo "RUNNING ON: $(hostname)"
echo "STARTED: $(date)"
echo '----------------------------'

eval "$(conda shell.bash hook)"
conda activate bigsnprenv_1911

Rscript prs.r --input ../test/public-data3.bed --summary ../test/public-data3-sumstats.txt --summarycols ../test/public-data3.json --threads 32 --output ../test/output/public-data3

echo '----------------------------'
echo "COMPLETED: $(date)"
echo '----------------------------'




