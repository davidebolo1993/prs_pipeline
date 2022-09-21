#!/bin/bash
#SBATCH --job-name=multi_bgen
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=450000M
#SBATCH --output=multi_bgen.log

echo '----------------------------'
echo "RUNNING ON: $(hostname)"
echo "STARTED: $(date)"
echo '----------------------------'

eval "$(conda shell.bash hook)"
conda activate bigsnprenv_1911

Rscript multi_bgen.r /group/soranzo/bigsnprtest_ukbiobk/test_prs_fede/testdir allbgen.rds 32

echo '----------------------------'
echo "COMPLETED: $(date)"
echo '----------------------------'

