#!/bin/bash
#SBATCH --job-name=multi_bed
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=100000M
#SBATCH --output=multi_bed.log

echo '----------------------------'
echo "RUNNING ON: $(hostname)"
echo "STARTED: $(date)"
echo '----------------------------'

module load plink/2.00_20211217
#make list of *fam files with ls *.fam > famlist.txt

plink2 --bfile <bedfiles.prefix> --merge-list <famlist.txt> --make-bed --out <outputfile.prefix> --threads <threads>

echo '----------------------------'
echo "COMPLETED: $(date)"
echo '----------------------------'

