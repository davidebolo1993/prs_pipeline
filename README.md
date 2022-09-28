# PRS pipeline


# Activate conda environment

A conda environment containining all the required packages to run the pipeline is available (bigsnprenv_1911).
Current versions for main packages (bigsnpr/bigstatsr) are reported below.

```
bigsnpr: 1.10.8
bigstatsr: 1.5.6
```

# Clone the pipeline

```bash
git clone https://github.com/davidebolo1993/prs_pipeline
cd prs_pipeline
```

# Run the script

The pipeline consists in a single r script based on [ldpred2](https://privefl.github.io/bigsnpr/articles/LDpred2.html).

```bash

#srun --nodes=1 --tasks-per-node=1 --partition cpu-interactive --cpus-per-task 8 --pty /bin/bash #this is to enter a computational node

conda activate bigsnprenv_1911
Rscript scripts/prs.r --help

#During startup - Warning message:
#Setting LC_CTYPE failed, using "C" 
#Loading required package: bigstatsr
#Attaching package: 'dplyr'

#The following objects are masked from 'package:data.table':

#    between, first, last

#The following objects are masked from 'package:stats':

#    filter, lag

#The following objects are masked from 'package:base':

#    intersect, setdiff, setequal, union

#R version 4.1.3 (2022-03-10)
#Platform: x86_64-conda-linux-gnu (64-bit)
#Running under: CentOS Linux 8 (Core)

#Matrix products: default
#BLAS/LAPACK: /project/alfredo/conda_envs/bigsnprenv_1911/lib/libopenblasp-r0.3.20.so

#locale:
# [1] LC_CTYPE=C                 LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] tools     stats     graphics  grDevices utils     datasets  methods  
#[8] base     

#other attached packages:
#[1] rjson_0.2.21      optparse_1.7.1    dplyr_1.0.9       fmsb_0.7.3       
#[5] stringr_1.4.0     magrittr_2.0.3    data.table_1.14.2 bigsnpr_1.10.8   
#[9] bigstatsr_1.5.6  

#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.8.3       bigsparser_0.6.1   pillar_1.7.0       compiler_4.1.3    
# [5] bigparallelr_0.3.2 rngtools_1.5.2     iterators_1.0.14   digest_0.6.29     
# [9] lifecycle_1.0.1    tibble_3.1.7       gtable_0.3.0       lattice_0.20-45   
#[13] doRNG_1.8.2        pkgconfig_2.0.3    rlang_1.0.2        Matrix_1.4-1      
#[17] foreach_1.5.2      DBI_1.1.2          cli_3.3.0          flock_0.7         
#[21] parallel_4.1.3     bigassertr_0.1.5   generics_0.1.2     vctrs_0.4.1       
#[25] getopt_1.20.3      grid_4.1.3         tidyselect_1.1.2   cowplot_1.1.1     
#[29] glue_1.6.2         R6_2.5.1           fansi_1.0.3        ggplot2_3.3.6     
#[33] purrr_0.3.4        scales_1.2.0       codetools_0.2-18   ellipsis_0.3.2    
#[37] assertthat_0.2.1   colorspace_2.0-3   utf8_1.2.2         stringi_1.7.6     
#[41] munsell_0.5.0      doParallel_1.0.17  crayon_1.5.1      

Usage: prs.r [options]

Options:
	-m MODEL, --model=MODEL
		model to use, either automatic or grid. Default to automatic

	-i INPUT, --input=INPUT
		.bed file (with companion .bim and .fam files) or (indexed) .bgen file [required]

	-s SUMMARY, --summary=SUMMARY
		GWAS summary stats or pre-calculated .tsv (with header) containing beta scores [required]

	--summarycols=SUMMARYCOLS
		.json file defining columns to use

	-p PHENOTYPE, --phenotype=PHENOTYPE
		.tsv file with phenotype (and also covariates, if any)

	--heritability=HERITABILITY
		heritability, if known

	-o OUTPUT, --output=OUTPUT
		output prefix [required]

	--threads=THREADS
		computing threads [1]

	--correlation=CORRELATION
		the correlation matrix provided as a pre-computed .rds object

	-h, --help
		Show this help message and exit


#run with test data

Rscript prs.r --input test/public-data3.bed --summary test/public-data3-sumstats.txt --summarycols test/public-data3.json --threads 8 --output test/output/public-data3
```

## Input files

### Model

Available models are "grid" and "automatic". Default to "automatic".
Detailed informations on each model are available in the [ldpred2 tutorial introduction](https://privefl.github.io/bigsnpr/articles/LDpred2.html) and the [paper](https://doi.org/10.1093/bioinformatics/btaa1029).

### Bed/bim/fam or Bgen

The pipeline accepts either a unique .bed file (with matching .bim/.fam) or a unique (indexed) .bgen file.
Further informations on these file formats are available [here](https://www.cog-genomics.org/plink/2.0/formats).
Examples are provided in the [test](test/) folder.

#### Build sinlge .bed/.bim/.fam or single .bgen from multiple inputs

Companion scripts are released in the [scripts](scripts/) folder to generate a single input .bed/.bim/.fam or .bgen if not available already (from chromosome-specific .bed/.bgen files).
For .bgen files, use [this script](scripts/multi_bgen.r). Positional arguments are a directory containing (indexed) .bgen files, the output, merged, .rds object and the number of computational threads.
The script also generates companion .bk, .bgen and .bgi files - the (indexed) .bgen files is actually empty but can be given as input to the pipeline that will load the informative .rds/.bk files instead.
A sbatch script that run the above mentioned R script is also available - adjust the parameters accordingly to the input files.
For .bed files, those can be merged using [plink2](https://www.cog-genomics.org/plink/2.0/). A dedicated sbatch script to be used as template is available in the [script](scripts/) folder.


#### Run the pipeline in parallel on multiple .bgen inputs

The pipeline can be run on .bgen files in parallel on a cluster - one node for each chromosome using [this script](scripts/bgen_parallel.sbatch). .bgen files must be provided in file-of-file-names format (.fofn) and the other parameters in the scipt should be set accordingly to the experimental setting. The sbatch script generates one output folder for each chromosome. Results for each chromosome can be summed up using this simple [bash script](scripts/sum_prs.sh) that produces a `all.prs.tsv` table in the output folder (the same used above).

### Summary

A standard summary statics file or a file with pre-computed scores can be given as input.
Summary statistics file format from GWAS catalog is described [here](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format#:~:text=Summary%20statistics%20are%20defined%20as,row%20for%20each%20variant%20analysed.).
PGS scoring file format from PGS catalog is described [here](https://www.pgscatalog.org/downloads/#:~:text=Formatted%20Files,-Format%3A%202.0&text=Each%20scoring%20file%20(variant%20information,gz%20).).
An example summary statistic file is available in the [test](test/) folder.

### Json

While the other inputs to the pipeline are standard files, the input .json file is used to handle summary statistics and matrices with pre-computed scores generated with different methods - so that the pipeline actually knows which column refers to which value.
The input .json file follows standard [json](https://en.wikipedia.org/wiki/JSON) specification, where, for each "key":"value" pair, users must specify:
- the name of the column as it is used within the pipeline ("key") - do not change this
- the name of the corresponding column in the sequencing summary file ("value") - adjust accordingly
- "n_eff" can be either a column in the summary file or a single mumeric value
- "is_beta_precomp" must be logical ("TRUE", "FALSE") and indicates whether to consider the input sequencing summary a true sequencing summary or a tabular file with precomputed beta scores (from PGS catalog, for instance)
Some of the values in the .json file may be missing - for instance, it's possible to specify chromosome and position instead of rsid or rsid alone (or both) but at least one of those combinations must be present in order to match the input with the information in .bed/.bim/.fam or .bgen format.
An example is provided in the [test](test/) folder.

### Phenotype

Table containing the phenotype as a third column (can be either categorical or continuous) and possible covariates afterwards. First and second columns can be anything (family and individual ID, for instance).
An example is provided in the [test](test/) folder.
 
### Heritability

Value of the heritability, if known. It will be calculated otherwise.

### Correlation matrix

If a correlation matrix for matching summary statistics and .bed/.bim/.fam files is available (because stored by a previous run of the pipeline), skip the computation of the correlation matrix

## Output files

The pipeline stores the correlation matrix and the ld as R objects (.rds format) - named with the output prefix followed by `.corr.rds` and `.ld.rds` respectively. The correlation matrix can be provided as input to skip the computation of the correlation matrix itself - which is time consuming.

### Automatic model

When using the automatic model without a table of phenotypes as input, the pipeline returns beta scores and predictions from the automatic model - named with the output prefix followed by `.auto.beta_scores.tsv` and `.auto.prs.tsv` respectively
When using the automatic model with a table of phenotypes as input, the pipeline further tests the generated predictions on the provided phenotype. Performances are stored as summary files - named with the output prefix followed by `.auto.summary.tsv` and `.auto.summary.rds` respectively

### Grid model

Grid model can only be run with a table of phenotypes as input. 
The pipeline returns beta scores and predictions from the grid model - named with the output prefix followed by `.auto.grid_scores.tsv` and `grid.prs.tsv` respectively - as well as performances of tested models as a tab-separated value file -  named with the output prefix followed by `.grid.allsummary.tsv` - and those of the best model as an additional R object - named with the output prefix followed by `.grid.bestsummary.rds`.

## Sbatch script

A [sbatch script](scripts/prs.sbatch) is provided to run the pipeline using an entire computational node. Parameters in the sbatch file should be adjusted coherently.
