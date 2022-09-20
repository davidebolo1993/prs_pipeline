# PRS pipeline


# Activate conda environment

A conda environment containining all the required packages to run the pipeline is available (bigsnprenv_1911).
Current versions for main packages (bigsnpr/bigstatsr) are reported below.

```
bigsnpr: 1.10.8
bigstatsr: 1.5.6
```

# Run the pipeline

The pipeline consists in a single r script:

```bash
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
		GWAS summary stats or pre-calculated tsv (with header) containing beta scores [required]

	--summarycols=SUMMARYCOLS
		.json file defining columns to use, as in test/sample.json

	-p PHENOTYPE, --phenotype=PHENOTYPE
		tsv file with phenotype (and also covariates, if any) as in test/pheno.tsv

	--heritability=HERITABILITY
		heritability, if known

	-o OUTPUT, --output=OUTPUT
		output prefix [required]

	--threads=THREADS
		computational threads [1]

	--correlation=CORRELATION
		the correlation matrix provided as a pre-computed .rds object [1]

	-h, --help
		Show this help message and exit

```

