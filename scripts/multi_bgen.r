library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(data.table)
library(magrittr)
library(stringr)
library(fmsb)
library(dplyr)
library(optparse)
library(tools)
library(rjson)
library(crayon)
sessionInfo()

prepare_multi_bgen<-function(x,rds_file,threads) {
  
  backing_file<-str_replace(file.path(rds_file), ".rds", "")
  fake_bgen<-str_replace(file.path(rds_file), ".rds", ".bgen")
  fake_bgi<-str_replace(file.path(rds_file), ".rds", ".bgen.bgi")

  if (file.exists(rds_file)) {

    now<-Sys.time()
    stop(red('[',now,'][Error] output file already exists. Nothing to do.'))

  }

  snp_bgi<-list()
  bgens<-list.files(inputdir, pattern=".bgen$", full.names=T)

  if (length(bgens) == 0) {

    now<-Sys.time()
    stop(red('[',now,'][Error] No .bgen files found in input folder. Nothing to do.'))

  }

  for (i in c(1:length(bgens))) {

    bgi_file<-file.path(str_replace(bgens[i], ".bgen", ".bgen.bgi"))

    if (! file.exists(bgi_file)) {

        now<-Sys.time()
        stop(red('[',now,'][Error] missing .bgi for file ',bgens[i]))

    }

    now<-Sys.time()
    message('[',now,'][Message] reading .bgi file ',bgi_file)
    bgi<-snp_readBGI(bgi_file,snp_id=NULL)
    snps_ids<-paste(bgi$chromosome, bgi$position, bgi$allele1, bgi$allele2, sep="_")
    snp_bgi[[i]]<-snps_ids
  
  }

  now<-Sys.time()
  message('[',now,'][Message] done')
  message('[',now,'][Message] creating .rds object')
  snp_readBGEN(bgens, backingfile=backing_file, list_snp_id=snp_bgi, ncores=threads, read_as ="dosage")
  now<-Sys.time()
  message('[',now,'][Message] done')

  file.create(fake_bgen)
  file.create(fake_bgi) 

}

args<-commandArgs(trailingOnly=TRUE)
inputdir<-file.path(args[1])
outputfile<-file.path(args[2])
threads<-as.numeric(args[3])

if (file_ext(outputfile) != "rds") {

    now<-Sys.time()
    stop(red('[',now,'][Error] output file must be a .rds object'))

}

outdir<-dirname(file.path(outputfile))
dir.create(outdir, recursive=TRUE,showWarning=FALSE)
prepare_multi_bgen(inputdir,outputfile,threads)