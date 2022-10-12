#!/usr/bin/env Rscript

# how to see installed packages in a specific conda environment being previously activated

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

# run functions before running the pipeline

check_model<-function(x) {

  if (!x %in% c('automatic', 'grid')) {

    now<-Sys.time()
    stop(red('[',now,'][Error] specified model can be either automatic or grid'))

  }

}

check_input<-function(x) {

  if (is.null(x)) {

  now<-Sys.time()
  stop(red('[',now,'][Error] a .bed/.bgen file must be provided as input'))

  } else {

    if (! file.exists(file.path(x))) {

      now<-Sys.time()
      stop(red('[',now,'][Error] the .bed/.bgen file must exist'))

    } else {

      if (file_ext(x) == "bed") {

        now<-Sys.time()
        message('[',now,'][Message] detected .bed input') 
        
        is_bgen<-FALSE
        bim_f<-str_replace(x, ".bed", ".bim")
        fam_f<-str_replace(x, ".bed", ".fam")

        if ((! file.exists(file.path(bim_f))) | (! file.exists(file.path(fam_f)))) {

          now<-Sys.time()
          stop(red('[',now,'][Error] .bim and .fam files must exist in .bed file location'))

        }

      } else if (file_ext(x) == "bgen") {

        now<-Sys.time()
        message('[',now,'][Message] detected .bgen input')
        is_bgen<-TRUE

        if (!file.exists(paste0(x, ".bgi"))) {

          now<-Sys.time()
          stop(red('[',now,'][Error] .bgen file must have a companion .bgi index. Use bgenix with .bgen file to create a proper index'))
        }

      } else {

        now<-Sys.time()
        stop(red('[',now,'][Error] unknown input file extension'))

      }

    }

  }

  is_bgen

}


check_output<-function(x) {

  if (is.null(x$output)) {

    now<-Sys.time()
    stop('[',now,'][Error] output prefix must be specified') 

  } else {

    #x$output<-normalizePath(file.path(x$output))
    out_dir<-dirname(x$output)
    dir.create(out_dir, recursive=TRUE,showWarning=FALSE)
  }

}

check_summary<-function(x) {

  if (is.null(x)) {

    now<-Sys.time()
    stop(red('[',now,'][Error] summary statistics or tsv with pre-computed beta scores must be provided'))

  } else {
    
    if (! file.exists(file.path(x))) {

      now<-Sys.time()
      stop(red('[',now,'][Error] summary statistics or tsv with pre-computed beta scores file must exist'))

    } 
  } 
}


# function for the internal train - maybe it will be uncommented in the future
# check_train<-function(x) {

#  train_test<-FALSE

#  if (!is.null(x)) { #requested train-test

#    val_type<-as.numeric(x)

#    if (is.na(val_type) | val_type > 100) { #incorrect number

#      now<-Sys.time()
#      stop('[',now,'][Error] invalid percentage for training') 
          
#    } else { #is valid number

#      train_test<-TRUE

#    }

#  }

#  train_test

#}


load_summary<-function(x, cols, threads) {

  stat_file<-file.path(x)
  cols_file<-file.path(cols)
  jfile<-fromJSON(file = cols_file)
  names_cols<-names(jfile)[c(1:length(jfile)-1)]
  select_values<-as.character(jfile)[c(1:length(jfile)-1)]
  precomp_betas<-type.convert(as.character(jfile)[length(jfile)],as.is=TRUE)
  
  if (class(jfile$n_eff) == "numeric") {

    n_val<-jfile$n_eff
    jfile<-jfile[which(names(jfile) != "n_eff")]
    select_values<-as.character(jfile)[c(1:length(jfile)-1)]
    names_cols<-names(jfile)[c(1:length(jfile)-1)]

  }

  sumstats<-data.frame(bigreadr::fread2(stat_file, nThread=threads, select=select_values))
  colnames(sumstats)<-names_cols

  if (! "n_eff" %in% colnames(sumstats)) {

    sumstats$n_eff<-n_val

  }

  if (!isTRUE(precomp_betas)) { 

    now<-Sys.time()
    message('[',now,'][Message] assuming summary stats file')

  } else {

    now<-Sys.time()
    message('[',now,'][Message] assuming pre-computed beta scores')

  }

  list(sumstats,precomp_betas)

}

load_bed<-function(x, threads) {

  now<-Sys.time()
  message('[',now,'][Message] reading .bed/.bim/.fam files')

  bed_file<-file.path(x)
  backing_file<-str_replace(bed_file, ".bed", "") #cut extension out
  bk_file<-file.path(str_replace(bed_file, ".bed", ".bk"))
  rds_file<-file.path(str_replace(bed_file, ".bed", ".rds"))

  if (!file.exists(bk_file)) {

    snp_readBed2(bed_file, backingfile=backing_file, ncores=threads) #this generate a .rds object that can be loaded into the environment

  }

  obj.bigSNP <- snp_attach(rds_file)
  obj.bigSNP

}


load_bgen<-function(x,threads) {

  now<-Sys.time()
  message('[',now,'][Message] reading .bgen file')
  
  bgen_file<-file.path(x)
  backing_file<-str_replace(bgen_file, ".bgen", "") #cut extension out
  bk_file<-file.path(str_replace(bgen_file, ".bgen", ".bk"))
  rds_file<-file.path(str_replace(bgen_file, ".bgen", ".rds"))
  bgi_file<-file.path(str_replace(bgen_file, ".bgen", ".bgen.bgi"))

  if (!file.exists(bk_file)) {

    bgi<-snp_readBGI(bgi_file,snp_id=NULL)
    snps_ids<-list(paste(bgi$chromosome, bgi$position, bgi$allele1, bgi$allele2, sep="_")) #do we want this or an external table?
    snp_readBGEN(bgen_file, backingfile=backing_file, list_snp_id=snps_ids, ncores=threads, read_as ="dosage")

  }

  obj.bigSNP <- snp_attach(rds_file)
  obj.bigSNP

}

check_phenotype<-function(x) {

  if (is.null(x$phenotype)) {

    now<-Sys.time()
    message(yellow('[',now,'][Warning] missing table of phenotypes/covariates. Model coherced to automatic even when grid is specified'))
    x$model<-"automatic"

    #if (! is.null(x$train)) {

      #now<-Sys.time()
      #stop('[',now,'][Error] cannot peform train/test evaluation with automatic model') 
    
    #}

  } else {

    if (! file.exists(file.path(x$phenotype))) {

      now<-Sys.time()
      stop(red('[',now,'][Error] if provided, table of phenotypes/covariates must exist'))

    }
  }

  x$model
}

check_predictions<-function(x,pheno) {

  sub_pheno<-pheno[,c(3:ncol(pheno))]

  if (class(sub_pheno) == "numeric") { #this is a vector

    sub_pheno<-cbind(sub_pheno)
    colnames(sub_pheno)<-"V3"

  }

  sub_val<-data.frame(cbind(sub_pheno,x))
  params<-list()

  if (length(unique(sub_val$V3)) == 2) { #binary trait  #maybe change to length(unique(subpheno$V3))

    res<-summary(glm(V3 ~ ., data=sub_val, family="binomial"))
    params$score<-res$coef["x",3]
    params$prt<-res$coef["x",4]
    params$r2<-1-(res$deviance/res$null.deviance) #calculate (adjusted) r2

  } else {
    
    res<-summary(lm(V3 ~ ., data=sub_val))
    params$score<-res$coef["x",3]
    params$prt<-res$coef["x",4]
    params$r2<-res$adj.r.squared
  
  }

  list(res,params)

}

get_complement<-function(nuc) {

  x<-"A"
  if (nuc=="A") x<-"T"
  if (nuc=="T") x<-"A"
  if (nuc=="C") x<-"G"
  if (nuc=="G") x<-"C"
  x

}

recover_missing_cols<-function(sumstats,map) {

  av_cols<-colnames(sumstats)
  #error if a1 not there ?   

  if ("rsid" %in% colnames(sumstats)) {

    sub_sumstats<-sumstats[which(sumstats$rsid  %in% map$rsid),]
    sub_map<-map[which(map$rsid  %in% sumstats$rsid),]
    sub_sumstats <- sub_sumstats[match(sub_map$rsid,sub_sumstats$rsid),]

    if (! "chr" %in% colnames(sumstats)) { #derive from rsid
      
      sub_sumstats$chr <- sub_map$chr
    
    }

    if (! "a0" %in% colnames(sumstats)) { #derive from rsid and a1

      a0<-rep(".", nrow(sub_map))
      
      for (i in 1:nrow(sub_map)) {

        a0m<-sub_map$a0[i]
        a1m<-sub_map$a1[i]

        if (sub_sumstats$a1[i] %in% c(a0m,a1m)) { #they appear as is, no need to reverse

          a0[i]<-ifelse(sub_sumstats$a1[i] == a0m, a1m, a0m)

        } else {

          a0[i]<-ifelse(get_complement(sub_sumstats$a1[i]) == a0m, get_complement(a1m), get_complement(a0m))

        }
      
      }

      sub_sumstats$a0<-a0

    }

    if (! "pos" %in% colnames(sumstats)) { #derive from rsid

      sub_sumstats$pos<-sub_map$pos
    
    }

  } else { #we do not have rsid but we have chr pos

    #create fake id
    map$fake_id<-paste(map$chr,map$pos,sep="_")
    sumstats$fake_id<-paste(sumstats$chr,sumstats$pos,sep="_")

    #exclude those that are duplicated
    sumstats<-sumstats[!(duplicated(sumstats$fake_id) | duplicated(sumstats$fake_id, fromLast=T)),]
    map<-map[!(duplicated(map$fake_id) | duplicated(map$fake_id, fromLast=T)),]
    
    #use fake_id instead of rsid
    sub_sumstats<-sumstats[which(sumstats$fake_id  %in% map$fake_id),]
    sub_map<-map[which(map$fake_id  %in% sumstats$fake_id),]
    sub_sumstats <- sub_sumstats[match(sub_map$fake_id,sub_sumstats$fake_id),]

    sub_sumstats$rsid<-sub_map$rsid

    if (! "a0" %in% colnames(sumstats)) { #derive from rsid and a1

      a0<-rep(".", nrow(sub_map))
      
      for (i in 1:nrow(sub_map)) {

        a0m<-sub_map$a0[i]
        a1m<-sub_map$a1[i]

        if (sub_sumstats$a1[i] %in% c(a0m,a1m)) { #they appear as is, no need to reverse

          a0[i]<-ifelse(sub_sumstats$a1[i] == a0m, a1m, a0m)

        } else {

          a0[i]<-ifelse(get_complement(sub_sumstats$a1[i]) == a0m, get_complement(a1m), get_complement(a0m))

        }
      
      }

      sub_sumstats$a0<-a0

    }

  }

  sub_sumstats

}

check_correlation<-function(x) {

  has_corr<-TRUE

  if (is.null(x)) {

    has_corr<-FALSE

  } else {

    if (! file.exists(file.path(x))) {

      now<-Sys.time()
      stop(red('[',now,'][Error] if provided, correlation matrix must exist'))

    } else {

        if (file_ext(file.path(x)) != "rds") {

            now<-Sys.time()
            stop(red('[',now,'][Error] if provided, correlation matrix must be provided in .rds file'))
        }

        if (! file.exists(str_replace(file.path(x), ".rds", ".sbk"))) {

            now<-Sys.time()
            stop(red('[',now,'][Error] if provided, correlation matrix must have a supporting .sbk file - same name, different extension - in the same directory'))

        }

        if (! file.exists(str_replace(file.path(x), ".corr.rds", ".ld.rds"))) {

            now<-Sys.time()
            stop(red('[',now,'][Error] if provided, correlation matrix must have a matchind .ld.rds file - different name, same extension - in the same directory'))

        }

    }

  }

  has_corr
}

check_index<-function(x) {

  has_index<-FALSE

  if (is.null(x)) {

    now<-Sys.time()
    message(yellow('[',now,'][Warning] no index provided. If an external correlation matrix is not provided, will use all the individuals for generating one'))

  } else {

    if (! file.exists(file.path(x))) {

      now<-Sys.time()
      stop(red('[',now,'][Error] if provided, index table must exists'))

    } else {

      has_index<-TRUE

    }
  }

  has_index

}


# start the prs script 
# create the list with all necessary files
options(warn = -1)

option_list = list(
  make_option(c('-m', '--model'), action='store', type='character', help='model to use, either automatic or grid. Default to automatic', default='automatic'),
  make_option(c('-i', '--input'), action='store', type='character', help='.bed file (with companion .bim and .fam files) or (indexed) .bgen file [required]'),
  make_option(c('-s', '--summary'), action='store', type='character', help='GWAS summary stats or pre-calculated .tsv (with header) containing beta scores [required]'),
  make_option(c('--summarycols'), action='store', type='character', help='.json file defining columns to use'),
  make_option(c('-p', '--phenotype'), action='store', type='character', help='.tsv file with phenotype (and also covariates, if any)'),
  make_option(c('--heritability'), action='store', type='numeric', help='heritability, if known'),
  make_option(c('-o', '--output'), action='store', type='character', help='output prefix [required]'),
  make_option(c('--threads'), action='store', type='numeric', help='computing threads [1]', default=1),
  make_option(c('--correlation'), action='store', type='character', help='the correlation matrix provided as a pre-computed .rds object'),
  make_option(c('--index'), action='store', type='character', help='indexes to limit correlation matrix calculation to')
  #make_option(c('-t', '--train'), action='store', type='character', help='train percentage for training-testing - internal validation'),
  #make_option(c('-x', '--external'), action='store', type='character', help='external validation set. Comma-separated .bed (or .bgen) and .pheno tsv. .bed should have accompanying .bim and .fam')
)

opt = parse_args(OptionParser(option_list=option_list))

#quickly check all input. This functions do not load anything in r

now<-Sys.time()
message('[',now,'][Message] performing pre-flight checks') 

check_model(opt$model)
is_bgen<-check_input(opt$input)
check_summary(opt$summary)
opt$model<-check_phenotype(opt)
check_output(opt)
has_corr<-check_correlation(opt$correlation)
#train_test<-check_train(opt$train)
#external_validation<-check_external(opt$external)
has_index<-check_index(opt$index)

now<-Sys.time()
message('[',now,'][Message] done') 
message('[',now,'][Message] loading data') 
message('[',now,'][Message] reading  .bed/.bgen')

if (!is_bgen) {

  obj.bigSNP<-load_bed(opt$input, threads=opt$threads)

} else {

  obj.bigSNP<-load_bgen(opt$input,threads=opt$threads)

}

G <- tryCatch({
  snp_fastImputeSimple(obj.bigSNP$genotypes)
  }, error = function(e) {
  obj.bigSNP$genotypes
})


if (is_bgen) {

  obj.bigSNP$fam <- snp_fake(n = nrow(G), m = 1)$fam

}

CHR <- as.integer(obj.bigSNP$map$chromosome) #this is somehow necessary for .bgen files, not for bed
POS <- obj.bigSNP$map$physical.pos

now<-Sys.time()
message('[',now,'][Message] done')

if (!is.null(opt$phenotype)) {

  message('[',now,'][Message] reading table of phenotypes')

  pheno<-data.frame(fread(opt$phenotype))
  colnames(pheno)<-paste0("V", c(1:ncol(pheno))) #be sure to have proper name for pheno columns

  now<-Sys.time()
  message('[',now,'][Message] done')
    
}

#not tested yet
#if (train_test) {

  #now<-Sys.time()
  #message('[',now,'][Warning] provided a percentage for training. Model coherced to grid even when automatic is specified')
  #opt$model<-"grid" 
  #n_val<-round(nrow(G)/100*as.numeric(opt$train))
  #ind.val<-sample(nrow(G), n_val)
  #ind.test<-setdiff(rows_along(G), ind.val)
  
#} else {

ind.test<-ind.val<-c(1:nrow(G))

#}

message('[',now,'][Message] reading summary statistics')

stats<-load_summary(opt$summary, opt$summarycols, opt$threads)
sumstats<-stats[[1]]
beta_is_precomp<-stats[[2]]

now<-Sys.time()
message('[',now,'][Message] done')
message('[',now,'][Message] matching variants between genotype data and summary statistics - or previously computed beta scores')

if (is_bgen) {

  map<-data.frame(obj.bigSNP$map)
  map$chromosome<-as.numeric(map$chromosome) #this is character otherwise
  map<-map[c("chromosome", "rsid", "physical.pos", "allele1", "allele2")]

} else { #is .bed

  map<-obj.bigSNP$map
  map<-map[c("chromosome", "marker.ID", "physical.pos", "allele1", "allele2")]

}

colnames(map)<-c("chr", "rsid", "pos", "a1", "a0")

#take care of incorrect rsid (?) how frequently does this happen?

rsid_map<-map$rsid
idsl<-strsplit(rsid_map, "_")
ids_<-sapply(idsl,"[[",1)
map$rsid<-ids_

#check columns and add missing fields in sumstats/pre-computed betas. Rsid is necessary

sumstats<-recover_missing_cols(sumstats,map)
sumstats$chr<-as.numeric(sumstats$chr)
df_beta<- snp_match(sumstats, map,join_by_pos = FALSE) #match by rsid for the time being

now<-Sys.time()
message('[',now,'][Message] done')

if (!beta_is_precomp) {

  if (!has_corr) { 

    now<-Sys.time()
    message('[',now,'][Message] computing correlation between variants')
    POS2 <- snp_asGeneticPos(CHR, POS, dir = dirname(opt$output), ncores = opt$threads)
    
    if (has_index) {
    
      now<-Sys.time()
      message('[',now,'][Message] loading indexes to subset correlation matrix')
      individual_idx<-as.integer(fread(opt$index)$V1)
    
    } else {
    
      individual_idx<-c(1:nrow(G))

    }

    for (chr in 1:22) {

        # print(chr)
        ind.chr <- which(df_beta$chr == chr) #subsample G

        if (length(ind.chr) == 0) next

        ## indices in 'G'
        ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
        corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000, ind.row=individual_idx, infos.pos = POS2[ind.chr2], ncores = opt$threads)

        if ((chr == 1) | (length(unique(df_beta$chr)) == 1)) { #also if we have a single chromosome
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, backingfile=file.path(paste0(opt$output, ".corr")), compact = TRUE)
        } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
        }
    }

    now<-Sys.time()
    message('[',now,'][Message] done')
    message('[',now,'][Message] storing correlation matrix and ld values to files')

    saveRDS(corr, file=file.path(paste0(opt$output, ".corr.rds")))
    saveRDS(ld, file=file.path(paste0(opt$output, ".ld.rds")))

  } else {

    now<-Sys.time()
    message('[',now,'][Message] loading precomputed correlation values')
    corr<-readRDS(file.path(opt$correlation))
    ld<-readRDS(file.path(str_replace(opt$correlation, ".corr.rds", ".ld.rds")))
  
  }

  now<-Sys.time()
  message('[',now,'][Message] done')

}

if (is.null(opt$heritability)) {

    now<-Sys.time()
    message('[',now,'][Message] estimating h2')

    ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, sample_size = df_beta$n_eff, blocks = NULL))
    h2_est <- ldsc[["h2"]]

    now<-Sys.time()
    message('[',now,'][Message] done')

  } else {

    h2_est<-opt$heritability

}

if (!beta_is_precomp) {

  if (opt$model == "automatic") {

    beta_tab_file<-file.path(paste0(opt$output,'.auto.beta_scores.tsv'))
    prs_tab_file<-file.path(paste0(opt$output,'.auto.prs.tsv'))
    
    summary_file<-file.path(paste0(opt$output,'.auto.summary.rds')) #use only if pheno provided
    summarytab_file<-file.path(paste0(opt$output,'.auto.summary.tsv')) #use only if pheno provided

    #summary_file2<-file.path(paste0(opt$output,'.auto.summary.external.rds')) #use only if ext
    #summarytab_file2<-file.path(paste0(opt$output,'.auto.summary.external.tsv')) #use only if ext

    now<-Sys.time()
    message('[',now,'][Message] running automatic model')

    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                  vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
                  allow_jump_sign = FALSE, shrink_corr = 0.95,
                  ncores = opt$threads)
    range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
    keep <- (range > (0.9 * quantile(range, 0.9)))
    beta<-rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
    pred <- big_prodVec(G, beta, ind.col = df_beta[["_NUM_ID_"]])

  } else { #is grid

    beta_tab_file<-file.path(paste0(opt$output,'.grid.beta_scores.tsv'))
    prs_tab_file<-file.path(paste0(opt$output,'.grid.prs.tsv'))
    
    #summary_file<-file.path(paste0(opt$output,'.grid.summary.rds'))
    #summarytab_file<-file.path(paste0(opt$output,'.grid.summary.tsv'))
    
    rank_file<-file.path(paste0(opt$output,'.grid.allsummary.tsv'))
    best_model_file<-file.path(paste0(opt$output,'.grid.bestsummary.rds'))

    #summary_file2<-file.path(paste0(opt$output,'.grid.summary.external.rds')) #use only if ext
    #summarytab_file2<-file.path(paste0(opt$output,'.grid.summary.external.tsv')) #use only if ext

    now<-Sys.time()
    message('[',now,'][Message] running grid model')

    h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
    p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
    params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
    beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = opt$threads)
    pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
    sub_pheno<-pheno[ind.val,c(3:ncol(pheno))] #it may be equal to pheno if no split

    if (length(unique(pheno$V3)) == 2) { #binary trait 

      res <- apply(pred_grid[ind.val, ], 2, function(x) {
        if (all(is.na(x))) {
          return(NA)
        } else {
          sub_val<-cbind(sub_pheno,x)
          colnames(sub_val)[1]<-"V3"
          summary(glm(V3 ~ ., data=data.frame(sub_val), family="binomial"))
        }
      })

      params$score<-do.call(c,lapply(res, function(x) x$coef["x",3]))
      params$prt<-do.call(c,lapply(res, function(x) x$coef["x",4]))
      params$r2<-1-(do.call(c,lapply(res, function(x) x$deviance)/do.call(c,lapply(res, function(x) x$null.deviance)))) #calculate (adjusted) r2

    } else { #quantitative/continuous trait

      res <- apply(pred_grid[ind.val, ], 2, function(x) {
        if (all(is.na(x))) {
          return(NA)
        } else {
          sub_val<-cbind(sub_pheno,x)
          colnames(sub_val)[1]<-"V3"
          summary(lm(V3 ~ ., data=data.frame(sub_val)))
        }
      })

      params$score<-do.call(c,lapply(res, function(x) x$coef["x",3]))
      params$prt<-do.call(c,lapply(res, function(x) x$coef["x",4]))
      params$r2<-do.call(c,lapply(res, function(x) x$adj.r.squared))

    }

    rank<- params %>%
      mutate(id = row_number()) %>%
      # filter(sparse) %>% not used. But leave here for future reference  
      arrange(desc(score))

    rank$snps<-nrow(df_beta)

    best_beta_grid_id<- rank %>%
      slice(1) %>%
      pull(id)    

    summary_best_beta_grid<-res[[best_beta_grid_id]]

    now<-Sys.time()
    message('[',now,'][Message] storing tested and best model summary to .tsv and .rds object')

    fwrite(rank, file=rank_file, col.names=T, row.names=F, sep="\t")
    saveRDS(summary_best_beta_grid, file = best_model_file)
    
    beta<-beta_grid[,best_beta_grid_id]

    #if (train_test) {

      #pred_test<-big_prodVec(G, beta, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
      #pred_train<-big_prodVec(G, beta, ind.row = ind.val, ind.col = df_beta[["_NUM_ID_"]])
      #pred<-c(pred_train,pred_test)
      #type_pred<-c(rep("TRAIN", length(pred_train)), rep("TEST", length(pred_test)))     
    
    #} else { #no validation requested or file with external validation given

      pred<-big_prodVec(G, beta, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]]) #here ind.test==ind.val

    #}

  }

} else { #betas as pre-computed values

  now<-Sys.time()
  message('[',now,'][Message] using pre-computed beta scores')

  beta_tab_file<-file.path(paste0(opt$output,'.precomp.beta_scores.tsv'))
  prs_tab_file<-file.path(paste0(opt$output,'.precomp.prs.tsv'))
  
  summary_file<-file.path(paste0(opt$output,'.precomp.summary.rds')) #use only if pheno provided
  summarytab_file<-file.path(paste0(opt$output,'.precomp.summary.tsv')) #use only if pheno provided
  
  #summary_file2<-file.path(paste0(opt$output,'.precomp.summary.external.rds')) #use only if ext
  #summarytab_file2<-file.path(paste0(opt$output,'.precomp.summary.external.tsv')) #use only if ext

  beta<-df_beta$beta
  pred <- big_prodVec(G, beta, ind.col = df_beta[["_NUM_ID_"]])

}

beta_tab<-data.frame(chr=df_beta$chr, rsid=df_beta$rsid,  pos=df_beta$pos, a0=df_beta$a0, a1=df_beta$a1, beta=beta)

#if (train_test) {

  #prs_tab<-data.frame(FID=obj.bigSNP$fam$family.ID,IID=obj.bigSNP$fam$sample.ID,PRED=pred, TYPE=type_pred)

#} else {

prs_tab<-data.frame(FID=obj.bigSNP$fam$family.ID,IID=obj.bigSNP$fam$sample.ID,PRED=pred, TYPE="ALL")

#}

now<-Sys.time()
message('[',now,'][Message] done')

now<-Sys.time()
message('[',now,'][Message] storing beta scores and predictions to .tsv files')

fwrite(beta_tab, file=beta_tab_file,col.names=T, row.names=F, sep="\t")
fwrite(prs_tab, file=prs_tab_file, col.names=T, row.names=F, sep="\t")

now<-Sys.time()
message('[',now,'][Message] done')

if (!is.null(opt$phenotype) & (opt$model != "grid" | beta_is_precomp)) {

  now<-Sys.time()
  message('[',now,'][Message] evaluating model on provided phenotype')

  res_params<-check_predictions(pred,pheno)
  res<-res_params[[1]]
  params<-res_params[[2]]
  pheno_tab<-data.frame(score=params$score, prt=params$prt, r2=params$r2, snps=nrow(df_beta))

  now<-Sys.time()
  message('[',now,'][Message] done')

  now<-Sys.time()
  message('[',now,'][Message] storing performances to .rds and .tsv files')

  saveRDS(res, file = summary_file)
  fwrite(pheno_tab, file=summarytab_file, sep="\t", col.names=TRUE, row.names=F)

  now<-Sys.time()
  message('[',now,'][Message] done')

}