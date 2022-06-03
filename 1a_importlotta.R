#!/usr/bin/env Rscript
#import instrument
library(data.table)
library(GagnonMR)
library(tidyverse)
library(tictoc)
library(furrr)

gwasvcf::set_bcftools()
gwasvcf::set_plink()
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
setwd("/mnt/sda/gagelo01/Projects/Dysbiose_project")

met_lab <-fread("/home/couchr02/Mendel_Commun/Christian/GWAS/proteome/Lotta/metabolites_label.txt")
colnames(met_lab) <- c("namefile", "abbreviation", "fullname", "class")

#convert label
traduction <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]


newrow<- data.table( id  = paste0("met-2-", 1:met_lab[,.N]),
                     trait = met_lab$abbreviation,
                     group_name = "Public", year = 2021,
                     author = "Lotta A, Luca",
                     consortium = "Fenland, EPIC-Norfolk (Metabolon), INTERVAL (Metabolon), INTERVAL (Nightingale), Kettunen et al. 2016 (Nightingale), Draisma et al. 2015 (Biocates p150), Shin et al.2014 (Twins-UK, Metabolon), Shin et al. 2014 (KORA, Metbolon)",
                                    sex = "Males and Females", population = "European", unit = "SD",
                     nsnp = "~10000000", sample_size = ">10000" ,initial_build = "HG19/GRCh37", category = "Metabolites",
                     pmid = 33414548, sd = 1, 
                     note = met_lab$fullname,
                     ncase = NA, ncontrol = NA)


df_index <- rbind(df_index, newrow)

lotta_ID <- df_index[pmid == "33414548", ]$id
format_wrapper <- function(lotta_id, traduction, df_index, met_lab) {
  
  namefile <- met_lab[abbreviation == df_index[id == lotta_id & pmid == "33414548", trait], namefile  ]
  data <- fread(paste0("/home/couchr02/Mendel_Commun/Christian/GWAS/proteome/Lotta/Results/", namefile, ".txt.gz"))
  
    GagnonMR::formattovcf_createindex(all_out = data,
                                    snp_col = "rsid",
                                    outcome_name = df_index[id == lotta_id, trait],
                                    beta_col = "Beta",
                                    se_col = "SE",
                                    pval_col = "Pvalue_MA",
                                    eaf_col = "Freq1_MA",
                                    effect_allele_col = "Allele1",
                                    other_allele_col =  "Allele2",
                                    ncase_col = NULL,
                                    ncontrol_col = NULL,
                                    samplesize_col = "Weight_MA",
                                    chr_col = "chr",
                                    pos_col = "pos",
                                    units = "SD",
                                    traduction = traduction,
                                    out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf",
                                    df_index = df_index,
                                    group_name = "public",
                                    year = 2021,
                                    author = "Lotta A, Luca",
                                    consortium = "Fenland, EPIC-Norfolk (Metabolon), INTERVAL (Metabolon), INTERVAL (Nightingale), Kettunen et al. 2016 (Nightingale), Draisma et al. 2015 (Biocates p150), Shin et al.2014 (Twins-UK, Metabolon), Shin et al. 2014 (KORA, Metbolon)",
                                    sex = "Males and Females",
                                    population = "Mix",
                                    initial_build = "HG19/GRCh37",
                                    category = "Metabolites",
                                    pmid = 33414548,
                                    note = df_index[id == lotta_id, note],
                                    should_create_id = FALSE,
                                    ID = lotta_id)
}


options(future.globals.maxSize= 1e10)
plan(multisession, workers = 8)

df_index_copy <- df_index
traduction_copy <- traduction
met_lab_copy <- met_lab
tic()

k <- list.files("/mnt/sda/gagelo01/Vcffile/Server_vcf")
l <- k[grepl("met-2-", k)]
n <- as.numeric(gsub("met-2-", "", l))
index <- which(!(1:173 %in% n))

future_map(lotta_ID[index], function(x) {format_wrapper(lotta_id = x, traduction =  traduction_copy, df_index =  df_index_copy, met_lab = met_lab_copy)},
           .options = furrr_options(seed = TRUE))

fwrite(df_index, "/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

message("this script finished without errors")

toc()

                     
                     
                     