#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)

setwd("/mnt/sda/gagelo01/Projects/small_MR_exploration/fg_BMI")

fg <- fread("Data/Raw/FG_combined_1000G_density_formatted_21-03-29.txt.gz")
setnames(fg, c("a1", "a2"), c("other_allele", "effect_allele"))

traduction = fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]


GagnonMR::formattovcf_createindex(all_out = fg,
                                  snp_col = "rsid",
                                  outcome_name = "Fasting_Glucose",
                                  beta_col = "beta",
                                  se_col = "se",
                                  pval_col = "p-value",
                                  eaf_col = NULL,
                                  effect_allele_col = "effect_allele",
                                  other_allele_col =  "other_allele",
                                  ncase_col = NULL,
                                  ncontrol_col = NULL,
                                  samplesize_col = "n",
                                  chr_col = NULL,
                                  pos_col = NULL,
                                  units = "natural logarithm transformed FI measured in pmol/L",
                                  traduction = traduction,
                                  out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf",
                                  df_index = fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt"),
                                  group_name = "public",
                                  year = 2021,
                                  author = "Lagou V",
                                  consortium = "MAGIC",
                                  sex = "Males and Females",
                                  population = "European",
                                  initial_build = "HG19/GRCh37",
                                  category = "Trait",
                                  pmid = 33402679,
                                  note = NA,
                                  should_create_id = TRUE,
                                  ID = NA )


message("This script finished without error")



