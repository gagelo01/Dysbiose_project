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
ao <-  fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao[grepl("glomerular", tolower(trait))][population == "European"]
setwd("/mnt/sda/gagelo01/Projects/Dysbiose_project")


wuttke <- fread("Data/Raw/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz")
stanzick <- fread("Data/Raw/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz")

traduction <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]


GagnonMR::formattovcf_createindex(all_out = wuttke,
                                  snp_col = "RSID",
                                  outcome_name = "wuttke_eGFR",
                                  beta_col = "Effect",
                                  se_col = "StdErr",
                                  pval_col = "P-value",
                                  eaf_col = "Freq1",
                                  effect_allele_col = "Allele1",
                                  other_allele_col =  "Allele2",
                                  ncase_col = NULL,
                                  ncontrol_col = NULL,
                                  samplesize_col = "n_total_sum",
                                  chr_col = "Chr",
                                  pos_col = "Pos_b37",
                                  units = "log(eGFR)",
                                  traduction = traduction,
                                  out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf",
                                  df_index = df_index,
                                  group_name = "public",
                                  year = 2019,
                                  author = "Wuttke Matthias",
                                  consortium = "CKDGEN",
                                  sex = "Males and Females",
                                  population = "European",
                                  initial_build = "HG19/GRCh37",
                                  category = "Trait",
                                  pmid = 31152163,
                                  note = "",
                                  should_create_id = TRUE,
                                  ID = NA)

df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")


GagnonMR::formattovcf_createindex(all_out = stanzick,
                                  snp_col = "RSID",
                                  outcome_name = "Stanzick_eGFR",
                                  beta_col = "Effect",
                                  se_col = "StdErr",
                                  pval_col = "P.value",
                                  eaf_col = "Freq1",
                                  effect_allele_col = "Allele1",
                                  other_allele_col =  "Allele2",
                                  ncase_col = NULL,
                                  ncontrol_col = NULL,
                                  samplesize_col = "n",
                                  chr_col = "chr",
                                  pos_col = "pos",
                                  units = "log(eGFR)",
                                  traduction = traduction,
                                  out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf",
                                  df_index = df_index,
                                  group_name = "public",
                                  year = 2021,
                                  author = "Stanzick Kira",
                                  consortium = "CKDGEN UKBbiobank",
                                  sex = "Males and Females",
                                  population = "European",
                                  initial_build = "HG19/GRCh37",
                                  category = "Trait",
                                  pmid = 34272381,
                                  note = "",
                                  should_create_id = TRUE,
                                  ID = NA)

message("This script finished without errors")