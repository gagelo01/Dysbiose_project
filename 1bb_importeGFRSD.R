#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)


gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
df_index[trait == "Stanzick_eGFR", ]



inst <- get_inst(vcffile = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-12-1/trait-12-1.vcf.gz")
sd <- coloc:::sdY.est(vbeta = (inst$se.exposure)^2,
                      maf = inst[, ifelse(eaf.exposure > 0.5,1-eaf.exposure,eaf.exposure)],
                      n = inst$samplesize.exposure)

sd

newrow <- df_index[trait == "Stanzick_eGFR", ]
newrow[, id := "trait-12-2"]
newrow[, note := "standardised with coloc"]
df_index <- rbind(df_index, newrow)

stanzick <- fread("Data/Raw/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz")

stanzick[, Effect := Effect / sd]
stanzick[, StdErr := StdErr / sd]

traduction <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]

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
                                  consortium = "CKDGEN",
                                  sex = "Males and Females",
                                  population = "European",
                                  initial_build = "HG19/GRCh37",
                                  category = "Trait",
                                  pmid = 31152163,
                                  note = "Standardised with coloc",
                                  should_create_id = FALSE,
                                  ID = "trait-12-2")
fwrite(df_index,"/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

######

inst <- get_inst(vcffile = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-11-1/trait-11-1.vcf.gz")
sd <- coloc:::sdY.est(vbeta = (inst$se.exposure)^2,
                      maf = inst[, ifelse(eaf.exposure > 0.5,1-eaf.exposure,eaf.exposure)],
                      n = inst$samplesize.exposure)

sd

newrow <- df_index[trait == "wuttke_eGFR", ]
newrow[, id := "trait-11-2"]
newrow[, note := "standardised with coloc"]
df_index <- rbind(df_index, newrow)

wuttke <- fread("Data/Raw/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt.gz")
wuttke[, Effect := Effect / sd]
wuttke[, StdErr := StdErr / sd]


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
                                  year = 2021,
                                  author = "Wuttke Matthias",
                                  consortium = "CKDGEN UKBbiobank",
                                  sex = "Males and Females",
                                  population = "European",
                                  initial_build = "HG19/GRCh37",
                                  category = "Trait",
                                  pmid = 34272381,
                                  note = "Standardise with coloc",
                                  should_create_id = FALSE,
                                  ID = "trait-11-2")

fwrite(df_index,"/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

message("This script finished without errors")