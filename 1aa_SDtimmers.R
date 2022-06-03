#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)


gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
df_index[id == "trait-7-1", ]

############onlt change this section
ancient_id<- "trait-7-1"
new_id<-"trait-7-2"
###################

inst <- get_inst(vcffile = paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ancient_id, "/", ancient_id, ".vcf.gz"))
sd <- coloc:::sdY.est(vbeta = (inst$se.exposure)^2,
                      maf = inst[, ifelse(eaf.exposure > 0.5,1-eaf.exposure,eaf.exposure)],
                      n = inst$samplesize.exposure)

sd

newrow <- df_index[id == ancient_id, ]
newrow[, id := new_id]
newrow[, note := "standardised with coloc"]
newrow[, unit := "SD"]
df_index <- rbind(df_index, newrow)

vcf <- VariantAnnotation::readVcf(paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ancient_id, "/", ancient_id, ".vcf.gz"))
tsmr <- vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(.) %>% as.data.table(.)
tsmr[, beta.exposure := beta.exposure / sd]
tsmr[, se.exposure := se.exposure / sd]

GagnonMR::formattovcf_createindex(all_out = tsmr,
                                  snp_col = "SNP",
                                  outcome_name = newrow[,trait],
                                  beta_col = "beta.exposure",
                                  se_col = "se.exposure",
                                  pval_col = "pval.exposure",
                                  eaf_col = "eaf.exposure",
                                  effect_allele_col = "effect_allele.exposure",
                                  other_allele_col =  "other_allele.exposure",
                                  ncase_col = NULL,
                                  ncontrol_col = NULL,
                                  samplesize_col = "samplesize.exposure",
                                  chr_col = "chr.exposure",
                                  pos_col = "pos.exposure",
                                  units = newrow$unit,
                                  traduction = NULL,
                                  out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf",
                                  df_index = df_index,
                                  group_name = newrow$group_name,
                                  year = newrow$year,
                                  author = newrow$author,
                                  consortium = newrow$consortium,
                                  sex = newrow$sex,
                                  population = newrow$population,
                                  initial_build = newrow$initial_build,
                                  category = newrow$category,
                                  pmid = newrow$pmid,
                                  note = newrow$note,
                                  should_create_id = FALSE,
                                  ID = new_id)
fwrite(df_index,"/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
message("this script finished without errors")
