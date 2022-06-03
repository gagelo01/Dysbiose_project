#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)

gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/")]
vec_id <-ao_small[pmid == 30224653, id]

##
traduction = fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]


for(i in 1:length(vec_id)) {

df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
vcffile <- paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", vec_id[i], "/", vec_id[i], ".vcf.gz")
inst<-gwasvcf::query_gwas(vcffile, chrompos = "1:20000-1200000") %>% 
  gwasglue::gwasvcf_to_TwoSampleMR(.) %>% 
  as.data.table(.)
sd <- coloc:::sdY.est(vbeta = (inst$se.exposure)^2,
                      maf = inst[, ifelse(eaf.exposure > 0.5,1-eaf.exposure,eaf.exposure)],
                      n = inst$samplesize.exposure)

##
tsmr <- VariantAnnotation::readVcf(vcffile) %>% 
  gwasglue::gwasvcf_to_TwoSampleMR(.) %>% 
  as.data.table(.)

tsmr[,beta.exposure := beta.exposure / sd]
tsmr[,se.exposure := se.exposure / sd]
tsmr[, chr.exposure := as.integer(chr.exposure) ]

newrow <- ao_small[id == vec_id[i]]
newrow[, id := paste0("trait-13-", i) ]
newrow[, trait := gsub(" ", "_", trait)]
newrow[,unit := "SD"]
setnames(newrow, "build",  "initial_build")
newrow[,category := "Trait"]
newrow[,note:= paste0("standardised with coloc initially in mmHg with SD ", sd)]
colinclude <- colnames(df_index)
newrow <- newrow[, ..colinclude]
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
df_index <- rbind(df_index, newrow)


GagnonMR::formattovcf_createindex(all_out = tsmr,
                                  snp_col = "SNP",
                                  outcome_name = newrow[, trait],
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
                                  units = "SD",
                                  traduction = traduction,
                                  out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf",
                                  df_index = df_index,
                                  group_name = "public",
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
                                  ID = newrow$id )


fwrite(df_index, "/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
}

message("This script finished without error")
