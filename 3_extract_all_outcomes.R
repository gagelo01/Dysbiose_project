#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)

setwd("/mnt/sda/gagelo01/Projects/Dysbiose_project")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/")]
ao_small[pmid == 32203549, unit := "SD"]

inst_clump <- fread( "Data/Modified/inst_clump.txt")


ID_mrbase_out <-  c("ieu-b-109", "ieu-b-110", "ieu-b-111", "ukb-b-19953","ukb-b-12141")
outcomes_mrbase <- paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", ID_mrbase_out, "/", ID_mrbase_out, ".vcf.gz")
ID_server_out <- c("trait-2-2", "trait-6-1", "trait-7-2", "trait-2-4", "trait-12-2", 
                   "dis-2-1", "dis-3-1", "dis-4-1", "dis-5-1", "dis-6-1", "dis-7-1", "dis-8-1",
                   "trait-13-1", "trait-13-2")
outcomes_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_out, "/", ID_server_out, ".vcf.gz")


options(future.globals.maxSize= 5e9)
plan(multisession, workers = 9, gc = TRUE) #I should try using multicore


# test <- gwasvcf::query_gwas(vcf = c(outcomes_mrbase, outcomes_server)[1], rsid = unique(inst_clump$SNP),
#                             proxies = "yes", bfile = ldref)
# test %>% 
#   gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>%
#   data.table::as.data.table(.)

outcome_all <- future_map(as.list(c(outcomes_mrbase, outcomes_server)), function(x, rsiid = unique(inst_clump$SNP)) {
 
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
   res <- gwasvcf::query_gwas(vcf = x, rsid = rsiid, proxies = "yes", bfile = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
                              tag_r2 = 0.8) %>% 
    gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>%
    data.table::as.data.table(.)
  return(res)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)

add_n <- c("ieu-b-109", "ieu-b-110", "ieu-b-111") 
for(i in 1:length(add_n)) {
  outcome_all[outcome %in% add_n[i], samplesize.outcome := ao_small[id %in% add_n[i],][,.(sample_size)]]
}

outcome_all[,outcome := outcome %>% ifelse(grepl("UKB-b-19953|UKB-b-12141", .), tolower(.), .)]
k<-ao_small[id %in% ID_mrbase_out, .(id, unit)]
setnames(k, "id", "trait")
theunits <- rbind(df_index[id %in% ID_server_out,.(trait, unit)], k)
setnames(theunits, "unit", "units.outcome")
mirge <- merge(outcome_all, theunits, by.x = "outcome", by.y = "trait")
mirge[outcome %in% "Stanzick_eGFR", units.outcome := "SD"]
fwrite(mirge, "Data/Modified/outcome_all")

print("this script finished without errors")