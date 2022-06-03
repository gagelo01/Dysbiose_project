#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)
library(GagnonMR)



########################## descriptive statistics   ###################################
inst_clump <- fread( "Data/Modified/inst_clump.txt")
outcome_all <- fread( "Data/Modified/outcome_all")
harm_all <- fread( "Data/Modified/harm_all.txt")
harm_all <- harm_all[exposure != "willer_LDL-cholesterol:ieu-a-300",]
primary_df_full <- fread( "/home/gagelo01/workspace/Projects/Dysbiose_project/Data/Modified/Primary/primary_df")
primary_df_full <- merge(primary_df_full, distinct(harm_all[,.(exposure, study)]), by = "exposure", all.x = TRUE)
primary_df_full[, lci := b - se*1.96]
primary_df_full[, uci := b + se*1.96]
primary_df <- primary_df_full[exposure != "willer_LDL-cholesterol:ieu-a-300" ]
primary_df[, exposure_outcome := paste0(exposure, "_", outcome)]


harm_all[, .N, by = "exposure_outcome"][order(N)]
list_sensitivity <- readRDS( "Data/Modified/Sensitivity/list_sensitivity")
veclog <- readRDS( "Data/Modified/Sensitivity/veclog")
MVMR_dat <- fread("/home/gagelo01/workspace/Projects/Dysbiose_project/Data/Modified/Sensitivity/MVMR")

#####
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ID_mrbase_out <-  c("ieu-b-109", "ieu-b-110", "ieu-b-111", "ukb-b-19953","ukb-b-12141")
ID_server_out <- c("trait-2-2", "trait-6-1", "trait-7-2", "trait-2-4", "trait-12-2", 
                   "dis-2-1", "dis-3-1", "dis-4-1", "dis-5-1", "dis-6-1", "dis-7-1", "dis-8-1",
                   "trait-13-1", "trait-13-2", "trait-12-2")


k <- df_index[id %in% ID_server_out,]
m <- ao[id %in% ID_mrbase_out, ]

dataset <-rbind(k,m,fill = TRUE )
dataset <- dataset[,.(trait,group_name, year, author,  consortium, sex, population, unit, nsnp, sample_size, ncase, ncontrol, pmid, note)]
######
return_format_data<-function(data) {
  return(data[, paste0(round(exp(b), digits = 2), " 95% CI=", round(exp(lci), digits = 2), "-",  round(exp(uci), digits = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}
return_format_data_noexp <-function(data) {
  return(data[, paste0(round(b, digits = 2), " 95% CI=", round(lci, digits = 2), "-",  round(uci, digits = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}
#included exposure
harm_all$exposure %>% unique %>% length
harm_all$outcome %>% unique %>% length
harm_all$exposure_outcome %>% unique %>% length


#Abstract
harm_all[study %in% c("sanna", "kettunen", "framingham", "lotta"), length(unique(exposure))]
harm_all[study %in% c("kurilshikov", "ruhlemann"),length(unique(exposure))]
sum(veclog)
harm_all[study != "willer", ]$exposure_outcome %>% unique %>% length
dt_sen <-rbindlist(list_sensitivity[veclog])
dt_sen[method == "Inverse variance weighted", ][abs(b)<0.1, .N]

##INtro

#Results
#para 1
harm_all[, length(unique(exposure))]
dataset[population != "European", ] 
harm_all[study %in% c("sanna", "kettunen", "framingham", "lotta"), length(unique(exposure))]
harm_all[study %in% c("kurilshikov", "ruhlemann"),length(unique(exposure))]
harm_all[,length(unique(exposure))]

#para 2
harm_all[study %in% c("sanna", "kettunen", "framingham", "lotta"), length(unique(exposure_outcome))]
primary_df[study %in% c("sanna", "kettunen", "framingham", "lotta"), mean(abs(b))]
primary_df[study %in% c("sanna", "kettunen", "framingham", "lotta") & pval < 0.05, .N]
primary_df[study %in% c("sanna", "kettunen", "framingham", "lotta") & fdr < 0.05, .N]
primary_df_full[exposure == "willer_LDL-cholesterol:ieu-a-300" & 
                  outcome %in% c("van_der_Harst_CAD", "Deelen_longevity")] %>% return_format_data(.)

#para 3
harm_all[study %in% c("kurilshikov", "ruhlemann"), length(unique(exposure))]
primary_df[study %in% c("kurilshikov", "ruhlemann"), length(unique(exposure_outcome))]
primary_df[study %in% c("kurilshikov", "ruhlemann"), mean(abs(b))]
primary_df[study %in% c("kurilshikov", "ruhlemann") & pval < 0.05, .N]
primary_df[study %in% c("kurilshikov", "ruhlemann") & fdr < 0.05, .N]

#Para 4
primary_df[pval < 0.05,.N]
length(list_sensitivity)
sum(veclog)
harm_all[exposure_outcome  %in% names(list_sensitivity)[veclog], all(steiger_dir)]
harm_all[exposure_outcome  %in% names(list_sensitivity)[veclog], any(is_in_pleiotropic_region)]

dtsen <- rbindlist(list_sensitivity[veclog])
dtsen[method == "Inverse variance weighted", ][abs(b)<0.1, .N]

dtsen[method == "Inverse variance weighted" & abs(b)>0.2,] %>%
  return_format_data(.)


#para7
{MVMR_dat[,b := as.numeric(b)]
MVMR_dat[, exposure_outcome := paste0(exposure, "_", outcome)]
MVMR_dat_split <- split(MVMR_dat, MVMR_dat$exposure_outcome)
list_dat <- vector(mode = "list", length(MVMR_dat_split))
for(i in 1:length(MVMR_dat_split)) {
dat <- MVMR_dat_split[[i]]
reference <- dat[MVMR == "no counfounder", b]
dat[, comparison := apply(.SD, 1, function(x) as.numeric(x[names(x)=="b"])/reference), ]
list_dat[[i]] <- dat
}

datfull <-rbindlist(list_dat)}

datfull[comparison < 0.6 & !(outcome == "ukb-b-19953" & MVMR == "with BMI"),][,.(exposure,outcome,b,comparison, MVMR)]
datfull[comparison < 0.6 & !(outcome == "ukb-b-19953" & MVMR == "with BMI"),][exposure == "lotta_Serotonin" & outcome == "van_der_Harst_CAD"] %>% 
  return_format_data(.)
datfull[comparison < 0.6 & !(outcome == "ukb-b-19953" & MVMR == "with BMI"),][exposure == "kurilshikov_order.Lactobacillales.id.1800"] %>% 
  return_format_data(.)


######Discussion
primary_df[fdr<0.05,]
dtsen[method == "Inverse variance weighted" & abs(b)>0.2,]
dtsen[method == "Inverse variance weighted" & abs(b)>0.1,]

#comparison with other studies
primary_df[exposure %in% "framingham_trimethylamine_N_oxide" & outcome %in% "van_der_Harst_CAD",] %>% 
  return_format_data(.)

#Method
primary_df[,length(unique(exposure_outcome))]
suppose <- length(unique(harm_all$outcome)) * harm_all[study != "willer", length(unique(exposure))]

suppose - primary_df[,length(unique(exposure_outcome))]  
