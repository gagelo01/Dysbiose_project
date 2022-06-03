#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)
library(GagnonMR)

setwd("/mnt/sda/gagelo01/Projects/Dysbiose_project")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
inst_clump <- fread( "Data/Modified/inst_clump.txt")
outcome_all <- fread( "Data/Modified/outcome_all")

arguments <- tidyr::crossing(exposure_name = unique(inst_clump$exposure), outcome_name = unique(outcome_all$outcome))
list_arguments <- split(arguments, 1:nrow(arguments))
harm_all <- map(list_arguments, function(x) {
  harm <- harmonise_data(inst_clump[exposure == x$exposure_name,], outcome_dat = outcome_all[outcome == x$outcome_name,], 
                         action = 1) #there all coded on the forward strand is rhee coded on the forward strand I assume that yes, because build 36 is coded on the forward strand
  return(harm)
}) %>% rbindlist(., fill = TRUE)

# harm_all1 <- harmonise_data(inst_clump[!(study %in%  c("framingham", "IBD")), ], outcome_all, action = 1)
# harm_all2 <- harmonise_data(inst_clump[study %in%  c("framingham", "IBD")], outcome_all, action = 2)
# harm_all <- rbindlist(list(harm_all1, harm_all2), fill = TRUE)
setDT(harm_all)
harm_all[, exposure_outcome := paste0(exposure,"_", outcome)]

##remove every association with less than 3 genetic instruments
harm_all<-harm_all[mr_keep==TRUE,]
to_remove <- harm_all[, .N, by ="exposure_outcome"][N<3,]$exposure_outcome
harm_all <- harm_all[!(exposure_outcome %in% to_remove),]
# ####remove ruhlemann exposure that are in kurilshikov
# harm_all[, taxa := tolower(exposure)]
# harm_all[study == c("kurilshikov"),taxa := taxa %>%gsub("kurilshikov_","",.)%>% gsub("\\.id..*", "", .) %>% gsub("\\.", "-", .) %>% gsub("--", "-", .)]
# harm_all[study == c("ruhlemann"), taxa := gsub("[\\(\\)]", "", regmatches(exposure, gregexpr("\\(.*?\\)", exposure)))]
# ku<-harm_all[study == c("kurilshikov"),taxa]
# ru<-harm_all[study == c("ruhlemann"), taxa ]
# to_remove <- harm_all[study == c("ruhlemann") & taxa %in% ru[(ru %in% ku)], unique(exposure)]
# harm_all <- harm_all[!(exposure %in% to_remove),]
###primary analysis
harm_all_split <- split(harm_all, harm_all$exposure_outcome)
perform_primary_analysis <- function(dat) {
  if (nrow(dat) > 3) {
    primary_analysis <- "mr_ivw" #multiplicative random effects IVW estimate. The standard error is corrected for under dispersion
  }
  else {
    primary_analysis <- "mr_ivw" #"mr_ivw_fe"
  }
  res <- TwoSampleMR::mr(dat, method_list = c(primary_analysis))
  return(res)
}
primary_df<-map(harm_all_split, perform_primary_analysis) %>% rbindlist

setDT(primary_df)
willer_primary <- primary_df[exposure == "willer_LDL-cholesterol:ieu-a-300"]
primary_df<- primary_df[exposure != "willer_LDL-cholesterol:ieu-a-300"]


#adjust for multiple  testing
primary_df$fdr <- p.adjust(primary_df$pval, method = "fdr")#OK so nothing comes out
primary_df[order(fdr)]

primary_df[fdr < 0.05,]
primary_df[,min(fdr)]

##steiger filtering
harm_all <- steiger_filtering(harm_all)
setDT(harm_all)
harm_all_steiger <- harm_all[harm_all$steiger_dir,]
setDT(harm_all_steiger)
harm_all_steiger_split <- split(harm_all_steiger, harm_all_steiger$exposure_outcome)

primary_df_steiger <-map(harm_all_steiger_split, perform_primary_analysis) %>% rbindlist

#remove SNPs in the ABO, APOE, HLA-A gene region
to_exclude = c("APOE", "ABO", "HLA-A")
window = 2e+06
gencode <- fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt")
list <- vector(mode = "list", length = length(to_exclude))
for (i in 1:length(to_exclude)) {
  bon <- gencode[gene_name == to_exclude[i], ]
  list[[i]] <- data.frame(chr = bon[1, ]$chr, start = min(bon$start) - 
                            window/2, end = max(bon$end) + window/2, gene_name = bon[1, 
                            ]$gene_name)
}
region_df <- rbindlist(list)
harm_all[, is_in_pleiotropic_region := FALSE]
for (i in 1:nrow(region_df)) {
  harm_all[(chr.exposure == region_df[i, ]$chr) & 
           (pos.exposure >= region_df[i, ]$start) & (pos.exposure <= 
                                                       region_df[i, ]$end), is_in_pleiotropic_region := TRUE]
}
harm_all_notpleio <- harm_all[is_in_pleiotropic_region == FALSE, ]
setDT(harm_all_notpleio)
harm_all_notpleio_split <- split(harm_all_notpleio, harm_all_notpleio$exposure_outcome)

primary_df_notpleio <- map(harm_all_notpleio_split, perform_primary_analysis) %>% rbindlist

fwrite(primary_df_steiger, "Data/Modified/primary_df_steiger.txt")
fwrite(primary_df_notpleio, "Data/Modified/primary_df_notpleio.txt")
fwrite(harm_all, "Data/Modified/harm_all.txt")
###robust MR methods
setDT(harm_all)
primary_df[, exposure_outcome := paste0(exposure,"_", outcome)]
exposure_sign <- primary_df[pval < 0.05, ]$exposure_outcome
exposure_sign <- exposure_sign[!grepl("willer_LDL", exposure_sign)]
harm_all[, exposure_outcome := paste0(exposure,"_", outcome)]
harm_all[,id.exposure := exposure]
harm_all[,id.outcome := outcome]
list_harm_sign <- split(harm_all[exposure_outcome %in% exposure_sign], harm_all[exposure_outcome %in% exposure_sign]$exposure_outcome)
list_sensitivity <- lapply(list_harm_sign, function(x) GagnonMR::all_mr_methods(dat = x, Primary_mr = "random_underdispersion"))

list_sensitivity <- map(list_sensitivity, function(x) x[!(method %in% 
                                                            c("Weighted mode", "MR Egger", "Robust adjusted profile score (RAPS)")),])

egger_intercept_test<- map_dbl(list_harm_sign, function(x) mr_pleiotropy_test(x) %>%
                                 as.data.table(.) %>%
                                 .$pval) < 0.05 #FALSE == intercept (average pleitropy) not significantly different from zero 
df_egger <- data.frame(exposure_outcome = names(egger_intercept_test),
                       intercept_test = unlist(egger_intercept_test))

map(list_sensitivity, function(x) x[,exposure_outcome := paste0(exposure, "_", outcome)])
list_sensitivity <- map(list_sensitivity, function(x) {merge(x,df_egger, by = "exposure_outcome")})

veclog <- vector(mode = "logical", length = length(list_sensitivity))
toinclude <- FALSE #Do we want to include those that could not be evaluated#
for(i in 1:length(list_sensitivity)) {
  dat <- list_sensitivity[[i]]
  outi <- dat[method %in%  c("MR-PRESSO (Outlier-corrected)", "IVW radial"), ]
  outlier_robust_val <- if(nrow(outi) == 0) {toinclude} else{outi[,pval < 0.05]}
  conta <- dat[method == "Contamination mixture", ]
  conta_val <- if(nrow(conta) == 0) {toinclude} else{conta[,pval < 0.05]}
  
  stat <- outlier_robust_val &
    dat[method ==  "Weighted median", pval < 0.05] &
    dat[1,intercept_test == FALSE] &
    conta_val
  
  veclog[i]<- stat
  names(veclog)[i]<-dat$exposure_outcome %>% unique
}
 
saveRDS(veclog, "Data/Modified/Sensitivity/veclog")
saveRDS(list_sensitivity, "Data/Modified/Sensitivity/list_sensitivity")
#power analysis
list_harm<- split(harm_all, harm_all$exposure_outcome)

vec_power <- map_dbl(list_harm, function(x) power.calculator_fromdat_vishner(x, effect = 0.1 ) )

primary_df[, exposure_outcome := paste0(exposure, "_", outcome)]
primary_df$power <- vec_power[match(primary_df$exposure_outcome, names(vec_power))]
primary_df <- rbindlist(list(primary_df, willer_primary), fill = TRUE)
fwrite(primary_df, "/home/gagelo01/workspace/Projects/Dysbiose_project/Data/Modified/Primary/primary_df")

#MVMR correcting for BMI and alcohol intake frequency

harm_all_sen <- harm_all[exposure_outcome %in% names(veclog[veclog == TRUE]),] #half (3 out of 6) sensitivity analyses with p-value <0.05
id_vec <- c("ukb-b-19953", "ukb-b-5779") #BMI alcohol in UKB
list_harm_sen <- split(harm_all_sen,  harm_all_sen$exposure_outcome)

sensitivity_mvmr <- function(dat, id) {
  the_clump <- inst_clump[exposure == unique(dat$exposure),]
  confound <- gwasvcf::query_gwas(vcf = paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", id, "/", id,".vcf.gz"), 
                                  rsid = dat$SNP,  proxies = "yes", bfile = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs") %>%
    gwasglue::gwasvcf_to_TwoSampleMR(., "exposure")
  
  inst_mvmr <- rbind(the_clump, confound, fill = TRUE)
  df_instrument <-  GagnonMR::prepare_for_mvmr(exposure_dat = inst_mvmr, d1 = inst_mvmr, clump_r2 = 0.01, clump_kb = 10000,
                                               harmonise_strictness = 1, pval_threshold = 1, clump_exp = NULL, should_clump = TRUE)
  
  outcome_names <- dat[,c("SNP", colnames(dat)[grepl("outcome", colnames(dat))])]
  outcome <- dat[,.SD, .SDcols = outcome_names]
  exposure_outcome_harmonized <- mv_harmonise_data(exposure_dat = df_instrument,
                                                   outcome_dat = outcome,
                                                   harmonise_strictness = 1)
  mvmr_results <- GagnonMR::mv_multiple_MendelianRandomization(exposure_outcome_harmonized = exposure_outcome_harmonized, 
                                                               only_IVW = TRUE)
  
  mvmr_results
  res<-mvmr_results[exposure == unique(dat$exposure),]
  return(res)
}

list_counfound_BMI <- map(list_harm_sen, function(x) sensitivity_mvmr(dat = x, id = id_vec[1]))
list_counfound_alcohol <- map(list_harm_sen, function(x) sensitivity_mvmr(dat = x, id = id_vec[2]))


no_counfound <- lapply(list_sensitivity[veclog], function(x) x[1,]) %>% rbindlist(.)
withbmi <- rbindlist(list_counfound_BMI)
withalcohol <- rbindlist(list_counfound_alcohol)

no_counfound$MVMR <- "no counfounder"
no_counfound$method <- "IVW"
# setnames(no_counfound, c("exposure", "outcome", "b", "se", "lci","uci", "pval"),
         # c("Exposure", "Outcome",  "Estimate", "StdError", "CILower",  "CIUpper",  "Pvalue"))
withbmi$MVMR <- "with BMI"
withalcohol$MVMR <- "with alcohol intake frequency"

MVMR<-rbindlist(list(no_counfound[,c("method", "exposure", "outcome",  "b", "se", "lci",  "uci",  "pval", "MVMR", "nsnp")],
                     withbmi, withalcohol), fill = TRUE)


MVMR<-MVMR[order(exposure, outcome),]

fwrite(MVMR, "Data/Modified/Sensitivity/MVMR")



message("This script finished without errors")

# 
# ########################## descriptive statistics   ###################################
# harm_all <- fread( "Data/Modified/harm_all.txt")
# harm_all <- harm_all[exposure != "willer_LDL-cholesterol:ieu-a-300",]
# primary_df <- fread( "/home/gagelo01/workspace/Projects/Dysbiose_project/Data/Modified/Primary/primary_df")
# primary_df <- primary_df[exposure != "willer_LDL-cholesterol:ieu-a-300"  ]
# primary_df[, exposure_outcome := paste0(exposure, "_", outcome)]
# harm_all[, .N, by = "exposure_outcome"][order(N)]
# list_sensitivity <- readRDS( "Data/Modified/Sensitivity/list_sensitivity")
# vec5<- sapply(list_sensitivity, function(x) sum(x[, pval] < 0.05 )) # est-ce qu'il y a
# #included exposure
# harm_all$exposure %>% unique %>% length
# harm_all$outcome %>% unique %>% length
# harm_all$exposure_outcome %>% unique %>% length
# #included microbial metabolites
# harm_all[study %in% c("sanna", "kettunen", "framingham"), length(unique(exposure))]
# harm_all[study %in% c("sanna", "kettunen", "framingham"), length(unique(exposure_outcome))]
# primary_df[exposure_outcome %in% harm_all[study %in% c("sanna", "kettunen", "framingham"), ]$exposure_outcome, ][,sum(pval < 0.05)]
# 
# primary_df_full[exposure == "willer_LDL-cholesterol:ieu-a-300" & outcome %in% c("van_der_Harst_CAD", "Timmers_parental_lifespan") , ][,.(outcome,b, b-1.96*se, b+1.96*se, pval)]
# 
# #included microbial taxa abundance
# harm_all[study %in% c("kurilshikov", "ruhlemann"), length(unique(exposure))]
# harm_all[study %in% c("kurilshikov", "ruhlemann"), length(unique(exposure_outcome))]
# primary_df[exposure_outcome %in% harm_all[study %in% c("kurilshikov", "ruhlemann"), ]$exposure_outcome, ][,sum(pval < 0.05)]
# 
# #general
# primary_df[,length(unique(exposure_outcome))]
# primary_df[,sum(pval < 0.05)]
# primary_df[which.min(fdr),]
# 
# #power
# mean(primary_df$power)
# sd(primary_df$power)
# 
# #sensitivity
# list_sensitivity <- readRDS( "Data/Modified/Sensitivity/list_sensitivity")
# length(list_sensitivity)
# sum(vec5>=4)
# list_sensitivity[vec5 >= 4]
# list_sensitivity$kettunen_Ile_NAFLD[method == "Inverse variance weighted", c(exp(b), exp(lci), exp(uci))]
# list_sensitivity$kettunen_Ile_NAFLD
# 
# ##mvmr
# MVMR<-fread( "Data/Modified/Sensitivity/MVMR")
# MVMR[, max(Pvalue), by = c("Exposure", "Outcome")]
# 
# #positive findings
# ok <- MVMR[, max(Pvalue), by = c("Exposure", "Outcome")]
# dt_sensitivity <- rbindlist(list_sensitivity, fill = TRUE)
# dt_sensitivity[exposure == ok[V1< 0.05,]$Exposure[1] & outcome == ok[V1< 0.05,]$Outcome[1],][method == "Inverse variance weighted",c( exposure, exp(b), exp(lci), exp(uci))]
# dt_sensitivity[exposure == ok[V1< 0.05,]$Exposure[2] & outcome == ok[V1< 0.05,]$Outcome[2],][method == "Inverse variance weighted",c( exposure, exp(b), exp(lci), exp(uci))]
# dt_sensitivity[exposure == ok[V1< 0.05,]$Exposure[3] & outcome == ok[V1< 0.05,]$Outcome[3],][method == "Inverse variance weighted",c( exposure, exp(b), exp(lci), exp(uci))]
# 
# 
# dt_sensitivity[exposure == ok[V1> 0.05 & V1 < 0.1,]$Exposure[1] & outcome == ok[V1> 0.05 & V1 < 0.1,]$Outcome[1],][method == "Inverse variance weighted",c( exposure, exp(b), exp(lci), exp(uci))]
# dt_sensitivity[exposure == ok[V1> 0.05 & V1 < 0.1,]$Exposure[2] & outcome == ok[V1> 0.05 & V1 < 0.1,]$Outcome[2],][method == "Inverse variance weighted",c( exposure, b, lci, uci)]
# 
