#!/usr/bin/env Rscript
library(data.table)
library(readxl)
library(TwoSampleMR)
library(tidyverse)
library(writexl)
library("xlsx")
library(GagnonMR)

setwd("/home/gagelo01/workspace/Projects/Dysbiose_project")
harm_all <- fread( "Data/Modified/harm_all.txt")
dataset_info<- read_excel("Data/Modified/informationondataused.xlsx")
expinf <- readxl::read_excel(path = "Data/Modified/exposures_cohort_info.xlsx") 
setDT(expinf)
########Supplementary table titles and description
dt_title <- data.table(title = paste0("Supplementary Table ", 1:7),
                       caption = c("Description of the datasets used.",
                                   "Description of the metabolites selection rational.",
                                   "Summary of the instrument selection criteria for each data source. LD R2, pvalue treshold and Fstatistics.",
                                   "Harmonised dataset for each exposure outcome.",
                                   "Primary MR results and associated statistics.",
                                    "Robust MR results and other sensitivity analyses results.",
                                   "Multivariable MR results and associated statistics."))
########Supplementary table with information
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao<-fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small<-ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/")]
ID_server_out <- c("trait-2-2", "trait-6-1", "trait-7-2", "trait-2-4", "trait-12-2", 
                   "dis-2-1", "dis-3-1", "dis-4-1", "dis-5-1", "dis-6-1", "dis-7-1", "dis-8-1",
                   "trait-13-1", "trait-13-2")
ID_mrbase_out <-  c("ieu-b-109", "ieu-b-110", "ieu-b-111", "ukb-b-19953","ukb-b-12141")

outcome_table <- rbind(ao[id %in% ID_mrbase_out, ], df_index[id %in% ID_server_out], fill = TRUE)
outcome_table[, author == "Richardson, Tom", unit := "SD"]
dis <- c("ukb-b-12141", "dis-2-1", "dis-3-1", "dis-4-1", "dis-5-1", "dis-6-1", "dis-7-1", "dis-8-1")
longlife <- c("trait-6-1", "trait-7-2") 
outcome_table[, category := id %>% ifelse(. %in% dis, "Disease", .) %>% 
                ifelse(. %in% longlife, "Mortality", .) %>%
                ifelse(!(. %in% c("Disease", "Mortality")), "Metabolic risk factor", .)]
outcome_table[,variable := "outcome"]
outcome_table <- outcome_table[,.(variable,category, trait, group_name, year, author, consortium, 
                                  sex, population, unit, nsnp,
                                  sample_size, pmid, note, ncase, ncontrol)]
outcome_table <- outcome_table[order(category, trait)] 

setDT(expinf)
expinf[trait == "Microbial relative abundance", trait := paste0(trait, " (", tolower(author), ")")]
expinf[, c("sex", "unit") := .("Males and Females", "SD")]
expinf[, category := ifelse(author %in% c("Kurilshikov", "Ruhlemann"), "Taxa abundance", "metabolites")]
lotta_tobind <- df_index[pmid == "33414548",][1,]
lotta_tobind[, trait := c("serotonin, leucine, isoleucine, valine, kynurenine")]
lotta_tobind[,category := "metabolites"]
lotta_tobind[, note := NULL]
lotta_tobind <- lotta_tobind[,.(trait, group_name, year, author, consortium, 
                sex, population, unit, nsnp, category,
                sample_size, pmid, ncase, ncontrol)]
exp <- rbind(lotta_tobind, expinf, fill = TRUE)
exp[,variable := "exposure"]

suptab1 <- rbind(outcome_table, exp, fill = TRUE)
suptab1[, c("group_name", "id", "initial_build", "sd", "nsnp") := NULL]
dattrait <- data.frame(trait = suptab1[,trait],
           url = c("https://datashare.is.ed.ac.uk/handle/10283/3203",
                   "https://ctg.cncr.nl/software/summary_statistics",
                   "http://diagram-consortium.org/downloads.html",
                   "https://www.megastroke.org/download.html",
                   "https://www.ebi.ac.uk/gwas/publications/34841290",
                   "https://gwas.mrcieu.ac.uk/files/ukb-b-12141/ukb-b-12141.vcf.gz",
                   "http://ckdgen.imbi.uni-freiburg.de/",
                   "https://data.mendeley.com/datasets/gbbsrpx6bs/1",
                   "https://gwas.mrcieu.ac.uk/files/ukb-b-19953/ukb-b-19953.vcf.gz",
                   "https://magicinvestigators.org/downloads/",
                   "https://magicinvestigators.org/downloads/",
                   "https://www.ebi.ac.uk/gwas/publications/32203549",
                   "https://www.ebi.ac.uk/gwas/publications/32203549",
                   "https://ckdgen.imbi.uni-freiburg.de/",
                   "https://www.ebi.ac.uk/gwas/publications/30224653",
                   "https://www.ebi.ac.uk/gwas/publications/30224653",
                   "https://www.ebi.ac.uk/gwas/publications/32203549",
                   "https://www.longevitygenomics.org/downloads",
                   "https://datashare.ed.ac.uk/handle/10283/3209",
                   "https://omicscience.org/apps/crossplatform/",
                   "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0350-x/MediaObjects/41588_2019_350_MOESM1_ESM.pdf",
                   "https://www.cell.com/action/showFullTableHTML?isHtml=true&tableId=tbl2&pii=S1550-4131%2813%2900257-X",
                   "http://www.computationalmedicine.fi/data#NMR_GWAS",
                   "https://www.nature.com/articles/s41588-020-00747-1",
                   "https://mibiogen.gcc.rug.nl/"))
                   
suptab1 <- cbind(suptab1, url = dattrait$url)
suptab1[, note := gsub("standardised with coloc", "standardised with coloc:::sdY.est function", note)]    

suptab1 <- merge(suptab1, dataset_info, by= "trait")
#Supplementary Table 2
suptab2<- readxl::read_excel(path = "Data/Modified/tables_metabolites.xlsx") %>% 
  as.data.table

#Supplementary Table 3
suptab3<- readxl::read_excel(path = "Data/Modified/instrumentselection.xlsx") %>% 
  as.data.table
setnames(suptab3, "min F-statistics", "F.statistics")
suptab3[, F.statistics := 10]
#Supplementary Table 4
harm_all <- fread( "Data/Modified/harm_all_clean.txt")
harm_all[,c("id.exposure", "id.outcome", "outcome", "exposure") := NULL]
setnames(harm_all, c("exposure_clean", "outcome_clean", "Category"), c("exposure", "outcome", "exposure_category"))
suptab4 <- harm_all
#Supplementary Table 5
primary_df <- fread( "Data/Modified/primary_df_clean.txt")
harm_all <- fread( "Data/Modified/harm_all_clean.txt")
harm_all <- distinct(harm_all[, .(exposure, outcome, fstat.exposure)])
harm_all[ , fstat.exposure := mean(fstat.exposure), by = c("exposure", "outcome")]
harm_all
primary_df <- merge(
  primary_df,
  distinct(harm_all),
  by = c("exposure", "outcome"))
primary_df[,c("id.exposure", "id.outcome", "outcome", "exposure") := NULL]
setnames(primary_df, c("exposure_clean", "outcome_clean", "Category"), c("exposure", "outcome", "exposure_category"))
primary_df <- primary_df[, .(exposure_category, outcome_category, exposure, outcome, method, nsnp,
b,se,pval,fdr, power, fstat.exposure) ]          
primary_df <- primary_df[order(exposure_category, outcome_category, exposure, outcome),]
suptab5 <- primary_df
suptab5[,power:=NULL]
#Supplementary table 6
list_sensitivity <- readRDS( "Data/Modified/Sensitivity/list_sensitivity")
veclog <- readRDS( "Data/Modified/Sensitivity/veclog")
suptab6 <- rbindlist(list_sensitivity)

#Supplementary Table 7
MVMR_dat <- fread("Data/Modified/MVMR_clean.txt")
MVMR_dat[,c( "outcome", "exposure") := NULL]
setnames(MVMR_dat, c("exposure_clean", "outcome_clean", "Category", "F_stastistics"),
         c("exposure", "outcome", "exposure_category", "conditional_F_stastistics"))
MVMR_dat<- MVMR_dat[, .(exposure_category, outcome_category, exposure, outcome, method, nsnp,
           b,se,lci,uci,pval,MVMR, nsnp, cochranQ,cochranQpval,conditional_F_stastistics)]
MVMR_dat <- MVMR_dat[order(exposure_category, outcome_category, exposure, outcome),]
suptab7<-MVMR_dat





writexl::write_xlsx(x = list("Tables captions and titles" = dt_title,
                             "Supplementary Table 1" = suptab1,
                             "Supplementary Table 2" = suptab2, 
                             "Supplementary Table 3" = suptab3, 
                             "Supplementary Table 4" = suptab4, 
                             "Supplementary Table 5" = suptab5, 
                             "Supplementary Table 6" = suptab6,
                             "Supplementary Table 7" = suptab7),
                    path = "Results/supplementary_tables.xlsx")











