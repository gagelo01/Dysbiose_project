#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)
library(GagnonMR)
library(readxl)

setwd("/mnt/sda/gagelo01/Projects/Dysbiose_project")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
#
inst_sanna <- fread("Data/Raw/inst_sanna")
inst_sanna[, units.exposure := "SD"]

rsid_traduction <-  fread("/mnt/sda/couchr02/rsids/rsid_pos_86millions.txt")
rsid_traduction[,chr_pos := paste0(chr, "_", pos)]
all_out <- rsid_traduction[rsid %in% inst_sanna$SNP]
mirge <- merge(inst_sanna, all_out[,c(1,3,2) ], by.x = "SNP", by.y = "rsid", all.x = TRUE)
setnames(mirge, c("pos", "chr"), c("pos.exposure", "chr.exposure"))
inst_sanna <- mirge

###"framingham"
traduction = fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]
framingham_pheno <- c("indole_3_propionate", "trimethylamine_N_oxide") #carnitine is also
format_framingham <- function(framingham_pheno) {
  instrument <- read_excel("Data/Raw/mmc2.xlsx", sheet = framingham_pheno)
  setDT(instrument)
  instrument <- instrument[Chr != "X", ]
  instrument[, Chr := as.integer(Chr)]
    all_out <- merge(instrument, traduction, by.x = c("rsID"), by.y = c("rsid"), all = FALSE) 
    setnames(all_out, c("min_all", "maj_all"), c("effect_allele","other_allele"))
    all_out[,c("Chr", "PhysPos") := NULL] #I think it is GRCH36 so better to remove it
    all_out <- all_out[(effect_allele == a0 | effect_allele == a1) & (other_allele == a0 | other_allele == a1) & a0 != a1  & effect_allele != other_allele, ] #because low number removed, coded on the forward strand
    all_out[effect_allele == a0, beta := beta*-1]
    all_out[effect_allele == a0, effect_allele := a1]
    all_out[other_allele == a1, other_allele := a0] 
    all_out[, eaf := ifelse(EUR < 0.5, MAF, 1-MAF)]
    
  instrument <- all_out
  
  instrument[,Phenotype := framingham_pheno]
  instrument[,units := "SD"]
  instrument[,samplesize := 2076]
  inst <- TwoSampleMR::format_data(
    instrument,
    type = "exposure",
    phenotype_col = "Phenotype",
    snp_col = "rsID",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele", #I verified and this is the effect allele
    other_allele_col = "other_allele", #I verified and this is the reference allele
    eaf_col = "eaf",
    pval_col = "pval",
    units_col = "units",
    samplesize_col = "samplesize",
    chr_col = "chr",
    pos_col = "position",
    )
  setDT(inst)
  return(inst)
} #do not estimate eaf, cause old GWAS

list_framingham <- pmap(data.frame(framingham_pheno = framingham_pheno), format_framingham)
inst_framingham <- rbindlist(list_framingham, fill = TRUE)
inst_framingham <- inst_framingham[pval.exposure < 1*10^-5, ]

##kettunen
kettunen_pheno <- c("Ace")
format_kettunen <- function(kettunen_pheno) {
  instrument <- fread(paste0("Data/Raw/Summary_statistics_MAGNETIC_kettunen_", kettunen_pheno, ".txt.gz"))
  setnames(instrument, "p-value", "p_value")
  instrument <- instrument[p_value < 1*10^-5, ]
  instrument[,Phenotype := kettunen_pheno]
  instrument[,units := "SD"]
  
  inst_kett<- format_data(
    instrument,
    type = "exposure",
    phenotype_col = "Phenotype",
    snp_col = "ID",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "EA",
    other_allele_col = "NEA",
    pval_col = "p_value",
    units_col = "units",
    ncase_col = "ncase",
    ncontrol_col = "ncontrol",
    samplesize_col = "n_samples",
    chr_col = "chromosome",
    pos_col = "position")
  
  setDT(inst_kett)
  return(inst_kett)
}

list_kettunen <- pmap(data.frame(kettunen_pheno = kettunen_pheno), format_kettunen)
inst_kettunen <- rbindlist(list_kettunen, fill = TRUE)


###Lotta
vec_id <- df_index[pmid == "33414548" ,][ trait %in% c("Ile", "Leu", "Val") | note %in% c("Serotonin", "Kynurenine"), ]$id

vec_vcffile <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", vec_id, "/", vec_id, ".vcf.gz")

inst_lotta <- map(as.list(vec_vcffile), function(x) get_inst(vcffile = x,pval = 1e-6, clump = FALSE, r2 = 0.01, kb = 10000)) %>% 
  rbindlist(., fill = TRUE)
inst_lotta[,units.exposure := "SD"]
#######ruhlemann
ruhlemann <-read_excel("Data/Raw/41588_2020_747_MOESM3_ESM.xlsx", sheet = 4, range = "A3:CN10003")
setDT(ruhlemann)
ruhlemann <- ruhlemann[META.P < 1*10^-6, ]

####change exposure of ruhlemann
key_exposure <-read_excel("Data/Raw/41588_2020_747_MOESM3_ESM.xlsx", sheet = 1, range = "A4:C283")
setDT(key_exposure)
setnames(key_exposure, "Feature name", "feature")

key_exposure[, Taxonomy := mgsub::mgsub(Taxonomy,pattern = c("K_", "P_", "C_", "O_", "F_", "G_"),
                                        c("kingdom-", "phylum-", "class-",  "order-",   "family-",  "genus-"))]
key_exposure[, feature := sub("_", "-", feature)]
ruhlemann[, feature := sub("NB_", "", feature)]
ruhlemann[, feature := sub("_", "-", feature)]
test <- merge(ruhlemann ,key_exposure[,c("feature", "Taxonomy")], by = "feature", all.x = TRUE)
test[, Taxonomy %>% unique %>% length, by = "feature"]$V1 %>% table
test[, feature %>% unique %>% length, by = "Taxonomy"]$V1 %>% table
test$Taxonomy %>% unique(.) %>% length(.) #Here I make a decision. I want my exposure to be "feature", but i will include Taxonomy
test[ , Taxonomy := tolower(Taxonomy)]
test[ , feature_taxonomy := paste0(feature, "(", Taxonomy, ")")]
ruhlemann <- test
##remove expousre ruhlemann already evaluated in kuri
ku_bac <- read_excel("Data/Raw/41588_2020_763_MOESM3_ESM.xlsx", sheet = 3, skip =2)
setDT(ku_bac)
ku <- ku_bac[Quant==TRUE, Taxon.name]
ku <- unique(ku)
ku<-ku%>% gsub("\\.id..*", "", .)
ku<-tolower(ku)
ku <- ku %>% gsub("\\.", "-", .) %>% gsub("--", "-", .)
ruhlemann <- ruhlemann[!(Taxonomy %in% ku), ]
#add rsid and eaf
mirge <- merge(ruhlemann, traduction, by.x = c("chrom", "pos"), by.y = c("chr", "position"), all = FALSE) #genome assembly hg19 (GRCh37).
all_out <- mirge
all_out[,units := "SD"]

setnames(all_out, c("A1", "A2"), c("other_allele", "effect_allele"))
all_out <- all_out[(effect_allele == a0 | effect_allele == a1) & (other_allele == a0 | other_allele == a1) & a0 != a1  & effect_allele != other_allele, ] #because low number removed, coded on the forward strand
all_out[effect_allele == a0, META.BETA := META.BETA*-1]
all_out[effect_allele == a0, effect_allele := a1]
all_out[other_allele == a1, other_allele := a0] 
all_out[, eaf := ifelse(EUR < 0.5, META.MAF, 1-META.MAF)]


##format
ruhlemann_inst  <- TwoSampleMR::format_data(
  all_out,
  type = "exposure",
  phenotype_col = "feature_taxonomy",
  snp_col = "rsid",
  beta_col = "META.BETA",
  se_col = "META.se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele", #See Supplementary Tables 2
  other_allele_col = "other_allele", #See Supplementary Tables 2
  pval_col = "META.P",
  units_col = "units",
  samplesize_col = "META.N",
  id_col = "id",
  chr_col = "chrom",
  pos_col = "pos") #

setDT(ruhlemann_inst)

inst_ruhlemann <- ruhlemann_inst

#kurilshikov
kurilshikov <- fread("Data/Raw/MBG.allHits.p1e4.txt")
kurilshikov <- kurilshikov[P.weightedSumZ < 1e6,]

inst_kurilshikov <- TwoSampleMR::format_data(
  kurilshikov,
  type = "exposure",
  phenotype_col = "bac",
  snp_col = "rsID",
  beta_col = "beta",
  se_col = "SE",
  eaf_col = ,
  effect_allele_col = "eff.allele",
  other_allele_col = "ref.allele",
  pval_col = "P.weightedSumZ",
  samplesize_col = "N",
  chr_col = "chr",
  pos_col = "bp")

setDT(inst_kurilshikov)
inst_kurilshikov[,units.exposure := "SD"]
all_out <- merge(inst_kurilshikov, traduction, by.x = c("chr.exposure", "pos.exposure"), by.y = c("chr", "position"), all = FALSE)
setDT(all_out)                                   
all_out <- all_out[(effect_allele.exposure == a0 | effect_allele.exposure == a1) & (other_allele.exposure == a0 | other_allele.exposure == a1)
                   &  a0 != a1 & effect_allele.exposure != other_allele.exposure, ]
all_out <- all_out[chr.exposure %in% 1:22, ]
all_out[, `:=`(chr.exposure, as.integer(chr.exposure))]
all_out[effect_allele.exposure == a0, `:=`(beta.exposure, beta.exposure * -1)]
all_out[effect_allele.exposure == a0, `:=`(effect_allele.exposure, a1)]
all_out[other_allele.exposure == a1, `:=`(other_allele.exposure, a0)]
all_out[, `:=`(eaf.exposure, EUR)]

all_out[, c("rsid", "a0", "a1","EUR") := NULL]
all_out <- all_out[pval.exposure < 1e-6,]
inst_kurilshikov <- all_out

#####change exposure and id column
list_inst <- list(inst_sanna, inst_framingham, inst_kettunen, inst_lotta, inst_ruhlemann, inst_kurilshikov)
study_name <- c( "sanna", "framingham", "kettunen", "lotta", "ruhlemann", "kurilshikov")
names(list_inst) <- study_name
for(i in 1:length(list_inst)) {
  inst <- list_inst[[i]]
  inst$exposure <- paste0(study_name[i], "_", inst$exposure)
  inst$id.exposure <- inst$exposure
  list_inst[[i]] <- inst
}

dt_inst_noselect <- rbindlist(list_inst, fill = TRUE)
fwrite(dt_inst_noselect, "Data/Modified/dt_inst_noselect.txt")

######add ldl cholesterol
##add ldl of GLGC consortium
ldl <- GagnonMR::get_inst("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/ieu-a-300/ieu-a-300.vcf.gz",pval = 5e-8, r2 = 0.01) 
ldl[, exposure := "willer_LDL-cholesterol:ieu-a-300" ]
ldl[,units.exposure :=  "SD"]
ldl[, data_source.exposure := NULL]
ldl[,id.exposure := exposure]

###rbind by names
dt_inst <- rbindlist(c(list_inst, list(ldl)), use.names = TRUE, fill = TRUE)
dt_inst[is.na(units.exposure),]
dt_inst[is.na(samplesize.exposure),]
dt_inst[,c("chrompos", "CHROM", "POS", "N_ALLELES", "N_CHR", "allele1", "freq1", "allele2", "freq2") := NULL]

####rsq
dt_inst <- add_rsq(dt_inst)
setDT(dt_inst)
dt_inst[is.na(rsq.exposure),] #parfait


#add column study and microbiote
dt_inst <- separate(dt_inst, col = "exposure", into = c("study", "microbiote"), sep = "_", remove = FALSE)

fwrite(dt_inst, "Data/Modified/all_inst.txt")

all_inst <- dt_inst

all_inst_split <- split(all_inst, all_inst$exposure)

inst_clump <- map(all_inst_split, function(x) { 
rsiid <- x %>% dplyr::select(rsid=SNP, pval=pval.exposure, id = id.exposure) %>%
  ieugwasr::ld_clump(., clump_r2 = 0.01, plink_bin=genetics.binaRies::get_plink_binary(), bfile=ldref) %>%
  {.$rsid}
return(x[SNP %in% rsiid]) }) %>% rbindlist(., fill = TRUE)

#fstat
inst_clump_split <- split(inst_clump, inst_clump$exposure)
fstat_fromdat <- function(dat) {
  k <- nrow(dat)
  n <- mean(dat$samplesize.exposure)
  r2sum <- sum(dat$rsq.exposure)
  Fstat <- ((n - k - 1)/k) * (r2sum/(1 - r2sum))
  return(Fstat)
}
vec_fstat <-lapply(inst_clump_split, fstat_fromdat)
inst_clump$fstat.exposure <- vec_fstat[match(inst_clump$exposure, names(vec_fstat))]
inst_clump[is.na(fstat.exposure)|fstat.exposure<0,] #parfait

#
exposure_to_remove <- inst_clump[, mean(fstat.exposure) > 10, by = "exposure"][V1 == FALSE,]$exposure #select mean fstat > 10
inst_clump <- inst_clump[!(exposure %in% exposure_to_remove), ]
exposure_to_remove <- inst_clump[, .N >= 3, by = "exposure"][V1 == FALSE, ]$exposure
inst_clump <- inst_clump[!(exposure %in% exposure_to_remove), ]

fwrite(inst_clump, "Data/Modified/inst_clump.txt")

print("script finished without errors")
