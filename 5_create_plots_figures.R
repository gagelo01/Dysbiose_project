library(ckbplotr) #it has to loaded at the beginning of the session
library(data.table)
library(readxl)
library(TwoSampleMR)
library(tidyverse)
library(GagnonMR)
library(ggpubr)
library(ggforestplot)
library(ggforce)

#creating report and plots
#clean primary_df
primary_df_full <- fread("/home/gagelo01/workspace/Projects/Dysbiose_project/Data/Modified/Primary/primary_df")
inst_clump <- fread( "Data/Modified/inst_clump.txt")


cleanify <- function(primary_df_full) {
primary_df_full$exposure_clean <- as.character(primary_df_full$exposure) %>% 
  ifelse(. == "framingham_betaine", "Betaine", .) %>%
  ifelse(. == "framingham_carnitine", "Carnitine", .) %>%
  ifelse(. == "framingham_choline", "Choline", .) %>%
  ifelse(. == "framingham_indole_3_propionate", "Indole-3-propionate", .) %>%
  ifelse(. == "kettunen_Ace", "Acetate", .) %>%
  ifelse(. == "lotta_Ile", "Isoleucine", .) %>%
  ifelse(. == "lotta_Leu", "Leucine", .) %>%
  ifelse(. == "lotta_Serotonin", "Serotonin", .) %>%
  ifelse(. == "lotta_Val" , "Valine", .) %>%
  ifelse(. == "lotta_PEA", "Phenylethylamine", . ) %>%
  ifelse(. == "lotta_Glu", "Glutamate", . ) %>% 
  ifelse(. == "lotta_Asp", "Aspartate", .) %>%
  ifelse(. == "lotta_C0", "Carnitine", .) %>%
  ifelse(. == "lotta_Kynurenine", "Kynurenine", .) %>% 
  ifelse(. == "framingham_trimethylamine_N_oxide", "Trimethylamine N-oxide (TMAO)", .) %>%
  ifelse(. == "sanna_fecal_propionate_levels", "Fecal propionate", .) %>%
  ifelse(. == "sanna_PWY-5022", "Pathway PWY-5022", .) %>%
  ifelse(. == "willer_LDL-cholesterol:ieu-a-300", "LDL cholesterol", .) 


# primary_df_full$outcome_clean <- sapply(strsplit(primary_df_full$outcome, "_"), function(x) x[length(x)])
primary_df_full$outcome_clean <- primary_df_full$outcome %>% 
  ifelse(. == "osteoporosis", "Osteoporosis", .) %>%
  ifelse(. == "NAFLD", "Non-alcoholic fatty liver disease", .) %>%
  ifelse(. == "van_der_Harst_CAD" , "Coronary artery disease", .) %>%
  ifelse(. == "Malik_Stroke", "Ischemic stroke", .) %>%
  ifelse(. == "Mahajan_Type2diabetes", "Type 2 diabetes", .) %>%
  ifelse(. == "Jansen_Alzheimer", "Alzheimer's disease", .) %>%
  ifelse(. == "Timmers_parental_lifespan", "Parental lifespan", .) %>%
  ifelse(. == "Wuttke_Chronic_kidney" , "Chronic kidney disease", .) %>%
  ifelse(. == "Deelen_longevity", "Human longevity", .) %>%
  ifelse(. == "Howard_Depression", "Depression", .) %>%
  ifelse(. == "diastolic_blood_pressure", "Diastolic blood pressure", .) %>%
  ifelse(. == "Fasting_Glucose_standardised", "Fasting glucose", .) %>%
  ifelse(. == "Fasting_Insulin", "Fasting insulin", .) %>%
  ifelse(. == "ieu-b-109" , "HDL cholesterol", .) %>%
  ifelse(. == "ieu-b-110", "LDL cholesterol", .) %>%
  ifelse(. == "ieu-b-111", "Triglycerides", .) %>%
  ifelse(. == "Stanzick_eGFR", "Estimated glomerular filtration rate", .) %>%
  ifelse(. == "systolic_blood_pressure" , "Systolic blood pressure", .) %>%
  ifelse(. == "ukb-b-12141", "Osteoporosis", .) %>%
  ifelse(. == "ukb-b-19953" , "Body mass index", .)

primary_df_full[, outcome_category := ifelse(outcome_clean %in% c("Coronary artery disease", "Ischemic stroke", 
                                                                  "Chronic kidney disease",
                                                                  "Type 2 diabetes", "Non-alcoholic fatty liver disease",
                                                     "Osteoporosis", "Depression", "Alzheimer's disease", "Human longevity", 
                                                     "Parental lifespan"), "Disease", "Cardiometabolic trait")]
                                                     

primary_df_full[, outcome_clean := factor(outcome_clean,
            levels = c("Coronary artery disease", "Ischemic stroke", "Chronic kidney disease", "Type 2 diabetes", "Non-alcoholic fatty liver disease",
                      "Osteoporosis", "Depression", "Alzheimer's disease", "Human longevity", "Parental lifespan", "HDL cholesterol", 
                      "LDL cholesterol","Triglycerides",  "Systolic blood pressure", "Diastolic blood pressure","Fasting glucose", 
                      "Fasting insulin", "Estimated glomerular filtration rate", "Body mass index"))]
#Category
primary_df_full$Category <- primary_df_full$exposure %>%
  ifelse(. == "IBD_IBD", "Dysbiotic disease", . ) %>%
  ifelse(grepl("fecal_propionate_levels",.), "Fecal Metabolites", .) %>%
  ifelse(grepl("PWY-", .) & grepl("sanna", .), "Microbial Pathway", .) %>%
  ifelse(grepl("framingham", .) | grepl("kettunen|lotta", .),  "Plasma Metabolites", .) %>%
  ifelse(grepl("dutch", .) & grepl("PWY", .), "Microbial Pathway", .) %>%
  ifelse(grepl("dutch", .) & !grepl("PWY", .), "Microbe Abundance", .) %>%
  ifelse(grepl("fin_gut", .) |  grepl("ruhlemann",.) | grepl("kurilshikov",.), "Microbe Abundance", .) %>%
  ifelse(grepl("willer_LDL", .), "Positive control", .)

all(primary_df_full$Category %in% c("Dysbiotic disease", "Fecal Metabolites", "Plasma Metabolites", "Microbe Abundance", "Microbial Pathway")) #if TRUE great

setDT(primary_df_full)
primary_df_full[,Category := factor(Category, levels = c("Dysbiotic disease", "Plasma Metabolites","Fecal Metabolites", "Microbial Pathway", "Microbe Abundance", "Positive control"))]
primary_df_full<- primary_df_full[order(Category,outcome, exposure)]

primary_df_full[Category == "Microbe Abundance", ]$exposure_clean <- sapply(strsplit(primary_df_full[Category == "Microbe Abundance", ]$exposure_clean, "_"),
                                                                           function(x) paste0(x[2:length(x)], collapse = "_"))

primary_df_full[Category == "Microbe Abundance", exposure_clean := gsub("(", " (", exposure_clean, fixed = TRUE)]    
return(primary_df_full)
}
primary_df_full <- cleanify(primary_df_full = primary_df_full)
fwrite(primary_df_full, "Data/Modified/primary_df_clean.txt")
#This document builds mostly on https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html

#figure 3 -> Create a forest plot of all IVW with IBD I want it to look like in Iyas and all.
#https://neilstats.github.io/ckbplotr/articles/make_forest_plot.html#forest-plot-with-row-labels

###a better way I think. All exposures. Outcome as title. One forest plot per outcome.
# primary_df <- primary_df_full[exposure == "IBD_IBD",]
# resultsA <- data.frame(variable = paste0(primary_df$outcome_clean),
#                          estimate = primary_df$b,
#                          stderr =  primary_df$se,
#                          n = primary_df$nsnp,
#                        P_value = formatC(primary_df$pval, format = "e", digits = 1))
# 
# 
# setDT(resultsA)
# resultsA <- resultsA[order(-estimate)]
# m<-make_forest_plot(panels         = list(resultsA),
#                  col.key        = "variable",
#                  panel.headings = "Primary MR analysis for inflammatory bowel disease \n and 10 health outcomes",
#                  exponentiate   = TRUE,
#                  nullval    = 1,
#                  pointsize        = 2,
#                  col.left         = c("n"),
#                  col.left.heading = c("n SNP"),
#                  col.right = "P_value",
#                  col.right.heading = c("OR (95% CI)", "P-value"),
#                  printplot = TRUE,
#                  xlab =  "Effect of IBD on chronic diseases and longevity",
#                  colour = "blue",
#                  cicolour = "darkred")
# 
# # ggsave("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/IBD_IVW_forest_plot.png",
# #        width=436/72,height=244/72,units="in",scale=1,
# #        device = "png")
# 
# b<- m[[1]]
# a <- readRDS( file = "Analysis/LD_score/Results/Forest_plot_IBD_Rg_benoitversion.rdata")
# 
# ggarrange(a, b, 
#           labels = c("A", "B"),
#           ncol = 1, nrow = 2)
# 
# ggsave("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Figures/2_IBD_twopanel_plot.tiff",
#        width=700/72,height=730/72,units="in",scale=1,
#        device = "tiff")
###huge code for other stuff
# list_primary_forest <- vector(mode = "list", length = length(unique(primary_df_full$outcome)))
# for(i in 1:length(unique(primary_df_full$outcome))) {
#   primary_df <- primary_df_full[outcome == unique(outcome)[i]]
# 
#   resultsA <- data.frame(variable = paste0(primary_df$exposure),
#                          estimate = primary_df$b,
#                          stderr =  primary_df$se,
#                          n = primary_df$nsnp)
# 
#   if(primary_df[1,outcome] == "Timmers_parental_lifespan") {
#     exponentiate <- FALSE
#     col.right.heading <- "Beta (95% CI)"
#     xlab <- "Difference in lifespan years caused by gut dysbiosis"
#     nullval <- 0
#   } else {
#     exponentiate <- TRUE
#     col.right.heading <- "HR (95% CI)"
#     xlab <- "Hazard ratio of having a disease caused by microbiome"
#     nullval <- NULL
#   }
# 
#   list_primary_forest[[i]] <-
#     make_forest_plot(panels         = list(resultsA),
#                      col.key        = "variable",
#                      panel.headings = paste0("Primarary MR analysis for ", primary_df[1,outcome]),
#                      exponentiate   = exponentiate,
#                      nullval    = nullval,
#                      pointsize        = 2,
#                      col.left         = c("n"),
#                      col.left.space   = c(0.02),
#                      col.left.heading = c("n SNP"),
#                      col.right.heading = col.right.heading,
#                      printplot = FALSE,
#                      xlab = xlab)
# 
#   print(list_primary_forest[[i]])
#   ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Primary/Forest_plot/",
#                 i, "_", unique(primary_df_full$outcome)[i], "_forest_plot.png" ),
#          width=12,height=12,units="in",scale=1,
#          device = "png")
# 
# }
# 
# saveRDS(list_primary_forest, "/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Primary/Forest_plot/list_primary_forest")
# 

#Balloon_plot figure 4 - 5
primary_df_full[, zscore := b/se]
primary_df_full[, Category := factor(Category, levels = c("Positive control", "Dysbiotic disease", "Fecal Metabolites",  "Microbial Pathway", "Plasma Metabolites",  "Microbe Abundance"  )),]
primary_df_full<- primary_df_full[order(Category),]
primary_df_full[exposure_clean == outcome_clean, c( "b", "se", "pval", "zscore") := NA]
dat_plot <- data.frame(exposure = factor(primary_df_full$outcome_clean, levels = unique(primary_df_full$outcome_clean)), #C,est mélangeant, mias puisque je veux les exposures sur l'ave des X j'inverse
                       outcome = factor(primary_df_full$exposure_clean, levels = unique(primary_df_full$exposure_clean)),
                       pval = as.numeric(primary_df_full$pval),
                       z_score = as.numeric(primary_df_full$zscore),
                       beta_score = as.numeric(primary_df_full$b),
                       Category = primary_df_full$Category,
                       outcome_category = primary_df_full$outcome_category)
                       #Category = factor(primary_df_full$Category, levels = sort(unique(primary_df_full$Category))))

dat_plot <- dat_plot[dat_plot$Category != "Dysbiotic disease",]
setDT(dat_plot)
dat_plot <- dat_plot[order(Category),]
###########
plot_balloon <- function(dat_plot, bonferroni_threshold = 0.05) {
  
  
  dat_plot$pval = as.numeric(dat_plot$pval)
  dat_plot$log10_pval = -log10(dat_plot$pval)
  dat_plot$shape_point = sapply(dat_plot$pval, FUN = function(x) {ifelse( x < bonferroni_threshold, "rond", "Non-significant")})
  
  dat_plot_rond = dat_plot
  dat_plot_rond$beta_score[which(dat_plot_rond$pval >= bonferroni_threshold)] = NA
  dat_plot_rond$log10_pval[which(dat_plot_rond$pval >= bonferroni_threshold)] = NA
  dat_plot_rond$shape_point = sapply(dat_plot_rond$pval, FUN = function(x) {ifelse( x < bonferroni_threshold, "rond", NA)})
  
  
  dat_plot_croix = dat_plot
  dat_plot_croix$beta_score[which(dat_plot_croix$pval < bonferroni_threshold)] = NA
  dat_plot_croix$log10_pval[which(dat_plot_croix$pval < bonferroni_threshold)] = NA
  dat_plot_croix$shape_point = sapply(dat_plot_croix$pval, FUN = function(x) {ifelse( x >= bonferroni_threshold, "Non-significant", NA)})
  
  balloon_plot = ggplot2::ggplot() +
    ggplot2::geom_point(data = dat_plot_croix, ggplot2::aes(x = exposure, y = outcome, shape = factor(shape_point)), size = 2, color = "gray20") +
    ggplot2::geom_point(data = dat_plot_rond, ggplot2::aes(x = exposure, y = outcome, size = log10_pval, color = beta_score)) +
    facet_grid(. ~  outcome_category, scales= "free_x") +
    ggplot2::scale_color_gradient2(name = "Beta",
                                   low = scales::muted("#5884E5"),
                                   mid = "white",
                                   high = scales::muted("#9E131E"),
                                   midpoint = 0
    ) +
    ggplot2::scale_shape_manual(name = "", values = c(4,1)) +
    ggplot2::scale_size(name = expression(-Log[10](P)), range = c(4,7)) +
    # ggplot2::coord_fixed(clip = "off", ratio = 1) +
    ggplot2::guides(size = guide_legend(order = 1),
                    shape = guide_legend(order = 2)) +
    ggplot2::theme(
      panel.grid.major.y = element_line(size = 0.25, colour = "gray60"),
      panel.grid.major.x = element_line(size = 0.25, colour = "gray60"),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(t = 2, r = 0.5, b = 0.5, l = 0.5, "cm"),
      legend.position = "right",
      legend.text = element_text(
        color = "gray20",
        size = 10,
        margin = margin(l = 0.2, r = 0.2)
      ),
      legend.spacing.y = unit(0.1, 'cm'),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size = unit(0.8, "cm"),
      axis.title = element_blank(),
      axis.line = element_line(size = 0.5, colour = "gray20"),
      axis.ticks = element_line(size = 0.5, colour = "gray20"),
      axis.text.y = element_text(
        size = 10,
        colour = "gray20"
      ),
      axis.text.x = element_text(
        angle = 60,
        size = 8,
        hjust = 1,
        face = "plain",
        colour = "gray20"
      ))
  
  print(balloon_plot)
  
}

#############
plot_balloon(dat_plot[dat_plot$Category != "Microbe Abundance",])

ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Figures/2_balloon_plot_nomicrobre.png"),
       width=638/72,height=458/72,units="in",scale=1,
       device = "png")


dat_plot2 <- dat_plot[dat_plot$Category %in% c("Microbe Abundance", "Positive control"), ]
plot_balloon(dat_plot2)
ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Figures/3_balloon_plot_onlymicrobre.png"),
       width=750/72,height=1400/72,units="in",scale=1, dpi=700,
       device = "png")

# dat_plot2 <- dat_plot[dat_plot$Category == "Microbe Abundance", ]
# 
# kk <-unique(dat_plot2$outcome)
# split_kk <- split(kk, ceiling(seq_along(kk)/(length(kk)/2+1)))
# 
# dat_plot2 <- dat_plot[(dat_plot$Category == "Microbe Abundance" & 
#                          dat_plot$outcome %in% split_kk[[1]]) |
#                         (dat_plot$Category ==  "Positive control"),]
# 
# 
# plot_balloon(dat_plot2)
# 
# ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Figures/3b_balloon_plot_onlymicrobre.png"),
#        width=750/72,height=900/72,units="in",scale=1, dpi=700,
#        device = "png")
# 
# dat_plot3 <- dat_plot[(dat_plot$Category == "Microbe Abundance" & 
#                          dat_plot$outcome %in% split_kk[[2]]) |
#                         (dat_plot$Category ==  "Positive control"),]
# 
# 
# plot_balloon(dat_plot3)
# 
# ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Figures/3b_balloon_plot_onlymicrobre.png"),
#        width=750/72,height=900/72,units="in",scale=1, dpi=700,
#        device = "png")

#######Sensitivity forest
#sensitivity
list_sensitivity <- readRDS( "Data/Modified/Sensitivity/list_sensitivity")
veclog <- readRDS( "Data/Modified/Sensitivity/veclog")

# x <- list_sensitivity[[1]]
# x[method %in% c("Weighted mode", "Weighted median"), any(pval<0.05)] & 
#   x[method %in% c("IVW radial", "Weighted median"), any(pval<0.05)]
# vec_which<-sapply(list_sensitivity, function(x) x[method %in% c("Weighted mode", "Weighted median"), pval]
# df_compare_traits_groups


length(list_sensitivity)
sum(veclog)
dt_sen <- rbindlist(list_sensitivity[veclog], fill = TRUE)
dt_sen <- cleanify(dt_sen)

dt_sen <- dt_sen[order(exposure_clean, outcome_clean), ]
max_chr_exp <-dt_sen[, max(nchar(exposure_clean))]
dt_sen[, outcome_clean := as.character(outcome_clean)]
max_chr_out <-dt_sen[, max(nchar(outcome_clean))]

#exposure_really_clean
# dt_sen$exposure_really_clean <- 
#   dt_sen[, paste0(exposure_clean, paste(rep("-", max_chr_exp - nchar(exposure_clean)),collapse = "")), by = seq_len(nrow(dt_sen))]$V1
# 
# 
# dt_sen$outcome_really_clean <- 
#   dt_sen[, paste0(paste(rep("-", max_chr_out - nchar(outcome_clean)),collapse = ""), outcome_clean), by = seq_len(nrow(dt_sen))]$V1
# 

dt_sen[, exposure_really_clean := paste0(exposure_clean, "  -->  ")]
dt_sen[, outcome_really_clean := outcome_clean]

dt_sen[, method := factor(method, levels = unique(method))]
dt_sen[, name := paste0(exposure_really_clean,  outcome_really_clean)]
dt_sen[, unique(nchar(name))]
dt_sen[, Category_other := ifelse(Category == "Microbe Abundance","Microbe Abundance", "Microbiota Associated Metabolites")]
dt_sen[method == "Contamination mixture", se := ((b - lci) + (uci - b))/(2*1.96) ]


my_forest_plot <- function(data) {
forestplot(
  df = data,
  name = name,
  se = se,
  estimate = b,
  pvalue = pval,
  psignif = 0.05,
  xlab = "Effect size (SD or log(OR)) per 1-SD\nincrease in gut microbiota features",
  ci = 0.95,
  colour = method,
  xlim = dt_sen[method != "MR Egger", round(c(min(lci),max(uci)), digits = 1)]
) + theme(legend.position = "right") +
  theme(text = element_text(size = 10)) +
  scale_x_continuous(breaks = c(-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4)) +
  ggforce::facet_col(
    # facets = ~outcome_category,
    facets = ~Category_other,
    scales = "free_y",
    space = "free"
  )
}

my_forest_plot(data= dt_sen)
ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Figures/4_forest_sensitivity.png"),
       width=732/72,height=437/72,units="in",scale=1,
       device = "png")

# Histogram of power figure 6 
primary_df <- primary_df_full
primary_df <- primary_df[exposure != "willer_LDL-cholesterol:ieu-a-300"]
under.8<-primary_df[, mean(power), by = outcome][V1 <0.8, ]$outcome
primary_df[, category := ifelse(outcome %in% under.8, outcome, "other")]

sum(primary_df$power<0.8)/length(primary_df$power)

ggplot(primary_df, aes(x=power)) +
  geom_histogram(binwidth = 0.1) +
  scale_x_continuous(breaks =seq(0.2, 1, by = 0.2)) +
  theme_classic() +
  annotate(geom = "text", x = 0.8, y = 300, 
           label = paste0(round(sum(primary_df$power>0.8)/length(primary_df$power)*100,1), "% \n over 0.8"),
           hjust = 0, vjust = 1, size = 5) +
  geom_vline(xintercept=0.8, linetype="dotted", colour = "red") +
  ggtitle(paste0("Power to find a 0.1 effect for all ", nrow(primary_df)," associations" ))

ggsave("Results/Figures/5_histogram_of_power.png",
       width=460/72,height=326/72, units="in", scale=1, device = "png")

#MVMR forest plot figure 7

MVMR_dat <- fread("/home/gagelo01/workspace/Projects/Dysbiose_project/Data/Modified/Sensitivity/MVMR")
MVMR_dat[, MVMR := MVMR %>% ifelse(. == "no counfounder", "no confounder", .)]
MVMR_dat <- cleanify(MVMR_dat)
MVMR_dat <- MVMR_dat[order(exposure_clean, outcome_clean), ]
max_chr_exp <-MVMR_dat[, max(nchar(exposure_clean))]
MVMR_dat[, outcome_clean := as.character(outcome_clean)]
max_chr_out <-MVMR_dat[, max(nchar(outcome_clean))]
MVMR_dat[,outcome_clean := as.character(outcome_clean)]
MVMR_dat[, exposure_really_clean := paste0(exposure_clean, "  -->  ")]
MVMR_dat[, outcome_really_clean := outcome_clean]

MVMR_dat[, method := factor(MVMR, levels = c("no confounder", "with BMI", "with alcohol intake frequency"))]
MVMR_dat[, name := paste0(exposure_really_clean,  outcome_really_clean)]
MVMR_dat[, Category_other := ifelse(Category == "Microbe Abundance","Microbe Abundance", "Microbiota Associated Metabolites")]

my_forest_plot(MVMR_dat)

ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Figures/MVMR_forest_plot.png"),
       width=730/72,height=430/72,units="in",scale=1,
       device = "png")

fwrite(MVMR_dat, "Data/Modified/MVMR_clean.txt")
# list_MVMR <- list("1a_MVMR_forest_plot" = MVMR_dat[outcome_category == "Disease",], 
#                   "1b_MVMR_forest_plot" = MVMR_dat[outcome_category == "Risk_factor",])
# 
# for(i in 1:length(list_MVMR)) {
# 
# MVMR <- list_MVMR[[i]]
# resultsA <- data.frame(variable = as.character(1:nrow(MVMR)),
#                        estimate = round(MVMR$b, digits =3),
#                        lci =  round(MVMR$lci, digits = 3),
#                        uci =  round(MVMR$uci, digits = 3),
#                        P_value = formatC(MVMR$pval, format = "e", digits = 1),
#                        n = MVMR$nsnp)
# 
# mylabels <- data.frame(heading1 = "",
#                        heading2 = paste0(MVMR$exposure_clean, " with ", MVMR$outcome_clean),
#                        heading3 = MVMR$MVMR,
#                        variable = as.character(1:nrow(MVMR)))
# 
# make_forest_plot(panels = list(resultsA),
#                  col.key = "variable",
#                  panel.headings = "IVW estimates with and without correcting for confounder",
#                  row.labels = mylabels,
#                  exponentiate = TRUE,
#                  pointsize = 2, 
#                  rows = unique(mylabels$heading1),
#                  col.stderr = NULL,
#                  col.lci = "lci",
#                  col.uci = "uci",
#                  col.left         = c("n"),
#                  col.left.heading = c("n SNP"),
#                  col.right = "P_value",
#                  col.right.heading = c("OR (95% CI)", "P-value"),
#                  xlab = "IVW estimates",
#                  colour = "blue",
#                  cicolour = "darkred",
#                  blankrows = c(0,1,0,0),
#                  xlim	= c(0.65, 1.8)
# )
# 
# ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Supplementary_material/", names(list_MVMR)[i],".png"),
#        width=710/72,height=576/72, units="in", scale=1,
#        device = "png")
# 
# }

#10-11 -> For all IVW estimate with p < 0.05, perform sensitivity analysis and present the scatter plot
list_sensitivity <- readRDS( "/home/gagelo01/workspace/Projects/Dysbiose_project/Data/Modified/Sensitivity/list_sensitivity")
list_sensitivity <- list_sensitivity[sapply(list_sensitivity, function(x) sum(x$pval < 0.05)) > 3]
list_sensitivity <- lapply(list_sensitivity, cleanify)
list_sensitivity <- lapply(list_sensitivity, function(x)  x[, c("exposure", "outcome", "exposure_copy", "outcome_copy") := .(exposure_clean, outcome_clean, exposure, outcome)])

harm_all <- fread( "Data/Modified/harm_all.txt")
harm_all <- cleanify(harm_all)
fwrite(harm_all, "Data/Modified/harm_all_clean.txt")
harm_all[, id.exposure := exposure]
harm_all[, id.outcome := outcome]
list_harm <- split(harm_all, harm_all$exposure_outcome)

list_harm<-list_harm[!sapply(list_harm, function(x) { x[1,"SNP"] == 0})]
lapply(list_harm, setDT)
list_harm <- lapply(list_harm, cleanify)
list_harm <- lapply(list_harm, function(x)  x[, c("exposure", "outcome") := .(exposure_clean, outcome_clean)])

for(j in 1:length(list_sensitivity)) {
  vecto <- rep(NA, length(list_harm))
for(i in 1:length(list_harm)) {
  expi <- list_harm[[i]][1, "exposure"] == list_sensitivity[[j]][1,"exposure"]
  outci <- list_harm[[i]][1, "outcome"]  == list_sensitivity[[j]][1,"outcome"]
  vecto[i] <-  expi & outci

}

index <- which(vecto)
li_sen <- list_sensitivity[[j]]

p <- mr_scatter_plot(li_sen, list_harm[[index]])
ggplot_object <- p[[1]] + theme_classic() +
 theme(legend.position="top")
ggplot_object
ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Supplementary_material/Scatter_plot/",
       j, "_", list_sensitivity[[j]][1,"exposure_copy"], "__" ,list_sensitivity[[j]][1,"outcome_copy"],"_", "scatter_plot.png"),
       width=538/72,height=335/72, units="in", scale=1,
       device = "png")

saveRDS(ggplot_object, file = paste0("Results/Supplementary_material/ggplot_object/", j, "_scatter.rdata"))
}


#8-9 -> For all IVW estimate with p < 0.05 and nSNPs > 1, perform  a forest plot that look like Iyas and all.-----

list_sensitivity_forest <- vector(mode = "list", length = length(list_sensitivity))

for(i in 1:length(list_sensitivity)) {
  li_sen <- list_sensitivity[[i]]
  if(li_sen[method == "Robust adjusted profile score (RAPS)", uci - lci] > 8) {
    li_sen[method == "Robust adjusted profile score (RAPS)", uci := b + 4]
    li_sen[method == "Robust adjusted profile score (RAPS)", lci := b - 4]
  }
  
  if(li_sen[1,"outcome"] == "Parental lifespan") {
    exponentiate <- FALSE
    col.right.heading <- c("Beta (95% CI)", "P-value")
    xlab <- paste0("1 year change in ", li_sen[1,]$outcome_clean, " per 1-SD increase of ", li_sen[1,]$exposure_clean)
    nullval <- 0
  } else {
    exponentiate <- TRUE
    col.right.heading <- c("OR (95% CI)", "P-value")
    xlab <- paste0("Odds ratio of ", li_sen[1,]$outcome_clean, " per 1-SD increase of ", li_sen[1,]$exposure_clean)
    nullval <- 1
  }
  
  
  
  resultsA <- data.frame(variable = li_sen$method,
                         estimate = round(li_sen$b, digits =3),
                         lci =  round(li_sen$lci, digits = 3),
                         uci =  round(li_sen$uci, digits = 3),
                         n = li_sen$nsnp,
                         pval = formatC(li_sen$pval, format = "e", digits = 1))
  
  mylabels <- data.frame(heading1 = li_sen$type_of_test,
                        heading2 = li_sen$method,
                         heading3 = as.character(NA),
                         variable = li_sen$method)
  
  list_sensitivity_forest[[i]] <-
    make_forest_plot(panels = list(resultsA),
                     col.key = "variable",
                     row.labels = mylabels,
                     exponentiate = exponentiate,
                     pointsize = 2, 
                     rows = unique(mylabels$heading1),
                     col.stderr = NULL,
                     col.lci = "lci",
                     col.uci = "uci",
                     col.left         = c("n"),
                     col.left.heading = c("n SNP"),
                     col.right.heading = col.right.heading,
                     col.right = "pval",
                     xlab = xlab,
                     nullval = nullval,
                     blankrows = c(0,0,0,0),
                     panel.names = "")
  
  ggplot_object <- list_sensitivity_forest[[i]]
  ggplot_object[[1]] 
  ggsave(paste0("/home/gagelo01/workspace/Projects/Dysbiose_project/Results/Supplementary_material/Forest_plot/",
                i, "_", li_sen[1, "exposure_copy"], "__", li_sen[1, "outcome_copy"], "_forest_plot.png" ),
         width=659/72,height=357/72,units="in",scale=1, device = "png")
  
  saveRDS(ggplot_object, file = paste0("Results/Supplementary_material/ggplot_object/", i, "_forest.rdata"))
  
}

###two panel supplementary figure

for(i in 1:length(list_sensitivity)) {
a  <- readRDS(file = paste0("Results/Supplementary_material/ggplot_object/", i, "_scatter.rdata"))
b <- readRDS(file = paste0("Results/Supplementary_material/ggplot_object/", i, "_forest.rdata"))

ggarrange(a, b[[1]], 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

ggsave(paste0("Results/Supplementary_material/TwoPanel/", i, "_scatter_forest.png"),
       width=490/72,height=475/72,units="in",scale=1,
       device = "png")
}

# ##forest plot between IBD and diseases
# harm_all_clean <- cleanify(harm_all)
# harm_all_clean[exposure_clean == "IBD_IBD", exposure_clean := "Inflammatory bowel disease"]
# harm_all_clean[,exposure := exposure_clean]
# harm_all_clean[,outcome := outcome_clean]
# singlesnp <- mr_singlesnp(dat = harm_all_clean[exposure == "Inflammatory bowel disease"], single_method = "mr_wald_ratio", all_method = "mr_ivw")
# setDT(singlesnp)
# index <- unique(singlesnp[,outcome])
# for(i in 1:length(index)) {
#   s<-mr_forest_plot(singlesnp_results = singlesnp[outcome == index[i]], exponentiate = TRUE) 
#   s[[1]] +
#     theme_classic() +
#     theme(legend.position="none") +
#     theme(axis.text.y = element_text(size=6),
#           axis.title=element_text(size=10))
#   
#   titre<-gsub(" ", "_", index[i])
#   ggsave(paste0("Results/Supplementary_material/Singlesnp/singlesnp_IBD_", titre, ".png"),
#          width=459/72,height=594/72,units="in",scale=1, device = "png")
#   
# }


