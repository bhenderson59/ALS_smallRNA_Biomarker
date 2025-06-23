###Analysis of droplet digital PCR data from Crestwood-HudsonAlpha ALS cohort, NEALS cohort, and Smith Family Clinic Parkinson's disease cohort###

#Requesting interactive session and loading into R
srun -c 10 -p interactive --pty /bin/bash
cd /cluster/home/bhenderson/ALS_Biomarker_Paper/ddPCR_analysis
module load cluster/R/4.1.1
R --vanilla

#Loading required R packages
library(glmnet)
library(broom)
library(dplyr)
library(ggplot2)
library(pROC)

#Creating the R workspace
wrkspc_img <- "ALS_biomarker_ddPCR_anlaysis_03_12_2025.Rvar"
load("ALS_biomarker_ddPCR_anlaysis_03_12_2025.Rvar")

#Loading in sample metadata and ddPCR expression data for each dataset
ALL_metadata <- read.csv("/cluster/home/bhenderson/ALS_biomarkers/miR_206_ddPCR_LASSO/regressions/ALL_SAMPLE_METADATA.csv")
NEALS_miR_206_ddPCR <- read.csv("/cluster/home/bhenderson/ALS_biomarkers/miR_206_ddPCR_LASSO/neals_miR_206_ddPCR.csv")
CHA_ALS_miR_206_ddPCR <- read.csv("/cluster/home/bhenderson/ALS_biomarkers/miR_206_ddPCR_LASSO/cha_als_miR_206_ddPCR.csv")
SFCHAPD_miR_206_ddPCR <- read.csv("/cluster/home/bhenderson/ALS_biomarkers/miR_206_ddPCR_LASSO/sfchapd_miR_206_ddPCR.csv")

#Adjusting ALS status from 0 or 1 to simply No or Yes
NEALS_miR_206_ddPCR$ALS_Status <- ifelse(NEALS_miR_206_ddPCR$ALS_Status == 1, "Yes", "No")
CHA_ALS_miR_206_ddPCR$ALS_Status <- ifelse(CHA_ALS_miR_206_ddPCR$ALS_Status == 1, "Yes", "No")
SFCHAPD_miR_206_ddPCR$ALS_Status <- ifelse(SFCHAPD_miR_206_ddPCR$ALS_Status == 1, "Yes", "No")

#Setting appropriate factor levels for each dataset
NEALS_miR_206_ddPCR$ALS_Status <- factor(NEALS_miR_206_ddPCR$ALS_Status, levels = c("No", "Yes"))
CHA_ALS_miR_206_ddPCR$ALS_Status <- factor(CHA_ALS_miR_206_ddPCR$ALS_Status, levels = c("No", "Yes"))
SFCHAPD_miR_206_ddPCR$ALS_Status <- factor(SFCHAPD_miR_206_ddPCR$ALS_Status, levels = c("No", "Yes"))

#Computing the mean and standard deviation for controls only so that I can then z-score for the CHA-ALS cohort
CHA_ALS_data_mean <- mean(CHA_ALS_miR_206_ddPCR$miR_206_AVG[CHA_ALS_miR_206_ddPCR$DISEASE == "Healthy Control"], na.rm = TRUE)
CHA_ALS_data_sd <- sd(CHA_ALS_miR_206_ddPCR$miR_206_AVG[CHA_ALS_miR_206_ddPCR$DISEASE == "Healthy Control"], na.rm = TRUE)
CHA_ALS_miR_206_ddPCR$miR_206_z <- (CHA_ALS_miR_206_ddPCR$miR_206_AVG - CHA_ALS_data_mean) / CHA_ALS_data_sd

#Computing the mean and standard deviation for controls only so that I can then z-score for the SFCHAPD cohort
SFCHAPD_data_mean <- mean(SFCHAPD_miR_206_ddPCR$miR_206_AVG[SFCHAPD_miR_206_ddPCR$DISEASE == "Healthy Control"], na.rm = TRUE)
SFCHAPD_data_sd <- sd(SFCHAPD_miR_206_ddPCR$miR_206_AVG[SFCHAPD_miR_206_ddPCR$DISEASE == "Healthy Control"], na.rm = TRUE)
SFCHAPD_miR_206_ddPCR$miR_206_z <- (SFCHAPD_miR_206_ddPCR$miR_206_AVG - SFCHAPD_data_mean) / SFCHAPD_data_sd

#Computing the mean and standard deviation for controls only so that I can then z-score for the NEALS cohort
NEALS_206_mean <- mean(NEALS_miR_206_ddPCR$miR_206_AVG[NEALS_miR_206_ddPCR$DISEASE == "Healthy Control"], na.rm = TRUE)
NEALS_206_sd <- sd(NEALS_miR_206_ddPCR$miR_206_AVG[NEALS_miR_206_ddPCR$DISEASE == "Healthy Control"], na.rm = TRUE)
NEALS_miR_206_ddPCR$miR_206_z <- (NEALS_miR_206_ddPCR$miR_206_AVG - NEALS_206_mean) / NEALS_206_sd

save.image(wrkspc_img)

#Merging the z-scored data
all_data_with_z <- bind_rows(NEALS_miR_206_ddPCR, CHA_ALS_miR_206_ddPCR, SFCHAPD_miR_206_ddPCR)

#Merging the z-scored data frame with the metadata data frame
merged_data <- merge(all_data_with_z, ALL_metadata, by = "Sample_ID")

#Splitting the data into controls only. The effect of disease status will preclude the identification of any covariates driving differences in miR-206 expression. To accurately removed these effects, we will look in controls only.
merged_data_only_CTL <- merged_data %>% filter(DISEASE == "Healthy Control")

#Labelling the controls by collection site
merged_data_only_CTL <- merged_data_only_CTL %>%
  mutate(
    SOURCE_ID = case_when(
      STUDY_ID == "CHA_ALS" ~ "Crestwood",
      STUDY_ID == "SFCHAPD" ~ "SFC",
      TRUE ~ "NEALS"  # Default case for all other Study_IDs
    ))

save.image(wrkspc_img)

#Running a simple linear regression to determine if there is an effect by sex on miR-206 expression
CTL_only_sex_lm <- lm(miR_206_z ~ SEX, data = merged_data_only_CTL)
summary(CTL_only_sex_lm)
CTL_only_sex_p_value <- summary(CTL_only_sex_lm)$coefficients[2, 4]

save.image(wrkspc_img)

#Plotting the CTL only z-scored data by sex, with the p-value added
source("/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/analysisScripts/jitterBoxplotXs.R")
source("/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/analysisScripts/changeColorAlpha.R")

benColPal <- c("#543005","#8C510A","#BF812D","#DFC27D","#F6E8C3","#C7EAE5","#80CDC1","#35978F","#01665E","#003C30")

plotVals <- merged_data_only_CTL[,"miR_206_z"]

sexFac <- merged_data_only_CTL$SEX

colKeySex <- rep(changeColorAlpha(benColPal[2],newAlpha=250),
                 length(sexFac))
colKeySex[which(sexFac=="Male")] <- changeColorAlpha(benColPal[8],newAlpha=250)

globalMin <- min(plotVals)
globalMax <- max(plotVals)

globalYlim <- c(globalMin,globalMax)

jitterXs <- jitterBoxplotsXs(plotVals,
                             plotFactor=factor(sexFac,levels = c("Female", "Male")),
                             maxJitter=0.2)

groupMedians <- tapply(plotVals,
                       INDEX=sexFac,
                       FUN=median)

CTL_sex_z <- paste("CTL_sex_miR206_z_boxplot.pdf",sep="_")

pdf(CTL_sex_z,width=1.5,height=2)

par(mai=c(0.4,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=jitterXs,
     y=plotVals,
     pch=1,
     col=colKeySex,
     cex=1,
     xlab="",
     xaxt='n',
     xlim= c(0.7, 2.3),
     ylab="miR-206 z-score",
     ylim=globalYlim)

lines(x=c(0.7,1.3),y=rep(groupMedians[1],2),col="black",lwd=1.5)
lines(x=c(1.7,2.3),y=rep(groupMedians[2],2),col="black",lwd=1.5)

axis(side=1,at=c(1,2),labels = FALSE)

par(xpd=NA)

plotLims <- par('usr')

text(x=c(1,2),y=c(plotLims[3],plotLims[3]),
     labels=c("Female","Male"),
     pos=1,offset=1,
     cex=0.8)

dev.off()

save.image(wrkspc_img)

#Running a simple linear regression to determine if there is an effect by source site

CTL_only_source_lm <- lm(miR_206_z ~ SOURCE_ID, data = merged_data_only_CTL)
summary(CTL_only_source_lm)

save.image(wrkspc_img)
#Plotting the CTL only z-scored data by source 

sourceFac <- merged_data_only_CTL$SOURCE_ID

colKeySource <- rep(changeColorAlpha(benColPal[1],newAlpha=250),
                 length(sourceFac))
colKeySource[which(sourceFac=="NEALS")] <- changeColorAlpha(benColPal[6],newAlpha=250)
colKeySource[which(sourceFac=="SFC")] <- changeColorAlpha(benColPal[4],newAlpha=250)

globalMin <- min(plotVals)
globalMax <- max(plotVals)

globalYlim <- c(globalMin,globalMax)

jitterXs <- jitterBoxplotsXs(plotVals,
                             plotFactor=factor(sourceFac,levels = c("Crestwood", "NEALS", "SFC")),
                             maxJitter=0.2)

groupMedians <- tapply(plotVals,
                       INDEX=sourceFac,
                       FUN=median)

CTL_source_z <- paste("CTL_source_miR206_z_boxplot.pdf",sep="_")

pdf(CTL_source_z,width=1.5,height=2)

par(mai=c(0.4,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=jitterXs,
     y=plotVals,
     pch=1,
     col=colKeySource,
     cex=1,
     xlab="",
     xaxt='n',
     yaxt='n',
     xlim= c(0.7, 3.3),
     ylab="",
     ylim=globalYlim,
     bty='n')

lines(x=c(0.7,1.3),y=rep(groupMedians[1],2),col="black",lwd=1.5)
lines(x=c(1.7,2.3),y=rep(groupMedians[2],2),col="black",lwd=1.5)
lines(x=c(2.7,3.3),y=rep(groupMedians[3],2),col="black",lwd=1.5)

axis(side=1,at=c(1,2,3),labels = FALSE, mgp = c(2,4,0))

par(xpd=NA)

plotLims <- par('usr')

text(x=c(1,2,3),y=c(plotLims[3],plotLims[3],plotLims[3]),
     labels=c("CHA","NEALS","SFC"),
     pos=1,offset=1,
     srt = 45,
     cex=0.8)

dev.off()

save.image(wrkspc_img)

load("/cluster/home/bhenderson/ALS_Biomarker_Paper/ddPCR_analysis/ALS_biomarker_ddPCR_anlaysis_03_12_2025.Rvar")

#Running a simple linear regression to see the effects of miR-206z by age

CTL_only_AGE_model <- lm(miR_206_z ~ AGE_COLLECT, data = merged_data_only_CTL)
summary(CTL_only_AGE_model)

#Plotting the CTL only z-scored data by age

source_colors <- c("Crestwood" = benColPal[1], 
                   "SFC" = benColPal[4],
                   "NEALS" = benColPal[6])

miR_206z_CTL_age <- ggplot(merged_data_only_CTL, aes(x = AGE_COLLECT, y = miR_206_z, color = SOURCE_ID)) +
  geom_point(size = 3, alpha = 0.8, shape = 16)+
  scale_color_manual(values = source_colors)+
  geom_smooth(method = "lm", color = "#e31a1c", linewidth = 1.5, se = FALSE) + 
  annotate("text", x = max(merged_data_only_CTL$AGE_COLLECT) - 5, y = max(merged_data_only_CTL$miR_206_z) - 1, label = "")+
  labs(title = "z-scored miR-206 by age at collection",
       x = "Age (years)",
       y = "miR-206 z-score" )+
  theme_minimal(base_size = 14, base_family = "sans") +  # Change font family here
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(color = "black"),
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )

pdf("miR206z_CTL_by_Age.pdf", width = 8, height = 6)
print(miR_206z_CTL_age)
dev.off()

save.image(wrkspc_img)

###At this point there are no significant effects in the controls based on sex, source, or age. Now we will plot the miR-206 z-scored values for all samples, and then create a linear model that wil predict ALS status.

#Labelling the all samples by collection site
merged_data <- merged_data %>%
  mutate(
    SOURCE_ID = case_when(
      STUDY_ID == "CHA_ALS" ~ "Crestwood",
      STUDY_ID == "SFCHAPD" ~ "SFC",
      TRUE ~ "NEALS"  # Default case for all other Study_IDs
    ))

#Running a simple linear regression to determine if there is an effect by sex on miR-206 expression
all_sample_sex_lm <- lm(miR_206_z ~ SEX, data = merged_data)
summary(all_sample_sex_lm)
summary_all_sample_sex_lm <- summary(all_sample_sex_lm)
summary_all_sample_sex_lm_results_table <- summary_all_sample_sex_lm$coefficients
results_summary_all_sample_sex_lm_tidy <- tidy(all_sample_sex_lm)
write.csv(summary_all_sample_sex_lm_results_table, "all_sample_sex_lm_results.csv", row.names = FALSE)
write.csv(results_summary_all_sample_sex_lm_tidy, "all_sample_sex_lm_results_tidy.csv", row.names = FALSE)
all_sample_sex_p_value <- summary(all_sample_sex_lm)$coefficients[2, 4]

#Plotting out all sample miR-206 z-score
all_sample_plotVals <- merged_data[,"miR_206_z"]

all_sample_sexFac <- merged_data$SEX

all_sample_colKeySex <- rep(changeColorAlpha(benColPal[2],newAlpha=250),
                 length(all_sample_sexFac))
all_sample_colKeySex[which(all_sample_sexFac=="Male")] <- changeColorAlpha(benColPal[8],newAlpha=250)

all_sample_globalMin <- min(all_sample_plotVals)
all_sample_globalMax <- max(all_sample_plotVals)

all_sample_globalYlim <- c(all_sample_globalMin,all_sample_globalMax)

all_sample_jitterXs <- jitterBoxplotsXs(all_sample_plotVals,
                             plotFactor=factor(all_sample_sexFac,levels = c("Female", "Male")),
                             maxJitter=0.2)

all_sample_groupMedians <- tapply(all_sample_plotVals,
                       INDEX=all_sample_sexFac,
                       FUN=median)

all_sample_sex_z <- paste("all_sample_sex_miR206_z_boxplot.pdf",sep="_")

pdf(all_sample_sex_z,width=1.5,height=2)

par(mai=c(0.4,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=all_sample_jitterXs,
     y=all_sample_plotVals,
     pch=1,
     col=all_sample_colKeySex,
     cex=1,
     xlab="",
     xaxt='n',
     xlim= c(0.7, 2.3),
     ylab="miR-206 z-score",
     ylim=all_sample_globalYlim)

lines(x=c(0.7,1.3),y=rep(all_sample_groupMedians[1],2),col="black",lwd=1.5)
lines(x=c(1.7,2.3),y=rep(all_sample_groupMedians[2],2),col="black",lwd=1.5)

axis(side=1,at=c(1,2),labels = FALSE)

par(xpd=NA)

plotLims <- par('usr')

text(x=c(1,2),y=c(plotLims[3],plotLims[3]),
     labels=c("Female","Male"),
     pos=1,offset=1,
     cex=0.8)

dev.off()

save.image(wrkspc_img)

#Running a simple linear regression to determine if there is an effect by source site on all samples

all_sample_source_lm <- lm(miR_206_z ~ SOURCE_ID, data = merged_data)
summary(all_sample_source_lm)

save.image(wrkspc_img)
#Plotting the CTL only z-scored data by source 

all_sample_sourceFac <- merged_data$SOURCE_ID

all_sample_colKeySource <- rep(changeColorAlpha(benColPal[1],newAlpha=250),
                    length(all_sample_sourceFac))
all_sample_colKeySource[which(sourceFac=="NEALS")] <- changeColorAlpha(benColPal[6],newAlpha=250)
all_sample_colKeySource[which(sourceFac=="SFC")] <- changeColorAlpha(benColPal[4],newAlpha=250)

all_sample_source_jitterXs <- jitterBoxplotsXs(all_sample_plotVals,
                             plotFactor=factor(all_sample_sourceFac,levels = c("Crestwood", "NEALS", "SFC")),
                             maxJitter=0.2)

all_sample_source_groupMedians <- tapply(all_sample_plotVals,
                       INDEX=all_sample_sourceFac,
                       FUN=median)

all_sample_source_z <- paste("all_sample_source_miR206_z_boxplot.pdf",sep="_")

pdf(all_sample_source_z,width=1.5,height=2)

par(mai=c(0.4,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=all_sample_source_jitterXs,
     y=all_sample_plotVals,
     pch=1,
     col=all_sample_colKeySource,
     cex=1,
     xlab="",
     xaxt='n',
     xlim= c(0.7, 3.3),
     ylab="miR-206 z-score",
     ylim=all_sample_globalYlim)

lines(x=c(0.7,1.3),y=rep(all_sample_source_groupMedians[1],2),col="black",lwd=1.5)
lines(x=c(1.7,2.3),y=rep(all_sample_source_groupMedians[2],2),col="black",lwd=1.5)
lines(x=c(2.7,3.3),y=rep(all_sample_source_groupMedians[3],2),col="black",lwd=1.5)

axis(side=1,at=c(1,2,3),labels = FALSE, mgp = c(2,4,0))

par(xpd=NA)

plotLims <- par('usr')

text(x=c(1,2,3),y=c(plotLims[3],plotLims[3],plotLims[3]),
     labels=c("CHA","NEALS","SFC"),
     pos=1,offset=1,
     srt = 45,
     cex=0.8)

dev.off()

save.image(wrkspc_img)

#Running a simple linear regression to see the effects of miR-206z by age on all samples

all_sample_AGE_model <- lm(miR_206_z ~ AGE_COLLECT, data = merged_data)
summary(all_sample_AGE_model)

#Plotting the CTL only z-scored data by age

source_colors <- c("Crestwood" = benColPal[1], 
                   "SFC" = benColPal[4],
                   "NEALS" = benColPal[6])

miR_206z_all_age <- ggplot(merged_data, aes(x = AGE_COLLECT, y = miR_206_z, color = SOURCE_ID)) +
  geom_point(size = 3, alpha = 0.8, shape = 16)+
  scale_color_manual(values = source_colors)+
  geom_smooth(method = "lm", color = "#e31a1c", linewidth = 1.5, se = FALSE) + 
  annotate("text", x = max(merged_data$AGE_COLLECT) - 5, y = max(merged_data$miR_206_z) - 1, label = "")+
  labs(title = "z-scored miR-206 by age at collection",
       x = "Age (years)",
       y = "miR-206 z-score" )+
  theme_minimal(base_size = 14, base_family = "sans") +  # Change font family here
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(color = "black"),
    plot.margin = margin(10, 10, 10, 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  )

pdf("miR206z_all_by_Age.pdf", width = 8, height = 6)
print(miR_206z_all_age)
dev.off()

save.image(wrkspc_img)

#Will run analysis first without any correction for age

#Need to split out the Crestwood data to train on
merged_data_only_CHA <- merged_data %>% filter(SOURCE_ID == "Crestwood")

#Fit logistic regression model using z-scored predictor
logistic_model <- glm(ALS_Status ~ miR_206_z, data = merged_data_only_CHA, family = binomial)

#Display model summary
summary(logistic_model)

# Ensure ALS_Status in test data is numeric for pROC
merged_data_only_CHA$ALS_Status <- as.numeric(merged_data_only_CHA$ALS_Status)

#Calculate ROC for Crestwood training data
CHA_only_roc_obj <- roc(merged_data_only_CHA$ALS_Status, merged_data_only_CHA$miR_206_z)

#Print AUC value
CHA_only_auc_value <- auc(CHA_only_roc_obj)
print(paste("AUC:", CHA_only_auc_value))

save.image(wrkspc_img)
#Determining the optimal threshold for specificity and sensitivity
mask_optimal_CHA_only<- min(which(CHA_only_roc_obj$specificities == max(CHA_only_roc_obj$specificities)))
optimal_thresh <- CHA_only_roc_obj$thresholds[mask_optimal_CHA_only]

#Getting the specificity and sensitivity for the optimal threshold in the Crestwood cohort
CHA_sens_spec <- coords(CHA_only_roc_obj, x = optimal_thresh, input = "threshold", ret = c("sensitivity", "specificity"))
opt_spec <- as.numeric(opt_thresh["specificity"])
opt_sens <- as.numeric(opt_thresh["sensitivity"])
opt_thr <- as.numeric(opt_thresh["threshold"])

# Define diagonal line limits based on actual ROC data range
min_spec <- min(1 - roc_df$specificity)
max_spec <- max(1 - roc_df$specificity)
min_sens <- min(roc_df$sensitivity)
max_sens <- max(roc_df$sensitivity)

#Plotting the ROC-AUC for the Crestwood only cohort
CHA_roc_plot <- ggplot() +
  geom_line(aes(x = 1 - rev(CHA_only_roc_obj$specificities), y = rev(CHA_only_roc_obj$sensitivities)),color = "#1f78b4", linewidth = 1.5) +                # ROC curve
  geom_segment(data = NULL, aes(x = min(1 - CHA_only_roc_obj$specificities), 
                                y = min(CHA_only_roc_obj$sensitivities), 
                                xend = max(1 - CHA_only_roc_obj$specificities), 
                                yend = max(CHA_only_roc_obj$sensitivities)), 
               linetype = "dashed", color = "gray")  +
  geom_point(data = data.frame(x = 1 - CHA_sens_spec$specificity, y = CHA_sens_spec$sensitivity), aes(x = x, y = y), 
             color = "red", size = 4) +
  labs(title = "", 
       x = "1 - Specificity", 
       y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  ) +
  annotate("text", x = 0.6, y = 0.1, 
           label = paste("AUC =", round(CHA_only_auc_value, 2)), 
           size = 6, hjust = 0, fontface = "bold")

pdf("CHA_miR_206_roc.pdf", width = 8, height = 6)
print(CHA_roc_plot)
dev.off()

#Plotting the z-scored miR-206 data for the Crestwood cohort
CHA_plotVals <- merged_data_only_CHA[,"miR_206_z"]

CHA_disFac <- merged_data_only_CHA$DISEASE

colKeyDisease <- rep(changeColorAlpha(benColPal[5],newAlpha=250),
                 length(CHA_disFac))
colKeyDisease[which(CHA_disFac=="ALS")] <- changeColorAlpha(benColPal[9],newAlpha=250)

CHA_globalMin <- min(CHA_plotVals)
CHA_globalMax <- max(CHA_plotVals)

CHA_globalYlim <- c(CHA_globalMin,CHA_globalMax)

CHA_jitterXs <- jitterBoxplotsXs(CHA_plotVals,
                             plotFactor=factor(CHA_disFac,levels = c("Healthy Control", "ALS")),
                             maxJitter=0.2)

CHA_groupMedians <- tapply(CHA_plotVals,
                       INDEX=CHA_disFac,
                       FUN=median)

CHA_dis_z <- paste("CHA_miR206_z_boxplot.pdf",sep="_")

pdf(CHA_dis_z,width=1.5,height=2)

par(mai=c(0.4,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=CHA_jitterXs,
     y=CHA_plotVals,
     pch=1,
     col=colKeyDisease,
     cex=1,
     xlab="",
     xaxt='n',
     xlim= c(0.7, 2.3),
     ylab="miR-206 z-score",
     ylim=NEALS_globalYlim)

lines(x=c(0.7,1.3),y=rep(CHA_groupMedians[2],2),col="black",lwd=1.5)
lines(x=c(1.7,2.3),y=rep(CHA_groupMedians[1],2),col="black",lwd=1.5)

x_limits <- par("usr")[1:2] 
x_start <- x_limits[1] + 0.1
x_end <- x_limits[2] - 0.1 
segments(x0 = x_start, x1 = x_end, 
         y0 = optimal_thresh, y1 = optimal_thresh, 
         col = "black", lwd = 1.5, lty = 2)

axis(side=1,at=c(1,2),labels = FALSE)

par(xpd=NA)

plotLims <- par('usr')

text(x=c(1,2),y=c(plotLims[3],plotLims[3]),
     labels=c("CTL","ALS"),
     pos=1,offset=1,
     cex=0.8)

dev.off()

save.image(wrkspc_img)

#Separate out the NEALS data to predict on at first
merged_data_only_NEALS <- merged_data %>% filter(SOURCE_ID == "NEALS")

#Generate predicted probabilities for ALS in the NEALS cohort
merged_data_only_NEALS$Predicted_Prob <- predict(logistic_model, newdata = merged_data_only_NEALS, type = "response")

#Convert probabilities to binary classification (threshold = 0.5)
merged_data_only_NEALS$Predicted_Class <- ifelse(merged_data_only_NEALS$Predicted_Prob > 0.5, 1, 0)

#Display first few predictions
head(merged_data_only_NEALS)

save.image(wrkspc_img)

#Compute and plot the ROC curve for the NEALS only cohort
# Ensure ALS_Status in test data is numeric for pROC
merged_data_only_NEALS$ALS_Status <- as.numeric(merged_data_only_NEALS$ALS_Status)

# Compute ROC curve
NEALS_only_roc_obj <- roc(merged_data_only_NEALS$ALS_Status, merged_data_only_NEALS$miR_206_z)

# Print AUC value
NEALS_only_auc_value <- auc(NEALS_only_roc_obj)
print(paste("AUC:", NEALS_only_auc_value))

# Save results to CSV
write.csv(merged_data_only_NEALS, "NEALS_only_predicted_results_zscored.csv", row.names = FALSE)

#Getting the specificity and sensitivity for the optimal threshold in the NEALS cohort
NEALS_sens_spec <- coords(NEALS_only_roc_obj, x = optimal_thresh, input = "threshold", ret = c("sensitivity", "specificity"))

#Plotting the ROC-AUC for the NEALS only cohort
NEALS_roc_plot <- ggplot() +
  geom_line(aes(x = 1 - rev(NEALS_only_roc_obj$specificities), y = rev(NEALS_only_roc_obj$sensitivities)),color = "#1f78b4", linewidth = 1.5) +                # ROC curve
  geom_segment(data = NULL, aes(x = min(1 - NEALS_only_roc_obj$specificities), 
                                y = min(NEALS_only_roc_obj$sensitivities), 
                                xend = max(1 - NEALS_only_roc_obj$specificities), 
                                yend = max(NEALS_only_roc_obj$sensitivities)), 
               linetype = "dashed", color = "gray")  +
  labs(title = "", 
       x = "1 - Specificity", 
       y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  ) +
  annotate("text", x = 0.6, y = 0.1, 
           label = paste("AUC =", round(NEALS_only_auc_value, 2)), 
           size = 6, hjust = 0, fontface = "bold")

pdf("NEALS_miR_206_roc.pdf", width = 8, height = 6)
print(NEALS_roc_plot)
dev.off()

save.image(wrkspc_img)

#Getting the specificity and sensitivity for the optimal threshold in the NEALS cohort
NEALS_sens_spec <- coords(NEALS_only_roc_obj, x = optimal_thresh, input = "threshold", ret = c("sensitivity", "specificity"))

#Plotting the z-scored miR-206 data for the NEALS cohort
NEALS_plotVals <- merged_data_only_NEALS[,"miR_206_z"]

NEALS_disFac <- merged_data_only_NEALS$DISEASE

colKeyDisease <- rep(changeColorAlpha(benColPal[5],newAlpha=250),
                     length(NEALS_disFac))
colKeyDisease[which(NEALS_disFac=="ALS")] <- changeColorAlpha(benColPal[9],newAlpha=250)

NEALS_globalMin <- min(NEALS_plotVals)
NEALS_globalMax <- max(NEALS_plotVals)

NEALS_globalYlim <- c(NEALS_globalMin,NEALS_globalMax)

NEALS_jitterXs <- jitterBoxplotsXs(NEALS_plotVals,
                                 plotFactor=factor(NEALS_disFac,levels = c("Healthy Control", "ALS")),
                                 maxJitter=0.2)

NEALS_groupMedians <- tapply(NEALS_plotVals,
                           INDEX=NEALS_disFac,
                           FUN=median)

NEALS_dis_z <- paste("NEALS_miR206_z_boxplot.pdf",sep="_")

pdf(NEALS_dis_z,width=1.5,height=2)

par(mai=c(0.4,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=NEALS_jitterXs,
     y=NEALS_plotVals,
     pch=1,
     col=colKeyDisease,
     cex=1,
     xlab="",
     xaxt='n',
     xlim= c(0.7, 2.3),
     ylab="miR-206 z-score",
     ylim=NEALS_globalYlim)

lines(x=c(0.7,1.3),y=rep(NEALS_groupMedians[2],2),col="black",lwd=1.5)
lines(x=c(1.7,2.3),y=rep(NEALS_groupMedians[1],2),col="black",lwd=1.5)

x_limits <- par("usr")[1:2] 
x_start <- x_limits[1] + 0.1
x_end <- x_limits[2] - 0.1 
segments(x0 = x_start, x1 = x_end, 
         y0 = optimal_thresh, y1 = optimal_thresh, 
         col = "black", lwd = 1.5, lty = 2)

axis(side=1,at=c(1,2),labels = FALSE)

par(xpd=NA)

plotLims <- par('usr')

text(x=c(1,2),y=c(plotLims[3],plotLims[3]),
     labels=c("CTL","ALS"),
     pos=1,offset=1,
     cex=0.8)

dev.off()


#Separate out the NEALS data and SFCHA data to predict on
merged_data_NEALS_SFCHA <- merged_data %>% filter(SOURCE_ID %in% c("NEALS", "SFC"))

#Generate predicted probabilities for ALS in the NEALS cohort
merged_data_NEALS_SFCHA$Predicted_Prob <- predict(logistic_model, newdata = merged_data_NEALS_SFCHA, type = "response")

#Convert probabilities to binary classification (threshold = 0.5)
merged_data_NEALS_SFCHA$Predicted_Class <- ifelse(merged_data_NEALS_SFCHA$Predicted_Prob > 0.5, 1, 0)

#Display first few predictions
head(merged_data_NEALS_SFCHA)

save.image(wrkspc_img)

#Compute and plot the ROC curve for the NEALS and SFCHA cohort
# Ensure ALS_Status in test data is numeric for pROC
merged_data_NEALS_SFCHA$ALS_Status <- as.numeric(merged_data_NEALS_SFCHA$ALS_Status)

# Compute ROC curve
NEALS_SFCHA_roc_obj <- roc(merged_data_NEALS_SFCHA$ALS_Status, merged_data_NEALS_SFCHA$miR_206_z)

# Print AUC value
NEALS_SFCHA_auc_value <- auc(NEALS_SFCHA_roc_obj)
print(paste("AUC:", NEALS_SFCHA_auc_value))

# Save results to CSV
write.csv(merged_data_NEALS_SFCHA, "NEALS_SFCHA_predicted_results_zscored.csv", row.names = FALSE)

save.image(wrkspc_img)

#Getting the specificity and sensitivity for the optimal threshold in the NEALS cohort
NEALS_SFCHA_sens_spec <- coords(NEALS_SFCHA_roc_obj, x = optimal_thresh, input = "threshold", ret = c("sensitivity", "specificity"))

#Plotting the ROC-AUC for the NEALS and SFCHA cohort
NEALS_SFCHA_roc_plot <- ggplot() +
  geom_line(aes(x = 1 - rev(NEALS_SFCHA_roc_obj$specificities), y = rev(NEALS_SFCHA_roc_obj$sensitivities)),color = "#1f78b4", linewidth = 1.5) +                # ROC curve
  geom_segment(data = NULL, aes(x = min(1 - NEALS_SFCHA_roc_obj$specificities), 
                                y = min(NEALS_SFCHA_roc_obj$sensitivities), 
                                xend = max(1 - NEALS_SFCHA_roc_obj$specificities), 
                                yend = max(NEALS_SFCHA_roc_obj$sensitivities)), 
               linetype = "dashed", color = "gray")  +
  labs(title = "", 
       x = "1 - Specificity", 
       y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  ) +
  annotate("text", x = 0.6, y = 0.1, 
           label = paste("AUC =", round(NEALS_SFCHA_auc_value, 2)), 
           size = 6, hjust = 0, fontface = "bold")

pdf("NEALS_SFCHA_miR_206_roc.pdf", width = 8, height = 6)
print(NEALS_SFCHA_roc_plot)
dev.off()

save.image(wrkspc_img)

#Plotting the z-scored miR-206 data for the NEALS and SFCHA cohorts

NEALS_SFCHA_plotVals <- merged_data_NEALS_SFCHA[,"miR_206_z"]

NEALS_disFac <- merged_data_NEALS_SFCHA$DISEASE

colKeyDisease <- rep(changeColorAlpha(benColPal[5],newAlpha=250),
                     length(NEALS_disFac))
colKeyDisease[which(NEALS_disFac=="ALS")] <- changeColorAlpha(benColPal[9],newAlpha=250)

NEALS_globalMin <- min(NEALS_plotVals)
NEALS_globalMax <- max(NEALS_plotVals)

NEALS_globalYlim <- c(NEALS_globalMin,NEALS_globalMax)

NEALS_jitterXs <- jitterBoxplotsXs(NEALS_plotVals,
                                 plotFactor=factor(NEALS_disFac,levels = c("Healthy Control", "ALS")),
                                 maxJitter=0.2)

NEALS_groupMedians <- tapply(NEALS_plotVals,
                           INDEX=NEALS_disFac,
                           FUN=median)

NEALS_dis_z <- paste("NEALS_miR206_z_boxplot.pdf",sep="_")

pdf(NEALS_dis_z,width=1.5,height=2)

par(mai=c(0.4,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=NEALS_jitterXs,
     y=NEALS_plotVals,
     pch=1,
     col=colKeyDisease,
     cex=1,
     xlab="",
     xaxt='n',
     xlim= c(0.7, 2.3),
     ylab="miR-206 z-score",
     ylim=NEALS_globalYlim)

lines(x=c(0.7,1.3),y=rep(NEALS_groupMedians[2],2),col="black",lwd=1.5)
lines(x=c(1.7,2.3),y=rep(NEALS_groupMedians[1],2),col="black",lwd=1.5)

x_limits <- par("usr")[1:2] 
x_start <- x_limits[1] + 0.1
x_end <- x_limits[2] - 0.1 
segments(x0 = x_start, x1 = x_end, 
         y0 = optimal_thresh, y1 = optimal_thresh, 
         col = "black", lwd = 1.5, lty = 2)

axis(side=1,at=c(1,2),labels = FALSE)

par(xpd=NA)

plotLims <- par('usr')

text(x=c(1,2),y=c(plotLims[3],plotLims[3]),
     labels=c("CTL","ALS"),
     pos=1,offset=1,
     cex=0.8)

dev.off()


#Separate out the SFC data to predict on
merged_data_only_SFC <- merged_data %>% filter(SOURCE_ID == "SFC")

#Generate predicted probabilities for ALS in the NEALS cohort
merged_data_only_SFC$Predicted_Prob <- predict(logistic_model, newdata = merged_data_only_SFC, type = "response")

#Convert probabilities to binary classification (threshold = 0.5)
merged_data_only_SFC$Predicted_Class <- ifelse(merged_data_only_SFC$Predicted_Prob > 0.5, 1, 0)

#Display first few predictions
head(merged_data_only_SFC)

save.image(wrkspc_img)

#Cannot plot the ROC curve for the SFC only cohort because there are no ALS cases. What I can do is plot the predicted probabilites in the cohort
SFC_only_predicted_prob <- ggplot(merged_data_only_SFC, aes(x = Predicted_Prob, fill = DISEASE)) +
  geom_density(alpha = 0.5) +  
  labs(title = "Predicted Probabilities for CTL vs. PD",
       x = "Predicted Probability",
       y = "Density") +
  scale_fill_manual(values = c("blue", "red")) +  # Adjust colors if needed
  theme_minimal()


pdf("SFC_only_predicted_prob.pdf", width = 8, height = 6)
print(SFC_only_predicted_prob)
dev.off()

save.image(wrkspc_img)


#Plotting the z-scored miR-206 data for the SFC cohort
SFC_plotVals <- merged_data_only_SFC[,"miR_206_z"]

SFC_disFac <- merged_data_only_SFC$DISEASE

colKeyDisease <- rep(changeColorAlpha(benColPal[5],newAlpha=250),
                     length(SFC_disFac))
colKeyDisease[which(SFC_disFac=="PD")] <- changeColorAlpha(benColPal[7],newAlpha=250)

SFC_globalMin <- min(SFC_plotVals)
SFC_globalMax <- max(SFC_plotVals)

SFC_globalYlim <- c(SFC_globalMin,SFC_globalMax)

SFC_jitterXs <- jitterBoxplotsXs(SFC_plotVals,
                                   plotFactor=factor(SFC_disFac,levels = c("Healthy Control", "PD")),
                                   maxJitter=0.2)

SFC_groupMedians <- tapply(SFC_plotVals,
                             INDEX=SFC_disFac,
                             FUN=median)

SFC_dis_z <- paste("SFC_miR206_z_boxplot.pdf",sep="_")

pdf(SFC_dis_z,width=1.5,height=2)

par(mai=c(0.4,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=SFC_jitterXs,
     y=SFC_plotVals,
     pch=1,
     col=colKeyDisease,
     cex=1,
     xlab="",
     xaxt='n',
     xlim= c(0.7, 2.3),
     ylab="miR-206 z-score",
     ylim=NEALS_globalYlim)

lines(x=c(0.7,1.3),y=rep(SFC_groupMedians[2],2),col="black",lwd=1.5)
lines(x=c(1.7,2.3),y=rep(SFC_groupMedians[1],2),col="black",lwd=1.5)

x_limits <- par("usr")[1:2] 
x_start <- x_limits[1] + 0.1
x_end <- x_limits[2] - 0.1 
segments(x0 = x_start, x1 = x_end, 
         y0 = optimal_thresh, y1 = optimal_thresh, 
         col = "black", lwd = 1.5, lty = 2)

axis(side=1,at=c(1,2),labels = FALSE)

par(xpd=NA)

plotLims <- par('usr')

text(x=c(1,2),y=c(plotLims[3],plotLims[3]),
     labels=c("CTL","PD"),
     pos=1,offset=1,
     cex=0.8)

dev.off()

save.image(wrkspc_img)

#Including sex in the logistic regression model to see if it has any affect on the results

# Ensure ALS_Status in test data is numeric for pROC
merged_data_only_CHA$ALS_Status <- as.(merged_data_only_CHA$ALS_Status)

#Fit logistic regression model using z-scored predictor
logistic_model_with_sex <- glm(ALS_Status ~ miR_206_z + SEX, data = merged_data_only_CHA, family = binomial)

#Display model summary
summary(logistic_model_with_sex)

#Predicted probabilities with sex in the model
with_sex_CHA_only_predicted_probabilities <- predict(logistic_model_with_sex, type = "response")

#Calculate ROC for Crestwood training data with sex
with_sex_CHA_only_roc_obj <- roc(merged_data_only_CHA$ALS_Status, with_sex_CHA_only_predicted_probabilities)

#Print AUC value
with_sex_CHA_only_auc_value <- auc(with_sex_CHA_only_roc_obj)
print(paste("AUC:", with_sex_CHA_only_auc_value))


#Plotting the ROC-AUC for the NEALS only cohort
with_sex_CHA_roc_plot <- ggplot() +
  geom_line(aes(x = 1 - rev(with_sex_CHA_only_roc_obj$specificities), y = rev(with_sex_CHA_only_roc_obj$sensitivities)),color = "#1f78b4", linewidth = 1.5) +                # ROC curve
  geom_segment(data = NULL, aes(x = min(1 - with_sex_CHA_only_roc_obj$specificities), 
                                y = min(with_sex_CHA_only_roc_obj$sensitivities), 
                                xend = max(1 - with_sex_CHA_only_roc_obj$specificities), 
                                yend = max(with_sex_CHA_only_roc_obj$sensitivities)), 
               linetype = "dashed", color = "gray")  +
  geom_point(data = data.frame(x = 1 - with_sex_CHA_sens_spec$specificity, y = with_sex_CHA_sens_spec$sensitivity), aes(x = x, y = y), 
             color = "red", size = 4) +
  labs(title = "", 
       x = "1 - Specificity", 
       y = "Sensitivity") +
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    text = element_text(color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14)
  ) +
  annotate("text", x = 0.6, y = 0.1, 
           label = paste("AUC =", round(with_sex_CHA_only_auc_value, 2)), 
           size = 6, hjust = 0, fontface = "bold")

pdf("with_sex_CHA_miR_206_roc.pdf", width = 8, height = 6)
print(with_sex_CHA_roc_plot)
dev.off()

save.image(wrkspc_img)
