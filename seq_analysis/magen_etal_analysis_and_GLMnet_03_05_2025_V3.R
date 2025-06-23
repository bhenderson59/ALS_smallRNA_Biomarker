#srun -c 40 -p interactive --pty /bin/bash
#cd /cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/Magen_etal_ALS_small_RNA
#module load cluster/R/4.4.1
#R --vanilla

source("/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/analysisScripts/explore_PCs.R")
source("/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/analysisScripts/jitterBoxplotXs.R")
source("/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/analysisScripts/changeColorAlpha.R")

require(glmnet)
require(matrixStats)
require(doMC)
library(BiocManager)
library(GenomicRanges)
library(sva)
#BiocManager::install("sva")

#set how many cores here. Make sure it's number of cores requested to slurm -2

nCores <- 38
registerDoMC(cores = nCores)

wrkspc_img <- "magenEtAl_Analysis_03_05_2025.rd"
#load(wrkspc_img)

#bring in raw cnts mat
raw_cnts_mat_File <- "/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/Magen_etal_ALS_small_RNA/magen_etal_combined_counts_mat.tsv"
raw_cnts_Df <- read.table(raw_cnts_mat_File,sep="\t",header=TRUE)

raw_cnts_mat <- as.matrix(raw_cnts_Df[,2:ncol(raw_cnts_Df)])
rownames(raw_cnts_mat) <- as.character(raw_cnts_Df[,1])

head(raw_cnts_mat[,1:10])
tail(raw_cnts_mat[,1:10])


#bring in sample key metadata
sample_key_file <- "/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/Magen_etal_ALS_small_RNA/magen_etal_combined_sample_key.tsv"
sample_key_df <- read.table(sample_key_file,sep="\t",header=TRUE)

head(sample_key_df)

sample_key_df$Disease <- factor(sample_key_df$status,levels=c("cntrl","ALS"))

#make batch a factor
batchFac <- paste("batch",sample_key_df$Batch,sep="_")
batchFac <- as.factor(batchFac)

sample_key_df$batchFac <- batchFac

#plot samples in order of total counts to look for min cnt cutoff
col_key_ALS_status <- rep("gray",nrow(sample_key_df))
col_key_ALS_status[which(sample_key_df$Disease=="ALS")] <- "red"

cnts_colSums <- colSums(raw_cnts_mat)

orderByCnts <- order(cnts_colSums)

min(cnts_colSums)

#png(file="Magen_etal_totMiRNA_cnts.png",width=10,height=2.5,units="in",res=1000,type="cairo")
pdf(file="Magen_etal_totMiRNA_cnts.pdf",width=10,height=2.5)
par(bty="l")
par(mai=c(0.8,0.7,0.1,0.1))
par(mgp=c(2.6,0.6,0))
barplot(cnts_colSums[orderByCnts],
        cex.axis=0.6,
        xlab="",
        ylab="aligned miRNA reads",
        las=2,
        cex.names=0.3,
        cex.lab=0.8,
        border=NA,
        col=col_key_ALS_status[orderByCnts])
dev.off()

#make a plot of fraction ALS vs tot counts threshold
sortTotCnts <- cnts_colSums[orderByCnts]
ALS_N <- rep(0,nrow(sample_key_df))
ALS_N[which(sample_key_df$Disease=="ALS")] <- 1

sortALS_N <- ALS_N[orderByCnts]

fracALS <- (sum(sortALS_N)-cumsum(sortALS_N))/c(nrow(sample_key_df):1)

pdf(file="Magen_etal_miRNACntsVFracALS.pdf",width=4,height=2.5)
par(bty="l")
par(mai=c(0.5,0.5,0.03,0.1))
par(mgp=c(2.6,0.6,0))
plot(x=sortTotCnts,
     y=fracALS,
     xlab="aligned miRNA reads",
     ylab="Fraction ALS samples",
     type="l",
     col="red",
     lwd=1.5,
     cex.axis=0.7,
     cex.lab=0.8,
     log='x')

dev.off()



#look at mir_486-5p as a blood cell contaiminant
mask_486_3p <- grep("hsa-miR-486-3p",as.character(rownames(raw_cnts_mat)))

#get cpm
cpm_mat <- t(apply(raw_cnts_mat,1,function(row_x){
	cpm_row_x <- 1e6*row_x/cnts_colSums
}))

dim(cpm_mat)

tot_486_3p_cpm <- cpm_mat[mask_486_3p,]

orderBy486Cpm <- order(tot_486_3p_cpm)

pdf(file="Magen_etal_miR_486_3p_CPM_cnts.pdf",width=10,height=2.5)
par(bty="l")
par(mai=c(0.8,0.7,0.1,0.1))
par(mgp=c(2.6,0.6,0))
barplot(tot_486_3p_cpm[orderBy486Cpm],
        cex.axis=0.6,
        xlab="",
        ylab="miR 486-3p CPM",
        las=2,
        cex.names=0.3,
        cex.lab=0.8,
        border=NA,
        col=col_key_ALS_status[orderBy486Cpm])
dev.off()


#log transforming cpm +1 because VST gives different results depending on
#run with all samples first

minNonZeroCpm <- min(as.numeric(cpm_mat)[which(as.numeric(cpm_mat)>0)])

logCPMmat <- t(apply(cpm_mat,1,function(x){log(x+minNonZeroCpm)}))
dim(logCPMmat)

pca_log_cpm <- prcomp(t(logCPMmat))
head(pca_ord_vsd$x)

colnames(sample_key_df)

#add log cnts_colSums and log(486-3p cpm)

sample_key_df$logTotCnts <- log(cnts_colSums)

rownames(logCPMmat)[mask_486_3p]
sample_key_df$log486cpm <- logCPMmat[mask_486_3p,]

colnames(sample_key_df)

selectColnamesForPCexplore <- c("batchFac",
                                "sex",
                                "Disease",
                                "Onset",
                                "Riluzole",
                                "logTotCnts",
                                "log486cpm")

pc_explore_mat_logCpm <- explore_PCs(pca_log_cpm,cov_df=sample_key_df,
                                  test_colnames=selectColnamesForPCexplore,
                                    out_fig_dir="PC_corr_plots")

write.csv(data.frame(covar=rownames(pc_explore_mat_logCpm),pc_explore_mat_logCpm),
          file="Magen_et_al_plasma_logCPM_PC_explore_res.csv",row.names=F,quote=F)



#place cutoffs on total cnts and miR-486-3p

#a lot of very under-sequenced libs

#tot cnts
totCntsQs <- quantile(cnts_colSums)

#there an inflection point in plot fo tot cnts vs frac ALS at 1e6
maskPassCnts <- which(cnts_colSums>=5e5)


length(maskPassCnts)

tapply(sample_key_df$Disease[maskPassCnts],
       INDEX=sample_key_df$Disease[maskPassCnts],
       length)


#are log(cnts_colSums) and log(tot_486_3p_cpm) correlated

lm486vsTot <- lm(log(tot_486_3p_cpm) ~ log(cnts_colSums))
summary(lm486vsTot)

#very negatively correlated

#miR-486-3p
cpm486Qs <- quantile(tot_486_3p_cpm)

#seems to be an inflectin around the 75th percentile
maskPass486 <- which(tot_486_3p_cpm<=cpm486Qs[4])
length(maskPass486)

tapply(sample_key_df$Disease[maskPass486],
       INDEX=sample_key_df$Disease[maskPass486],
       length)

#now want to keep only batches that have samples passing from both ALS and cntrl
#AND have at least 1 of each sex

maskPassCntsAnd486 <- intersect(maskPassCnts,maskPass486)

uniqPassCntsand486Batches <- 
  as.character(unique(sample_key_df$batchFac[maskPassCntsAnd486]))


batchCatCnts <- lapply(uniqPassCntsand486Batches,function(batchX){
  subMaskX <- which(sample_key_df$batchFac[maskPassCntsAnd486] == batchX)
  
  diseaseX <- sample_key_df$Disease[maskPassCntsAnd486][subMaskX]
  sexX <- sample_key_df$sex[maskPassCntsAnd486][subMaskX]
  
  nALS <- length(which(diseaseX=="ALS"))
  nCtl <- length(which(diseaseX=="cntrl"))
  
  nM <- length(which(sexX=="M"))
  nF <- length(which(sexX=="F"))
  
  outVec <- c(nALS,nCtl,nM,nF)
  names(outVec) <- c("ALS","cntrl","M","F")
  return(outVec)
})


subMaskBatchHasBal <- which(unlist(lapply(batchCatCnts,function(x){
  if(is.null(x)){
    outVal <- 0
  }
  else{
    outVal <- min(x)
  }
  return(outVal)
}))>0)

batchCatCnts[subMaskBatchHasBal]

passingBalBatches <- uniqPassCntsand486Batches[subMaskBatchHasBal]

maskBatchBal <- which(sample_key_df$batchFac %in% passingBalBatches)


#keep samples passing tot cnts, 486 cnts, and balanced batches

maskPassAll <- intersect(maskBatchBal,
                         intersect(maskPassCnts,maskPass486))

tapply(sample_key_df$Disease[maskPassAll],
       INDEX=sample_key_df$Disease[maskPassAll],
       length)


#make pass all version of cnts mat, cpm mat, and sample key

pass_cnts_mat <- raw_cnts_mat[,maskPassAll]

pass_cpm_mat <- cpm_mat[,maskPassAll]


pass_sample_key <- sample_key_df[maskPassAll,]

pass_sample_key <- droplevels(pass_sample_key)

#filter smallRNA by median counts

medianMirCnts <- rowMedians(pass_cnts_mat)

medianMirCnts[grep("miR-206",rownames(pass_cnts_mat))]


pdf(file="Magen_et_al_plasma_medianMiRCnts_Hist.pdf",
    width=5,
    height=2.5)

par(bty="l")
par(mai=c(0.5,0.5,0.3,0.2))
par(mgp=c(1.4,0.6,0))

histOuts <- hist(log10(medianMirCnts),
                 cex.axis=0.6,
                 xlab="binned log10(cnts)",
                 ylab="n",
                 main ="small RNA median counts",
                 cex.lab=0.8,
                 breaks=25)

dev.off()

maskPassMedianCnts <- which(medianMirCnts >=50)


#find miRNAs that are all zero in any batch (these cannot be batch corrected)
minMaxMiRnaCntByBatch <- apply(pass_cnts_mat,1,function(rowX){
  cntsPerBatchX <- tapply(rowX,
                          INDEX=pass_sample_key$batchFac,
                          FUN=max)
  return(min(as.numeric(cntsPerBatchX)))
})

maskPassMinMaxBatch <- which(minMaxMiRnaCntByBatch>0)

#keep only miRNAs with at least 1 count in every sample
rowMinsPassCnts <- rowMins(pass_cnts_mat)

maskPassNoZeroCnts <- which(rowMinsPassCnts>0)

#maskPassAllMiRNA <- intersect(maskPassMedianCnts,maskPassMinMaxBatch)
maskPassAllMiRNA <- intersect(maskPassMedianCnts,maskPassNoZeroCnts)

pass_mirna_filt_cpm_mat <- pass_cpm_mat[maskPassAllMiRNA,]

pass_mirna_filt_logCPMmat <- apply(pass_mirna_filt_cpm_mat,2,log)

dim(pass_mirna_filt_cpm_mat)
dim(pass_mirna_filt_logCPMmat)

which(!is.finite(pass_mirna_filt_logCPMmat))

rowSDsPass_mirna_filt_logCPMmat <- rowSds(pass_mirna_filt_logCPMmat)
quantile(rowSDsPass_mirna_filt_logCPMmat)

#use ComBat batch normalization
dim(pass_mirna_filt_logCPMmat)

head(pass_sample_key)
dim(pass_sample_key)

modMatrixForComBat <- model.matrix(~.,data=pass_sample_key[,c(3,7)])

dim(modMatrixForComBat)

batCorrPassMirnaFiltLogCPMmat <- ComBat(dat=pass_mirna_filt_logCPMmat,
                                        batch=pass_sample_key$batchFac,
                                        mod=modMatrixForComBat,
                                        par.prior=FALSE,
                                        prior.plots=FALSE)


dim(batCorrPassMirnaFiltLogCPMmat)

head(batCorrPassMirnaFiltLogCPMmat[,1:20])


#re-run PCA with passing samples, passing miRNA, and batch corrected only

pca_pass_log_cpm <- prcomp(t(batCorrPassMirnaFiltLogCPMmat))

head(pca_pass_log_cpm$x[,1:10])


pc_explore_mat_PASS_logCpm <- explore_PCs(pca_pass_log_cpm,cov_df=pass_sample_key,
                                     test_colnames=selectColnamesForPCexplore,
                                     out_fig_dir="PC_PASSING_corr_plots")

write.csv(data.frame(covar=rownames(pc_explore_mat_PASS_logCpm),pc_explore_mat_PASS_logCpm),
          file="Magen_et_al_plasma_PASS_logCPM_PC_explore_res.csv",row.names=F,quote=F)


save.image(wrkspc_img)


#why not use residuals from DESeq2? From Michael Love himself:
#"There isn’t really a residual in a GLM that is similar to the residual in a linear model. I’d recommend working with the VST data if you want to assess the variation in the data before and after controlling for variables."

#make Df with variable to be controlled (in lm)
colnames(pass_sample_key)

#what is balance between sex in cases and controls?

tapply(pass_sample_key$sex,
       INDEX=pass_sample_key$Disease,
       function(x){tapply(x,x,length)})

controlledColnames <- c("sex",
                        "logTotCnts",
                        "log486cpm")

# controlledColnames <- c("sex")

maskKeptSampleCols <- which(colnames(pass_sample_key) %in% controlledColnames)

keptMetaLmDf <- pass_sample_key[,maskKeptSampleCols,drop=FALSE]
dim(keptMetaLmDf)
head(keptMetaLmDf)

keptLmDf <- data.frame(logCpm=I(t(batCorrPassMirnaFiltLogCPMmat)),keptMetaLmDf)


lmCovarFormula <- paste("logCpm",
                        paste(controlledColnames,collapse="+"),
                        sep="~")

lmCovarFormula <- as.formula(lmCovarFormula)

lmKeptVsdVsCovar <- lm(formula = lmCovarFormula,
                       data=keptLmDf)


#first get summaries and make a nice output table
summLmKeptVsdVsCovar <- summary(lmKeptVsdVsCovar)

length(summLmKeptVsdVsCovar)

listSummCoeffs <- lapply(summLmKeptVsdVsCovar,coefficients)
names(listSummCoeffs) <- rownames(batCorrPassMirnaFiltLogCPMmat)

#make results wide with each row a separate small RNA

listSummResRows <- lapply(listSummCoeffs,function(matX){
  #don't report intercept
  keepMatX <- matX[c(2:nrow(matX)),]
  
  #as.numeric() on the transpose t() of keepMatX makes a vector by rows
  allResX <- as.numeric(t(keepMatX))
  
  tempCoeffColnames <- colnames(keepMatX)
  tempCoeffColnames[4] <- "pValue"
  
  colnamesX <- paste(rep(rownames(keepMatX),each=4),tempCoeffColnames,sep="_")
  outMatX <- matrix(allResX,nrow=1)
  colnames(outMatX) <- colnamesX
  return(outMatX)
})

allSummResMat <- do.call(rbind,listSummResRows)

# allSummResDf <- data.frame(small_rna_row_data_df[maskKeepMedianCnts,],
#                            allSummResMat)

allSummResDf <- data.frame(miRNA=rownames(batCorrPassMirnaFiltLogCPMmat),
                           allSummResMat)


write.csv(allSummResDf,file="Magen_LM_regress_Summary_selected_Covars.csv",
          row.names=FALSE,
          quote=FALSE)

#now get residuals from the lm for each smallRNA

lmResids <- lmKeptVsdVsCovar$residuals
dim(lmResids)
lmResids[1:5,1:5]

#check if resids have 0 mean and unit variance
residMeans <- colMeans(lmResids)
quantile(residMeans)

residSDs <- colSds(lmResids)
head(residSDs[order(residSDs,decreasing=TRUE)],10)
quantile(residSDs)

residSDs[grep("miR-206",colnames(lmResids))]


#they do not have unit variance
#scale them
scaledResids <- scale(lmResids)
dim(scaledResids)
quantile(colMeans(scaledResids))
quantile(colSds(scaledResids))

save.image(wrkspc_img)

#make heatmap of scaled resids to look for structure
#try with all
head(pass_sample_key)

colByDisease <- rep("blue",nrow(scaledResids))
colByDisease[which(pass_sample_key$Disease == "ALS")] <- "red"

png(file="Magen_scaled_Covar_Resids_All_heatmap.png",
    width=10,
    height=8,
    units="in",
    res=1000,
    type="cairo")

heatmap(scaledResids,
        labRow=pass_sample_key$Disease,
        labCol=colnames(scaledResids),
        margins = c(10,5),
        RowSideColors=colByDisease)
dev.off()



#lm resids to see predictors of disease

dataLmResidsVsDiseaseDf <- data.frame(resids=I(scaledResids),
                                      Disease=pass_sample_key$Disease)

lmDiseaseVsResids <- lm(resids ~ Disease,data=dataLmResidsVsDiseaseDf)

summLmDiseaseVsResids <- summary(lmDiseaseVsResids)
length(summLmDiseaseVsResids)

listDiseaseVsResids <- lapply(summLmDiseaseVsResids,function(summX){
  coeffMatX <- coefficients(summX)
  outMatX <- coeffMatX[2,,drop=FALSE]
  colnamesX <- colnames(coeffMatX)
  colnamesX[4] <- "pValue"
  colnamesX <- paste(rownames(outMatX),colnamesX,sep="_")
  colnames(outMatX) <- colnamesX
  return(outMatX)
})

diseaseVsResidsLmResMat <- do.call(rbind,listDiseaseVsResids)
dim(diseaseVsResidsLmResMat)

diseaseVsResidsLmResDf <- data.frame(miRNA=colnames(scaledResids),
                                     medianMirCnts=medianMirCnts[maskPassAllMiRNA],
                                     diseaseVsResidsLmResMat)

head(diseaseVsResidsLmResDf)

write.csv(diseaseVsResidsLmResDf,
          file="Magen_LM_regress_Summary_Resids_ALS_vs_cntrl.csv",
          row.names=FALSE,
          quote=FALSE)

diseaseVsResidsLmResDf[grep("miR-206",diseaseVsResidsLmResDf$miRNA),]


save.image(wrkspc_img)


#now run LASSO on scaled resids with leave one out cross validation
set.seed(1990)

head(pass_sample_key)

ALSRespNumeric <- rep(0,nrow(pass_sample_key))
ALSRespNumeric[which(pass_sample_key$Disease == "ALS")] <- 1

#weight samples incidence
ALSweight <- length(which(ALSRespNumeric==0))/length(which(ALSRespNumeric==1))

weightsALS <- rep(1,length(ALSRespNumeric))
weightsALS[which(ALSRespNumeric==0)] <- ALSweight

#run glmnet with cross
cv_glmnet_als <- cv.glmnet(x= scaledResids,
                           y= ALSRespNumeric,
                           folds=nrow(scaledResids),
                           #folds=10,
                           maxit=1e5,
                           family="binomial",
                           keep=TRUE,
                           type.measure="deviance",
                           relax=FALSE,
                           weights=weightsALS,
                           parallel=TRUE)


png(filename = "Magen_etal_ALS_LmResids_lasso_plot.png",
    width=5,height=3,units="in",res=1000,type="cairo")
par(mai=c(0.6,0.6,0.3,0.2))
par(bg="white")
par(bty='l')
par(mgp=c(1.7,0.6,0))

plot(cv_glmnet_als)
dev.off()

#get coeffs of lambda.min and lambda.1se
cv_glmnet_als$lambda.1se
cv_glmnet_als$lambda.min

cv_glmnet_als$lambda

maskPicked <- 7

temp_coeff <- as.matrix(coefficients(cv_glmnet_als,s= cv_glmnet_als$lambda[maskPicked]))
raw_coeff <- temp_coeff[which(abs(temp_coeff[,1])>0),]
lamdaPickedlasso_coeff_df <- data.frame(feature=names(raw_coeff)[2:length(raw_coeff)],
                                        coeff=raw_coeff[2:length(raw_coeff)],
                                        pickedLambda=cv_glmnet_als$lambda[maskPicked])

write.csv(lamdaPickedlasso_coeff_df,file="Magen_etal_ALS_AllSamples_resids_lambda_picked_features.csv",row.names=F,quote=F)

#try small lambda
save.image(wrkspc_img)

##
source("/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/analysisScripts/jitterBoxplotXs.R")
source("/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/analysisScripts/changeColorAlpha.R")

benColPal <- c("#8C510A","#D8B365","#F6E8C3","#C7EAE5","#5AB4AC","#01665E")



plotVals <- scaledResids[,which(colnames(scaledResids)=="hsa-miR-206")]

alsFac <- pass_sample_key$Disease

colKeyAls <- rep(changeColorAlpha(benColPal[3],newAlpha=250),
                 length(alsFac))
colKeyAls[which(alsFac=="ALS")] <- changeColorAlpha(benColPal[6],newAlpha=250)

globalMin <- min(plotVals)
globalMax <- max(plotVals)

globalYlim <- c(globalMin,globalMax)

jitterXs <- jitterBoxplotsXs(plotVals,
                             plotFactor=alsFac,
                             maxJitter=0.2)

groupMedians <- tapply(plotVals,
                       INDEX=alsFac,
                       FUN=median)

fileX <- paste("magen_etal_MiR206_resids_boxplot.pdf",sep="_")

pdf(fileX,width=1.5,height=2)

par(mai=c(0.4,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=jitterXs,
     y=plotVals,
     pch=1,
     col=colKeyAls,
     cex=1,
     xlab="",
     xlim=c(0.7,2.3),
     xaxt='n',
     ylab="miR-206 residuals",
     ylim=globalYlim)

lines(x=c(0.7,1.3),y=rep(groupMedians[1],2),col="black",lwd=1.5)
lines(x=c(1.7,2.3),y=rep(groupMedians[2],2),col="black",lwd=1.5)

axis(side=1,at=c(1,2),labels = FALSE)

par(xpd=NA)

plotLims <- par('usr')

text(x=c(1,2),y=c(plotLims[3],plotLims[3]),
     labels=c("Cntrl","ALS"),
     pos=1,offset=1,
     cex=0.8)

dev.off()


save.image(wrkspc_img)
