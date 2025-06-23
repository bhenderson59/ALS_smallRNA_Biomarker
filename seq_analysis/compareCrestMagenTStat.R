#srun -c 4 -p interactive --pty /bin/bash
#cd /cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/compareCrestMagenSummStats
#module load cluster/R/4.4.1
#R --vanilla


library(GenomicRanges)

#bring in Magen ALS vs cntrl regression summary stats

magenSumStatFile <- 
  "/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/Magen_etal_ALS_small_RNA/Magen_LM_regress_Summary_Resids_ALS_vs_cntrl.csv"

magenSumStatDF <- read.csv(magenSumStatFile)

#bring in Crestwood ALS vs cntrl regression summary stats

creatSumStatFile <- 
  "/cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/CrestwoodSeqAnalysisClean/Crestwood_LM_regress_Summary_Resids_ALS_vs_cntrl.csv"

crestSumStatDF <- read.csv(creatSumStatFile)

#convert to comparable miRNA names, find matches and order

#magen
head(magenSumStatDF)

cleanMagenNames <- gsub("-","_",magenSumStatDF$miRNA)
head(cleanMagenNames)

#crestwood

head(crestSumStatDF)

cleanCrestNames <- sub("_[^_]*$", "", crestSumStatDF$miRNA)

head(cleanCrestNames)

#find names common to both

commonMiRNANames <- sort(intersect(cleanMagenNames,cleanCrestNames))

magenToCommonHits <- findMatches(commonMiRNANames,cleanMagenNames)

crestToCommonHits <- findMatches(commonMiRNANames,cleanCrestNames)


commonSummStatsDf <- 
  data.frame(miRNA=commonMiRNANames,
             magenTstat=magenSumStatDF$DiseaseALS_t.value[subjectHits(magenToCommonHits)],
             magenPval= magenSumStatDF$DiseaseALS_pValue[subjectHits(magenToCommonHits)],
             crestTstat=crestSumStatDF$DiseaseALS_t.value[subjectHits(crestToCommonHits)],
             crestPval= crestSumStatDF$DiseaseALS_pValue[subjectHits(crestToCommonHits)])


globalMin <- min(c(commonSummStatsDf$magenTstat,
                   commonSummStatsDf$crestTstat))

globalMax <- max(c(commonSummStatsDf$magenTstat,
                   commonSummStatsDf$crestTstat))

globalPlotLims <- c(globalMin,globalMax)


benBigColPal <- c("#543005","#8C510A","#BF812D","#DFC27D","#F6E8C3","#C7EAE5",
                  "#80CDC1","#35978F","#01665E","#003C30")


fileX <- paste("CrestVsMagenTstatPlot.pdf",sep="_")
pdf(fileX,width=2,height=2)

par(mai=c(0.43,0.43,0.1,0.07))
par(mgp=c(1.2,0.4,0))
par(bty="l")
par(cex.axis=0.7)
par(cex.lab=0.8)
par(xpd=NA)

plot(x=commonSummStatsDf$magenTstat,
     y=commonSummStatsDf$crestTstat,
     xlim=globalPlotLims,
     ylim=globalPlotLims,
     xlab=expression("Magen et al. " * italic(t) * "-stat."),
     ylab=expression("This Study " * italic(t) * "-stat."),
     col=benBigColPal[9],
     type='p',
     pch=1,
     cex=0.8)

mask206 <- which(commonMiRNANames=="hsa_miR_206")

points(x=commonSummStatsDf$magenTstat[mask206],
       y=commonSummStatsDf$crestTstat[mask206],
       pch=16,
       cex=0.8,
       col=benBigColPal[9])

dev.off()

dim(commonSummStatsDf)
#187   5
#
commonSummStatsDf[which(commonSummStatsDf$magenTstat>3),]

commonSummStatsDf[intersect(which(commonSummStatsDf$magenTstat>1),
                            which(commonSummStatsDf$crestTstat>1)),]



