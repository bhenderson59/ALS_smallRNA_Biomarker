#srun -c 4 -p interactive --pty /bin/bash
#cd /cluster/home/broberts/Ben_miRNA_work/ALSbiomarker/CrestTotCntsPrettyPlot
#module load cluster/R/4.4.1
#R --vanilla

require(matrixStats)
library(GenomicRanges)

#load in raw counts and name the columns from filenames that match the sample metadata
raw_cnts_base_dir <- "/cluster/home/bhenderson/ALS_biomarkers/Sequencing/ALS_SeqCenter_Groups1-3_pipeline/new_pipeline/smallRNA_alignments"

raw_cnts_files <- list.files(raw_cnts_base_dir,pattern="sRNA_feature_cnts.tsv",recursive=T)

short_file_bases <- unlist(lapply(strsplit(raw_cnts_files,split="_S[0-9]{1,2}_R1_001_trim"),function(split_x){
  #short_file_bases <- unlist(lapply(strsplit(raw_cnts_files,split="\\.consol"),function(split_x){
  return(split_x[1])
}))

short_file_bases <- gsub("-","_",short_file_bases)

raw_cnts_list <- lapply(raw_cnts_files,function(file_x){
  print(file_x)
  path_x <- paste(raw_cnts_base_dir,file_x,sep="/")
  raw_cnts_df_x <- read.table(path_x,sep="\t",header=T)
})

head(raw_cnts_list[[1]],20)

raw_cnts_mat <- do.call(cbind,lapply(raw_cnts_list,function(df_x){
  return(as.numeric(df_x$counts))
}))

small_rna_row_data_df <- raw_cnts_list[[1]][,1:5]

rownames(raw_cnts_mat) <- as.character(small_rna_row_data_df[,1])

colnames(raw_cnts_mat) <- short_file_bases

tail(raw_cnts_mat)

#bring in sample key metadata
sample_key_file <- "/cluster/home/bhenderson/ALS_biomarkers/analysis/CHA_ALS_seqcenter_sRNA_biomarkers_Group1-3_analysis/ALS_seqcenter_Groups1-3_biomarker_metadata.csv"
sample_key_df <- read.csv(sample_key_file)

head(sample_key_df)

sample_key_df$file_base <- short_file_bases
rownames(sample_key_df) <- short_file_bases
sample_key_df$Disease <- factor(sample_key_df$Disease,levels=c("CTL","ALS"))

cbind(sample_key_df$Sample_ID,sample_key_df$file_base)

#sample_key_df is ordered like raw_cnts_mat

#keep only W for Race because other categories are not powered
#convert "w", "M","E" (Egyptian) to "W" in Race column
sample_key_df$Race <- 
  as.factor(gsub("w|E|M","W",as.character(sample_key_df$Race)))

tapply(sample_key_df$Race,sample_key_df$Race,length)

maskKeepW <- which(sample_key_df$Race == "W")

sample_key_df <- sample_key_df[maskKeepW,]

raw_cnts_mat <- raw_cnts_mat[,maskKeepW]


#make batch a factor
batchFac <- sample_key_df$Library_Prep_Group
batchFac <- as.factor(batchFac)

sample_key_df$batchFac <- batchFac

head(sample_key_df)

#plot samples in order of total counts to look for min cnt cutoff

benBigColPal <- c("#543005","#8C510A","#BF812D","#DFC27D","#F6E8C3","#C7EAE5",
                  "#80CDC1","#35978F","#01665E","#003C30")



col_key_ALS_status <- rep("#F6E8C3",nrow(sample_key_df))
col_key_ALS_status[which(sample_key_df$Disease=="ALS")] <- "#01665E"

cnts_colSums <- colSums(raw_cnts_mat)

orderByCnts <- order(cnts_colSums)

min(cnts_colSums)

pdf(file="Crestwood_totMiRNA_cntsPRETTY.pdf",width=6.5,height=1.25)
par(bty="l")
par(mai=c(0.1,0.7,0.05,0.01))
par(mgp=c(2.5,0.7,0))
barplot(cnts_colSums[orderByCnts],
        cex.axis=0.7,
        xlab="",
        ylab="small RNA reads",
        las=2,
        cex.names=0.3,
        cex.lab=0.8,
        border=NA,
        col=col_key_ALS_status[orderByCnts])
dev.off()
