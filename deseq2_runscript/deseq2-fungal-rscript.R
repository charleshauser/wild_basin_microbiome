# DESeq2 script for WBCRC fungal microbiome count data (no bacteria)

# clear workspace
rm(list = ls())

library("DESeq2")

setwd("~/Desktop/deseq2_runscript")

# SET UP COUNTS MATRIX
# read counts file (no taxonomy)
df <- read.table("combined_onlyfungi_noCL", sep = ",", row.names = "OTU_ID",
				 header = TRUE)
head(df, 1)

# sort df by col name
df_ordered <- df[ , order(names(df))]
head(df_ordered, 1)

# convert df to numeric matrix
countData <- data.matrix(df_ordered)
dim(countData) # ADD DIMENSIONS
head(countData, 2)

# pre-filter
# remove rows which only 1 count or none
countData <- countData[rowSums(countData) > 1, ]
head(countData, 1)

# convert 0 -> 1 to avoid log errors
countData[countData == 0] <- 1
head(countData)

# SET UP METADATA MATRIX
# read in sample metadata (includes sampleID, scientificName, commonName, 
# yearCollected, and soilFraction)
colData <- read.csv("sample_metadata", sep = ",", row.names = "sampleID", 
					header = TRUE)
head(colData)

# order by row name so column of count matrix and rows of column data are in
# the same order
colData <- colData[order(names(colData)), ]
head(colData)

#verify all sample names are present in both files
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

# CONSTRUCT OBJECT AND RUN PIPELINE
# construct dds object using soilFraction as experiment design
dds <- DESeqDataSetFromMatrix(countData = countData,
							  colData = colData,
							  design = ~ soilFraction)
dds
# run pipeline
dds <- DESeq(dds)

####################### ANALYSIS #######################

# PCA ANALYSIS

# log transform
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
	 pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)

# assess overall similarity between samples
# "use the R function dist to calculate the Euclidean distance between samples.
#  To ensure we have a roughly equal contribution from all genes, we use it on
#  the rlog-transformed data. We need to transpose the matrix of values using t, 
#  because the dist function expects the different samples to be rows of its 
#  argument, and different dimensions (here, genes) to be columns"
sampleDists <- dist( t( assay(rld) ) )
sampleDists

library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists)
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png("allSx_heatmap.png")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# poisson Distance dissimilarity heatmap
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
png("poissonDissimilarity_heatmap.png")
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)
dev.off()


# PCA plot
png("PCA_allSx.png")
plotPCA(rld, intgroup = c("yearCollected", "soilFraction"))
dev.off()

(pcaData <- plotPCA(rld, intgroup = c("yearCollected", "soilFraction"), returnData=TRUE))
percentVar <- round(100 * attr(pcaData, "percentVar"))
library("ggplot2")

# add additonal emtadata to colData file for use in PCA descrimination
ggplot(pcaData, aes(PC1, PC2, color=yearCollected,shape=soilFraction)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()


# MSDS PLot
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=yearCollected)) + geom_point(size=3) +
    coord_fixed()

mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds)))
png("mdsPoisData.png")
# FOLLOWING COMMAND DOES NOT WORK; "CONDITION" NOT DEFINED FOR MICROBIOME DATA
# ggplot(mdsPois, aes(X1,X2,color=condition)) + geom_point(size=3) + 
#        coord_fixed()


# DIFFERENTIAL TESTING
# For differential testing we recommend the DESeq function applied to raw 
# counts (NOT log transformed)
dds <- DESeq(dds)
res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)

# compare any 2 levels contrast (variable, numerator, denominator)
# FIXED ERROR HERE: variable must be design of dds; numerator and denominator
# must be from resultsNames(dds) with variable substring removed
# EX: variable = hello, resultsNames(dds) = [hello1, hello2], 
#     contrast = c(hello, 1, 2)
resCe_Ctl <-results(dds, contrast=c("soilFraction", "E", "B"))
sum(resCe_Ctl$padj < 0.1, na.rm=TRUE)#53
resSig <- subset(resCe_Ctl, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])

write.csv(as.data.frame(resCe_Ctl),file="Ce_vs_CTL_results.csv")

sessionInfo()