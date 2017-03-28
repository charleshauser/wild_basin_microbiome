# GET SIGNIFICANT OTUs
#----------------------
# load results
results <- read.table("Ce_vs_CTL_results.csv", sep = ",", header = TRUE)

# drop all insignificant rows (padj > 0.05)
significant_results <- results[results$padj <= 0.05 & !is.na(results$padj), ]
head(significant_results, 1)

# make sure OTU colum is named "OTU_ID"
colnames(significant_results)[1] <- "OTU_ID"
head(significant_results, 1)

# get enriched endo (significant results where log2FoldChange > 0)
enriched_endo <- significant_results[significant_results$log2FoldChange > 0, ]
head(enriched_endo, 1)
dim(enriched_endo)

# get enriched bulk (significant results where log2FoldChange < 0)
enriched_bulk <- significant_results[significant_results$log2FoldChange < 0, ]
head(enriched_bulk, 1)
dim(enriched_bulk)

# ADDING PHYLOGENY
#------------------
# get table with OTU_ID and ConsensusLineage
otu_table <- read.table("combined_onlyfungi.csv", sep = ",", header = TRUE)
head(otu_table, 1)
only_cl <- otu_table[ , c("OTU_ID", "ConsensusLineage")]
head(only_cl, 1)

# merge to add 
merged_endo <- merge(only_cl, enriched_endo, by = "OTU_ID")
head(merged_endo, 1)
dim(merged_endo)

merged_bulk <- merge(only_cl, enriched_bulk, by = "OTU_ID")
head(merged_bulk)
dim(merged_bulk)

# make sure no enriched OTUs lost in merge
all(enriched_endo$OTU_ID %in% merged_endo$OTU_ID)
all(enriched_bulk$OTU_ID %in% merged_bulk$OTU_ID)

# make sure column names are appropriate
colnames(merged_endo)
colnames(merged_bulk) 

# write files
write.csv(as.data.frame(merged_endo), file = "enriched_endo.csv", row.names = FALSE)
write.csv(as.data.frame(merged_bulk), file = "enriched_bulk.csv", row.names = FALSE)

sessionInfo()