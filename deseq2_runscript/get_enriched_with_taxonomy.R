# get table with OTU_ID and ConsensusLineage
otu_table <- read.table("~/Dropbox/fungal_microbiome/combined/combined_onlyfungi.csv", sep = ",", header = TRUE)
head(otu_table, 1)
only_cl <- otu_table[ , c("OTU_ID", "ConsensusLineage")]
head(only_cl, 1)

# get enriched endo otus in table
# MAKE SURE OTU_ID COLUMN IS ACTUALLY NAMED "OTU_ID"
enriched_endo <- read.table("enriched_endo.csv", sep = ",", header = TRUE)
head(enriched_endo, 1)
colnames(enriched_endo)[1] <- "OTU_ID"
head(enriched_endo, 1)

# get enriched bulk otus in table
# MAKE SURE OTU_ID COLUMN IS ACTUALLY NAMED "OTU_ID"
enriched_bulk <- read.table("enriched_bulk.csv", sep = ",", header = TRUE)
head(enriched_bulk, 1)
colnames(enriched_bulk)[1] <- "OTU_ID"
head(enriched_bulk, 1)

# MERGE
merged <- merge(only_cl, enriched_bulk, by = "OTU_ID")
head(merged, 1)

# make sure no enriched OTUs lost in merge
all(enriched_endo$OTU_ID %in% merged$OTU_ID)

# make sure column names are appropriate
colnames(merged) 

# write files
write.csv(as.data.frame(merged), file = "enriched_endo_CL.csv", row.names = FALSE)