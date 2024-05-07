setwd("~/PPI-analysis")

# Load required libraries.
source("Functions for PPI analysis.R")
loadLibraries()
rm(loadLibraries)

# Import list of PPIs.
all_PPIs <- as.data.frame(read_xlsx("Data/PPI characterisations (Zhang et al., 2024).xlsx"))

# Filter for PPIs involving R-genes.
NLRs <- as.data.frame(read_xlsx("Data/Arabidopsis NLRs.xlsx"))
NLR_PPIs <- data.frame()

for (row in 1:nrow(all_PPIs)) {
  if (all_PPIs[row, "GeneID1"] %in% NLRs$Gene | all_PPIs[row, "GeneID2"] %in% NLRs$Gene) {
    NLR_PPIs <- rbind(NLR_PPIs, all_PPIs[row,])
  }
}

# Create a heatmap for the number of interactions per gene.
for (gene in unique(NLR_PPIs$GeneID1, NLR_PPIs$GeneID2)) {
  
}