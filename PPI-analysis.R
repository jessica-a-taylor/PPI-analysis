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

# Identify interaction hotspots.
interaction_counts <- data.frame(Gene = unique(c(NLR_PPIs$GeneID1, NLR_PPIs$GeneID2)),
                                 Counts = rep(0, times = length(unique(c(NLR_PPIs$GeneID1, NLR_PPIs$GeneID2)))))

for (gene in interaction_counts$Gene) {
  interaction_counts[interaction_counts$Gene==gene, "Counts"] <- length(c(which(NLR_PPIs$GeneID1==gene), 
                                                                          which(NLR_PPIs$GeneID2==gene)))
}

AT2G01023_network <- NLRs[c(which(NLRs$Gene %in% NLR_PPIs[c(which(NLR_PPIs$GeneID1=="AT2G01023"), 
                                                            which(NLR_PPIs$GeneID2=="AT2G01023")),"GeneID1"]),
                            which(NLRs$Gene %in% NLR_PPIs[c(which(NLR_PPIs$GeneID1=="AT2G01023"), 
                                                            which(NLR_PPIs$GeneID2=="AT2G01023")),"GeneID2"])),
                          c(1:4,7)]

AT3G41762_network <- NLRs[c(which(NLRs$Gene %in% NLR_PPIs[c(which(NLR_PPIs$GeneID1=="AT3G41762"), 
                                                            which(NLR_PPIs$GeneID2=="AT3G41762")),"GeneID1"]),
                            which(NLRs$Gene %in% NLR_PPIs[c(which(NLR_PPIs$GeneID1=="AT3G41762"), 
                                                            which(NLR_PPIs$GeneID2=="AT3G41762")),"GeneID2"])),
                          c(1:4,7)]

# Compare the response of interacting genes to flg22.
flg_exp <- as.data.frame(read_xlsx("Data/flg-induced expression.xlsx"))
AT2G01023_network_flg <- as.matrix(flg_exp[which(flg_exp$AGI %in% AT2G01023_network$Gene),c(2:6)])
rownames(AT2G01023_network_flg) <- flg_exp[which(flg_exp$AGI %in% AT2G01023_network$Gene),1]
heatmap(AT2G01023_network_flg, Colv = NA)

AT3G41762_network_flg <- as.matrix(flg_exp[which(flg_exp$AGI %in% AT3G41762_network$Gene),c(2:6)])
rownames(AT3G41762_network_flg) <- flg_exp[which(flg_exp$AGI %in% AT3G41762_network$Gene),1]
heatmap(AT3G41762_network_flg, Colv = NA)

#########################

# Filter for PPIs involving late genes.
late_up <- as.data.frame(read_csv("Data/Flg/Late_up_genes.csv"))
lateUp_PPIs <- data.frame()

for (row in 1:nrow(all_PPIs)) {
  if (all_PPIs[row, "GeneID1"] %in% late_up$Gene | all_PPIs[row, "GeneID2"] %in% late_up$Gene) {
    lateUp_PPIs <- rbind(lateUp_PPIs, all_PPIs[row,])
  }
}

# Identify interaction hotspots.
interaction_counts <- data.frame(Gene = unique(c(lateUp_PPIs$GeneID1, lateUp_PPIs$GeneID2)),
                                 Counts = rep(0, times = length(unique(c(lateUp_PPIs$GeneID1, lateUp_PPIs$GeneID2)))))

for (gene in interaction_counts$Gene) {
  interaction_counts[interaction_counts$Gene==gene, "Counts"] <- length(c(which(lateUp_PPIs$GeneID1==gene), 
                                                                          which(lateUp_PPIs$GeneID2==gene)))
}

# Produce a venn diagram showing the overlap in PPI networks between gene sets.

