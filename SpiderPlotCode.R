#Alina Ulezko Antonova
#May 14th, 2023
#Pathway Analysis Spider Plot

#Library Load
library(readxl)
library(triwise)
library(limma)
library(Biobase)
data(vandelaar)
library(Seurat)

#List of top DEGs calculated with Biotouring
k <- read_excel("/storage1/fs1/mcolonna/Active/Alina/LabCollabs/Simone/analysis/Gene_list.xlsx")


########## Macs annotation on Seurat object ######
macs <- readRDS("/storage1/fs1/mcolonna/Active/Simone/Cp_scRNAseq/Analysis_Combined_7.23.23/macs_v2_072323.rds")
# with resolution 0.5, cluster correspondence is as follows:
# 6 microglia
# 5 cd9 
# 8 mitotic
# 7,9 interferon
# 0, 1  mhc II
# 4 ccr2
# 10 meningeal
# 2, 3 cd163

#Annotate
Idents(macs) <- macs@meta.data$RNA_snn_res.0.5
macs$name <- plyr::mapvalues(
  x = macs$RNA_snn_res.0.5, 
  from = c("6", "5", "8", "7", "9", "0", "1", "4", "10", "2", "3"), 
    to = c("microglia", "cd9", "mitotic", 
           "ifn", "ifn", "mhcii", "mhcii", "ccr2", "meningeal", "cd163", "cd163")
)

#Take average expression of each gene in count Matrix per cluster, in the 3 macs clusters of interest
Idents(macs) <- macs@meta.data$name
macs_s <- subset(macs, idents =c("cd9", "mhcii", "cd163"))
macs_avg <- AverageExpression(macs_s, assay ="RNA", return.seurat = T)
macs_avg_counts <- macs_avg@assays$RNA@counts
macs_avg_counts_df <- as.data.frame(macs_avg_counts)
colnames(macs_avg_counts) <- c("MHCII", "CD163", "CD9")
barycoords = transformBarycentric(as.matrix(macs_avg_counts))


genes_to_label <- c("F13a1",
"H2-Eb1",
"Ctsd",
"Cd163",
"Cd74",
"Cd9",
"Fcrls",
"H2-Ab1",
"Trem2",
"Mrc1",
"H2-Aa",
"Hexb",
"Pf4",
"Apoe",
"Stab1",
"Itgb5",
"Tgfbr1",
"Cd63")
k <- plotDotplot(barycoords, Goi = genes_to_label, showlabels = T)
ggsave("/storage1/fs1/mcolonna/Active/Simone/Cp_scRNAseq/Analysis_Combined_7.23.23/dotplot_highlighted.png",k, dpi=300, width=8, height=8)
ggsave("/storage1/fs1/mcolonna/Active/Simone/Cp_scRNAseq/Analysis_Combined_7.23.23/dotplot_highlighted.svg",k, dpi=300, width=8, height=8)
ggsave("/storage1/fs1/mcolonna/Active/Simone/Cp_scRNAseq/Analysis_Combined_7.23.23/dotplot_highlighted.pdf",k, dpi=300, width=8, height=8)


