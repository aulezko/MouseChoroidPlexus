#Alina Ulezko Antonova
#May 14th, 2023
#Pathway Analysis Radar Chart

#Library Load
library(readxl)
library(pcaExplorer) #for topGOtable, v.2.20.2
library(GeneTonic) #for shake_topGOtableResult, v.1.6.4
library(fmsb) #for radarchart, v.0.7.1

#List of top DEGs calculated with Biotouring
k <- read_excel("/storage1/fs1/mcolonna/Active/Alina/LabCollabs/Simone/analysis/Gene_list.xlsx")

#Extract gene symbols of top DEGs per cluster that were higher than 0.1 log2FC enriched 
#and then perform pathway enrichment analysis with pcaExplorer, using the GO database
#MHC II
genes_mhcii <- k$Gene[k$Cluster == "MHCII" & k$`Log2(FC)` > 0.1]
topgoDE_mhcii_vs_others <-
  pcaExplorer::topGOtable(genes_mhcii,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Mm.eg.db",
                          geneID = "symbol",
                          topTablerows = 10000)
res_enrich_mhcii <- shake_topGOtableResult(topgoDE_mhcii_vs_others)
#CD9
genes_cd9 <- k$Gene[k$Cluster == "CD9" & k$`Log2(FC)` > 0.1]
topgoDE_cd9_vs_others <-
  pcaExplorer::topGOtable(genes_cd9,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Mm.eg.db",
                          geneID = "symbol",
                          topTablerows = 10000)
res_enrich_cd9 <- shake_topGOtableResult(topgoDE_cd9_vs_others)
#CD163
genes_cd163 <- k$Gene[k$Cluster == "CD163" & k$`Log2(FC)` > 0.1]
topgoDE_cd163_vs_others <-
  pcaExplorer::topGOtable(genes_cd163,
                          bg_ids,
                          ontology = "BP",
                          mapping = "org.Mm.eg.db",
                          geneID = "symbol",
                          topTablerows = 10000)
res_enrich_cd163 <- shake_topGOtableResult(topgoDE_cd163_vs_others)

#IN EXCEL: From these tables, extract the p value of top enriched pathways 
#(cross-comparing these results with those of metascape manually), and then order 
#the columns such that the order of the pathways matches the order of the clusters in the rows. 
#Then, transform these p values into -log10(pvalue).

#Make the radar plot:
p <- read_excel("/storage1/fs1/mcolonna/Active/Alina/LabCollabs/Simone/analysis/pathways/pathways_radar.xlsx")
#make clusters be rownames
names <- rownames(p)
p$...1 <- NULL
rownames(p) <- names
#Scale pathways by cluster, so that maximum is 2 (arbitrary)
# Define a scaling function
scale_max_2 <- function(x) {
  (x / max(x)) * 2
}
# Apply the scaling function to each column
p1 <- as.data.frame(apply(p, 2, scale_max_2))
#Add the 2 rows that are needed for the radarchart() function to work. The first row are 0s and the second row are 3s. 
#These numbers were chosen based on a balanced visualization of the area occupied by the radars
new_row <- as.data.frame(t(rep(0, ncol(p1))))
# Set the names of the new row to match the names of df
names(new_row) <- names(p1)
# Insert the new row at the top of the dataframe
p1 <- rbind(new_row, p1)
new_row <- as.data.frame(t(rep(3, ncol(p1))))
# Set the names of the new row to match the names of df
names(new_row) <- names(p1)
# Insert the new row at the top of the dataframe
p1 <- rbind(new_row, p1)
#see how a simple radarchart would look
radarchart(p1)

#customize the chart and save as pdf
dev.off()
pdf("radar_new.pdf")
radarchart( p1  , axistype=1, 
            #custom polygon
            pcol=c("#DE3D26",#r
                   "#4166B0", #b
                   "#159747") , pfcol=c("#DE3D2633",
                                        "#4166B033",
                                        "#15974733") , 
            plwd=3, #width of line
            plty=1, # continuous line
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,4,0.5), cglwd=0.8,
            #custom labels
            vlcex=0.8)
dev.off()


