install.packages("BiocManager")
install.packages("gprofiler2")
BiocManager::install("readxl") #to install the “readxl” package from Bioconductor
BiocManager::install("dendextend") #to install the “dendextend” package from Bioconductor
BiocManager::install("ggplot2") #to install the “ggplot2” package from Bioconductor
BiocManager::install("dplyr") #to install the “dplyr” package from Bioconductor
BiocManager::install("DESeq2") #to install the “DESeq2” package from Bioconductor
BiocManager::install("pheatmap") #to install the “pheatmap” package from Bioconductor
BiocManager::install("RColorBrewer") #to install the “RColorBrewer” package from Bioconductor
BiocManager::install("EnhancedVolcano") #to install the “EnhancedVolcano” package from Bioconductor
BiocManager::install("writexl") #to install the “writexl” package from Bioconductor

getwd()
setwd("C:/Users/Documents/Kofoworola")

countTable <- read.table("GSE147507_RawReadCounts_Human.tsv", header = T, sep = "\t",
                         row.names = 1)
colnames(countTable) #To print all column names
countTable <- countTable[, 1:6]
sampleInfo <- read_xlsx("sampleTable.xlsx")
rownames("sample")

# Data quality control
ddsMat <- DESeqDataSetFromMatrix(countData = countTable, colData = sampleinfo, design= ~virus)
ddsMat<- ddsMat[rowSums(counts(ddsMat))>1,]
rld <- rlog(ddsMat, blind = F)
sampleDist <- dist(t(assay(rld))) #Pairwise distance calculation
sampleMatrix <- as.matrix(sampleDist) #Converted to matrix
rownames(sampleMatrix)<-rld$sample #Changing row names of “sampleMatrix” to the values stored in rld$sample
colnames(sampleMatrix = NULL) #Setting the column names of “sampleMatrix” to NULL

# To identify trends, clusters or outliers in data
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
png("sampleHeatmap_NHBE.png", width=23, height=20, units="cm", res=500)
pheatmap(SampleMatrix, clustering_distance_rows =sampleDist, clustering_distance_cols=sampleDist, color=colors)
dev.off()

#To determine the principal components of your data
data <- plotPCA(rld, intgroup=c("virus"), returnData=T)
percentVar <- round(100*attr(data,"percentVar")) #Round values for easier interpretation
#Create plot
png("PCAplot_NHBE.png", width = 22, height = 20, units = "cm", res = 500)
ggplot(data, aes(PC1, PC2, shape=virus, color=virus)) +geom_point(size=3)+
  scale_color_manual(values = c("mock"="#9EA4BD", " SARS_CoV_2"="#424861"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ", percentVar[2],"% variance")) +
  theme(
    text = element_text(size = 20), # Set the overall font size to 20
    axis.title = element_text(size = 20), # Set the axis title font size to 20
    axis.text = element_text(size = 20), # Set the axis text font size to 20
    legend.text = element_text(size = 20), # Set the legend text font size to 20
    legend.title = element_text(size = 20) # Set the legend title font size to 20
  )
dev.off()




#To identify differentially expressed genes
dds <- DESeq(ddsMat)
resultsNames(dds)
resNHBE <- as.data.frame(results(dds, name = " virus_SARS_CoV_2_vs_mock"))
resNHBE <- filter(resNHBE, padj < 0.05) # Filter for significantly expressed genes only.
#Volcano plot to visualize the results of differential expression analysis.
png("Volcano Plot_NHBE.png")
EnhancedVolcano(resNHBE,
                lab = rownames(resNHBE),
                boxedLabels = F,
                x = "log2FoldChange", y = "padj", pCutoff = 0.05, FCcutoff = 0.5,
                col = c("grey", "grey", "grey", "#630436"), colAlpha = 0.5, labSize = 4,
                drawConnectors = F,
                title = "volcano plot of infected NHBE cells",
                subtitle = "differentially expressed genes") + coord_flip()
dev.off()
#splitting the gene into up and down regulated ones.
resNHBEUp <- rownames(filter(resNHBE, log2FoldChange > 0.5))
resNHBEDown <- rownames(filter(resNHBE, log2FoldChange < -0.5))
#determine overrepresented gene sets among the upregulated genes.
gsea <- gost(resNHBEUp, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaNHBEUp <- filter(gsea$result, term_size < 500, term_size > 10)
write_xlsx(gseaNHBEUp, path = "gseaNHBEup.xlsx") #Save result in excel sheet
#determine overrepresented gene sets among the downregulated genes.
gsea <- gost(resNHBEDown, organism = "hsapiens", sources = c("GO"), evcodes = T)
gseaNHBEDown <- filter(gsea$result, term_size < 500, term_size > 10)
write_xlsx(gseaNHBEDown, path = "gseaNHBEdown.xlsx") #save result in excel sheet
#display the 20 most significant, up-regulated GO terms.
ggplot(gseaNHBEUp[1:20,], aes(x = reorder(term_name, -log10(p_value)),
                                y = -log10(p_value))) +
  geom_bar(stat = "identity", fill = "orchid") + coord_flip() +
  xlab("GO terms") + ylab("-log10(p-Value)") + theme_minimal()
ggsave("BarPlot_upregulated.png", width = 15, height = 20, units = "cm", dpi = 500)
#display the 20 most significant, up-regulated GO terms.
ggplot(gseaNHBEDown[1:20,], aes(x = reorder(term_name, -log10(p_value)),
                              y = -log10(p_value))) +
  geom_bar(stat = "identity", fill = "orchid") + coord_flip() +
  xlab("GO terms") + ylab("-log10(p-Value)") + theme_minimal()
ggsave("BarPlot_downregulated.png", width = 15, height = 20, units = "cm", dpi = 500)