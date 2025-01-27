#RNA-seq analysis of 16HBEs exposed to 8.5µM CBD, CBDQ, or vehicle control for 12 or 24 hours

#setting working directory and loading packages
setwd("/Users/charlotte/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/CBD and D8 Project/RNA-Seq/CBD D8 RNA Seq/My Analysis")

library(DESeq2)
library(apeglm)
library(pheatmap)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(tibble)
library(dplyr)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(ggfortify)

#Loading in data and metadata
setwd("/Users/charlotte/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/CBD and D8 Project/RNA-Seq/CBD D8 RNA Seq/My Analysis")
cts <- read.csv('raw_counts_cbd_cbdq12.csv', header = TRUE, sep = ",")
coldata <- read.csv('coldata_3groups.csv', header = TRUE, sep = ",")

#checking that order of samples is the same in cts and coldata
head(cts)
coldata

#constructing a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition, tidy = TRUE)
dds

#pre-filtering genes
smallestGroupSize <- 3 #change this to the minimum number of samples you had in any of your treatment groups
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize #this requires genes to have 3 or more samples with counts of 10 or more. 
dds <- dds[keep,]

#specifying the reference level (control)
dds$condition <- relevel(dds$condition, ref = "VC")

#differential expression analysis for all groups combined 
dds <- DESeq(dds)

#transformed values
vsd <- vst(dds, blind=FALSE) #variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #regularized log transformation
ntd <- normTransform(dds) #normalized counts transformation
head(assay(vsd), 3)

#PCA of all groups
tiff(filename = "CBD_CBDQ_PCA.tiff", width = 1000, height = 2000, res = 300)
plotPCA(vsd, intgroup=c("condition")) #run this line only to visualize in R without saving
dev.off()

pcaData <- plotPCA(vsd, intgroup = "condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$condition <- factor(pcaData$condition, levels=c("VC", "CBD", "CBDQ"), labels=c("VC", "CBD", "CBDQ"))

#formatted PCA
tiff(filename = "CBD_CBDq_12hr_PCA.tiff", width = 1500, height = 1500, res = 300)
pca <- ggplot(pcaData, aes(PC1, PC2))
pca + geom_point(size=3, aes(color=condition, fill=condition, group=condition)) + 
  ggtitle("Principle Component Analysis 12 hours") + 
  xlab(paste0("PC1: ", percentVar[1],"% varience")) + 
  ylab(paste0("PC2: ", percentVar[2],"% varience")) + 
  theme(legend.title=element_blank()) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("VC" = "lightblue", "CBD" = "lightgreen", "CBDQ" = "tan1")) +
  coord_fixed() +
  ylim(-6, 4)
dev.off()

#Getting human gene symbols
library("org.Hs.eg.db") #downloading human database, need to change if different species
columns(org.Hs.eg.db)

#Extracting results for each separate comparison
#comparison 1 (CBD vs VC)
res1 <- results(dds, contrast=c("condition","CBD","VC")) #reference or control should be second
resOrdered1 <- res1[order(res1$padj),] #ordering by adjusted p-value
res.df1 <- as.data.frame(resOrdered1)
res.df1$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res.df1), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res.df1, file="CBDvsVC.csv") #writing results file for this specific comparison
res05_1 <- results(dds, contrast=c("condition","CBD","VC"), alpha=0.05)
summary(res05_1) #use this output to see an overview of results

#comparison 2 (CBDQ vs VC)
res2 <- results(dds, contrast=c("condition","CBDQ","VC")) #reference or control should be second
resOrdered2 <- res2[order(res2$padj),] #ordering by adjusted p-value
res.df2 <- as.data.frame(resOrdered2)
res.df2$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res.df2), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res.df2, file="CBDQvsVC.csv") #writing results file for this specific comparison
res05_2 <- results(dds, contrast=c("condition","CBDQ","VC"), alpha=0.05)
summary(res05_2) #use this output to see an overview of results

#comparison 3 (CBDQ vs CBD)
res3 <- results(dds, contrast=c("condition","CBDQ","CBD")) #reference or control should be second
resOrdered3 <- res3[order(res3$padj),] #ordering by adjusted p-value
res.df3 <- as.data.frame(resOrdered3)
res.df3$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res.df3), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res.df3, file="CBDQvsCBD.csv") #writing results file for this specific comparison
res05_3 <- results(dds, contrast=c("condition","CBDQ","CBD"), alpha=0.05)
summary(res05_3) #use this output to see an overview of results

#Making results file of only significant genes for each comparison
#Then making a vector of significant gene probe IDs to use in heatmap
res_sig1 <- res1[which(res1$padj < 0.05 & abs(res1$log2FoldChange) >= 1), ] 
sig_CBDvsVC <- rownames(res_sig1)
res_sig1$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_sig1), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(res_sig1, file="CBDvsVC_sig.csv")
summary(res_sig1)

res_sig2 <- res2[which(res2$padj < 0.05 & abs(res2$log2FoldChange) >= 1), ] 
sig_CBDQvsVC <- rownames(res_sig2)
res_sig2$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_sig2), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(res_sig2, file="CBDQvsVC_sig.csv")
summary(res_sig2)

res_sig3 <- res3[which(res3$padj < 0.05 & abs(res3$log2FoldChange) >= 1), ] 
sig_CBDQvsCBD <- rownames(res_sig3)
res_sig3$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res_sig3), keytype = "ENSEMBL", column = "SYMBOL")
write.csv(res_sig3, file="CBDQvsCBD_sig.csv")
summary(res_sig3)

#Combining vectors of significant genes for each comparison and then removing duplicate probe IDs
all_sig_genes <- c(sig_CBDvsVC, sig_CBDQvsVC, sig_CBDQvsCBD)
all_sig_genes_unique <- unique(all_sig_genes)

#Making matrix of normalized counts to use in heatmap
normcounts <- counts(dds, normalized=TRUE)
head(normcounts)

#filtering normalized counts by vector of significant genes previously made
norm_sigonly <- subset(normcounts, rownames(normcounts) %in% all_sig_genes_unique)

#Setting up to make heatmap
hmtable <- norm_sigonly
head(hmtable)
dim(hmtable)

#scaling data to obtain z-scores
hmtable_transposed <- t(hmtable)
head(hmtable_transposed)
hmtable_scaled <- scale(hmtable_transposed)
hmtable <- t(hmtable_scaled)

#plotting and saving heatmap 
tiff(filename = "Sig_Genes_AllGroups_HM.tiff", width = 2500, height = 2000, res = 300)
Heatmap(hmtable,
        #column_dend_reorder = FALSE,
        heatmap_legend_param = list(title = "Z-score"),
        column_names_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        top_annotation = HeatmapAnnotation(Group = coldata$condition, simple_anno_size = unit(.4, "cm"),
                                           show_annotation_name = FALSE,                                  
                                           col = list(Group = c("VC" = "lightblue", "CBD" = "lightgreen", "CBDQ" = "orange")),
                                           annotation_legend_param = list(Group = list(labels = c("VC", "CBD", "CBDQ"), at = c("VC", "CBD", "CBDQ")))))
dev.off()

#making Venn diagram
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")

library("ggVennDiagram")

#using the vectors of significant genes made previously
venn_list <- list(A = sig_CBDvsVC, 
                  B = sig_CBDQvsVC, 
                  C = sig_CBDQvsCBD)

tiff(filename = "CBD_CBDQ_Venn.tiff", width = 2000, height = 1000, res = 300)
ggVennDiagram(venn_list, category.names = c("CBD vs. VC", "CBDQ vs. VC", "CBDQ vs. CBD")) +
  scale_x_continuous(expand = expansion(mult = .1)) #expanding x axis so labels aren't cut off
dev.off()

#to get gene lists from certain intersections in Venn
intersect <- intersect(sig_CBDQvsVC, sig_CBDQvsCBD) #getting overlapping probe IDs
intersect <- as.data.frame(intersect) #making it a dataframe
colnames(intersect) = c("GeneID")
intersect$GeneName <- mapIds(org.Hs.eg.db, keys = intersect$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
intersect$GeneName
