#RNA-Seq Analysis of Vaped CBD in HBECs

#setting working directory and loading packages
setwd("/Users/charlotte/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Lab/CBD and D8 Project/RNA-Seq/Vaped CBD RNA-Seq/Love Data/My Analysis")

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
cts <- read.csv('rawcounts.csv', header = TRUE, sep = ",")
coldata <- read.csv('coldata.csv', header = TRUE, sep = ",")

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
dds$condition <- relevel(dds$condition, ref = "Apical_CBD")

#differential expression analysis for all groups combined 
dds <- DESeq(dds)

#transformed values
vsd <- vst(dds, blind=FALSE) #variance stabilizing transformation
rld <- rlog(dds, blind=FALSE) #regularized log transformation
ntd <- normTransform(dds) #normalized counts transformation
write.csv(ntd, file="VapedCBD_normalizedcounts.csv")
head(assay(vsd), 3)

#PCA of all groups
tiff(filename = "VapedCBD_PCA.tiff", width = 2000, height = 2000, res = 300)
nudge <- position_nudge(y = 2)
plotPCA(vsd, intgroup=c("condition")) +
  geom_label(aes(label = coldata$sample), position = nudge, size = 2.5) +
  coord_fixed(xlim = c(-25,20))
dev.off()

pcaData <- plotPCA(vsd, intgroup = "condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$condition <- factor(pcaData$condition, levels=c("Apical_CBD", "CBD_Juice_20Puff", "CBD_Oil_20Puff"), labels=c("Apical_CBD", "CBD_Juice_20Puff", "CBD_Oil_20Puff"))

#formatted PCA
tiff(filename = "VapedCBD_PCA.tiff", width = 1500, height = 1500, res = 300)
pca <- ggplot(pcaData, aes(PC1, PC2))
pca + geom_point(size=3, aes(color=condition, fill=condition, group=condition)) + 
  ggtitle("Principle Component Analysis 12 hours") + 
  xlab(paste0("PC1: ", percentVar[1],"% varience")) + 
  ylab(paste0("PC2: ", percentVar[2],"% varience")) + 
  theme(legend.title=element_blank()) + 
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("lightskyblue", "plum2", "peachpuff"), labels = c("Non-vaped CBD Juice", "20 Puffs CBD Juice", "20 Puffs CBD Distillate")) +
  coord_fixed() 
dev.off()

#Getting human gene symbols
library("org.Hs.eg.db") #downloading human database, need to change if different species
columns(org.Hs.eg.db)

#Extracting results for each separate comparison
#comparison 1 (20 Puff Juice vs Apical)
res1 <- results(dds, contrast=c("condition","CBD_Juice_20Puff","Apical_CBD")) #reference or control should be second
resOrdered1 <- res1[order(res1$padj),] #ordering by adjusted p-value
res.df1 <- as.data.frame(resOrdered1)
res.df1$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res.df1), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res.df1, file="20PuffJuicevsApical.csv") #writing results file for this specific comparison
res05_1 <- results(dds, contrast=c("condition","CBD_Juice_20Puff","Apical_CBD"), alpha=0.05)
summary(res05_1) #use this output to see an overview of results

#comparison 2 (20 Puff Oil vs Apical)
res2 <- results(dds, contrast=c("condition","CBD_Oil_20Puff","Apical_CBD")) #reference or control should be second
resOrdered2 <- res2[order(res2$padj),] #ordering by adjusted p-value
res.df2 <- as.data.frame(resOrdered2)
res.df2$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res.df2), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res.df2, file="20PuffOilvsApical.csv") #writing results file for this specific comparison
res05_2 <- results(dds, contrast=c("condition","CBD_Oil_20Puff","Apical_CBD"), alpha=0.05)
summary(res05_2) #use this output to see an overview of results

#comparison 3 (20 Puff Oil vs 20 Puff Juice)
res3 <- results(dds, contrast=c("condition","CBD_Oil_20Puff","CBD_Juice_20Puff")) #reference or control should be second
resOrdered3 <- res3[order(res3$padj),] #ordering by adjusted p-value
res.df3 <- as.data.frame(resOrdered3)
res.df3$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res.df3), keytype = "ENSEMBL", column = "SYMBOL") #adding gene names
write.csv(res.df3, file="20PuffOilvs20PuffJuice.csv") #writing results file for this specific comparison
res05_3 <- results(dds, contrast=c("condition","CBD_Oil_20Puff","CBD_Juice_20Puff"), alpha=0.05)
summary(res05_3) #use this output to see an overview of results

#Making results file of only significant genes for each comparison
#Then making a vector of significant gene probe IDs to use in heatmap
res_sig1 <- res1[which(res1$padj < 0.05 & abs(res1$log2FoldChange) >= 1), ] 
summary(res_sig1)
write.csv(res_sig1, file="JuicevsApical_sig.csv")
sig_JuicevsApical <- rownames(res_sig1) 

res_sig2 <- res2[which(res2$padj < 0.05 & abs(res2$log2FoldChange) >= 1), ] 
summary(res_sig2)
write.csv(res_sig2, file="OilvsApical_sig.csv")
sig_OilvsApical <- rownames(res_sig2)

res_sig3 <- res3[which(res3$padj < 0.05 & abs(res3$log2FoldChange) >= 1), ] 
summary(res_sig3)
write.csv(res_sig3, file="OilvsJuice_sig.csv")
sig_OilvsJuice <- rownames(res_sig3) 

#Combining vectors of significant genes for each comparison and then removing duplicate probe IDs
all_sig_genes <- c(sig_JuicevsApical, sig_OilvsApical, sig_OilvsJuice)
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
donors <- c("37P", "41N", "41P", "69N")
donor_colors <- setNames(RColorBrewer::brewer.pal(4, "Set2"), donors)

tiff(filename = "Sig_Genes_AllGroups_HM.tiff", width = 2500, height = 2000, res = 300)
Heatmap(hmtable,
        #column_dend_reorder = FALSE,
        heatmap_legend_param = list(title = "Z-score"),
        column_names_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        top_annotation = HeatmapAnnotation(
                  Group = coldata$condition, 
                  Donor = coldata$donor,
                  Sex = coldata$sex,
                  simple_anno_size = unit(.4, "cm"),
                  show_annotation_name = FALSE,  
                  col = list(Group = c("Apical_CBD" = "lightskyblue", "CBD_Juice_20Puff" = "plum2", "CBD_Oil_20Puff" = "peachpuff"), 
                              Donor = donor_colors,
                              Sex = c("F" = "lightpink", "M" = "lightblue")),
                  annotation_legend_param = list(Group = list(labels = c("Apical_CBD", "CBD_Juice_20Puff", "CBD_Oil_20Puff"), at = c("Apical_CBD", "CBD_Juice_20Puff", "CBD_Oil_20Puff")))))
dev.off()


#extracting gene lists for specific clusters
HM <- Heatmap(hmtable, km=4)  #in first heatmap, there we 3 main sub-clusters, split it into these 3 to obtain gene IDs for each clusster

tiff(filename = "HM_clusters.tiff", width = 2500, height = 2000, res = 300)
HM
dev.off()

HM <- draw(HM)  #Show the heatmap

#formatted heatmap with clusters
tiff(filename = "Vaped CBD 4 clusters HM.tiff", width = 2500, height = 2000, res = 300)
Heatmap(hmtable,
        km = 4,
        heatmap_legend_param = list(title = "Z-score"),
        column_names_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 6),
        show_row_names = FALSE,
        row_title = NULL,
        top_annotation = HeatmapAnnotation(
          Group = coldata$name, 
          simple_anno_size = unit(.4, "cm"),
          show_annotation_name = FALSE,  
          col = list(Group = c("Non-vaped CBD Juice" = "lightskyblue", "20 Puffs CBD Juice" = "plum2", "20 Puffs CBD Distillate" = "brown")), 
          annotation_legend_param = list(Group = list(labels = c("Non-vaped CBD Juice", "20 Puffs CBD Juice", "20 Puffs CBD Distillate"), at = c("Non-vaped CBD Juice", "20 Puffs CBD Juice", "20 Puffs CBD Distillate")))))
dev.off()

r.dend <- row_dend(HM)  #extracting row dendrogram
rcl.list <- row_order(HM)  #returns the order (numbered rownames) of each cluster specified above

lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters

library(magrittr) # needed to load the pipe function '%%'

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(GeneID = rownames(hmtable[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

clu_df$GeneName <- mapIds(org.Hs.eg.db, keys = clu_df$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
write.table(clu_df, file= "HM_gene_clusters.txt", sep="\t", quote=F, row.names=FALSE)


#Making Venn diagram
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")

library("ggVennDiagram")

#ggVennDiagram works for up to 7 comparisons, just add more comparisons to the list function below
#using the vectors of significant genes made previously
venn_list <- list("20 Puffs CBD Juice vs. Non-vaped CBD Juice" = sig_JuicevsApical, 
                  "20 Puffs CBD Dilstillate vs. Non-vaped CBD Juice" = sig_OilvsApical, 
                  "20 Puffs CBD Dilstillate vs. 20 Puffs CBD Juice" = sig_OilvsJuice)

tiff(filename = "Sig_Venn.tiff", width = 2500, height = 2000, res = 300)

ggVennDiagram(venn_list, category.names = c("Juice vs. Apical", "Oil vs. Apical", "Oil vs. Juice")) +
  scale_x_continuous(expand = expansion(mult = .1)) + #expanding x axis so labels aren't cut off 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

dev.off()

#making euler plot
library(gplots)
venn(venn_list)
my_venn <- venn(venn_list)
class(my_venn)
library(pryr)
names(attributes(my_venn))
venn_intersections <- attr(x = my_venn, "intersections")
vennDF <- print(as.data.frame(do.call(cbind, venn_intersections)))
write.csv(vennDF, file ="EulerIntersections_VapedCBD.csv")

library(eulerr)

fit1 <- euler(venn_list)
fit1
fit2 <- euler(venn_list, shape = "ellipse")
fit2
plot(fit2)

tiff(filename = "Sig_Euler.tiff", width = 2000, height = 2000, res = 300)
plot(fit2, 
     legend = list(side = "bottom"),
     quantities = list(type = c("counts")))
dev.off()

#getting gene lists from intersections 
intersect <- intersect(sig_JuicevsApical, sig_OilvsApical) #getting overlapping probe IDs
JuiceOnly <- setdiff(sig_JuicevsApical, sig_OilvsApical)
OilOnly <- setdiff(sig_OilvsApical, sig_JuicevsApical)

intersect <- as.data.frame(intersect) #making it a dataframe
colnames(intersect) = c("GeneID")
intersect$GeneName <- mapIds(org.Hs.eg.db, keys = intersect$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
intersect$GeneName
write.csv(intersect,file = "JuiceandOil_intersection.csv")

JuiceOnly <- as.data.frame(JuiceOnly) #making it a dataframe
colnames(JuiceOnly) = c("GeneID")
JuiceOnly$GeneName <- mapIds(org.Hs.eg.db, keys = JuiceOnly$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
JuiceOnly$GeneName
write.csv(JuiceOnly, file ="JuiceGenesOnly.csv")

OilOnly <- as.data.frame(OilOnly) #making it a dataframe
colnames(OilOnly) = c("GeneID")
OilOnly$GeneName <- mapIds(org.Hs.eg.db, keys = OilOnly$GeneID, keytype = "ENSEMBL", column = "SYMBOL") #adding on gene names
OilOnly$GeneName
write.csv(OilOnly,file = "OilGenesOnly.csv")

