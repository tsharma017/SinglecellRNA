### Let's get all the package that will be necessary for our analysis
library(caret)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library()
# In this one we will do analysis of SHH vs Group3 -> Let's get the data
df <- read.csv('/Users/tejsharma/Downloads/Ct.csv')
Ann <- read.csv('/Users/tejsharma/Downloads/Sra.csv')
# for annotation converting that to factor
Ann$treatment <- factor(Ann$treatment)
rownames(Ann) <- Ann[,1]
Ann <- Ann[,-2]
#### Naming the Rows of df
#rownames(df) <- df[,-1]
#df=df[,-1]
df <- df_numeric <- df[sapply(df, is.numeric)]
#####No we gonna see the data distribution using boxplot
boxplot(df_numeric,outline=FALSE,col="cornflowerblue") ### data seems to have skewness on its one extremities
###Next thing we gonna do is take the top 10000 genes that has highest variability
tdf_numeric <- t(df)
SDs=apply(tdf,2,sd )
topPreds=order(SDs,decreasing = TRUE)[1:10000]
df_filtered=t(tdf[,topPreds])
df_filtered <- as.data.frame(df_filtered)
integer_cols <- sapply(df_filtered, is.integer)
df_filtered[!integer_cols] <- lapply(df_filtered[!integer_cols], as.integer)
####### done majority of filtering######################33333


###### We gonna work on data for further analysis#########333
dds <- DESeqDataSetFromMatrix(countData = df_filtered, colData = Ann, design = ~treatment)
dds$GROUP <- factor(dds$GROUP, levels = c('GROUP3', 'GROUP4'))
keeps <- rowSums(counts(dds)) >= 5 # we finna take genes with count total 5
dds <- dds[keeps,]
dds <- DESeq(dds)
plotDispEsts(dds)
## we  gonna make R_object Dataframe
deseq_result <- results(dds)
deseq_result <- as.data.frame(deseq_result)

#### Let's try to visualize this one on PCA plot
#first we gonna perform variance stabalization using VST
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup=c('GROUP'))


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
#set a color scheme
colors <- colorRampPalette(rev(brewer.pal(9,'Blues')))(255)


#generate the heatmap
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors) # In this plot we basically plotted distance between two samples

# Now plotting log transformation of top 50 genes
top_hits <- deseq_result[order(deseq_result$padj),][1:50,]
top_hits <- row.names(top_hits)
rld <- rlog(dds, blind = FALSE)


#pheatmap(assay(rld)[top_hits,], cluster_rows = FALSE,
#       show_rownames = FALSE, cluster_cols = FALSE)
pheatmap(assay(rld)[top_hits,],show_rownames = FALSE)
##Heatmap of Z scores. we will use the top 10 genes
normalized_counts <- counts(dds, normalized = TRUE)
cal_z_score <- function(x) (x-mean(x))/sd(x)
zscore_all <- t(apply(normalized_counts,1, cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap(zscore_subset,show_rownames = FALSE)

###########################################

## MA plot
plotMA(dds,ylim=c(-2,2))


#remove the noise
resLFc <- lfcShrink(dds,coef = 'GROUP_GROUP3_vs_SHH',
                    type='apeglm')


#Now generating MA plot after error removal
plotMA(resLFc,ylim=c(-2,2))

####################################
# Volcano plot
# Change resLFC to a dataframe
resLFc <- as.data.frame(resLFc)

# Adding label to the genes
resLFc$diffexpressed <- 'No'
resLFc$diffexpressed[resLFc$log2FoldChange > 1.5  & resLFc$padj < 0.05] <- 'UP'
resLFc$diffexpressed[resLFc$log2FoldChange < -1.5 & resLFc$padj < 0.05] <- 'DOWN'

# Annotate gene names for upregulated and downregulated genes
resLFc$delabel <- ifelse(abs(resLFc$log2FoldChange) > 1.5 & resLFc$padj < 0.05,
                         rownames(resLFc),
                         NA)

ggplot(data = resLFc, aes(x = log2FoldChange, y = -log10(padj),
                          col = diffexpressed, label = delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c('blue', 'black', 'red')) +
  theme(text = element_text(size = 20))
#####################################################################################################################################################
#############################################################################################################33
##############################################################################
#Comparing expressed genes
up_regulated <- rownames(resLFc[resLFc$diffexpressed=='UP',])
down_regulated <- rownames(resLFc[resLFc$diffexpressed=='DOWN',])
boxplot(df["ENSG00000118271.12",c("DMB006","RCMB38", "RCMB45",'RCMB51')])
boxplot(df["ENSG00000118271.12",c("MB002","MB511H", "RCMB40",'RCMB28')])
boxplot(c(13,156, 4744,30220))
boxplot(c(314,202,159,22))
