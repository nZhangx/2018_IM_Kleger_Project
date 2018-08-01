#R script

library(ensembldb)
library("EnsDb.Hsapiens.v86")
edb <- EnsDb.Hsapiens.v86
tx2gene <- transcripts(edb, columns = c(listColumns(edb , "tx"), "gene_name"), return.type = "DataFrame")
tx2gene <- tx2gene[,c(1,7)]

nm <- c("DP_1","DP_2", "EHB_1","EHB_2","L_1","L_2")
file <- file.path(".",paste0(nm,"quant.sf"))
all(file.exists(file))
names(file) <- nm

library(tximport)

#txi <- tximport(file, type = "salmon", tx2gene = tx2gene)

txi <- tximport(file, type = "salmon", abundanceCol = "TPM" , importer = function(x) read.csv2(x, sep = '\t', header = TRUE), tx2gene = tx2gene, countsCol = "NumReads", ignoreTxVersion = TRUE)


head(txi$counts)

df <- read.csv2("expDesign.csv", sep = '\t', header = T)
rownames(df) <- df$run
#txi$length <- txi$length[!is.na(txi$length)]
txi$length[is.na(txi$length)] <- 1
#txi$length[txi$length == 0] <- 1
#txi$counts[is.na(txi$counts)] <- 0

###### OPTIONAL TRIAL : normalize by quantile ######
 library(preprocessCore)
# txi$counts <- normalize.quantiles(txi$counts)




library(DESeq2)
df$batch <- as.factor(df$batch)
ddxi <- DESeqDataSetFromTximport(txi, colData = df, design =  ~ condition + batch)

##DESEQ analysis
png("dispest.png", width=1200, height = 900) #Plot dispersion 
plotDispEsts(dds)
dev.off()



####Prefiltering#####

row_sum_cutoff = 30
pvalue_cutoff = 0.05
lfc_cutoff = 1

keep <- rowSums(counts(ddxi)) >= row_sum_cutoff
ddxi <- ddxi[keep,]

dds <- DESeq(ddxi)

library(limma)
library(sva)

#dds <- removeBatchEffect(dds)

#PCA plot
vst.vals <-vst(dds, blind = FALSE)
vst.mat <- assay(vst.vals)
pca <- plotPCA(vst.vals, intgroup = colnames(vst.vals))
pdf("deseq_vst_pca_batched2.pdf")
pca
dev.off()


vst.norm <- removeBatchEffect(as.data.frame(vst.vals),batch = df$batch)
pca <- plotPCA(vst.vals, intgroup = colnames(vst.vals))
png("deseq_vst_pca_.png", width=1200, height = 900)
pca
dev.off()
######## individual analysis #########
# res <- results(dds, name="condition_celltype")
# res <- results(dds, contrast=c("condition","treated","untreated"))
# resnames <- resultsNames(dds)


res_DP_EHB <- results(dds, contrast = c("condition", "DP", "EHB"))
res_DP_EHB <- na.omit(res_DP_EHB)
res_DP_EHB.sig <- res_DP_EHB[res_DP_EHB$padj <= pvalue_cutoff,]
res_DP_EHB.sig <- res_DP_EHB.sig[abs(res_DP_EHB.sig$log2FoldChange) >= lfc_cutoff,]


res_EHB_L <- results(dds, contrast = c("condition", "L", "EHB"))
res_EHB_L <- na.omit(res_EHB_L)
res_EHB_L.sig <- res_EHB_L[res_EHB_L$padj <= pvalue_cutoff,]
res_EHB_L.sig <- res_EHB_L.sig[abs(res_EHB_L.sig$log2FoldChange) >= lfc_cutoff,]

res_DP_L <- results(dds, contrast = c("condition", "L", "DP"))
res_DP_L <- na.omit(res_DP_L)
res_DP_L.sig <- res_DP_L[res_DP_L$padj <= pvalue_cutoff,]
res_DP_L.sig <- res_DP_L.sig[abs(res_DP_L.sig$log2FoldChange) >= lfc_cutoff,]


# Extrat norm counts and plot PCA
norm.counts <- counts(dds, normalized=TRUE)
write.table(norm.counts, "deseq_normalized_counts.txt", sep="\t", quote=F, col.names = NA)

# Claculate regularized logarithms

#vst.vals <-vst(dds, blind = FALSE)
vst.vals <-vst(dds, blind = FALSE)
vst.mat <- assay(vst.vals)

vst.norm <- removeBatchEffect(vst.mat,df$batch)
#pca plot from normalized
library(factoextra)
res.pca <- prcomp(vst.norm, scale = TRUE)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

library("vsn")

ddsBlind <- estimateDispersions( dds )
vsd <- varianceStabilizingTransformation( ddsBlind )

#par(mfrow=c(1,2))
notAllZero = (rowSums(counts(dds))>0)
pdf("log2+1meansd.pdf")
meanSdPlot(log2(counts(dds)[notAllZero, ] + 1))
dev.off()
pdf("vsdmeansd.pdf")
meanSdPlot(vst.vals[notAllZero, ])
dev.off()


pdf("rawQT.pdf", width=1200, height = 900)
meanSdPlot(txi$counts[notAllZero, ])
dev.off()

#rlog
rld <- rlog(dds, blind = FALSE, fitType='mean')

# In sparseTest(counts(object, normalized = TRUE), 0.9, 100, 0.1) :
#   the rlog assumes that data is close to a negative binomial distribution, an assumption
# which is sometimes not compatible with datasets where many genes have many zero counts
# despite a few very large counts.
# In this data, for 16.2% of genes with a sum of normalized counts above 100, it was the case 
# that a single sample's normalized count made up more than 90% of the sum over all samples.
# the threshold for this warning is 10% of genes. See plotSparsity(dds) for a visualization of this.
# We recommend instead using the varianceStabilizingTransformation or shifted log (see vignette).

######## Blind ###########
# blind means the variable in the design is taken or not into account of the
# analysis. TRUE would be completely unsupervised transfromation
##########################


#analysis plot
library("dplyr")
library("ggplot2")
library("magrittr")

dds <- estimateSizeFactors(dds)

df1 <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vst.vals)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df1)[1:2] <- c("x", "y")  
df1 <- data.frame(df1)

png("variance_plotQT.png", width=1200, height = 900)
ggplot(df1, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

dev.off()

##### SAMPLE DISTANCE #####
sampleDists <- dist(t(assay(vst.vals)))

library("pheatmap")
library("RColorBrewer")




#####
pca <- plotPCA(vst.vals, intgroup = "condition")
png("deseq_vst_pcaQT.png", width=1200, height = 900)
pca
dev.off()


### PLOT for top expressed gene #########

#compare DP and EHB
topGene <- rownames(res_DP_EHB)[which.min(res$res_DP_EHB)]
plotCounts(dds, gene = topGene, intgroup= "condition")

####### Gene cluster ############
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vst.vals)), decreasing = TRUE), 50)

mat  <- assay(vst.vals)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vst.vals)[, "condition"])
rownames(anno) <- colnames(vst.vals)
png("geneFilterHeatmapQT50.png", width=1200, height = 900)
pheatmap(mat, annotation_col = anno)
dev.off()



##### MORE HEATMAP #######
cdsFull <- estimateSizeFactors( dds )
cdsFull <- estimateDispersions( cdsFull )
cdsFullBlind = estimateDispersions( cdsFull, method = "blind" )

vsdFull = varianceStabilizingTransformation( cdsFullBlind )
library("RColorBrewer")
library("gplots")
select = order(rowMeans(counts(cdsFull)), decreasing=TRUE)[1:50] #more select yield better clade
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
png("vsdHM.png",width=1200, height = 900)
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
dev.off()
#untransformed data
png("ddsHM.png",width=1200, height = 900)
heatmap.2(counts(cdsFull)[select,], col = hmcol, trace="none", margin=c(10,6))
dev.off()




#Account for batch effect
#library(sva)

