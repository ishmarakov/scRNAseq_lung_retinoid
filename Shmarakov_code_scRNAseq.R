knitr::opts_chunk$set(echo = TRUE)

library(Seurat)
library(Matrix)
library(fossil) 
library(dplyr)
library(plyr)
library(liger)
library(ggplot2)
library(cowplot)

# Setting folder location for saving output files. This is also the same location as input data.
mydir <- "C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/"
setwd("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/") 

# Objects to save.
Rda.sparse.path <- paste0(mydir, "UV-cells_subsample.Rda")
Rda.path <- paste0(mydir, "UV-cells_nobatchcorrect.Rda")
Rda.Seurat3.path <- paste0(mydir, "UV-cells_Seurat3.Rda")
Rda.liger.path <- paste0(mydir, "UV-cells_liger.Rda")

# Read in all four input expression matrices
celseq.data <- Read10X(data.dir = ("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/WI001/WI001_cellranger_count_outs/filtered_feature_bc_matrix"))
celseq2.data <- Read10X(data.dir = "C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/WI02/WI02_Analysis/WI02_Analysis_Extracted/WI02_cellranger_count_outs/filtered_feature_bc_matrix");

dim(celseq.data)
dim(celseq2.data)

# Creating and setting up Seurat objects for each dataset with the following 6 steps.
# 1. CreateSeuratObject
# 2. subset
# 3. NormalizeData
# 4. FindVariableFeatures
# 5. ScaleData 
# 6. Update @meta.data slot in Seurat object with tech column (celseq, celseq2)
# Look at the distributions of number of genes per cell before and after FilterCells.

# creating seurat object celseq (from sample WI01) and excluding samples with less than 5 cells and less than 200 detected genes
celseq <- CreateSeuratObject(counts = celseq.data, min.cells = 5, min.features = 200)

# calculating the mitochondrial transcript percentage for each cell:
mito.genes <- grep(pattern = "^MT-", x = rownames(x = celseq), value = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = celseq, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = celseq, slot = 'counts'));
celseq[['percent.mito']] <- percent.mito;

# calculating the ribosomal transcript percentage for each cell:
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = celseq), value = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = celseq, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = celseq, slot = 'counts'));
celseq[['percent.ribo']] <- percent.ribo;

# ploting as violin plots
pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/VlnPlot.WI01.pdf"), width = 13, height = 6);
vln <- VlnPlot(object = celseq, features = c("percent.mito", "percent.ribo"), ncol = 2);
print(vln);
dev.off();

pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/VlnPlot.nCount.25Kmax.WI01.pdf"), width = 10, height = 10)
vln <- VlnPlot(object = celseq, features = "nCount_RNA", y.max=25000)
print(vln)
dev.off();

pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/VlnPlot.nFeature.WI01.pdf"), width = 10, height = 10)
vln <- VlnPlot(object = celseq, features = "nFeature_RNA")
print(vln)
dev.off()

# using Seurat's FeatureScatter function to create scatterplots of the relationships among QC variables
pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Scatter1.WI01.pdf"), width = 8, height = 6);
scatter <- FeatureScatter(object = celseq, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size=0.1)
print(scatter);
dev.off();

pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Scatter2.WI01.pdf"), width = 8, height = 6);
scatter <- FeatureScatter(object = celseq, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size=0.1)
print(scatter);
dev.off();

pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Scatter3.WI01.pdf"), width = 8, height = 6);
scatter <- FeatureScatter(object = celseq, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1)
print(scatter);
dev.off();

# filtering the cells to remove debris, dead cells, and probable doublets
min <- min(celseq@meta.data$nFeature_RNA);
m <- median(celseq@meta.data$nFeature_RNA)
max <- max(celseq@meta.data$nFeature_RNA)    
s <- sd(celseq@meta.data$nFeature_RNA)
min1 <- min(celseq@meta.data$nCount_RNA)
max1 <- max(celseq@meta.data$nCount_RNA)
m1 <- mean(celseq@meta.data$nCount_RNA)
s1 <- sd(celseq@meta.data$nCount_RNA)
Count93 <- quantile(celseq@meta.data$nCount_RNA, 0.93) # calculate value in the 93rd percentile
print(paste("Feature stats:",min,m,max,s));
print(paste("UMI stats:",min1,m1,max1,s1,Count93));

# filtering cells to exclude cells with less than 200 genes and more than 10 % mitochondrial genes 
celseq <- subset(x = celseq, subset = nFeature_RNA > 200  & nCount_RNA < Count93 & percent.mito < 0.1)
pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/VlnPlot.nFeature.after.subset.WI01.pdf"), width = 10, height = 10)
vln <- VlnPlot(object = celseq, features = "nFeature_RNA")
print(vln)
dev.off()

# normalizing and scaling data
celseq <- NormalizeData(celseq, normalization= "LogNormalize", scale.factor = 10000)
celseq <- FindVariableFeatures(celseq, selection.method = "vst", nfeatures = 2000)
celseq <- ScaleData(celseq)
celseq[["tech"]] <- "celseq"

# celseq2 (sample WI02)
# creating seurat object celseq2 (from sample WI02) and excluding samples with less than 5 cells and less than 200 detected genes
celseq2 <- CreateSeuratObject(counts = celseq2.data, min.cells = 5, min.features = 200)

# calculating the mitochondrial transcript percentage for each cell:
mito.genes <- grep(pattern = "^MT-", x = rownames(x = celseq2), value = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = celseq2, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = celseq2, slot = 'counts'));
celseq2[['percent.mito']] <- percent.mito;

# calculating the ribosomal transcript percentage for each cell:
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = celseq2), value = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = celseq2, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = celseq2, slot = 'counts'));
celseq2[['percent.ribo']] <- percent.ribo;

# ploting as violin plots
pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/VlnPlot.WI02.pdf"), width = 13, height = 6);
vln <- VlnPlot(object = celseq2, features = c("percent.mito", "percent.ribo"), ncol = 2);
print(vln);
dev.off();

pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/VlnPlot.nCount.40Kmax.WI02.pdf"), width = 10, height = 10)
vln <- VlnPlot(object = celseq2, features = "nCount_RNA", y.max=40000)
print(vln)
dev.off();

pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/VlnPlot.nFeature.WI02.pdf"), width = 10, height = 10)
vln <- VlnPlot(object = celseq2, features = "nFeature_RNA")
print(vln)
dev.off()

# using Seurat's FeatureScatter function to create scatterplots of the relationships among QC variables
pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Scatter1.WI02.pdf"), width = 8, height = 6);
scatter <- FeatureScatter(object = celseq2, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size=0.1)
print(scatter);
dev.off();

pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Scatter2.WI02.pdf"), width = 8, height = 6);
scatter <- FeatureScatter(object = celseq2, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size=0.1)
print(scatter);
dev.off();

pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Scatter3.WI02.pdf"), width = 8, height = 6);
scatter <- FeatureScatter(object = celseq2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1)
print(scatter);
dev.off();

# filtering the cells to remove debris, dead cells, and probable doublets
min <- min(celseq2@meta.data$nFeature_RNA);
m <- median(celseq2@meta.data$nFeature_RNA)
max <- max(celseq2@meta.data$nFeature_RNA)    
s <- sd(celseq2@meta.data$nFeature_RNA)
min1 <- min(celseq2@meta.data$nCount_RNA)
max1 <- max(celseq2@meta.data$nCount_RNA)
m1 <- mean(celseq2@meta.data$nCount_RNA)
s1 <- sd(celseq2@meta.data$nCount_RNA)
Count93 <- quantile(celseq2@meta.data$nCount_RNA, 0.93) # calculate value in the 93rd percentile
print(paste("Feature stats:",min,m,max,s));
print(paste("UMI stats:",min1,m1,max1,s1,Count93));

# filtering cells to exclude cells with less than 200 genes and more than 10 % mitochondrial genes 
celseq2 <- subset(x = celseq2, subset = nFeature_RNA > 200  & nCount_RNA < Count93 & percent.mito < 0.1)
pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/VlnPlot.nFeature.after.subset.WI02.pdf"), width = 10, height = 10)
vln <- VlnPlot(object = celseq2, features = "nFeature_RNA")
print(vln)
dev.off()

# normalizing and scaling data
celseq2 <- NormalizeData(celseq2, normalization= "LogNormalize", scale.factor = 10000)
celseq2 <- FindVariableFeatures(celseq2, selection.method = "vst", nfeatures = 2000)
celseq2 <- ScaleData(celseq2)
celseq2[["tech"]] <- "celseq2"

# Save the sub-sampled Seurat objects
save(celseq, celseq2, file = Rda.sparse.path)

load(Rda.sparse.path)

# Merging Seurat objects. Original sample identities are stored in gcdata[["tech"]].
# Cell names will now have the format tech_cellID (smartseq2_cell1...)
add.cell.ids <- c("celseq", "celseq2")
gcdata <- merge(x = celseq, y = list(celseq2), add.cell.ids = add.cell.ids, merge.data = FALSE)
Idents(gcdata) <- "tech"  # use identity based on sample identity

# Looking at how the number of genes per cell varies across the different technologies.
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/VlnPlot.nFeature_RNA.merged1.jpg"), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(gcdata, "nFeature_RNA", group.by = "tech")
print(vln)
dev.off();

#Normalizing and scaling the merged data. 
gcdata <- NormalizeData(gcdata, normalization.method = "LogNormalize", scale.factor = 10000)
var.genes <- SelectIntegrationFeatures(SplitObject(gcdata, split.by = "tech"), nfeatures = 2000, verbose = TRUE, fvf.nfeatures = 2000, selection.method = "vst")
VariableFeatures(gcdata) <- var.genes
gcdata <- ScaleData(gcdata, features = VariableFeatures(gcdata))


# Doing PCA on data including only the variable genes.
gcdata <- RunPCA(gcdata, features = VariableFeatures(gcdata), npcs = 100, ndims.print = 1:15, nfeatures.print = 15)

# Coloring the PC biplot by the scRNA-seq technology.
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/PC_biplot.jpg"), width = 10, height = 8, units="in", res=300);
p1 <- DimPlot(gcdata, reduction = "pca", dims = c(1, 2), group.by = "tech")
print(plot_grid(p1));
dev.off();

##Generating an elbow plot of principal component standard deviations
elbow <- ElbowPlot(object = gcdata)
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/PCA.elbow.jpg"), width = 6, height = 8, units="in", res=300);
print(elbow);
dev.off();

##using a bootstrapping technique called Jackstraw analysis to estimate a p-value for each component, print out a plot, and save the p-values to a file
gcdata <- JackStraw(object = gcdata, num.replicate = 100, dims=30); # takes around 4 minutes
gcdata <- ScoreJackStraw(object = gcdata, dims = 1:30)
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/PCA.jackstraw.jpg"), width = 10, height = 6, units="in", res=300);
js <- JackStrawPlot(object = gcdata, dims = 1:30)
print(js);
dev.off();
pc.pval <- gcdata@reductions$pca@jackstraw@overall.p.values; # get p-value for each PC
write.table(pc.pval, file=sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/PCA.jackstraw.scores.xls", date), quote=FALSE, sep='\t', col.names=TRUE);

# Clustering the cells using the first twenty principal components.
gcdata <- FindNeighbors(gcdata, reduction = "pca", dims = 1:15, k.param = 20)

gcdata <- FindClusters(gcdata, resolution = 0.2, algorithm = 1, random.seed = 100)

# Creating a UMAP visualization. 
gcdata <- RunUMAP(gcdata, dims = 1:15, reduction = "pca", n.neighbors = 15, min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)  

# Visualizing the Louvain clustering
# the clustering is stored in @meta.data in column seurat_clusters and the technology is stored in the column tech. 
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Louvain_clustering.jpg"), width = 10, height = 8, units="in", res=300);
p2 <- DimPlot(gcdata, reduction = "umap", group.by = "seurat_clusters")
print(plot_grid(p2));
dev.off();

# Visualizing the batches on the UMAP
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/batches_UMAP.jpg"), width = 10, height = 8, units="in", res=300);
p3 <- DimPlot(gcdata, reduction = "umap", group.by = "tech")
print(plot_grid(p3));
dev.off();

# Adjusted rand index test for overlap between technology and cluster labelings. 
# This goes between 0 (completely dissimilar clustering) to 1 (identical clustering). 
# The adjustment corrects for chance grouping between cluster elements.
ari <- dplyr::select(gcdata[[]], tech, seurat_clusters)
ari$tech <- plyr::mapvalues(ari$tech, from = c("celseq", "celseq2"), to = c(0, 1))
adj.rand.index(as.numeric(ari$tech), as.numeric(ari$seurat_clusters))

# Saving current progress.
save(gcdata, file = Rda.path)
# To load the data, run the following command.
# load(Rda.path)

# Batch correction: canonical correlation analysis (CCA) + mutual nearest neighbors (MNN) using Seurat v3
# The first piece of code will identify variable genes that are highly variable in at least 2/4 datasets. We will use these variable genes in our batch correction.
ob.list <- list(celseq, celseq2)

# Identifying anchors on the 2 datasets, commonly shared variable genes across samples, and integrate samples.
gcdata.anchors <- FindIntegrationAnchors(object.list = ob.list, anchor.features = 2000, dims = 1:15)
gcdata <- IntegrateData(anchorset = gcdata.anchors, dims = 1:15)

DefaultAssay(gcdata) <- "integrated"

# Running the standard workflow for visualization and clustering.
# The integrated data object only stores the commonly shared variable genes.
gcdata <- ScaleData(gcdata, do.center = T, do.scale = F)

gcdata <- RunPCA(gcdata, npcs = 40, ndims.print = 1:5, nfeatures.print = 5)

jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Integrated_batches_UMAP.jpg"), width = 10, height = 8, units="in", res=300);
p4 <- DimPlot(gcdata, dims = c(1, 2), reduction = "pca", split.by = "tech")
print(plot_grid(p4));
dev.off();

# Clustering. Choosing the dimensional reduction type to use and the number of aligned canonical correlation vectors to use.
gcdata <- FindNeighbors(gcdata, reduction = "pca", dims = 1:15, k.param = 20)

#Clustering with different resolutions
gcdata <- FindClusters(gcdata, resolution = 0.2, algorithm = 1, random.seed = 100)
gcdata <- FindClusters(gcdata, resolution = 0.3, algorithm = 1, random.seed = 100)
gcdata <- FindClusters(gcdata, resolution = 0.4, algorithm = 1, random.seed = 100)
gcdata <- FindClusters(gcdata, resolution = 0.5, algorithm = 1, random.seed = 100)
gcdata <- FindClusters(gcdata, resolution = 0.6, algorithm = 1, random.seed = 100)
gcdata <- FindClusters(gcdata, resolution = 0.7, algorithm = 1, random.seed = 100)
gcdata <- FindClusters(gcdata, resolution = 0.8, algorithm = 1, random.seed = 100)
gcdata <- FindClusters(gcdata, resolution = 0.9, algorithm = 1, random.seed = 100)


# Building clustering tree to know cluster resolution value
library(clustree)
head(gcdata[[]])
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Clustree.jpg"), width = 10, height = 8, units="in", res=300);
p41 <- clustree(gcdata)
print(plot_grid(p41));
dev.off();

gcdata <- FindClusters(gcdata, resolution = 0.2, algorithm = 1, random.seed = 100)

# UMAP. Choosing the dimensional reduction type to use and the number of aligned canonical correlation vectors to use.
gcdata <- RunUMAP(gcdata, dims = 1:15, reduction = "pca", n.neighbors = 15, min.dist = 0.5, spread = 1, metric = "euclidean", seed.use = 1)  

## After data integration, using the original expression data in all visualization and DE tests.
DefaultAssay(gcdata) <- "RNA"  

# Visualizing the Louvain clustering and the batches on the UMAP. 
# The clustering is stored in @meta.data in column seurat_clusters 
# and the technology is stored in the column tech. Remember you can also use DimPlot.
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Intergated_UMAP_clusters.jpg"), width = 16, height = 8, units="in", res=300);
p5 <- DimPlot(gcdata, reduction = "umap", group.by = "seurat_clusters")
p6 <- DimPlot(gcdata, reduction = "umap", group.by = "tech")
print(plot_grid(p5, p6));
dev.off();

# Looking to see how the adjusted rand index changed compared to using no batch correction.
ari <- dplyr::select(gcdata[[]], tech, seurat_clusters)
ari$tech <- plyr::mapvalues(ari$tech, from = c("celseq", "celseq2"), to = c(0, 1))
adj.rand.index(as.numeric(ari$tech), as.numeric(ari$seurat_clusters))

# Identifying conserved marker genes across the batches. Differential gene expression is done across each batch, and the p-values are combined.
markers <- FindConservedMarkers(gcdata, ident.1 = 0, grouping.var = "tech", assay = "RNA", print.bar = T)
head(markers)
write.table(markers, "C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Conserved.Markers.xls", na = "NA", append = FALSE, col.names = TRUE, row.names = TRUE, sep = "/", quote = TRUE)

# Visualizing the expression of the first 5 marker genes on UMAP across the different batches using DoHeatmap.
gcdata <- ScaleData(gcdata, features = rownames(gcdata), do.center = T, do.scale = F)
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Intergated_Heatmap.jpg"), width = 10, height = 8, units="in", res=300);
p7 <- DoHeatmap(gcdata, features = rownames(markers)[1:5], group.by = "tech", disp.max = 3)
print(plot_grid(p7));
dev.off();

# Markers for clusters
genes <- c("Col1a1", "Epcam", "Pecam1", "Ptprc")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Intergated_Feature_plot.jpg"), width = 10, height = 8, units="in", res=300);
p8 <- FeaturePlot(gcdata, genes, ncol = 2)
print(plot_grid(p8));
dev.off();

genes <- c("Lrat", "Rbp1", "Rbp4", "Cd36")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Intergated_Ret_genes_plot.jpg"), width = 10, height = 8, units="in", res=300);
p9 <- FeaturePlot(gcdata, genes, ncol = 2)
print(plot_grid(p9));
dev.off();

# performing DEG analysis on all clusters simultaneously using the default differential expression test (Wilcoxon), then save the results to a file.
DEGs <- FindAllMarkers(object=gcdata, logfc.threshold=1, min.diff.pct=.2);
write.table(DEGs, file=sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/DEGs.Wilcox.xls"), quote=FALSE, sep="/", row.names=FALSE);

# choosing the top 10 DEGs in each cluster, and print them to a heatmap using DoHeatmap and a red/white/blue color scheme
top10 <- DEGs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
pdf(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/heatmap.pdf"), height=20, width=15);
DoHeatmap(gcdata, features=top10$gene, slot="scale.data", disp.min=-2, disp.max=2, group.by="ident", group.bar=TRUE) + scale_fill_gradientn(colors = c("blue", "white", "red")) + theme(axis.text.y = element_text(size = 10));
dev.off();

# Building violin plots
library(Seurat)
library(patchwork)
library(ggplot2)

## removing the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

features<- c("Col1a1", "Pdgfra", "Ptgis", "Sftpa1", "Sftpb", "Epcam", "Cdh5", "Pecam1", "Plvap", "Ptprc", "Il1b", "Lsp1")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Intergated_Violin_plot_markers.jpg"), width = 10, height = 16, units="in", res=300);
p10 <- StackedVlnPlot(obj = gcdata, features = features)
print(plot_grid(p10));
dev.off();

features<- c("Lrat", "Rbp1", "Rbp4", "Lpl", "Gpihbp1", "Cd36", "Scarb1", "Vldlr", "Ldlr", "Lrp1", "Apoe", "Dgat1", "Dgat2")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Intergated_Violin_plot_Retinoid_Lipid_markers.jpg"), width = 10, height = 16, units="in", res=300);
p31 <- StackedVlnPlot(obj = gcdata, features = features)
print(plot_grid(p31));
dev.off();

## extracting the information on cell numbers in clusters
library(data.table)
library(magrittr)

# extract meta data
md <- gcdata@meta.data %>% as.data.table
# the resulting md object has one "row" per cell

# count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md[, .N, by = c("seurat_clusters")]

saveRDS(gcdata, "C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/gcdata.rds")

## Subsetting fibroblast clusters only in CF.only
CFonly<- subset(x=gcdata, idents = c('0','1','2'))

DefaultAssay(CFonly) <- "RNA"

saveRDS(CFonly, "C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/CF.only.rds")

CFonly <- readRDS(file = "C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/CFonly.rds")

## Fibroblast genes from Hurskainen
features <- c("Mfap4", "Macf1", "Npnt", "Limch1", "Gyg", "Tgfbi", "Aspn", "Htra1", "Agt", "Pdlim3", "Dcn", "Col14a1", "Mfap5", "Col1a1", "Col1a2", "Crip1", "Acta2", "Myh11", "Tagln", "Igfbp5", "Pdzd2", "Higd1b", "Gucy1a1", "Postn", "Gucy1b1", "Actc1", "Actg2", "Tnnt2", "Des")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Hurskainen.jpg"), width = 10, height = 4, units="in", res=300);
p12 <- DotPlot(obj = CFonly, features = features) + theme(axis.text.x = element_text(angle = 90, hjust=1))
print(plot_grid(p12));
dev.off();

## Fibroblast genes from Liu
features <- c("Tcf21", "Plin2", "Fgf10", "G0s2", "Gyg", "Macf1", "Wnt2", "Col13a1", "Acta2", "Myh11", "Tagln", "Pdgfra", "Tgfbi", "Hhip", "Enpp2", "Wnt5a", "Pdgfrb", "Higd1b", "Cox4i2", "Notch3", "Ebf1", "Gucy1a3", "Pdzd2", "Postn", "Agtr2", "Prss35", "Igfbp7", "Fbln5", "Ptn", "Heyl", "Fstl1", "Tm4sf1", "Hmmr", "Mki67", "Pcna", "Top2a", "Spc25", "Cdca3", "Ccnb2", "Hist1h2ap", "Upk3b", "Wt1", "Msln", "Upk1b", "Aldh1a2", "Cpe", "Gm12840", "Aqp1")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Liu.jpg"), width = 10, height = 4, units="in", res=300);
p13 <- DotPlot(obj = CFonly, features = features) + theme(axis.text.x = element_text(angle = 90, hjust=1))
print(plot_grid(p13));
dev.off();

## Fibroblast genes combined from Xie
features <- c("Itga8", "Cxcl14", "Npnt", "Hsd11b1", "Tcf21", "Mfap4", "Spon1", "Limch1", "Cdh11", "Col13a1", "Pi16", "Mmp3", "Clec3b", "Cygb", "Dcn", "Rbp4", "Gsn", "Col14a1", "Dpep1", "Fbln1", "Hhip", "Aspn", "Mustn1", "Enpp2", "Igfbp3", "Grem2", "Tagln", "Lum", "Bmp5", "Acta2", "Ear2", "Mrc1", "Ccl6", "Plet1", "Ctss", "Mpeg1", "Atp6v0d2", "Chil3", "Abcg1", "Krt79", "Cybb", "Lpl", "Plin2", "Higd1b", "Cox4i2", "Notch3", "Fam162b", "Postn", "Col8a1", "Lmcd1", "Tmem178", "Hbegf", "Lipg", "Pdgfrb", "Myh11", "Vsnl1")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Xie.jpg"), width = 10, height = 4, units="in", res=300);
p14 <- DotPlot(obj = CFonly, features = features) + theme(axis.text.x = element_text(angle = 90, hjust=1))
print(plot_grid(p14));
dev.off();

## Fibroblast markers combined from Travaglini
features <- c("Col1a1", "Col1a2", "Pdgfra", "Bsg", "Eln", "Bgn", "Tagln", "Acta2", "Itga8", "Vegfd", "Gpm6b", "Tbx2", "Spint2", "Fgfr4", "Gpc3", "Dkk3", "Aoc3", "Nkd1", "Rspo1", "Ptgis", "Serpinf1", "Pi16", "Mfap5", "Igfbp4", "C3", "Pdgfrl", "Sfrp2", "Pdlim4", "Thy1", "Scara5", "Il32", "Apoe", "Fst", "Plin2", "Aspn", "Wif1", "Fgf18", "Scx", "Lgr6", "Myh11", "Cnn1", "Actg2", "Cox4i2", "Rergl", "Kcna5")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Travaglini.jpg"), width = 10, height = 4, units="in", res=300);
p15 <- DotPlot(obj = CFonly, features = features) + theme(axis.text.x = element_text(angle = 90, hjust=1))
print(plot_grid(p15));
dev.off();

## Fibroblast markers from Tsukui
features <- c("Npnt", "Ces1d", "Slc7a10", "Pi16", "Ccl11", "Il33", "Adh7", "Dcn", "Hhip", "Aspn", "Fgf18", "Pdgfrb", "Higd1b", "Cox4i2")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Tsukui.jpg"), width = 10, height = 4, units="in", res=300);
p16 <- DotPlot(obj = CFonly, features = features) + theme(axis.text.x = element_text(angle = 90, hjust=1))
print(plot_grid(p16));
dev.off();

## Fibroblast genes Grouped
features <- c("Slc7a10", "Tagln", "Fgfr4", "Ces1d", "Limch1", "Col13a1", "Pdgfra", "Npnt", "Bsg", "Itga8", "Vegfd", "Tbx2", "Ptgis", "Mfap4", "Spon1", "Cdh11", "Apoe", "Macf1", "Tcf21", "Rbp1", "Lrp1", "G0s2", "Gyg", "Wnt2", "Pdzd2", "Plin2", "Lrat", "Lpl", "Col1a1", "Col1a2", "Igfbp7", "Bgn", "Gpc3", "Eln", "C3", "Mfap5", "Fstl1", "Dcn","Dkk3", "Cygb", "Mmp3", "Rbp4", "Heyl", "Col14a1", "Tm4sf1", "Adh7", "Pi16", "Serpinf1", "Crip1", "Htra1", "Pdgfrb", "Hhip", "Acta2", "Tgfbi", "Aspn")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Fibroblast_Genes_Grouped.jpg"), width = 10, height = 4, units="in", res=300);
p17 <- DotPlot(obj = CFonly, features = features) + theme(axis.text.x = element_text(angle = 90, hjust=1))
print(plot_grid(p17));
dev.off();

## Fibroblast cluster markers
features <- c("Fos", "Ccn1", "Atf3", "Ier2", "Egr1", "Zfp36", "Ccn2", "Rhob", "1200007C13Rik", "Nfkbiz", "Sat1", "Nr4a1", "Ppp1r15a", "Ier3", "Klf4", "Klf2", "Btg2", "Fosb", "Socs3", "Jun", "Mt1", "Junb", "Ces1d", "Slc7a10", "Tagln", "Gm14964", "Cdkn2c", "Hopx", "Igfbp3", "Cdkn1a", "Maff", "Gadd45g", "Mt2", "Dcn", "Rbp4", "C4b", "Mmp3", "Cygb", "Igfbp6", "Nbl1", "Col3a1", "C3", "Fxyd6",  "Col14a1",  "Cpxm1",  "Mfap5",  "Mmp2", "Fbln1",  "Pmepa1", "Gas6", "Man2a1", "Dpt",  "Col1a1", "Gas1", "Sfrp1",  "Eln", "Adamts2", "Ucp2", "Igfbp4", "Cd9", "Ly6e")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Fibroblast_Cluster_Markers.jpg"), width = 10, height = 4, units="in", res=300);
p18 <- DotPlot(obj = CFonly, features = features) + theme(axis.text.x = element_text(angle = 90, hjust=1))
print(plot_grid(p18));
dev.off();

## Fibroblast cluster refined
features <- c("Ccn1", "Ccn2", "Rhob", "Ier2", "Jun", "Junb", "Jund", "Fos", "Fosb", "Egr1", "Atf3", "Ppp1r15a", "Zfp36", "Socs3", "Btg2", "Cebpd", "Dusp1", "Sat1", "Klf6", "Nr4a1", "Tagln", "Pdzd2", "Hopx", "Igfbp3", "Cdkn2c", "Slc7a10", "Ces1d", "Gm14964", "Pmepa1", "Cav1", "Sfrp1", "Nbl1", "Gas1", "Gas6", "Gadd45g", "Cdkn1a", "Klf2", "Klf4", "Igfbp4", "Igfbp6", "Col1a1", "Col3a1", "Col14a1", "Eln", "Fbln1", "Dpt", "Dcn", "Mfap5", "Mmp2", "Mmp3", "Timp1", "Cpxm1", "Man2a1", "Cygb", "Mt1", "Mt2", "Fxyd6", "Fxyd5", "Rbp4", "Pltp", "C3", "C4b")
jpeg(sprintf("C:/Users/igors/Documents/Igor_docs_2019/Lung/SingleCellTranscriptomics/Seurat_Analysis/Batch_corrected/Fibroblast_Cluster_Markers_Refined.jpg"), width = 10, height = 4, units="in", res=300);
p18 <- DotPlot(obj = CFonly, features = features) + theme(axis.text.x = element_text(angle = 90, hjust=1))
print(plot_grid(p18));
dev.off();



