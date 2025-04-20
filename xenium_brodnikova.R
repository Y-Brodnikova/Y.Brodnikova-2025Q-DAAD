library(Voyager)
library(SFEData) 
library(SingleCellExperiment)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(ggplot2)
library(ggforce)
library(stringr)
library(scater) 
library(scuttle)
library(scales)
library(BiocParallel)
library(BiocSingular)
library(bluster)
library(scran)
library(patchwork)
library(RBioFormats)
library(fs)
library(sf)
library(arrow)
library(dplyr)
library(tidyr)
library(BiocNeighbors)

library(tiff)
library(rjson)
library(Matrix)
library(SFEData)
library(RBioFormats)
library(ExperimentHub)



theme_set(theme_bw())
options(timeout = Inf)

base <- "."
fn <- file.path(base, "Xenium_V1_FFPE_TgCRND8_5_7_months_outs")

#Read Xenium outputs 
(sfe <- readXenium(fn, row.names = "symbol", add_molecules = TRUE))

#There are 58682 cells in this dataset

#Plot the tissue
plotGeometry(sfe, colGeometryName = "cellSeg")+
  ggtitle("Plot the cell") +
  theme(plot.title = element_text(hjust = 0.5)) 

plotGeometry(sfe, colGeometryName = "nucSeg")+
  ggtitle("Plot the nuclei") +
  theme(plot.title = element_text(hjust = 0.5)) 

plotCellBin2D(sfe, binwidth = 25, hex = TRUE)+
  ggtitle("Plot cell density in space") +
  theme(plot.title = element_text(hjust = 0.5))

plotTxBin2D(sfe, binwidth = 20, flip = TRUE, hex = TRUE)+
  ggtitle("Plot total transcript density") +
  theme(plot.title = element_text(hjust = 0.5))

#So, there are different areas: high cell density and lower cell density. They also have a very different spatial structure. 
#A region with a higher cell density also has a higher transcript density, but in a sparse region there are regions with 
#a higher transcript density but not a higher cell density.

#Plot the images
getImg(sfe)
plotImage(sfe, image_id = "morphology_focus", show_axes = TRUE,
          dark = TRUE) +
  ggtitle("Plot the images") +
  theme(plot.title = element_text(hjust = 0.5))


plotImage(sfe, image_id = "morphology_focus", channel = 1L, show_axes = TRUE, dark = TRUE,
          normalize_channels = TRUE, palette = viridis_pal(option = "H")(255)) +
  ggtitle("Plot the images") +
  theme(plot.title = element_text(hjust = 0.5))

bbox1 <- c(xmin = 3500, xmax = 4500, ymin = -1500, ymax = -250)
bbox2 <- c(xmin = 3500, xmax = 4500, ymin = -2450, ymax = -2050)

bboxes_sf <- c(st_as_sfc(st_bbox(bbox1)), st_as_sfc(st_bbox(bbox2)))
plotImage(sfe, image_id = "morphology_focus", channel = 1L, show_axes = TRUE, dark = TRUE, 
          normalize_channels = TRUE, palette = viridis_pal(option = "H")(255)) +
  geom_sf(data = bboxes_sf, fill = NA, color = "white", linewidth = 1) +
  ggtitle("Plot the images") +
  theme(plot.title = element_text(hjust = 0.5))

plotImage(sfe, image_id = "morphology_focus", channel = 1L,
          bbox = bbox1, normalize_channels = TRUE, palette = viridis_pal(option = "H")(255))

plotImage(sfe, image_id = "morphology_focus", channel = 1L,
          bbox = bbox2, normalize_channels = TRUE, palette = viridis_pal(option = "H")(255))

plotGeometry(sfe, colGeometryName = "cellSeg", fill = FALSE, dark = TRUE,
             image_id = "morphology_focus", channel = 1L, bbox = bbox1, 
             normalize_channels = TRUE, palette = viridis_pal(option = "H")(255))
plotGeometry(sfe, colGeometryName = "nucSeg", fill = FALSE, dark = TRUE,
             image_id = "morphology_focus", channel = 1L, bbox = bbox1, 
             normalize_channels = TRUE, palette = viridis_pal(option = "H")(255))
plotGeometry(sfe, colGeometryName = c("cellSeg", "nucSeg"), fill = FALSE, 
             dark = TRUE, image_id = "morphology_focus", channel = 1L, bbox = bbox1, 
             normalize_channels = TRUE, palette = viridis_pal(option = "H")(255))

#Quality control

#Negative controls

names(rowData(sfe))
unique(rowData(sfe)$Type)

is_blank <- rowData(sfe)$Type == "Unassigned Codeword"
sum(is_blank)
is_neg <- rowData(sfe)$Type == "Negative Control Probe"
sum(is_neg)
is_neg2 <- rowData(sfe)$Type == "Negative Control Codeword"
sum(is_neg2)
is_any_neg <- is_blank | is_neg | is_neg2
sfe <- addPerCellQCMetrics(sfe, subsets = list(unassigned = is_blank,
                                               negProbe = is_neg,
                                               negCodeword = is_neg2,
                                               any_neg = is_any_neg))
names(colData(sfe))
cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_percent$")]
plotColDataHistogram(sfe, cols_use, bins = 100) +
  ggtitle("Proportion of transcript counts coming from any negative control") +
  theme(plot.title = element_text(hjust = 0.5))

#The histogram is dominated by the bin at zero and there are some extreme outliers too few to be seen 
#but evident from the scale of the x axis. 

plotColDataHistogram(sfe, cols_use, bins = 100) + 
  scale_x_log10() + 
  annotation_logticks(sides = "b") +
  ggtitle("Histogram only for cells with at least 1 count from a negative control") +
  theme(plot.title = element_text(hjust = 0.5))

#For the most part cells have less than 10% of spots assigned to the negative controls.

cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_sum$")]
plotColDataHistogram(sfe, cols_use, bins = 20, ncol = 2) +
  scale_x_continuous(breaks = scales::breaks_extended(Q = c(1,2,5))) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  ggtitle("Distribution of the number of negative control counts per cell") +
  theme(plot.title = element_text(hjust = 0.5))

#The counts are low, mostly zero, but there are some cells with up to 3 counts of all types aggregated.
#Here we have 5,070 out of 58,682 cells (about 8.64%) that have at least one count from a negative control. 
#However, only 0.047% of all transcripts come from the negative control, indicating a low false positive rate overall.
  
neg_features <- rownames(sfe)[is_any_neg]
tx <- arrow::open_dataset("E:/Alzheimers/TgCRND8/Xenium_V1_FFPE_TgCRND8_5_7_months_outs/transcripts.parquet")
tx_filtered <- tx |>
  filter(feature_name %in% neg_features) |>
  select(x_location, y_location) |>
  collect()

names(colData(sfe))

#Where are cells with higher proportion of negative control spots located in space?

data("ditto_colors")
sfe$more_neg <- sfe$subsets_any_neg_percent > 5
plotSpatialFeature(sfe, "subsets_any_neg_sum", colGeometryName = "cellSeg",
                   show_axes = TRUE) +
  geom_sf(data = SpatialFeatureExperiment::centroids(sfe)[sfe$more_neg,],
          size = 1.5, color = ditto_colors[1]) +
  theme(legend.position = "bottom") + 
  ggtitle("Cell location with a higher proportion of negative control spots in space") +
  theme(plot.title = element_text(hjust = 0.5))

#Cells that have any spots of negative control and that have a higher proportion of negative control 
#appear randomly distributed in space.

#Where are negative control spots located?
tx_filtered |>
  ggplot(aes(x_location, -y_location)) +
  geom_bin2d(binwidth = 25) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(x = NULL, y = NULL) + 
  ggtitle("Location of negative control spots") +
  theme(plot.title = element_text(hjust = 0.5))

#Generally there are more negative control spots in the region with higher cell and transcript density, 
#but especially regions with high cell density, which is not surprising. There don’t visually seem to be regions 
#with more negative control spots not accounted for by cell and transcript density.

rm(tx)

#Cells

n_panel <- sum(!is_any_neg)
sfe$nCounts <- colSums(counts(sfe)[!is_any_neg,])
sfe$nGenes <- colSums(counts(sfe)[!is_any_neg,] > 0)
colData(sfe)$nCounts_normed <- sfe$nCounts/n_panel
colData(sfe)$nGenes_normed <- sfe$nGenes/n_panel

plotColDataHistogram(sfe, c("nCounts", "nGenes")) + 
  ggtitle("Total number of transcript molecules detected from genes\nand the number of detected genes per total number of genes probed") +
  theme(plot.title = element_text(hjust = 0.5))

plotColDataHistogram(sfe, c("nCounts_normed", "nGenes_normed"))

#The number of detected genes has discrete values, which excludes negative controls

plotSpatialFeature(sfe, "nCounts", colGeometryName = "cellSeg") + 
  ggtitle("Spatial structure of nCounts") +
  theme(plot.title = element_text(hjust = 0.5))

plotSpatialFeature(sfe, "nGenes", colGeometryName = "cellSeg") + 
  ggtitle("Spatial structure of nGenes") +
  theme(plot.title = element_text(hjust = 0.5))

#These plots clearly show that nCounts and nGenes are spatially structured and may be biologically relevant.
#However, the presence of reading artifacts is noted. Some cells and genes are seen to be distant from all others

plotColData(sfe, x="nCounts", y="nGenes", bins = 70) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  ggtitle("Relation between nCounts and nGenes") +
  theme(plot.title = element_text(hjust = 0.5))

plotColDataHistogram(sfe, c("cell_area", "nucleus_area"), scales = "free_y") + 
  ggtitle("Distribution of cell area") +
  theme(plot.title = element_text(hjust = 0.5))

#There’s a very long tail. The nuclei are much smaller than the cells

plotSpatialFeature(sfe, "cell_area", colGeometryName = "cellSeg", show_axes = TRUE) + 
  ggtitle("Cell area distributed in space") +
  theme(plot.title = element_text(hjust = 0.5))

#The largest cells are clearly structured, which may be a biological reason

plotSpatialFeature(sfe, "nucleus_area", colGeometryName = "nucSeg", show_axes = TRUE) + 
  ggtitle("Nuclei area distributed in space") +
  theme(plot.title = element_text(hjust = 0.5))

#Zoom into smaller regions to see the nature of the very large cells

bbox3 <- c(xmin = 3500, xmax = 4500, ymin = -1500, ymax = -250)
plotGeometry(sfe, colGeometryName = c("cellSeg", "nucSeg"), fill = FALSE, 
             dark = TRUE, image_id = "morphology_focus", 
             channel = 1L, bbox = bbox3, normalize_channels = TRUE) + 
  ggtitle("Dentate gyrus") +
  theme(plot.title = element_text(hjust = 0.5))

#Large cells in this case are not an artifact, because this is a dentate gyrus characterized by high lithin density 
#and predominant expression of DGsg with little or no expression of DGpo and/or DGmo.

#There may also be some cells that do not have nuclei, that may have nuclei in a different z-plane, or that may be 
#false positives, as well as some false negatives in areas with fluorescence and seemingly nuclei but no cell segmentation.

colData(sfe)$prop_nuc <- sfe$nucleus_area / sfe$cell_area
summary(sfe$prop_nuc)

plotColDataHistogram(sfe, "prop_nuc") + 
  ggtitle("Distribution of the fraction of the cell in the z-plane occupied by the nucleus") +
  theme(plot.title = element_text(hjust = 0.5))

#The NA’s are cells that don’t have nuclei detected.There are no cells with fully occupied nuclei.
#In general, the nucleus occupies about 10% of the cell area 

plotSpatialFeature(sfe, "prop_nuc", colGeometryName = "cellSeg") + 
  ggtitle("Nuclei proportion in space") +
  theme(plot.title = element_text(hjust = 0.5))

#Cells in some histological regions have larger proportions occupied by the nuclei

#Does cell area relate to nCounts, nucleus area, proportion of area taken up by nucleus, and so on?

cols_use <- c("nCounts", "nGenes", "cell_area", "nucleus_area", "prop_nuc")
df <- colData(sfe)[,cols_use] |> as.data.frame()

ggplot(df) +
  geom_bin2d(aes(x = .panel_x, y = .panel_y), bins = 30) +
  geom_autodensity(fill = "gray90", color = "gray70", linewidth = 0.2) +
  facet_matrix(vars(tidyselect::everything()), layer.diag = 2) +
  scale_fill_distiller(palette = "Blues", direction = 1) + 
  ggtitle("Correlation of cell area with the number of cells, nucleus area, \nand the proportion of the area occupied by the nucleus") +
  theme(plot.title = element_text(hjust = 0.5))

#So, in general, nCounts, nGenes, cell area and nucleus area, and fraction of cell area are positively correlated with each other

sfe$has_nuclei <- !is.na(sfe$nucleus_area)

plotSpatialFeature(sfe, "has_nuclei", colGeometryName = "cellSeg",
                   show_axes = TRUE) +
  geom_sf(data = SpatialFeatureExperiment::centroids(sfe)[!sfe$has_nuclei,], 
          color = ditto_colors[1], size = 0.3) + 
  ggtitle("Cells containing a nucleus") +
  theme(plot.title = element_text(hjust = 0.5))

#So, all cells have a nucleus

#From what has already been investigated, cells with negative benchmarks and extremely large cells may not be suspicious, 
#but cells far outside the tissue will be suspicious. For some types of spatial neighborhood graphs, such as the k-nearest 
#neighbor graph, these cells will also affect the spatial analysis. Therefore, it is necessary to remove these cells. 

g <- findKNN(spatialCoords(sfe)[,1:2], k = 5, BNPARAM = AnnoyParam())
max_dist <- rowMaxs(g$distance)
data.frame(max_dist = max_dist) |> 
  ggplot(aes(max_dist)) +
  geom_histogram(bins = 100) +
  scale_y_continuous(transform = "log1p") +
  scale_x_continuous(breaks = breaks_pretty()) +
  annotation_logticks(sides = "l") + 
  ggtitle("Maximum distance between the neighbors") +
  theme(plot.title = element_text(hjust = 0.5))

min_dist <- rowMins(g$distance)
data.frame(min_dist = min_dist) |> 
  ggplot(aes(min_dist)) +
  geom_histogram(bins = 100) +
  scale_y_continuous(transform = "log1p") +
  scale_x_continuous(breaks = breaks_pretty()) +
  annotation_logticks(sides = "l") + 
  ggtitle("Minimum distance between the neighbors") +
  theme(plot.title = element_text(hjust = 0.5))

#There’s a long tail, which must be the small clusters of cells away from the tissue.

sfe$main_tissue <- !(max_dist > 95 | min_dist > 75)
plotSpatialFeature(sfe, "main_tissue", colGeometryName = "cellSeg",
                   show_axes = TRUE) +
  geom_sf(data = SpatialFeatureExperiment::centroids(sfe)[!sfe$main_tissue,], 
          color = ditto_colors[1], size = 0.3) + 
  ggtitle("Clusters of cells far from the tissue") +
  theme(plot.title = element_text(hjust = 0.5))

plotColDataHistogram(sfe, c("nCounts", "cell_area")) +
  scale_x_log10() +
  annotation_logticks(sides = "b") + 
  ggtitle("Histogram of cell clusters") +
  theme(plot.title = element_text(hjust = 0.5))

sfe <- sfe[,sfe$main_tissue & sfe$nCounts > 10]
sfe <- sfe[rowSums(counts(sfe)) > 0,]
dim(sfe)

#In total, 292 cells were removed

#Genes
rowData(sfe)$means <- rowMeans(counts(sfe))
rowData(sfe)$vars <- rowVars(counts(sfe))

rowData(sfe)$is_neg <- rowData(sfe)$Type != "Gene Expression"
plotRowData(sfe, x = "means", y = "is_neg") +
  scale_y_log10() +
  annotation_logticks(sides = "b") + 
  ggtitle("Comparison of gene expression and negative controls in cells ") +
  theme(plot.title = element_text(hjust = 0.5))

plotRowData(sfe, x="means", y="vars", color_by = "is_neg") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks() +
  coord_equal() +
  labs(color = "Negative control") + 
  ggtitle("Comparison of gene expression and negative controls in cells ") +
  theme(plot.title = element_text(hjust = 0.5))

# Negative controls and real genes form mostly separate clusters. Negative controls stick close to the line, 
#while real genes are overdispersed. 

#Spatial autocorrelation of QC metrics

colGraph(sfe, "knn5") <- findSpatialNeighbors(sfe, method = "knearneigh", 
                                              dist_type = "idw", k = 5, 
                                              style = "W")
sfe$nucleus_area <- ifelse(is.na(sfe$nucleus_area), 0, sfe$nucleus_area)
sfe <- colDataMoransI(sfe, c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                      colGraphName = "knn5")
colFeatureData(sfe)[c("nCounts", "nGenes", "cell_area", "nucleus_area"),]

#Global Moran’s I indicates positive spatial autocorrelation. (nCounts = 0.41, nGenes = 0.51, cell_area = 0.46, nucleus_area = 0.21)

sfe <- colDataUnivariate(sfe, type = "localmoran", 
                         features = c("nCounts", "nGenes", "cell_area", 
                                      "nucleus_area"),
                         colGraphName = "knn5", BPPARAM = SnowParam(2))
plotLocalResult(sfe, "localmoran",
                features = c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                colGeometryName = "centroids", scattermore = TRUE,
                divergent = TRUE, diverge_center = 0, pointsize = 1.5,
                ncol = 1)

localResultAttrs(sfe, "localmoran", "nCounts")

plotLocalResult(sfe, "localmoran", attribute = "pysal",
                features = c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                colGeometryName = "centroids", scattermore = TRUE, pointsize = 1.5,
                ncol = 1, size = 3)

bbox4 <- c(xmin = 3500, xmax = 5500, ymin = -2450, ymax = -25)

plotLocalResult(sfe, "localmoran",
                features = c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                colGeometryName = "cellSeg", divergent = TRUE, diverge_center = 0,
                ncol = 2, bbox = bbox4)

plotLocalResult(sfe, "localmoran", attribute = "pysal",
                features = c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                colGeometryName = "cellSeg", bbox = bbox4, ncol = 2)

plotLocalResult(sfe, "localmoran", attribute = "-log10p_adj", 
                features = c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                colGeometryName = "cellSeg", divergent = TRUE, 
                diverge_center = -log10(0.05),
                ncol = 2, bbox = bbox4)

#Here in this dense region, for these metrics, local Moran’s I is generally significant.
#It seems that locally we may have a mixture of positive and negative spatial autocorrelation.

sfe <- colDataUnivariate(sfe, "moran.plot", "nCounts", colGraphName = "knn5")
moranPlot(sfe, "nCounts", binned = TRUE, plot_influential = FALSE, bins = c(150, 80)) + 
  ggtitle("Moran plot for nCounts") +
  theme(plot.title = element_text(hjust = 0.5)) 

#There are no obvious clusters here based on the 2D histogram showing point density in the scatter plot. 

#Moran’s I

sfe <- logNormCounts(sfe, size.factor = sfe$cell_area)
system.time(
  sfe <- runMoransI(sfe, colGraphName = "knn5", BPPARAM = SnowParam(2))
)
plotRowData(sfe, x = "moran_sample01", y = "is_neg") + 
  ggtitle("Moran’s I for negative controls") +
  theme(plot.title = element_text(hjust = 0.5)) 

#Generally the negative controls are tightly clustered around 0, while the real genes have positive 
#Moran’s I, which means there is generally no technical artifact spatial trend.

#What are the genes with the highest Moran’s I?

top_moran <- rownames(sfe)[order(rowData(sfe)$moran_sample01, decreasing = TRUE)[1:4]]
plotSpatialFeature(sfe, top_moran, colGeometryName = "centroids",
                   scattermore = TRUE, pointsize = 1.5, ncol = 2) +
  plot_annotation(title = "Genes with the highest Moran’s I") &
  theme(plot.title = element_text(hjust = 0.5))

top_moran2 <- rownames(sfe)[order(rowData(sfe)$moran_sample01, decreasing = TRUE)[5:8]]
plotSpatialFeature(sfe, top_moran2, colGeometryName = "centroids",
                   scattermore = TRUE, pointsize = 1.5, ncol = 2) +
  plot_annotation(title = "Genes with the highest Moran’s I") &
  theme(plot.title = element_text(hjust = 0.5))

#How does Moran’s I relate to gene expression level?

plotRowData(sfe, x = "means", y = "moran_sample01", color_by = "is_neg") +
  scale_x_log10() + annotation_logticks(sides = "b")

#Very highly expressed genes have higher Moran’s I, but there are some less expressed genes with higher Moran’s I as well.

sfe <- sfe[rowData(sfe)$Type == "Gene Expression",]

#Non-spatial dimension reduction and clustering

set.seed(29)
sfe <- runPCA(sfe, ncomponents = 30, scale = TRUE, BSPARAM = IrlbaParam(), ntop = Inf)
ElbowPlot(sfe, ndims = 30) +
  plot_annotation(title = "PC elbow plot") &
  theme(plot.title = element_text(hjust = 0.5))

plotDimLoadings(sfe, dims = 1:8) +
  plot_annotation(title = "Plot top PC loadings of genes") &
  theme(plot.title = element_text(hjust = 0.5))

spatialReducedDim(sfe, "PCA", 4, colGeometryName = "centroids", divergent = TRUE,
                  diverge_center = 0, ncol = 1, scattermore = TRUE, pointsize = 1.5) +
  plot_annotation(title = "PCA") &
  theme(plot.title = element_text(hjust = 0.5))

spatialReducedDim(sfe, "PCA", components = 5:8, colGeometryName = "centroids", 
                  divergent = TRUE, diverge_center = 0, 
                  ncol = 1, scattermore = TRUE, pointsize = 1.5) +
  plot_annotation(title = "PCA") &
  theme(plot.title = element_text(hjust = 0.5))

#Non-spatial clustering and locating the clusters in space

colData(sfe)$cluster <- clusterRows(reducedDim(sfe, "PCA")[,1:15],
                                    BLUSPARAM = SNNGraphParam(
                                      cluster.fun = "leiden",
                                      cluster.args = list(
                                        resolution_parameter = 0.5,
                                        objective_function = "modularity")))

plotPCA(sfe, ncomponents = 4, colour_by = "cluster", scattermore = TRUE) +
  plot_annotation(title = "PCA") &
  theme(plot.title = element_text(hjust = 0.5))

plotSpatialFeature(sfe, "cluster", colGeometryName = "centroids", scattermore = TRUE,
                   pointsize = 1.2, size = 3) +
  plot_annotation(title = "Plot the location of the clusters in space") &
  theme(plot.title = element_text(hjust = 0.5))

#Since the principal components are spatially structured, the clusters found from the PCs are spatially structured.

#Differential expression

#Cluster marker genes are found with Wilcoxon rank sum test

markers <- findMarkers(sfe, groups = colData(sfe)$cluster,
                       test.type = "wilcox", pval.type = "all", direction = "up")
markers[[6]]

genes_use <- vapply(markers, function(x) rownames(x)[1], FUN.VALUE = character(1))
plotExpression(sfe, genes_use, x = "cluster", point_fun = function(...) list()) +
  plot_annotation(title = "Significant markers for each cluster") &
  theme(plot.title = element_text(hjust = 0.5))

genes_use2 <- unique(unlist(lapply(markers, function(x) rownames(x)[1:5])))
plotGroupedHeatmap(sfe, genes_use2, group = "cluster", colour = scales::viridis_pal()(100))

#Local spatial statistics of marker genes

plotSpatialFeature(sfe, genes_use, colGeometryName = "centroids", ncol = 3,
                   pointsize = 0, scattermore = TRUE) +
  plot_annotation(title = "Significant markers genes in space as a reference") &
  theme(plot.title = element_text(hjust = 0.5))

#Global Moran’s I of these marker genes

setNames(rowData(sfe)[genes_use, "moran_sample01"], genes_use)

#All these marker genes have positive spatial autocorrelation, but some stronger than others.

#Local Moran’s I of these marker genes 

sfe <- runUnivariate(sfe, "localmoran", features = genes_use, colGraphName = "knn5",
                     BPPARAM = SnowParam(2))
plotLocalResult(sfe, "localmoran", features = genes_use, 
                colGeometryName = "centroids", ncol = 3, divergent = TRUE,
                diverge_center = 0, scattermore = TRUE, pointsize = 0)

