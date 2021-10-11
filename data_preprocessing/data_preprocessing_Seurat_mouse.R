####################################
# sc.mm.gam.gender                 #
# data preprocessing and analysis  #
####################################

# Load packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

library(future)
options(future.globals.maxSize = 100 * 1024 * 1024^2)
plan("multiprocess", workers = 4)

### Creat Seurat Object
message( "==>Reading 10x data<==" )
##download GSE129087
##rawdata  contains : barcodes.tsv,genes.tsv,matrix.mtx
samples <- c("R1","R2","R3")

mat <- Read10X(data.dir = "rawdata/", gene.column = 1)
obj <- CreateSeuratObject(counts = mat , assay = "RNA" )
obj@meta.data$orig.ident <- gsub("[A-Z]+-","",rownames(obj@meta.data))
obj@meta.data$orig.ident <- factor(obj@meta.data$orig.ident)
obj <- SetIdent(obj, value = "orig.ident")

color.sample <- scales::hue_pal()(length(levels(obj@meta.data$orig.ident)))
names(color.sample) <- levels(obj@meta.data$orig.ident)
obj@misc[["color.sample"]] <- color.sample


##add Feature Data
ref_name_file <- "name_list.xls"
fdata <- read.table(ref_name_file, row.names = 1, stringsAsFactors = F, sep = "\t")
colnames(fdata) <- c("merge_name", "name", "type")
fdata$merge_name <- fdata$name
fdata$merge_name[fdata$merge_name == "-"] <- rownames(fdata)[fdata$merge_name == "-"]
index <- c(which(duplicated(fdata$merge_name, fromLast=T)), which(duplicated(fdata$merge_name, fromLast=F)))
fdata$merge_name[index] <- paste0(fdata$merge_name[index], " (", rownames(fdata)[index], ")")
fdata <- AddUnderscore(fdata)
obj@misc[["fdata"]] <- fdata

##add p data
obj@misc[["pdata"]] <- FetchData(obj, c("orig.ident", "nFeature_RNA", "nCount_RNA"))

##add mito.percent
StatFeatures <- function(object, features = NULL, col.name = NULL, stat_pct = FALSE, assay = NULL, add_to_pdata = FALSE){
		if ( file.exists(features[1]) )
				features <- readLines(con = features[1])
		features <- FindFeaturesID(object = object, features = features, unlist = FALSE)
		if ( is.null(assay) )
				assay <- DefaultAssay(object = object)
		metadata <- Matrix::colSums(x = GetAssayData(object = object, slot = "counts", assay = assay)[unlist(features), , drop = FALSE])
		if (stat_pct)
				metadata <- metadata / object@meta.data[[paste0("nCount_", assay)]] * 100

		if (!is.null(x = col.name)) {
				object@misc[[col.name]] <- features
				if (add_to_pdata)
						object@misc$pdata[[col.name]] <- metadata
				object <- AddMetaData(object = object, metadata = metadata, col.name = col.name)
				return(object)
		} else {
				return(metadata)
		}
}

mito_list <- c("MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5","MT-ND6")
obj <- StatFeatures(obj, mito_list, col.name = "percent.mito", stat_pct = T, add_to_pdata = T)


###Filter Cells
cells.use <- Cells(obj)
#nCount_RNA: [ -.inf,20000 ]
standard <- c(-.inf,20000)
cells.use <- obj@meta.data %>% tibble::rownames_to_column(var = "cells") %>% filter(.data[["nCount_RNA"]] >= standard[1] & .data[["nCount_RNA"]] <= standard[2] & cells %in% cells.use) %>% select(cells) %>% unlist()
#nFeature_RNA: [ 500,5000 ]
standard <- c(500,5000)
cells.use <- obj@meta.data %>% tibble::rownames_to_column(var = "cells") %>% filter(.data[["nFeature_RNA"]] >= standard[1] & .data[["nFeature_RNA"]] <= standard[2] & cells %in% cells.use) %>% select(cells) %>% unlist()
#percent.mito: [ -.inf,20 ]
standard <- c(-.inf,20)
cells.use <- obj@meta.data %>% tibble::rownames_to_column(var = "cells") %>% filter(.data[["percent.mito"]] >= standard[1] & .data[["percent.mito"]] <= standard[2] & cells %in% cells.use) %>% select(cells) %>% unlist()

##Do Filt
obj <- obj[,cells.use]

##Cell Cycle
obj <- DoCellCycleScoring(obj)

### Normalization Data
normalization.method <- "LogNormalize"
scale.factor <- 10000
nfeatures <- 2000

DefaultAssay(obj) <- "RNA"

vars.regress <- c("nCount_RNA","percent.mito","CC.Difference")
obj@misc$vars.regress <- vars.regress

obj <- NormalizeData(obj, normalization.method = normalization.method, scale.factor = scale.factor)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
obj <- ScaleData(object = obj, vars.to.regress = obj@misc$vars.regress, features = rownames(obj))


### Integration
message( "==> Do Integration <==" )
pca_dims <- 50               # number of PCA dimensions to compute and use in tSNE, UMAP
dims <- seq(pca_dims)
anchor.features <- 3000
k.filter <- min(200, ceiling(min(sapply(object.list, ncol))/2))

old.assay <- DefaultAssay(obj)
anchors <- FindIntegrationAnchors(object.list = object.list, dims = dims, normalization.method = normalization.method, anchor.features = anchor.features, k.filter = k.filter)
integrated <- IntegrateData(anchorset = anchors, dims = dims, normalization.method = normalization.method)
integrated <- ScaleData(integrated, verbose = FALSE)
integrated@misc <- object@misc
integrated[[old.assay]] <- object[[old.assay]]
integrated@reductions <- object@reductions


### tSNE & UMAP 
pc.num <- 50
sig.PCs <- seq(pc.num)

obj <- integrated
assay <- DefaultAssay(obj)  # integrated
#run PCA
obj <- RunPCA(object = obj, assay = assay, npcs = pc.num, features = VariableFeatures(obj), verbose = FALSE )
obj[[paste0("pca_", assay)]] <- obj[["pca"]]
#run tSNE
obj <- RunTSNE(obj, dims = sig.PCs)
obj[[paste0("tsne_", assay)]] <- obj[["tsne"]]
#run UMAP
obj <- RunUMAP(obj, dims = sig.PCs, umap.method = "uwot")
obj[[paste0("umap_", assay)]] <- obj[["umap"]]


### Find clusters
clustering_resolution <- 0.5   # Resolution parameter for Seurat clustering	

reduction <- "pca"
dims <- seq(obj[["pca"]])
Idents(obj) <- "seurat_clusters"

obj <- FindNeighbors(obj, reduction = reduction, dims = dims, force.recalc = TRUE)
obj <- FindClusters(obj, resolution = clustering_resolution , temp.file.location = getwd())

color.cluster <- rainbow(length(levels(obj@meta.data$seurat_clusters)))
names(color.cluster) <- levels(obj@meta.data$seurat_clusters)
obj@misc[["color.cluster"]] <- color.cluster

###Find maker genes
min.pct <- 0.25    #minimum percent of gene-expressed cells 
logfc <- 0.25      #log Fold change
pvalue <- 0.01     #pvalue threshold

Idents(object) <- "seurat_clusters"
object.markers <- FindAllMarkers(object = object, only.pos = TRUE,min.pct = min.pct, logfc.threshold = logfc,return.thresh = pvalue, pseudocount.use = 0 )

#Find top5 marker
top_num <- 5
top <- object.markers %>% group_by( cluster ) %>%
		arrange(desc(avg_logFC), p_val, p_val_adj, .by_group = TRUE) %>% filter(1:n() <= top_num)


### Save data object 
message( "==>Output obj.Rda<==" )
DefaultAssay(obj) <- "RNA"
save(obj, file = "obj.Rda")
save(obj.markers, file = "markers.Rda")

message( "==>All Done!<==" )

