####################################
# sc.mm.gam.gender                 #
# Figures - revision               #
####################################

#load load 1.data_pre obj file 
load("obj.Rda")

###Draw tSNE
#Figure 3a  cus_1
p1 <- DimPlot(obj, reduction = "tsne", group.by = "orig.ident", label = TRUE, cols = obj@misc$color.sample)
ggsave(p1, filename = "tSNE_3a.1.pdf", width = 10, height = 10, limitsize = F)

p2 <- DimPlot(obj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE, cols = obj@misc$color.cluster)
ggsave(p2, filename = "tSNE_3a.2.pdf", width = 10, height = 10, limitsize = F)

#Figure 3c  sup_2
relation <- c("Sensitive","Non-sensitive","Non-sensitive","Non-sensitive","Non-sensitive","Non-sensitive","Non-sensitive","Sensitive","Sensitive")
names(relation) <- c("0","2","3","4","5","6","7","9","10")
color <- c("#00F4F4","#FE0000")
names(color) <- c("Sensitive","Non-sensitive")

sub_obj <- subset(obj,idents=names(relation))
sub_obj@meta.data$type <- relation[as.vector(sub_obj@meta.data$seurat_clusters)]

p1 <- DimPlot(sub_obj, reduction = "tsne", group.by = "type", label = TRUE, cols = color)
ggsave(p1, filename = "tSNE_3c.pdf", width = 10, height = 10, limitsize = F)

#Figure 3d  sup_2
marker <- c("Prg4","Sparcl1","Inhba","Col11a1","Mfap4","Sfrp2")
color <- c("#3871B3","#FB0005")
names(color) <- c("Sensitive","Non-sensitive")
p2 <- VlnPlot(object = sub_obj, features = marker , group.by = "type", pt.size = 0, combine = TRUE,nCol = 3, cols = color)
ggsave(p2, filename = "Vlnplot_3d.pdf", width = 10, height = 7, limitsize = F)


#function 
drawtSNE <- function(data, group.by = "seurat_clusters" , cols = NULL, outfile = NULL){

	p <- ggplot(data = data, mapping = aes(x = "tSNE_1" , y = "tSNE_2")) + 
		geom_point(aes_string(color = group.by )) +
		scale_color_manual(colors = cols, na.value = "lightgrey") +
		theme_classic() +
		theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())

	if( ! is.null(outfile)){
		ggsave(p, filename = outfile, width = 10, height = 10, limitsize = F)
	}else{
		return(p)
	}
}

#Figure 4a
infile <- "Fig4a.data.xls"
data <- read.table(infile, row.names = 1, stringsAsFactors = F, sep = "\t", header = T)

color <- c("#C6C5E8","#8052d8","#FB73E","#F0B142","#B02F96","#F14149","#28C580","#1CE8BF","#33ccff","#F3854F")
names(color) <- c("B cells","DC ","Endothelial","fibroblast","Macrophages","macrophages/Monocytes","Monocytes","MSC","pericyte","T cells")
drawtSNE(data,group.by = "Celltype", cols = color, outfile = "tSNE_4a.pdf")

#Figure 4b
infile <- "Fig4bd.data.xls"
data <- read.table(infile, row.names = 1, stringsAsFactors = F, sep = "\t", header = T)
use.data <- subset(data, Cluster %in% c("Fib1","Fib2","Fib3","Fib4"))

color1 <- c("#F0B142","#0099ff","#F14149","#28C580")
names(color1) <- c("Fib1","Fib2","Fib3","Fib4")
drawtSNE(use.data,group.by = "Cluster", cols = color1, outfile = "tSNE_4b.pdf")

color2 <- c("#0033cc","#F14149")
names(color2) <- c("Fib a","Fib b")
relation <- c("Fib a","Fib a","Fib b","Fib b")
names(relation) <- c("Fib1","Fib2","Fib3","Fib4")

use.data$newType <- relation[as.vector(use.data$Cluster)]

drawtSNE(use.data,group.by = "newType", cols = color2, outfile = "tSNE_4d.pdf")

message( "==>All Done!<==" )
