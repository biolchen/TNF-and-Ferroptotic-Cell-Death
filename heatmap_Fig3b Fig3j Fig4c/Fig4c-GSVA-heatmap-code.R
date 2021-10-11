library(pheatmap)
data = read.table("Fig4c-GSVA-heatmap.txt",header=T,row.names=1,sep="\t",check.names = F)
mycol<-colorRampPalette(c("navy", "white", "firebrick3"))(100)

p<-pheatmap(data,
            scale = "row",
            border ="NA",
            color = mycol,
            cluster_rows=T,cluster_cols=T,
            show_rownames=T,show_colnames=T,
            cellheight =12,cellwidth = 30,
            fontsize=10,fontsize_col=15,
            angle_col = 45,
            
)