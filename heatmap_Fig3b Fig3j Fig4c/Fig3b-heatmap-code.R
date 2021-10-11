library(pheatmap)  
data = read.table("fig3b-heatmap.txt",header=T,row.names=1,sep="\t")
exp = data[,1:9] 
classify=read.table("fig3b-heatmap-classify.txt",head=T,row.names=1,sep="\t") 
mycol<-colorRampPalette(c(#1860a8, #f0f0f0, #c03030),bias=1)(30)
p<-pheatmap(exp,  scale = "row",
             border =NA,
             color = mycol,
             cluster_rows=F,cluster_cols=F,
             show_rownames=T,show_colnames=T,
             cellheight = 7,cellwidth = 20,
             fontsize=7,fontsize_col=8,
             angle_col = 45,  
             treeheight_row=25, treeheight_col=20,
             gaps_col=3,
             gaps_row=22,
             annotation_row=classify
             
)
p