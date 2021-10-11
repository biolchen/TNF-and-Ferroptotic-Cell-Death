library(pheatmap)
data = read.table("Fig3j-heatmap.txt",header=T,row.names=1,sep="\t",check.names = F)
dim(data)
head(data)
exp = data[,1:2]
dim(exp)
head(exp)
mycol<-colorRampPalette(c("navy", "white", "firebrick3"))(100)

data2=log2(exp+0.0001)
p<-pheatmap(data2,  
            border ="grey",
            color = mycol,
            cluster_rows=T,cluster_cols=F,
            show_rownames=T,show_colnames=T,
            cellheight = 16,cellwidth = 26,
            fontsize=16,fontsize_col=18,
)
p
