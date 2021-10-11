library(ggplot2)

all_pval <- read.table("pvalues.txt",head=T,sep="\t",check.name=0,comment.char = "",quote="")
all_means <- read.table("means.txt",head=T,sep="\t",check.name=0,comment.char = "",quote="")
all_pval <- all_pval[,-c(1:11)]
all_means <- all_means[,-c(1:11)]
intr_pairs <- all_pval$interacting_pair

selected_rows <- c("TNF_TNFRSF1A","TNF_TNFRSF1B","TGFB1_TGFbeta receptor1","TGFB1_TGFbeta receptor2","IL1 receptor_IL1A","IL1 receptor_IL1B","IL6 receptor_IL6","IL7 receptor_IL7","IL10 receptor_IL10","IL17 receptor AC_IL17A","IL15 receptor_IL15","PDGFA_PDGFRA","FGF1_FGFR1")
selected_columns <- c("macrophages|macrophages","macrophages|non-sensitive","macrophages|sensitive","non-sensitive|macrophages","non-sensitive|non-sensitive","non-sensitive|sensitive","sensitive|macrophages","sensitive|non-sensitive","sensitive|sensitive")

sel_pval <- all_pval[match(selected_rows, intr_pairs), selected_columns]
sel_means <- all_means[match(selected_rows, intr_pairs), selected_columns]

df_names <- expand.grid(selected_rows, selected_columns)
pval <- unlist(sel_pval)
pval[pval==0] <- 0.0009
pr = unlist(as.data.frame(sel_means))
pr[pr==0] <- 1

plot.data <- cbind(df_names,pval,log2(pr))
colnames(plot.data) <- c('pair', 'clusters', 'pvalue', 'mean')

my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

p <- ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette)+
    labs(x = "", y = "")+ theme_bw() +
    theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    theme(axis.text.x = element_text(size=16,angle=90,hjust=1),axis.text.y = element_text(size=16),
    legend.title = element_text(size = 16))

ggsave(p,filename = "ligand.pdf", width = 14, height = 10)
