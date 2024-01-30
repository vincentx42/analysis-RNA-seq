####Screen DEGs (differentially expressed genes)####
# check distribution by boxplot
boxplot(dat[1,]~group_list) 

bp=function(g){         
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}

bp(dat[1,]) 
bp(dat[2,])
bp(dat[3,])
bp(dat[4,])
dim(dat)

library(limma)
# expression matrix
exprSet <-  dat

# design group matrix
design <- model.matrix(~0 + factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(exprSet)

# contrast matrix
contrast.matrix <- makeContrasts(paste0(unique(group_list)[2:1], collapse = '-'),
                                        levels = design)

# if manually set contrast>>
# contrast.matrix <- makeContrasts("Vemurafenib-Control", levels = design)
contrast.matrix

# screen out degs
deg = function(exprSet, design, contrast.matrix){
  fit1 <- lmFit(exprSet, design)
  fit2 <- contrasts.fit(fit1, contrast.matrix)
  fit3 <- eBayes(fit2)
  
  tempoutput <- topTable(fit3, coef = 1, n=Inf)
  nrDEG <- na.omit(tempoutput)
  # write.csv(nrDEG, "limma_notrend.results.csv", quote = F)

  return(nrDEG)
}

deg = deg(exprSet, design, contrast.matrix)
head(deg)
save(deg, file = './data/deg.Rdata')
write.csv(deg, file = './results/deg.csv')

# check some genes
exprSet['CD36', ]


####Visualization####
# Volcano plot and heatmap

## Volcano plot
nrDEG = deg
head(nrDEG)
attach(nrDEG)
plot(logFC, -log10(P.Value))

library(ggpubr)
df = nrDEG
df$v = -log10(P.Value) 
ggscatter(df, x = "logFC", y = "v", size = 0.5)

df$regulated = ifelse(df$P.Value > 0.01, 'stable', 
            ifelse(df$logFC > 2, 'up', 
                    ifelse(df$logFC < -2, 'down', 'stable') )
)

table(df$regulated)

df$name = rownames(df)
head(df)

ggscatter(df, x = "logFC", y = "v", size=0.5, color = 'regulated')
ggscatter(df, x = "logFC", y = "v", color = "regulated", size = 0.5,
          label = "name", repel = T,
          # label.select = rownames(df)[df$regulated != 'stable'] ,
          label.select = head(rownames(deg)), 
          palette = c("#00AFBB", "#E7B800", "#FC4E07"))
ggsave(filename = './results/volcano.png')


ggscatter(df, x = "AveExpr", y = "logFC", size = 0.2)
df$p_c = ifelse(df$P.Value < 0.001, 'p<0.001',
                ifelse(df$P.Value < 0.01,'0.001<p<0.01', 'p>0.01'))
table(df$p_c )
ggscatter(df,x = "AveExpr", y = "logFC", color = "p_c",size=0.2, 
          palette = c("green", "red", "black") )
ggsave('./results/MA.png')


## Heatmap 
if(T){ 
  load(file = './data/step1-2-output.Rdata')
  dat[1:4,1:4]
  table(group_list)
  x=deg$logFC 
  names(x)=rownames(deg) 
  cg=c(names(head(sort(x),100)),
       names(tail(sort(x),100)))
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
  n=t(scale(t(dat[cg,])))
  
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) 
  pheatmap(n,show_colnames =F,
           show_rownames = F,
           cluster_cols = F, 
           annotation_col=ac,filename = 'heatmap_top200_DEG.png')
  
}