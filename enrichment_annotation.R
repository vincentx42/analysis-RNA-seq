####Functional and pathway enrichment (GO/KEGG)####
rm(list = ls())
load(file = './data/deg.Rdata')
head(deg)

threshold = 1.5
deg$group = ifelse(deg$P.Value>0.05,'stable',
               ifelse(deg$logFC > threshold,'up',
                       ifelse(deg$logFC < -threshold,'down','stable') )
)

table(deg$group)
deg$symbol=rownames(deg)
head(deg)


library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG = deg
head(DEG)

DEG = merge(DEG, df, by.y = 'SYMBOL',by.x = 'symbol')
head(DEG)
save(DEG, file = './data/anno_DEG.Rdata')

gene_up = DEG[DEG$group == 'up', 'ENTREZID'] 
gene_down = DEG[DEG$group == 'down', 'ENTREZID'] 
gene_diff = c(gene_up, gene_down)
gene_all = as.character(DEG[, 'ENTREZID'])

data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList = DEG$logFC
names(geneList) = DEG$ENTREZID
geneList=sort(geneList,decreasing = T)


## KEGG pathway analysis
kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red", guide = 'none') + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}

if(T){
  ###   over-representation test
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff = 0.9)
  head(kk.up)[,1:6]
  dotplot(kk.up );ggsave('./results/kk.up.dotplot.png')
  
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  dotplot(kk.down );ggsave('./results/kk.down.dotplot.png')
  
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  dotplot(kk.diff );ggsave('./results/kk.diff.dotplot.png')
  
  # transfer to dataframe
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  
  
  g_kegg = kegg_plot(up_kegg, down_kegg)
  print(g_kegg)
  ggsave(g_kegg, filename = './results/kegg_up_down.png')
  
  ###  GSEA 
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 120,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  
  down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = './results/kegg_up_down_gsea.png')
  
  
}

### GO database analysis 
g_list=list(gene_up = gene_up,
            gene_down = gene_down,
            gene_diff = gene_diff)

if(T){
  go_enrich_results <- lapply(g_list, function(gene) {
    lapply( c('BP','MF','CC') , function(ont) {
      cat(paste('Now process ',ont ))
      ego <- enrichGO(gene          = gene,
                      universe      = gene_all,
                      OrgDb         = org.Hs.eg.db,
                      ont           = ont ,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.99,
                      qvalueCutoff  = 0.99,
                      readable      = TRUE)
      
      print( head(ego) )
      return(ego)
    })
  })
  save(go_enrich_results,file = './data/go_enrich_results.Rdata')
  
}

load(file = './data/go_enrich_results.Rdata')

n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
for (i in 1:3){
  for (j in 1:3){
    fn=paste0('./results/dotplot_',n1[i],'_',n2[j],'.png')
    cat(paste0(fn,'\n'))
    png(fn,res=150,width = 1080)
    print(dotplot(go_enrich_results[[i]][[j]]))
    dev.off()
  }
}
