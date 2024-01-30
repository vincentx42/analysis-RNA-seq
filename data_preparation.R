####Download data from GEO####
library(GEOquery)
options(stringsAsFactors = F)
# check file integrity 
f='./data/GSE42872_eSet.Rdata'
if(!file.exists(f)){
  gset <- getGEO('GSE42872', destdir = "./data",
                 AnnotGPL = F,     ## annotation
                 getGPL = F)       ## platform
  save(gset,file = f)   ## save to local folder
}

# load the Rdata
load('./data/GSE42872_eSet.Rdata')
class(gset)  
length(gset)  
class(gset[[1]])

# Extract information
a = gset[[1]]
str(a)
dat = exprs(a)
dim(dat)

dat[1:4,1:4] 
boxplot(dat,las=2)
pd = pData(a) ## clinical information

library(stringr)
group_list = str_split(pd$title, ' ', simplify = T)[,4]
table(group_list)

# another way to load expression matrix
dat_2 = read.table('./data/GSE42872_series_matrix.txt.gz',
                   sep = '\t',
                   quote = '',
                   fill = T,
                   comment.char = '!')
class(dat_2)
str(dat_2)


####Tranfer gene ID####
if(F){
  library(GEOquery)
  gpl <- getGEO('GPL6244', destdir="./data")
  colnames(Table(gpl))  
  head(Table(gpl)[,c(1,15)])
  probe2gene=Table(gpl)[,c(1,15)]
  head(probe2gene)
  library(stringr)  
  save(probe2gene,file='probe2gene.Rdata')
}
#
# load(file='probe2gene.Rdata')
# ids=probe2gene 


BiocManager::install("hugene10sttranscriptcluster.db")
library(hugene10sttranscriptcluster.db)
#extract relationship between probe_id and gene symbol
ids = toTable(hugene10sttranscriptclusterSYMBOL) 
head(ids) 
summary(ids)
colnames(ids)=c('probe_id','symbol')  
ids = ids[ids$symbol != '', ]
ids = ids[ids$probe_id %in%  rownames(dat), ]

dat = dat[ids$probe_id, ] 

ids$median = apply(dat, 1, median) 
ids=ids[order(ids$symbol,ids$median,decreasing = T), ]
ids=ids[!duplicated(ids$symbol),]
dat=dat[ids$probe_id,] 
rownames(dat)=ids$symbol
save(dat, group_list, file = './data/step1-2-output.Rdata')


####Step3 check your expression matrix####
rm(list = ls())
options(stringsAsFactors = F)
load(file = './data/step1-2-output.Rdata')

# Expression level check for common genes
dat['ZZZ3', ]
boxplot(dat[, 1])

# Distribution check
library(reshape2)
exprSet_L = melt(exprSet)
colnames(exprSet_L) = c('probe', 'sample', 'value')
exprSet_L$group = rep(group_list, each=nrow(exprSet))

ggplot(exprSet_L, aes(x=sample, y=value, fill=group))+
  geom_boxplot()

ggplot(exprSet_L, aes(value, fill=group))+
  geom_histogram()


# Initial clustering
exprSet = dat
colnames(exprSet) = paste(group_list, 1:6, sep='')
nodePar <- list(lab.cex = 0.6, pch=c(NA, 19),
                cex=0.7, col='blue')
hc <- hclust(dist(t(exprSet)))
par(mar=c(5,5,5,10))
plot(as.dendrogram(hc), nodePar=nodePar, horiz=T)

# PCA
dat=t(dat)
dat=as.data.frame(dat)
dat=cbind(dat,group_list) 
library("FactoMineR")

dat.pca <- PCA(dat[, -ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
ggsave('./results/all_samples_PCA.png')

# Heatmap
cg = names(tail(sort(apply(exprSet, 1, sd)), 1000))
library(pheatmap)
pheatmap(exprSet[cg, ], show_colnames=F, show_rownames=F) 

n=t(scale(t(exprSet[cg, ]))) # 'scale' can normalise the data
n[n>2]=2 
n[n< -2]= -2
n[1:4, 1:4]
pheatmap(n, show_colnames=F, show_rownames=F)

ac=data.frame(g=group_list)
rownames(ac)=colnames(n)

pheatmap(n, show_colnames=F, show_rownames=F,
         annotation_col=ac, filename = './results/heatmap_top1000_sd.png')

