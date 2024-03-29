rm(list=ls()) 
rm(ids)


load("D:/DEG/2018-11-29day1/result-r-data1/GSE8021_raw_exprSet_IDO.Rdata")

load(file='GSE17708_raw_exprSet.Rdata')
#exprSet=log(exprSet)
exprSet=raw_exprSet
library(hgu133plus2.db)
library(hgu133plus2.db)
hgu133a2.db
source("http://bioconductor.org/biocLite.R");
biocLite(hgu133a2.db )

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hgu133a2.db", version = "3.8")
library("hgu133a2.db")
ids=toTable(mouse4302SYMBOL)
ids=toTable(hgu133a2SYMBOL)
#length(unique(ids$symbol))
#tail(sort(table(ids$symbol)))
#table(sort(table(ids$symbol)))
#plot(table(sort(table(ids$symbol))))
#判断前面一个向量内的元素是否在后面一个向量中，返回布尔值
exprSet=raw_exprSet
table(rownames(exprSet) %in% ids$probe_id)
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)

ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
exprSet[1:5,1:5]

jimmy <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}

new_exprSet <- jimmy(exprSet,ids)

save(new_exprSet,ids,group_list,
     file='GSE9106_new_exprSet.Rdata')


load(file='GSE8021_new_exprSet.Rdata')



#
new_exprSet=diffgene_exprset[hub$x,]
colnames(new_exprSet)=group_list
write.csv(new_exprSet,"93sample_93hubgene.csv")
#






exprSet=new_exprSet
if(T){
  
  library(reshape2)
  exprSet_L=melt(exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  
  
  #group_list=as.character(datTraits[,2])
  
  
  exprSet_L$group=rep(group_list,each=nrow(exprSet))
  head(exprSet_L)
  ### ggplot2
  library(ggplot2)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
  print(p)
  
  p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p) 
  ## hclust
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),
                  cex = 0.7, col = "blue")
  hc=hclust(dist(t(exprSet)))
  par(mar=c(5,5,5,10))
  png('hclust.png',res=120)
  plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
  dev.off()
  
  ## PCA
  
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_list
  png('pca11.png',res=120)
  autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')+theme_bw()
  dev.off()
}