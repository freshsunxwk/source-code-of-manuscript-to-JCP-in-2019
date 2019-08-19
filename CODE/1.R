rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE53986', destdir=".",
                 AnnotGPL = F,
                 getGPL = F)
  save(eSet,file='GSE53986_eSet.Rdata')
}
load('GSE3325_eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
phe=pData(b)
library(stringr)
group_list= str_split(as.character(phe$title),',',simplify = T)[,2]
save(raw_exprSet,group_list,
     file='GSE3325_raw_exprSet.Rdata')

library("affy")
library("limma")
library("affyPLM")
library("RColorBrewer")
library("pheatmap")
library("affycoretools")
library("simpleaffy")
library("BiocGenerics", lib.loc="D:/R-3.5.0/library")

rm(list=ls())
if(F){
  library(GEOquery)
  eSet <- getGEO('GSE9106', destdir=".",
                 AnnotGPL = T,
                 getGPL = F)
  save(eSet,file='GSE9106_eSet.Rdata')
}
load('GSE17708_eSet.Rdata')

keytypes("mouse4302.db")
keytypes(org.Hs.eg.db)
GPL5356 <- read.delim("D:/DEG/2018-11-29day1/result-r-data1/GPL5356.csv")

GPL=GPL[,c(1,3)]

colnames(GPL)=c("probe_id","symbol")  
ids=GPL

data_RMA[1]

2986---2----94

b = eSet[[1]]
raw_exprSet=exprs(b) 

phe=pData(b)


raw_exprSet=exprs(data_mas5)
raw_exprSet=exprs(data_RMA)
raw_exprSet[1:4,1:4]
pData(data_RMA)=phe


phe
library(stringr)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,3]
save(raw_exprSet,group_list,phe,ids,
     file='GSE8021_raw_exprSet_IDO.Rdata')

group_list= str_split(as.character(phe$characteristics_ch1),';',simplify = T)[,1]

phe$title
library(stringr)
group_list= str_split(as.character(phe$title),' ',simplify = T)[,11]
group_list=paste0('group',group_list,'h')
save(raw_exprSet,group_list,
     file='GSE9106_raw_exprSet.Rdata')

load(file='GSE17708_raw_exprSet.Rdata')
library(hgu133plus2.db)
eg2probe=toTable(hgu133plus2SYMBOL)
eg2probe[eg2probe$symbol=='TRAF4',]
raw_exprSet[1:4,1:4]
exprSet=log2(raw_exprSet)
dat=data.frame(gene= exprSet['211899_s_at',] ,
               mut= group_list)
head(dat)
if(require('ggpubr')){
  library(ggpubr)
  # google search : ggpubr boxplot add p-value
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  p <- ggboxplot(dat, x = "mut", y = "gene",
                 color = "mut", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means(method = "t.test")
}

if(require('ggstatsplot')){
  library(ggstatsplot)
  ggbetweenstats(data = dat, x = mut,  y = gene)
}


if(require('ggplot2')){
  library(ggplot2)
  ggplot(dat,aes(x=mut,y=gene))+
    geom_boxplot()+
    theme_bw()
}