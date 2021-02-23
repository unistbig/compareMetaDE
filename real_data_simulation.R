#filtering step은 없어도 될 것.
#--------------------importing data into R--------------------------------------------------------#
library(MetaDE)
Sys.setenv(JAVA_HOME='~/java/jdk1.8.0_261')
Sys.setenv(CLASSPATH="~/java/jdk1.8.0_261/jre/lib/ext")
library(tidyverse)
library(magrittr)
library(MetaQC)

args = commandArgs(trailingOnly = T)
if (length(args) < 4)
  stop("Must set 4 indexs", call. = F)
input_dir=args[1]
res_dir=args[2]
true_gene_threshold=as.numeric(args[3])
k = as.integer(args[4])




load(file.path(res_dir,paste0('truegenes_thr_',true_gene_threshold ,'.Rdata')))



setwd(input_dir)
dir.create(file.path(res_dir,'Meta_Res_real'),showWarnings = F)

study.names<-c("Welsh","Yu","Lapointe","Varambally","Singh","Wallace","Nanni","Dhanasekaran", "Tomlins")
prostate.raw<-MetaDE.Read(study.names,skip=rep(1,9),via="txt",matched=T,log=F)

library(magrittr)
study.good<-c("Welsh","Yu","Lapointe","Varambally","Singh","Wallace")
study.bad<-c("Nanni","Dhanasekaran", "Tomlins")
rep_ind<-1
study.now<-list()
for(i in 1:(length(study.good)-1)){
  for(j in (i+1):length(study.good)){
    study.now[[rep_ind]]<-c(study.good[c(i,j)])
    rep_ind%<>%+1
  }
}


prostate.resampled<-prostate.raw
for(i in 1:length(prostate.resampled)){
  prostate.resampled[[i]]$y<-sample(prostate.resampled[[i]]$y, size=length(prostate.resampled[[i]]$y))
}
# save(prostate.resampled, file=paste0(res_dir,'prostate.resampled.Rdata'))
#load('prostate.resampled.Rdata')

study.names<-c("Welsh","Yu","Lapointe","Varambally","Singh","Wallace","Nanni","Dhanasekaran", "Tomlins")


data.QC.raw<-list()
for(i in 1:length(study.names)){
  if(study.names[[i]] %in% study.now[[k]]){
    data.QC.raw[[i]]<-prostate.raw[[study.names[i]]]
  }else{
    data.QC.raw[[i]]<-prostate.resampled[[study.names[i]]]
  }
}

names(data.QC.raw)<-study.names
data.QC.merged<-MetaDE.merge(data.QC.raw)
data.QC.filtered<-MetaDE.filter(data.QC.merged,c(0.3,0.3))#for studyall
ndata=length(data.QC.filtered)
data.QC.filtered<-MetaDE:::MetaDE.impute(data.QC.filtered,y=0.3)
weight=rep(1,length(data.QC.filtered))
for (wi in 1:length(weight)) {
  weight[wi] = length(data.QC.filtered[[wi]]$y)
}
set.seed(1234)
r_num=sample(1:ndata, 15 ,replace=T)[i]
MetaDE.Res<-list()
MetaDE.Res.ind<-list()
MetaDE.Res.rankProd<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",ndata),meta.method=c("rankProd"),asymptotic=F, nperm = 300, metade.perm = F, na.rm = F)

for(ind.tail in c('high','low')){
  tailed=switch (ind.tail,
                 "low" = "left",
                 "high" = "right"
  )
  
  ind.Res1<-ind.analysis(data.QC.filtered,ind.method=rep("modt",ndata),nperm=300,tail=ind.tail,miss.tol = 0.3)
  ind.Res2<-ind.cal.ES(data.QC.filtered,paired=rep(F,ndata),nperm=300,miss.tol=0.3)
  
  MetaDE.Res.ind[['ind.modt']]<-ind.Res1
  MetaDE.Res.ind[['ind.ES']]<-ind.Res2
  MetaDE.Res.ind$x<-data.QC.filtered
  MetaDE.Res.ind$index<-k
  
  MetaDE.Res[["Fisher"]]<-MetaDE.pvalue(ind.Res1,meta.method=c("Fisher"),asymptotic = T)
  MetaDE.Res[["Stouffer"]]=MetaDE.pvalue(ind.Res1,meta.method=c("Stouffer"),asymptotic = T)
  MetaDE.Res[["wFisher"]]=MetaDE.pvalue(ind.Res1,meta.method=c("wFisher_sora"),asymptotic = T, weight = weight)
  MetaDE.Res[["ordmeta"]]=MetaDE.pvalue(ind.Res1,meta.method=c("ordmeta"),asymptotic = T)
  MetaDE.Res[["roP_2"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=2,asymptotic = T)
  MetaDE.Res[["roP_4"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=2,asymptotic = T)
  MetaDE.Res[["roP_6"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=6,asymptotic = T)
  MetaDE.Res[["roP_random"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=r_num,asymptotic = T)
  MetaDE.Res[["REM"]]<-MetaDE.ES(ind.Res2, meta.method = "REM", tail=ind.tail)
  MetaDE.Res[["FEM"]]<-MetaDE.ES(ind.Res2, meta.method = "FEM", tail=ind.tail)
  MetaDE.Res[["rankProd"]]<-MetaDE.Res.rankProd
  
  save(MetaDE.Res,file=file.path(res_dir,'Meta_Res_real',paste0('Real_' ,ndata ,'study_' ,k,'_rep_',tailed,'_tailed_Meta_Res.RData')))
  save(MetaDE.Res.ind,file=file.path(res_dir,'Meta_Res_real',paste0('Real_' ,ndata ,'study_' ,k,'_rep_',tailed,'_tailed_Meta_Res_ind.RData')))
}


