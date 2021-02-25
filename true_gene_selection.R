#filtering step은 없어도 될 것.
#--------------------importing data into R--------------------------------------------------------#
library(MetaDE)
Sys.setenv(JAVA_HOME='~/java/jdk1.8.0_261')
Sys.setenv(CLASSPATH="~/java/jdk1.8.0_261/jre/lib/ext")
library(tidyverse)
library(magrittr)
library(MetaQC)

args = commandArgs(trailingOnly = T)
if (length(args) < 2)
  stop("Must set 2 indexs", call. = F)
input_dir=args[1]
res_dir=args[2]
true_gene_threshold=as.numeric(args[3])

setwd(input_dir)
dir.create(file.path(res_dir,'Meta_Res_real'),showWarnings = F)

study.names<-c("Welsh","Yu","Lapointe","Varambally","Singh","Wallace","Nanni","Dhanasekaran", "Tomlins")
prostate.raw<-MetaDE.Read(study.names,skip=rep(1,9),via="txt",matched=T,log=F)

# ##################################################################################
# #9study
prostate.merged<-MetaDE.merge(prostate.raw)
dim(prostate.merged[[1]][[1]])
prostate.filtered<-MetaDE.filter(prostate.merged,c(0.3,0.3))

ndata=length(prostate.filtered)

prostate.filtered<-MetaDE:::MetaDE.impute(prostate.filtered,y=0.3)

weight=rep(1,length(prostate.filtered))
for (wi in 1:length(weight)) {
  weight[wi] = length(prostate.filtered[[wi]]$y)
}
MetaDE.Res<-list()
MetaDE.Res.ind<-list()

MetaDE.Res.rankProd<-MetaDE.rawdata(prostate.filtered,ind.method=rep("modt",ndata),meta.method=c("rankProd"),asymptotic=F, nperm = 300, metade.perm = F, na.rm=F)

for(ind.tail in c('high', 'low')){
  tailed=switch (ind.tail,
                 "low" = "left",
                 "high" = "right"
  )
  #ndata=length(study.names)



  # start<-Sys.time()
  ind.Res1<-ind.analysis(prostate.filtered,ind.method=rep("modt",ndata),nperm=300,tail=ind.tail,miss.tol = 0.3)
  ind.Res2<-ind.cal.ES(prostate.filtered,paired=rep(F,ndata),nperm=300,miss.tol=0.3)

  MetaDE.Res.ind[['ind.modt']]<-ind.Res1
  MetaDE.Res.ind[['ind.ES']]<-ind.Res2
  MetaDE.Res[["Fisher"]]<-MetaDE.pvalue(ind.Res1,meta.method=c("Fisher"),asymptotic = T)
  MetaDE.Res[["Stouffer"]]=MetaDE.pvalue(ind.Res1,meta.method=c("Stouffer"),asymptotic = T)
  MetaDE.Res[["wFisher"]]=MetaDE.pvalue(ind.Res1,meta.method=c("wFisher_sora"),asymptotic = T, weight = weight)
  MetaDE.Res[["ordmeta"]]=MetaDE.pvalue(ind.Res1,meta.method=c("ordmeta"),asymptotic = T)
  MetaDE.Res[["roP_2"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=2,asymptotic = T)
  MetaDE.Res[["roP_4"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=4,asymptotic = T)
  MetaDE.Res[["roP_6"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=6,asymptotic = T)
  MetaDE.Res[["REM"]]<-MetaDE.ES(ind.Res2, meta.method = "REM", tail=ind.tail)
  MetaDE.Res[["FEM"]]<-MetaDE.ES(ind.Res2, meta.method = "FEM", tail=ind.tail)
  MetaDE.Res[["rankProd"]]<-MetaDE.Res.rankProd



  save(MetaDE.Res,file=file.path(res_dir,'Meta_Res_real',paste0('Real_9study',tailed,'_tailed_Meta_Res.RData')))
  save(MetaDE.Res.ind,file=file.path(res_dir,'Meta_Res_real',paste0('Real_9study',tailed,'_tailed_Meta_Res_ind.RData')))
}

# ##################################################################################
# #6 study
study.names<-c("Welsh","Yu","Lapointe","Varambally","Singh","Wallace")
data.QC.raw<-list()
for(i in 1:length(study.names)){
  data.QC.raw[[i]]<-prostate.raw[[study.names[[i]]]]
}

names(data.QC.raw)<-study.names
data.QC.merged<-MetaDE.merge(data.QC.raw)
data.QC.filtered<-MetaDE.filter(data.QC.merged,c(0.2,0.2))#for study.good



MetaDE.Res<-list()
MetaDE.Res.ind<-list()
ndata=length(data.QC.filtered)
data.QC.filtered<-MetaDE:::MetaDE.impute(data.QC.filtered,y=0.3)
weight=rep(1,length(data.QC.filtered))
for (wi in 1:length(weight)) {
  weight[wi] = length(data.QC.filtered[[wi]]$y)
}

MetaDE.Res.rankProd<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",ndata),meta.method=c("rankProd"),asymptotic=F, nperm = 300, metade.perm = F, na.rm = F)
for(ind.tail in c('high','low')){
  tailed=switch (ind.tail,
                 "low" = "left",
                 "high" = "right"
  )


  ndata=length(data.QC.filtered)


  # start<-Sys.time()
  ind.Res1<-ind.analysis(data.QC.filtered,ind.method=rep("modt",ndata),nperm=300,tail=ind.tail)
  ind.Res2<-ind.cal.ES(data.QC.filtered,paired=rep(F,ndata),nperm=300,miss.tol=0.3)

  MetaDE.Res.ind[['ind.modt']]<-ind.Res1
  MetaDE.Res.ind[['ind.ES']]<-ind.Res2
  MetaDE.Res[["Fisher"]]<-MetaDE.pvalue(ind.Res1,meta.method=c("Fisher"),asymptotic = T)
  MetaDE.Res[["Stouffer"]]=MetaDE.pvalue(ind.Res1,meta.method=c("Stouffer"),asymptotic = T)
  MetaDE.Res[["wFisher"]]=MetaDE.pvalue(ind.Res1,meta.method=c("wFisher_sora"),asymptotic = T, weight = weight)
  MetaDE.Res[["ordmeta"]]=MetaDE.pvalue(ind.Res1,meta.method=c("ordmeta"),asymptotic = T)
  MetaDE.Res[["roP_2"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=2,asymptotic = T)
  MetaDE.Res[["roP_4"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=4,asymptotic = T)
  MetaDE.Res[["roP_6"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=6,asymptotic = T)
  MetaDE.Res[["REM"]]<-MetaDE.ES(ind.Res2, meta.method = "REM", tail=ind.tail)
  MetaDE.Res[["FEM"]]<-MetaDE.ES(ind.Res2, meta.method = "FEM", tail=ind.tail)
  MetaDE.Res[["rankProd"]]<-MetaDE.Res.rankProd

  save(MetaDE.Res,file=file.path(res_dir,'Meta_Res_real',paste0('Real_6study',tailed,'_tailed_Meta_Res.RData')))
  save(MetaDE.Res.ind,file=file.path(res_dir,'Meta_Res_real',paste0('Real_6study',tailed,'_tailed_Meta_Res_ind.RData')))
}
# 
# 
# ##################################################################################
# #true gene selection
truegenes<-list()
for(ind.tail in c('high','low')){
  tailed=switch (ind.tail,
                 "low" = "left",
                 "high" = "right"
  )
  load(file.path(res_dir,'Meta_Res_real',paste0('Real_6study',tailed,'_tailed_Meta_Res.RData')))
  FDR_tot<-c()
  for(i in 1:length(MetaDE.Res)){
    if(names(MetaDE.Res)[i]=="rankProd"){
      if(ind.tail=='high'){
        rp_temp<-as.matrix(MetaDE.Res[[i]]$FDR.down)
        colnames(rp_temp)="rankProd"
        FDR_tot%<>%cbind(rp_temp)
      }else if(ind.tail=='low'){
        rp_temp<-as.matrix(MetaDE.Res[[i]]$FDR.up)
        colnames(rp_temp)="rankProd"
        FDR_tot%<>%cbind(rp_temp)
      }else{
        stop('unexpected case')
      }
    }else if(names(MetaDE.Res)[i] %in% c("FEM","REM")){
      FDR_tot%<>%cbind(MetaDE.Res[[i]]$FDR)
    }else{
      res_temp<-MetaDE.Res[[i]]$meta.analysis$FDR
      if(str_detect(names(MetaDE.Res)[i],"roP")){
        colnames(res_temp)=names(MetaDE.Res)[i]
      }
      FDR_tot%<>%cbind(res_temp)
    }
  }
  truede<-c()
  falsede<-c()
  for(i in 1:length(MetaDE.Res)){
    FDR_temp<-FDR_tot
    FDR_temp<-FDR_temp[,-i]
    FDR_temp<-FDR_temp[complete.cases(FDR_temp),]
    max_FDR<-apply(FDR_temp, 1, FUN=function(x){x[order(x, decreasing=T)[3]]})
    min_FDR<-apply(FDR_temp, 1, FUN=function(x){x[order(x)[3]]})

    truede%<>%union(rownames(FDR_temp)[max_FDR<true_gene_threshold])
    falsede%<>%union(rownames(FDR_temp)[min_FDR>true_gene_threshold])
  }
  truegenes[[ind.tail]]$truede<-truede
  truegenes[[ind.tail]]$falsede<-falsede
}
save(truegenes, file=file.path(res_dir,'truegenes_thr_0.01.Rdata'))

