library(MetaDE)
library(metapro)
# Sys.setenv(JAVA_HOME='~/java/jdk1.8.0_261')
# Sys.setenv(CLASSPATH="~/java/jdk1.8.0_261/jre/lib/ext")
library(MetaQC)

args = commandArgs(trailingOnly = T)
if (length(args) < 6)
  stop("Must set 6 indexs", call. = F)
i = as.integer(args[1])
n.diffexp = as.integer(args[2])
total_study = as.integer(args[3])
associated_study = as.integer(args[4])
input_dir=args[5]
res_dir=args[6]
dir.create(file.path(res_dir,'Meta_Res'))



load(file.path(input_dir,'voom',paste0(total_study,'total_studies_',associated_study,'associated-studies_',n.diffexp,'de_',i,'_rep_voom.RData')))
MetaDE.Res<-list()
MetaDE.Res.ind<-list()
ndata=length(simulated_studies)
data.QC.filtered<-simulated_studies
data.QC.filtered<-MetaDE:::MetaDE.impute(data.QC.filtered,y=0.3)
weight=rep(1,length(data.QC.filtered))
for (wi in 1:length(weight)) {
  weight[wi] = length(data.QC.filtered[[wi]]$y)
}
MetaDE.Res.rankProd<-MetaDE.rawdata(data.QC.filtered,ind.method=rep("modt",ndata),meta.method=c("rankProd"),asymptotic=F, nperm = 300, metade.perm = F, na.rm=F)

set.seed(1234)
r_num=sample(1:10, 10)[i%%10]*2
for(ind.tail in c('high', 'low')){
  tailed=switch (ind.tail,
                 "low" = "left",
                 "high" = "right"
  )
  
  ind.Res1<-ind.analysis(data.QC.filtered,ind.method=rep("modt",ndata),nperm=300,tail=ind.tail)
  ind.Res2<-ind.cal.ES(data.QC.filtered,paired=rep(F,ndata),nperm=300,miss.tol=0.3)
  
  MetaDE.Res.ind[['ind.modt']]<-ind.Res1
  MetaDE.Res.ind[['ind.ES']]<-ind.Res2
  MetaDE.Res[["Fisher"]]<-MetaDE.pvalue(ind.Res1,meta.method=c("Fisher"),asymptotic = T)
  MetaDE.Res[["Stouffer"]]=MetaDE.pvalue(ind.Res1,meta.method=c("Stouffer"),asymptotic = T)
  MetaDE.Res[["wFisher"]]=MetaDE.pvalue(ind.Res1,meta.method=c("wFisher_sora"),asymptotic = T)
  MetaDE.Res[["ordmeta"]]=MetaDE.pvalue(ind.Res1,meta.method=c("ordmeta"),asymptotic = T)
  MetaDE.Res[["roP_2"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=2,asymptotic = T)
  MetaDE.Res[["roP_4"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=4,asymptotic = T)
  MetaDE.Res[["roP_6"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=6,asymptotic = T)
  MetaDE.Res[["roP_random"]]=MetaDE.pvalue(ind.Res1,meta.method=c("roP"), rth=r_num,asymptotic = T)
  MetaDE.Res[["REM"]]<-MetaDE.ES(ind.Res2, meta.method = "REM")
  MetaDE.Res[["FEM"]]<-MetaDE.ES(ind.Res2, meta.method = "FEM")
  MetaDE.Res[["rankProd"]]<-MetaDE.Res.rankProd
  
  
  
  
  save(MetaDE.Res,file=file.path(res_dir,'Meta_Res',paste0(total_study,'total_studies_',associated_study,'associated-studies_',n.diffexp,'de_',i,'_rep_',tailed,'_tailed_Meta_Res.RData')))
  save(MetaDE.Res.ind,file=file.path(res_dir,'Meta_Res',paste0(total_study,'total_studies_',associated_study,'associated-studies_',n.diffexp,'de_',i,'_rep_',tailed,'_tailed_Meta_Res_ind.RData')))
}