library(ROCR)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(colortools)

args = commandArgs(trailingOnly = T)
if (length(args) < 4)
  stop("Must set 4 indexs", call. = F)
rep.end = as.integer(args[1])
#n.diffexp = as.integer(args[2])
total_study = as.integer(args[2])
#associated_study = as.integer(args[4])
input_dir=args[3]
res_dir=args[4]


nvar=1000
n.diffexps=c(100,300,600)
associated_studies=c(2,5,10)

colorset<-brewer.pal(8, "Set1")

select_color = function(method_color)
{
  answer = switch(method_color,
                  "rankProd"=colorset[1],
                  "roP_2"=colorset[2],
                  "roP_4"=adjacent(colorset[2],plot=F)[2],
                  "roP_6"=adjacent(colorset[2],plot=F)[3],
                  "Stouffer"=colorset[3],
                  "REM"=colorset[4],
                  "FEM"=colorset[5],
                  "Fisher"=colorset[6],
                  "wFisher"= colorset[7],
                  "ordmeta"=colorset[8],
                  "roP_union"=sequential(colorset[2])[20]
      
  )
  return(answer)
}



figure.dir<-file.path(res_dir,'Simulation_analysis_plot')
dir.create(figure.dir,showWarnings = F)


# methods_used=c("rankProd", "roP_2", "roP_4", "roP_6", "roP_union", "Stouffer", "REM", "FEM", "Fisher", "wFisher", "ordmeta")
methods_used=c("rankProd", "roP_union", "Stouffer", "REM", "FEM", "Fisher", "wFisher", "ordmeta")
{
  tpr = tfdr = auc = NDE = METHOD = COLOR = REPEAT = nDE_Factor = NDE_datasets =NULL
  for(n.diffexp in n.diffexps){
    for(associated_study in associated_studies){
      for(ind.tail in c('high', 'low')){
        tailed=switch (ind.tail,
                       "low" = "left",
                       "high" = "right"
        )
        for(i in 1:rep.end)
        {
          tpr_temp=c()
          tfdr_temp=c()
          auc_temp=c()
          load(file.path(res_dir,'Meta_Res',paste0(total_study,'total_studies_',associated_study,'associated-studies_',n.diffexp,'de_',i,'_rep_',tailed,'_tailed_Meta_Res.RData')))
          
          FDR_tot<-c()
          for(k in 1:length(MetaDE.Res)){
            if(names(MetaDE.Res)[k]=="rankProd"){
              if(ind.tail=='high'){
                rp_temp<-as.matrix(MetaDE.Res[[k]]$FDR.down)
                colnames(rp_temp)="rankProd"
                FDR_tot%<>%cbind(rp_temp)
              }else if(ind.tail=='low'){
                rp_temp<-as.matrix(MetaDE.Res[[k]]$FDR.up)
                colnames(rp_temp)="rankProd"
                FDR_tot%<>%cbind(rp_temp)
              }else{
                stop('unexpected case')
              }
            }else if(names(MetaDE.Res)[k] %in% c("FEM","REM")){
              FDR_tot%<>%cbind(MetaDE.Res[[k]]$FDR)
            }else{
              res_temp<-MetaDE.Res[[k]]$meta.analysis$FDR
              if(str_detect(names(MetaDE.Res)[k],"roP")){
                colnames(res_temp)=names(MetaDE.Res)[k]
              }
              FDR_tot%<>%cbind(res_temp)
            }
          }
          
          col_FDR<-FDR_tot%>%colnames()
          col_FDR%<>%gsub(pattern='wFisher_sora', replacement = 'wFisher')
          col_FDR%<>%gsub(pattern='ordmeta_sora', replacement = 'ordmeta')
          col_FDR%<>%gsub(pattern='roP_random', replacement = 'roP_union')
          colnames(FDR_tot)<-col_FDR
          FDR<-FDR_tot
          
          
          if(n.diffexp > 0){
            if(n.diffexp < nrow(FDR)){
              if(tailed=='left'){
                TrueGene = paste('gene',1:(n.diffexp/2),sep="")
                FalseGene = paste('gene', ((n.diffexp/2)+1):(nrow(FDR)),sep="")
              }else if(tailed=='right'){
                TrueGene = paste('gene',((n.diffexp/2)+1):n.diffexp,sep="")
                FalseGene = paste('gene', c(1:(n.diffexp/2),(n.diffexp+1):(nrow(FDR))),sep="")
              }else{
                stop('tail should be left or right')
              }
            }else if(n.diffexp == nrow(FDR)){
              TrueGene = paste('gene',1:nrow(FDR),sep="")
              FalseGene = ''
            }else{
              stop("n.diffexp cannot exceed number of total genes")
            }
          }
          
          
          for(p in colnames(FDR)){
            FDRm<-FDR[,p]
            FDRm<-FDRm[!is.na(FDRm)]
            GeneName = names(FDRm)
            indexTrue = which(GeneName%in%TrueGene)
            indexFalse = which(GeneName%in%FalseGene)
            
            tpr_temp = append(tpr_temp, length(which(FDRm[indexTrue]<0.05))/length(FDRm[indexTrue]))
            if(length(which(FDRm<0.05))<=5){
              tfdr_temp=append(tfdr_temp, NaN)
            }else{
              tfdr_temp = append(tfdr_temp, length(setdiff(which(FDRm<0.05),indexTrue))/length(which(FDRm<0.05)))
            }
            
            # tfdr_temp = append(tfdr_temp, length(setdiff(which(FDRm<0.05),indexTrue))/length(which(FDRm<0.05)))
            
            label = rep(0, length(FDRm))
            label[indexTrue] = 1
            pred = prediction(predictions = 1-FDRm, labels = label)
            auc_temp = append(auc_temp, performance(pred, 'auc')@y.values[[1]][1])
            REPEAT = append(REPEAT,i)
            NDE = append(NDE, paste("pDE = ",round(n.diffexp*100/nvar,2),"%",sep=""))
            #UPPROP = append(UPPROP, paste("upDE = ",round(prop*100,2),"%",sep=""))
            METHOD = append(METHOD, p)
            COLOR = append(COLOR, select_color(p))
            NDE_datasets = append(NDE_datasets, associated_study)
          }
          tpr = append(tpr, tpr_temp)
          tfdr = append(tfdr, tfdr_temp)
          auc = append(auc, auc_temp)
          
        }
      }
    }
  }
  res = data.frame(Methods = METHOD, Repeat=REPEAT , nDE=NDE, nDEds = NDE_datasets, TPR = tpr, trueFDR = tfdr, AUC = auc, Color = COLOR)
  #res[which(res$nDE=="pDE = 60%"),]
  res$Color=factor(res$Color)
  res$nDEds=paste(res$nDEds, "associated studies", sep=' ')
  res$nDEds=factor(res$nDEds, levels=unique(res$nDEds))
  
  res2 = melt(res, measure.vars=c("AUC","TPR","trueFDR"))
  
  
  
  default_order=c("rankProd", "roP_2", "roP_4", "roP_6", "roP_union", "Stouffer", "REM", "FEM", "Fisher", "wFisher", "ordmeta")
  unique_order<-as.vector(unique(res2$Methods))
  res2$Methods<-factor(res2$Methods, levels=unique(res2$Methods)[match(default_order,unique(res2$Methods))])
  #default_order=as.vector(unique(res2$Methods))
  
  
  order_used=default_order[default_order %in% methods_used]
  order_used<-order_used[order_used %in% unique_order]
  axis_order=order_used
  
  # sub.res = res2[res2$variable == "AUC",]
  # sub.res = res2[res2$variable == "TPR",]
  # sub.res = res2[res2$variable == "trueFDR",]
  
  for(restype in c('AUC', 'TPR', 'trueFDR')){
    sub.res = res2[res2$variable == restype,]
    pd = position_dodge(width=0.0)
    
    gbase = ggplot(sub.res, aes(y=value, x=Methods, color=Methods))+geom_boxplot(position=pd,outlier.shape=8)+
      facet_grid(nDE~nDEds, scales='free')+
      scale_x_discrete(limits = axis_order)+
      theme(axis.text.x=element_text(angle=90, hjust=1))+
      scale_colour_manual(name = "Methods",
                          # labels = default_order,
                          values = sapply(default_order,FUN=select_color))
    
    gline = gbase
    
    tt=paste('Simulation /',restype ,'/', nvar ,'genes / 50% upDE / 20 datasets')
    tt=gsub(pattern = 'upDE', replacement = 'Bal', x = tt)
    # print(gline+aes(x=Methods)+labs(x='', y=restype)+ggtitle(tt))
    print(gline+aes(x=Methods)+labs(x='', y=restype))
    
    
    
    
    figurename=gsub(pattern = " / ", replacement = "_", x = tt)
    figurename=gsub(pattern=" = ",replacement="_",x=figurename,fixed=T)
    figurename=gsub(pattern="%",replacement="percent",x=figurename,fixed=T)
    figurename=paste(figurename,".pdf",sep = "")
    ggsave(file=file.path(figure.dir,figurename),width = 10,height = 8)
    
    dev.off()
  }
}
