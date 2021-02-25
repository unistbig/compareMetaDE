library(ROCR)
library(reshape2)
library(ggplot2)
library(stringr)
library(magrittr)
library(RColorBrewer)
library(colortools)
args = commandArgs(trailingOnly = T)
if (length(args) < 4)
  stop("Must set 4 indexs", call. = F)
input_dir1=args[1]
input_dir2=args[2]
res_dir=args[3]
true_gene_threshold=as.numeric(args[4])

figure.dir<-file.path(res_dir,'Real_study_analysis_plot')

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

rep.start=1
rep.end=15
# methods_used=c("rankProd", "roP_2", "roP_4", "roP_6", "roP_union", "Stouffer", "REM", "FEM", "Fisher", "wFisher", "ordmeta")
methods_used=c("rankProd", "roP_2", "roP_4", "roP_6", "Stouffer", "REM", "FEM", "Fisher", "wFisher", "ordmeta")
{
  tpr = tfdr = auc = NDE = METHOD = COLOR = REPEAT = nDE_Factor = NDE_datasets =NULL
  
  load(file.path(res_dir,paste0('truegenes_thr_',true_gene_threshold,'.Rdata')))
  for(ind.tail in c('high', 'low')){
    tailed=switch (ind.tail,
                   "low" = "left",
                   "high" = "right"
    )
    
    for(i in rep.start:rep.end)
    {
      
      load(file.path(input_dir1,paste0('Real_9study_' ,i,'_rep_',tailed,'_tailed_Meta_Res.RData')))
      tpr_temp=c()
      tfdr_temp=c()
      auc_temp=c()
      
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
      
      TrueGene<-truegenes[[ind.tail]]$truede
      FalseGene<-truegenes[[ind.tail]]$falsede
      
      
      for(p in colnames(FDR)){
        FDRm<-FDR[,p]
        FDRm<-FDRm[!is.na(FDRm)]
        GeneName = names(FDRm)
        indexTrue = which(GeneName%in%TrueGene)
        indexFalse = which(GeneName%in%FalseGene)
        FDRm<-FDRm[c(indexTrue,indexFalse)]
        indexTrue=1:length(indexTrue)
        indexFalse=(length(indexTrue)+1):(length(indexTrue)+length(indexFalse))
        
        
        
        tpr_temp = append(tpr_temp, length(which(FDRm[indexTrue]<0.05))/length(FDRm[indexTrue]))
        if(length(which(FDRm<0.05))<=5){
          tfdr_temp=append(tfdr_temp, NaN)
        }else{
          tfdr_temp = append(tfdr_temp, length(setdiff(which(FDRm<0.05),indexTrue))/length(which(FDRm<0.05)))
        }
        
        label = rep(0, length(FDRm))
        label[indexTrue] = 1
        pred = prediction(predictions = 1-FDRm, labels = label)
        auc_temp = append(auc_temp, performance(pred, 'auc')@y.values[[1]][1])
        REPEAT = append(REPEAT,i)
        METHOD = append(METHOD, p)
        COLOR = append(COLOR, select_color(p))
      }
      tpr = append(tpr, tpr_temp)
      tfdr = append(tfdr, tfdr_temp)
      auc = append(auc, auc_temp)
      
    }
  }
  
  res = data.frame(Methods = METHOD, Repeat=REPEAT , TPR = tpr, trueFDR = tfdr, AUC = auc, Color = COLOR)
  
  res$Color=factor(res$Color)
  res2 = melt(res, measure.vars=c("AUC","TPR","trueFDR"))
  default_order=c("rankProd", "roP_2", "roP_4", "roP_6", "roP_union", "Stouffer", "REM", "FEM", "Fisher", "wFisher", "ordmeta")
  unique_order<-as.vector(unique(res2$Methods))
  res2$Methods<-factor(res2$Methods, levels=unique(res2$Methods)[match(default_order,unique(res2$Methods))])
  
  
  order_used=default_order[default_order %in% methods_used]
  order_used<-order_used[order_used %in% unique_order]
  axis_order=order_used
  
  for(restype in c('AUC', 'TPR', 'trueFDR')){
    sub.res = res2[res2$variable == restype,]
    pd = position_dodge(width=0.0)
    
    gbase = ggplot(sub.res, aes(y=value, x=Methods, color=Methods))+geom_boxplot(position=pd,outlier.shape=8)+
      scale_x_discrete(limits = axis_order)+
      theme(axis.text.x=element_text(angle=90, hjust=1))+
      scale_colour_manual(name = "Methods",
                          # labels = default_order,
                          values = sapply(default_order,FUN=select_color))
    gline = gbase
    plot(gline)
    
    tt=paste('Real data /',restype ,'/ 7 study permutation / 9 studies')
    tt=gsub(pattern = 'upDE', replacement = 'Bal', x = tt)
    print(gline+aes(x=Methods)+labs(x='', y=restype))
    
    
    
    
    figurename=gsub(pattern = " / ", replacement = "_", x = tt)
    figurename=gsub(pattern=" = ",replacement="_",x=figurename,fixed=T)
    figurename=gsub(pattern="%",replacement="percent",x=figurename,fixed=T)
    figurename=paste(figurename,".pdf",sep = "")
    ggsave(file=paste(figure.dir,"/",figurename,sep=""),width = 10,height = 8)
    dev.off()
  }
}
