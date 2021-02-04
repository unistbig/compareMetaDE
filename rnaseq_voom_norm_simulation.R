library(limma)
library(SimSeq)
library(tidyverse)
library(Biobase)
library(magrittr)


rtnorm <- function(n, mean, sd, a = -Inf, b = Inf){
  qnorm(runif(n, pnorm(a, mean, sd), pnorm(b, mean, sd)), mean, sd)
}


getDisp = function(mean, mean.condition, disp.condition)
{
  pool = disp.condition[which(mean.condition>(mean-20) & mean.condition<(mean+20))]
  if(length(pool)==0){value = disp.condition[which.min(abs(mean.condition-mean))]
  }else{value = sample(pool,1)}
  return(value)
}

getDisp2 = function(mean, mean.condition, disp.condition)
{
  pool = disp.condition[which(mean.condition>(mean-20) & mean.condition<(mean+20))]
  if(length(pool)==0){
    pool = disp.condition[which(abs(mean.condition-mean)<sort(abs(mean.condition-mean))[20])]
    value=sample(pool,1)
  }else{value = sample(pool,1)}
  return(value)
}


generateDatasetParameter = function(){
  #  Effect size 1.2/1.5~
  data(kidney, package='SimSeq') # Load TCGA KIRC RNA-seq data. 144 samples (72 cancer and matched normal, respectively)
  k_count = kidney$counts # RNA-seq read count data
  index.cancer=(1:72)*2 # cancer sample index
  index.normal=index.cancer-1 # normal sample index
  
  k_count= k_count[,c(index.cancer, index.normal)] # Arrange samples for convinience
  
  # Get count mean and dispersion using edgeR package
  
  # Mean and dispersion values are obtained separately from cancer and normal samples when different dispersion is assummed between two sample types.
  # Mean and dispersion values from normal samples
  dge.normal=DGEList(counts=k_count[,73:144], group = factor(rep(2,72)))
  dge.normal=calcNormFactors(dge.normal)
  dge.normal=estimateCommonDisp(dge.normal)
  dge.normal=estimateTagwiseDisp(dge.normal)
  disp.normal = dge.normal$tagwise.dispersion # Dispersion
  mean.normal = apply(k_count[,73:144],1,mean)
  
  
  # Mean and dispersion values from cancer samples
  dge.cancer=DGEList(counts=k_count[,1:72], group=factor(rep(1,72)))
  dge.cancer=calcNormFactors(dge.cancer)
  dge.cancer=estimateCommonDisp(dge.cancer)
  dge.cancer=estimateTagwiseDisp(dge.cancer)
  disp.cancer = dge.cancer$tagwise.dispersion
  mean.cancer = apply(k_count[,1:72],1,mean)
  
  
  # Gene filtering: Genes having small read count (<10) are filtered
  k_mean.total = apply(k_count,1,mean)
  k_index.filter = which(k_mean.total < 10)
  k_mean.total = k_mean.total[-k_index.filter]
  disp.normal = disp.normal[-k_index.filter]
  disp.cancer = disp.cancer[-k_index.filter]
  mean.normal = mean.normal[-k_index.filter]
  mean.cancer = mean.cancer[-k_index.filter]
  
  
  # Mean and dispersion values obtained using all samples when same dispersion is assummed between two sample types..
  k_dge.total = DGEList(counts = k_count, group = factor(c(rep(1,72),rep(2,72))))
  k_dge.total = calcNormFactors(k_dge.total)
  k_dge.total = estimateCommonDisp(k_dge.total)
  k_dge.total = estimateTagwiseDisp(k_dge.total)
  k_disp.total = k_dge.total$tagwise.dispersion
  k_disp.total = k_disp.total[-k_index.filter]
  
  
  dataset.parameters=list(k_count=k_count, disp.normal=disp.normal, mean.normal=mean.normal, disp.cancer=disp.cancer, mean.cancer=mean.cancer, k_mean.total=k_mean.total, k_index.filter=k_index.filter, k_disp.total=k_disp.total)
  
  return(dataset.parameters)
}

library(limma)
library(edgeR)
dataset.parameters<-generateDatasetParameter()

dataset='hello_wolrd'
simul.data='KIRC'

data_gen=function(dataset.parameters, fraction.upregulated=0.5, dat_num=20, diff.dat_num=2, n.var=1000, n.diffexp=600, random_sampling=FALSE){
  
  n=dat_num
  
  s<-sample(c(3:20),size=n*1000,replace=T)
  s_ratio<-2^abs(rtnorm(n*1000, mean=0, sd=1, a=-log2(3), b=log2(3)))
  
  
  s_1=s
  s_2=s*s_ratio
  
  s_1_temp<-(s-min(s))*((20-3)/(max(s*s_ratio)-min(s)))
  s_2_temp<-(s*s_ratio-min(s))*((20-3)/(max(s*s_ratio)-min(s)))
  
  
  adjust_ind<-which(s*s_ratio>20)
  s_1[adjust_ind]=(s_1_temp[adjust_ind]+3)
  s_2[adjust_ind]=(s_2_temp[adjust_ind]+3)
  s_1 %<>% round()
  s_2 %<>% round()
  
  s_1_temp=c()
  s_2_temp=c()
  
  
  sampled<-sample(n*1000,size=n)
  s_1_temp<-s_1[sampled]
  s_2_temp<-s_2[sampled]
  
  #hist(rtnorm(n*1000, mean=5.5, sd=1, a=0, b=11))
  sampler<-floor(rtnorm(1, mean=n/2, sd=1, a=0, b=(n+1)))
  sampler[sampler==(n+1)]=n
  sample_ratio<-sample(n,sampler)
  s_con_use=c(s_1_temp[sample_ratio],s_2_temp[-sample_ratio])
  s_tum_use=c(s_2_temp[sample_ratio],s_1_temp[-sample_ratio])
  
  s1=s_con_use
  s2=s_tum_use
  k_count=dataset.parameters$k_count
  disp.normal=dataset.parameters$disp.normal
  disp.cancer=dataset.parameters$disp.cancer
  mean.normal=dataset.parameters$mean.normal
  mean.cancer=dataset.parameters$mean.cancer
  k_mean.total=dataset.parameters$k_mean.total
  k_index.filter=dataset.parameters$k_index.filter
  k_disp.total=dataset.parameters$k_disp.total
  
  
  
  random.index = sample(1:length(k_mean.total), size=n.var)
  sample.mean1 = k_mean.total[random.index]
  
  sample.mean2 = sample.mean1
  
  if(n.diffexp!=0){
    upindex = 1:round(n.diffexp*fraction.upregulated)
    dnindex = round(n.diffexp*fraction.upregulated+1):n.diffexp
    

    factor1 = 1.3+rexp(n = length(upindex), rate=1) #fold change 1.3~
    factor2 = 1.3+rexp(n = length(dnindex), rate=1)
    
    sample.mean2[upindex] = sample.mean2[upindex]*factor1
    sample.mean2[dnindex] = sample.mean2[dnindex]/factor2
  }
  
  mean.condition1=mean.normal
  mean.condition2=mean.cancer
  mean.total=k_mean.total
  disp.condition1=disp.normal
  disp.condition2=disp.cancer
  disp.total=k_disp.total
  
  sample.disp1 = disp.total[random.index]
  sample.disp2 = sample.disp1
  
  
  random_diffexp=c()
  for(j in seq_len(n.var)){
    random_diffexp%<>%rbind(sample(dat_num,diff.dat_num))
  }
  
  dat_tot=list()
  L=3^runif(dat_num, min=-1, max=1)
  for(dat_ind in seq_len(dat_num)){
    counts = matrix(nrow=n.var, ncol = s1[dat_ind]+s2[dat_ind])
    if(random_sampling==TRUE){
      rand1=runif(s1[dat_ind],min=0.7,max=1.3)
      rand2=runif(s2[dat_ind],min=0.7,max=1.3)
    }else{
      rand1=rep(1,s1[dat_ind])
      rand2=rep(1,s2[dat_ind])
    }
    #####################################################################
    L_temp=L[dat_ind]
    sample.mean1_temp = sample.mean1*L
    sample.mean2_temp = sample.mean2*L
    sample.disp1_temp = sapply(sample.mean1_temp, FUN = getDisp2, mean.condition = mean.total, disp.condition = disp.total, simplify = T, USE.NAMES = F)
    sample.disp2_temp = sapply(sample.mean2_temp, FUN = getDisp2, mean.condition = mean.total, disp.condition = disp.total, simplify = T, USE.NAMES = F)
    #####################################################################
    
    for(i in 1:n.var)
    {
      if(dat_ind %in% random_diffexp[i,]){
        counts[i,1:s1[dat_ind]] = sapply(rand1, FUN = function(x) rnbinom(1, 1/sample.disp1_temp[i], mu=sample.mean1_temp[i]*x))
        counts[i,(s1[dat_ind]+1):(s1[dat_ind]+s2[dat_ind])] = sapply(rand2, FUN = function(x) rnbinom(1, 1/sample.disp2_temp[i], mu=sample.mean2_temp[i]*x))
      }else{
        counts[i,1:s1[dat_ind]] = sapply(rand1, FUN = function(x) rnbinom(1, 1/sample.disp1_temp[i], mu=sample.mean1_temp[i]*x))
        counts[i,(s1[dat_ind]+1):(s1[dat_ind]+s2[dat_ind])] = sapply(rand2, FUN = function(x) rnbinom(1, 1/sample.disp1_temp[i], mu=sample.mean1_temp[i]*x))
      }
    }
    
    
    rownames(counts)=paste0('gene',seq_len(nrow(counts)))
    colnames(counts)=c(paste0('Normal.',seq_len(s1[dat_ind])),paste0('Cancer.',seq_len(s2[dat_ind])))
    dat_tot[[paste0('Dataset_',dat_ind)]][['x']]=counts
    dat_tot[[paste0('Dataset_',dat_ind)]][['y']]=c(rep(0,s1[dat_ind]),rep(1,s2[dat_ind]))
    dat_tot[[paste0('Dataset_',dat_ind)]][['k']]=L_temp
  }
  return(dat_tot)
}

args = commandArgs(trailingOnly = T)
if (length(args) < 3)
  stop("Must set 3 indexs", call. = F)
i = as.integer(args[1])
n.diffexp = as.integer(args[2])
diff.dat_num = as.integer(args[3])
res_dir=args[4]

dir.create(file.path(res_dir,'count'),showWarnings = F)
simulated_studies<-data_gen(dataset.parameters, n.diffexp = n.diffexp, diff.dat_num = diff.dat_num)
save(simulated_studies, file=file.path(res_dir,'count',paste0(diff.dat_num,'diffds_',n.diffexp,'de_',i,'_rep_count.RData')))

for(k in 1:length(simulated_studies)){
  nf <- edgeR::calcNormFactors(simulated_studies[[k]]$x, method = 'TMM')
  voom.data <- limma::voom(simulated_studies[[k]]$x, design = model.matrix(~factor(simulated_studies[[k]]$y)), lib.size = colSums(simulated_studies[[k]]$x) * nf)
  voom.data$genes <- rownames(simulated_studies[[k]]$x)
  simulated_studies[[k]]$count<-simulated_studies[[k]]$x
  simulated_studies[[k]]$x<-voom.data$E
}
dir.create(file.path(res_dir,'voom'),showWarnings = F)
save(simulated_studies, file=file.path(res_dir,'voom',paste0(diff.dat_num,'diffds_',n.diffexp,'de_',i,'_rep_voom.RData')))

