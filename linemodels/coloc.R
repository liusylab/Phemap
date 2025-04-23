library(tidyverse)
library(data.table)
library(coloc)
args<-commandArgs(T)
pheno<-args[1]

file_list<-list.files(paste0('./coloc/data/mydata/',pheno,'/'),pattern = paste0("^",pheno,"_.*_subdata_1M.txt$"))
coloc_sum1<-data.frame()
for(i in file_list){
  print(i)
  name<-sub('_subdata_1M.txt','',i)
  mydata<-fread(paste0('./coloc/data/mydata/',pheno,'/',i),header = T)
  mydata<-mydata[,c(1:11)] #CHR	BP	SNP	A1	A2	freq	b	se	p	Direction	N
  mydata$p<-ifelse(mydata$p<=1e-300,1e-300,mydata$p)
  EAdata<-fread(paste0('./coloc/data/East_Asian_hg38/',pheno,'/',name,'_subdata_1M_hg38.txt'),header = T)
  colnames(EAdata)[2:3]<-c('CHR','BP')
  EAdata$CHR[EAdata$CHR=='X']<-23
  EAdata$EA_P<-ifelse(EAdata$EA_P<=1e-300,1e-300,EAdata$EA_P)
  EAdata$EA_A1<-toupper(EAdata$EA_A1)
  EAdata$EA_A2<-toupper(EAdata$EA_A2)
  merge_data<-merge(mydata,EAdata,by=c('CHR','BP'))
  for(i in 1:nrow(merge_data)){
    if(merge_data$A1[i]==merge_data$EA_A1[i] & merge_data$A2[i]==merge_data$EA_A2[i]){
      merge_data$align[i]<-'YES'
    }else if(merge_data$A1[i]==merge_data$EA_A2[i] & merge_data$A2[i]==merge_data$EA_A1[i]){
      merge_data$EA_BETA[i]<-(-merge_data$EA_BETA[i])
      merge_data$EA_FRQ[i]<-1-merge_data$EA_FRQ[i]
      merge_data$EA_A1[i]<-merge_data$A1[i]
      merge_data$EA_A2[i]<-merge_data$A2[i]
      merge_data$align[i]<-'YES'
    }else{
      merge_data$align[i]<-'NO'
    }
  }
  merge_data$MAF<-ifelse(merge_data$freq>0.5,1-merge_data$freq,merge_data$freq)
  merge_data$EA_MAF<-ifelse(merge_data$EA_FRQ>0.5,1-merge_data$EA_FRQ,merge_data$EA_FRQ)
  colnames(merge_data)[3]<-'my_SNP'
  merge_data$SNP<-paste0('chr',merge_data$CHR,'_',merge_data$BP,'_',merge_data$A1,'_',merge_data$A2)
  merge_data<-merge_data[,c(24,1:6,22,7:15,23,16:21)]

  n_list<-which(duplicated(merge_data$SNP))
  if(length(n_list)==0){print('no duplicate')}else{
  snp_list<-merge_data[n_list,]$SNP
  a2<-subset(merge_data,SNP%in%snp_list)
  a2<-a2[a2$my_SNP==a2$EA_SNP,]
  if(nrow(a2)!=1){a2<-a2[a2$align=='YES',]}
  merge_data<-subset(merge_data,!(SNP %in% snp_list))
  merge_data<-rbind(merge_data,a2)
  merge_data<-merge_data[order(merge_data$SNP),]
  }

  fwrite(merge_data,paste0('./coloc/data/mergedata/',pheno,'/',name,'_mergedata.txt'),col.names = T,row.names = F,sep = '\t',quote = F)
  
  my.res<-coloc.abf(dataset1 = list(snp=merge_data$SNP,pvalues=merge_data$p,MAF=merge_data$MAF,beta=merge_data$b,se=merge_data$se,N=merge_data$N,type='quant'),
                    dataset2 = list(snp=merge_data$SNP,pvalues=merge_data$EA_P,MAF=merge_data$EA_MAF,beta=merge_data$EA_BETA,se=merge_data$EA_SE,N=merge_data$EA_N,type='quant'),
                    p1=1e-04,p2=1e-04,p12=5e-06)
  coloc_sum<-as.data.frame(t(my.res$summary))
  coloc_sum$pheno<-pheno
  coloc_sum$SNP<-sub(".*_(rs\\d+).*", "\\1", name)
  fwrite(coloc_sum,paste0('./coloc/output/',pheno,'/',name,'_coloc_summary.txt'),col.names = T,row.names = F,sep = '\t',quote = F)
  coloc_sum1<-rbind(coloc_sum1,coloc_sum)
  if(my.res[[1]][6]>=0.8){
    results<-my.res$results[order(my.res$results$SNP.PP.H4,decreasing = T),]
    cs<-cumsum(results$SNP.PP.H4)
    cred_set<-results[1:which(cs>0.95)[1],]
    fwrite(results,paste0('./coloc/output/',pheno,'/',name,'_coloc_results.txt'),col.names = T,row.names = F,sep = '\t',quote = F)
    fwrite(cred_set,paste0('./coloc/output/',pheno,'/',name,'_coloc_95_cred_result.txt'),col.names = T,row.names = F,sep = '\t',quote = F)
  }
}
fwrite(coloc_sum1,paste0('./coloc/output/',pheno,'/',pheno,'_coloc_summary_all.txt'),col.names = T,row.names = F,sep = '\t',quote = F)
