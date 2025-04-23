library(data.table)
library(dplyr)
library(readxl)
library(MRPRESSO)
library(TwoSampleMR)
mr_dir<-'./MR/phemap_BBJ/results'
results_dir<-'./MR/phemap_BBJ/MRPRESSO_results'
pheno<-'Abortion'
Twosamplemr_results<-read.csv(paste0(mr_dir,'/',pheno,'/res_GCTA_',pheno,'_1.csv'),sep=',')
Twosamplemr_results<-subset(Twosamplemr_results,exposure!='exposure')
print(head(Twosamplemr_results))
heterogeneity_results<-Twosamplemr_results[!duplicated(Twosamplemr_results$exposure),]
    for(expos_name in heterogeneity_results$exposure){
  mydata_raw<-read.table(paste0(mr_dir,'/',pheno,'/mydata_',pheno,'_',expos_name,'.txt'),header = T,sep='\t')
  if(length(mydata_raw$SNP)<4){next}else{
  PRESSO<-mr_presso(BetaOutcome ='beta.outcome', BetaExposure = 'beta.exposure', SdOutcome ='se.outcome', SdExposure = 'se.exposure', 
                     OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = mydata_raw, NbDistribution = 5000,  
                     SignifThreshold = 0.05)
  presso_result<-PRESSO$`Main MR results`
  presso_result$Exposure<-expos_name
  presso_result$Outcome<-pheno
  presso_result$Outliers_snp<-length(PRESSO$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
  presso_result$Distortion_p<-PRESSO$`MR-PRESSO results`$`Distortion Test`$Pvalue
  write.table(presso_result,paste0(results_dir,'/',pheno,'_PRESSO/',pheno,'_as_outcome_presso.csv'),quote = F,append = T,sep=',',row.names = F)
  #remove Outliers TwosampleMR analysis
  Outliers_index<-PRESSO$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  print(Outliers_index)
  if(length(Outliers_index)== 0){next}else{if(Outliers_index=='No significant outliers'){next}else{
  mydata_rm<-mydata_raw[-Outliers_index,]
  write.table(mydata_rm,file = paste0(results_dir,'/',pheno,'_PRESSO/mydata_',expos_name,'_',pheno,'.txt'),quote = F,sep = '\t')
  res <- mr(mydata_rm)
  if(length(res)==0){next}else{
    write.table(res,file=paste0(results_dir,'/',pheno,'_PRESSO/',pheno,'_as_outcome_TwosampleMR.csv'),quote = F,append = T,sep=',',row.names = F)
    if(res$nsnp[1]>1){
      het <- mr_heterogeneity(mydata_rm)
      colnames(het)[5]<-'method_het'
      if(length(het$Q_pval)==2&het$Q_pval[1]>het$Q_pval[2]){het<-het[2,]}else{het<-het[1,]}
      res1<-merge(res,het,all.x = T,by=c('id.exposure','id.outcome','outcome','exposure'))
      pleio <- mr_pleiotropy_test(mydata_rm)
      colnames(pleio)[6:7]<-c('pleio_se','pleio_p')
      res1<-merge(res1,pleio,all.x = T,by=c('id.exposure','id.outcome','outcome','exposure'))
      write.table(res1,file=paste0(results_dir,'/',pheno,'_PRESSO/',pheno,'_as_outcome_TwosampleMR_1.csv'),quote = F,append = T,sep = ',',row.names = F)
      pdf(file = paste0(results_dir,'/',pheno,'_PRESSO/','/het_pleio_',expos_name,'_',pheno,'.pdf'))
      single <- mr_leaveoneout(mydata_rm)
      print(mr_leaveoneout_plot(single))
      print(mr_scatter_plot(res,mydata_rm))
      res_single <- mr_singlesnp(mydata_rm)
      print(mr_forest_plot(res_single))
      print(mr_funnel_plot(res_single))
      dev.off()
    }else{write.table(res,file=paste0(results_dir,'/',pheno,'_PRESSO/',pheno,'_as_outcome_TwosampleMR_1.csv'),quote = F,append = T,sep = ',',row.names = F)
      }}}}}}