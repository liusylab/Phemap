#This is a R script with 122 gestational phenotypes as exposures and Abortion in the Biobank Japan as the outcome.

library(TwoSampleMR)
library(data.table)
library(R.utils)
library(readxl)
pheno='Abortion'

expos_path<-'./meta_analysis/phemap_220108/GCTA_input_file_230313/'
outcome_path<-'./MR/BBJ_summary/BBJ/'
results_path<-'./MR/phemap_BBJ/results_20240317/'
IVfile<-list.files('./meta_analysis/phemap_220108/GCTA_file_230314/',pattern = '_IV_clean.txt$')
phenotype<-c('maternal_age')
for(expos in phenotype){
  name=expos
  cc <- paste0('./meta_analysis/phemap_220108/GCTA_file_230314/',name,'_IV_clean.txt')   
  IV_exposure<-read.table(cc,header=T,sep='\t',stringsAsFactors = F)
  gwasresults_exposure<-fread(paste0(expos_path,name,'.meta.ma'),header=T,sep='\t',stringsAsFactors = F)
  gwasresults_exposure<-subset(gwasresults_exposure,SNP %in% IV_exposure$SNP)
  gwasresults_exposure$p<-as.numeric(gwasresults_exposure$p)
  gwasresults_exposure$phenotype<-name
  exp_dat <- format_data(
    gwasresults_exposure,
    phenotype_col = 'phenotype',
    type='exposure',
    snp_col = 'SNP',
    beta_col = 'beta',
    se_col = 'se',
    effect_allele_col ='effect',
    other_allele_col = 'other',
    eaf_col = 'eaf',
    pval_col = 'p')
  outname<-pheno
    if(name==outname){next}else{
      outcc <- paste0(outcome_path,outname,'.gwas.summary')
      gwasresults_outcome<-fread(outcc, header=TRUE, stringsAsFactors=FALSE, sep='\t')
      gwasresults_outcome$phenotype<-outname
      gwasresults_outcome<-subset(gwasresults_outcome,p!='NA')
      snps<-subset(gwasresults_outcome,gwasresults_outcome$SNP %in% exp_dat$SNP)
      if(length(snps$SNP)==0){next}else{
        out_dat <- format_data(
          dat=gwasresults_outcome,
          type = 'outcome',
          snps = exp_dat$SNP,
          header = TRUE,
          phenotype_col = 'phenotype',
          snp_col = 'SNP',
          beta_col = 'beta',
          se_col = 'se',
          effect_allele_col = 'A1',
          other_allele_col = 'A2',
          eaf_col = 'eaf',
          pval_col = 'p')
        
        #harmonise_data
        mydata <- harmonise_data(
          exposure_dat=exp_dat,
          outcome_dat=out_dat,
          action= 1)
        write.table(mydata,file = paste0(results_path,outname,'/mydata_',outname,'_',name,'.txt'),quote = F,sep = '\t')
        mydata_steiger<-subset(mydata,is.na(beta.outcome)==F)
	steiger<-mr_steiger(p_exp=mydata_steiger$pval.exposure,p_out=mydata_steiger$pval.outcome,gwasresults_exposure$N[1],gwasresults_outcome$N[1],r_exp = mydata_steiger$beta.exposure,r_out = mydata_steiger$beta.outcome)
        s<-data.frame(steiger[1:12])
        s$exposure<-name
        s$outcome<-outname
        write.table(s,file=paste0(results_path,outname,'/steiger_',outname,'.csv'),quote = F,append = T,sep=',',row.names = F)
        #run mr
        res <- mr(mydata)
        if(length(res)==0){next}else{
          write.table(res,file=paste0(results_path,outname,'/res_GCTA_',outname,'.csv'),quote = F,append = T,sep=',',row.names = F)
          if(res$nsnp[1]>1){
            het <- mr_heterogeneity(mydata)
            colnames(het)[5]<-'method_het'
            if(length(het$Q_pval)==2&het$Q_pval[1]>het$Q_pval[2]){het<-het[2,]}else{het<-het[1,]}
            res1<-merge(res,het,all.x = T,by=c('id.exposure','id.outcome','outcome','exposure'))
            pleio <- mr_pleiotropy_test(mydata)
            colnames(pleio)[6:7]<-c('pleio_se','pleio_p')
            res1<-merge(res1,pleio,all.x = T,by=c('id.exposure','id.outcome','outcome','exposure'))
            write.table(res1,file=paste0(results_path,outname,'/res_GCTA_',outname,'_1','.csv'),quote = F,append = T,sep = ',',row.names = F)
            pdf(file = paste0(results_path,outname,'/het_pleio_',outname,'_',name,'.pdf'))
            single <- mr_leaveoneout(mydata)
            print(mr_leaveoneout_plot(single))
            print(mr_scatter_plot(res,mydata))
            res_single <- mr_singlesnp(mydata)
            print(mr_forest_plot(res_single))
            print(mr_funnel_plot(res_single))
            dev.off()
          }else{write.table(res,file=paste0(results_path,outname,'/res_GCTA_',outname,'_1','.csv'),quote = F,append = T,sep = ',',row.names = F)}}}}}