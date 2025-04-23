#Due to unstable internet connectivity, we directly wrote out the functions used and modified them according to our specific needs. 
#The reference website is: https://marinalearning.netlify.app/2021/03/22/setting-up-multivariable-mendelian-randomization-analysis/

library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(dplyr)
library(TwoSampleMR)
library(MVMR)
library(data.table)
library(plinkbinr)
library(ieugwasr)
#Defining function
ld_clump_local_YT=function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin) {
  #debug:
  # dat= data.frame(rsid=exposure$SNP, 
  #                     pval=exposure$pval.exposure, 
  #                     id=exposure$id.exposure)
  # .....................
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", 
                    "sh")
 
    fn <- tempfile(tmpdir='./MVMR/bin/tempfile')
    
    write.table(data.frame(SNP = dat[["rsid"]], P = dat[["pval"]]), 
                file = fn, row.names = F, col.names = T, quote = F)

    fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ",
                   shQuote(bfile, type = shell), " --clump ",
                   shQuote(fn,type = shell), " --clump-p1 ", clump_p, " --clump-r2 ",
                   clump_r2, " --clump-kb ", clump_kb, " --out ",
                   shQuote(fn,type = shell))

    system(fun2)
    
    #if ######.clumped exists
    if(file.exists(paste(fn, ".clumped", sep = ""))){
      res <- read.table(paste(fn, ".clumped", sep = ""), header = T)
      y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
      if (nrow(y) > 0) {
        message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), 
                " variants due to LD with other variants or absence from LD reference panel")
      }
      unlink(paste(fn, "*", sep = ""))
      return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
    }else{ #does not exists clumped data
      return(dat)
    }
  }

ld_clump_YT=function (dat = NULL, clump_kb = 1000, clump_r2 = 0.1, clump_p = 0.99, 
            pop = "EAS", access_token = NULL, bfile = NULL, plink_bin = NULL){
    stopifnot("rsid" %in% names(dat))
    stopifnot(is.data.frame(dat))
    if (is.null(bfile)) {
      message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
    }
    if (!"pval" %in% names(dat)) {
      if ("p" %in% names(dat)) {
        warning("No 'pval' column found in dat object. Using 'p' column.")
        dat[["pval"]] <- dat[["p"]]
      }
      else {
        warning("No 'pval' column found in dat object. Setting p-values for all SNPs to clump_p parameter.")
        dat[["pval"]] <- clump_p
      }
    }
    if (!"id" %in% names(dat)) {
      dat$id <- random_string(1)
    }
    if (is.null(bfile)) {
      access_token = check_access_token()
    }
    ids <- unique(dat[["id"]])
    res <- list()
    for (i in 1:length(ids)) {
      x <- subset(dat, dat[["id"]] == ids[i])
      if (nrow(x) == 1) {
        message("Only one SNP for ", ids[i])
        res[[i]] <- x
      }
      else {
        message("Clumping ", ids[i], ", ", nrow(x), " variants, using ", 
                pop, " population reference")
        if (is.null(bfile)) {
          res[[i]] <- ld_clump_api(x, clump_kb = clump_kb, 
                                   clump_r2 = clump_r2, clump_p = clump_p, pop = pop, 
                                   access_token = access_token)
        }
        else {
          res[[i]] <- ld_clump_local_YT(x, clump_kb = clump_kb, 
                                     clump_r2 = clump_r2, clump_p = clump_p, bfile = bfile, 
                                     plink_bin = plink_bin)
        }
      }
    }
    res <- dplyr::bind_rows(res)
    return(res)
  }

get_mv_exposures<-function(tophits_list,full_gwas_list){
  exposures<-bind_rows(tophits_list)
  
  temp<-exposures
  temp$id.exposure<-1
  temp<-ld_clump_YT(tibble(rsid=temp$SNP,pval=temp$pval.exposure,id=temp$id.exposure),
  plink_bin='./plink_linux_x86_64_20231211/plink',
  bfile='./MR/MVMR/bin/tempfile/chr1-22.CGP.beagle52.filter_MAF.0.001.addrsid.filt',
  clump_kb=1000,clump_r2=0.1)
  colnames(temp)[1]<-'SNP'
  exposures<-filter(exposures,SNP %in% temp$SNP)
  
  for(i in 1:length(full_gwas_list)){
    full_gwas_list[[i]]<-full_gwas_list[[i]] %>% filter(SNP %in% exposures$SNP)
  }
  
  d1<-bind_rows(full_gwas_list)%>%distinct()
  id_exposure<-unique(d1$id.outcome)
  tmp_exposure<-d1%>%filter(id.outcome==id_exposure[1])%>%convert_outcome_to_exposure()
  tmp_outcome<-d1%>%filter(id.outcome!= id_exposure[1])
  
  d<-harmonise_data(exposure_dat = tmp_exposure,outcome_dat = tmp_outcome,action = 2)
  snps_not_in_all<-d%>%
    count(SNP)%>%
    filter(n<length(tophits_list)-1)%>% #why?
    pull(SNP)
  d<-filter(d,!SNP%in%snps_not_in_all) 
  
  dh1x<-d%>%filter(id.outcome==id.outcome[1])%>%
    select(SNP,contains('exposure'))
  dh2x<-d%>%select(SNP,contains('outcome'))
  names(dh2x)<-gsub('outcome','exposure',names(dh2x))
  exposure_dat<-bind_rows(dh1x,dh2x)%>%
    select(-c('samplesize.exposure','mr_keep.exposure','pval_origin.exposure'))%>%
    distinct()
  
  return(exposure_dat)
}

make_mvmr_input<-function(exposure_dat,outcome.id.mrbase=NULL,outcome.data=NULL){
  if(!is.null(outcome.id.mrbase)){
    outcome_dat<-extract_outcome_data(snps =unique(exposure_dat$SNP),
                                      outcomes = outcome.id.mrbase)
  }else if(!is.null(outcome.data)){
    outcome_dat<-outcome.data%>%filter(SNP %in% exposure_dat$SNP)
  }
  exposure_dat<-exposure_dat%>%mutate(id.exposure=exposure)
  outcome_harmonised<-mv_harmonise_data(exposure_dat,outcome_dat)
  exposures_order<-colnames(outcome_harmonised$exposure_beta)
  
  no_exp<-dim(outcome_harmonised$exposure_beta)[2] #count exposures
  colnames(outcome_harmonised$exposure_beta)<-paste0('betaX',1:no_exp)
  colnames(outcome_harmonised$exposure_se)<-paste0('seX',1:no_exp)
  
  XGs<-left_join(as.data.frame(outcome_harmonised$exposure_beta) %>% rownames_to_column('SNP'),
                 as.data.frame(outcome_harmonised$exposure_se) %>% rownames_to_column('SNP'),
                 by='SNP')
  
  YG<-data.frame(beta.outcome=outcome_harmonised$outcome_beta,
                 se.outcome=outcome_harmonised$outcome_se) %>%
    mutate(SNP=XGs$SNP)
  
  return(list(YG=YG,XGs=XGs,exposures=exposures_order))
}

tidy_mvmr_output<-function(mvmr_res){
  mvmr_res%>%
    as.data.frame()%>%
    rownames_to_column('exposure')%>%
    rename(b=Estimate,se='Std. Error',pval='Pr(>|t|)')%>%
    TwoSampleMR::generate_odds_ratios()
}

#start

args<-commandArgs(T)
pheno<-args[1]
  print(pheno)
  #1st:tidy up exposure GWAS data.
    gwas_exposure_format1<-vroom(paste0('./Taiwan/data/GCTA_input_data/',pheno,'.ma'),
                                 col_select = c('SNP','beta','se','effect','other','eaf','p'))%>%
      format_data(.,type='exposure',
                  snp_col='SNP',
                  beta_col='beta',
                  se_col='se',
                  effect_allele_col='effect',
                  other_allele_col='other',
                  eaf_col='eaf',
                  pval_col='p')%>%
      mutate(exposure=paste0(pheno,'_Taiwan'))
  gwas_exposure_format2<-vroom(paste0('./meta_GWAS_summary_statistics_230313/',pheno,'.meta.format.gz'),
                              col_select = c('SNP','b','se','A1','A2','freq','p') )%>%
    format_data(.,type='exposure',
                snp_col='SNP',
                beta_col='b',
                se_col='se',
                effect_allele_col='A1',
                other_allele_col='A2',
                eaf_col='freq',
                pval_col='p')%>%
    mutate(exposure=paste0(pheno,'_nipt'))
  
  #2nd:tidy up the tophits data
  external_IV<-read.table(paste0('./Taiwan/MVMR/IV_Taiwan/',pheno,'_GCTA_IV.txt'),header = T)
  external_IV<-external_IV[,c(1:3)] #CHR SNP BP
  external_IV<-merge(external_IV,gwas_exposure_format1,by='SNP',all.x=T)
  external_IV<-subset(external_IV,select=c('exposure','SNP','beta.exposure','se.exposure','effect_allele.exposure','other_allele.exposure',
                                           'eaf.exposure','pval.exposure','mr_keep.exposure','pval_origin.exposure','id.exposure'))
  nipt_IV<-read.table(paste0('./GCTA_file_230314/',pheno,'_IV_clean.txt'),header = T)
  nipt_IV<-nipt_IV[,c(1:3)]
  nipt_IV<-merge(nipt_IV,gwas_exposure_format2,by='SNP',all.x = T)
  nipt_IV<-subset(nipt_IV,select = c('exposure','SNP','beta.exposure','se.exposure','effect_allele.exposure','other_allele.exposure',
                                     'eaf.exposure','pval.exposure','mr_keep.exposure','pval_origin.exposure','id.exposure'))
  
  #3rd:load full GWAS for exposure, but in .outcome data
  colnames(gwas_exposure_format1)<-c('SNP','beta.outcome','se.outcome','effect_allele.outcome','other_allele.outcome','eaf.outcome',
                                     'pval.outcome','outcome','mr_keep.outcome','pval_origin.outcome','id.outcome')
  colnames(gwas_exposure_format2)<-c('SNP','beta.outcome','se.outcome','effect_allele.outcome','other_allele.outcome','eaf.outcome',
                                     'pval.outcome','outcome','mr_keep.outcome','pval_origin.outcome','id.outcome')
  
  #4rd:check list objects for exposures
  tophits_list<-list(external_IV,nipt_IV)
  full_gwas_list<-list(gwas_exposure_format1,gwas_exposure_format2)
  
  #5rd:multivariable MR with TwosampleMR
  exposure_dat<-get_mv_exposures(tophits_list,full_gwas_list)
  print(head(exposure_dat$SNP))
  
  #outcome data list
  outcome_list<-list.files('./MR/BBJ_summary/BBJ/',pattern = '.gwas.summary$')
  for(out in outcome_list){
    out_pheno<-sub('.gwas.summary','',out)
    print(out_pheno)
    outcome_dat<-fread(paste0('./MR/BBJ_summary/BBJ/',out),header = T)
    outcome_dat$phenotype<-out_pheno
    outcome_dat<-subset(outcome_dat,p!='NA')
    outcome_dat<-subset(outcome_dat,SNP!="")
    outcome_dat<-format_data(
      dat = outcome_dat,
      type = 'outcome',
      snps = exposure_dat$SNP,
      header = TRUE,
      phenotype_col = 'phenotype',
      snp_col = 'SNP',
      beta_col = 'beta',
      se_col = 'se',
      effect_allele_col = 'A1',
      other_allele_col = 'A2',
      eaf_col = 'eaf',
      pval_col = 'p'
    )
    #2smr
    mvdat<-mv_harmonise_data(exposure_dat,outcome_dat)
    res_bmi<-mv_multiple(mvdat)
    #print(res_bmi$results)
    results_2smr<-res_bmi$result %>%
      split_outcome()%>%
      separate(outcome,'outcome',sep = "[(]")%>%
      mutate(outcome=stringr::str_trim(outcome))%>%
      generate_odds_ratios()%>%
      select(-id.exposure,-id.outcome)
    fwrite(results_2smr,paste0('./Taiwan/MVMR/results_0816/',pheno,'/',pheno,'_',out_pheno,'_2smr_results.txt'),col.names = T,row.names = F,sep = '\t',quote = F)
    
    #MVMR
    mvmr_input<-make_mvmr_input(exposure_dat = exposure_dat,outcome.data = outcome_dat,outcome.id.mrbase = NULL)
    mvmr_out<-format_mvmr(BXGs = mvmr_input$XGs %>% select(contains('beta')), #exposure beta
                          BYG = mvmr_input$YG$beta.outcome, #outcome beta
                          seBXGs = mvmr_input$XGs %>% select(contains('se')), #exposure se
                          seBYG = mvmr_input$YG$se.outcome, #outcome se
                          RSID = mvmr_input$XGs$SNP) #SNPs
    #head(mvmr_out)
    mvmr_res<-ivw_mvmr(r_input = mvmr_out)
    #print(mvmr_res)
    result_mvmr<-mvmr_res%>%
      tidy_mvmr_output()%>%
      mutate(exposure=mvmr_input$exposures,
             outcome=out_pheno)%>%
      select(exposure,outcome,everything())
    fwrite(result_mvmr,paste0('./Taiwan/MVMR/results_0816/',pheno,'/',pheno,'_',out_pheno,'_mvmr_results.txt'),col.names = T,row.names = F,sep = '\t',quote = F)
    sres<-strength_mvmr(r_input=mvmr_out,gencov=0)
    sres$exposure<-pheno
    sres$outcome<-out_pheno
    fwrite(sres,paste0('./Taiwan/MVMR/results_0816/',pheno,'/',pheno,'_',out_pheno,'_mvmr_sres.txt'),col.names = T,row.names = T,sep = '\t',quote = F)
    pres<-pleiotropy_mvmr(r_input=mvmr_out,gencov=0)
    pres1<-data.frame(Q_statistic=pres$Qstat[1],pval=pres$Qpval[1],exposure=pheno,outcome=out_pheno)
    fwrite(pres1,paste0('./Taiwan/MVMR/results_0816/',pheno,'/',pheno,'_',out_pheno,'_mvmr_pres.txt'),col.names = T,row.names = F,sep = '\t',quote = F)
  }
  
