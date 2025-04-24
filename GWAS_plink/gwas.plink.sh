#/usr/bin/sh
#---
#Perform genome wide association analysis using plink
#Step1 perform inverse rank transformation (quantile transformation) for phenotypic data
#Step2 generate plink shell
#---

outdir=./output
mkdir -p $outdir


script1=./2022106_GWASPIP/bin/get.pheno.py
script2=./2022106_GWASPIP/bin/normalized_phenotyp.r
script3=./2022106_GWASPIP/bin/g.worksh_1.py
R=`which Rscript`

allpheno_table=./input/coding.txt

allsample=./input/pheno.table
vcflist=./input/NIPT.imputation.vcf.bychr.list
covar9=./input/plink.PC1-10.bmi.xCov.txt
hweinfo=./NIPT_all_snp.info.gz
imputevcf=./hg38.all.stitch.bed.sorted.dbsnp.gwas.clinvar.snpeff.vcf.gz

for pheno in SBP DBP GDM GT
do
    mkdir -p $outdir/$pheno/input/
    mkdir -p $outdir/$pheno/output/plink/
    mkdir -p $outdir/$pheno/output/plink_merge/
    mkdir -p $outdir/$pheno/bin

#---
#Step1 clean and normalize phenotypes
#---
    case $pheno in
        SBP|DBP
           model='linear';
	   python $script1 $allpheno_table $allsample $pheno >$outdir/$pheno/input/${pheno}_pheno.table
            $R $script2 $outdir/$pheno/input/${pheno}_pheno.table Y -9 -9 $outdir/$pheno/input/;
		    phenotable=$outdir/$pheno/input/${pheno}_pheno.table.rmout.qtrans;
           echo 'Linear regression done';;
        GDM|GT
           model='logistic';
	   python $script1 $allpheno_table $allsample $pheno >$outdir/$pheno/input/${pheno}_pheno.table
           phenotable=$outdir/$pheno/input/${pheno}_pheno.table;
           echo 'Logistic regression done';;
    esac

#---
#Step 2 generate plink shell
#---
    case $pheno in
        SBP|DBP
            python $script3 $vcflist $outdir/$pheno/output/plink/ $covar9 $pheno $phenotable $hweinfo linear > $outdir/$pheno/bin/step1.plink.work.sh;;
        GDM|GT
            python $script3 $vcflist $outdir/$pheno/output/plink/ $covar9 $pheno $phenotable $hweinfo logistic > $outdir/$pheno/bin/step1.plink.work.sh;;
    esac

done
