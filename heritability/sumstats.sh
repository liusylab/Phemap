#!/bin/sh
pheno_path=./heritability/LDSC/ldsc_input_file
result_path=./heritability/LDSC/sumstats_gz
mkdir -p ${result_path}/shell
n=1
for pheno_file in `ls ${pheno_path}/*.txt`
do
pheno_name=${pheno_file##*/}
name=${pheno_name%%.txt*}
echo ${name}
n_row=`less ./heritability/LDSC/bin/phemap_sample_number.txt | grep -wn ${name} | awk -F ":" '{print $1}'`
num=`less ./heritability/LDSC/bin/phemap_sample_number.txt | awk 'NR=='${n_row}'{print $2}'`
echo ${num}
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10240
#SBATCH --partition=cpu
#SBATCH --job-name=${name}
#SBATCH --output=${result_path}/shell/slurm_${name}.sh
echo \"process will start at :\"
date
munge_sumstats.py --sumstats ${pheno_path}/${name}.txt --maf-min 0.01 --N ${num} --out ${result_path}/${name}
echo \"process end at : \"
date" > ${result_path}/shell/sbatch_${name}.txt
let n+=1
done