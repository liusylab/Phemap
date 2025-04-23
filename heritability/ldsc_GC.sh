#!/bin/sh
pheno_path=./heritability/LDSC/sumstats_gz
result_path=./heritability/LDSC/genetic_correlation
mkdir -p ${result_path}/shell
mkdir -p ${result_path}/result
n=1
for trait1_file in `ls ${pheno_path}/*.sumstats.gz`
do
trait1_name=${trait1_file##*/}
trait1=${trait1_name%%.sumstats.gz*}
for trait2_file in `ls ${pheno_path}/*.sumstats.gz`
do
trait2_name=${trait2_file##*/}
trait2=${trait2_name%%.sumstats.gz*}
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10240
#SBATCH --partition=cpu
#SBATCH --job-name=${trait1}_${trait2}
#SBATCH --output=${result_path}/shell/slurm_${trait1}_${trait2}.sh
echo \"process will start at :\"
date
ldsc.py --rg ${trait1_file},${trait2_file} --ref-ld-chr ./GCTA_ldscore/10k_ldscore/chr --w-ld-chr ./GCTA_ldscore/10k_ldscore/chr --out ${result_path}/result/${trait1}.gz.sumstats.gz_${trait2}.gz.sumstats.gz
echo \"process end at : \"
date" > ${result_path}/shell/sbatch_${trait1}_${trait2}.txt
let n+=1
echo $n
done
done