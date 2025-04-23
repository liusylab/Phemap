#!/bin/sh
pheno_path=./heritability/LDSC/sumstats_gz
result_path=./heritability/LDSC/heritability
mkdir -p ${result_path}/shell
mkdir -p ${result_path}/result
n=1
for trait_file in `ls ${pheno_path}/*.sumstats.gz`
do
trait_name=${trait_file##*/}
trait=${trait_name%%.sumstats.gz*}
echo ${trait}
echo "#!/bin/sh
#SBATCH --nodes=1                  
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10240
#SBATCH --partition=cpu
#SBATCH --job-name=${trait}
#SBATCH --output=${result_path}/shell/slurm_${trait}.sh
echo \"process will start at :\"
date
ldsc.py --h2 ${trait_file} --ref-ld-chr ./GCTA_ldscore/10k_ldscore/chr --w-ld-chr ./GCTA_ldscore/10k_ldscore/chr --out ${result_path}/result/${trait}_h2
echo \"process end at : \"
date" > ${result_path}/shell/sbatch_${trait}.txt
let n+=1
echo $n
done