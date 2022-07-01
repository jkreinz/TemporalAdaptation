#!/bin/bash
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=80
#SBATCH --time=12:00:00
#SBATCH --job-name freebayes_parallel

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=1

module load NiaEnv
module load CCEnv
module load nixpkgs/16.09
module load gcc/7.3.0
module load parallel/20160722
module load freebayes/1.2.0

#run on multiple nodes
HOSTS=$(scontrol show hostnames $SLURM_NODELIST | tr '\n' ,)

#generate 100KB chunks across the genome
#regionsfile=$(/home/w/wrighste/kreinerj/freebayes-1.3.2/scripts/fasta_generate_regions.py /home/w/wrighste/kreinerj/amaranth/AMATA_finishedtohypo_renamedfinal.fasta 100000)


# iterate over regions using gnu parallel to dispatch jobs
cat /home/w/wrighste/kreinerj/amaranth/100kbregions | parallel --env OMP_NUM_THREADS,PATH,LD_LIBRARY_PATH --joblog slurm-$SLURM_JOBID.log -k -j $SLURM_NTASKS_PER_NODE -S $HOSTS --wd $PWD "freebayes -f /home/w/wrighste/kreinerj/amaranth/AMATA_finishedtohypo_renamedfinal.fasta --use-best-n-alleles 4 /scratch/w/wrighste/kreinerj/bams/*dd.bam > /scratch/w/wrighste/kreinerj/region_vcfs/{}.vcf " --region {}


#cat /home/w/wrighste/kreinerj/amaranth/listofregions.bed | parallel --verbose -j 8 "java -jar gatk \
        #        -T HaplotypeCaller \
        #        -R /home/w/wrighste/kreinerj/amaranth/AMATA_finishedtohypo_renamedfinal.fasta \
        #        -ploidy 1 \
        #        -ERC GVCF \
        #        -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
        #        -L {}
#        -I {}.sort.bam -o $src/SNPs/{}.gvcf"


#join seperately called regions, taking only first header and removing any duplicates
# cat *vcf| python2.7 /scratch/w/wrighste/kreinerj/software/freebayes/vcflib/scripts/vcffirstheader \
        #    | /scratch/w/wrighste/kreinerj/software/freebayes/vcflib/bin/vcfstreamsort -w 1000 | /scratch/w/wrighste/kreinerj/software/freebayes/vcflib/bin/vcfuniq >> allsites.vcf #
