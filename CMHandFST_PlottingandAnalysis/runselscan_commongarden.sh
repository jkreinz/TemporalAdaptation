#go from haps to vcf
for i in {1..16}
do
paste -d" " <(awk '{print $1":"$3"_"$4"_"$5}' ../Scaffold_${i}_commongarden_phased_flipped.haps) <(cut -d" " -f2- ../Scaffold_${i}_commongarden_phased_flipped.haps) > Scaffold_${i}.hap
bgzip Scaffold_${i}.hap
cp ../Scaffold_${i}_commongarden_phased.vcf.sample ./Scaffold_${i}.samples
/ohta1/julia.kreiner/software/bcftools/bcftools convert --hapsample2vcf Scaffold_${i} > Scaffold_${i}.vcf
done

#format recombination map
for i in {1..16}; do grep -v "#" Scaffold_${i}.vcf | awk '{print $1 "\t" $2}' > Scaffold_${i}.sites; done
for i in {1..16}
do 
Rscript /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/ldhat/LDhat_workflow/amaranth/interpolate.R Scaffold_${i}.sites /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/ldhat/LDhat_workflow/amaranth/Scaffold_${i}.map Scaffold_${i}_selscan.map
done

#match map and vcf positions
for i in {1..16}
do
awk '{print $1 "\t" $4}' Scaffold_${i}_selscan.map > Scaffold_${i}_selscan_map.positions
#make environment specific vcfs
vcftools --vcf Scaffold_${i}.vcf --keep /ohta1/julia.kreiner/waterhemp/commongarden/regional_vcfs/bychrom/byenv_analysis/ag_samps --recode --out ag_Scaffold_${i} --positions Scaffold_${i}_selscan_map.positions &
vcftools --vcf Scaffold_${i}.vcf --keep /ohta1/julia.kreiner/waterhemp/commongarden/regional_vcfs/bychrom/byenv_analysis/nat_samps --recode --out nat_Scaffold_${i} --positions Scaffold_${i}_selscan_map.positions
done

#run selscan
seq 1 16 | parallel -j 4 "selscan --xpehh --vcf-ref nat_Scaffold_{}.recode.vcf --vcf ag_Scaffold_{}.recode.vcf --map Scaffold_{}_selscan.map --out XPEHH_nat_ag_Scaffold_{} --threads 10"
seq 1 16 | parallel -j 4 "selscan --xpnsl --vcf-ref nat_Scaffold_{}.recode.vcf --vcf ag_Scaffold_{}.recode.vcf --out xpnsl_nat_ag_Scaffold_{} --threads 10"

