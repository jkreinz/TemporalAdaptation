for i in {1..100}
do

slim temporal_slim_sqrtfitness.script #run a single iteration of simulation to test correlation between estimated s and true s with spatial/temporal sampling scheme
rm *.vcf.gz

for i in `ls *.vcf | sed 's/_sample.vcf//g' | sed 's/t//g'`
do
bgzip t${i}_sample.vcf && tabix -p vcf t${i}_sample.vcf.gz  #individual vcfs outputted from SLIM from each temporal/spatial sample
done

vcf-merge *.vcf.gz > allsamples.vcf 

cat <(grep "#" allsamples.vcf) <(grep 'MT=2' allsamples.vcf) | grep -v "MULTIALLELIC" > m2_allsamples.vcf #drop multiallelic sites

sed -e 's/\t\./\t0\|0/g' m2_allsamples.vcf > m2_allsamples_homoref.vcf #change ./. notation to 0/0
grep -v "#" m2_allsamples_homoref.vcf | awk '{print $8}' | awk 'split($1,a,";") {print a[9]}' | sed 's/S=//g' > selection_coefficients.txt #extract true s coefficients

vcftools --vcf m2_allsamples_homoref.vcf --012 --out m2_allsamples_nomulti

Rscript getcorr_estimatetrues.R  >> correlation_coefficients_sqrt.out #estimate s coefficients & perform correlation
done
