#use samtools to get fragment size
while read files; do java -jar ~/software/picard.jar CollectInsertSizeMetrics I=${files}.final.sorted_rmdup.rescaled.sorted.bam H=${files}.hist O=${files}.out; done < listoffiles
grep "insert size average" *stats > allsamp_insertlengths.txt
awk -F":" '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' allsamp_fraglengths.txt > allsamp_fraglengths_final.txt

#calculate allelic bias for two sets of alleles

awk -f ABforhets.awk reference_genos_n91_noheader.vcf > reference_AB_bygeno.txt
awk -f ABforhets.awk alternate_genos_n63_noheader.vcf > alternate_AB_bygeno.txt #will merge estimates for reference and alternate ag-alleles (remember, both can be ag since sometimes ref allele is more common in ag)

grep -v "#" herbarium_fullfilt_ancestryinformSNPs.recode.vcf > herbarium_fullfilt_ancestryinformSNPs_noheader.vcf #remove header
awk -f /ohta1/julia.kreiner/waterhemp/commongarden/regional_vcfs/bychrom/plink/CMH/ABforhets.awk herbarium_fullfilt_ancestryinformSNPs_noheader.vcf

#where ABforhets.awk script calculates the allelic bias for each genotype in the vcf (ref count over total read count)
BEGIN {OFS = "\t"}
{ printf $1"\t"$2"\t"$4"\t"$5"\t"
        for (i=10; i<=NF; ++i){
                split($i,a,":")
                if (toupper(a[1]) == "0/1" ) {
                        split(a[3],b,","); printf b[1]/(b[1]+b[2])"\t"}
                else
                        printf "NA" "\t"
                                }
        printf "\n"
}
