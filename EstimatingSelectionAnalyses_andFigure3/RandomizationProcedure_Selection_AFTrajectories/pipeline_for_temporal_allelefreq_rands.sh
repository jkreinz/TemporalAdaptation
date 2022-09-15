                                                                                                                                                                                                                                                                                                                                                                                                                                                                         1,1           Top
#########################
#this pipeline uses a combination of VCFtools, PLINK, and bash to get an expectation of allele frequency change based on a starting (contemporary) allele frequency distribution

#the basic workflow is as follows:
#1) get rounded allele freqs for focal SNPs in environment of interest
#2) get rounded allele freqs for all sites across the genome
### make sure we are tracking the ag allele (those at a greater freq in ag than nat) by flipping to track the ref when alt is less freq in ag
#3) randomly sample sites across the genome to match the same starting (focal) allele frequency distribution and sample size
#4) for each random sampling procedure, extract genotypes from historical samples at those sites

##############
#step 1
##############
#get list of individuals
awk '{print $2}' ../../byenv_analysis/ag.fam > ag.inds

#get freqs for focal sites (n=154) within those individuals (here agricultural inds)
vcftools --vcf CMH_clumped_outliers_cg.vcf --freq --out ag_herbmatched --keep ag.inds --positions herb_cmh_FDR10p_clumped_matched_refalt.pos

#get frequency of observed focal alleles in ag by subsetting (taking ref or alt based on which is more common in ag)
plink2 --bfile CMH_clumped_outliers_cg --freq --keep nat.fam --nonfounders --out nat_clump_outliers --allow-extra-chr #get allele frequency at all sites in nat
plink2 --bfile CMH_clumped_outliers_cg --freq --keep ag.fam --nonfounders --out ag_clump_outliers --allow-extra-chr #get allele frequency at all sites in nat
paste <( cut -f2-5 ag_clump_outliers.afreq) <(cut -f3-5 nat_clump_outliers.afreq) > ag_nat_cmh_clump_freqs.txt
#check which top CMH outlier sites have higher freq in nat compared to ag
awk '$4 < $7' ag_nat_cmh_clump_freqs.txt > allelestoflip_clump_agassoc.txt #ones where reference is the ag allele

cat <(grep -Ff <(cut -f1 allelestoflip_clump_agassoc.txt | cut -d"[" -f1 | tr ":" "\t") ag_herbmatched.frq | awk '{print $1 "\t" $2 "\t" $5}' | tr ":" "\t" ) <(grep -vFf <(cut -f1 allelestoflip_clump_agassoc.txt | cut -d"[" -f1 | tr ":" "\t") ag_herbmatched.frq | awk '{print $1 "\t" $2 "\t" $6}' | tr ":" "\t")  | grep -v "CHROM" > ag_freq_contemporary_n154.frq


#bin agricultural allele frequencies in agricultural environments (round to the nearest 100th) for matching with permuted set
R
library(data.table)
frq<-fread("ag_freq_contemporary_n154.frq",header=F)
names(frq)[4]<-"MAF"
mean(frq$MAF) #check the mean allele freq
write.table(unique(round(frq$MAF,digits=2)),"bins_agAFs_n154.txt",quote=F,row.name=F)
write.table(table((round(frq$MAF,digits=2))), "bins_agAFs_n154_occurence.txt",quote=F,row.name=F,col.name=F)
q()
n

##############
#step 2
##############
#for generating a set of random loci to estimate AF change over, we first need to make sure they are present in both herbarium and contemporary data
#first copy herb/historical files over
cp /ohta2/julia.kreiner/herbarium/femaleref/redupped_rescaled/regional_vcfs/bychrom/allchroms_herblastrerun_QUALDP_biSNP.fam .
cp /ohta2/julia.kreiner/herbarium/femaleref/redupped_rescaled/regional_vcfs/bychrom/allchroms_herblastrerun_QUALDP_biSNP.bim .
cp /ohta2/julia.kreiner/herbarium/femaleref/redupped_rescaled/regional_vcfs/bychrom/allchroms_herblastrerun_QUALDP_biSNP.bed .

#fix snp names to match contemporary set (plink based naming)
awk 'split($2,a,"_") {print $1 "\t" a[1]"_"a[2]"[b37]"$5","$6 "\t" 0 "\t" $4 "\t"$5 "\t" $6}' /ohta2/julia.kreiner/herbarium/femaleref/redupped_rescaled/regional_vcfs/bychrom/allchroms_herblastrerun_QUALDP_biSNP.bim > allchroms_herblastrerun_QUALDP_biSNP_mod.bim
mv allchroms_herblastrerun_QUALDP_biSNP.bed allchroms_herblastrerun_QUALDP_biSNP_mod.bed
mv allchroms_herblastrerun_QUALDP_biSNP.fam  allchroms_herblastrerun_QUALDP_biSNP_mod.fam

#make sure we're looking at alleles that are present in both herb and cg datasets
join -j2 <( sort -k2,2 allchroms_herblastrerun_QUALDP_biSNP_mod.bim) <( sort -k2,2 ../commongarden_finalfilt.bim) | awk '{print $1}' > cg_herb_alleles_overlappingpos_ag.txt #1.66M alleles in common

#get ag observed AFs across the genome
mkdir AF_perms/ag_analysis_redo
plink --bfile ../commongarden_finalfilt --freq --keep ../../byenv_analysis/ag.fam --out AF_perms/ag_analysis_redo/CG_ag --allow-extra-chr --nonfounders --extract cg_herb_alleles_overlappingpos_ag.txt
#get nat observed AFs across the genome
plink --bfile ../commongarden_finalfilt --freq --keep ../../byenv_analysis/nat.fam --out AF_perms/ag_analysis_redo/CG_nat --allow-extra-chr --nonfounders --extract cg_herb_alleles_overlappingpos_ag.txt

#constrain to sites where freq of alt greater in ag compared to nat, and for sites > 0.5, freq of ref greater in nat compared to ag (since ref would the be higher freq in ag > nat)
paste AF_perms/ag_analysis_redo/CG_ag.frq AF_perms/ag_analysis_redo/CG_nat.frq | awk '$5 > $11 {print $1 "\t" $2 "\t" $5}' > CG_ag.frq.constrained #alternate (genome-wide) ag alleles
paste AF_perms/ag_analysis_redo/CG_ag.frq AF_perms/ag_analysis_redo/CG_nat.frq | awk '$5< $11 {print $1 "\t" $2 "\t" $5}' > CG_ag.frq.constrained.flipped #reference (genome-wide) ag alleles

#round genome-wide frequencies in ag and nat to the nearest hundredth for matching
awk '{$3 = sprintf("%.2f",$3)} {print $1 "\t" $2 "\t" $3}' CG_ag.frq.constrained | tail -n+2 > AF_perms/ag_analysis_redo/CG_ag.frq.constrained.binned
awk '{$3 = sprintf("%.2f",$3)} {print $1 "\t" $2 "\t" $3}' CG_ag.frq.constrained.flipped | tail -n+2 > AF_perms/ag_analysis_redo/CG_ag.frq.constrained.flipped.binned

#produce seperate files for each binned AF (for ease later)
#remove header
tail -n+2 bins_agAFs_n154.txt > bins_agAFs_n154_vec.txt
while read freq
do
grep -F $freq AF_perms/ag_analysis_redo/CG_ag.frq.constrained.binned > AF_perms/ag_analysis_redo/${freq}_loci.txt
done < bins_agAFs_n154_vec.txt

#again for SNPS with freq > 0.5 (in which case track reference allele)
sort bins_agAFs_n154_vec.txt | awk '$1 > .49' > bins_agAFs_n154_vec_uppers.txt
awk '$3 = sprintf("%.2f",1-$3) {print $1 "\t" $2 "\t" $3}' AF_perms/ag_analysis_redo/CG_ag.frq.constrained.flipped.binned > AF_perms/ag_analysis_redo/CG_ag.frq.constrained.flipped.binned.f
while read freq
do
grep -F $freq AF_perms/ag_analysis_redo/CG_ag.frq.constrained.flipped.binned.f > AF_perms/ag_analysis_redo/${freq}_loci.txt
done < bins_agAFs_n154_vec_uppers.txt

#produce a 012 file for each frequency class from common garden data
#first remake plink files with just matching alleles, and just ag samples
plink --bfile ../commongarden_finalfilt --keep ../../byenv_analysis/ag.fam --out CG_herbmatchedsnps_genomewide --allow-extra-chr --nonfounders --extract cg_herb_alleles_overlappingpos_ag.txt --make-bed

#now extract 012 in agricultural samples for each frequency class
while read freq
do
awk '{print $2}' AF_perms/ag_analysis_redo/${freq}_loci.txt > AF_perms/ag_analysis_redo/intermed.pos
plink --bfile CG_herbmatchedsnps_genomewide --out AF_perms/ag_analysis_redo/cg_${freq}  --allow-extra-chr --nonfounders --recodeA --extract  AF_perms/ag_analysis_redo/intermed.pos
done < bins_agAFs_n154_vec.txt

#transpose 012
cat bins_agAFs_n154_vec.txt | parallel -j 20 "bash transposeraw_ag.sh {}"

#flip 0 and 2's for AFs > 0.5 (since want to count the number of reference alleles, not alternate)
while read freq
do
cat AF_perms/ag_analysis_redo/cg_${freq}.012 | sed 's/ 2/ 3/g' | sed 's/ 0/ 2/g' | sed 's/ 3/ 0/g' > AF_perms/ag_analysis_redo/cg_${freq}_f.012
done < bins_agAFs_n154_vec_uppers.txt

#cut snpnames for lower AFs (i.e. alternate ag alleles) - we'll use this below to track alleles from each ref or alt set
awk '$1 < 0.51' bins_agAFs_n154_vec.txt > bins_agAFs_n154_vec_lowers.txt
while read freq
do
cat AF_perms/ag_analysis_redo/cg_${freq}.012 > AF_perms/ag_analysis_redo/cg_${freq}_f.012
done < bins_agAFs_n154_vec_lowers.txt

##############
#step 3
##############
#get permed itts, resampling ag-frequency matched alleles in agricultural environments
mkdir AF_perms/ag_analysis_redo/resamples
seq 1 1000 | parallel -j 20 "bash resample_ag.sh {}" #keep track of whether an allele was sampled from the flipped or unflipped set


##############
#step 4
##############
#match resampled snps to herbarium data
#have to keep in mind the flipping procedure to make sure we track the correct allele (alt or ref)

#for each 1000 randomiziations...
for itt in {1..1000}
do
#0.48 is highest freq of alt (ag) allele in ag environments 
grep -F "0.48" -B100 AF_perms/ag_analysis_redo/resamples/cg_${itt}_out.txt | awk '{print $1 "\t" $2}' | sed -e 's/[0-9]*\.[0-9]*\t//g' | awk '{print $1}' | sed 's/..$//'  > AF_perms/ag_analysis_redo/resamples/lower_matches_${itt}.txt
#0.52 is lowest freq of alt (ref) allele in ag environments
grep -F "0.52" -A100 AF_perms/ag_analysis_redo/resamples/cg_${itt}_out.txt | awk '{print $1 "\t" $2}' | sed -e 's/[0-9]*\.[0-9]*\t//g' | awk '{print $1}' | sed 's/..$//' > AF_perms/ag_analysis_redo/resamples/upper_matches_${itt}.txt

#get genos in historical sequence dataset
plink --bfile allchroms_herblastrerun_QUALDP_biSNP_mod  --allow-extra-chr --nonfounders --recodeA --extract AF_perms/ag_analysis_redo/resamples/lower_matches_${itt}.txt --out AF_perms/ag_analysis_redo/resamples/lowermatches_${itt} &
plink --bfile allchroms_herblastrerun_QUALDP_biSNP_mod  --allow-extra-chr --nonfounders --recodeA --extract AF_perms/ag_analysis_redo/resamples/upper_matches_${itt}.txt --out AF_perms/ag_analysis_redo/resamples/uppermatches_${itt}

#transpose and merge lower (alt) and upper (ref) sets, both which are tracking the ag > nat allele
cat <(cut -d" " -f7- AF_perms/ag_analysis_redo/resamples/lowermatches_${itt}.raw | \
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' ) <(cut -d" " -f7- AF_perms/ag_analysis_redo/resamples/uppermatches_${itt}.raw | \
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' | sed 's/ 2/ 3/g' | sed 's/ 0/ 2/g' | sed 's/ 3/ 0/g' ) | sort -Vk1,1 > AF_perms/ag_analysis_redo/resamples/herb_ordered_${itt}_matches.012

#sort contepmorary snps to match historical order
cut -d$'\t' -f2 AF_perms/ag_analysis_redo/resamples/cg_${itt}_out.txt | sort -Vk1,1 > AF_perms/ag_analysis_redo/resamples/cg_ordered_${itt}_matches.012

done



##############
#calculating AF stats
##############
cd AF_perms/ag_analysis_redo/resamples/
mkdir results
cd results

#get list of inds
awk '{print $2}' /ohta1/julia.kreiner/waterhemp/commongarden/regional_vcfs/bychrom/plink/CMH/allchroms_herblastrerun_QUALDP_biSNP_mod.fam > herb.ind
awk '{print $2}' /ohta1/julia.kreiner/waterhemp/commongarden/regional_vcfs/bychrom/plink/CMH/CG_herbmatchedsnps_genomewide.fam > cg.ind
cp ../../../nat_analysis_redo/resamples/results/AFperm_statistics.R .
cp /ohta1/julia.kreiner/waterhemp/commongarden/regional_vcfs/bychrom/plink/CMH/AF_perms/resamples/results/3waymerged_sampleinfo.txt .

#R script for producing summary stats on each perm
seq 1 1000 | parallel -j 40 "Rscript AFperm_statistics.R {}"
cat summarystat_*.txt > allsummaries_agperms.txt


#quick print to screen of overall results across perms
R
sumz<-read.table("allsummaries_agperms.txt")
names(sumz)<-c("percentin", "meanchange", "joint_s", "joint_s_pre", "joint_s_post", "medianchange", "meanslope", "medianslope", "meansloperr", "mediansloperr", "meanpval", "medianpval")
quantile(sumz$meanchange, probs=c(0.025,0.975))
quantile(sumz$joint_s, probs=c(0.025,0.975))
quantile(sumz$joint_s_pre, probs=c(0.025,0.975))
quantile(sumz$joint_s_post, probs=c(0.025,0.975))
q()
