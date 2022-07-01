#estimate CMH genome-wide
plink --allow-extra-chr --bfile commongarden_finalfilt_maf01 --make-founders --mh --out commongarden_finalfilt_maf01 --within commongarden_finalfilt_maf01.fam

#run clumping algorithm on CMH results
plink --allow-extra-chr --bfile commongarden_finalfilt_maf01 --clump ag_v_nat_withinpairs_maf01.cmh --clump-kb 1000 --clump-r2 .25 --exclude commongarden_finalfilt_maf01.dupvar2 --make-founders --out ag_v_nat_withinpairs_maf01_1mb_r25

#estimate Fst genome-wide
vcftools --gzvcf ../allchroms_commongarden_dustm_hardfilt_AC0_QualDP_biall_AB_sorted.vcf.gz --weir-fst-pop ag_samps --weir-fst-pop nat_samps --keep ag_samps --keep nat_samps --out ag_v_nat

