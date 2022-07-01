#get list of positions by scaf from map
for i in {1..16}
do
awk '{print $1}' Scaf${i}.recomb > Scaf${i}_map.pos
awk -v scaf="Scaffold_${i}" '($1 == scaf)' ../commongarden_and_PNAS_snps_missing2.bim > commongarden_and_PNAS_snps_missing2.scaf${i}
awk -F"\t" 'NR==FNR{e[$1]=1;next};e[$4]' Scaf${i}_map.pos commongarden_and_PNAS_snps_missing2.scaf${i} | awk '{print $2}' > scaf${i}.plinkpos
done

#comparing between homogenous allopatric pops (N versus K)
plink --bfile ../commongarden_and_PNAS_snps_missing2 --fst case-control --allow-extra-chr --allow-no-sex
#keep only sites that have higher fst then the mean
#awk -F"\t" '$5 > 0.1977654' plink.fst |grep -v "nan" | awk '{print $2}' > ancestry_informative_sites.txt
awk -F"\t" '$5 > 0.4' plink.fst |grep -v "nan" | awk '{print $2}' > ancestry_informative_sites.txt

#get AFs from ref populations
for i in {1..16}
do
cd Scaffold_${i}
grep "Scaffold_${i}:" ../ancestry_informative_sites.txt > scaf${i}.plinkpos
plink --bfile ../../commongarden_and_PNAS_snps_missing2 --keep ../K.fam --out varrudisAF_mergedSNPs_Scaf${i} --freq --allow-extra-chr --chr Scaffold_${i} --extract scaf${i}.plinkpos
plink --bfile ../../commongarden_and_PNAS_snps_missing2 --keep ../N.fam --out vartubercAF_mergedSNPs_Scaf${i} --freq --allow-extra-chr --chr Scaffold_${i} --extract scaf${i}.plinkpos
awk '{print $5}' varrudisAF_mergedSNPs_Scaf${i}.frq |  tail -n+2 | awk '{$1 = sprintf("%.3e", $1)} 1' > rudis_scaf${i}.AAF
awk '{print $5}' vartubercAF_mergedSNPs_Scaf${i}.frq | tail -n+2 | awk '{$1 = sprintf("%.3e", $1)} 1' > tuber_scaf${i}.AAF
cd ../
done

#extract just freqs
for i in {1..16}
do
awk '{print $5}' varrudisAF_mergedSNPs_Scaf${i}.frq | tail -n+2 > varrudisAF_mergedSNPs_Scaf${i}.AAfrq
awk '{print $5}' vartubercAF_mergedSNPs_Scaf${i}.frq | tail -n+2 > vartubercAF_mergedSNPs_Scaf${i}.AAfrq
done


#sci format
for i in {1..16}
do
awk '{$1 = sprintf("%.3e", $1)} 1' vartubercAF_mergedSNPs_Scaf${i}.AAfrq > tuber_scaf${i}.AAF
awk '{$1 = sprintf("%.3e", $1)} 1' varrudisAF_mergedSNPs_Scaf${i}.AAfrq > rudis_scaf${i}.AAF
done

#get 012 output of matched positions, by population
#format missing data accordingly
for i in {1..16}
do
cd Scaffold_${i}
        while read pop
        do
        plink --bfile ../../commongarden_and_PNAS_snps_missing2 --recodeA --chr Scaffold_${i} --allow-extra-chr --out commongarden_and_PNAS_snps_missing_Scaf${i}_${pop} --extract scaf${i}.plinkpos --keep ../${pop}.fam --allow-no-sex
	cut -d" " -f7- commongarden_and_PNAS_snps_missing_Scaf${i}_${pop}.raw | tail -n+2 | sed 's/NA/-1/g' | tr " " "\t" > commongarden_and_PNAS_snps_missing_Scaf${i}_${pop}.geno
        done < ../listoffams
cd ../
done


#fix the positional info
for i in {1..9}
do
pop=`head -n1 listoffams`
cd Scaffold_${i}
	cut -d" " -f7- commongarden_and_PNAS_snps_missing_Scaf${i}_${pop}.raw | head -n+1 | sed 's/Scaffold_.://g' |  sed 's/_.....//g' | tr " " "\n" > commongarden_and_PNAS_snps_missing_Scaf${i}.pos
cd ..
done

#and again
for i in {10..16}
do
pop=`head -n1 listoffams`
cd Scaffold_${i}
cut -d" " -f7- commongarden_and_PNAS_snps_missing_Scaf${i}_${pop}.raw | head -n+1 | sed 's/Scaffold_..://g' |  sed 's/_.....//g' | tr " " "\n" > commongarden_and_PNAS_snps_missing_Scaf${i}.pos
cd ../
done

#impute genetic map for lamp
for i in {1..16}
do
cd Scaffold_${i}
cut -d":" -f2 scaf${i}.plinkpos | cut -d"_" -f1 > scaf${i}.justpos
Rscript /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/imputemap.R /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/Scaffold_${i}_finalformat_recombrate.txt scaf${i}.justpos Scaf${i}.recomb

awk '{print $3}' Scaf${i}.recomb > Scaf${i}.recombrate
cd ..
done

#format configure files, esp with pop specific mean ancestry estimates
paste <(cut -f1 ancestry_bypop.forauto ) <(cut -f2- ancestry_bypop.forauto | awk -F"," '{printf " %.2f,%.2f\n", $1,$2}') > ancestry_bypop_rounded.forauto #round ancestry estimates, too many decimals was throwing errors

for i in {1..16}
do
cd Scaffold_${i}
        while read pop
        do
        anc=`grep $pop ../ancestry_bypop_rounded.forauto | awk '{print $2}'`
        cat ../config.test | sed "s/Pontoon/$pop/g" | sed "s/0.8528583,0.1471417/$anc/g" | sed "s/scaf10/scaf$i/g" | sed "s/Scaf10/Scaf$i/g" > ${pop}_scaf${i}.config
        done < ../listoffams
cd ../
done

#run lamp
while read pop
do
seq 1 16 | parallel -j 16 "cd Scaffold_{} && ~/software/lamp.x86_64-unknown-linux-gnu.REL_2_5/bin/lamp ${pop}_scaf{}.config"
done < listoffams


###
#use R to plot and get avrgs


