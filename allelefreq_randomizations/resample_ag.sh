itt=$1

rm AF_perms/ag_analysis_redo/resamples/cg_${itt}_out.txt

while read freq occur
do
paste <(echo $freq) <(shuf -n $occur AF_perms/ag_analysis_redo/cg_${freq}_f.012) >> AF_perms/ag_analysis_redo/resamples/cg_${itt}_out.txt
done < bins_agAFs_n154_occurence.txt

cut -d" " -f-2 AF_perms/ag_analysis_redo/resamples/cg_${itt}_out.txt | awk 'split($1,a,"_") {print a[1]"_"a[2]}' > AF_perms/ag_analysis_redo/resamples/cg_${itt}.snppos
