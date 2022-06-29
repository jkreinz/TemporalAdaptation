freq=$1

cut -d" " -f7- AF_perms/ag_analysis_redo/cg_${freq}.raw | \
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
}' > AF_perms/ag_analysis_redo/cg_${freq}.012
