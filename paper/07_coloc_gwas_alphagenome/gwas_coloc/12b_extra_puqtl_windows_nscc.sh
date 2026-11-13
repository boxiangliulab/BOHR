#!/bin/bash
# LD-window GWAS coloc, step 12b (NSCC): LD windows of the 122 puQTL leads of the non-significant MRAS units (MRAS panel)
# that are not in Figure 4's lead table. Same recipe as Figure 4's 01_ldwin_nscc.pbs: PLINK 1.9 --r2 on the MAF >= 0.05
# panel bfile, the tissue's donors, r2 >= 0.1 within 1 Mb; kept per lead: n partners and leftmost/rightmost partner.
set -eo pipefail
P=/home/users/nus/e1327989/bin/plink; R=/home/users/nus/e1327989/YuhaoDong/GTEx/LD_ref_v9_plink19; W=/scratch/users/nus/e1327989/fig6_ldwin
mkdir -p $W/extra_win/out; cd $W/extra_win
tail -n +2 $W/puqtl_leads_missing.tsv | awk -F'\t' '{print > ("in." $1 "." $3 ".txt")}' ; for f in in.*.txt; do cut -f2 $f > $f.snps; done
echo -e "tissue\tsnp\tpos\tn_partners\tmin_partner_pos\tmax_partner_pos" > $W/extra_win/puqtl_extra_windows.tsv
for f in in.*.txt; do x=${f#in.}; x=${x%.txt}; T=${x%.chr*}; C=${x##*.}
  $P --bfile $R/bfile/$C --keep-allele-order --keep $R/donors/$T.txt --ld-snp-list $f.snps --r2 --ld-window 99999 \
     --ld-window-kb 1000 --ld-window-r2 0.1 --threads 4 --memory 8000 --out out/$x > out/$x.stdout 2>&1
  awk -v T=$T 'NR>1 { s=$3; p=$5; if (!(s in n)) { n[s]=0; lo[s]=p; hi[s]=p; bp[s]=$2 } if ($6 != s) n[s]++; if (p<lo[s]) lo[s]=p; if (p>hi[s]) hi[s]=p }
       END { for (s in n) print T "\t" s "\t" bp[s] "\t" n[s] "\t" lo[s] "\t" hi[s] }' out/$x.ld >> $W/extra_win/puqtl_extra_windows.tsv
  rm -f out/$x.ld
done
echo "extra windows: $(( $(wc -l < $W/extra_win/puqtl_extra_windows.tsv) - 1 )) of $(( $(wc -l < $W/puqtl_leads_missing.tsv) - 1 ))"
