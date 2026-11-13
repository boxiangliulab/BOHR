#!/usr/bin/env python3
"""2026-10-09: MAGIC 2021 European fasting glucose / fasting insulin / 2-h glucose, from the GWAS
Catalog harmonised files (GRCh38: hm_chrom, hm_pos), converted with convert_gwas_catalog() copied verbatim from
new_GWAS/scripts/prepare_table_s4_gwas.py (the converter used for MAGIC HbA1c, GCST90002244), so the prepared files have
exactly the format and QC of the others. Output keeps the analysis file names (the content is now the European result)."""
import gzip, math, os, sys
W = "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/gwas_eur_20261009"

def convert_gwas_catalog(src, out, name, use_log_or=False):
    n_w = n_s = 0
    with gzip.open(src, 'rt') as fin, open(out, 'w') as fout:
        header = fin.readline().rstrip('\n').split('\t')
        col = {h: i for i, h in enumerate(header)}
        fout.write("SNPID_hg19 SNPID_hg38 BETA SE Pvalue A1FREQ R2 ALTfreq\n")
        for line in fin:
            f = line.rstrip('\n').split('\t')
            if len(f) < len(header):
                continue
            chrom = f[col['hm_chrom']]; pos = f[col['hm_pos']]
            if chrom in ('NA', '') or pos in ('NA', ''):
                n_s += 1; continue
            pval = f[col['p_value']]
            if pval in ('NA', ''):
                n_s += 1; continue
            if use_log_or:
                raise SystemExit("not used here")
            else:
                bs = f[col['hm_beta']]; ss = f[col['standard_error']]
                if bs in ('NA', '') or ss in ('NA', ''):
                    n_s += 1; continue
                try:
                    beta = float(bs); se = float(ss)
                except Exception:
                    n_s += 1; continue
                if se <= 0:
                    n_s += 1; continue
            a1 = f[col.get('hm_effect_allele_frequency', -1)] if 'hm_effect_allele_frequency' in col else 'NA'
            if a1 in ('', 'NA'):
                a1 = 'NA'
            fout.write(f"NA {chrom}:{pos} {beta:.6g} {se:.6g} {pval} {a1} NA NA\n")
            n_w += 1
    print(f"  {name}: {n_w} variants, {n_s} skipped")

JOBS = [("MAGIC_FG_EUR_GCST90002232.h.tsv.gz", "Fasting_Glucose_Chen_2021_Nature_Genetics_hg38.txt"),
        ("MAGIC_FI_EUR_GCST90002238.h.tsv.gz", "Fasting_Insulin_Chen_2021_Nature_Genetics_hg38.txt"),
        ("MAGIC_2hGlu_EUR_GCST90002227.h.tsv.gz", "2hGlu_Chen_2021_Nature_Genetics_hg38.txt")]
for src, name in JOBS:
    op = os.path.join(W, "prepared", name)
    if os.path.exists(op): print(name, "exists"); continue
    convert_gwas_catalog(os.path.join(W, "raw", src), op, name)
