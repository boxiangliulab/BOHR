"""Figure 4 panel e, LD-window coloc, sensitivity "eQTL leads by recursive clumping" (author 2026-10-12),
step s1: the GTEx v11 significant eQTL variants (eQTLs.signif_pairs: p below the gene's nominal threshold) of every eGene
x tissue of the Fig. 4e analysis set (gene_base matched on the version-stripped ID).
Output: results/coloc_eqtl_ldwin/sens_clump/sig/<tissue>.tsv.gz (gene_base, chr, variant_id, pval_nominal).
Usage: python s1_sig_pairs.py <tissue>"""
import sys, os
import pandas as pd, pyarrow.parquet as pq

O = "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_ldwin"
D = "/home/e1327989/project/GTEx/Rawdata/gtex_v11/GTEx_Analysis_v11_eQTL"
tissue = sys.argv[1]
out = os.path.join(O, "sens_clump", "sig", tissue + ".tsv.gz")
if os.path.exists(out):
    sys.exit(0)
P = pd.read_csv(os.path.join(O, "pairs.tsv.gz"), sep="\t", usecols=["tissue", "gene_base", "chr", "egene"])
G = P[(P.tissue == tissue) & (P.egene)].drop_duplicates(["gene_base", "chr"])
t = pq.read_table(os.path.join(D, f"{tissue}.v11.eQTLs.signif_pairs.parquet"), columns=["phenotype_id", "variant_id", "pval_nominal"]).to_pandas()
t["gene_base"] = t.phenotype_id.str.replace(r"\..*$", "", regex=True)
t = t[t.gene_base.isin(G.gene_base)].merge(G[["gene_base", "chr"]], on="gene_base")
t = t[t.variant_id.str.split("_").str[0] == t.chr]
res = t[["gene_base", "chr", "variant_id", "pval_nominal"]].sort_values(["chr", "gene_base", "pval_nominal"])
tmp = out + ".tmp.gz"
res.to_csv(tmp, sep="\t", index=False, compression="gzip")
os.replace(tmp, out)
print(tissue, "eGenes", G.gene_base.nunique(), "with sig rows", res.gene_base.nunique(), "rows", len(res))
