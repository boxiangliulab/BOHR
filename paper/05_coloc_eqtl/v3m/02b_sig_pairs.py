"""Mirror rule v3 coloc, step 2b: GTEx v11 significant eQTL variants of the eGenes that need clumping (new_inputs/clump_genes.tsv),
exactly as coloc_eqtl_ldwin/sens_clump/s1_sig_pairs.py (eQTLs.signif_pairs, same chromosome).
Output: results/coloc_eqtl_v3m/new_inputs/sig/<tissue>.tsv.gz (gene_base, chr, variant_id, pval_nominal)."""
import os
import pandas as pd, pyarrow.parquet as pq
O = "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_v3m/new_inputs"
D = "/home/e1327989/project/GTEx/Rawdata/gtex_v11/GTEx_Analysis_v11_eQTL"
CG = pd.read_csv(os.path.join(O, "clump_genes.tsv"), sep="\t")
for tissue, G in CG.groupby("tissue"):
    t = pq.read_table(os.path.join(D, f"{tissue}.v11.eQTLs.signif_pairs.parquet"), columns=["phenotype_id", "variant_id", "pval_nominal"]).to_pandas()
    t["gene_base"] = t.phenotype_id.str.replace(r"\..*$", "", regex=True)
    t = t[t.gene_base.isin(G.gene_base)].merge(G[["gene_base", "chr"]], on="gene_base")
    t = t[t.variant_id.str.split("_").str[0] == t.chr]
    res = t[["gene_base", "chr", "variant_id", "pval_nominal"]].sort_values(["chr", "gene_base", "pval_nominal"])
    res.to_csv(os.path.join(O, "sig", tissue + ".tsv.gz"), sep="\t", index=False, compression="gzip")
    print(tissue, "genes", G.gene_base.nunique(), "with sig rows", res.gene_base.nunique(), "rows", len(res))
