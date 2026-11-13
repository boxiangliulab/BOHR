"""Figure 4 panel e, LD-window coloc, sensitivity (eQTL leads by recursive clumping), step s4 (copy of 05_extract_eqtl.py;
only ROOT changes, to results/coloc_eqtl_ldwin/sens_clump):
cut the GTEx v11 eQTL all-pairs to what the coloc needs.
For one tissue: rows of the genes in results/coloc_eqtl_ldwin/sens_clump/eqtl_windows.tsv (gene matched on the version-stripped
Ensembl ID), variants inside the span of that gene's pair unions (LD windows of the puQTL signal and of each clumped eQTL lead; step s3).
Output: results/coloc_eqtl_ldwin/sens_clump/eqtl_subset/<tissue>.tsv.gz (gene_base, variant_id, pval_nominal, slope, slope_se, af).
Usage: python 01_extract_eqtl.py <tissue>"""
import sys, os
import pandas as pd, pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc

ROOT = "/projects/group/bliu/e1327989/GTEx/QTL_figure/rerun_v2/results/coloc_eqtl_ldwin/sens_clump"
D = "/home/e1327989/project/GTEx/Rawdata/gtex_v11/GTEx_Analysis_v11_eQTL_all_associations/GTEx_Analysis_v11_eQTL_all_associations"
tissue = sys.argv[1]
out = os.path.join(ROOT, "eqtl_subset", tissue + ".tsv.gz")
if os.path.exists(out):
    sys.exit(0)
W = pd.read_csv(os.path.join(ROOT, "eqtl_windows.tsv"), sep="\t")
W = W[W.tissue == tissue]
parts = []
for ch, Wc in W.groupby("chr"):
    t = pq.read_table(os.path.join(D, f"{tissue}.v11.allpairs.{ch}.parquet"),
                      columns=["gene_id", "variant_id", "pval_nominal", "slope", "slope_se", "af"])
    base = pc.list_element(pc.split_pattern(t["gene_id"], pattern="."), 0)
    t = t.append_column("gene_base", base).filter(pc.is_in(base, value_set=pa.array(Wc.gene_base.unique())))
    x = t.drop(["gene_id"]).to_pandas()
    x["pos"] = x.variant_id.str.split("_").str[1].astype("int64")
    x = x.merge(Wc[["gene_base", "start", "end"]], on="gene_base")
    parts.append(x[(x.pos >= x.start) & (x.pos <= x.end)][["gene_base", "variant_id", "pval_nominal", "slope", "slope_se", "af"]])
    del t, x
res = pd.concat(parts, ignore_index=True)
tmp = out + ".tmp.gz"
res.to_csv(tmp, sep="\t", index=False, compression="gzip")
os.replace(tmp, out)
print(tissue, "rows", len(res), "genes", res.gene_base.nunique(), "of", W.gene_base.nunique())
