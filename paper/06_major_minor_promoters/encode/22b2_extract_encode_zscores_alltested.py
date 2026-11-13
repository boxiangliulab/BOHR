#!/usr/bin/env python3
"""22b2_extract_encode_zscores_alltested.py -- (copy of 22b for the all-tested-promoter TSS elements of 25a) pull ENCODE v4 core-collection Z-scores (H3K4me3, H3K27ac, DNase) for the
element sets of 22a (data/zscore_accessions_alltested.txt) out of the SCREEN matrices (TATA_Enhancer/encode_v4_zscores/).
Each archive holds: human_cCRE_<MARK>_core_collection_zscore.npy (rows = GRCh38-cCREs.Rows.bed order, columns =
<mark>-biosamples-metadata-core-collection-columns.txt order). The .npy is extracted once and read with a memory map.
Output per mark: data/zscores_<mark>.tsv.gz (acc + one column per biosample) and data/zscores_<mark>_columns.tsv."""
import gzip, os, subprocess, sys
import numpy as np

E = "/home/e1327989/project/GTEx/TATA_Enhancer/encode_v4_zscores"
W = "/home/e1327989/project/GTEx/TATA_Enhancer/Figure5_rerun_v3"
MARKS = {"h3k4me3": "H3K4me3", "h3k27ac": "H3K27ac", "dnase": "DNase"}
acc = [l.strip() for l in open(os.path.join(W, "data/zscore_accessions_alltested.txt")) if l.strip()]
want = set(acc)

for m, M in MARKS.items():
    tgz = os.path.join(E, f"human_{m}_core_collection_zscores.tar.gz")
    out = os.path.join(W, f"data/zscores_{m}_alltested.tsv.gz")
    if not os.path.exists(tgz) or os.path.exists(out):
        print(f"{m}: skip (archive missing or output exists)"); continue
    d = os.path.join(E, m); os.makedirs(d, exist_ok=True)
    npy = os.path.join(d, f"human_cCRE_{M}_core_collection_zscore.npy")
    if not os.path.exists(npy):
        subprocess.run(["tar", "-xzf", tgz, "-C", d], check=True)
    rows = os.path.join(d, "GRCh38-cCREs.Rows.bed")
    cols = [f for f in os.listdir(d) if f.endswith("-biosamples-metadata-core-collection-columns.txt")][0]
    idx, names = [], []
    with open(rows) as fh:
        for i, line in enumerate(fh):
            a = line.split("\t")[4]
            if a in want:
                idx.append(i); names.append(a)
    X = np.load(npy, mmap_mode="r")
    meta = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, cols))]
    assert X.shape[1] == len(meta), (X.shape, len(meta))
    sub = np.asarray(X[np.array(idx), :], dtype=np.float32)
    with gzip.open(out, "wt") as fo:
        fo.write("acc\t" + "\t".join(f"c{j}" for j in range(sub.shape[1])) + "\n")
        for a, r in zip(names, sub):
            fo.write(a + "\t" + "\t".join(f"{v:.3f}" for v in r) + "\n")
    with open(os.path.join(W, f"data/zscores_{m}_alltested_columns.tsv"), "w") as fo:
        fo.write("col\texperiment_id\tfile_id\tbiosample_value\tassay\tcelltype_value\tcelltype\tlifestage\tsampletype\torgan_tissue\tcollection\n")
        for j, r in enumerate(meta):
            fo.write(f"c{j}\t" + "\t".join(r) + "\n")
    print(f"{m}: matrix {X.shape}, extracted {len(idx)} of {len(want)} accessions, {sub.shape[1]} biosamples; "
          f"value range {np.nanmin(sub):.2f}..{np.nanmax(sub):.2f}, NaN {np.isnan(sub).sum()}")
