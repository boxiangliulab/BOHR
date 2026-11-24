#!/usr/bin/env python3
"""187b_extract_encode_zscores_alltested_v3c.py -- ENCODE Z-scores of the TSS elements of 187a. Copy of
Figure5_rerun_v3/scripts/22b2_extract_encode_zscores_alltested.py. Changes, and nothing else:
  1. input accession list = <B1>/data/zscore_accessions_alltested_v3c.txt (from 187a, units v3c) instead of
     Figure5_rerun_v3/data/zscore_accessions_alltested.txt;
  2. outputs = <B1>/data/zscores_<mark>_alltested_v3c.tsv.gz and _columns.tsv (new files; old ones untouched);
  3. the skip test. 22b2 skipped a mark when human_<mark>_core_collection_zscores.tar.gz was missing. Those archives were moved
     to NSCC on 2026-10-08 (encode_v4_zscores/OFFLOADED_TO_NSCC.md), so 22b2 would now silently skip all three marks. This copy
     skips only when the extracted .npy is missing (the .npy, Rows.bed and the columns file are still on the lab disk) and never
     un-tars. The reading of the .npy, the row matching and the output format are those of 22b2."""
import gzip, os, sys
import numpy as np

E = "/home/e1327989/project/GTEx/TATA_Enhancer/encode_v4_zscores"
B1 = "/home/e1327989/project/GTEx/TATA_Enhancer/Figure5_rerun_v4/b1_tss_elements_v3c_20261113"
MARKS = {"h3k4me3": "H3K4me3", "h3k27ac": "H3K27ac", "dnase": "DNase"}
acc = [l.strip() for l in open(os.path.join(B1, "data/zscore_accessions_alltested_v3c.txt")) if l.strip()]
want = set(acc)

for m, M in MARKS.items():
    d = os.path.join(E, m)
    npy = os.path.join(d, f"human_cCRE_{M}_core_collection_zscore.npy")
    out = os.path.join(B1, f"data/zscores_{m}_alltested_v3c.tsv.gz")
    if not os.path.exists(npy) or os.path.exists(out):          # changed test (was: archive missing or output exists)
        print(f"{m}: SKIP (npy missing or output exists)"); sys.exit(1)
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
    with open(os.path.join(B1, f"data/zscores_{m}_alltested_v3c_columns.tsv"), "w") as fo:
        fo.write("col\texperiment_id\tfile_id\tbiosample_value\tassay\tcelltype_value\tcelltype\tlifestage\tsampletype\torgan_tissue\tcollection\n")
        for j, r in enumerate(meta):
            fo.write(f"c{j}\t" + "\t".join(r) + "\n")
    print(f"{m}: matrix {X.shape}, extracted {len(idx)} of {len(want)} accessions, {sub.shape[1]} biosamples; "
          f"value range {np.nanmin(sub):.2f}..{np.nanmax(sub):.2f}, NaN {np.isnan(sub).sum()}")
