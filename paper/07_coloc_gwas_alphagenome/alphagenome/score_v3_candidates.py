"""
Score the v3 MRAS candidate set (Figure 6, 2026-10-03) on the FIXED promoter features used for the original 95
variants: the same COMMON CAGE peak set (predicted_cage_peaks_COMMON.tsv, Adrenal-defined) and the same COMMON PRO-cap
peak set (predicted_procap_peaks_COMMON.tsv), read from the original tables and NOT re-called. The scoring is identical
to MRAS_perbase_analysis/scripts/mras_score_variants.py:
- per peak, sum the sense(+) tracks of the scope;
- signed log2((ALT + eps)/(REF + eps)), eps = 1e-3;
- P45112 share = pk3 / (pk3 + pk4).
Variants: the original 95 (alphagenome_results/perbase_MRAS) plus the 75 add-ons (alphagenome_results/perbase_MRAS_v3add).
Sanity check: the recomputed values for the 95 must equal the original score_CAGE_peaks.tsv.
Outputs go to coloc_final_nogate_v3/alphagenome_addon/tables/. Nothing in the original tables is modified.
"""
import numpy as np, pandas as pd, glob, os, sys, warnings
warnings.filterwarnings('ignore')
ORIG = "/projects/group/bliu/e1327989/alphagenome_results/perbase_MRAS"
ADD  = "/projects/group/bliu/e1327989/alphagenome_results/perbase_MRAS_v3add"
TAB  = "/projects/group/bliu/e1327989/GTEx/alphagenome/MRAS_perbase_analysis/tables"
A3   = "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/alphagenome_addon"
OUT  = f"{A3}/tables"; os.makedirs(OUT, exist_ok=True)
EPS = 1e-3; REGION = (138347300, 138349500)
WET = {'chr3_138348480_G_A_b38': 'v480', 'chr3_138348784_C_T_b38': 'v784', 'chr3_138349284_G_A_b38': 'c284'}

files = {os.path.basename(f)[:-4]: f for f in glob.glob(f"{ORIG}/per_snp/*.npz")}
n_orig = len(files)
ORIG_SET = set(files)   # source is decided by membership, NOT by path prefix ('perbase_MRAS_v3add' starts with 'perbase_MRAS')
for f in glob.glob(f"{ADD}/per_snp/*.npz"):
    s = os.path.basename(f)[:-4]
    if s not in files: files[s] = f
print(f"variants: original {n_orig}, add-on {len(files) - n_orig}, total {len(files)}")
cand = pd.read_csv(f"{A3}/MRAS_v3_candidate_set_LDstandard_PLINK.tsv", sep="\t")
todo = set(pd.read_csv(f"{A3}/MRAS_v3_new_candidates_to_score_ALL75.txt", header=None)[0])
missing = sorted(todo - set(files)); print(f"add-on variants still missing: {len(missing)}", missing[:10])

meta = {m: pd.read_parquet(f"{ORIG}/meta/{m}_tracks.parquet").reset_index(drop=True) for m in ['CAGE', 'PROCAP']}
for m in ['CAGE', 'PROCAP']:
    mp = f"{ADD}/meta/{m}_tracks.parquet"
    if os.path.exists(mp):
        ma = pd.read_parquet(mp).reset_index(drop=True)
        assert list(ma['name']) == list(meta[m]['name']), f"{m} track order differs between runs"
pos = np.load(files[next(iter(files))])['pos_CAGE']
CAGEpk = pd.read_csv(f"{TAB}/predicted_cage_peaks_COMMON.tsv", sep="\t")
PCP = pd.read_csv(f"{TAB}/predicted_procap_peaks_COMMON.tsv", sep="\t")

def scope_idx(mod, scope):
    m = meta[mod]; msk = (m.strand == '+'); bn = m.biosample_name.astype(str)
    if scope == 'hpa': msk = msk & bn.str.lower().str.contains('pituitar|adrenal')
    elif scope == 'hek293': msk = msk & (bn == 'HEK293')
    elif scope == 'cellline': msk = msk & (m.biosample_type.astype(str) == 'cell_line')
    elif scope == 'K562': msk = msk & (bn == 'K562')
    elif scope == 'A673': msk = msk & (bn == 'A673')
    return m.index[msk].tolist()
SCOPES = {'CAGE': ['all', 'hpa', 'hek293', 'cellline'], 'PROCAP': ['all', 'K562', 'A673']}
def sumiv(sig, a, b): w = (pos >= a) & (pos <= b); return float(sig[w].sum())

rows = {'CAGE': [], 'PROCAP': []}; srows = []
pk3 = CAGEpk[CAGEpk.peak_id == 'pk3_P45112+83'].iloc[0]; pk4 = CAGEpk[CAGEpk.peak_id == 'pk4_P45113+13'].iloc[0]
for s, f in sorted(files.items()):
    d = np.load(f); assert np.array_equal(d['pos_CAGE'], pos), s
    base = dict(snp=s, lab=WET.get(s, ''), pos=int(s.split('_')[1]), source='original95' if s in ORIG_SET else 'v3_addon')
    for mod, peaks in [('CAGE', CAGEpk), ('PROCAP', PCP)]:
        for sc in SCOPES[mod]:
            ii = scope_idx(mod, sc)
            if not ii: continue
            ref = d[f'ref_{mod}'][:, ii].sum(1); alt = d[f'alt_{mod}'][:, ii].sum(1)
            for _, pk in peaks.iterrows():
                R = sumiv(ref, pk.start, pk.end); A = sumiv(alt, pk.start, pk.end)
                rows[mod].append({**base, 'scope': sc, 'peak_id': pk.peak_id, 'nearest_prom': pk.nearest_prom,
                                  'ref_sum': round(R, 3), 'alt_sum': round(A, 3), 'diff': round(A - R, 3),
                                  'log2fc': round(float(np.log2((A + EPS) / (R + EPS))), 4)})
            if mod == 'CAGE':
                R3 = sumiv(ref, pk3.start, pk3.end); A3_ = sumiv(alt, pk3.start, pk3.end)
                R4 = sumiv(ref, pk4.start, pk4.end); A4 = sumiv(alt, pk4.start, pk4.end)
                rs = R3 / (R3 + R4 + EPS); as_ = A3_ / (A3_ + A4 + EPS)
                srows.append({**base, 'scope': sc, 'ref_share112': round(rs, 4), 'alt_share112': round(as_, 4), 'delta_share112': round(as_ - rs, 4)})
C = pd.DataFrame(rows['CAGE']); PR = pd.DataFrame(rows['PROCAP']); SH = pd.DataFrame(srows)
for t in (C, PR, SH):
    t['in_v3_candidate_set'] = t.snp.isin(set(cand.snp)); t['in_original95'] = t.source == 'original95'
C.to_csv(f"{OUT}/score_CAGE_peaks_v3.tsv", sep="\t", index=False)
PR.to_csv(f"{OUT}/score_PROCAP_peaks_v3.tsv", sep="\t", index=False)
SH.to_csv(f"{OUT}/score_CAGE_share_v3.tsv", sep="\t", index=False)

# sanity check against the original table (95 variants)
o = pd.read_csv(f"{TAB}/score_CAGE_peaks.tsv", sep="\t")
mrg = o.merge(C[C.source == 'original95'], on=['snp', 'scope', 'peak_id'], suffixes=('_orig', '_new'))
dmax = float((mrg.log2fc_orig - mrg.log2fc_new).abs().max())
print(f"sanity: {len(mrg)} rows matched to the original table; max |log2fc diff| = {dmax:.2e}")

# ranking (scope all): max |log2fc| over pk3/pk4, the original metric
piv = C[(C.scope == 'all') & (C.peak_id.isin(['pk3_P45112+83', 'pk4_P45113+13']))].pivot_table(
    index=['snp', 'lab', 'pos', 'source', 'in_v3_candidate_set'], columns='peak_id', values='log2fc').reset_index()
piv['maxabs'] = piv[['pk3_P45112+83', 'pk4_P45113+13']].abs().max(axis=1)
piv = piv.merge(SH[SH.scope == 'all'][['snp', 'delta_share112']], on='snp')
for scope_name, sub in [('all scored (95 + add-on)', piv), ('v3 candidate set (LD standard)', piv[piv.in_v3_candidate_set])]:
    sub = sub.sort_values('maxabs', ascending=False).reset_index(drop=True)
    r = sub.index[sub.lab == 'v784']
    print(f"\n=== {scope_name}: n = {len(sub)}; v784 rank by max|log2FC| = {int(r[0]) + 1 if len(r) else 'NA'} ===")
    print(sub.head(8)[['snp', 'lab', 'source', 'pk3_P45112+83', 'pk4_P45113+13', 'maxabs', 'delta_share112']].to_string(index=False))
piv.sort_values('maxabs', ascending=False).to_csv(f"{OUT}/ranking_v3_candidates.tsv", sep="\t", index=False)
print(f"\nsource counts: {piv.source.value_counts().to_dict()}")
print(f"add-on variants with |log2FC| > 0.02 at pk3/pk4: {int(((piv.source == 'v3_addon') & (piv.maxabs > 0.02)).sum())}")
print(piv[(piv.source == 'v3_addon')].sort_values('maxabs', ascending=False).head(8)[['snp', 'pk3_P45112+83', 'pk4_P45113+13', 'maxabs', 'delta_share112']].to_string(index=False))
