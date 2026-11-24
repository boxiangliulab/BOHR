"""
Scores the 170 variants of the earlier AlphaGenome runs and the one new MRAS candidate of the final version (r5e80)
(chr3_138437434_C_T_b38, per-base run alphagenome_results/perbase_MRAS_r5e80add) on the same fixed COMMON CAGE and
PRO-cap peaks with the formulas of score_v3_candidates.py (sense(+) tracks of the scope summed per peak; signed
log2((ALT + eps)/(REF + eps)), eps = 1e-3; P45112 share = pk3 / (pk3 + pk4)).
Differences from score_v3_candidates.py: a third npz directory; source label 'r5e80_addon'; candidate flag = the 72-variant
set of 20_candidates.R (checks/candidates_r5e80_vs_paper.tsv, status != 'dropped in r5e80'); in_v3_candidate_set is kept.
Checks: all 170 earlier variants reproduce alphagenome_addon/tables/score_*_v3.tsv exactly; CAGE 'all' = 273 sense tracks.
Outputs: this module's alphagenome/tables/. Nothing existing is modified.
"""
import numpy as np, pandas as pd, glob, os, warnings
warnings.filterwarnings('ignore')
ORIG = "/projects/group/bliu/e1327989/alphagenome_results/perbase_MRAS"
ADD  = "/projects/group/bliu/e1327989/alphagenome_results/perbase_MRAS_v3add"
R5   = "/projects/group/bliu/e1327989/alphagenome_results/perbase_MRAS_r5e80add"
TAB  = "/projects/group/bliu/e1327989/GTEx/alphagenome/MRAS_perbase_analysis/tables"
A3   = "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/alphagenome_addon"
MOD  = "/projects/group/bliu/e1327989/GTEx/new_GWAS/coloc_final_nogate_v3/r5e80_20261119/supp_S7_S9"
OUT  = f"{MOD}/alphagenome/tables"; os.makedirs(OUT, exist_ok=True)
NEW  = "chr3_138437434_C_T_b38"
EPS = 1e-3
WET = {'chr3_138348480_G_A_b38': 'v480', 'chr3_138348784_C_T_b38': 'v784', 'chr3_138349284_G_A_b38': 'c284'}

files, src = {}, {}
for d, lab in [(ORIG, 'original95'), (ADD, 'v3_addon'), (R5, 'r5e80_addon')]:
    for f in sorted(glob.glob(f"{d}/per_snp/*.npz")):
        s = os.path.basename(f)[:-4]
        if s not in files: files[s] = f; src[s] = lab
print("variants by source:", pd.Series(src).value_counts().to_dict())
assert NEW in files and src[NEW] == 'r5e80_addon', "new variant not scored yet"
assert len(files) == 171

meta = {m: pd.read_parquet(f"{ORIG}/meta/{m}_tracks.parquet").reset_index(drop=True) for m in ['CAGE', 'PROCAP']}
for d in (ADD, R5):
    for m in ['CAGE', 'PROCAP']:
        ma = pd.read_parquet(f"{d}/meta/{m}_tracks.parquet").reset_index(drop=True)
        assert list(ma['name']) == list(meta[m]['name']), f"{m} track order differs ({d})"
        assert list(ma['strand']) == list(meta[m]['strand']), f"{m} strand differs ({d})"
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
ntr = {(mod, sc): len(scope_idx(mod, sc)) for mod in SCOPES for sc in SCOPES[mod]}
print("sense tracks per scope:", ntr); assert ntr[('CAGE', 'all')] == 273 and ntr[('PROCAP', 'all')] == 6
def sumiv(sig, a, b): w = (pos >= a) & (pos <= b); return float(sig[w].sum())

rows = {'CAGE': [], 'PROCAP': []}; srows = []
pk3 = CAGEpk[CAGEpk.peak_id == 'pk3_P45112+83'].iloc[0]; pk4 = CAGEpk[CAGEpk.peak_id == 'pk4_P45113+13'].iloc[0]
for s, f in sorted(files.items()):
    d = np.load(f); assert np.array_equal(d['pos_CAGE'], pos), s
    base = dict(snp=s, lab=WET.get(s, ''), pos=int(s.split('_')[1]), source=src[s])
    for mod, peaks in [('CAGE', CAGEpk), ('PROCAP', PCP)]:
        for sc in SCOPES[mod]:
            ii = scope_idx(mod, sc)
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

# candidate flags: the earlier set (161, as in the v3 tables) and the set of the final version (72)
paper = set(pd.read_csv(f"{A3}/MRAS_v3_candidate_set_LDstandard_PLINK.tsv", sep="\t").snp)
cr = pd.read_csv(f"{MOD}/checks/candidates_r5e80_vs_paper.tsv", sep="\t")
r5 = set(cr.loc[cr.status != 'dropped in r5e80', 'variant']); assert len(r5) == 72 and NEW in r5
assert r5 <= set(files), "an r5e80 candidate is not scored"
for t in (C, PR, SH):
    t['in_v3_candidate_set'] = t.snp.isin(paper); t['in_original95'] = t.source == 'original95'
    t['in_r5e80_candidate_set'] = t.snp.isin(r5)

# reproduction check: the 170 earlier variants equal the v3 tables in every column
for name, t in [('score_CAGE_peaks', C), ('score_PROCAP_peaks', PR), ('score_CAGE_share', SH)]:
    o = pd.read_csv(f"{A3}/tables/{name}_v3.tsv", sep="\t", keep_default_na=False, na_values=[])
    n = t[t.source != 'r5e80_addon'].drop(columns='in_r5e80_candidate_set').reset_index(drop=True)
    o = o.astype(str); n = n.astype(str)
    assert list(o.columns) == list(n.columns), name
    key = [c for c in ['snp', 'scope', 'peak_id'] if c in o.columns]
    o = o.sort_values(key).reset_index(drop=True); n = n.sort_values(key).reset_index(drop=True)
    ok = o.equals(n); print(f"reproduce {name}_v3.tsv ({len(o)} rows, 170 variants): {'identical' if ok else 'DIFFERENT'}")
    assert ok, name

C.to_csv(f"{OUT}/score_CAGE_peaks_r5e80.tsv", sep="\t", index=False)
PR.to_csv(f"{OUT}/score_PROCAP_peaks_r5e80.tsv", sep="\t", index=False)
SH.to_csv(f"{OUT}/score_CAGE_share_r5e80.tsv", sep="\t", index=False)

piv = C[(C.scope == 'all') & (C.peak_id.isin(['pk3_P45112+83', 'pk4_P45113+13']))].pivot_table(
    index=['snp', 'lab', 'pos', 'source', 'in_r5e80_candidate_set'], columns='peak_id', values='log2fc').reset_index()
piv['maxabs'] = piv[['pk3_P45112+83', 'pk4_P45113+13']].abs().max(axis=1)
piv = piv.merge(SH[SH.scope == 'all'][['snp', 'delta_share112']], on='snp')
sub = piv[piv.in_r5e80_candidate_set].sort_values('maxabs', ascending=False).reset_index(drop=True)
sub['rank_of_72'] = sub.maxabs.rank(ascending=False, method='min').astype(int)
sub.to_csv(f"{OUT}/ranking_r5e80_candidates.tsv", sep="\t", index=False)
print(f"\nr5e80 candidates: {len(sub)}; top 6:"); print(sub.head(6).to_string(index=False))
print("\nnew variant:"); print(sub[sub.snp == NEW].to_string(index=False))
print("\nnew variant, all scopes (CAGE pk3/pk4, share; PRO-cap):")
print(C[(C.snp == NEW) & C.peak_id.isin(['pk3_P45112+83', 'pk4_P45113+13'])][['scope', 'peak_id', 'ref_sum', 'alt_sum', 'log2fc']].to_string(index=False))
print(SH[SH.snp == NEW][['scope', 'delta_share112']].to_string(index=False))
print(PR[PR.snp == NEW][['scope', 'peak_id', 'log2fc']].to_string(index=False))
print("max |log2FC| over all CAGE peaks and scopes, new variant:", C[C.snp == NEW].log2fc.abs().max())
print("done")
