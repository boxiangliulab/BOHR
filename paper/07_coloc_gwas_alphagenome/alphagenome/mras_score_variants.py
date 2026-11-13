"""
MASTER variant scoring on FIXED promoter features (the correct, promoter-focused analysis).
 - CAGE  : per-peak on COMMON peak set (Adrenal-defined, predicted_cage_peaks_COMMON.tsv);
           scopes all/HPA/HEK293/cellline; + promoter-share pk3(P45112) vs pk4(P45113).
 - PROCAP: per-peak on its OWN common peak set (called from PROCAP all sense+ REF); scopes all/K562/A673.
 - ATAC / DNASE : WHOLE promoter region [138347300,138349500] (broad accessibility, peaks not separable);
           scopes all / adrenal; REF_SUM/ALT_SUM/DIFF/log2FC only.
All from promoter-fixed per-base windows (same window for every variant). Sense(+) for CAGE/PROCAP.
Outputs tables/score_*.tsv.
"""
import numpy as np, pandas as pd, glob, os, warnings
warnings.filterwarnings('ignore')
from scipy.signal import find_peaks
CB="/projects/group/bliu/e1327989/alphagenome_results/perbase_MRAS"
TAB="/projects/group/bliu/e1327989/GTEx/alphagenome/MRAS_perbase_analysis/tables"
TSS={'P45111':138347648,'P45112':138348579,'P45113':138348850,'P45114':138349146}
EPS=1e-3; REGION=(138347300,138349500)
WET={'chr3_138348480_G_A_b38':'v480','chr3_138348784_C_T_b38':'v784','chr3_138349284_G_A_b38':'c284'}
snps=[os.path.basename(f)[:-4] for f in glob.glob(f"{CB}/per_snp/*.npz")]
pos=np.load(f"{CB}/per_snp/{snps[0]}.npz")['pos_CAGE']
# leader info
BASE="/projects/group/bliu/e1327989/GTEx/TATA_Enhancer/new_lrt_fixed_analysis/selected_cutoffs/mu0.70_p0.35"
P=pd.concat([pd.read_csv(f"{BASE}/qtl_{t}_proxy.txt",sep="\t",usecols=lambda c:c in ['promoter','snp','is_leader']) for t in ['mu','p0','p1']])
P=P[P.promoter.str.startswith('PromoterID_4511')]
leadmap=P[P.is_leader==True].groupby('snp').promoter.apply(lambda x:','.join(sorted(set(p.replace('PromoterID_','P') for p in x)))).to_dict()

def meta_of(mod): return pd.read_parquet(f"{CB}/meta/{mod}_tracks.parquet").reset_index(drop=True)
def near(p): k=min(TSS,key=lambda t:abs(TSS[t]-p)); return k,int(p-TSS[k])

# ---- PROCAP own common peaks (call on all sense+ REF mean) ----
def canon_ref(mod):
    acc=np.zeros((len(pos),len(meta_of(mod))))
    for s in snps: acc+=np.load(f"{CB}/per_snp/{s}.npz")[f'ref_{mod}']
    return acc/len(snps)
def call_peaks(y):
    reg=(pos>=REGION[0])&(pos<=REGION[1]); x=pos[reg]; ys=np.convolve(y[reg],np.ones(9)/9,mode='same')
    if ys.max()<=0: return []
    pk,_=find_peaks(ys,prominence=ys.max()*0.008,distance=12,height=ys.max()*0.05)
    if not len(pk): return []
    pk=np.sort(pk); cl=[[pk[0]]]
    for s in pk[1:]:
        (cl[-1].append(s) if x[s]-x[cl[-1][-1]]<=100 else cl.append([s]))
    out=[]
    for ci,c in enumerate(cl):
        lo_s,hi_s=c[0],c[-1]; thr=ys[c].max()*0.20
        lcap=0 if ci==0 else int((cl[ci-1][-1]+lo_s)//2); rcap=len(ys)-1 if ci==len(cl)-1 else int((hi_s+cl[ci+1][0])//2)
        l=lo_s; r=hi_s
        while l>lcap and ys[l-1]>=thr: l-=1
        while r<rcap and ys[r+1]>=thr: r+=1
        smt=c[int(np.argmax(ys[c]))]; nm,dd=near(int(x[smt]))
        out.append(dict(start=int(x[l]),end=int(x[r]),summit=int(x[smt]),nearest_prom=nm,dist=dd))
    return out

pm=meta_of('PROCAP'); pref=canon_ref('PROCAP')[:,pm.index[pm.strand=='+']].mean(1)
PCP=pd.DataFrame(call_peaks(pref))
PCP['peak_id']=[f"pp{i+1}_{r.nearest_prom}{r.dist:+d}" for i,r in PCP.iterrows()]
PCP.to_csv(f"{TAB}/predicted_procap_peaks_COMMON.tsv",sep="\t",index=False)
CAGEpk=pd.read_csv(f"{TAB}/predicted_cage_peaks_COMMON.tsv",sep="\t")
print("PROCAP common peaks:\n",PCP[['peak_id','start','end','summit','nearest_prom','dist']].to_string(index=False))

# ---- scope track indices ----
def scope_idx(mod,scope):
    m=meta_of(mod); msk=(m.strand=='+') if scope not in ('atac_all','dnase_all','region_all') and 'strand' in m.columns and (m.strand=='+').any() else pd.Series(True,index=m.index)
    bn=m.biosample_name.astype(str)
    if scope=='hpa': msk=msk&bn.str.lower().str.contains('pituitar|adrenal')
    elif scope=='hek293': msk=msk&(bn=='HEK293')
    elif scope=='cellline': msk=msk&(m.biosample_type.astype(str)=='cell_line')
    elif scope=='K562': msk=msk&(bn=='K562')
    elif scope=='A673': msk=msk&(bn=='A673')
    elif scope=='adrenal': msk=msk&bn.str.lower().str.contains('adrenal')
    return m.index[msk].tolist()

SCOPES={'CAGE':['all','hpa','hek293','cellline'],'PROCAP':['all','K562','A673'],
        'ATAC':['all','adrenal'],'DNASE':['all','adrenal']}
def sumiv(sig,a,b): w=(pos>=a)&(pos<=b); return float(sig[w].sum())

# ---- peak-level: CAGE + PROCAP ----
for mod,peaks in [('CAGE',CAGEpk),('PROCAP',PCP)]:
    idxc={sc:scope_idx(mod,sc) for sc in SCOPES[mod]}
    rows=[]
    for s in snps:
        d=np.load(f"{CB}/per_snp/{s}.npz"); base=dict(snp=s,lab=WET.get(s,''),pos=int(s.split('_')[1]),leads=leadmap.get(s,''))
        for sc in SCOPES[mod]:
            ii=idxc[sc]
            if not ii: continue
            ref=d[f'ref_{mod}'][:,ii].sum(1); alt=d[f'alt_{mod}'][:,ii].sum(1)
            for _,pk in peaks.iterrows():
                R=sumiv(ref,pk.start,pk.end); A=sumiv(alt,pk.start,pk.end)
                rows.append({**base,'scope':sc,'peak_id':pk.peak_id,'nearest_prom':pk.nearest_prom,
                             'ref_sum':round(R,3),'alt_sum':round(A,3),'diff':round(A-R,3),
                             'log2fc':round(float(np.log2((A+EPS)/(R+EPS))),4),'active':round(max(R,A),3)})
    pd.DataFrame(rows).to_csv(f"{TAB}/score_{mod}_peaks.tsv",sep="\t",index=False)
    print(f"wrote score_{mod}_peaks.tsv ({len(rows)} rows)")

# ---- CAGE promoter-share pk3(112) vs pk4(113) ----
pk3=CAGEpk[CAGEpk.peak_id.str.contains('P45112')].sort_values('nsub' if 'nsub' in CAGEpk else 'start').iloc[-1] if (CAGEpk.peak_id.str.contains('P45112')).any() else None
pk3=CAGEpk[CAGEpk.peak_id=='pk3_P45112+83'].iloc[0]; pk4=CAGEpk[CAGEpk.peak_id=='pk4_P45113+13'].iloc[0]
idxc={sc:scope_idx('CAGE',sc) for sc in SCOPES['CAGE']}
srows=[]
for s in snps:
    d=np.load(f"{CB}/per_snp/{s}.npz"); base=dict(snp=s,lab=WET.get(s,''),pos=int(s.split('_')[1]),leads=leadmap.get(s,''))
    for sc in SCOPES['CAGE']:
        ii=idxc[sc]
        if not ii: continue
        ref=d['ref_CAGE'][:,ii].sum(1); alt=d['alt_CAGE'][:,ii].sum(1)
        R3=sumiv(ref,pk3.start,pk3.end); A3=sumiv(alt,pk3.start,pk3.end)
        R4=sumiv(ref,pk4.start,pk4.end); A4=sumiv(alt,pk4.start,pk4.end)
        rs=R3/(R3+R4+EPS); as_=A3/(A3+A4+EPS)
        srows.append({**base,'scope':sc,'ref_share112':round(rs,4),'alt_share112':round(as_,4),
                      'delta_share112':round(as_-rs,4)})
SH=pd.DataFrame(srows); SH.to_csv(f"{TAB}/score_CAGE_share.tsv",sep="\t",index=False)
print("wrote score_CAGE_share.tsv")

# ---- ATAC/DNASE whole-region ----
for mod in ['ATAC','DNASE']:
    idxc={sc:scope_idx(mod,sc) for sc in SCOPES[mod]}
    rows=[]
    for s in snps:
        d=np.load(f"{CB}/per_snp/{s}.npz"); base=dict(snp=s,lab=WET.get(s,''),pos=int(s.split('_')[1]),leads=leadmap.get(s,''))
        for sc in SCOPES[mod]:
            ii=idxc[sc]
            if not ii: continue
            ref=d[f'ref_{mod}'][:,ii].sum(1); alt=d[f'alt_{mod}'][:,ii].sum(1)
            R=sumiv(ref,*REGION); A=sumiv(alt,*REGION)
            rows.append({**base,'scope':sc,'ref_sum':round(R,3),'alt_sum':round(A,3),'diff':round(A-R,3),
                         'log2fc':round(float(np.log2((A+EPS)/(R+EPS))),4)})
    pd.DataFrame(rows).to_csv(f"{TAB}/score_{mod}_region.tsv",sep="\t",index=False)
    print(f"wrote score_{mod}_region.tsv")

# ---- sanity: top CAGE 112/113 effects (all-tissue) + share ----
c=pd.read_csv(f"{TAB}/score_CAGE_peaks.tsv",sep="\t")
main=c[(c.scope=='all')&(c.peak_id.isin(['pk3_P45112+83','pk4_P45113+13']))].pivot_table(index=['snp','lab','pos','leads'],columns='peak_id',values='log2fc').reset_index()
main['maxabs']=main[['pk3_P45112+83','pk4_P45113+13']].abs().max(axis=1)
print("\n=== CAGE all-tissue: top10 by |effect on 112/113 peaks| ===")
print(main.sort_values('maxabs',ascending=False).head(10).to_string(index=False))
print("\n=== v480/v784/c284 CAGE share (all-tissue) ===")
print(SH[(SH.scope=='all')&(SH.lab!='')].to_string(index=False))
