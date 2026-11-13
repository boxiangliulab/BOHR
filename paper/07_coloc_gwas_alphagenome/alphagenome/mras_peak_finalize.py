"""
Finalize the COMMON predicted-CAGE peak set for MRAS = the Adrenal_Gland-called peaks (user choice,
visually most appropriate). Save as the shared peak set used by ALL variants, and draw a merged
confirmation figure (all sense+ mean REF + 4 tracks overlaid, normalized) with the common peaks boxed.
"""
import numpy as np, pandas as pd, glob, os, warnings
warnings.filterwarnings('ignore')
from scipy.signal import find_peaks
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
for fp in glob.glob("/home/e1327989/.local/share/mamba/envs/gtexfig/fonts/*.ttf"):
    try: font_manager.fontManager.addfont(fp)
    except Exception: pass
plt.rcParams.update({"font.family":"Arial","font.size":7,"pdf.fonttype":42,"ps.fonttype":42})

CB="/projects/group/bliu/e1327989/alphagenome_results/perbase_MRAS"
TAB="/projects/group/bliu/e1327989/GTEx/alphagenome/MRAS_perbase_analysis/tables"
FIG="/projects/group/bliu/e1327989/GTEx/alphagenome/MRAS_perbase_analysis/figures"
EX=pd.read_csv(f"{TAB}/mras_gencode_v47_exons.tsv",sep="\t")
TSS={'P45111':138347648,'P45112':138348579,'P45113':138348850,'P45114':138349146}
REGION=(138347300,138349500)
SMOOTH=9; PROM_FRAC=0.008; DIST=12; CLUSTER_GAP=100; EXPAND_FRAC=0.20; FLOOR=3.0

meta=pd.read_parquet(f"{CB}/meta/CAGE_tracks.parquet").reset_index(drop=True)
snps=[os.path.basename(f)[:-4] for f in glob.glob(f"{CB}/per_snp/*.npz")]
pos=np.load(f"{CB}/per_snp/{snps[0]}.npz")['pos_CAGE']
acc=np.zeros((len(pos),len(meta)))
for s in snps: acc+=np.load(f"{CB}/per_snp/{s}.npz")['ref_CAGE']
canon=acc/len(snps); reg=(pos>=REGION[0])&(pos<=REGION[1]); x=pos[reg]
def grp(kind):
    m=meta.strand=='+'
    if kind=='adrenal': m&=meta.biosample_name.astype(str).str.lower().str.contains('adrenal')
    elif kind=='pituitary': m&=meta.biosample_name.astype(str).str.lower().str.contains('pituitar')
    elif kind=='hek293': m&=meta.biosample_name.astype(str)=='HEK293'
    return meta.index[m].tolist()
def near(p): k=min(TSS,key=lambda t:abs(TSS[t]-p)); return k,int(p-TSS[k])
def call_peaks(y):
    ys=np.convolve(y,np.ones(SMOOTH)/SMOOTH,mode='same')
    pk,_=find_peaks(ys,prominence=ys.max()*PROM_FRAC,distance=DIST,height=FLOOR); pk=np.sort(pk)
    cl=[[pk[0]]]
    for s in pk[1:]:
        (cl[-1].append(s) if x[s]-x[cl[-1][-1]]<=CLUSTER_GAP else cl.append([s]))
    out=[]
    for ci,c in enumerate(cl):
        lo_s,hi_s=c[0],c[-1]; thr=ys[c].max()*EXPAND_FRAC
        lcap=0 if ci==0 else int((cl[ci-1][-1]+lo_s)//2); rcap=len(ys)-1 if ci==len(cl)-1 else int((hi_s+cl[ci+1][0])//2)
        l=lo_s; r=hi_s
        while l>lcap and ys[l-1]>=thr: l-=1
        while r<rcap and ys[r+1]>=thr: r+=1
        smt=c[int(np.argmax(ys[c]))]; nm,dd=near(int(x[smt]))
        out.append(dict(start=int(x[l]),end=int(x[r]),summit=int(x[smt]),nearest_prom=nm,dist=dd,nsub=len(c)))
    return out,ys

# COMMON peak set = Adrenal
adr=canon[:,grp('adrenal')].mean(1)[reg]
common,_=call_peaks(adr)
CP=pd.DataFrame(common)
CP['peak_id']=[f"pk{i+1}_{r.nearest_prom}{r.dist:+d}" for i,r in CP.iterrows()]
CP=CP[['peak_id','start','end','summit','nearest_prom','dist','nsub']]
CP.to_csv(f"{TAB}/predicted_cage_peaks_COMMON.tsv",sep="\t",index=False)
print("COMMON peak set (from Adrenal):"); print(CP.to_string(index=False))

# ---- merged confirmation figure ----
fig,axes=plt.subplots(2,1,figsize=(9,4.2),sharex=True,gridspec_kw={'height_ratios':[1,0.45]})
ax=axes[0]
allmean=canon[:,grp('all')].mean(1)[reg]
allmean_n=allmean/allmean.max()
ax.fill_between(x,0,allmean_n,color="#c9d7ec",lw=0,step='mid',label="all sense+ (merged) REF",zorder=1)
for kind,c in [('adrenal','#d73027'),('pituitary','#1a9850'),('hek293','#7b3294')]:
    yy=canon[:,grp(kind)].mean(1)[reg]; ax.plot(x,yy/yy.max(),color=c,lw=0.8,label=kind,zorder=3)
for _,pkd in CP.iterrows():
    ax.axvspan(pkd.start,pkd.end,color="#F4A582",alpha=0.35,zorder=0)
    ax.text((pkd.start+pkd.end)/2,1.02,f"{pkd.peak_id}\n{pkd.start}-{pkd.end}",fontsize=5,ha='center',va='bottom',color='#7b241c')
for nm,t in TSS.items(): ax.axvline(t,color='#3366cc',ls='--',lw=0.6)
ax.set_ylabel("REF CAGE (normalized to own max)",fontsize=7); ax.set_ylim(0,1.18)
ax.legend(fontsize=6,frameon=False,ncol=4,loc='upper right')
ax.set_title("MRAS common predicted-CAGE peak set (from Adrenal) — confirmation on merged + per-track REF",fontsize=8)
axe=axes[1]; yv=0
for _,t in EX[EX.feature=='transcript'].sort_values('start').iterrows():
    if t.end<REGION[0] or t.start>REGION[1]: continue
    exs=EX[(EX.feature=='exon')&(EX.tx==t.tx)]; exs=exs[(exs.end>=REGION[0])&(exs.start<=REGION[1])]
    if not len(exs): continue
    axe.plot([max(t.start,REGION[0]),min(t.end,REGION[1])],[yv,yv],color='0.5',lw=0.5)
    for _,e in exs.iterrows(): axe.add_patch(plt.Rectangle((e.start,yv-0.35),max(e.end-e.start,6),0.7,color='#2c3e50'))
    yv-=1
for _,pkd in CP.iterrows(): axe.axvspan(pkd.start,pkd.end,color="#F4A582",alpha=0.35)
axe.set_ylim(yv+0.4,0.8); axe.set_yticks([]); axe.set_ylabel("MRAS\nexons",fontsize=6)
for nm,t in TSS.items(): axe.axvline(t,color='#3366cc',ls='--',lw=0.5); axe.text(t,0.9,nm,rotation=90,fontsize=5,color='#3366cc',ha='center',va='bottom')
axe.set_xlim(*REGION); axe.set_xlabel("chr3 position (orange = COMMON peak; blue dashed = promoter TSS)",fontsize=7)
fig.tight_layout()
for e in ("pdf","png"): fig.savefig(f"{FIG}/fig_peak_common.{e}",dpi=300,bbox_inches="tight")
print("wrote fig_peak_common + predicted_cage_peaks_COMMON.tsv")
