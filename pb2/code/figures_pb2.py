#!/usr/bin/env python3
"""PB2 figures: manhattan, host-group heatmap, 3D structure, temporal."""
import csv, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser

DOMAINS={'PB2':[(1,110,'PB1-binding'),(250,500,'cap-binding region'),(530,759,'627/NLS domain')]}
def load(p): return list(csv.DictReader(open(p),delimiter='\t'))

rows=load('results/PB2.enrichment.tsv')
pos=np.array([int(r['pos']) for r in rows]); q=np.array([max(float(r['q']),1e-300) for r in rows])
y=-np.log10(q); fm=np.array([float(r['freq_m']) for r in rows])
g2=np.array([r['g2_pass']=='True' for r in rows])
fig,ax=plt.subplots(figsize=(13,4.2))
sizes=4+300*fm
ax.scatter(pos[~g2],y[~g2],s=sizes[~g2],c='#7aa6c2',alpha=0.45,lw=0,label='tested sites')
ax.scatter(pos[g2],y[g2],s=sizes[g2]+10,c='#c22f2f',alpha=0.9,lw=0,label='G2 candidates')
for p0,p1,name in DOMAINS['PB2']:
    ax.axvspan(p0,p1,color='k',alpha=0.03)
    ax.text((p0+p1)/2,2,name,ha='center',fontsize=7,rotation=90,color='gray')
for mk,lab in [(627,'E627K'),(701,'D701N'),(591,'Q591K'),(271,'T271A')]:
    r=next((r for r in rows if int(r['pos'])==mk),None)
    if r and lab.endswith(r['alt']):
        yv=-np.log10(max(float(r['q']),1e-300))
        ax.annotate(lab,(mk,yv),textcoords='offset points',xytext=(6,6),fontsize=8,color='#114',
                    arrowprops=dict(arrowstyle='-',lw=0.4,color='#889'))
ax.axhline(-np.log10(0.05),ls='--',c='k',lw=0.7,alpha=0.6)
ax.set_xlabel('PB2 position (Gs/Gd/96 numbering)'); ax.set_ylabel('-log10(FDR q)')
ax.set_title('PB2: mammalian-enriched substitutions, H5N1 2021-2026 (point size = mammalian frequency)')
ax.legend(loc='upper left',fontsize=7,frameon=False)
plt.tight_layout(); plt.savefig('results/fig_manhattan_PB2.png',dpi=180); plt.close()
print('manhattan done')

gf=load('results/PB2.groupfreq.tsv')
groups=sorted(set(r['group'] for r in gf))
sites=sorted(set((int(r['pos']),r['alt']) for r in gf))
M=np.zeros((len(sites),len(groups)))
for r in gf: M[sites.index((int(r['pos']),r['alt']))][groups.index(r['group'])]=float(r['freq'])
fig,ax=plt.subplots(figsize=(7,max(2.2,0.35*len(sites)+1.5)))
im=ax.imshow(M,aspect='auto',cmap='Reds',vmin=0,vmax=max(0.2,M.max()))
ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups,rotation=35,ha='right',fontsize=8)
ax.set_yticks(range(len(sites))); ax.set_yticklabels([f'{p}{a}' for p,a in sites],fontsize=8)
for i in range(len(sites)):
    for j in range(len(groups)):
        ax.text(j,i,f'{M[i,j]:.2f}',ha='center',va='center',fontsize=6.5,
                color='white' if M[i,j]>0.6*max(0.2,M.max()) else 'black')
ax.set_title('PB2 candidate frequencies by mammalian host group')
plt.colorbar(im,label='frequency',shrink=0.8)
plt.tight_layout(); plt.savefig('results/fig_heatmap_PB2.png',dpi=180); plt.close()
print('heatmap done')

# temporal: mammalian isolates and E627K frequency by year
man=load('parsed/manifest_seg1.tsv')
years={}
for r in man:
    if r['year'] and r['gene']=='PB2':
        y_=r['year']; years.setdefault(y_,[0,0]); years[y_][0]+=1
cand627={r['accession'] for r in [] }  # placeholder
# E627K carriers: need sequences; use enrichment m counts per year approximation via manifest is not possible; skip per-year marker, plot counts only
ys=sorted(years); tot=[years[y_][0] for y_ in ys]
fig,ax=plt.subplots(figsize=(7,3))
ax.bar(ys,tot,color='#4a80b4')
ax.set_title('PB2 segment records by collection year (H5N1, all hosts)')
ax.set_xlabel('year'); ax.set_ylabel('records')
plt.tight_layout(); plt.savefig('results/fig_temporal_PB2.png',dpi=180); plt.close()
print('temporal done')

# 3D
m=MMCIFParser(QUIET=True).get_structure('x','structures/8R1J.cif')[0]
def ca(ch):
    return {r.id[1]:r['CA'].coord for r in m[ch] if r.id[0]==' ' and 'CA' in r}
C={'PB1':ca('B'),'PA':ca('A'),'PB2':ca('C'),'ANP32B':ca('G')}
cands=[r for r in load('results/cands_PB2.tsv')]
fig=plt.figure(figsize=(13,5.6))
for vi,(elev,azim,ttl) in enumerate([(20,-60,'view 1'),(80,-60,'top view')]):
    ax=fig.add_subplot(1,2,vi+1,projection='3d')
    for name,col,alpha in [('PB1','#4a80b4',0.25),('PA','#54a06a',0.25),('PB2','#999999',0.2),('ANP32B','#e08a3c',0.8)]:
        pts=np.array(list(C[name].values()))
        ax.scatter(pts[:,0],pts[:,1],pts[:,2],s=3,c=col,alpha=alpha,lw=0,label=name)
    for r in cands:
        p0=int(r['pos'])
        if p0 in C['PB2']:
            p=C['PB2'][p0]
            ax.scatter([p[0]],[p[1]],[p[2]],s=60,c='#c22f2f',edgecolor='k',lw=0.4,depthshade=False)
            ax.text(p[0],p[1],p[2],f" {r['avian_cons']}{p0}{r['alt']}",fontsize=6.5)
    ax.view_init(elev=elev,azim=azim); ax.set_title(ttl,fontsize=8)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    if vi==0: ax.legend(fontsize=7,frameon=False,loc='upper left')
fig.suptitle('H5N1 polymerase + human ANP32B (PDB 8R1J, CA trace) with PB2 mammalian-enriched candidates',fontsize=9)
plt.tight_layout(); plt.savefig('results/fig_structure3d.png',dpi=180); plt.close()
print('structure3d done')
