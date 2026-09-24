#!/usr/bin/env python3
"""Map PB2 candidate residues onto 8R1J (H5N1 polymerase dimer + human ANP32B).
8R1J: copy1 PA=A/PB1=B/PB2=C, copy2 PA=D/PB1=E/PB2=F, ANP32B=G.
8R1L: ANP32B=D, PB2=C (partial 41-248), PA=A(202-716), PB1=B.
Reports min CA-CA distance to ANP32B and to PB1 catalytic Asp anchors (305,306,445-447),
plus distance to PB2 residue 627 (host-adaptation hub). Unresolved -> THIN.
Usage: structure_map_pb2.py <cands_tsv> <out_tsv>"""
import sys, csv, warnings
warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
import numpy as np

ANCHOR_POS=[305,306,445,446,447]
def chain_ca(model,ch):
    out={}
    if ch not in model: return out
    for res in model[ch]:
        if res.id[0]!=' ': continue
        if 'CA' in res: out[res.id[1]]=res['CA'].coord
    return out
def min_dist(p,cmap):
    if not cmap: return np.nan
    arr=np.array(list(cmap.values()))
    return float(np.sqrt(((arr-p)**2).sum(axis=1)).min())

def main(cands,out_tsv):
    m1=MMCIFParser(QUIET=True).get_structure('j','structures/8R1J.cif')[0]
    ca1={c:chain_ca(m1,c) for c in 'ABCDEFG'}
    m2=MMCIFParser(QUIET=True).get_structure('l','structures/8R1L.cif')[0]
    ca2={c:chain_ca(m2,c) for c in 'ABCD'}
    rows=[]
    for r in csv.DictReader(open(cands),delimiter='\t'):
        pos=int(r['pos'])
        hit=None
        for ch in ('C','F'):
            if pos in ca1[ch]: hit=('8R1J',ch); break
        if not hit and pos in ca2['C']: hit=('8R1L','C')
        if hit:
            pdb,ch=hit
            p=(ca1 if pdb=='8R1J' else ca2)[ch][pos]
            if pdb=='8R1J':
                pb1ch='B' if ch=='C' else 'E'
                anp=ca1['G']; pb1=ca1[pb1ch]
            else:
                anp=ca2['D']; pb1=ca2['B']
            anchors=[pb1[a] for a in ANCHOR_POS if a in pb1]
            d_anp=min_dist(p,anp)
            d_cat=min((np.linalg.norm(p-a) for a in anchors),default=np.nan)
            d627=np.linalg.norm(p-ca1['C'][627]) if 627 in ca1['C'] else np.nan
            status=f'resolved-{pdb}-chain{ch}'
        else:
            d_anp=d_cat=d627=np.nan; status='UNRESOLVED-THIN'
        rows.append({**r,'struct_status':status,
            'dist_ANP32B_A':f'{d_anp:.1f}' if d_anp==d_anp else 'NA',
            'dist_catalytic_A':f'{d_cat:.1f}' if d_cat==d_cat else 'NA',
            'dist_PB2_627_A':f'{d627:.1f}' if d627==d627 else 'NA'})
    keys=list(rows[0].keys())
    with open(out_tsv,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys,delimiter='\t'); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r['pos'],r['alt'],r['struct_status'],'ANP32B',r['dist_ANP32B_A'],'cat',r['dist_catalytic_A'],'627',r['dist_PB2_627_A'])
if __name__=='__main__':
    main(sys.argv[1],sys.argv[2])
