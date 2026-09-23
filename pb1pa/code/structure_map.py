#!/usr/bin/env python3
"""Map candidate residues onto 8R1J (H5N1 polymerase dimer + human ANP32B).
Reports: resolved?, CA coords, min distance to ANP32B, to PB1 catalytic Asp, to PB2 chain.
Outputs: structure_map.tsv + 3D figures."""
import sys, csv, warnings
warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
import numpy as np

CHAINS = {'PA': 'A', 'PB1': 'B', 'PB2': 'C', 'ANP32B': 'G'}

def load_coords(cif='structures/8R1J.cif'):
    s = MMCIFParser(QUIET=True).get_structure('x', cif)
    model = s[0]
    coords = {}
    for name, ch in CHAINS.items():
        chain = model[ch]
        cmap = {}
        for res in chain:
            if res.id[0] != ' ': continue
            if 'CA' in res:
                cmap[res.id[1]] = res['CA'].coord
        coords[name] = cmap
    return coords

def min_dist(p, cmap):
    if not cmap: return np.nan
    arr = np.array(list(cmap.values()))
    return float(np.sqrt(((arr - p) ** 2).sum(axis=1)).min())

def main(candidates_tsv, out_tsv):
    coords = load_coords()
    # catalytic anchors: find Asp in PB1 motifs (D305 region 300-315, SDD 440-455)
    pb1 = coords['PB1']
    anchors = [p for p in [305, 306, 445, 446, 447] if p in pb1]
    anp = coords['ANP32B']
    pb2 = coords['PB2']
    rows = []
    for r in csv.DictReader(open(candidates_tsv), delimiter='\t'):
        prot, pos = r['protein'], int(r['pos'])
        ch = 'PA' if prot == 'PA' else 'PB1' if prot == 'PB1' else None
        if ch is None: continue
        cmap = coords[ch]
        if pos in cmap:
            p = cmap[pos]
            d_anp = min_dist(p, anp)
            d_cat = min((np.linalg.norm(p - pb1[a]) for a in anchors), default=np.nan)
            d_pb2 = min_dist(p, pb2)
            status = 'resolved'
        else:
            d_anp = d_cat = d_pb2 = np.nan
            status = 'UNRESOLVED-thin'
        rows.append({**r, 'struct_status': status,
                     'dist_ANP32B_A': f'{d_anp:.1f}' if d_anp == d_anp else 'NA',
                     'dist_catalytic_A': f'{d_cat:.1f}' if d_cat == d_cat else 'NA',
                     'dist_PB2_A': f'{d_pb2:.1f}' if d_pb2 == d_pb2 else 'NA'})
    keys = list(rows[0].keys()) if rows else []
    with open(out_tsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
        w.writeheader(); w.writerows(rows)
    print(f'mapped {len(rows)} candidates; unresolved={sum(1 for r in rows if r["struct_status"]!="resolved")}')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
