#!/usr/bin/env python3
"""Map candidate residues onto 8R1J (H5N1 polymerase dimer + human ANP32B).
Both polymerase copies are checked: copy1 PA=A/PB1=B/PB2=C, copy2 PA=D/PB1=E/PB2=F;
ANP32B=G. A residue unresolved in copy1 falls back to copy2 and the chain used is
recorded. Reports min CA-CA distance to ANP32B, PB1 catalytic Asp anchors, PB2 chain.
Outputs: structure_map.tsv."""
import sys, csv, warnings
warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser
import numpy as np

COPY1 = {'PA': 'A', 'PB1': 'B', 'PB2': 'C'}
COPY2 = {'PA': 'D', 'PB1': 'E', 'PB2': 'F'}
ANP = 'G'
ANCHOR_POS = [305, 306, 445, 446, 447]

def chain_ca(model, ch):
    out = {}
    for res in model[ch]:
        if res.id[0] != ' ':
            continue
        if 'CA' in res:
            out[res.id[1]] = res['CA'].coord
    return out

def min_dist(p, cmap):
    if not cmap:
        return np.nan
    arr = np.array(list(cmap.values()))
    return float(np.sqrt(((arr - p) ** 2).sum(axis=1)).min())

def main(candidates_tsv, out_tsv, cif='structures/8R1J.cif'):
    model = MMCIFParser(QUIET=True).get_structure('x', cif)[0]
    ca = {c: chain_ca(model, c) for c in 'ABCDEFG'}
    anp = ca[ANP]
    rows = []
    for r in csv.DictReader(open(candidates_tsv), delimiter='\t'):
        prot, pos = r['protein'], int(r['pos'])
        if prot not in ('PA', 'PB1'):
            continue
        hit = None
        for copy in (COPY1, COPY2):
            ch = copy[prot]
            if pos in ca[ch]:
                hit = (copy, ch)
                break
        if hit:
            copy, ch = hit
            p = ca[ch][pos]
            pb1 = ca[copy['PB1']]
            pb2 = ca[copy['PB2']]
            anchors = [pb1[a] for a in ANCHOR_POS if a in pb1]
            d_anp = min_dist(p, anp)
            d_cat = min((np.linalg.norm(p - a) for a in anchors), default=np.nan)
            d_pb2 = min_dist(p, pb2)
            status = f'resolved-8R1J-chain{ch}'
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
        w.writeheader()
        w.writerows(rows)
    print(f'wrote {out_tsv} ({len(rows)} rows)')

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
