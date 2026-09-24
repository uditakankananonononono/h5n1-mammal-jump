#!/usr/bin/env python3
"""Direct residue composition at curated marker positions (no test-row dependence)."""
import sys, csv
from collections import Counter
from Bio import SeqIO
sys.path.insert(0, 'code')
from classify_hosts import classify

JOBS = {'PB1': ('parsed/prot_seg2.fa', 'parsed/manifest_seg2.tsv', 'refs/PB1_ref.fa'),
        'PA': ('parsed/prot_seg3.fa', 'parsed/manifest_seg3.tsv', 'refs/PA_ref.fa'),
        'PB1-F2': ('parsed/prot_seg2.fa', 'parsed/manifest_seg2.tsv', 'refs/PB1F2_ref.fa'),
        'PB2': ('parsed/prot_seg1.fa', 'parsed/manifest_seg1.tsv', 'refs/PB2_ref.fa')}

def main():
    markers = [m for m in csv.DictReader(open('code/markers.tsv'), delimiter='\t')]
    out = {}
    for prot, (fa, man, reffa) in JOBS.items():
        poss = {int(m['position']) for m in markers if m['protein'] == prot}
        if not poss: continue
        ref = str(next(SeqIO.parse(reffa, 'fasta')).seq)
        meta = {}
        for row in csv.DictReader(open(man), delimiter='\t'):
            cls = classify(row['host'], row['strain'])
            if not row['year']: continue
            meta.setdefault(row['accession'], cls)
        mc = {p: Counter() for p in poss}; ac = {p: Counter() for p in poss}
        for rec in SeqIO.parse(fa, 'fasta'):
            acc = rec.id.split('|')[0]; gene = rec.id.split('|')[1]
            if gene != prot or acc not in meta: continue
            seq = str(rec.seq).replace('*','').replace('X','')
            if len(seq) != len(ref): continue
            cls = meta[acc]
            if cls not in ('mammal','avian'): continue
            tgt = mc if cls=='mammal' else ac
            for p in poss:
                if p <= len(seq):
                    tgt[p][seq[p-1]] += 1
        for p in poss:
            out[(prot,p)] = (mc[p], ac[p])
    w = csv.writer(open('results/marker_composition.tsv','w'), delimiter='\t')
    w.writerow(['protein','position','marker_alt','tier','alt_m','alt_av','tot_m','tot_av','freq_m','freq_a','avian_consensus','mammal_consensus'])
    for m in markers:
        key = (m['protein'], int(m['position']))
        mc, ac = out.get(key, (Counter(), Counter()))
        tm, ta = sum(mc.values()), sum(ac.values())
        am, aa = mc[m['alt_aa']], ac[m['alt_aa']]
        w.writerow([m['protein'], m['position'], m['alt_aa'], m['tier'], am, aa, tm, ta,
                    f'{am/tm:.5f}' if tm else 'NA', f'{aa/ta:.5f}' if ta else 'NA',
                    ac.most_common(1)[0][0] if ac else '-', mc.most_common(1)[0][0] if mc else '-'])
        print(f"{m['protein']:7s} {m['position']:>4s}{m['alt_aa']} [{m['tier'][:1]}] m={am}/{tm} ({am/tm:.4f}) av={aa}/{ta} ({aa/ta:.5f}) avian_cons={ac.most_common(1)[0][0] if ac else '-'}")
if __name__ == '__main__':
    main()
