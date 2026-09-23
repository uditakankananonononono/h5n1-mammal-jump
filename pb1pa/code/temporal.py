#!/usr/bin/env python3
"""Within-cattle monthly emergence of candidate PA alleles (+ PB2 627K reference)."""
import sys, csv, re
from collections import Counter, defaultdict
from Bio import SeqIO
sys.path.insert(0, 'code')
from classify_hosts import classify

ALLELES = {'PA': [(68,'S'),(94,'V'),(137,'R'),(219,'I'),(362,'R'),(432,'I'),(655,'F'),(277,'P'),(113,'R')],
           'PB2': [(627,'K'),(701,'N')]}
JOBS = {'PA': ('parsed/prot_seg3.fa', 'parsed/manifest_seg3.tsv', 'refs/PA_ref.fa', 716),
        'PB2': ('parsed/prot_seg1.fa', 'parsed/manifest_seg1.tsv', 'refs/PB2_ref.fa', 759)}

def month_key(cdate):
    MON = {'Jan':'01','Feb':'02','Mar':'03','Apr':'04','May':'05','Jun':'06','Jul':'07','Aug':'08','Sep':'09','Oct':'10','Nov':'11','Dec':'12'}
    m = re.match(r'^(\d{2})-([A-Za-z]{3})-(\d{4})$', cdate or '')
    if m: return f'{m.group(3)}-{MON[m.group(2)]}'
    m = re.match(r'^([A-Za-z]{3})-(\d{4})$', cdate or '')
    if m: return f'{m.group(2)}-{MON[m.group(1)]}'
    m = re.match(r'^(\d{4})-(\d{2})-(\d{2})$', cdate or '')
    if m: return f'{m.group(1)}-{m.group(2)}'
    return None
MON = {m: i+1 for i, m in enumerate(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])}
def month_sort(k):
    if not k: return '9999-99'
    p = k.split('-')
    if len(p) == 1: return p[0] + '-06'
    mon = p[1]
    mm = f'{MON.get(mon, 6):02d}' if mon.isalpha() else mon
    return f'{p[0]}-{mm}'

w = csv.writer(open('results/temporal_cattle.tsv', 'w'), delimiter='\t')
w.writerow(['protein', 'pos', 'alt', 'month', 'n_cattle', 'n_alt', 'freq'])
for prot, (fa, man, reffa, L) in JOBS.items():
    meta = {}
    for row in csv.DictReader(open(man), delimiter='\t'):
        cls = classify(row['host'], row['strain'])
        if cls != 'mammal': continue
        h = (row['host'] or '').lower(); s = (row['strain'] or '').lower()
        if not re.search(r'\bbos\b|bovine|cattle|\bcow\b|dairy', h or s): continue
        mk = month_key(row['collection_date'])
        if not mk: continue
        meta[row['accession']] = mk
    tot = Counter(); hit = Counter()
    for rec in SeqIO.parse(fa, 'fasta'):
        acc = rec.id.split('|')[0]; gene = rec.id.split('|')[1]
        if gene != prot or acc not in meta: continue
        seq = str(rec.seq).replace('*','').replace('X','')
        if len(seq) != L: continue
        mk = meta[acc]
        tot[mk] += 1
        for pos, alt in ALLELES[prot]:
            if seq[pos-1] == alt: hit[(pos, alt, mk)] += 1
    for pos, alt in ALLELES[prot]:
        for mk in sorted(tot, key=month_sort):
            if tot[mk] < 10: continue
            w.writerow([prot, pos, alt, mk, tot[mk], hit[(pos, alt, mk)], f'{hit[(pos,alt,mk)]/tot[mk]:.4f}'])
print('wrote results/temporal_cattle.tsv')
