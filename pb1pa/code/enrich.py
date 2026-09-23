#!/usr/bin/env python3
"""Position-resolved mammalian-enrichment analysis for one protein.
Usage: enrich.py <protein> <fasta> <manifest> <ref_fasta> <out_prefix> [--include-lab]
Reference numbering = 8R1J construct sequence (PB1 757aa / PA 716aa).
"""
import sys, csv, re, hashlib
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from scipy.stats import fisher_exact
import numpy as np

sys.path.insert(0, 'code')
from classify_hosts import classify

HOST_GROUPS = [
    ('human',      [r'homo sapiens', r'\bhuman\b']),
    ('cattle',     [r'\bbos\b', r'bovine', r'cattle', r'\bcow\b', r'dairy']),
    ('felids',     [r'felis', r'\bcat\b', r'feline', r'lynx', r'puma', r'panthera', r'leopard', r'tiger', r'\blion\b', r'acinonyx', r'cheetah']),
    ('mustelids',  [r'neovison', r'\bmink\b', r'mustela', r'marten', r'martes', r'fisher', r'gulo', r'wolverine', r'badger', r'otter', r'lutra', r'lontra', r'enhydra']),
    ('pinnipeds',  [r'seal', r'phoca', r'halichoerus', r'mirounga', r'pusa', r'otaria', r'sea lion', r'zalophus', r'eumetopias', r'arctocephalus', r'callorhinus', r'odobenus', r'walrus', r'hydrurga', r'leptonychotes', r'lobodon']),
]
def host_group(host, strain):
    h = (host or '').lower(); st = (strain or '').lower()
    m = re.match(r'^a/([^/]+)/', st); text = h if h else (m.group(1) if m else '')
    for g, pats in HOST_GROUPS:
        for p in pats:
            if re.search(p, text): return g
    return 'other_mammal'

def load_ref(path):
    r = next(SeqIO.parse(path, 'fasta'))
    return str(r.seq)

def map_to_ref(seq, ref, aligner):
    """Return dict ref_pos(1-based)->aa for this sequence."""
    if len(seq) == len(ref):
        return {i+1: c for i, c in enumerate(seq)}
    aln = aligner.align(ref, seq)[0]
    blocks = aln.aligned
    out = {}
    for (rs, re_), (qs, qe) in zip(blocks[0], blocks[1]):
        for k in range(re_ - rs):
            out[rs + k + 1] = seq[qs + k]
    return out

def main():
    protein, fasta, manifest, ref_fa, out_prefix = sys.argv[1:6]
    include_lab = '--include-lab' in sys.argv
    ref = load_ref(ref_fa)
    L = len(ref)
    # manifest: accession -> (class, group, year)
    meta = {}
    for row in csv.DictReader(open(manifest), delimiter='\t'):
        cls = classify(row['host'], row['strain'])
        if cls == 'mammal_lab' and include_lab: cls = 'mammal'
        if not row['year']: continue  # frame: known 2021+ only
        meta.setdefault(row['accession'], (cls, host_group(row['host'], row['strain']), row['year']))
    aligner = PairwiseAligner()
    aligner.mode = 'global'; aligner.match_score = 2; aligner.mismatch_score = -1
    aligner.open_gap_score = -5; aligner.extend_gap_score = -1
    m_counts = defaultdict(Counter)  # pos -> residue -> n (mammal)
    a_counts = defaultdict(Counter)
    m_groups = defaultdict(lambda: defaultdict(Counter))  # group -> pos -> residue -> n
    n_m = n_a = 0
    m_iso_group = Counter()
    seen_iso = set()
    def norm_strain(s):
        s = (s or '').lower()
        s = re.sub(r'-(original|passaged?|l\d+|egg|mdck|reassortant.*|clone.*)$', '', s)
        return s
    hap_m = defaultdict(Counter); hap_a = defaultdict(Counter)
    seen_hap = set()
    ident_samples = []
    for rec in SeqIO.parse(fasta, 'fasta'):
        acc = rec.id.split('|')[0]
        gene = rec.id.split('|')[1] if '|' in rec.id else ''
        if gene != protein: continue
        if acc not in meta: continue
        cls, grp, year = meta[acc]
        if cls not in ('mammal', 'avian'): continue
        seq = str(rec.seq).replace('*', '').replace('X', '')
        if len(seq) < L * 0.9: continue
        mp = map_to_ref(seq, ref, aligner)
        if len(ident_samples) < 200 and len(seq) == len(ref):
            ident_samples.append(sum(1 for i in range(L) if seq[i] == ref[i]) / L)
        md5 = hashlib.md5(seq.encode()).hexdigest()
        parts = rec.description.split('|')
        strain_field = parts[4] if len(parts) > 4 else ''
        iso_key = (norm_strain(strain_field), md5)
        if iso_key in seen_iso: continue
        seen_iso.add(iso_key)
        is_new_hap = md5 not in seen_hap
        seen_hap.add(md5)
        if cls == 'mammal':
            n_m += 1; m_iso_group[grp] += 1
            for p, c in mp.items():
                m_counts[p][c] += 1; m_groups[grp][p][c] += 1
                if is_new_hap: hap_m[p][c] += 1
        else:
            n_a += 1
            for p, c in mp.items():
                a_counts[p][c] += 1
                if is_new_hap: hap_a[p][c] += 1
    print(f'{protein}: n_mammal={n_m} n_avian={n_a} groups={dict(m_iso_group)}')
    if ident_samples:
        print(f'{protein}: median identity to ref (same-length sample) = {np.median(ident_samples):.3f}')
    # enrichment tests
    rows = []
    for p in range(1, L + 1):
        if not a_counts[p]: continue
        av_cons = a_counts[p].most_common(1)[0][0]
        res_set = set(m_counts[p]) | set(a_counts[p])
        for r in res_set:
            if r == av_cons or r == '-': continue
            a = m_counts[p][r]; c = a_counts[p][r]
            if a == 0 and c == 0: continue
            b = n_m - a; d = n_a - c
            if a + c < 3: continue  # minimum information
            or_ = ((a + .5) * (d + .5)) / ((b + .5) * (c + .5))
            pval = fisher_exact([[a, b], [c, d]], alternative='greater')[1]
            rows.append(dict(protein=protein, pos=p, avian_cons=av_cons, alt=r,
                             m=a, m_other=b, av=c, av_other=d,
                             freq_m=a / n_m, freq_a=c / n_a, OR=or_, p=pval))
    # BH FDR
    rows.sort(key=lambda x: x['p'])
    N = len(rows)
    for i, r in enumerate(rows):
        r['q'] = min(1.0, r['p'] * N / (i + 1))
    # monotone fix
    qs = [r['q'] for r in rows]
    for i in range(N - 2, -1, -1):
        rows[i]['q'] = min(rows[i]['q'], rows[i + 1]['q'])
    # G2 flags
    for r in rows:
        r['g2_pass'] = (r['q'] < 0.05 and r['OR'] >= 5 and r['freq_m'] >= 0.05 and r['freq_m'] >= 3 * r['freq_a'])
        # haplotype view
        hm = hap_m[r['pos']][r['alt']]; ha = hap_a[r['pos']][r['alt']]
        r['hap_m'] = hm; r['hap_a'] = ha
    with open(f'{out_prefix}.enrichment.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ['protein'], delimiter='\t')
        w.writeheader(); w.writerows(rows)
    # group frequencies for candidates
    cand = [r for r in rows if r['g2_pass']]
    with open(f'{out_prefix}.groupfreq.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['protein', 'pos', 'alt', 'group', 'freq', 'n_group'])
        for r in cand:
            for g in m_iso_group:
                n_g = m_iso_group[g]
                cnt = m_groups[g][r['pos']][r['alt']]
                w.writerow([protein, r['pos'], r['alt'], g, f'{cnt / n_g:.4f}', n_g])
    # jackknife for candidates
    groups = [g for g, n in m_iso_group.items() if n >= 10]
    jk = {}
    for r in cand:
        keep = 0; tot = 0
        for g in groups:
            a2 = r['m'] - m_groups[g][r['pos']][r['alt']]
            n_m2 = n_m - m_iso_group[g]
            if n_m2 < 20: continue
            b2 = n_m2 - a2
            pval = fisher_exact([[a2, b2], [r['av'], r['av_other']]], alternative='greater')[1]
            tot += 1
            if pval < 0.01: keep += 1
        jk[(r['pos'], r['alt'])] = (keep, tot)
    with open(f'{out_prefix}.jackknife.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['protein', 'pos', 'alt', 'jk_keep', 'jk_total', 'jk_frac', 'robust'])
        for r in cand:
            k, t = jk[(r['pos'], r['alt'])]
            w.writerow([protein, r['pos'], r['alt'], k, t, f'{k / t:.3f}' if t else 'NA', (k / t >= 0.8) if t else False])
    print(f'{protein}: tested={N} g2_candidates={len(cand)}')

if __name__ == '__main__':
    main()
