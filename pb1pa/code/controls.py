#!/usr/bin/env python3
"""G1 positive-control validation: curated markers vs enrichment outputs."""
import csv, sys

RES = {'PB2': 'results/PB2_control.enrichment.tsv', 'PB1': 'results/PB1.enrichment.tsv',
       'PA': 'results/PA.enrichment.tsv', 'PB1-F2': 'results/PB1F2.enrichment.tsv'}

def main(out='results/controls.tsv'):
    found = {}
    for prot, path in RES.items():
        try:
            for r in csv.DictReader(open(path), delimiter='\t'):
                found[(prot, int(r['pos']), r['alt'])] = r
        except FileNotFoundError:
            print('MISSING', path)
    rows = []
    for m in csv.DictReader(open('code/markers.tsv'), delimiter='\t'):
        key = (m['protein'], int(m['position']), m['alt_aa'])
        r = found.get(key)
        if r is None:
            status = 'absent_in_data'
            row = {**m, 'm_count': 0, 'av_count': 0, 'freq_m': 0, 'freq_a': 0, 'OR': '', 'q': '', 'status': status}
        else:
            mc = int(r['m']); avc = int(r['av'])
            testable = mc >= 3
            recovered = testable and float(r['q']) < 0.05 and float(r['freq_m']) > float(r['freq_a'])
            status = 'recovered' if recovered else ('testable_not_significant' if testable else 'untestable_rare')
            row = {**m, 'm_count': mc, 'av_count': avc, 'freq_m': r['freq_m'], 'freq_a': r['freq_a'], 'OR': r['OR'], 'q': r['q'], 'status': status}
        rows.append(row)
    keys = list(rows[0].keys())
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
        w.writeheader(); w.writerows(rows)
    t1 = [r for r in rows if r['tier'] == '1_validated' and r['protein'] != 'PB2']
    t1_test = [r for r in t1 if r['status'] != 'untestable_rare' and r['status'] != 'absent_in_data']
    t1_rec = [r for r in t1_test if r['status'] == 'recovered']
    tech = [r for r in rows if r['protein'] == 'PB2'][0]
    print(f"TECHNICAL (PB2 E627K): {tech['status']} q={tech['q']}")
    print(f"TIER1 PB1/PA: {len(t1_rec)}/{len(t1_test)} recovered (gate: >=2)")
    for r in rows:
        print(f"  {r['protein']:7s} {r['position']:>4s}{r['alt_aa']} [{r['tier']}] {r['status']:26s} m={r['m_count']} av={r['av_count']}")

if __name__ == '__main__':
    main()
