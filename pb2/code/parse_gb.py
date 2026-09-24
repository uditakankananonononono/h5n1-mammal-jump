#!/usr/bin/env python3
"""Parse H5N1 genbank segments -> manifest TSV + protein FASTA.
Fields: accession, segment, strain, host, collection_date, country, gene, length_aa
QC: collection year 2021+, CDS >= 95% canonical length.
"""
import glob, re, csv, sys
from Bio import SeqIO

CANON = {'PB2':759, 'PB1':757, 'PA':716}
MINFRAC = 0.95

def parse_date(s):
    """Return (year, month or None) from genbank collection_date formats."""
    if not s: return (None, None)
    m = re.match(r'^(\d{4})$', s.strip())
    if m: return (int(m.group(1)), None)
    m = re.match(r'^(\d{2})-([A-Za-z]{3})-(\d{4})$', s.strip())
    if m: return (int(m.group(3)), m.group(2))
    m = re.match(r'^([A-Za-z]{3})-(\d{4})$', s.strip())
    if m: return (int(m.group(2)), m.group(1))
    m = re.match(r'^(\d{4})-(\d{2})-(\d{2})$', s.strip())
    if m: return (int(m.group(1)), m.group(2))
    return (None, None)

def run(seg, gz_path, man_out, fa_out):
    recs = 0; kept = 0
    man = open(man_out, 'w', newline='')
    w = csv.writer(man, delimiter='\t')
    w.writerow(['accession','segment','strain','host','collection_date','year','country','gene','protein_len','seq_md5','ncbi_title'])
    fa = open(fa_out, 'w')
    import hashlib
    import itertools
    def records():
        for f in sorted(glob.glob(gz_path)):
            with open(f, errors='replace') as fh:
                yield from SeqIO.parse(fh, 'genbank')
    if True:
        for r in records():
            if 'H5N1' not in r.description: continue
            recs += 1
            src = next((f for f in r.features if f.type=='source'), None)
            q = src.qualifiers if src else {}
            host = (q.get('host') or [''])[0].strip()
            strain = (q.get('strain') or q.get('isolate') or [''])[0].strip()
            cdate = (q.get('collection_date') or [''])[0].strip()
            country = (q.get('country') or [''])[0].strip()
            year, mon = parse_date(cdate)
            for f in r.features:
                if f.type != 'CDS': continue
                gene = (f.qualifiers.get('gene') or [''])[0].strip()
                prod = (f.qualifiers.get('product') or [''])[0].strip()
                tr = (f.qualifiers.get('translation') or [''])[0].strip()
                if not tr: continue
                key = gene if gene in CANON else None
                if not key:
                    if prod.startswith('polymerase PB2') or prod=='PB2': key='PB2'
                    elif prod.startswith('polymerase PB1') or prod=='PB1': key='PB1'
                    elif prod.startswith('polymerase PA') or prod=='PA' or prod=='polymerase acidic protein': key='PA'
                    elif 'PB1-F2' in prod or 'PB1-F2' in gene: key='PB1-F2'
                if not key: continue
                if key in CANON and len(tr) < CANON[key]*MINFRAC: continue
                if year is not None and year < 2021: continue
                md5 = hashlib.md5(tr.encode()).hexdigest()
                w.writerow([r.id.split('.')[0], seg, strain, host, cdate, year or '', country, key, len(tr), md5, r.description[:120]])
                fa.write(f'>{r.id.split(".")[0]}|{key}|{host[:40]}|{cdate}|{strain[:60]}\n{tr}\n')
                kept += 1
    man.close(); fa.close()
    print(f'seg{seg}: records={recs} kept_CDS={kept}')

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
