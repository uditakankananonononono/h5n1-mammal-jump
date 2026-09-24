#!/usr/bin/env python3
"""Parse GenBank flatfile dump: per accession -> collection_date, host, taxon id."""
import re, json
recs={}
acc=None; cur=None; src=False
for line in open('data/pb2.gb', errors='replace'):
    if line.startswith('LOCUS'):
        parts=line.split()
        acc=parts[1] if len(parts)>1 else None
        cur={'collection_date':'','host':'','taxon':''}; src=False
    elif acc is not None:
        if line.startswith('VERSION'):
            m=re.search(r'VERSION\s+(\S+)', line)
            if m: cur['accession']=m.group(1)
        elif line.startswith('  source'):
            src=True
        elif src:
            if line.startswith('                     /'):
                m=re.match(r'\s+/(\w+)=?"?([^"]*)"?\s*$', line)
                if m:
                    k,v=m.group(1),m.group(2)
                    if k=='collection_date': cur['collection_date']=v
                    elif k=='host': cur['host']=v
                    elif k=='db_xref' and v.startswith('taxon:'): cur['taxon']=v[6:]
            elif not line.startswith('                     '):
                src=False
    if line.startswith('//') and acc is not None:
        if 'accession' not in cur: cur['accession']=acc
        recs[cur['accession']]=cur; acc=None
json.dump(recs, open('results/gb_meta.json','w'))
print('parsed', len(recs))
