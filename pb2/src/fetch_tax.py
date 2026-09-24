#!/usr/bin/env python3
import json, re, time, urllib.request, urllib.parse, os
meta=json.load(open('results/gb_meta.json'))
taxids=sorted({m['taxon'] for m in meta.values() if m.get('taxon')}, key=int)
lin=json.load(open('results/taxlineage.json')) if os.path.exists('results/taxlineage.json') else {}
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
B=200
for i in range(0,len(taxids),B):
    batch=[t for t in taxids[i:i+B] if t not in lin]
    if not batch: continue
    url=E+'?'+urllib.parse.urlencode({'db':'taxonomy','id':','.join(batch),'retmode':'xml'})
    for attempt in range(4):
        try:
            x=urllib.request.urlopen(url,timeout=60).read().decode('utf-8','replace'); break
        except Exception as e:
            print('retry',e); time.sleep(3)
    for block in re.findall(r'<Taxon>.*?</Taxon>', x, re.S):
        tid=re.search(r'<TaxId>(\d+)</TaxId>', block).group(1)
        name=re.search(r'<ScientificName>([^<]*)</ScientificName>', block).group(1)
        ids=re.findall(r'<TaxId>(\d+)</TaxId>', block)
        lin[tid]={'name':name,'lineage_ids':ids[1:]}
    print('tax batch',i,'done', flush=True); time.sleep(0.5)
json.dump(lin, open('results/taxlineage.json','w'))
print('total lineages:', len(lin))
