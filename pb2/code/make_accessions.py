#!/usr/bin/env python3
"""Accession inventory from raw genbank (all records) + parsed manifests (used records)."""
import glob, re, csv, sys
for seg in ['1','2','3']:
    accs = []
    for f in sorted(glob.glob(f'data/raw/seg{seg}.*.gb')):
        accs += re.findall(r'^VERSION\s+(\S+)', open(f, errors='replace').read(), re.M)
    open(f'manifests/accessions_seg{seg}_raw.txt','w').write('\n'.join(sorted(set(accs)))+'\n')
    print(f'seg{seg}: {len(set(accs))} raw accessions')
