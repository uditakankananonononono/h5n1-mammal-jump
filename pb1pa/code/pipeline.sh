#!/bin/bash
# pb1pa end-to-end pipeline (reproducibility entry point)
set -e
cd "$(dirname "$0")/.."
python3 parse_gb.py 1 'data/raw/seg1.*.gb' parsed/manifest_seg1.tsv parsed/prot_seg1.fa
python3 parse_gb.py 2 'data/raw/seg2.*.gb' parsed/manifest_seg2.tsv parsed/prot_seg2.fa
python3 parse_gb.py 3 'data/raw/seg3.*.gb' parsed/manifest_seg3.tsv parsed/prot_seg3.fa
python3 code/classify_hosts.py parsed/manifest_seg1.tsv parsed/manifest_seg2.tsv parsed/manifest_seg3.tsv results/host_audit.tsv
python3 code/enrich.py PB2 parsed/prot_seg1.fa parsed/manifest_seg1.tsv refs/PB2_ref.fa results/PB2_control
python3 code/enrich.py PB1 parsed/prot_seg2.fa parsed/manifest_seg2.tsv refs/PB1_ref.fa results/PB1
python3 code/enrich.py PA  parsed/prot_seg3.fa parsed/manifest_seg3.tsv refs/PA_ref.fa results/PA
python3 code/enrich.py PB1-F2 parsed/prot_seg2.fa parsed/manifest_seg2.tsv refs/PB1F2_ref.fa results/PB1F2
echo DONE
