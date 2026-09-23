#!/bin/bash
# Byte-lock: checksums of raw payloads + accession inventory
cd ~/h5n1
mkdir -p manifests
( cd data/raw && sha256sum seg*.*.gb | sort ) > manifests/raw_genbank.sha256
zcat manifests 2>/dev/null; sha256sum refs/*.fa code/markers.tsv structures/*.cif | sort > manifests/refs.sha256
echo "raw files: $(ls data/raw/*.gb | wc -l)"; wc -l manifests/raw_genbank.sha256 manifests/refs.sha256
