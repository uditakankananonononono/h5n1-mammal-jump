#!/bin/bash
E=https://eutils.ncbi.nlm.nih.gov/entrez/eutils
TERM='"Influenza A virus"[Organism] AND "H5N1"[All Fields] AND "PB2"[Gene] AND ("2020/01/01"[PDAT] : "2026/12/31"[PDAT])'
curl -sS -G $E/esearch.fcgi --data-urlencode db=nuccore --data-urlencode "term=$TERM" \
  --data-urlencode usehistory=y --data-urlencode retmax=0 -o data/esearch.xml
WEB=$(grep -o '<WebEnv>[^<]*' data/esearch.xml | sed 's/<WebEnv>//')
QK=$(grep -o '<QueryKey>[^<]*' data/esearch.xml | sed 's/<QueryKey>//')
N=$(grep -o '<Count>[0-9]*' data/esearch.xml | head -1 | sed 's/<Count>//')
echo "count=$N"
fetch_range () { # $1 rettype $2 outfile $3 batchsize
  local N2=$N; : > "$2"
  for ((s=0; s<N2; s+=$3)); do
    for att in 1 2 3 4 5; do
      tmp=$(mktemp)
      if curl -sS --max-time 90 -G $E/efetch.fcgi --data-urlencode db=nuccore \
        --data-urlencode "query_key=$QK" --data-urlencode "WebEnv=$WEB" \
        --data-urlencode "retstart=$s" --data-urlencode "retmax=$3" \
        --data-urlencode "rettype=$1" --data-urlencode retmode=text -o "$tmp" && [ -s "$tmp" ]; then
        cat "$tmp" >> "$2"; rm "$tmp"; break
      fi
      rm -f "$tmp"; echo "retry $att $1 at $s"; sleep $((att*3))
    done
    echo "$1 $s done $(date +%T)"; sleep 0.34
  done
}
fetch_range fasta_cds_aa data/pb2_cds_aa.fna 400
fetch_range gb data/pb2.gb 500
echo FETCH_COMPLETE
