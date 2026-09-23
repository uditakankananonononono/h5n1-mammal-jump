# pb1pa data manifest - sources and regeneration

All checksums are SHA-256. Raw payloads are byte-locked in `data/raw_genbank.sha256` (36 files),
`data/refs.sha256` (reference sequences + structures), `data/parsed.sha256` (6 derived files).
Pipeline: `code/pipeline.sh` (fetch -> parse -> enrich -> controls -> figures -> structure map).

## Raw GenBank (36 batch files, seg{1,2,3}.<offset>.gb)
Source: NCBI E-utilities, db=nuccore.
- esearch (per segment N in 1,2,3): `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=H5N1%5BAll+Fields%5D+AND+%22segment+N%22%5BAll+Fields%5D+AND+2021%3A2026%5BPDAT%5D&retmax=0&usehistory=y`
- efetch batches of 2000 via WebEnv/query_key: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&query_key=1&WebEnv=<env>&retstart=<offset>&retmax=2000&rettype=gb&retmode=text`
- Record counts at fetch time: seg1 23,938 / seg2 23,913 / seg3 23,927 (2026-09-23). Counts grow as GenBank grows; the byte-locked files and `data/accessions_seg{1,2,3}_raw.txt` freeze the analyzed set.
- Accession lists: same esearch with `rettype=acc` (`code/make_accessions.py`).

## Parsed tables (derived; regenerate with `code/parse_gb.py <seg> 'data/raw/segN.*.gb' parsed/manifest_segN.tsv parsed/prot_segN.fa`)
- parsed/manifest_seg{1,2,3}.tsv, parsed/prot_seg{1,2,3}.fa - checksums in `data/parsed.sha256`. Re-parse of the locked raw bytes reproduces them byte-identically (verified twice, independent sandboxes).

## Reference sequences (checksums in `data/refs.sha256`)
- `refs/PB1_ref.fa`, `refs/PA_ref.fa`: 8R1J construct sequences (PB1 757 aa, PA 716 aa), from PDB 8R1J entity sequences - `https://www.rcsb.org/fasta/entry/8R1J`
- `refs/PB2_ref.fa` = `refs/pb2_ref_query.fa`: UniProt Q9Q0V1 (A/Goose/Guangdong/1/1996 PB2) - `https://rest.uniprot.org/uniprotkb/Q9Q0V1.fasta`
- `refs/PB1F2_ref.fa` = basis of `refs/pb1f2_query.fa`: UniProt P0C791 (PB1-F2) - `https://rest.uniprot.org/uniprotkb/P0C791.fasta`
- `code/markers.tsv`: curated literature marker panel (refs in paper).

## Structures (checksums in `data/refs.sha256`)
- `structures/8R1J.cif` - `https://files.rcsb.org/download/8R1J.cif`
- `structures/8R1L.cif` - `https://files.rcsb.org/download/8R1L.cif`
Chain mapping used by `code/structure_map.py`: 8R1J PA=A/D, PB1=B/E, PB2=C/F, ANP32B=G; 8R1L ANP32B=A, PA=B (construct spans PA 202-716), PB1=C, PB2=D.
