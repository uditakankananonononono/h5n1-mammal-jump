# LOCKED GATES - Project H5N1-MAMMAL-JUMP / Builder 1 slice (PB2)
Locked: 2026-09-24 08:10 IST, BEFORE touching any outcome data. Sign-off: B01 v3 builder
(agent-01M38M9G8VARPSWXTCCJJV9DZB, key SHA256:r6of2QA9CbmipjNWsmV2uaSBk2+JWc6ZKl5xchyAKS0).
Mirrors the pb1pa slice gates (locked 2026-09-23 13:30 IST) for the PB2 segment.

## Question
Which PB2 amino-acid substitutions are statistically enriched in mammalian (vs avian)
H5N1 isolates of the 2.3.4.4b era (2021-2026), which are novel vs literature-known, and
where do they sit on the H5N1 polymerase / human ANP32B complex?

## Prior-art verdict: CROWDED (proceed with differentiation)
Same closest works as pb1pa gates (Microbiol Spectr 2025 doi:10.1128/spectrum.00948-25;
SciDirect S2666517425000392; Gabriel 2015 Nat Commun ncomms8491) plus the classic PB2
host-determinant literature (Subbarao 1993 E627K; Steel 2009 D701N; Mehle 2006/2010
ANP32-dependence). Differentiation: (a) data current to Sept 2026 incl. US dairy-cattle
B3.13 and marine-mammal outbreaks; (b) effect-size + leave-one-host-group-out robustness,
not marker lookup; (c) structural mapping of ALL candidates onto PDB 8R1J/8R1L with
confidence tiers; (d) feeds builder 2's joint whole-complex map; (e) reusable auditor +
byte-locked reproducibility.

## Sampling frame
- NCBI nuccore: Influenza A H5N1 ("Influenza A virus"[Organism] AND "H5N1"[All Fields]),
  PB2 gene (segment 1), retrieved 2026-09-24 via esearch+efetch, collection date
  2021-01-01..2026-09-23 from the /collection_date qualifier (retrieval PDAT window
  2020-01-01.. to avoid edge loss), near-complete CDS (>=95% of canonical 759 aa /
  2277 nt).
- Host classification from the /host qualifier: Mammalia vs Aves by NCBI Taxonomy
  lineage (Mammalia taxid 40674, Aves taxid 8782). Unknown/experimental/lab strains
  excluded: records with no resolvable host lineage, and records whose /host is a
  laboratory model with no wild isolation evidence (Mus musculus) are EXCLUDED, listed
  by accession in the audit report.
- Avian background restricted to the same period.
- Dedup: identical protein sequences collapsed with isolate counts retained.

## G1 - Positive control (pipeline validation), locked
- Technical control: the pipeline MUST recover E627K as significantly mammal-enriched
  (Fisher exact, BH FDR<0.05, mammalian freq > avian freq).
- PB2 literature control: MUST recover >=2 of the 3 canonical experimentally validated
  PB2 mammalian markers (D701N [Steel 2009 PLoS Pathog 10.1371/journal.ppat.1000652],
  Q591K [Mehle 2006 PLoS Pathog; Yamada 2010 Nature], T271A [Bussey 2010 J Virol
  10.1128/JVI.01199-10]) among significantly enriched sites, PROVIDED each marker's
  variant allele is carried by >=3 mammalian isolates in the sample. Markers absent or
  below that count are excluded from the gate and reported as untestable - never
  re-fished.
- Literature Tier-2 list (recovery rate reported, NOT gated): PB2 E627V, K526R, D740N,
  A199S, D253N, M64T, E158G, D1956E(n/a excluded), G590S/Q591R (the "Q591K-like" pair
  counted separately). Final list frozen at data byte-lock.

## G2 - Novel-signal gate, locked
A position is reported as a NOVEL candidate only if ALL hold:
- Fisher exact BH-FDR < 0.05 (per alignment position, mammal vs avian)
- Haldane-corrected odds ratio >= 5
- mammalian frequency >= 5% AND >= 3x avian frequency
- robust: remains FDR<0.1 in >=80% of leave-one-host-group-out jackknife runs
  (host groups: human, cattle, felids, mustelids, pinnipeds, other mammals)
If no position passes, the result is reported as a negative, with numbers.

## G3 - Data adequacy, locked
>=100 mammalian AND >=500 avian isolates after QC. Otherwise halt, report thin.

## G4 - Structural mapping, locked
Primary: PDB 8R1J / 8R1L (H5N1 polymerase + human ANP32B). Secondary: AlphaFold DB
full-length PB2 for regions absent from PDB. Every candidate position mapped; residues
in pLDDT<70 (AF) or unresolved (PDB) regions marked THIN. Distances to the ANP32B
interface and to the polymerase active site reported for every candidate.

## G5 - Reproducibility, locked
Accession lists + SHA256 of every raw file byte-locked at milestone 2; pipeline re-run
from locked bytes must reproduce identical mutation tables.

## Negative results
Preserved and reported verbatim; no re-fishing past a failed gate.
