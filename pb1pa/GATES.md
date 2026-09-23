# LOCKED GATES - Project H5N1-MAMMAL-JUMP / Builder 2 slice (PB1 + PA)
Locked: 2026-09-23 13:30 IST, BEFORE touching any outcome data. Sign-off: builder 2.

## Question
Which PB1 and PA amino-acid substitutions are statistically enriched in mammalian
(vs avian) H5N1 isolates of the 2.3.4.4b era (2021-2026), which are novel vs
literature-known, and where do they sit on the H5N1 polymerase / human ANP32B complex?

## Prior-art verdict: CROWDED (proceed with differentiation)
Closest works:
1. "Decoding non-human mammalian adaptive signatures of 2.3.4.4b H5N1..." Microbiol Spectr 2025. https://doi.org/10.1128/spectrum.00948-25
2. "Genomic signatures and host adaptation of H5N1 clade 2.3.4.4b" (2025). https://www.sciencedirect.com/science/article/pii/S2666517425000392
3. Gabriel et al. 2015, Nat Commun, experimental H5N1 polymerase screen. https://www.nature.com/articles/ncomms8491
Differentiation: (a) data current to Sept 2026 incl. US dairy-cattle B3.13 and marine-mammal
outbreaks; (b) effect-size + leave-one-host-group-out robustness gates, not marker lookup;
(c) structural mapping of ALL candidates onto PDB 8R1J/8R1L H5N1 polymerase-ANP32B with
confidence tiers; (d) joint PB1+PA map merged with builder 1's PB2 list -> whole-complex view;
(e) reusable auditor tool + byte-locked reproducibility.

## Sampling frame
- NCBI Virus/nuccore: Influenza A H5N1, collection date 2021-01-01..2026-09-23, segments
  2 (PB1), 3 (PA), and 1 (PB2, technical control only), near-complete CDS (>=95% of
  canonical length; PB1 757aa, PA 716aa, PB2 759aa).
- Host classification from /host qualifier: Mammalia vs Aves; unknown/experimental/lab
  strains excluded. Avian background restricted to same period.
- Dedup: identical protein sequences collapsed with isolate counts retained.

## G1 - Positive control (pipeline validation), locked
- Technical control: identical pipeline on PB2 MUST recover E627K as significantly
  mammal-enriched (Fisher exact, BH FDR<0.05, mammalian freq > avian freq).
- PB1/PA control: MUST recover >=2 of the 3 canonical experimentally validated PB1/PA
  mammalian markers (PA T97I [Gabriel 2015 Nat Commun; Virology 2014, PMID 25194918],
  PB1 N105S [Gabriel 2015], PB1-F2 N66S [Conenello 2007 PLoS Pathog, 10.1371/journal.ppat.0030141])
  among significantly enriched sites, PROVIDED each marker's variant allele is carried by
  >=3 mammalian isolates in the sample. Markers absent/below that count are excluded from
  the gate and reported as untestable - never re-fished.
- Literature Tier-2 list (recovery rate reported, NOT gated): PB1 D3V, L13P, D375N, D622G;
  PA K356R, N383D, S224P; PB1 L473V - final list with URLs frozen at data byte-lock.

## G2 - Novel-signal gate, locked
A position is reported as a NOVEL candidate only if ALL hold:
- Fisher exact BH-FDR < 0.05 (per segment, per alignment position, mammal vs avian)
- Haldane-corrected odds ratio >= 5
- mammalian frequency >= 5% AND >= 3x avian frequency
- robust: remains FDR<0.1 in >=80% of leave-one-host-group-out jackknife runs
  (host groups: human, cattle, felids, mustelids, pinnipeds, other mammals)
If no position passes, the result is reported as a negative, with numbers.

## G3 - Data adequacy, locked
>=100 mammalian AND >=500 avian isolates per segment after QC. Otherwise halt, report thin.

## G4 - Structural mapping, locked
Primary: PDB 8R1J / 8R1L (H5N1 polymerase + human ANP32B). Secondary: AlphaFold DB
full-length PB1/PA for regions absent from PDB. Every candidate position mapped;
residues in pLDDT<70 (AF) or unresolved (PDB) regions marked THIN. Distances to ANP32B
interface and polymerase active site reported for every candidate.

## G5 - Reproducibility, locked
Accession lists + SHA256 of every raw file byte-locked at milestone 2; pipeline re-run
from locked bytes must reproduce identical mutation tables.

## Negative results
Preserved and reported verbatim; no re-fishing past a failed gate.

<!-- Fleet-rule compliance: this file is re-committed standalone at head; gates (locked 2026-09-23 13:30 IST, sha256 e8d22d27...) preceded ALL results commits in this repo. -->
