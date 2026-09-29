# NMPFamsDB2 Pipeline: Fixes and Added Steps

2026-09-29 · Igor Tolstoy

## Summary

The NMPFamsDB2 pipeline now runs end to end on macOS and Linux with eight corrected steps and seven steps from the paper added as commands. The starting point was [PavlopoulosLab/NMPFamsDB2](https://github.com/PavlopoulosLab/NMPFamsDB2) at commit `1cdf36e`, which is byte-identical to the Zenodo code archive ([record 22253220](https://zenodo.org/records/22253220)), checked against the paper ([doi:10.1093/nar/gkag938](https://doi.org/10.1093/nar/gkag938)).

- **Corrected:** contig join order, low-complexity filter, Pfam hit parsing, join sort order, DIAMOND input file, DIAMOND coverage, alignment trimming guide, macOS multiprocessing in all three scripts, DIAMOND database naming.
- **Added from the paper:** AntiFam filter, LAST before DIAMOND, family profile generation, the 16-sequence cutoff, ColabFold prediction, Foldseek search, structural superfamilies.
- **Files changed:** `commands.md`, `aligner.py`, `trimmer.py`, `redundancy_removal.py`; added `Dockerfile` and this file.

The `Dockerfile` builds a linux/amd64 image with every tool below at the tested version (11.1 GB; all tools respond inside the image). BLAST+ 2.17.0 has no linux-aarch64 build on bioconda, so the image is pinned to amd64. The original repository is kept as the `upstream` remote.

## Environment

All tests ran on macOS arm64 in the conda env `nmpfams`; ColabFold has its own env `colabfold`. Two tools are newer than the paper's versions (HMMER 3.3.2, tantan 26).

| Tool | Version tested | Used for |
| --- | --- | --- |
| HMMER | 3.4 | Pfam and AntiFam search |
| tantan | 51 | low-complexity masking |
| LAST | 1654 | first reference search |
| DIAMOND | 2.2.8 | second reference search |
| MMseqs2 | 18.8cc5c | Linclust clustering |
| MAFFT | 7.526 | family alignment |
| HH-suite (hhfilter, hhmake) | 3.3.0 | redundancy filter |
| Foldseek | 10.941cd33 | structure search and clustering |
| ColabFold | 1.5.5 | structure prediction (options checked, not run) |
| BLAST+ | 2.17.0 | remote ClusteredNR search |
| Python | 3.11, pandas 3.0.6, scikit-learn 1.9.1 | the three scripts |

| Database | Version | Check |
| --- | --- | --- |
| Pfam-A | 38.2 (Jan 2026; paper used 37.0) | md5 `7ab3c4e215d0daaea3004e37c4e24f8a` matches EBI; 30,134 models, all with TC |
| AntiFam | 8.0 (Oct 2023) | 278 models; 278 GA, 1 TC |

## Databases and resources for a full run

The downloadable databases total about 130 GB compressed; the IMG/M protein inputs add an estimated 1–6.5 TB. Sizes are the servers' `Content-Length` unless marked as estimates.

| Database | Step | Download size |
| --- | --- | --- |
| Pfam-A 38.2 | Pfam filter | 0.42 GB (2.25 GB unpacked) |
| AntiFam 8.0 | AntiFam filter | 0.006 GB |
| IMG/M reference proteins, 539,040,456 | LAST + DIAMOND | ~85 GB FASTA (estimate) |
| IMG/M metagenome and metatranscriptome ORFs, 40.3 billion, + scaffold lengths | input | ~6.4 TB raw, ~1.0 TB after filtering (estimate) |
| Foldseek `CATH50` | structure search | 1.01 GB |
| Foldseek `PDB` (`pdb100`) | structure search | 2.33 GB |
| Foldseek `Alphafold/UniProt50` | structure search, no CATH/PDB hit | 122.56 GB (full `afdb`: 491.47 GB) |
| AlphaFold2 weights (`alphafold_params_colab_2022-12-06`) | ColabFold | 4.10 GB |
| NMPFamsDB2 results, Zenodo 17225887 (optional) | reuse, comparison | 26.22 GB |

The IMG/M estimates use 157.7 bytes of FASTA per protein, measured on 2 million NMPFamsDB2 proteins (mean 103.4 aa). The paper notes the filters enrich for smaller proteins, so raw ORFs are likely longer and the raw estimate low. The paper reports no CPU, GPU or runtime figures.

| Step | Resource estimate | Basis |
| --- | --- | --- |
| hmmsearch vs Pfam, 6.23 billion proteins | ~48,000 CPU-hours | measured 455 CPU-s for 7,000 proteins of mean 240 aa (2.7×10⁻⁴ CPU-s per residue), scaled to 103 aa |
| LAST, DIAMOND vs 539 million proteins | not measured | — |
| MMseqs2 Linclust, 1.74 billion proteins | not measured | — |
| MAFFT, trimming, hhfilter, 608,258 families | not measured | — |
| ColabFold, 156,608 models (Zenodo HQ + MQ + LQ counts) | GPU; not measured | 5 models × 3 recycles per family |

## Fixes

Each fix was tested by running the original command and the corrected one on the same input. Fixes 3 and 4 used the real HMMER 3.4 output of 7,000 test proteins against Pfam 38.2 (4,617 with a Pfam hit, 2,383 without).

| # | Step | Behaviour at `1cdf36e` | Change | Test result |
| --- | --- | --- | --- | --- |
| 1 | Remove genes near contig ends | `join` with unsorted `scaffolds_info.tsv` drops proteins with no warning | sort both inputs with `LC_ALL=C` | valid proteins kept: 1 of 2 before, 2 of 2 after |
| 2 | Low-complexity filter | proteins with a run of 10 `X` skip the ≥35 aa check; `sed 's/X//g'` also edits IDs | one rule: keep if ≥35 unmasked residues; sequence keeps its `X`; ID untouched | 18-residue protein (62 of 80 masked) dropped; `GaX01` ID kept |
| 3 | Collect Pfam hits | `awk -F '\t'` on space-delimited domtblout returns whole lines | skip `#` lines, split on whitespace, `sort -u` | proteins left after Pfam removal: 7,000 before, 2,383 after |
| 4 | Remove Pfam hits | unsorted protein table into `join -v` | `LC_ALL=C sort` before `join`; same sort order for the reference join | Pfam hits passed as novel: 4,581 before, 0 after |
| 5 | Collect DIAMOND hits | reads `output_table`; DIAMOND writes `matches.m8` | read `matches.m8` | reference hits found: 0 before, 2 of 2 after |
| 6 | DIAMOND coverage | `(end - start) / length` is one residue short | `(end - start + 1) / length` | alignment covering exactly 70%: 0.69 not a hit before, 0.70 hit after |
| 7 | `trimmer.py` on macOS | pool created at import; spawn workers re-run it | `if __name__ == "__main__":` guard | before: killed at 30 s, 324 RuntimeErrors; after: exit 0, 5 rows |
| 8 | Trimming guide | first sequence used (`first_is_representative = True`) | central sequence (minimum mean Hamming distance), as in the paper | outlier-first cluster: 4 columns kept before, 12 after |
| 9 | `aligner.py`, `redundancy_removal.py` on macOS | same as 7 | same guard | before: killed at 60 s (2,484 and 2,498 RuntimeErrors); after: 6 aligned, 5 kept |
| 10 | DIAMOND databases | all four reference files built as `-d nr`; the last one wins | one database per file, results merged | queries matched: 1 of 4 before, 4 of 4 after |

Fix 8 also changes what `hhfilter -M first` measures against. hhfilter uses the first sequence for match columns and coverage, which is now the central sequence. In the outlier test, the old order kept only the outlier; the new order keeps the 3 related sequences and drops the outlier (coverage 4/12).

## Masked residues (X)

Keeping the masked residues as `X` in the searched sequence finds more known proteins than deleting them. Deleting them joins the two flanks, so a homolog must align across a gap the length of the masked stretch. A protein that loses its hit this way is passed on as novel. The paper describes tantan masking but not deletion; the deletion comes from the repository's `sed 's/X//g'`.

**DIAMOND** — a 150 aa reference; query = the same protein mutated to the given identity, residues 51–100 masked; 50 random mutants per level; hit = ≥30% identity, ≥70% coverage of both proteins.

| Identity kept | Hit both ways | Hit only with X kept | Hit only when deleted | No hit |
| --- | --- | --- | --- | --- |
| 50% | 3 | 12 | 0 | 35 |
| 60% | 8 | 27 | 0 | 15 |
| 70% | 32 | 14 | 0 | 4 |
| 80% | 45 | 5 | 0 | 0 |
| 90% | 50 | 0 | 0 | 0 |

**hmmsearch `--cut_tc` against Pfam 38.2** — 1,000 sequences sampled with `hmmemit` from 200 random families (model length ≥120), a centred window masked. Counts are for the 875 sequences with a Pfam hit before masking.

| Masked residues | Hit both ways | Hit only with X kept | Hit only when deleted | No hit |
| --- | --- | --- | --- | --- |
| 10 | 747 | 67 | 7 | 54 |
| 20 | 564 | 137 | 3 | 171 |
| 40 | 438 | 8 | 22 | 407 |

Across the three window sizes, 212 hits were kept only with `X` and 32 only with deletion. The sequences were sampled from the models and the windows placed by position, so the size of the effect on real tantan output is not measured.

## Added steps from the paper

Seven steps from the paper's Methods were not in the repository and are now in `commands.md`. Five were run on test data; ColabFold and the database-scale Foldseek searches were checked against the tools' options and code but not run.

| Step | Command and parameters | Test |
| --- | --- | --- |
| AntiFam filter | `hmmsearch --cut_ga AntiFam.hmm`; hits removed together with Pfam hits | 3 of 3 AntiFam seed proteins removed; 0 of 2 Pfam-sampled proteins |
| LAST, then DIAMOND | `lastdb -p -c` over all four reference files; `lastal -f BlastTab+`; ≥30% identity, ≥70% coverage of both proteins; DIAMOND only on proteins with no LAST hit | 90%, 60% and 50% identity queries removed by LAST; random protein kept as novel |
| Family profiles | families with ≥100 members from `cluDB.tsv` → `aligner.py` (MAFFT) → `trimmer.py` → `redundancy_removal.py` (hhfilter `-id 95 -cov 70`) | 40-member family: 40 of 40 rows through all three scripts (threshold lowered to 30 for the test) |
| ≥16 sequences | count of sequences left after hhfilter; one A3M per family, central sequence first | A3M with 40 sequences, first sequence gapless; accepted by `hhmake` |
| ColabFold | `colabfold_batch --model-type alphafold2_ptm --num-models 5 --num-recycle 3 --rank ptm`, no `--templates`; rank-1 model kept; HQ pTM ≥ 0.7, MQ 0.5–0.7 | not run (needs AlphaFold2 weights and a GPU); file names and JSON keys (`ptm`, `plddt`) checked in ColabFold 1.5.5 `batch.py` |
| Foldseek search | `easy-search` of HQ/MQ models against `CATH50` and `PDB`; models with no hit then against `Alphafold/UniProt50` | run on 1UBQ, 1AAR, 1CRN against a local database; CATH, PDB, AlphaFoldDB not downloaded |
| Superfamilies | `foldseek easy-cluster -c 0.8 --cov-mode 0 --tmscore-threshold 0.5` | 1UBQ grouped with both 1AAR chains; 1CRN alone |

One bug in the new A3M writer was found and fixed during testing: closing the file after every line made each write truncate it, leaving 1 sequence instead of 40.

## Assumptions and open points

Four details are not specified in the paper and are choices made in `commands.md`; each is labelled there.

- **LAST command line:** `lastdb -p -c` and `lastal -f BlastTab+`; the paper gives only the 30% / 70% thresholds.
- **16 effective sequences:** taken as the number of sequences after hhfilter. The paper calls the cutoff "the minimum number of sequences per MSA" and gives no formula; NMPFamsDB v1 does not define it either. HH-suite's NEFF is a different, smaller number (2.7 for a 40-sequence test family).
- **Model selection:** ColabFold `--rank ptm`, rank-1 model kept; the paper selects by "highest average pTM" and "best average pLDDT".
- **Structure databases:** Foldseek prebuilt `CATH50` and `PDB` (`pdb100`) instead of the paper's CATH v4.4 and PDB assemblies (March 2024).

One observation on the paper's structural hit rule: an alignment TM-score above 0.5 can come from a short alignment between unrelated proteins. In the test, ubiquitin (1UBQ) against crambin (1CRN) scored 0.546, with query- and target-normalised TM-scores of 0.17 and 0.27 and an E-value of 1.2. The rule is implemented as the paper states it; an added E-value or normalised TM-score condition would remove such hits.

Not in `commands.md`: the paper's taxonomy (Kraken2, MMseqs2 taxonomy, Whokaryote, EukRep, geNomad), gene-neighbourhood and biome analyses.

## Test data and file locations (local)

| What | Location |
| --- | --- |
| Filter, Pfam-join and trimming tests (fixes 1–8) | `~/Work/protein_families/test/` |
| X deletion vs X kept, DIAMOND | `~/Work/protein_families/test/xstrip/` |
| X deletion vs X kept, hmmsearch on Pfam 38.2 | `~/Work/protein_families/test/pfam/` |
| Script guard tests (MAFFT, hhfilter) | `~/Work/protein_families/test/guards/` |
| DIAMOND database naming test | `~/Work/protein_families/test/makedb/` |
| AntiFam, LAST → DIAMOND, profiles, A3M tests | `~/Work/protein_families/test/missing/` |
| Foldseek tests (1UBQ, 1AAR, 1CRN) | `~/Work/protein_families/test/foldseek/` |
| Pfam 38.2, AntiFam 8.0, NMPFamsDB2 data | `/Volumes/T7/ref/pfam/`, `/Volumes/T7/ref/antifam/`, `/Volumes/T7/ref/nmpfamsdb2/` |
