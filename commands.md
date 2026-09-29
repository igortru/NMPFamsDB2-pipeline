# Novel Protein Discovery Pipeline

All data are retrieved from IMG/M.

## Input files

### Proteins

Proteins are stored as a TSV file (`protein.tsv`).

| Column | Description |
|---------|-------------|
| 1 | Protein ID |
| 2 | Protein sequence |
| 3 | Start position |
| 4 | End position |

### Scaffolds

Scaffold information is stored as a TSV file (`scaffolds_info.tsv`).

| Column | Description |
|---------|-------------|
| 1 | Scaffold ID |
| 2 | Scaffold length |

---

# Filtering

To minimize the inclusion of fragmented or incomplete protein predictions while retaining potentially functional proteins, a four step filtering strategy was applied. The filtering criteria included:

1. Protein length
2. Scaffold length
3. Gene position within the scaffold
4. Presence of low complexity regions

---

## Step 1. Filter by protein length (≥35 aa)

```bash
awk -F '\t' '{print $1"\t"$2"\t"length($2)"\t"$3"\t"$4}' protein.tsv \
| awk -F '\t' '$3>=35 {print $1"\t"$2"\t"$4"\t"$5}' \
> filtered_step1_proteins.tsv
```

---

## Step 2. Filter by scaffold length (≥500 bp)

Keep scaffolds that are at least 500 bp long.

```bash
awk -F '\t' '$2>=500 {print $1}' scaffolds_info.tsv \
| sort -k1,1 \
> sorted_scaffolds_morethan500.tsv
```

Prepare the protein table for joining.

```bash
sed 's/|/\t/g' filtered_step1_proteins.tsv \
| awk -F '\t' '{print $1"|"$2"\t"$3"\t"$4"\t"$5"\t"$6}' \
| sort -k1,1 \
> sorted_filtered_step1_proteins.tsv
```

Retain only proteins located on scaffolds ≥500 bp.

```bash
join -t $'\t' -1 1 -2 1 \
sorted_scaffolds_morethan500.tsv \
sorted_filtered_step1_proteins.tsv \
> filtered_step2_proteins.tsv
```

---

## Step 3. Remove genes near scaffold edges

Genes beginning within 10 nucleotides of the scaffold start or ending within 10 nucleotides of the scaffold end are removed.

`join` requires both inputs sorted on the join key; an unsorted `scaffolds_info.tsv` silently drops proteins.

```bash
LC_ALL=C sort -t $'\t' -k1,1 scaffolds_info.tsv > sorted_scaffolds_info.tsv
LC_ALL=C sort -t $'\t' -k1,1 filtered_step2_proteins.tsv \
| LC_ALL=C join -t $'\t' -1 1 -2 1 - sorted_scaffolds_info.tsv \
| awk -F '\t' '$4 > 10 && $5 < $6 - 10 {print $1"\t"$2"\t"$3}' \
> filtered_step3_proteins.tsv
```

Convert the filtered proteins to FASTA.

```bash
awk '{print ">"$1"|"$2"\n"$3}' filtered_step3_proteins.tsv \
> filtered_step3_proteins.fa
```

---

## Step 4. Remove proteins with extensive low complexity regions

Mask low complexity regions using **tantan**.

```bash
tantan -p -x X filtered_step3_proteins.fa \
> masked_low_complexity_region.fa
```

Convert the masked FASTA back to TSV.

```bash
awk '/^>/ {
    if (seq) print name"\t"seq;
    name=substr($0,2);
    seq="";
    next
}
{
    seq=seq $0
}
END {
    if (seq) print name"\t"seq
}' masked_low_complexity_region.fa \
> masked_low_complexity_region.tsv
```

Retain proteins with at least 35 unmasked amino acids. The masked sequence (with X) is kept as-is; stripping X joins the flanking segments and breaks alignments to reference proteins.

```bash
awk -F '\t' '{s=$2; n=gsub(/X/,"",s); if (length(s)>=35) print $1"\t"$2}' \
masked_low_complexity_region.tsv \
> final_filtered_proteins.tsv
```

Generate the final FASTA file.

```bash
awk '{print ">"$1"\n"$2}' final_filtered_proteins.tsv \
> final_filtered_proteins.fa
```

---

# Searching

A protein was classified as **novel** if it showed no detectable matches to any of:

* Pfam (v37)
* AntiFam (spurious-protein families)
* the IMG/M reference protein collection (LAST, then DIAMOND)

---

# Search against Pfam

To facilitate analysis of the initial filtered dataset (~6 billion proteins), the FASTA file was split into approximately 6,000 smaller files. Each chunk was searched independently against the Pfam HMM database using `hmmsearch`.

```bash
hmmsearch \
    --cut_tc \
    --cpu 45 \
    --domtblout results_chunkx.domtblout \
    -o results_chunkx.hmmout \
    Pfam.hmm \
    chunkx.fa
```

After all searches completed, the `*.domtblout` files were merged into a single results table. Every protein appearing in this table was considered a Pfam hit.

## domtblout format

| Column | Description |
|---------|-------------|
| 1 | Target name (protein ID) |
| 2 | Target accession |
| 3 | Target length |
| 4 | Query name (Pfam domain) |
| 5 | Query accession |
| 6 | Query length |
| 7 | E-value |
| 8 | Full sequence bit score |
| 9 | Full sequence bias |
| 10 | Domain number |
| 11 | Total domains |
| 12 | Conditional E-value |
| 13 | Independent E-value |
| 14 | Domain bit score |
| 15 | Domain bias |
| 16 | HMM start |
| 17 | HMM end |
| 18 | Alignment start |
| 19 | Alignment end |
| 20 | Envelope start |
| 21 | Envelope end |
| 22 | Alignment accuracy |
| 23 | Target description |

# Search against AntiFam

The same search is run against AntiFam. AntiFam models carry gathering thresholds (GA) only, so `--cut_ga` is used as recommended in the AntiFam release notes (`--cut_tc` fails: "TC bit thresholds unavailable").

```bash
hmmsearch \
    --cut_ga \
    --cpu 45 \
    --domtblout antifam_chunkx.domtblout \
    -o antifam_chunkx.hmmout \
    AntiFam.hmm \
    chunkx.fa
```

Merge the AntiFam `*.domtblout` files into `final_antifam.domtblout`.

## Collect Pfam and AntiFam hits

`domtblout` is space-delimited and contains `#` comment lines; split on whitespace and skip comments.

```bash
cat final_hmmresults.domtblout final_antifam.domtblout \
| grep -v '^#' \
| awk '{print $1}' \
| LC_ALL=C sort -u \
> protein_Pfam_hits.txt
```

## Remove Pfam and AntiFam hits

`final_filtered_proteins.tsv` is not sorted; `join` needs both inputs sorted in the same (C) order, otherwise Pfam hits leak through.

```bash
LC_ALL=C sort -t $'\t' -k1,1 final_filtered_proteins.tsv \
| LC_ALL=C join -t $'\t' -v 1 -1 1 -2 1 - protein_Pfam_hits.txt \
> pfam_novel_proteins.tsv
```

Generate FASTA.

```bash
awk '{print ">"$1"\n"$2}' pfam_novel_proteins.tsv \
> pfam_novel_proteins.fa
```

---

# Search against IMG/M reference proteins

Reference proteins are stored in four FASTA files:

* `bacteria.fa`
* `archaea.fa`
* `eukarya.fa`
* `viruses.fa`

Two sequential searches are used with the same hit criteria (≥30% identity, ≥70% query and subject coverage): LAST first, then DIAMOND on the proteins with no LAST hit.

## LAST search

The paper gives the thresholds but not the LAST command line; the options below (protein database `-p`, BLAST-like tabular output with query/subject lengths) are this repository's choice. `BlastTab+` columns 1–14 are the same as the DIAMOND `--outfmt 6` columns used below, so the same hit filter applies.

```bash
lastdb -p -c -P 45 refDB bacteria.fa archaea.fa eukarya.fa viruses.fa

lastal -P 45 -f BlastTab+ refDB pfam_novel_proteins.fa \
| grep -v '^#' \
> last_matches.tsv

awk -F '\t' '{
    if ($3>=30 &&
        (($8-$7+1)/$13)>=0.7 &&
        (($10-$9+1)/$14)>=0.7)
    print $1
}' last_matches.tsv \
| LC_ALL=C sort -u \
> last_hits_ids

LC_ALL=C join -t $'\t' -v 1 -1 1 -2 1 \
pfam_novel_proteins.tsv \
last_hits_ids \
> last_novel_proteins.tsv

awk '{print ">"$1"\n"$2}' last_novel_proteins.tsv \
> last_novel_proteins.fa
```

## Create DIAMOND databases

One database per reference file (a shared name would overwrite the previous one).

```bash
for db in bacteria archaea eukarya viruses; do
    diamond makedb \
        --in ${db}.fa \
        -d ${db}
done
```

## Search proteins

Search each database and merge the results.

```bash
for db in bacteria archaea eukarya viruses; do
    diamond blastp \
        -d ${db}.dmnd \
        -q last_novel_proteins.fa \
        -o matches_${db}.m8 \
        --outfmt 6 \
        qseqid sseqid pident length mismatch gapopen \
        qstart qend sstart send evalue bitscore qlen slen
done
cat matches_bacteria.m8 matches_archaea.m8 matches_eukarya.m8 matches_viruses.m8 \
> matches.m8
```

### Hit criteria

A protein is considered a reference match if it satisfies all of the following:

* Percent identity ≥30%
* Query coverage ≥70%
* Subject coverage ≥70%

Collect matching protein IDs.

```bash
awk '{
    if ($3>=30 &&
        (($8-$7+1)/$13)>=0.7 &&
        (($10-$9+1)/$14)>=0.7)
    print $1
}' matches.m8 \
| LC_ALL=C sort -u \
> reference_hits_ids
```

Remove reference hits (`last_novel_proteins.tsv` is already C-sorted).

```bash
LC_ALL=C join -t $'\t' -v 1 -1 1 -2 1 \
last_novel_proteins.tsv \
reference_hits_ids \
> novel_proteins.tsv
```

Generate the final FASTA.

```bash
awk '{print ">"$1"\n"$2}' novel_proteins.tsv \
> novel_proteins.fa
```

---

# Clustering

Novel proteins are clustered using **MMseqs2 Linclust**.

Create the sequence database.

```bash
mmseqs createdb novel_proteins.fa seqDB
```

Run Linclust.

```bash
mmseqs linclust seqDB cluDB tmp \
    --min-seq-id 0.3 \
    -c 0.8 \
    --cov-mode 0
```

Export cluster assignments.

```bash
mmseqs createtsv seqDB seqDB cluDB cluDB.tsv
```

---

# Profile generation (families with ≥100 members)

Build a family table (`cluster \t member \t sequence`) for clusters with at least 100 members. `cluDB.tsv` has `representative \t member`.

```bash
awk -F '\t' '{n[$1]++} END{for (c in n) if (n[c]>=100) print c}' cluDB.tsv \
| LC_ALL=C sort > families100.txt

LC_ALL=C sort -t $'\t' -k1,1 cluDB.tsv \
| LC_ALL=C join -t $'\t' - families100.txt \
| LC_ALL=C sort -t $'\t' -k2,2 \
| LC_ALL=C join -t $'\t' -1 2 -2 1 -o 1.1,1.2,2.2 - novel_proteins.tsv \
| LC_ALL=C sort -t $'\t' -k1,1 \
> families100.tsv
```

Split by family into chunks (all rows of a family in the same chunk), then run the three scripts in order. Each takes a text file listing chunk paths and writes TSV to stdout.

```bash
python aligner.py chunks.txt > aligned.tsv              # MAFFT per family
# re-chunk aligned.tsv by family -> aligned_chunks.txt
python trimmer.py aligned_chunks.txt > trimmed.tsv      # central-sequence column trimming
# re-chunk trimmed.tsv by family -> trimmed_chunks.txt
python redundancy_removal.py trimmed_chunks.txt > filtered.tsv   # hhfilter -id 95 -cov 70
```

## Keep families with ≥16 sequences after filtering

The paper retains "families with at least 16 effective sequences after filtering", described as "the minimum number of sequences per MSA". Here this is taken as the number of sequences left after `hhfilter`; the paper does not give a formula (HH-suite's `hhmake` NEFF is a different, smaller number). Each family is written as an A3M file (the last column of `filtered.tsv`), first sequence = central sequence.

```bash
mkdir -p a3m
awk -F '\t' 'NF>=5 {n[$1]++} END{for (c in n) if (n[c]>=16) print c}' filtered.tsv > families_neff16.txt
awk -F '\t' 'NR==FNR {keep[$1]=1; next}
             NF>=5 && ($1 in keep) {f="a3m/"$1".a3m"; if (f!=prev) {if (prev) close(prev); prev=f}
                                    print ">"$2"\n"$NF > f}' \
    families_neff16.txt filtered.tsv
```

---

# Structure prediction (ColabFold)

ColabFold implementation of AlphaFold2 in de novo mode (no `--templates`), family A3M as input MSA, five models, three recycles. Models are ranked by pTM; the rank-1 model is kept (the paper selects the model with the highest average pTM and best average pLDDT).

```bash
colabfold_batch \
    --model-type alphafold2_ptm \
    --num-models 5 \
    --num-recycle 3 \
    --rank ptm \
    a3m/ colabfold_out/
```

Collect the rank-1 pTM per family and classify: HQ pTM ≥ 0.7, MQ 0.5 ≤ pTM < 0.7, LQ pTM < 0.5. Only HQ and MQ models go to the structural search.

```bash
for j in colabfold_out/*_scores_rank_001_*.json; do
    fam=$(basename "$j" | sed 's/_scores_rank_001_.*//')
    python -c "import json,sys; d=json.load(open(sys.argv[1])); print(sys.argv[2], d['ptm'], sum(d['plddt'])/len(d['plddt']))" "$j" "$fam"
done > model_scores.tsv

awk '$2>=0.5 {print $1}' model_scores.tsv > models_hq_mq.txt

mkdir -p models
while read -r fam; do
    cp colabfold_out/${fam}_unrelaxed_rank_001_*.pdb models/${fam}.pdb
done < models_hq_mq.txt
```

---

# Structural homology search (Foldseek)

HQ/MQ models (rank-1 PDB files in `models/`) are searched against CATH and PDB. Foldseek's prebuilt `CATH50` and `PDB` databases are used here; the paper used CATH v4.4 and PDB assemblies (March 2024), so versions differ.

```bash
foldseek databases CATH50 cathDB tmp
foldseek databases PDB pdbDB tmp

for db in cathDB pdbDB; do
    foldseek easy-search models/ ${db} ${db}_hits.tsv tmp \
        --format-output query,target,alntmscore,qtmscore,ttmscore,qlen,tlen,evalue
done
```

Hit rule as described in the paper: alignment TM-score > 0.5; otherwise, if the query is shorter than the target, query-normalised TM-score ≥ 0.5; if the target is shorter, target-normalised TM-score ≥ 0.5.

```bash
cat cathDB_hits.tsv pdbDB_hits.tsv \
| awk -F '\t' '$3>0.5 || ($6<$7 && $4>=0.5) || ($7<$6 && $5>=0.5) {print $1}' \
| LC_ALL=C sort -u > structural_hits.txt
```

Models without a CATH/PDB hit are searched against AlphaFoldDB the same way.

```bash
foldseek databases Alphafold/UniProt50 afdb tmp
ls models/ | sed 's/\.pdb$//' | LC_ALL=C sort \
| LC_ALL=C comm -23 - structural_hits.txt > no_exp_hit.txt
# copy the models listed in no_exp_hit.txt to models_noexp/, then:
foldseek easy-search models_noexp/ afdb afdb_hits.tsv tmp \
    --format-output query,target,alntmscore,qtmscore,ttmscore,qlen,tlen,evalue
```

## Superfamilies

Foldseek clustering in bidirectional mode, 80% coverage, TM-score 0.5.

```bash
foldseek easy-cluster models/ superfamilies tmp \
    -c 0.8 \
    --cov-mode 0 \
    --tmscore-threshold 0.5
```
