"""Convert an IMG/VR nucleotide FASTA into the pipeline's input tables.

Writes
  scaffolds_info.tsv : scaffold ID \t length \t topology (circular | linear)
  protein.tsv        : protein ID \t sequence \t start \t end

IDs follow the IMG/M layout the pipeline expects (commands.md step 2 splits on "|"):
  scaffold = IMGVR|<uvig>        protein = IMGVR|<uvig>|<uvig>_<n>
<uvig> = first whitespace-separated token of the FASTA header, "|" replaced by "_";
with --first-field, only the part before the first "|" (use when that part alone is unique);
GVMAG genome bins share one UViG ID across contigs, so they get <uvig>__<contig>.

Genes are predicted with pyrodigal-gv (Prodigal with viral genetic-code models, meta mode).
Coordinates are 1-based, start < end on both strands; stop codon excluded from the protein.

Circular sequences are detected by a direct terminal repeat (DTR, >= --min-dtr bp; the
sequence end repeats its start, as in CheckV/geNomad). For these, the repeat copy is removed,
genes are predicted on the sequence joined to itself, and only complete genes starting in the
first copy are kept, so genes spanning the origin are predicted whole. A gene that wraps the
origin is written with start > end. The reported length is the length without the repeat.

usage: python imgvr_to_pipeline.py IMGVR_nucleotides.fna[.gz] OUTDIR [--threads N] [--min-len 500]
"""
import argparse
import gzip
import multiprocessing as mp
import os
import sys

import pyrodigal_gv

finder = None


def init_worker():
    global finder
    finder = pyrodigal_gv.ViralGeneFinder(meta=True)


def read_fasta(path, first_field=False):
    op = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                tok = line[1:].split()[0] if line[1:].strip() else ""
                if first_field:
                    f = tok.split("|")
                    # GVMAG = genome bin: one UViG ID over several contigs, the contig is field 3
                    name = f"{f[0]}__{f[2]}" if "_GVMAG-" in f[0] and len(f) > 2 else f[0]
                else:
                    name = tok.replace("|", "_")
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def dtr_length(seq, min_dtr):
    """Length of the longest direct terminal repeat (prefix == suffix), 0 if none >= min_dtr."""
    if len(seq) < 2 * min_dtr:
        return 0
    seed = seq[:min_dtr]
    p = seq.find(seed, len(seq) // 2)
    while p != -1:
        if seq[p:] == seq[:len(seq) - p]:
            return len(seq) - p
        p = seq.find(seed, p + 1)
    return 0


def predict(record, min_dtr):
    uvig, seq = record
    seq = seq.upper()
    dtr = dtr_length(seq, min_dtr) if min_dtr else 0
    rows = []
    if dtr:
        seq = seq[:len(seq) - dtr]
        L = len(seq)
        genes = [g for g in finder.find_genes((seq + seq).encode())
                 if g.begin <= L and not g.partial_begin and not g.partial_end]
    else:
        L = len(seq)
        genes = finder.find_genes(seq.encode())
    for n, gene in enumerate(genes, 1):
        prot = gene.translate(include_stop=False)
        end = gene.end - L if gene.end > L else gene.end
        rows.append(f"IMGVR|{uvig}|{uvig}_{n}\t{prot}\t{gene.begin}\t{end}\n")
    return uvig, L, "circular" if dtr else "linear", rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("outdir")
    ap.add_argument("--threads", type=int, default=os.cpu_count())
    ap.add_argument("--min-len", type=int, default=0,
                    help="skip gene prediction on sequences shorter than this (still listed in scaffolds_info.tsv)")
    ap.add_argument("--min-dtr", type=int, default=20,
                    help="minimum direct terminal repeat (bp) to call a sequence circular; 0 = treat all as linear")
    ap.add_argument("--first-field", action="store_true",
                    help='UViG ID = header text before the first "|" (fails on duplicate IDs)')
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    n_seq = n_prot = 0
    seen = set()
    with open(os.path.join(a.outdir, "scaffolds_info.tsv"), "w") as sc, \
         open(os.path.join(a.outdir, "protein.tsv"), "w") as pr, \
         mp.Pool(a.threads, initializer=init_worker) as pool:
        records = ((u, s) for u, s in read_fasta(a.fasta, a.first_field))
        for uvig, length, topology, rows in pool.imap(predict_or_skip(a.min_len, a.min_dtr), records, chunksize=64):
            if uvig in seen:
                sys.exit(f"duplicate UViG ID: {uvig}")
            seen.add(uvig)
            sc.write(f"IMGVR|{uvig}\t{length}\t{topology}\n")
            pr.writelines(rows)
            n_seq += 1
            n_prot += len(rows)
    print(f"sequences={n_seq} proteins={n_prot}", file=sys.stderr)


class predict_or_skip:
    def __init__(self, min_len, min_dtr):
        self.min_len = min_len
        self.min_dtr = min_dtr

    def __call__(self, record):
        if len(record[1]) < self.min_len:
            return record[0], len(record[1]), "linear", []
        return predict(record, self.min_dtr)


if __name__ == "__main__":
    main()
