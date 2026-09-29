"""Convert an IMG/VR nucleotide FASTA into the pipeline's input tables.

Writes
  scaffolds_info.tsv : scaffold ID \t length
  protein.tsv        : protein ID \t sequence \t start \t end

IDs follow the IMG/M layout the pipeline expects (commands.md step 2 splits on "|"):
  scaffold = IMGVR|<uvig>        protein = IMGVR|<uvig>|<uvig>_<n>
<uvig> = first whitespace-separated token of the FASTA header, "|" replaced by "_";
with --first-field, only the part before the first "|" (use when that part alone is unique).

Genes are predicted with pyrodigal-gv (Prodigal with viral genetic-code models, meta mode).
Coordinates are 1-based, start < end on both strands; stop codon excluded from the protein.

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
                name = tok.split("|")[0] if first_field else tok.replace("|", "_")
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def predict(record):
    uvig, seq = record
    rows = []
    for n, gene in enumerate(finder.find_genes(seq.encode()), 1):
        prot = gene.translate(include_stop=False)
        rows.append(f"IMGVR|{uvig}|{uvig}_{n}\t{prot}\t{gene.begin}\t{gene.end}\n")
    return uvig, len(seq), rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("outdir")
    ap.add_argument("--threads", type=int, default=os.cpu_count())
    ap.add_argument("--min-len", type=int, default=0,
                    help="skip gene prediction on sequences shorter than this (still listed in scaffolds_info.tsv)")
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
        for uvig, length, rows in pool.imap(predict_or_skip(a.min_len), records, chunksize=64):
            if uvig in seen:
                sys.exit(f"duplicate UViG ID: {uvig}")
            seen.add(uvig)
            sc.write(f"IMGVR|{uvig}\t{length}\n")
            pr.writelines(rows)
            n_seq += 1
            n_prot += len(rows)
    print(f"sequences={n_seq} proteins={n_prot}", file=sys.stderr)


class predict_or_skip:
    def __init__(self, min_len):
        self.min_len = min_len

    def __call__(self, record):
        if len(record[1]) < self.min_len:
            return record[0], len(record[1]), []
        return predict(record)


if __name__ == "__main__":
    main()
