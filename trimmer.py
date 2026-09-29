from sys import argv, stderr, stdout
from collections import OrderedDict as odict
import pandas as pd
import multiprocessing as mp
from sklearn.metrics import pairwise_distances
import numpy as np


# GLOBALS TO EDIT
number_of_files = 18
# False = use the central sequence (min mean Hamming distance), as described in the paper
first_is_representative = False


# aa globals definitions
alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
states = len(alpha_1)
alpha_3 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
           'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'GAP']
aa_1_N = {a: n for n, a in enumerate(alpha_1)}
aa_3_N = {a: n for n, a in enumerate(alpha_3)}
aa_N_1 = {n: a for n, a in enumerate(alpha_1)}
aa_1_3 = {a: b for a, b in zip(alpha_1, alpha_3)}
aa_3_1 = {b: a for a, b in zip(alpha_1, alpha_3)}


def parse_alignments_from_nmpf_tsv(tsv_input: str):
    """
    :param tsv_input:
    :return:
    """
    alignments = []
    # we expect the aligned tsv to have the following column format:
    # 0: cluster ID, 1: header, 2: sequence, 3: aligned sequence
    tsv = pd.read_csv(tsv_input, sep="\t", header=None)
    # get the cluster IDs from the data frame's column 0
    clusters = sorted(list(set(tsv[0])))
    # for each cluster, isolate the corresponding lines from the data frame and construct a headers and a sequences
    # list, then append them to a list of alignments as a
    # dictionary with the form: {headers: headers, sequences: sequences}
    for c in clusters:
        lines = tsv[tsv[0] == c]
        headers = []
        sequences = []
        for i in lines.index:
            line = lines.loc[i]
            header = "%s|%s" % (line[0], line[1])
            seq = line[3]
            headers.append(header)
            sequences.append(seq)
        alignments.append({"headers": headers, "sequences": sequences})
    return alignments


def aa_to_n(x: list):
    """
    :param x:
    :return:
    """
    x = np.array(x)
    if x.ndim == 0:
        x = x[None]
    return [[aa_1_N.get(a, states - 1) for a in y] for y in x]


def detect_gaps(seq: str):
    positions = []
    for i in range(len(seq)):
        if seq[i] == "-" or seq[i] == ".":
            positions.append(i)
    return positions


def remove_positions(seq: str, positions: list):
    new_seq = ""
    for i in range(len(seq)):
        if i not in positions:
            new_seq += seq[i]
    return new_seq


def identify_representative_sequence(alnobj: dict):
    sequences = alnobj["sequences"]
    idx_of_center_seq = pairwise_distances(np.array(aa_to_n(sequences)), metric="hamming").mean(0).argmin()
    return idx_of_center_seq


def trim_msa_based_on_representative(alnobj: dict, idx_of_center_seq: int = 0):
    """

    :param alnobj: the alignment object
    :param idx_of_center_seq: numerical identifier of the centroid sequence
    :return: trimmed_msa: a dictionary containing the trimmed_msa, where keys are headers, values are sequences
    """
    sequences = alnobj["sequences"]
    headers = alnobj["headers"]
    center_seq_header = headers[idx_of_center_seq]
    center_seq_sequence = sequences[idx_of_center_seq]
    pos_to_remove = detect_gaps(center_seq_sequence)
    center_seq_sequence = remove_positions(center_seq_sequence, pos_to_remove)

    trimmed_msa = odict({})
    trimmed_msa[center_seq_header] = center_seq_sequence

    for i in range(len(sequences)):
        if i != idx_of_center_seq:
            seq = remove_positions(sequences[i], pos_to_remove)
            head = headers[i]
            trimmed_msa[head] = seq

    return trimmed_msa


def convert_trimmed_msa_to_tsv(tsv_input: str, trimmed_msa: odict):
    tsv = pd.read_csv(tsv_input, sep="\t", header=None)

    for header, trimmed_seq in trimmed_msa.items():
        itms = header.split("|")
        cluster = itms[0]
        seq_header = "|".join(itms[1:])
        line = tsv[(tsv[0] == cluster) & (tsv[1] == seq_header)].reset_index()
        if len(line) != 0:
            row = line.loc[0]
            str_line = f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t{trimmed_seq}\n"
            stdout.write(str_line)


def process_msa(alnobj, first_is_representative):
    if first_is_representative is True:
        representative_index = 0
    else:
        representative_index = identify_representative_sequence(alnobj)
    trimmed = trim_msa_based_on_representative(alnobj, idx_of_center_seq=representative_index)
    return trimmed


def process_tsv(tsv_input, first_is_representative=True):
    try:
        msas = parse_alignments_from_nmpf_tsv(tsv_input)
        for alnobj in msas:
            print(f"Trimming alignment for {alnobj['headers'][0].split('|')[0]}", file=stderr)
            trimmed = process_msa(alnobj, first_is_representative=first_is_representative)
            convert_trimmed_msa_to_tsv(tsv_input, trimmed)
    except:
        print(f"Problem with chunk {tsv_input}", file=stderr)


# if __name__ == "__main__":
#     script, tsv_with_alignments = argv
#     process_tsv(tsv_with_alignments.rstrip())


def parallelize(chunks, number_of_files):
    print(chunks, file=stderr)
    pool = mp.Pool(number_of_files)
    proc = [pool.apply_async(process_tsv, args=(chunk, first_is_representative,)) for chunk in chunks]
    result = [p.get() for p in proc]
    for r in result:
        print(r, file=stderr)


# guard required for multiprocessing "spawn" start method (macOS/Windows default)
if __name__ == "__main__":
    script, list_of_chunks = argv

    with open(list_of_chunks, "r") as fl:
        chunks = []
        for i in fl:
            if i.rstrip() != "":
                chunks.append(i.rstrip())
    fl.close()



    for i in range(0, len(chunks), number_of_files):
        chunks_subset = chunks[i:i+number_of_files]
        parallelize(chunks_subset, number_of_files)
