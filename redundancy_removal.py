import subprocess
import pandas as pd
from sys import argv, stdout, stderr
import multiprocessing as mp

def run_hhfilter(fasta_file):
#output file
    filtered_file = fasta_file.replace(".fasta", "_filtered.a3m")
#hhfilter command
    command = "hhfilter -i "+fasta_file+" -o "+filtered_file+" -M first -id 95 -cov 70"
#write output 
    subprocess.call(command, shell=True, stderr=subprocess.DEVNULL)

    return filtered_file



def get_chunk_clusters(chunk):
    csv = pd.read_csv(chunk, sep="\t", header=None, usecols=[0])
    clusters = sorted(list(set(csv[0])))
    return clusters


def process_tsv(chunk, cluster, output_fasta):
    with open(chunk, "r") as cf, open(output_fasta, "w") as out:
        for line in cf:
            columns = line.rstrip().split("\t")
            if columns[0] == cluster:
                header = f"{columns[0]}|{columns[1]}"
                seq = columns[3]
                fasta = f">{header}\n{seq}\n"
                out.write(fasta)
                #print(fasta)
    out.close()


def parse_fasta(fasta_file):
    fas_fl = open(fasta_file, "r")
    fastas = fas_fl.read().split(">")
    fas_fl.close()
    fastas.pop(0)
    fasta_dict = {}
    for fas in fastas:
        lines = fas.split("\n")
        title = lines[0]
        seq = "".join(lines[1:])
        fasta_dict[title] = seq
    return fasta_dict


def append_alignment_to_tsv(input_tsv, aligned_sequences):
    inp_fl = open(input_tsv, "r")
    rows = []
    for line in inp_fl:
        columns = line.rstrip().split("\t")
        header = f"{columns[0]}|{columns[1]}"
        if header not in aligned_sequences.keys():
            aligned_seq = None
        else:
            aligned_seq = aligned_sequences[header]
        if aligned_seq is not None:
            rows.append(line.rstrip() + "\t" + aligned_seq)
    return rows


def process_chunk(chunk):
    clusters = get_chunk_clusters(chunk)
    #print(clusters, file=stderr)
    cluster_text = []
    #output_tsv = chunk + "_filtered.tsv"
    #flout = open(output_tsv, "w")
    for cluster in clusters:
        print("Running cluster " + cluster + " from chunk " + chunk, file=stderr)
        try:
            unaligned_fasta = f"cluster_{cluster}.fasta"
            aligned_fasta = f"cluster_{cluster}_filtered.a3m"
            process_tsv(chunk, cluster, unaligned_fasta)
            print(unaligned_fasta, file=stderr)
            run_hhfilter(unaligned_fasta)
            aligned_sequences = parse_fasta(aligned_fasta)
            rows = append_alignment_to_tsv(chunk, aligned_sequences)
            cluster_text.append("\n".join(rows))
            subprocess.call(f"rm cluster_{cluster}.fasta cluster_{cluster}_filtered.a3m", shell=True)
            stdout.write("\n".join(rows)+"\n")
            #flout.write("\n".join(cluster_text)+"\n")
        except:
            print("Error in cluster " + cluster + " from chunk " + chunk)
    #flout.close()


def parallelize(chunks, number_of_files):

    pool = mp.Pool(number_of_files)

    proc = [pool.apply_async(process_chunk, (chunk,)) for chunk in chunks]
    result = [p.get() for p in proc]
    for r in result:
        print(r, file=stderr)


# guard required for multiprocessing "spawn" start method (macOS/Windows default)
if __name__ == "__main__":
    script, list_of_chunks = argv

    with open(list_of_chunks, "r") as fl:
        chunks = [i.rstrip() for i in fl]
    fl.close()


    number_of_files = 40

    for i in range(0, len(chunks), number_of_files):
        chunks_subset = chunks[i:i+number_of_files]
        parallelize(chunks_subset, number_of_files)
