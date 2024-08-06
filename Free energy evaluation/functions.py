import pandas as pd
import numpy as np
from Bio import SeqIO
import re

def random_sequence_position_generator(length, number, bed_file):
    random_positions = []
    for index, row in bed_file.iterrows():
        chromosome = row.iloc[0]
        start = row.iloc[1]
        end = row.iloc[2]
        region = row[3]
        for i in range(number):
            random_start = np.random.randint(start, end - length)
            random_end = random_start + length
            random_positions.append([chromosome, random_start, random_end, region])
    random_df = pd.DataFrame(random_positions, columns=["Chromosome", "Start", "End", "Region"])
    return random_df

def get_sequences(positions_df, fasta_file):
    sequences = {}
    DNA_seqs = []
    for record in SeqIO.parse(fasta_file,'fasta'):
        sequences[record.id.lower()] = str(record.seq)
    for index,row in positions_df.iterrows():
        chrom = row.iloc[0].lower()
        begin = row.iloc[1] - 1
        end = row.iloc[2]
        if chrom in sequences:
            DNA_seqs.append(sequences[chrom][begin:end])
    positions_df['Sequence'] = DNA_seqs
    return positions_df

def write_fasta(file_name, sequences):
    with open(file_name, "w") as fasta_file:
        for seq_id, sequence in sequences.items():
            fasta_file.write(f">{seq_id}\n{sequence}\n")

def read_output(content):
    mfe_pattern = r'\((-?\d+\.\d+)\)'
    fete_pattern = r'\[(-?\d+\.\d+)\]'
    ed_pattern = r'ensemble diversity\s+(-?\d+\.\d+)'
    mfe = re.findall(mfe_pattern, content)
    fete = re.findall(fete_pattern, content)
    ed = re.findall(ed_pattern, content)
    array = [mfe[0], fete[0], ed[0]]
    return array

