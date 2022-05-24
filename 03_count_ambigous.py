from Bio import SeqIO
import os
import glob


paths = glob.glob('data/cdhit/segment_*.de.fasta')

tmp = paths[1]
print(tmp)


def cal_content(record):
    # length A C G T N other_IUPAC - ? invalid_nucleotides
    strain_name = record.id.split('|')[0]
    seq_len = len(record)
    record.seq = record.seq.upper()
    count_a = record.seq.count('A')
    count_c = record.seq.count('C')
    count_g = record.seq.count('G')
    count_t = record.seq.count('T')
    other_IUPAC = ['RYSWKMBDHV']
    count_other = sum([record.seq.count(x) for x in other_IUPAC])
    count_N = record.seq.count('N')
    count_gap = record.seq.count('-')
    count_q = record.seq.count('?')
    count_invalid = seq_len - count_a - count_c - count_g - count_t - count_other - count_N - count_gap - count_q
    return strain_name, seq_len, count_a, count_c, count_g, count_t, count_other, count_N, count_gap, count_q, count_invalid


for record in SeqIO.parse(tmp, 'fasta'):
    res = cal_content(record)
    if res[9] != 0 or res[10] != 0:
        print(*res, sep='\t')
