from Bio import SeqIO
import pandas as pd
import numpy as np
import random
import subprocess
import os
import re
from glob import glob
from io import StringIO
from Bio.Blast import NCBIXML


SEGMENT = 'MP'


out_file = 'data/blast_res/%s.xml' % SEGMENT


blast_res_l = []
with open(out_file) as f:
    blast_records = NCBIXML.parse(f)
    res = []
    for idx, record in enumerate(blast_records):
        print('************* QUERY %s INFORMATION *************' % idx)
        print('Query name: ', record.query)
        print('Query length: ', record.query_length)
        
        print('********** Result Description **********')
        for i, (des_, align_) in enumerate(zip(record.descriptions, record.alignments)):
            isl_name, isl_acc, date, seg_name = align_.hit_def.split('#')
            seg_acc = align_.hit_id
            res.append([record.query, seg_acc, isl_acc, isl_name, date, seg_name])
            if i < 5:
                print(des_.score, seg_acc, isl_acc, date, seg_name, isl_name, sep='\t')
        else:
            print('...\n')

        blast_res_l.append(pd.DataFrame(res, columns=['q_acc', 'seg_acc', 'isl_acc', 'isl_name', 'date', 'seg_name']))


blast_res = pd.concat(blast_res_l, axis=0, ignore_index=True)
print(blast_res)

date_range = sorted(blast_res['date'].unique())
print(date_range[:5], date_range[-5:])


blast_summary = blast_res.groupby(['q_acc', 'isl_name']).size() / len(blast_res_l)
print(blast_summary)

unique_sbj = blast_res.seg_acc.unique()
unique_que = blast_res.q_acc.unique()
union_ = np.union1d(unique_que, unique_sbj)

print(len(unique_que), len(unique_sbj), len(union_))
 

# 提取 blast hits 序列
hit_accs = union_.tolist()
cmd = "blastdbcmd -db /home/zeng/blastdb/fludb/fludb2 -entry '%s'" % ','.join(hit_accs)
blastdbcmd = subprocess.run(cmd, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env={'PATH':'/home/zeng/opt/ncbi-blast-2.11.0+/bin'})
# print(blastdbcmd.returncode)
# print(blastdbcmd.stderr)
# print(blastdbcmd.stdout)

hits_handle = StringIO()
hits_handle.write(blastdbcmd.stdout)
hits_handle.seek(0)
hits_seqs = list(SeqIO.parse(hits_handle, 'fasta'))

print(hits_seqs[0])

# SeqIO.write(hits_seqs, 'data/blast_res/blast_res_%s.fasta' % SEGMENT, 'fasta')

# with open(out_file + '_blast_summary.txt', 'w') as f:
#     f.write("%s\n" % '\t'.join(clade_strain_names_sampled))
#     f.write("%s\n" % '\t'.join(hit_names))
#     for k, v  in blast_summary.iteritems():
#         f.write('%s\t%s\n' % (k, v))
