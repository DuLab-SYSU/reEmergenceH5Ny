import os
import glob
import numpy as np
import pandas as pd
import subprocess
from io import StringIO
from Bio import SeqIO
import random
import matplotlib.pyplot as plt
from scipy.stats import gamma


#* crosstab for dataset1
df = pd.read_csv('./genotype2.csv', sep='\t')
df.index = pd.to_datetime(df.date)
df = df.to_period('M')

df2 = df.groupby([df.index, 'group']).size().unstack(1, fill_value=0)

group_order = ['1_1_1_1_1_1_1_1', '0_0_0_0_0_0_0_0', '0_0_5_0_0_0_0_0', '4_2_3_0_3_2_0_0', '2_2_1_1_2_1_1_1',
               '3_2_2_1_2_1_1_2', '1_1_1_1_2_1_1_1',
               '3_5_1_1_2_1_1_1', '5_4_4_0_4_0_3_3',
               '6_3_6_0_5_0_2_4', '3_1_2_1_2_1_1_1']


#* down-sampling by a exponential distribution
def sample_by_exp(x):
    if x.sum() <= 30:
        return x

    tmp1 = x[x==0]
    tmp2 = x[x!=0]
    rv = gamma(1, scale=1.)
    sample_p = rv.pdf(np.linspace(0.1, 3, len(tmp2)))
    tmp2 = (tmp2 * sample_p).round()
    return pd.concat([tmp1, tmp2]).sort_index()

df_sampled = df2.apply(sample_by_exp)


tmp = pd.crosstab(index=df.index, columns=df.group, values=df.Acc, aggfunc=lambda x: x.to_list())
res = []
for column_ in df2.columns:
    for index_ in df2.index:
        acc_l = tmp.loc[index_, column_]
        sample_size = round(df_sampled.loc[index_, column_])
        if sample_size:
            random.seed(123456)
            sampling = random.sample(acc_l, sample_size)
            res.extend(sampling)


print("Total number of sampling:", len(res))
print(' '.join(res[0:10]) + ' ... ' + ' '.join(res[-10:]))

with open('sampling_accs_of_dataset1.txt', 'w') as f:
    f.write('\n'.join(res))


# fig, (ax1, ax2) = plt.subplots(2, 1)
# df2[group_order].plot(kind='area', colormap='viridis', stacked=False, ax=ax1, legend=False)
# ax1.legend(loc='upper left', bbox_to_anchor=(1.05, 1), fontsize=8)
# df_sampled[group_order].plot(kind='area', colormap='viridis', stacked=False, ax=ax2, legend=False)
# plt.show()



# SEGMENT = 'NS'

# tmp = df_final[SEGMENT].dropna().tolist()
# query_str = ','.join(tmp)

# cmd = "blastdbcmd -db /home/zeng/blastdb/fludb/fludb2 -entry '%s'" % query_str
# blastdbcmd = subprocess.run(cmd, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env={'PATH':'/home/zeng/opt/ncbi-blast-2.11.0+/bin'})

# hits_handle = StringIO()
# hits_handle.write(blastdbcmd.stdout)
# hits_handle.seek(0)
# hits_seqs = list(SeqIO.parse(hits_handle, 'fasta'))
# for record in hits_seqs:
#     seq_acc = record.id
#     attr = record.description.split(' ', 1)[1]
#     isl_name, isl_acc, _, seg_name = attr.split('#')
#     host = df_final.loc[isl_acc, 'host2']
#     subregion = df_final.loc[isl_acc, 'Subregion']
#     subtype = df_final.loc[isl_acc, 'Subtype']
#     date = df_final.loc[isl_acc, 'Collection_Date']
#     lat = df_final.loc[isl_acc, 'Latitude']
#     lon = df_final.loc[isl_acc, 'Longitude']
#     record.id = '|'.join([seq_acc, isl_acc, isl_name, date, host, subregion, subtype, str(lat), str(lon)])
#     record.name = ''
#     record.description = ''
#     # break

# SeqIO.write(hits_seqs, 'data/ext_seqs/seqs_ext_%s.fasta' % SEGMENT, 'fasta')


# file_dir = 'data/ext_seqs/'
# files = glob.glob(file_dir + '*.align.fasta')
# print(files)

# file = files[3]
# records = list(SeqIO.parse(file, 'fasta'))
# for record in records:
#     seqacc, islacc, islname, date, host, location, subtype, lat, lon = record.description.replace(' ', '_').split('|')
#     record.id = '|'.join([seqacc, date, host, location, islacc, islname, subtype, lat, lon])
#     record.name = ''
#     record.description = ''

# SeqIO.write(records, file, 'fasta')
