from Bio import SeqIO
import pandas as pd


LPAIV_acc = ['EPI_ISL_661312']

# parse sequence file
seq_dict = {}
df = pd.DataFrame()
all_records = SeqIO.parse('data/all_nt.fasta', 'fasta')
for idx, record in enumerate(all_records):
    fileds = record.id.split('|')
    seg_acc = 'EPI' + fileds[0]
    isl_acc, isl_name, date, seg_names = fileds[2:]
    seq = str(record.seq)
    if isl_acc in LPAIV_acc:
        print(seg_acc, isl_acc, isl_name, date, seg_names, sep='\t')
        continue
    seq_dict[seg_acc] = seq
    df = df.append([[seg_acc, isl_acc, isl_name, date, seg_names]], ignore_index=True)


df.columns = ['seg_acc', 'isl_acc', 'isl_name', 'date', 'seg_name']
# df.to_csv('seq_info.csv', sep='\t', index=False)
df[['seg_acc', 'date']].to_csv('seq_date.csv', sep='\t', index=False)


# trans long table to wide table
df_de = df.drop_duplicates(subset=['isl_acc', 'seg_name'], keep='first', ignore_index=True)
df_trans = df_de.pivot(index='isl_acc', columns='seg_name', values='seg_acc')
isl_info = df[['isl_acc', 'isl_name', 'date']].drop_duplicates(keep='first', ignore_index=True)

# drop duplicate sequences
df_final = pd.merge(isl_info, df_trans, left_on='isl_acc', right_index=True)
df_tmp = df_final[['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']].applymap(lambda x: seq_dict.get(x))
df_tmp['isl_acc'] = df_final['isl_acc']
df_tmp['date'] = df_final['date']


# drop duplicate isolates with keep oldest isolate
df_final2 = df_tmp.sort_values(by=['isl_acc', 'date']).dropna(axis=0, how='any')
df_final3 = df_final2.drop_duplicates(subset=['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'], keep='first', ignore_index=True)
# print(df_final2[['isl_acc', 'date']])


df_final4 = df_final[df_final.isl_acc.isin(df_final3.isl_acc)].copy()



def parse_date(x):
    date_ = x.split('-')
    if len(date_) == 3 or len(date_) == 2:
        return int(date_[0]) + (int(date_[1]) - 1) / 12
    elif len(date_) == 1:
        return int(date_[0])

# trans date to num format
df_final4['date2'] = df_final4.date.map(parse_date)
# df_final4.to_csv('isl_info.csv', sep='\t', index=False)


# wirte sequences by segment
# for seg_name in ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']:
#     with open('data/by_seg/segment_%s.fasta' % seg_name, 'w') as f:
#         for acc in df_final4[seg_name]:
#             f.write('>%s\n%s\n' % (acc, seq_dict.get(acc)))
