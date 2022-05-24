from os import replace
import pandas as pd


seg_l = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']
df_l = []
for seg in seg_l:
    file_path = './data/05_group_res/group_%s.txt' % seg
    df = pd.read_csv(file_path, sep='\t', names=['group', 'acc'])
    df['seg_name'] = seg
    df_l.append(df)

df_cat = pd.concat(df_l, axis=0, ignore_index=True)

df_trans = df_cat.pivot(index='acc', columns='seg_name', values='group')
df_trans = df_trans.applymap(int)


df_trans = df_trans.reset_index()

df_trans['Acc'] = df_trans['acc'].str.split('|').str.get(0)
df_trans['date'] = df_trans['acc'].str.split('|').str.get(1)

dt = pd.to_datetime(df_trans.date)


df2 = df_trans.sort_values(by='date')
print(df2)


# TODO 片段 group 编号重排序
df2 = df2.fillna(9)

for seg in seg_l:
    tmp = df2[seg].value_counts()
    group_trans = {}
    for idx, (k, v) in enumerate(tmp.iteritems()):
        print(f'For segment {seg}, id {k} to id {idx}')
        group_trans[int(k)] = idx

    df2[seg] = df2[seg].map(group_trans)
    print(group_trans)

# for seg in seg_l:
#     tmp = df2.groupby(seg)['date'].apply(min)
#     tmp = tmp.sort_values()
#     group_order =  tmp.index.to_list()
#     group_trans = dict(zip(range(len(group_order)), group_order))
#     df2[seg] = df2[seg].map(group_trans)


df2['group'] = df2[['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']].apply(lambda x: '_'.join(list((map(str, x.tolist())))), axis=1)
print(df2)

print(df2['group'].value_counts())

# df2.to_csv('genotype3.csv', sep='\t', index=False)
