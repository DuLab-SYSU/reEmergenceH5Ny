import pandas as pd
import numpy as np
import glob
import os

meta_dir = os.path.expanduser('../data/')
files = glob.glob(meta_dir + 'gisaid_epi*.xls')


df_l = [pd.read_excel(file_path) for file_path in files]
df = pd.concat(df_l, axis=0)

seg_name_l = ['PB2 Segment_Id', 'PB1 Segment_Id', 'PA Segment_Id', 'HA Segment_Id', 'NP Segment_Id', 'NA Segment_Id', 'MP Segment_Id', 'NS Segment_Id']

df1 = df[['Isolate_Id', 'Isolate_Name', 'Subtype', 'Host', 'Location', 'Collection_Date', 'Domestic_Status']]
t = df1.Collection_Date.str.split('-', expand=True)
# print(t.columns)
t.columns = ['year', 'month', 'day']


def p(x):
    try:
        return x.split('|')[0]
    except AttributeError:
        return np.nan

df2 = df[seg_name_l].applymap(p)
print(df2)


df3 = pd.concat([df1, t, df2], axis=1)
print(df3)

df3.to_csv('data/all_aiv_info.csv', sep='\t', index=False)
