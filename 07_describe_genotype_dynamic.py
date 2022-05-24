import os
import glob
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import cm
from matplotlib.colors import to_hex, ListedColormap
import seaborn as sns
from sklearn.neighbors import KernelDensity
from datetime import timedelta



df = pd.read_csv('./genotype.csv', sep='\t', index_col=0)
df.index = pd.to_datetime(df.date)
df = df.to_period('M')

df2 = df.groupby([df.index, 'group']).size().unstack(1, fill_value=0)
# print(df2)


group_time_order = df.groupby(['group']).apply(lambda df_: df_.index.value_counts().index[0]).sort_values().index.values

group_order = ['1_1_1_1_1_1_1_1', '0_0_0_0_0_0_0_0', '4_2_3_0_3_2_0_0', '2_2_1_1_2_1_1_1', '3_2_2_1_2_1_1_2', '1_1_1_1_2_1_1_1',
               '0_0_5_0_0_0_0_0', '3_5_1_1_2_1_1_1', '5_4_4_0_4_0_3_3', '6_3_6_0_5_0_2_4', '3_1_2_1_2_1_1_1']

minor_group = ['0_0_5_0_0_0_0_0', '3_5_1_1_2_1_1_1', '5_4_4_0_4_0_3_3', '6_3_6_0_5_0_2_4', '3_1_2_1_2_1_1_1']


color_list = ['#BC3C29', '#FFDC91', '#0072B5', '#FFDC91', '#20854E', '#FFDC91', '#E18727', '#FFDC91', '#7876B1', '#FFDC91',  '#EE4C97']


color_map3 = dict(zip(group_order, color_list))
cmap =ListedColormap(color_list)


# TODO 统一x坐标轴

fig, (ax, ax1, ax2) = plt.subplots(3, 1, figsize=(12, 6))

ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
plt.setp(ax.get_xticklabels(), rotation=15)


for i, group in enumerate(group_order):
    tmp = df[df['group'] == group]['date']
    s = pd.to_datetime(tmp.min())
    e = pd.to_datetime(tmp.max()) + timedelta(7)

    ax.hlines(i, s, e, colors=color_map3[group], linewidth=4, alpha=0.8)
    ax.text(e, i, group, color=color_map3[group], fontsize=8, alpha=0.8)

ax.grid(axis='x')
ax.yaxis.set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

time_idx = sorted(df.index.unique())
groups_len = len(time_idx)
X_plot = np.linspace(0, groups_len, 1000)[:, np.newaxis]
time_obs = df.index.map(lambda x: time_idx.index(x))


# ************************ KDE scheme *******************************
# genotype = '4_2_3_1_3_2_1_2'
# X = time_obs[df['group'] == genotype].values[:, np.newaxis]
# kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(X)

# log_dens = kde.score_samples(X_plot)
# exp_dens = np.exp(log_dens)

# fig = plt.figure()
# ax = plt.subplot(211)

# ax.hist(time_obs[df['group'] == genotype].values, density=True, alpha=0.8)
# ax.fill(X_plot[:, 0], exp_dens, fc='#AAAAFF', alpha=0.8)
# plt.xticks(ticks=np.percentile(X_plot.squeeze(), [(i+0.5)*100/groups_len for i in range(groups_len)]), labels=time_idx, rotation=15)
# plt.ylabel('Density')

# ax = plt.subplot(212)
# xticks = np.percentile(X_plot.squeeze(), [(i+0.5)*100/groups_len for i in range(groups_len)])

# x2 = xticks[:, np.newaxis]
# exp2 = np.exp(kde.score_samples(x2))

# ax.bar(xticks, exp2 * len(X), alpha=0.8, width=1)
# plt.xticks(ticks=xticks, labels=time_idx, rotation=15)
# plt.ylabel('Frequency')
# plt.show()


# ************************ Group Dynamic *******************************
def genotype_obs(genotype):
    X = time_obs[df['group'] == genotype].values[:, np.newaxis]
    kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(X)
    log_dens = kde.score_samples(X_plot)
    return genotype, np.exp(log_dens) * len(X)


epi_year_count_kde = pd.DataFrame(dict([genotype_obs(x) for x in group_order]))

# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 4))

df2[group_order].plot.bar(ax=ax1, legend=False, colormap=cmap, stacked=True, rot=15, xlabel='', ylabel='# of genomes', alpha=0.8)

epi_year_pct_kde = epi_year_count_kde.div(epi_year_count_kde.sum(axis=1), axis=0)
epi_year_pct_kde[group_order].plot.area(colormap=cmap, ax=ax2, ylabel='Fraction of genomes', xlabel='Time', alpha=0.8)
_ = plt.xticks(ticks=np.percentile(np.linspace(0, 999, 1000), [(i+0.5)*100/groups_len for i in range(groups_len)]), labels=time_idx, rotation=15)
_ = plt.ylim(0, 1)
_ = plt.xlim(0, 1000)
plt.grid(False)
plt.legend(loc='upper left', bbox_to_anchor=(1.01, 3.35), fontsize=8)
plt.subplots_adjust(bottom=0.15, right=0.87, hspace=0.35)
plt.show()
