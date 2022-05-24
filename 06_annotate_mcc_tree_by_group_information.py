import os
import glob
import itertools
import pandas as pd
from ete3 import Tree, TreeStyle, NodeStyle, AttrFace, TextFace, RectFace
from matplotlib import cm
from matplotlib.colors import to_hex, ListedColormap
import dendropy
from cycler import cycler
import json


# * read genotype file
genotype = pd.read_csv('./genotype.csv', sep='\t', index_col=9)

acc_after_2020 = genotype.index.values

# * read metadata file
metainfo = pd.read_csv('./dataset1_metainfo.csv', sep='\t', index_col=0)
metainfo = metainfo.join(genotype['group'])

metainfo2 = pd.read_csv('./nextstrain/metadata1.tsv', sep='\t', index_col=2)
node_time = json.load(open('./nextstrain/node.json'))
node_muts = json.load(open('./nextstrain/node1.json'))
node_loc = json.load(open('./nextstrain/node2.json'))

# * define color schema
# * group for segment
color_hex = ['#BC3C29', '#E18727', '#0072B5', '#20854E', '#7876B1', '#6F991D', '#FFDC91', '#EE4C97']
color_map = dict(zip(range(8), color_hex))


# * group for strain
# 0_0_0_0_0_0_0_0    345
# 1_1_1_1_1_1_1_1     77
# 2_2_1_1_2_1_1_1     13
# 3_2_2_1_2_1_1_2     10
# 4_2_3_0_3_2_0_0      3
# 1_1_1_1_2_1_1_1      2
# 0_0_5_0_0_0_0_0      1
# 3_5_1_1_2_1_1_1      1
# 5_4_4_0_4_0_3_3      1
# 6_3_6_0_5_0_2_4      1
# 3_1_2_1_2_1_1_1      1
group_order = ['1_1_1_1_1_1_1_1', '0_0_0_0_0_0_0_0', '2_2_1_1_2_1_1_1',
               '3_2_2_1_2_1_1_2', '4_2_3_0_3_2_0_0', '1_1_1_1_2_1_1_1',
               '0_0_5_0_0_0_0_0', '3_5_1_1_2_1_1_1', '5_4_4_0_4_0_3_3',
               '6_3_6_0_5_0_2_4', '3_1_2_1_2_1_1_1']


color_country = {'Europe': '#BC3C29', 'JapanKorea': '#E18727', 'Africa': '#0072B5', 'China': '#20854E', 'WestAsia': '#7876B1', 'NorthAmerica': '#6F991D', 'SouthAsia': '#FFDC91'}



def layout(node):
    node.img_style['size'] = 0
    node.img_style['vt_line_width'] = 12
    node.img_style['hz_line_width'] = 12
    node.img_style['vt_line_color'] = color_country[node.region]
    node.img_style['hz_line_color'] = color_country[node.region]

    if node.is_leaf():
        isl_acc = node.seq_acc
        if isl_acc in acc_after_2020:
            for idx, seg_name in enumerate(['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']):
                group_i = int(genotype.loc[isl_acc, seg_name])
                # g = RectFace(width=72, height=6, fgcolor=color_map[group_i], bgcolor=color_map[group_i])
                g = RectFace(width=60, height=8, fgcolor=color_map[group_i], bgcolor=color_map[group_i])
                node.add_face(g, idx, 'aligned')


def main():
    # * read mcc tree & trans for ete3
    SEGMENT = 'HA'
    tree_to_visualization = Tree('./nextstrain/timetree.tree', format=1)
    tree_to_visualization.ladderize(1)

    # * annotate tree with metainfo
    for node in tree_to_visualization.traverse("postorder"):
        strain = node.name
        if node.is_leaf():
            seq_acc = metainfo2.loc[strain, 'strain_epi']
        # date = metainfo2.loc[strain, 'Collection_Date']
        # host = metainfo2.loc[strain, 'host2']
        region = node_loc['nodes'][strain]['region']
        date = node_time['nodes'][strain]['date']
        muts = ','.join(node_muts['nodes'][strain]['muts'])

        # group = genotype.loc[seq_acc, 'group']

        # leaf.name = seq_acc
        node.add_features(seq_acc=seq_acc, region=region, date=date, muts=muts)

    ts = TreeStyle()
    # * legend
    for i, (a, b) in enumerate(color_country.items()):
        g = RectFace(width=100, height=100, fgcolor='black', bgcolor=b)
        l = TextFace(str(a), fsize=48)
        l.margin_left = 50
        ts.legend.add_face(g, column=0)
        ts.legend.add_face(l, column=1)

    g = RectFace(width=100, height=100, fgcolor='white', bgcolor='white')
    ts.legend.add_face(g, column=2)

    for i, (a, b) in enumerate(color_map.items()):
        g = RectFace(width=100, height=100, fgcolor='black', bgcolor=b)
        l = TextFace(str(a), fsize=48)
        l.margin_left = 50
        ts.legend.add_face(g, column=3)
        ts.legend.add_face(l, column=4)
    ts.legend_position = 1
    return tree_to_visualization, ts



if __name__ == '__main__':
    tree, ts = main()

    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.scale = 60000
    ts.show_scale = False
    ts.margin_left = 50
    ts.margin_right = 50
    ts.margin_top = 20
    ts.margin_bottom = 20

    tree.show(tree_style=ts, layout=layout)
    # tree.render('fig1a.pdf', tree_style=ts, layout=layout, units='mm', h=45)

# tree.show(tree_style=ts, layout=layout)
# # tree.render('fig1a.pdf', tree_style=ts, layout=layout, units='mm', h=200)

# # ts.title.add_face(TextFace(tree_file, fsize=20), column=0)
# # tree.render(tree_file.replace('.tree', '.pdf'), tree_style=ts, layout=layout)
