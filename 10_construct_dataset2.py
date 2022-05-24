import os
import glob
from re import sub
from ete3 import Tree, TreeStyle, AttrFace, TextFace
from numpy.matrixlib import defmatrix
import pandas as pd
from collections import defaultdict
import subprocess
from Bio import SeqIO
from io import StringIO


# read dataset1 metainfo
dataset1 = pd.read_csv('dataset1_metainfo.csv', sep='\t', index_col=0)
dataset1_isl_accs = dataset1.index.tolist()
# read aiv database metainfo
all_aiv_info = pd.read_csv('data/all_aiv_info.csv', sep='\t', index_col=0)

genotype = pd.read_csv('genotype2.csv', sep='\t', index_col=9)


with open('accs_sample_after_2020.txt') as f:
    accs_after_2020 = [x.strip() for x in f.readlines()]
print(len(accs_after_2020), end='\t')
print(' '.join(accs_after_2020[0:3]) + ' ... ' + ' '.join(accs_after_2020[-3:]))

# read tree
SEGMENT = 'PA'
tree_file = 'data/07_blast_res/blast_res_%s.tree' % SEGMENT
tree = Tree(tree_file, format=0, quoted_node_names=True)

mid_node = tree.get_midpoint_outgroup()
tree.set_outgroup(mid_node)
tree.ladderize(1)

# read group infomation
acc2clstr = genotype[SEGMENT].astype(int).to_dict()
clstr2accs = defaultdict(list)
for k, v in acc2clstr.items():
    clstr2accs[v].append(k)

color_hex2 = ['#BC3C29', '#0072B5', '#E18727', '#20854E', '#7876B1', '#FFDC91', '#EE4C97', '#6F99AD']
color_map2 = {a: color_hex2[a] for a in clstr2accs.keys()}

# annotate leaf nodes
for leaf in tree.iter_leaves():
    seq_acc, attr = leaf.name.split(' ', 1)
    isl_name, isl_acc, date, seg_name = attr.split('#')
    leaf.name = isl_acc
    subtype = all_aiv_info.loc[isl_acc, 'Subtype']
    group = acc2clstr.get(isl_acc)
    if group != None:
        leaf.add_feature('group', group)
    else:
        leaf.add_feature('group', 'none')
    leaf.add_features(seq_acc=seq_acc, isl_name=isl_name, date=date, seg_name=seg_name, subtype=subtype)
    if isl_acc in dataset1_isl_accs:
        leaf.add_feature('is_ds1', True)
    else:
        leaf.add_feature('is_ds1', False)


#* extract ext seqs by MCRA
ext_accs = []
ext_accs_dict = defaultdict(list)
for group in clstr2accs.keys():
    node_l = [node_i for node_i in tree.get_monophyletic(values=[group], target_attr='group')]
    mcra = tree.get_common_ancestor(*node_l)
    mcra.img_style['bgcolor'] = color_map2[group]



#* extract ext seqs by monophyletic
# ext_accs = []
# ext_accs_dict = defaultdict(list)
# for group in clstr2accs.keys():
#     for node_i in tree.get_monophyletic(values=[group], target_attr='group'):
#         node_i.img_style['bgcolor'] = color_map2[group]
#         pp_node = node_i.up.up
#         pp_node.img_style['bgcolor'] = color_map2[group]
#         for leaf in pp_node.iter_leaf_names():
#             if leaf in dataset1_isl_accs:
#                 continue
#             elif leaf in ext_accs:
#                 continue
#             else:
#                 ext_accs_dict[group].append(leaf)
#                 ext_accs.append(leaf)
# # ext_accs
# print(len(ext_accs), end='\t')
# print(' '.join(ext_accs[0:3]) + ' ... ' + ' '.join(ext_accs[-3:]))
# for k, v in ext_accs_dict.items():
#     print(k, len(v))

#! extract ext metainfo
# meta_info_keep = all_aiv_info[all_aiv_info.index.isin(accs_after_2020 + ext_accs)]
# meta_info_keep.to_csv('ext_metainfo_%s.csv' % SEGMENT, sep='\t', index=True)

#! extract ext sequences
# HA_accs = meta_info_keep['%s Segment_Id' % SEGMENT].to_list()
# blastdbcmd = subprocess.run(
#     'blastdbcmd -db fludb/fludb2 -entry_batch -',
#     shell=True,
#     universal_newlines=True,
#     input='\n'.join(HA_accs),
#     stdout=subprocess.PIPE,
#     stderr=subprocess.DEVNULL
# )

# stdout = blastdbcmd.stdout
# records = list(SeqIO.parse(StringIO(stdout), 'fasta'))
# for record in records:
#     isl_acc = record.description.split('#')[1]
#     record.id = isl_acc
#     record.name = ''
#     record.description = ''

# SeqIO.write(records, 'ext_%s.fasta' % SEGMENT, 'fasta')


def layout2(leaf):
    if leaf.is_leaf():
        if leaf.name in dataset1_isl_accs:
            leaf.add_feature('is_ds1', True)
            leaf.img_style['vt_line_color'] = 'red'
            leaf.img_style['hz_line_color'] = 'red'
            leaf.img_style['fgcolor'] = 'red'
            leaf.img_style['size'] = 6
            leaf.img_style['vt_line_type'] = 0
            leaf.img_style['hz_line_type'] = 0
        elif leaf.subtype == 'A / H5N8':
            leaf.add_feature('is_ds1', False)
            leaf.img_style['vt_line_color'] = 'blue'
            leaf.img_style['hz_line_color'] = 'blue'
            leaf.img_style['fgcolor'] = 'blue'
            leaf.img_style['size'] = 6
            leaf.img_style['vt_line_type'] = 2
            leaf.img_style['hz_line_type'] = 2
            # long_str = AttrFace('name', fsize=15, fgcolor='blue')
            # leaf.add_face(long_str, 2, 'aligned')

        elif leaf.subtype != 'A / H5N8':
            leaf.img_style['fgcolor'] = 'green'
            leaf.img_style['size'] = 18
            leaf.img_style['vt_line_color'] = 'green'
            leaf.img_style['hz_line_color'] = 'green'
            leaf.img_style['fgcolor'] = 'green'
            leaf.img_style['vt_line_type'] = 1
            leaf.img_style['hz_line_type'] = 1
            # long_str = AttrFace('name', fsize=15, fgcolor='green')
            # leaf.add_face(long_str, 2, 'aligned')


ts = TreeStyle()
ts.show_leaf_name = False
# ts.show_branch_length = True
ts.show_branch_support = False
ts.scale = 40000
ts.show_scale = False
ts.draw_guiding_lines = True

ts.title.add_face(TextFace(tree_file, fsize=60), column=0)
tree.show(tree_style=ts, layout=layout2)
# tree.show(tree_style=ts)

# tree.render(tree_file.replace('.tree', 'test.pdf'), tree_style=ts, layout=layout2)
