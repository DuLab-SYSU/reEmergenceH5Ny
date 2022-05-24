import os
import re
import pandas as pd
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, RectFace, TextFace


with open('./ancestor_H9/annotated_tree.nexus') as f:
	raw_tree = f.read()


group_info = pd.read_csv('./ancestor_H9/clade_info2.csv', sep='\t')
taxa_l = group_info.query('Clade_info == "#003366"')['Name'].to_list()


matched = re.search(r'Tree tree1=(.+)\sEnd;', raw_tree)
string_ = matched.group(1)


string2_ = re.sub(r'\[&.*?\]', '',string_)

node_name_with_anno = re.findall(r'(\w+):.+?\[&(.*?)\]', string_)
node_name_with_anno = list(map(lambda x: (x[0], x[1].split('=')[1]), node_name_with_anno))
node_name_with_anno = dict(node_name_with_anno)

tree = Tree(string2_, format=1)

trunk_nodes = []
for group, df in group_info.groupby('Clade_info'):
	taxa_i = df['Name'].to_list()
	# print(group, taxa_i[:3])
	anc = tree.get_common_ancestor(taxa_i)
	if len(taxa_i) == len(anc):
		lineages = anc.get_ancestors()
		trunk_nodes.extend(lineages)
		trunk_nodes.append(anc)
		trunk_nodes.extend(anc.get_children())
		for node in anc.traverse():
			node.img_style['hz_line_color'] = group
			node.img_style['vt_line_color'] = group


trunk_nodes = list(set(trunk_nodes))

for node in tree.traverse():
	if node.name in node_name_with_anno:
		node.add_feature('mutations', eval(node_name_with_anno[node.name]))


for node in trunk_nodes:
	node.img_style['hz_line_color'] = 'red'
	if "mutations" in node.features:
		node.add_face(AttrFace('mutations', fsize=200, fgcolor='red'), 0, 'branch-top')


tree.ladderize(1)

ts = TreeStyle()
ts.show_leaf_name = False
ts.show_branch_length = False
ts.show_branch_support = False
ts.scale = 60000
ts.mode = 'r'
ts.show_scale = True
ts.margin_right = 2000
ts.extra_branch_line_type = 1
ts.extra_branch_line_color = 'gray'
ts.complete_branch_lines_when_necessary = True

tree.show(tree_style=ts)
