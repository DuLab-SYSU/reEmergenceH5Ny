import numpy as np
from functools import partial
from ete3 import Tree, faces, TreeStyle
import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib.colors import to_hex


def parse_tree(tree_file, format=0, quoted_node_names=True):
    t = Tree(tree_file, format=format, quoted_node_names=quoted_node_names)
    MPD_values = []
    for node in t.traverse():
        if 'mpd' in node.features:
            MPD_values.append(float(node.mpd))
    return t, sorted(MPD_values)


def find_clstr(node, MPD_threshold):
    if 'mpd' not in node.features:
        return False
    if (float(node.mpd) < MPD_threshold) and (int(node.n_leaf) >= 3) and (float(node.support) >= 0.7):
        return True
    else:
        return False


def main(t, MPD_values, pct_th):
    MPD_median = np.quantile(MPD_values, pct_th)
    find_clstr2 = partial(find_clstr, MPD_threshold=MPD_median)
    sub_trees = list(t.iter_leaves(is_leaf_fn=find_clstr2))
    return MPD_median, sub_trees


def tree_view(t, sub_trees, save=False):
    # viridis = cm.get_cmap('viridis', len(sub_trees))
    # viridis_hex = list(map(to_hex, viridis.colors))
    viridis_hex = [
        '#BC3C29', '#0072B5', '#E18727', '#20854E', '#7876B1', '#FFDC91',
        '#EE4C97', '#6F99AD'
    ]

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.scale = 100
    ts.show_scale = False

    for i, node in enumerate(sub_trees):
        group_color = viridis_hex[i]
        node.img_style['fgcolor'] = group_color
        node.img_style['hz_line_color'] = group_color
        node.img_style['vt_line_color'] = group_color

        longNameFace = faces.TextFace("Group %s" % i,
                                      fsize=6,
                                      fgcolor=group_color)
        longNameFace.margin_top = 1
        longNameFace.margin_right = 2
        longNameFace.margin_left = 2
        longNameFace.margin_bottom = 1
        # longNameFace.border.width = 1
        # node.add_face(longNameFace, column=0, position='branch-right')
        for des in node.iter_descendants():
            des.img_style['size'] = 0
            des.img_style['fgcolor'] = group_color
            des.img_style['hz_line_color'] = group_color
            des.img_style['vt_line_color'] = group_color

    if save:
        t.render(save, tree_style=ts)
    else:
        t.show(tree_style=ts)


if __name__ == "__main__":
    SEGMENT = 'HA'
    tree_path = 'data/04_beast_res/beast_res2/segment_%s_attr.newick' % SEGMENT

    t, MPD_values = parse_tree(tree_path)
    x, y, z = [], [], []
    for pct_th in np.linspace(0, 1, 1000, endpoint=False):
        MPD_median, sub_trees = main(t, MPD_values, pct_th)
        x.append(pct_th)
        y.append(len(sub_trees))
        z.append(MPD_median)

    for a, b, c in zip(x[-20:], z[-20:], y[-20:]):
        print(a, b, c)

    # fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
    # ax1.hist(MPD_values, alpha=0.7, bins=30)
    # ax1.set_xlabel(r'$MPD$')
    # ax1.set_ylabel(r'$Number of nodes$')
    # ax2.plot(x, z, alpha=0.7)
    # ax2.set_xlabel(r'$Quantile$')
    # ax2.set_ylabel(r'$MPD$')
    # ax3.scatter(x, y, alpha=0.7, s=2, marker='D')
    # ax3.plot(x, y, alpha=0.7)
    # ax3.set_xlabel(r'$Quantile$')
    # ax3.set_ylabel(r'$Number of clusters$')
    # plt.subplots_adjust(hspace=0.5)
    # plt.show()

    def clstr(pct_th):
        t, MPD_values = parse_tree(tree_path)

        MPD_median, sub_trees = main(t, MPD_values, pct_th)
        print('*' * 40)
        print("Quantile %s: " % SEGMENT, MPD_median)
        print("Number of sub-trees: ", len(sub_trees))

        tree_view(t, sub_trees)


    clstr(0.997)

    # th_dict = {'PB2': 0.998, 'PB1': 0.998, 'PA': 0.996, 'HA': 0.996, 'NP': 0.994, 'NA': 0.995, 'MP': 0.996, 'NS': 0.996}

    # for segment_name, pct_th in th_dict.items():
    #     t2, MPD_values = parse_tree('beast/beast_res2/segment_%s_attr.newick' % segment_name)
    #     MPD_median, sub_trees = main(t2, MPD_values, pct_th)
    #     print('*' * 20)
    #     print(segment_name)
    #     print("Quantile %s: " % pct_th, MPD_median)
    #     print("Number of sub-trees: ", len(sub_trees))
    #     tree_view(t2, sub_trees)

    #     # with open('clstr_%s' % segment_name, 'w') as f:
    #     #     sub_idx = 0
    #     #     for sub_tree in sub_trees:
    #     #         for leaf in sub_tree.get_leaf_names():
    #     #             f.write('%s\t%s\n' % (leaf, sub_idx))
    #     #         sub_idx += 1
