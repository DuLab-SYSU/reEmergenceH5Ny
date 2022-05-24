from ete3 import Tree
import numpy as np
import itertools
import matplotlib.pyplot as plt
import sys


def cal_dist_mat(tree):
    leaf_names = tree.get_leaf_names()
    num_leaves = len(leaf_names)
    dist_mat = np.zeros((num_leaves, num_leaves))
    for i, leaf_i in enumerate(leaf_names):
        for j, leaf_j in enumerate(leaf_names):
            if j > i:
                dist = tree.get_distance(leaf_i, leaf_j)
                dist_mat[i, j] = dist
                dist_mat[j, i] = dist
            if j == i:
                dist_mat[i, j] = 0
    return dist_mat, leaf_names


def cal_MDP(node, leaf_names, dist_mat):
    sub_leaf_names = node.get_leaf_names()
    leaf_combinations = list(itertools.combinations(sub_leaf_names, 2))
    n_pair = len(leaf_combinations)
    dist_l = [dist_mat[leaf_names.index(leaf_i), leaf_names.index(leaf_j)] for leaf_i, leaf_j in leaf_combinations]
    cum_dist = sum(dist_l)
    return cum_dist / n_pair


if __name__ == '__main__':
    SEGMENT, tree_path = sys.argv[1], sys.argv[2]

    t = Tree(tree_path, format=0, quoted_node_names=True)
    # t.set_outgroup(t.get_midpoint_outgroup())
    t.ladderize(1)

    dist_mat, leaf_names = cal_dist_mat(t)

    for node in t.traverse():
        if node.is_leaf():
            node.add_features(mpd=float(node.dist), n_leaf=1)
        else:
            mpd = cal_MDP(node, leaf_names, dist_mat)
            n_leaf = len(node.get_leaves())
            node.add_features(mpd=mpd, n_leaf=n_leaf)


    t.write(format=1, outfile='segment_%s_attr.newick' % SEGMENT, features=['mpd', 'n_leaf', 'support'])
