import dendropy
from dendropy.calculate import popgenstat
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from ete3 import Tree, TreeStyle
from itertools import cycle

# * 3 indenpend long mcmc run was merged to obtain a high support tree topology
# * 200,000,000 chain lenght with 40000 log frequency, so 5000 trees for each run
# * 5000 trees was burnin

# * read mcc tree & multiple sequences alignment

def read_data(tree_file, fasta_file):
    dataset1 = dendropy.DataSet()
    tns = dataset1.new_taxon_namespace(label='global')
    dataset1.attach_taxon_namespace(tns)
    dataset1.read(path=tree_file, schema='nexus', preserve_underscores=True)
    dataset1.read(path=fasta_file, schema='fasta', data_type='dna')

    tree = dataset1.tree_lists[0][0]
    # * reorder tree topology
    tree.ladderize(0)
    # * store phylogenetic distance matrix
    pdm = tree.phylogenetic_distance_matrix()
    return tree, pdm


# * calculate mean / median patristic distance for internal nodes
def cal_mpd(node, pdm, method='median'):
    pdc_l = []
    for i, t1 in enumerate(node.leaf_nodes()[:-1]):
        for t2 in node.leaf_nodes()[i+1:]:
            pdc_l.append(pdm.patristic_distance(t1.taxon, t2.taxon))
    if method == 'median':
        return np.median(pdc_l)
    else:
        return np.mean(pdc_l)


def cal_distribution(tree, pdm):
    # * calculate whole tree patristic distance distribution
    pdc_l = []
    for i, t1 in enumerate(tree.leaf_nodes()[:-1]):
        for t2 in tree.leaf_nodes()[i+1:]:
            weighted_patristic_distance = pdm.patristic_distance(t1.taxon, t2.taxon)
            pdc_l.append(weighted_patristic_distance)

    mpd_l = []
    for node in tree.preorder_internal_node_iter():
        mpd_i = cal_mpd(node, pdm)
        node.annotations.add_new(name='mpd', value=mpd_i)
        mpd_l.append(mpd_i)
        node.label = '%.8f' % float(node.annotations.get_value('posterior'))
    return pdc_l, mpd_l


def extract_subtrees(t_percentile, distribution, tree, pdm, log1=False, log2=True):
    # * two strategy here
    # * first,
    threshold = np.percentile(distribution, t_percentile)

    sub_trees = []
    para_group = []

    root = tree.seed_node
    stack = []
    stack.append(root)
    while stack:
        tmp = stack.pop()
        if tmp.is_leaf():
            para_group.append(tmp)
            continue

        posterior_value = float(tmp.annotations.get_value('posterior'))
        num_leaves = len(tmp.leaf_nodes())
        mpd = tmp.annotations.get_value('mpd')

        if mpd < threshold and posterior_value > 0.7 and num_leaves >= 2:
            if log1:
                print("Branch support: %.4f, Median patristic distance: %.6f, # of sequences %d" % (posterior_value, mpd, num_leaves))
            sub_trees.append(tmp)

        else:
            stack.extend(tmp.child_nodes())

    total = 0
    for sub_tree in sub_trees:
        total += len(sub_tree.leaf_nodes())

    if log2:
        print("Percentile threshold: %.2f, Distance threshold: %.6f, # of clusters: %d, # of sequences classfied: %d, # of sequnences unclassfied: %d" % (t_percentile, threshold, len(sub_trees), total, len(para_group)))
    return sub_trees, para_group


if __name__ == '__main__':
    
    fig = plt.figure()
    for idx, SEGMENT in enumerate(['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS']):

        tree_file = './data/03_beast_analysis/beast_db1/combined_%s.mcc.tree' % SEGMENT
        fasta_file = './data/01_sequence_of_dataset1/02_align_unified/db1_seq_by_acc_%s.fasta' % SEGMENT
        tree, pdm = read_data(tree_file, fasta_file)
        pdc_l, mpd_l = cal_distribution(tree, pdm)

        # * sensitive analysis
        x_l = []
        res = []
        res2 = []
        for i in np.arange(90, 100, 0.05):
            sub_trees, para_group = extract_subtrees(i, mpd_l, tree, pdm, log2=True)
            m = len(sub_trees)
            l = len(para_group)
            c_l = [len(sub_tree.leaf_nodes()) for sub_tree in sub_trees]
            c_l.extend([1]*l)
            n = sum(c_l)
            c = (1/n)*(sum([c_i*c_i for c_i in c_l]))
            mcs = c / n
            x_l.append(i)
            res.append(mcs)
            res2.append(m+l)

        ax = fig.add_subplot(2, 4, idx+1)
        ax.plot(x_l, res, 'k-', marker='d', label='Mean cluster size')
        ax.set_title(SEGMENT)

        ax2 = ax.twinx()
        ax2.plot(x_l, res2, 'r-', marker='D', label='# of clusters')

    plt.show()


    # * plot the distribution of whole tree distance and MPD
    # x1 = np.array(pdc_l)
    # x2 = np.array(mpd_l)
    # kde1 = stats.gaussian_kde(x1)
    # kde2 = stats.gaussian_kde(x2)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # # ax.plot(x1, np.zeros(x1.shape), 'b+', ms=10)
    # x_eval = np.linspace(0, 10, num=100)
    # ax.plot(x_eval, kde1(x_eval), 'k-', label="PDC dis")
    # ax.plot(x_eval, kde2(x_eval), 'r-', label="MPD dis")
    # plt.show()

    # sub_trees, para_group = extract_subtrees(99.2, mpd_l, tree, pdm, log=True)
    # tree.write(path='test.nwk', schema='newick', suppress_rooting=True, suppress_internal_node_labels=False)



    # def tree_visualization(tree, pdm, mpd_l, t, SEGMENT):
    #     sub_trees, para_group = extract_subtrees(t, mpd_l, tree, pdm, log1=True)
    #     nwk = tree.as_string(schema='newick', suppress_rooting=True, suppress_internal_node_labels=False)

    #     tree_to_visualization = Tree(nwk, format=0, quoted_node_names=True)

    #     color_list = ['#E64B35', '#4DBBD5', '#00A087', '#F39B7F', '#8491B4', '#7E6148', '#DC0000']
    #     color_generator = cycle(color_list)

    #     f = open('group_%s.txt' % SEGMENT, 'w')
    #     for i, sub_tree in enumerate(sub_trees):
    #         sub = [t.taxon.label for t in sub_tree.leaf_nodes()]
    #         node = tree_to_visualization.get_common_ancestor(sub)
    #         node.img_style['bgcolor'] = next(color_generator)
    #         for sub_i in sub:
    #             f.write("%s\t%s\n" % (i, sub_i))
    #     for para_i in para_group:
    #         i += 1
    #         f.write("%s\t%s\n" % (i, para_i.taxon.label))
    #     f.close()

    #     ts = TreeStyle()
    #     ts.scale = 1000

    #     tree_to_visualization.show(tree_style=ts)
        # tree_to_visualization.render('group_%s.pdf' % SEGMENT, tree_style=ts)


    # t_l = [99, 99, 99, 99.8, 98.6, 99.6, 99, 99.5]
    # seg_l = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'NS', 'MP']
    # for idx, (SEGMENT, t) in enumerate(zip(seg_l, t_l)):
    #     print("********* group results for %s ***********" % SEGMENT)
    #     tree_file = 'data/04_beast_res/beast_db1_longRun/combined_%s.mcc.tree' % SEGMENT
    #     fasta_file = 'data/02_align_unified/db1_seq_by_acc_%s.fasta' % SEGMENT
    #     tree, pdm = read_data(tree_file, fasta_file)
    #     pdc_l, mpd_l = cal_distribution(tree, pdm)

    #     tree_visualization(tree, pdm, mpd_l, t, SEGMENT)

    #     print('\n')
    # tree_visualization(tree, pdm)
    # tree_to_visualization.render('test.png', tree_style=ts)
