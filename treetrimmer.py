#!/usr/bin/env python
# PolyFastA.py
# Santiago Sanchez-Ramirez, University of Toronto, santiago.snchez@gmail.com

import collections
import argparse
import sys
import os
from itertools import product
import multiprocessing as mp
try:
    import dendropy as dp
    from dendropy.calculate import treecompare
except:
    print("Dendropy is not installed. Try: pip install [--user] dendropy")
    sys.exit(1)

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def parse_arguments(treetrimmer):
    # parse arguments
    
    parser = argparse.ArgumentParser(prog="treetrimmer.py",
    formatter_class=argparse.RawTextHelpFormatter,
    description="""
    """ + treetrimmer + """ reduces conflict in phylogenomic data sets by resolving
    paralogy using a phylogenetic approach. Tree subsets within gene clusters
    or orthogroups are compared to a species tree, and the subset with less
    conflict is selected.

    Weighted Robinson-Foulds distances are used as metrics to assess all
    possible resolved tree subsets.

    The program depends on the DendroPy phylogenetic library:
    https://dendropy.org

    \n""",
    epilog="""
    Example:

    python treetrimmer.py -p 10 -a myAlignments/ -af fasta -t myTrees/ -tf newick -st my_species_tree.tre -o output_dir/ 

    -a, -t, and -o refer to directory names, and -st refers to the species tree file.

    If the species tree labels are these:
    (((A,B),C),D);
    
    The expected gene tree labels should be:
    ((((A.g1,A.g5),B.g7),C.g45),(C.g32,D.g5));

    Where the \".\" dot delimits the species label.
    Note that this can be changed with the --tip_delim argument.

    Gene trees are expected to have branch lenghts.
    Branches/Nodes with low support can be collapsed into polytomies with the --support_limit argument.

    \n""")
    parser.add_argument(
    '--file_list', '-fl', metavar="PATH", type=str, required=True,
    help='a directory with alignment files.')
    parser.add_argument(
    '--aln_dir', '-ad', metavar="PATH", type=str, required=True,
    help='a directory with alignment files.')
    parser.add_argument(
    '--tree_dir', '-td', metavar="PATH", type=str, required=True,
    help='a directory with tree files.')
    parser.add_argument(
    '--species_tree', '-st', metavar="PATH", type=str, required=True,
    help='a species tree.')
    parser.add_argument(
    '--aln_format', '-af', type=str, default='fasta', choices=['fasta','nexus','phylip'],
    help='format of alignment files. (default: fasta)')
    parser.add_argument(
    '--tree_format', '-tf', type=str, default='newick', choices=['newick','nexus'],
    help='format of tree files. (default: newick)')
    parser.add_argument(
    '--species_tree_format', '-stf', type=str, default='newick', choices=['newick','nexus'],
    help='format of species tree file. (default: newick)')
    parser.add_argument(
    '--tip_delim', '-d', metavar="CHAR", type=str, default=".",
    help='character to delimit the species label in gene tree. (default: \'.\'')
    parser.add_argument(
    '--tip_delim_index', '-di', metavar="INT", type=int, default=0,
    help='zero-based index of the species label.')
    parser.add_argument(
    '--support_limit', '-sl', metavar="[INT | FLOAT]", type=str, default='0',
    help='zero-based index of the species label.')
    parser.add_argument(
    '--threads', '-p', metavar="INT", type=int,
    help='numer of processes to run in parallel. (default: all available)')
    parser.add_argument(
    '--silent', '-si', action="store_true", default=False,
    help='run in silent mode. (defalt: False)')
    parser.add_argument(
    '--out', '-o', metavar="DIR", type=str, default="treetrimmer_output",
    help='name of output directory. (default: treetrimmer_output)')

    args = parser.parse_args()
    return(args)

def split_multi_and_singlecopy(tree, delimiter=".", index=0):
    tn = tree.taxon_namespace
    labs = [ x.label for x in tn ]
    labs_spp = [ x.split(delimiter)[index] for x in labs ]
    spp_counts = dict(collections.Counter(labs_spp))
    spp_para = [ x for x,y in zip(spp_counts.keys(), spp_counts.values()) if y > 1 ]
    labs_para = [ [ labs[i] for i in range(len(labs)) if labs_spp[i] == p ] for p in spp_para ]
    labs_sing = [ labs[labs_spp.index(x)] for x,y in zip(spp_counts.keys(), spp_counts.values()) if y == 1 ] 
    return(labs_para, labs_sing)

def get_tip_edge(tree, lab):
    edge = tree.edges(filter_fn=lambda x: x.is_leaf() and x.head_node.taxon.label.startswith(lab))[0]
    return(edge)

def get_sequence(aln, lab, no_gaps=True):
    taxon_i = [ i for i in range(len(aln.taxon_namespace)) if aln.taxon_namespace[i].label.startswith(lab) ][0]
    if no_gaps:
        return("".join([ x for x in aln[taxon_i].symbols_as_list() if x != "-" ]))
    else:
        return(aln[taxon_i].symbols_as_string())

def resolve_inparalogs(tree, dat, inpara):
    keep_inpara = []
    seq_len = [ len(get_sequence(dat, i)) for i in [ i.label for i in dat.taxon_namespace ]]
    mean_seq_len = sum(seq_len)/len(seq_len)
    for group in inpara:
        seq_len = [ len(get_sequence(dat, i)) for i in group ]
        seq_len_diff = [ abs(mean_seq_len - x) for x in seq_len ]
        if not all([seq_len_diff[0] == i for i in seq_len_diff]):
            keep_inpara.append(group[ seq_len_diff.index(min(seq_len_diff)) ])
        else:
            tip_length = [ get_tip_edge(tree, i).length for i in group ]
            keep_inpara.append(group[ tip_length.index(min(tip_length)) ])
    return(keep_inpara)

def is_monophyletic(tree, taxa):
    mrca = tree.mrca(taxon_labels=taxa)
    leaf_names = [ l.taxon.label for l in mrca.leaf_nodes() ]
    return(collections.Counter(leaf_names) == collections.Counter(taxa))

def consolidate_taxon_namespace(tree1, tree2, delimiter=".", index=0):
    tns = dp.TaxonNamespace()
    t1 = dp.Tree.get(data=tree1.as_string(schema="newick"), schema="newick", taxon_namespace=tns)
    t2 = tree2.clone(2)
    for t in t2.taxon_namespace:
        t.label = t.label.split(delimiter)[index]
    t3 = dp.Tree.get(data=t2.as_string(schema="newick"), schema="newick", taxon_namespace=tns)
    return(t1,t3)

def collapse_low_branches(tree, thresh):
    tree_copy = tree.clone(2)
    for n in tree_copy.nodes(filter_fn=lambda x: not x.is_leaf() ):
        if n.label and int(n.label) < thresh:
            n.edge_length = None
        elif not n.label:
            n.edge_length = None
    tree_copy.collapse_unweighted_edges()
    return(tree_copy)

def reduce_alignment(dat, taxa):
    d2 = dat.clone(2)
    tx = []
    for taxon in d2.taxon_namespace:
        if taxon.label not in taxa:
            tx.append(taxon)
    d2.remove_sequences(tx)
    return(d2)

def select_combination_RF(tree1, tree2, dat, all_comb, singletons, inpara, thresh):
    RF = []
    trees_origlabs = []
    trees_spplabs = []
    taxon_name_combinations = []
    for cb in all_comb:
        tmp_lab = singletons + inpara + list(cb)
        t2 = collapse_low_branches(tree2, thresh).extract_tree_with_taxa_labels(labels=tmp_lab)
        tree_labels = [ x.label for x in t2.taxon_namespace ]
        t1, t3 = consolidate_taxon_namespace(tree1, t2)
        RF.append(treecompare.weighted_robinson_foulds_distance(t1, t3))
        trees_origlabs.append(t2)
        trees_spplabs.append(t3)
        taxon_name_combinations.append(tmp_lab)
    best_wRF = RF.index(min(RF))
    RF_best_score = RF[best_wRF]
    RF_best_tree = trees_origlabs[best_wRF]
    best_labels = taxon_name_combinations[best_wRF]
    best_dat = reduce_alignment(dat, best_labels)
    return(RF_best_score, RF_best_tree, best_dat)

def worker(args, spp_tree, data):
    file1 = args.aln_dir + "/" + data[0]
    file2 = args.tree_dir + "/" + data[1]
    matrix = dp.ProteinCharacterMatrix.get(path=file1, schema=args.aln_format)
    gene_tree = dp.Tree.get(path=file2, preserve_underscores=True, schema=args.tree_format)
    gene_tree.deroot()
    all_paralogs, singletons = split_multi_and_singlecopy(gene_tree)
    inparalogs = [ x for x in all_paralogs if is_monophyletic(gene_tree, x) ]
    if inparalogs:
        keep_inpara = resolve_inparalogs(gene_tree, matrix, inparalogs)
        paralogs = [ x for x in all_paralogs if x not in inparalogs ]
    else:
        keep_inpara = []
        paralogs = all_paralogs
    all_comb = list(product(*paralogs))
    RF_score, best_tree, r_matrix = select_combination_RF(spp_tree, gene_tree, matrix, all_comb, singletons, keep_inpara, int(args.support_limit))
    r_matrix.write(path=args.out+"/trees"+"/"+data[0], schema="fasta")
    best_tree.write(path=args.out+"/alignments"+"/"+data[1], schema="newick")
    if not args.silent:
        print("Done for",color.BOLD+data[0]+color.END,"and",color.BOLD+data[1]+color.END,", best RF score:",color.BOLD+str(round(RF_score,2))+color.END)

if __name__ == '__main__':
    treetrimmer = color.BOLD + color.GREEN + "TreeTrimmer" + color.END
    args = parse_arguments(treetrimmer)
    if not args.silent:
        print(treetrimmer,"Called arguments:",args,sep="\n")
    spp_tree = dp.Tree.get(path=args.species_tree, preserve_underscores=True, schema=args.species_tree_format)
    spp_tree.deroot()
    if not os.path.exists(args.out):
        os.makedirs(args.out+"/trees")
        os.makedirs(args.out+"/alignments")
    with open(args.file_list,'r') as f:
        files = f.read().splitlines()
        files = [ x.split("\t") for x in files ]
    all_args = ( (args, spp_tree, i) for i in files )
    if args.threads:
        with mp.Pool(args.threads) as pool:
             pool.starmap(worker, all_args)
    else:
        with mp.Pool() as pool:
            pool.starmap(worker, all_args)

