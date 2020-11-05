import argparse
import sys
import os
from ete3 import PhyloTree, NCBITaxa, SeqGroup, Tree
import re
from collections import defaultdict
from multiprocessing import Pool
import itertools


def get_species(name):
    return name.split('.')[0]

def iter_tree_paths():
    ''' Iterates over all tree files, yielding one tree path at a time '''
#    print '##'+'\t'.join('Seed (co-)orthologs type target_taxid target_species orthologs treefile'.split())
    for c, treepath in enumerate(open(tree_list)):
        treepath = str(treepath)
        treepath = treepath.strip()
        yield treepath

def process_tree(treepath):
    ''' processes a tree to extract orthology relationships between target taxid and the rest
     of species, organized by orthology type and species code '''
    treepath = str(treepath)
    treepath = treepath.rstrip()
    t = PhyloTree(treepath, sp_naming_function=get_species)
    treefile = os.path.basename(treepath)
    t.dist = 0

    outgroup = t.get_midpoint_outgroup()
    try:
        t.set_outgroup(outgroup)
        t.standardize()
    except:
        if args.pairs_table:
            if len(t) == 1:
                sys.stderr.write(treefile+'len(t) == 1'+'\n')
                return ([],[])   
                #return (['aa', 'aa'] ,[['aa', 'aa']])
            
            else:
                sys.stderr.write(treefile+'len(t) != 1'+'\n')
                l = t.get_leaf_names()
                r = l[0]
                t.set_outgroup(r)
                pass
                #return ([],[])
                #return  (['None', 'None'] ,[['None', 'None']])
        else:
            if len(t) == 1:
                sys.stderr.write(treefile+'len(t) == 1'+'\n')
                return []
            else:
                sys.stderr.write(treefile+'len(t) != 1'+'\n')
                return []

    
    names={}
    for leaf in t:
        try:
            sp = str(leaf.name.split('.')[0])
            leaf.taxid = str(sp)
            sci_name = ncbi.get_taxid_translator([sp])
            names[sp] = sci_name[int(sp)]

        except:
            names[sp] = ''
       
        if args.conv_table:
            try:
                good_name = "%s" %(conversion[leaf.name][0])
            except:
                good_name = leaf.name
            leaf.good_name = good_name

    node2content = t.get_cached_content()
    target_species = set([target_taxid])
    
    def is_sp_specific(_node):
        _species = set([_leaf.species for _leaf in  node2content[_node]])
        if not (_species - target_species):
            return True
        return False

    #traverse only target taxid leaves
    if collapse == 'yes':
        for n in t.get_leaves(is_leaf_fn=is_sp_specific):
            if n.children:
                for ch in n.get_children():
                    ch.detach()
                n.taxid = target_taxid
                n.name = "{%s}" %('|'.join([_lf.name for _lf in node2content[n]]))
                if args.conv_table:
                    n.good_name = "{%s}" %('|'.join([_lf.good_name for _lf in node2content[n]]))

    all_ortholgs_tree = []
    all_ortholgs_pairs = []     
    event_lines = []

    for ev in t.get_descendant_evol_events():
        if ev.etype == "S":
            source_seqs = ev.node.children[0]
            ortho_seqs = ev.node.children[1]

            if target_taxid:
                sp_1=set()
                for leaf in source_seqs:
                    sp_1.add(leaf.taxid)
                sp_2=set()
                for leaf in ortho_seqs:
                    sp_2.add(leaf.taxid)
                
                if str(target_taxid) in sp_1:
                    source_seqs, ortho_seqs = source_seqs, ortho_seqs
                elif str(target_taxid) in sp_2:
                    source_seqs, ortho_seqs = ortho_seqs, source_seqs
                else:
                    continue           

            if args.conv_table:
                co_orthologs = [leaf.good_name for leaf in source_seqs]
                co_orthologs.sort()
            else:
                co_orthologs = [leaf.name for leaf in source_seqs]
                co_orthologs.sort()

            orthologs = defaultdict(set)
            for leaf in ortho_seqs:
                sp = str(leaf.name.split('.')[0])
                if args.conv_table:
                    orthologs[sp].add(leaf.good_name)
                else:
                    orthologs[sp].add(leaf.name)  

            if len(source_seqs) == 1:
                _otype = "one-to-"
            else:
                _otype = "many-to-"
        
            for sp, orth in orthologs.items():
                if len(orth) == 1:
                    otype = _otype + "one"
                else:
                    otype = _otype + "many"
            
                event_lines.append('\t'.join([','.join(co_orthologs), otype , str(sp), ','.join(sorted(orth)), treefile, names[sp],'\n']))
                

            if args.pairs_table:
                
                source_seqs_names = []
                ortho_seqs_names= []

                for node in source_seqs:
                    for leaf in node:
                        if args.conv_table:
                            name = leaf.good_name
                        else: 
                            name = leaf.name
                        source_seqs_names.append(name)
            
                for node in ortho_seqs:
                    for leaf in node:
                        if args.conv_table:
                            name = leaf.good_name
                        else:
                            name = leaf.name
                        ortho_seqs_names.append(name)

                all_ortholgs_node = itertools.product(source_seqs_names, ortho_seqs_names)
                all_ortholgs_tree.append(all_ortholgs_node)

                for node in all_ortholgs_tree:
                    for pair in node:
                        all_ortholgs_pairs.append(pair)

                #return (event_lines, all_ortholgs_pairs)

    if args.pairs_table:
        return (event_lines, all_ortholgs_pairs)
    else:
        return (event_lines)
                
                


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--collapse', dest='collapse', default='yes', type=str, choices=["yes", "no"])
    parser.add_argument('--target', dest='target_taxid', type=str, required=False)
    parser.add_argument('--tree_list', dest='tree_list', type=str)
    parser.add_argument('--output', dest='out_table', default='orthologs.tsv', type=str)
    parser.add_argument('--output_pair', dest='pairs_table', type=str, required = False)
    parser.add_argument('--conversion', dest='conv_table', required = False)
    parser.add_argument('--cpu', dest='cpu', type=int, default=5)
    args = parser.parse_args()

    target_taxid = args.target_taxid
    collapse = args.collapse
    tree_list = args.tree_list
    CPU = args.cpu
    out_file = args.out_table

    conversion = {}
    if args.conv_table:
        for line in open(args.conv_table):
            code, realname, annotations = line.split("\t")
            conversion[code] = [realname, annotations]

    ncbi = NCBITaxa()
    
    pool = Pool(CPU)
    if args.pairs_table:
        with open(out_file, 'w') as OUT_t, open(args.pairs_table, 'w') as OUT_pair:
            OUT_t.write('##'+'\t'.join('Seed_(co-)orthologs type target_taxid  orthologs treefile target_species'.split())+'\n')
            for event_lines, pair_line in pool.imap_unordered(process_tree, iter_tree_paths()):
                for line in event_lines:
                    OUT_t.write(line)
                for line in pair_line:
                    OUT_pair.write('\t'.join(line)+'\n')

                    
    else:
        with open(out_file, 'w') as OUT:
            OUT.write('##'+'\t'.join('Seed_(co-)orthologs type target_taxid  orthologs treefile target_species'.split())+'\n')
            for event_lines in pool.imap_unordered(process_tree, iter_tree_paths()):
                try:
                    for line in event_lines:
                        OUT.write(line)
                except:
                    OUT.write('ERROR'+'\n')  