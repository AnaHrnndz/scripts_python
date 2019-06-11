import argparse
import sys
import os
from ete3 import PhyloTree, NCBITaxa, SeqGroup, Tree
import re
from collections import defaultdict
from multiprocessing import Pool

def get_species(name):
    return name.split('.')[0]

def iter_tree_paths():
    ''' Iterates over all tree files, yielding one tree path at a time '''
    print '##'+'\t'.join('Seed (co-)orthologs type target_taxid target_species orthologs '.split())
    for c, treepath in enumerate(open("/data/lamprey/lamprey/scripts_orthology/treepaths.txt")):
        treepath = str(treepath)
        treepath = treepath.strip()
        yield treepath

def process_tree(treepath):
    ''' processes a tree to extract orthology relationships between target taxid and the rest
     of species, organized by orthology type and species code '''
    treepath = str(treepath)
    treepath = treepath.rstrip()
    t = PhyloTree(treepath, sp_naming_function=get_species)
    # traverse all leaves in tree file and get taxid
    leaf_count = 0
    for leaf in t:
        leaf_count += 1
        tax = int(leaf.name.split(".", 1)[0])

        #get scientific name and convert taxid from int to str
        sci_name = names.get(tax)
        leaf.taxid = str(tax)

        #rename leaves names
        try:
            good_name = "%s" %(conversion[leaf.name][0])
        except:
            good_name = leaf.name
           
        
        good_name = re.sub("[ |\t,:)(;\n\]\[]+", "_", good_name)
        leaf.good_name = good_name
        

    #obtain cluster name from tree file path
    clus_name = os.path.split(treepath)[-1].replace(".fa.final_tree.nw", "")
    try:
        base_name = conversion[clus_name][0].replace('|', '_')
    except:
        base_name = clus_name[0]
    t.dist = 0

    #colapses plat specific
    node2content = t.get_cached_content()
    target_species = set([target_taxid])
    
    def is_sp_specific(_node):
        _species = set([_leaf.species for _leaf in  node2content[_node]])
        if not (_species - target_species):
            return True
        return False

    #traverse only lamprey leaves
    if collapse == 'yes':
        for n in t.get_leaves(is_leaf_fn=is_sp_specific):
            if n.children:
                for ch in n.get_children():
                    ch.detach()
                n.taxid = target_taxid
                n.name = "%s" %('|'.join([_lf.name for _lf in node2content[n]]))
                n.good_name = "{%s}" %('|'.join([_lf.good_name for _lf in node2content[n]]))
         
    #set outgroup
    outgroup = t.get_midpoint_outgroup()
    try:
        t.set_outgroup(outgroup)
    except:
        if len(t) == 1:
            return
        else:
            raise
    
    node2content = t.get_cached_content()
   
    event_lines = []
    for ev in t.get_descendant_evol_events():
        if ev.etype == "S":
            
            source_seqs = node2content[ev.node.children[0]]
            ortho_seqs = node2content[ev.node.children[1]]
            
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
           
            
            #co_orthologs is a list with lamprey seed in source_seqs
            co_orthologs = [leaf.good_name for leaf in source_seqs if leaf.taxid == str(target_taxid)]
            co_orthologs.sort()
            
            #orthologs is a list of all ortho_seqs names
            orthologs = defaultdict(set)
            for leaf in ortho_seqs:
                sp = int(leaf.taxid)
                orthologs[sp].add(leaf.good_name)  
            
            if len(co_orthologs) == 1:
                _otype = "one-to-"
            else:
                _otype = "many-to-"
        
            for sp, orth in orthologs.iteritems():
                if len(orth) == 1:
                    otype = _otype + "one"
                else:
                    otype = _otype + "many"
                    
                event_lines.append('\t'.join([','.join(co_orthologs), otype , str(sp), names[sp], ','.join(sorted(orth)),'\n']))
    return event_lines
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--collapse', dest='collapse', default='yes', type=str, choices=["yes", "no"])
    parser.add_argument('--target', dest='target_taxid', type=str)
    parser.add_argument('--cpu', dest='cpu', type=int, default=5)
    args = parser.parse_args()

    target_taxid = args.target_taxid
    collapse = args.collapse
    CPU = args.cpu

    # Set up taxonomy information (currently optimized for Lamprey analysis)
    conversion = {}
    taxids = set()
    rooting_dict = {'8364': 1, '13735': 1, '7918': 1, '7227': 4, '7574': 4, '8090': 1, '7719': 2,
                    '7739': 2, '6183': 4, '8839': 1, '7757': 0, '9606': 1, '7668': 3, '10224': 3,
                    '9315': 1, '7868': 1, '10090': 1, '6500': 4, '7955': 1, '6359': 4, '7740': 2,
                    '45351': 5, '6239': 4, '31033': 1, '7897': 1}

    for key in rooting_dict:
        taxids.add(int(key))

    for line in open("/data/spongilla/conversion.tsv"):
        code, realname, annotations = line.split("\t")
        conversion[code] = [realname, annotations]

    ncbi = NCBITaxa()
    names = ncbi.get_taxid_translator(taxids)
    names[target_taxid]=target_taxid

    pool = Pool(CPU)
    with open('orthology_events.tsv', 'w') as ORTHOLOGS:
        for event_lines in pool.imap_unordered(process_tree, iter_tree_paths()):
            for line in event_lines:
                print >>ORTHOLOGS, line
