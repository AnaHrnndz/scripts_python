import re
from ete3 import SeqGroup, Tree
import sys

alg_file = sys.argv[1] # in fasta format
tree_file = sys.argv[2] # in newick format

alg = SeqGroup(alg_file)
for k,v in alg.name2id.items():
    # converts ilegal newick chars from alg names. 
    # Comment this line if not necessary
    k = re.sub('[:,();]','_', k)
    alg.name2id[k] = v

tree = Tree(tree_file)
for leaf in tree:
    print (">%s\n%s" %(leaf.name, alg.get_seq(leaf.name)))