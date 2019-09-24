import re
from ete3 import SeqGroup, Tree
import sys


tree_file = sys.argv[1] # in newick format
original_fasta = SeqGroup(sys.argv[2])
pruned_fasta = open(sys.argv[3], 'w')
star_target = str(sys.argv[4])
end_target = str(sys.argv[5])


tree=Tree(tree_file)
R = tree.get_midpoint_outgroup()
tree.set_outgroup(R) 
    
name_list=[]  
for num, leaf in enumerate(tree):
    name_list.append(leaf.name)
    if star_target == leaf.name:
        star_pos=num
    if end_target == leaf.name:
        end_pos=num
        
pruned_list=name_list[star_pos:(end_pos+1)]
print pruned_list

#for ele in pruned_list:
#    print >>pruned_fasta,">%s\n%s"%(ele, original_fasta.get_seq(ele))
    
for ele in name_list:
    if ele not in pruned_list:
        print >>pruned_fasta,">%s\n%s"%(ele, original_fasta.get_seq(ele))

