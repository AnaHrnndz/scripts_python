from ete3 import  Tree, PhyloTree, TreeStyle, add_face_to_node, SeqMotifFace, TextFace, faces,  AttrFace
from ete3 import NCBITaxa
from ete3.treeview.svg_colors import random_color
import sys

ncbi=NCBITaxa()

def layout(node):
    node.img_style["size"] = 0
    
    
    #add boostrap in not-leaf nodes
    if not node.is_leaf():
        boostFace=faces.TextFace(node.support, fgcolor= "grey", fsize=34)
        add_face_to_node(boostFace, node, column=0, position="branch-top" )
    
    if node.is_leaf():
        
        #get taxid and name 
        taxid=(node.name.split('.')[0])
        seq_name=(node.name.split('.')[1])
        
        #add predicted name (eggnog-mapper)
        try:
            pred_name=predict_name[node.name]
            predNameFace = faces.TextFace(pred_name,fgcolor = "salmon" , fsize=34)
            predNameFace.margin_right=20
            predNameFace.margin_left=20
            add_face_to_node(predNameFace, node, column=2, position="branch-right" )
        except:
            predNameFace = faces.TextFace('--',fgcolor="salmon", fsize=34)
            add_face_to_node(predNameFace, node, column=2, position="branch-right")
        
        #add sequence name (seq_name)
        seqNameFace = faces.TextFace(seq_name,fgcolor = "grey" , fsize=34)
        add_face_to_node(seqNameFace, node, column=1, position="branch-right" )
        
        #get scientific name from taxid
        sp_name = ncbi.get_taxid_translator([(node.name.split('.')[0])])
        node.name=sp_name[int(taxid)]
        
        #add alignment
        seqFace = SeqMotifFace(node.sequence, gap_format="blank")
        add_face_to_node(seqFace, node, column=3, position="aligned")
        
        
        #add specific faces to nodes: Example metazoa nodes in blue
        #lin= ncbi.get_lineage(taxid)
        #metazoa
        #if int('33208') in lin:
        #    N = AttrFace("name", fsize=34, fgcolor="blue")
        #    N.margin_left = 20
        #    N.margin_right= 20
        #    faces.add_face_to_node(N, node, column=0)
        
        
        
#tree file
t=sys.argv[1]

#alg file
alg=sys.argv[2]

#eggnog-mapper table
predic_table=sys.argv[3]
out_img_name=sys.argv[4]

#open tree and tree rooting (midpoint)
tree_upp = PhyloTree(newick=t, alignment=alg, alg_format="fasta")
R = tree_upp.get_midpoint_outgroup()
tree_upp.set_outgroup(R)

#load table eggnog-table
predict_name={}
for line in open(predic_table):
    if line.rstrip() and not line.startswith("#"):
        prot_name=line.split("\t")[0]
        pred_name=line.split("\t")[5]
        predict_name[prot_name]=pred_name


ts=TreeStyle()
ts.show_leaf_name = False
ts.tree_width = 2000

ts.layout_fn = layout

tree_upp.render(out_img_name, tree_style=ts, h=9000, dpi=300, units ="px")