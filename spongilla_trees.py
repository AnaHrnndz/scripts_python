from ete3 import  Tree, PhyloTree, TreeStyle, add_face_to_node, SeqMotifFace, TextFace, faces,  AttrFace
from ete3 import NCBITaxa
from ete3.treeview.svg_colors import random_color
import sys

ncbi=NCBITaxa()

def layout(node):
    node.img_style["size"] = 0
    # add species name
    # add protein name
    # add original sequence name
    
    if node.is_leaf():
        
        taxid=(node.name.split('.')[0])        
        
        #add conversion name or seq_name for spongilla
        if taxid!='6055':
            original_name=orig_name[node.name]
            original_name=original_name.split('.')[1]
            oriNameFace = faces.TextFace(original_name, fsize=14)
            oriNameFace.margin_right=10 
            oriNameFace.margin_left=10
            add_face_to_node(oriNameFace, node, column=1,position="branch-right")
        
        elif taxid =='6055':
            #highlight spongilla seed names
            if node.name in seed_names:               
                oriNameFace = faces.TextFace((node.name), fsize=14)
                oriNameFace.margin_right=10 
                oriNameFace.margin_left=10
                oriNameFace.background.color="Aqua" 
                add_face_to_node(oriNameFace, node, column=1,position="branch-right")
            else:
                oriNameFace = faces.TextFace((node.name), fsize=14)
                add_face_to_node(oriNameFace, node, column=1,position="branch-right")
            
        #add predicted name (eggnog 4.5)
        try:
            pred_name=predict_name[node.name]
            predNameFace = faces.TextFace(pred_name,fgcolor = "purple" , fsize=14)
            add_face_to_node(predNameFace, node, column=2, position="branch-right" )
        except:
            predNameFace = faces.TextFace('---',fgcolor="purple", fsize=14)
            add_face_to_node(predNameFace, node, column=2, position="branch-right")
            
                   
        #add alignment          
        seqFace = SeqMotifFace(node.sequence, gap_format="blank")
        add_face_to_node(seqFace, node, column=0, position="aligned")
        
        #change node.name for species name 
        sp_name = ncbi.get_taxid_translator([(node.name.split('.')[0])])
        node.name=sp_name[int(taxid)]
        
        #add Jake colors
        if node.name.startswith("Homo"):
            # Add an static face that handles the node name
            N = AttrFace("name", fsize=14, fgcolor="red")
            faces.add_face_to_node(N, node, column=0)

        elif node.name.startswith("Spongilla"):
            N = AttrFace("name", fsize=14, fgcolor="green")
            faces.add_face_to_node(N, node, column=0)
            
        elif node.name.startswith("Sycon"):
            N = AttrFace("name", fsize=14, fgcolor="green")
            faces.add_face_to_node(N, node, column=0)

        elif  node.name.startswith("Amphimedon"):
            N = AttrFace("name", fsize=14, fgcolor="green")
            faces.add_face_to_node(N, node, column=0)

        elif node.name.startswith("Oscarella"):
            N = AttrFace("name", fsize=14, fgcolor="green")
            faces.add_face_to_node(N, node, column=0)

        elif  node.name.startswith("Gallus"):
            N = AttrFace("name", fsize=14, fgcolor="red")
            faces.add_face_to_node(N, node, column=0)

        elif  node.name.startswith("Branchiostoma"):
            N = AttrFace("name", fsize=14, fgcolor="red")
            faces.add_face_to_node(N, node, column=0)

        elif  node.name.startswith("Trichoplax"):
            N = AttrFace("name", fsize=14, fgcolor="orange")
            faces.add_face_to_node(N, node, column=0)

        elif node.name.startswith("Nematostella"):
            N = AttrFace("name", fsize=14, fgcolor="orange")
            faces.add_face_to_node(N, node, column=0)

        elif node.name.startswith("Hydra"):
            N = AttrFace("name", fsize=14, fgcolor="orange")
            faces.add_face_to_node(N, node, column=0)

        elif  node.name.startswith("Drosophila"):
            N = AttrFace("name", fsize=14, fgcolor="blue")
            faces.add_face_to_node(N, node, column=0)

        elif  node.name.startswith("Crassostrea"):
            N = AttrFace("name", fsize=14, fgcolor="blue")
            faces.add_face_to_node(N, node, column=0)

        else:
            N = AttrFace("name", fsize=14, fgcolor="black")
            faces.add_face_to_node(N, node, column=0)

t=sys.argv[1]
alg=sys.argv[2]
predic_table=sys.argv[3]
seed_table=sys.argv[4]
out_img_name=sys.argv[5]

tree_upp = PhyloTree(newick=t, alignment=alg, alg_format="fasta")
R = tree_upp.get_midpoint_outgroup()
tree_upp.set_outgroup(R)

ts=TreeStyle()
ts.show_leaf_name = False
ts.tree_width = 2000

ts.layout_fn = layout

orig_name={} 
for line in open("/data/spongilla/conversion.tsv"):
    if line.strip() and not line.startswith("#"):
        prot_name=line.split("\t")[0]
        ori_name=line.split("\t")[1]
        orig_name[prot_name]=ori_name

predict_name={}
for line in open(predic_table):
    if line.rstrip() and not line.startswith("#"):
        prot_name=line.split("\t")[0]
        pred_name=line.split("\t")[1]
        predict_name[prot_name]=pred_name
        
seed_names=[]
for line in open(seed_table):
    line=line.rstrip()
    seed_names.append(line)
        
#t_klf_522.show(tree_style=ts)
tree_upp.render(out_img_name, tree_style=ts, h=14000, w=8000,units ="px",dpi=1000)
