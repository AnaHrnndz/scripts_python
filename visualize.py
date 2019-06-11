import re
from string import lower
import hashlib
import time
import os

from ete_dev import (add_face_to_node, AttrFace, TextFace, TreeStyle,
                     SequenceFace, random_color, SeqMotifFace)

class TreeImage(object):
     def ly_basic(self, node):
        if node.is_leaf():
            node.img_style['size'] = 0
        else:
            node.img_style['size'] = 0
            node.img_style['shape'] = 'square'

        node.img_style['hz_line_width'] = 1
        node.img_style['vt_line_width'] = 1
                    
     def ly_leaf_names(self, node):
         if node.is_leaf():
             if node.species != 0 and self.ncbitaxa:
                 spF = TextFace("%s" %getattr(node, "spname", ""), fsize=8, fgcolor='#777777', ftype='Helvetica')
                 geneF = TextFace("(%s)" %node.name[:50], fsize=10, fgcolor='#444444', fstyle='italic', ftype='Helvetica')
                 
                 add_face_to_node(spF, node, column=0, position='branch-right')
                 add_face_to_node(geneF, node, column=1, position='branch-right')
                 
             else:
                 geneF = TextFace("%s" %node.name, fsize=8, fgcolor='#777777', ftype='Helvetica')
                 add_face_to_node(geneF, node, column=0, position='branch-right')

     def ly_supports(self, node):
         if not node.is_leaf():
             supFace = TextFace("%0.3g" %(node.support), fsize=7, fgcolor='indianred')
             add_face_to_node(supFace, node, column=0, position='branch-top')
                
     def ly_tax_labels(self, node):
         if node.is_leaf():
             c = self.LABEL_START_COL
             node_lin = getattr(node, "named_lineage", [])
             for tname in node_lin:
                  if tname.lower() in self.TRACKED_CLADES:
                       linF = self.LIN2FACE[tname.lower()]
                       # linF = TextFace(tname, fsize=10, fgcolor='white')
                       # linF.margin_left = 3
                       # linF.margin_right = 2
                       # linF.background.color = self.LIN2COLOR[tname.lower()]
                       add_face_to_node(linF, node, c, position='aligned')
                       c += 1
            
             for n in xrange(c, len(self.TRACKED_CLADES)):
                  add_face_to_node(TextFace('', fsize=10, fgcolor='slategrey'), node, c, position='aligned')
                  c+=1

     def ly_condensed_alg(self, node):
         if node.is_leaf():
             if 'sequence' in node.features:
                 seqFace = SeqMotifFace(node.sequence, [],
                                        intermotif_format="line",
                                        seqtail_format="line", scale_factor=1)
                 add_face_to_node(seqFace, node, 10, aligned=True)

     def ly_block_alg(self, node):
          if node.is_leaf():
             if 'sequence' in node.features:
                 # [10, 100, "[]", None, 10, "black", "rgradient:blue", "arial|8|white|domain Name"],
                 motifs = []
                 last_lt = None
                 for c, lt in enumerate(node.sequence):
                     if lt != '-':
                         if last_lt is None:
                             last_lt = c
                         if c+1 == len(node.sequence):
                             start, end = last_lt, c
                             motifs.append([start, end, "()", 0, 12, "slategrey", "slategrey", None])
                             last_lt = None

                     elif lt == '-':
                         if last_lt is not None:
                             start, end = last_lt, c-1
                             motifs.append([start, end, "()", 0, 12, "grey", "slategrey", None])
                             last_lt = None
                 seqFace = SeqMotifFace(node.sequence, motifs,
                                        intermotif_format="line",
                                        seqtail_format="line", scale_factor=1)
                 add_face_to_node(seqFace, node, 10, aligned=True)
    
     def __init__(self, tree, alg_type=None, ncbitaxa=True):
         self.tree = tree
         self.tree.dist = 0
         self.svg = None
         self.layouts = None
         self.treestyle = TreeStyle()
         self.alg_type = alg_type
         
         ts = self.treestyle
         ts.draw_aligned_faces_as_table = False
         ts.draw_guiding_lines = False
         ts.show_leaf_name = False
         ts.show_branch_support = False
         ts.tree_width = 250
         ts.layout_fn = [self.ly_basic, self.ly_leaf_names, self.ly_supports]
         self.ncbitaxa = ncbitaxa
         
         if ncbitaxa:
              ts.layout_fn.append(self.ly_tax_labels)

         if self.alg_type == 'condensed':
              ts.layout_fn.append(self.ly_condensed_alg)
         elif self.alg_type == 'block':
              ts.layout_fn.append(self.ly_block_alg)
              
         self.TRACKED_CLADES = set()         
         if ncbitaxa: 
              tree.set_species_naming_function(spname)
              annotate_tree_with_ncbi(tree)

              for lf in tree:
                   for lname in getattr(lf, 'named_lineage', [])[2:9]:
                        self.TRACKED_CLADES.add(lname.lower())

              self.TRACKED_CLADES.update(map(lower, ["Eukaryota", "Viridiplantae",  "Fungi", "birds", "Basidiomycota", "Ascomycota",
                                          "Alveolata", "Metazoa", "Stramenopiles", "Rhodophyta", "mammalia", "Primates",
                                          "Amoebozoa", "Crypthophyta", "Bacteria","Archaea", "Arthropoda", "rodentia", "laurastheria",
                                          "Alphaproteobacteria", "Betaproteobacteria", "Cyanobacteria",
                                          "Gammaproteobacteria",]))




              colors = random_color(num=len(self.TRACKED_CLADES), s=0.45)
              self.LIN2COLOR = dict([(ln, colors[i]) for i, ln in enumerate(sorted(self.TRACKED_CLADES))])
              self.LIN2FACE = {}
              for tname in self.TRACKED_CLADES:
                   linF = TextFace(tname, fsize=10, fgcolor='white', tight_text=False, ftype="Arial")
                   linF.margin_left = 3
                   linF.margin_right = 2
                   linF.inner_background.color = self.LIN2COLOR[tname.lower()]
                   self.LIN2FACE[tname.lower()] = linF              

              
         self.LABEL_START_COL = 100
         self.ALG_START_COL = len(self.TRACKED_CLADES)+11
         #tree.ladderize()
         self.svg, self.img_map = tree.render("%%return", tree_style=ts)

              
def annotate_tree_with_ncbi(tree):
    from ete_dev.ncbi_taxonomy import ncbiquery as ncbi
    ncbi.connect_database()
    ncbi.annotate_tree(tree, attr_name='species')
    
def spname2(name):
    m = re.search('\{([^}]+)\}', name)
    if m:
        return m.groups()[0]
    else:
        return name.split('|')[0].strip().replace('_', ' ')
        taxid = name.split('.', 1)[0]

        tax2name = ncbi.get_taxid_translator([taxid])
        if int(taxid) not in tax2name:
            print 'name', name        , taxid, tax2name
        return tax2name.get(int(taxid), taxid)

def spname(name):
    try:
        return int(name.split('.', 1)[0].strip())
    except Exception:
        return 0
    
