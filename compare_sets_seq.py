import sys
import re

eggnog_cls = "/scratch/plaza/soft/eggnog_clustering/results_eggnog5.0/arc_members.tsv"
mmseqs_cls = sys.argv[1]
compa_tab = sys.argv[2] 



#load eggnog clusters
e_clusters = {}
seq_to_eclus= {}
for line in open(eggnog_cls):
    fields = line.rstrip().split("\t")
    cluster = fields[1]
    members = list(map(str.strip, fields[4].split(",")))
    
    for m in members:
        seq_to_eclus[m] = cluster
    #supervised OG not included
    match_COG = re.match('arCOG[0-9]*', cluster)
    
    if match_COG:
        continue
    else:
        e_clusters[cluster] = set(members)

#load mmseqs/mcl clusters
m_clusters= {}
for line in open(mmseqs_cls):
    fields = line.rstrip().split("\t")
    cluster_name = fields[0].rstrip()
    members = list(map(str.strip, fields[1].split(",")))
    clu_type = list()
    if len(members) > 1:
        for m in members:
            try:
                clu_type.append(seq_to_eclus[m])
            except:
                clu_type.append('NO_EGG')

        if all(re.match( 'arCOG[0-9]*', typ) for typ in clu_type):
            continue
        else:
            m_clusters[cluster_name] = set(members)


#load collapsed table: all eggnog clusters that match with the same mmseqs/mcl 
# Ej: {MMSEQS_1654165: eggnog1, eggnog2, eggnog3}....
egg_same_mmseqs = {}
for line in open(compa_tab):
    fields = line.rstrip().split("\t")
    mmseqs_clus = fields[0]
    egg_clus_mem = list(map(str.strip, fields[1].split(",")))
    egg_clus_mem_set = set()
    for mem in egg_clus_mem:
        match_COG = re.match( 'arCOG[0-9]*', mem)
        if match_COG:
            continue
        else:
            egg_clus_mem_set.add(mem)
    egg_same_mmseqs[mmseqs_clus] = egg_clus_mem_set


#merge all seqs for each eggnog cluster with the same mmseqs match
egg_same_mmseqs_seqs = {}
for clus, egg_members in egg_same_mmseqs.items():
    if egg_members != None:
        seqs_egg_members = set()
        for mem in egg_members:
            for seq in e_clusters[mem]:
                seqs_egg_members.add(seq)
        egg_same_mmseqs_seqs[clus]=seqs_egg_members


#Compare all seqs from all eggnog cluster with the same mmseqs/mcl match , with all seqs that include that mmseqs/mcl match
OUT = open("compare_arq_extended_noarCOG.tsv", 'w')
for mm_clus, egg_mem in egg_same_mmseqs.items():
    if mm_clus != 'None':
        try: 
            mmseq_seqs = m_clusters[mm_clus]
        except:
            continue

        all_seqs = set(egg_same_mmseqs_seqs[mm_clus])
        egg_sp = set ()
        for seq in all_seqs:
            sp = seq.split('.')[0]
            egg_sp.add(sp)
        mmseq_seqs = m_clusters[mm_clus]
        mmseq_seqs_set = set(mmseq_seqs)
        mmseqs_sp = set()
        for seq in mmseq_seqs_set:
            sp = seq.split('.')[0]
            mmseqs_sp.add(sp)
        d = all_seqs.difference(mmseq_seqs_set)
        OUT.write('\t'.join(map(str,(mm_clus, len(egg_mem),  len(d), len(all_seqs),len(mmseq_seqs_set),len(egg_sp), len(mmseqs_sp),'\n'))))

OUT.close()


