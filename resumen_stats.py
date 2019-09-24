import sys
import re
import argparse


eggnog_expanded = "./eggnog_clustering/results_eggnog5.0/mcl_no_weights/euk/compare_euk_extended_nokogs.tsv"
eggnog_to_mmseqs =  "./eggnog_clustering/results_eggnog5.0/mcl_no_weights/euk/compare_to_euk.tsv"
mmseq_collapsed = "./eggnog_clustering/results_eggnog5.0/mcl_no_weights/euk/mcl_to_euk.tsv"
mmseqs_clus_members = "./eggnog_clustering/results_eggnog5.0/mcl_no_weights/clus_mem.tsv"
eggnog_cls = "/scratch/plaza/soft/eggnog_clustering/results_eggnog5.0/Euk_members.tsv"
hhblits_table = "/scratch/plaza/soft/eggnog_clustering/results_eggnog5.0/mcl_no_weights/hhblits_hits_1.abc"


def get_list(dict_query, column_query, condition, value):
    save_name = set()
    for key, info_d in dict_query.items():
        for info, real_val in info_d.items():
            if info == column_query:
                if condition == 'equal':
                    if int(real_val) == int(value):
                        save_name.add(key)
                elif condition == 'menor':
                    if int(real_val) < int(value):
                        save_name.add(key)
                elif condition == 'mayor':
                    if int(real_val) > int(value):
                        save_name.add(key)
            else:
                continue
    return list(save_name)


def get_info(list_clusters, egg_exp_dict, egg_to_mmseq, e_clusters, mmseq_mem_dict):

    type_clus = list_clusters[0]
    match_MMSEQ = re.match( 'MMSEQS_[0-9]*', type_clus )
    match_MCL = re.match( 'MCL_[0-9]*', type_clus)
    

    if match_MMSEQ or match_MCL:
        for mem in list_clusters:
            count = 0
            try:
                info_member = egg_exp_dict[mem]
                print (mem, '\t'.join(info_member.values()))
                related_mmseqs = set()
                for egg_clus, info_val in egg_to_mmseq.items():
                    for key, val in info_val.items():
                        if key == 'mmseqs_best_match':
                            if val == mem:
                                iter_name = mem.split('.')[0]+'_'+str(count)
                                info_egg = egg_to_mmseq[egg_clus]
                                related_mmseqs_egg = map(str.strip, info_egg['mmseq_related'].split(","))
                                for rel in related_mmseqs_egg:
                                    related_mmseqs.add(rel)
                                print ('\t', egg_clus, iter_name, '\t'.join(info_egg.values()))
                                count +=1
                print('EGGNOGS PRESENTS IN', len(related_mmseqs), 'MMSEQS/MCL CLUSTERS','\n')
                            
            except:
                clus_members = mmseq_mem_dict[mem]
                print(mem, 'CLUSTER NOT PRESENT IN EGGNOG BACT')
                print('MEMBERS:', ','.join(clus_members), '\n')

    else:
        for mem in list_clusters:
            info_member = egg_to_mmseq[mem]
            name_mmseqs_best = info_member['mmseqs_best_match']
            info_mmseqs = egg_exp_dict[name_mmseqs_best]
            e_members = e_clusters[mem]
            m_members = mmseq_mem_dict[name_mmseqs_best]
            print (mem, '\t'.join(info_member.values()), ','.join(e_members.difference(m_members)))
            print ('\t', name_mmseqs_best, '\t'.join(info_mmseqs.values()))
            related_mmseqs = map(str.strip, info_member['mmseq_related'].split(','))

            for related in related_mmseqs:
                if related != str(name_mmseqs_best):
                    try:
                        info_mmseqs_related = egg_exp_dict[related]
                        print ('\t', related, '\t'.join(info_mmseqs_related.values()), '\n')
                    except:
                         print('\t', related, 'CLUSTER NOT PRESENT IN EGGNOG BACT', '\n')

    return 




def get_info_seq(seq_names, e_clust, m_clust):
    print ('SEARCH SEQ', len(seq_names))
    for seq in seq_names:
        for e_key, e_mem in e_clust.items():
            if seq in e_mem:
                print(seq, e_key)
        for m_key, m_mem in m_clust.items():
            if seq in m_mem:
                print (seq, m_key)

            
def check_hhblits(mmseq_names, hhblits_dict ):
    print ('SEARCH HHLITS HITS:', mmseq_names)
    name_1 = mmseq_names[0]
    name_2 = mmseq_names[1]
    print(name_1)
    info_name_1 = hhblits_dict[name_1]
    for hit in info_name_1:
        if hit['hit'] == name_2:
            print(hit)
    
    print(name_2)
    info_name_2 = hhblits_dict[name_2]
    for hit in info_name_2:
        if hit['hit'] == name_1:
            print(hit)
            
    
    # name_2 = mmseq_names[1]
    # print(name_2)
    # info_name_2 = hhblits_dict[name_2]
    # for k, val in info_name_2.items():
        # if info_name_1['hit'] == mmseq_names[0]:
            # get_hit_2 = hhblits_dict[name_2]
            # print (get_hit_2)
            # break
    


eggnog_expanded_dict = {}
for line in open(eggnog_expanded, 'r'):
    info_dict = {}
    fields = line.rstrip().split("\t")
    mmseq_name = fields[0]
    info_dict['repeat'] = fields[1]
    info_dict['seq_miss'] = fields[2]
    info_dict['num_seq_egg'] = fields[3]
    info_dict['num_seq_mmseq']= fields[4]
    info_dict['num_sp_egg'] = fields[5]
    info_dict['num_sp_mmseq'] = fields[6]
    eggnog_expanded_dict[mmseq_name]=info_dict


eggnog_to_mmseqs_dict = {}
for line in open(eggnog_to_mmseqs, 'r'):
    info_dict = {}
    fields = line.rstrip().split("\t")
    egg_name = fields[0]
    if fields[1] == 'None':
        info_dict['mmseqs_best_match'] = '0'
    else:
        info_dict['mmseqs_best_match'] = fields[1]
    info_dict['num_seq_egg'] = fields[2]
    info_dict['num_seq_mmseq'] = fields[3]
    if fields[4] == 'None':
        info_dict['num_seq_miss'] = '0'
    else:
        info_dict['num_seq_miss'] = fields[4]
    info_dict['perct'] = fields[5]
    info_dict['num_splits'] = fields[6]
    info_dict['num_sp_egg'] = fields[7]
    info_dict['num_sp_mmseq'] = fields[8]
    info_dict['num_sp_miss'] = fields[9]
    try:
        info_dict['mmseq_related'] = fields[10]
    except:
        info_dict['mmseq_related'] = '-'
    eggnog_to_mmseqs_dict[egg_name] = info_dict

mmseq_coll_dict = {}

for line in open(mmseq_collapsed,'r'):
    info_dict = {}
    fields = line.rstrip().split("\t")
    mmseq_name = fields[0]
    info_dict['list_eggnogs'] = list(fields[1])
    if len(info_dict['list_eggnogs']) == 1:
        check_COG = info_dict['list_eggnogs']
        match_COG = re.match( 'COG[0-9]*', check_COG)
        if not match_COG:
            mmseq_coll_dict[mmseq_name] = info_dict
    else:
        mmseq_coll_dict[mmseq_name] = info_dict

mmseq_mem_dict = {}
seqs = set()
for line in open(mmseqs_clus_members,'r'):
    fields = line.rstrip().split("\t")
    mmseq_name = fields[0].rstrip()
    members = map(str.strip, fields[1].split(","))
    mmseq_mem_dict[mmseq_name] = set(members)

e_clusters = {}
for line in open(eggnog_cls):
    fields = line.rstrip().split("\t")
    cluster = fields[1]
    members = map(str.strip, fields[4].split(","))
    e_clusters[cluster] = set(members)

    #supervised OG not included
    # match_COG = re.match( 'COG[0-9]*', cluster)
    # if match_COG:
        # continue
    # else:
        #e_clusters[cluster] = set(members)

hhblits_dict = {}
for line in open(hhblits_table):
    fields = line.rstrip().split("\t")

    name_1 = fields[0]
    name_2 = fields[1]
    hit_dict = {}
    hit_dict['hit'] = name_2
    #hit_dict['eval'] = fields[2]
    #hit_dict['score'] = fields[3]
    #hit_dict['cov'] = fields[4]
    
        
    if name_1 in hhblits_dict:
        hhblits_dict[name_1].append(hit_dict)
    else:
        hhblits_dict[name_1] = list()
        hhblits_dict[name_1].append(hit_dict)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--query_dict', dest='q_dict', default = None, type=str, choices=["eggnog_expanded_dict", "eggnog_to_mmseqs_dict"], required = False)
    parser.add_argument('--query_col', dest='col_query', default = None, type=str, required = False)
    parser.add_argument('--condition', dest='condi', default = None, choices=["equal", "mayor", "menor"], type=str, required = False)
    parser.add_argument('--value', dest='user_val', default = None,  type=int, required = False)
    parser.add_argument('--cluster_name', dest = 'clus_query', type = str, required = False)
    parser.add_argument('--query_seq', dest = 'q_seq', type = str, required = False)
    parser.add_argument('--hhblits_hits', dest = 'hhblits', type = str, required = False)
    args = parser.parse_args()

    q_dict = args.q_dict
    col_query = args.col_query
    condi = args.condi
    user_val = args.user_val
    clus_query = args.clus_query
    q_seq = args.q_seq
    hhblits_hits = args.hhblits

    if hhblits_hits:
        mmseqs_names = list()
        try:
            mmseqs_names = list(map(str.strip, hhblits_hits.split(",")))
        except:
            print('NOT VALID MMSEQS NAMES')

        check_hhblits(mmseqs_names, hhblits_dict)
           

    elif q_seq:
        seq_names = list()
        try:
            seq_names = list(map(str.strip, q_seq.split(",")))    
        except:
            seq_names.append(str(q_seq))
        get_info_seq(seq_names, e_clusters, mmseq_mem_dict)

    else:
        if q_dict == 'eggnog_expanded_dict':
            names = get_list(eggnog_expanded_dict, col_query, condi, user_val)
            print (len(names), 'clusters in expanded eggnog')
    
        elif q_dict == 'eggnog_to_mmseqs_dict':
            names = get_list(eggnog_to_mmseqs_dict, col_query, condi, user_val)
            print (len(names), 'clusters')
        else: 
            names = list()
            try:
                names = list(map(str.strip, clus_query.split(",")))    
            except:
                names.append(str(clus_query))

        get_info(names, eggnog_expanded_dict, eggnog_to_mmseqs_dict, e_clusters, mmseq_mem_dict)
