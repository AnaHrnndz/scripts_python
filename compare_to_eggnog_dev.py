from collections import defaultdict
import hashlib
import numpy as np
from datetime import datetime
import argparse


def get_eggnog_splits(e_clusters, m_clusters, f_name):
    """
    e_clusters = {eclu_name: set(eclu_members), ...}
    m_clusters = {mclu_name: set(mclu_members), ...}
    """
    RESULTS_F = open(f_name, "w")

    # precalculate where does it appear each individual
    # sequence name (what clusters)
    # seq2ecluster = defaultdict(list)
    # for clu_name, members in e_clusters.iteritems():
    #     for seq in members:
    #         seq2ecluster[seq].append(clu_name)

    seq2mcluster = defaultdict(list)
    for clu_name, members in m_clusters.iteritems():
        for seq in members:
            seq2mcluster[seq].append(clu_name)

    print "members mapping to clusters done"

    # For each eggnog cluster, lets compare if it is split or fully
    # contained in an mmseqs cluster.
    for e_clu, e_members in e_clusters.iteritems():
        related_mmseqs_clusters = set()
        for member in e_members:
            related_mmseqs_clusters.update(seq2mcluster.get(member, []))

        best_match = None
        best_match_missing = None
        for m_clu in related_mmseqs_clusters:
            m_members = m_clusters[m_clu]
            missings = len((e_members - m_members) - genes_to_be_excluded)
            if best_match_missing is None or missings < best_match_missing: 
                best_match_missing = missings
		percent = (float(best_match_missing)/float((len(e_members))) )
                best_match = m_clu
            if missings == 0: 
		percent = (float(missings)/float((len(e_members))))
                break
        #print >>RESULTS_F, '\t'.join(map(str, [e_clu, best_match, len(e_members), len(m_members), best_match_missing]))
        #print >>RESULTS_F, '\t'.join(map(str, [e_clu, best_match, len(e_members), len(m_clusters[best_match]), best_match_missing]))
	print >>RESULTS_F, '\t'.join(map(str, [e_clu, best_match, len(e_members), len(m_clusters[best_match]), best_match_missing, percent, len(related_mmseqs_clusters), ','.join(related_mmseqs_clusters)]))
	RESULTS_F.flush()
    RESULTS_F.close()

def compared_clusters(e_clusters, m_clusters, f_name):
    """
    e_clusters = {eclu_name: set(eclu_members), ...}
    m_clusters = {mclu_name: set(mclu_members), ...}
    """
    RESULTS_F = open(f_name, "w")

    # precalculate where does it appear each individual
    # sequence name (what clusters)
    seq2ecluster = defaultdict(list)
    for clu_name, members in e_clusters.iteritems():
        for seq in members:
            seq2ecluster[seq].append(clu_name)

    seq2mcluster = defaultdict(list)
    for clu_name, members in m_clusters.iteritems():
        for seq in members:
            seq2mcluster[seq].append(clu_name)

    print "members mapping to clusters done"

    # For each MMseqs clusters, lets compare content only against the subset of
    # eggnog clusters that had at least one sequence in common (avoids brute
    # force double loop).

    # eggnog cluster contained or overlapping with mmseqs clusters 
    for m_clu, m_members in m_clusters.iteritems():
        related_eggnog_clusters = set()
        for member in m_members:
            member_e_clusters = seq2ecluster.get(member, [])
            related_eggnog_clusters.update(member_e_clusters)
            
        best_match = None
        e_missing_min = float("inf")
        for e_clu in related_eggnog_clusters:
            common = len(m_members & e_clusters[e_clu])
            m_missing = len(m_members - e_clusters[e_clu])
            e_missing = len(e_clusters[e_clu] - m_members)
            print "we are at", e_clu
            if e_missing < e_missing_min:
                best_match = (e_clu, common, m_missing, e_missing)
                e_missing_min = e_missing            
    
        if best_match != None:
            print "best_match", best_match
            e_clu, common, m_missing, e_missing = best_match
            m_missing_percent = len(m_members - e_clusters[e_clu])/float(len(m_members))
            e_missing_percent = len(e_clusters[e_clu] - m_members)/float(len(e_clusters[e_clu]))
            
            print >> RESULTS_F, "\t".join(map(str, (m_clu, len(related_eggnog_clusters), e_clu, len(m_members), len(e_clusters[e_clu]), common, m_missing, e_missing, m_missing_percent, e_missing_percent)))


parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='cls_tsv_file', type=str,
                   help='clusters tsv file')
parser.add_argument('-o', dest='comparison_file', type=str,
                   help='clusters fasta file')

args = parser.parse_args()

startTime = datetime.now()

eggnog_cls = "/scratch/pipeline_mmseqs/egg5/results/e5.luca_members.tsv"
mmseqs_cls = args.cls_tsv_file

m_clusters = defaultdict(list)

# getting all mmseqs clusters in a dict
for line in open(mmseqs_cls):
    fields = line.rstrip().split("\t")
    cluster_name = fields[0]
    members = map(str.strip, fields[1].split(","))
    if len(members) > 1:
        m_clusters[cluster_name] = set(members)
print '%d mmseqs clusters loaded' %(len(m_clusters))

# getting all eggnog clusters in a dict
e_clusters = {}
for line in open(eggnog_cls):
    fields = line.rstrip().split("\t")
    cluster = fields[1]
    members = map(str.strip, fields[4].split(","))
    e_clusters[cluster] = set(members)

#for line in open("/scratch/huerta/eggnog_clustering/genes_in_multiple_eggnog_clusters"):
#    genes_to_be_excluded.add(line.rstrip())    
genes_to_be_excluded=set()    
print '%d eggnog clusters loaded' %(len(e_clusters))


# converting to hashes to find identical clusters
m_hashes = {}
e_hashes = {}
for c1 in m_clusters.keys():
    c1_string = ','.join(sorted(m_clusters[c1]))
    m_hashes[hashlib.md5(c1_string).hexdigest()] = c1
for c2 in e_clusters.keys():
    c2_string = ','.join(sorted(e_clusters[c2]))
    e_hashes[hashlib.md5(c2_string).hexdigest()] = c2
print "converting to hashes done"

h1_set = set(m_hashes.keys())
h2_set = set(e_hashes.keys())

common = h1_set & h2_set
print "identical clusters", len(common)

#for h in common:
#    del m_clusters[m_hashes[h]]
#    del e_clusters[e_hashes[h]]

print 'common clusters deleted from set'
print len(m_clusters), len(e_clusters)

t1 = datetime.now()

get_eggnog_splits(e_clusters, m_clusters, args.comparison_file)


print "comparison done in ", datetime.now() - t1

print "script running for", datetime.now() - startTime
