from Bio import SeqIO
from collections import defaultdict

print('load fasta file')
sequences = SeqIO.index('/scratch/plaza/soft/eggnog_clustering/results_eggnog5.0/clean_fasta.faa', 'fasta')
print ("parsing fasta done... now writing fasta files") 
mcl_cluster = open ('./2QVH7_members.tsv', 'r')

clusters = defaultdict(set)
for line in mcl_cluster:
    fields = line.rstrip().split("\t")
    cluster = fields[0]
    members = fields[1].split(',')
    for sequence in members:
        sp = sequence.split('.')[0]
        clusters[sp].add(sequence)



for sp, members in clusters.items():
    out_fasta = open(sp+'.faa', 'w')
    for seq_name in members:
        try:
            print (">%s" %(seq_name), file = out_fasta)
            print (sequences[seq_name].seq, file =out_fasta)
        except:
            sys.stderr.write("sequence missing", seq_name, "in cluster", clu_name)
    out_fasta.close()
