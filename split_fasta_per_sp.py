from ete3 import SeqGroup

sp_mem = {}
in_fasta = SeqGroup('/home/plaza/research/dom_walk/raw/COG0484.faa')

for num, (name, seq, _) in enumerate(in_fasta):
    sp = name.split('.')[0]
    if sp not in sp_mem:
        sp_mem[sp] = []
    sp_mem[sp].append(name)
    
print ('writing fastas per sp')
for k, val in sp_mem.items():
    out_fasta = open('/home/plaza/research/dom_walk/analysis/fasta_per_sp/'+k+'.faa', 'w')
    for seq_name in val:
        print (">%s" %(seq_name), file = out_fasta)
        print (in_fasta.get_seq(seq_name), file =out_fasta)
    out_fasta.close()
