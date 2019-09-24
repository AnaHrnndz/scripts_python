import sys
import re
import random, string
#from seqgroup import SeqGroup
from ete3 import SeqGroup
from ete3 import NCBITaxa

ncbi=NCBITaxa()
seq_file= SeqGroup(sys.argv[1])
OUT_fasta = open(sys.argv[2], 'w')
out_table = open(sys.argv[3], 'w')
proc_type=str(sys.argv[4])


def proces_prote_check_taxid(seq_file):
	for seq_num, (name, seq, _) in enumerate(seq_file):
		try:
			taxid=name.split('.')[0]
			taxid2name=ncbi.get_taxid_translator([taxid])
			code=''.join(random.choices(string.ascii_letters + string.digits, k=5))
			code_name = '.'.join([taxid,code])
		#	name = re.sub('\|', '_', name)
		 #	name = re.sub('\:|\;|\*|,|\"|\[|\]|\(|\)|\/', '_', name)
		#	name = re.sub(r"\\", '_', name)
			seq = re.sub('\*|"|\+', '',seq)
			seq = re.sub('Z|B|J', 'X',seq)
			print (">%s\n%s" %(code_name, seq), file= OUT_fasta)
			print ("%s\t%s" %(code_name, name), file = out_table)
		except:
			print ('ERROR in', name)
	OUT_fasta.close()

def proces_prote_add_taxid(seq_file):
	for seq_num, (name, seq, _) in enumerate(seq_file):
		try:
			taxid=str(sys.argv[4])
			code=''.join(random.choices(string.ascii_letters + string.digits, k=5))
			code_name = '.'.join([name_id,code])
			#name = re.sub('\|', '_', name)
			#name = re.sub('\:|\;|\*|,|\"|\[|\]|\(|\)|\/', '_', name)
			#name = re.sub(r"\\", '_', name)
			#name = '%s.%s' %(taxid, name)
			seq = re.sub('\*|"|\+', '',seq)
			seq = re.sub('Z|B|J', 'X',seq)
			print(">%s\n%s" %(code_name, seq), file = OUT_fasta)
			print ("%s\t%s" %(code_name, name), file = out_table)
		except:
			print ('ERROR in', name)
	OUT_fasta.close()

    
def proces_metag(seq_file):
	for seq_num, (name, seq, _) in enumerate(seq_file):
		try:
			code_name=''.join(random.choices(string.ascii_letters + string.digits, k=5))
			#name = name.split(' ')[0]
			#name = re.sub('\|', '_', name)
			#name = re.sub('\:|\;|\*|,|\"|\[|\]|\(|\)|\/', '_', name)
			#name = re.sub(r"\\", '_', name)
			seq = re.sub('\*|"|\+', '',seq)
			seq = re.sub('Z|B|J', 'X',seq)
			print >>OUT_fasta, '>%s\n%s'%(code_name, seq)
			print ("%s\t%s" %(code_name, name), file = out_table)
		except:
			print ('ERROR in', name)
	OUT_fasta.close()


if proc_type=='metag':
    	proces_metag(seq_file)
elif proc_type=='check_taxid':
	proces_prote_check_taxid(seq_file)
elif proc_type=='add_taxid':
	proces_prote_add_taxid(seq_file)
else:
	print ('ERROR in processing type')
