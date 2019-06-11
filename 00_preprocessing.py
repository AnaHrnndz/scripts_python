import sys
import re
#from seqgroup import SeqGroup
from ete3 import SeqGroup
from ete3 import NCBITaxa

ncbi=NCBITaxa()
seq_file= SeqGroup(sys.argv[1])
OUT = open(sys.argv[2], 'w')
proc_type=str(sys.argv[3])


def proces_prote_check_taxid(seq_file):
	for seq_num, (name, seq, _) in enumerate(seq_file):
		try:
			name_id=name.split('.')[0]
			taxid2name=ncbi.get_taxid_translator([name_id])
			name = re.sub('\|', '_', name)
			name = re.sub('\:|\;|\*|,|\"|\[|\]|\(|\)|\/', '_', name)
			name = re.sub(r"\\", '_', name)
			seq = re.sub('\*|"|\+', '',seq)
			seq = re.sub('Z|B|J', 'X',seq)
			print (">%s\n%s" %(name, seq), file= OUT)
		except:
			print ('ERROR in', name)
	OUT.close()

def proces_prote_add_taxid(seq_file):
	for seq_num, (name, seq, _) in enumerate(seq_file):
		try:
			taxid=str(sys.argv[4])
			name = re.sub('\|', '_', name)
			name = re.sub('\:|\;|\*|,|\"|\[|\]|\(|\)|\/', '_', name)
			name = re.sub(r"\\", '_', name)
			name = '%s.%s' %(taxid, name)
			seq = re.sub('\*|"|\+', '',seq)
			seq = re.sub('Z|B|J', 'X',seq)
			print(">%s\n%s" %(name, seq), file = OUT)
		except:
			print ('ERROR in', name)
	OUT.close()

    
def proces_metag(seq_file):
	for seq_num, (name, seq, _) in enumerate(seq_file):
		try:
			name = name.split(' ')[0]
			name = re.sub('\|', '_', name)
			name = re.sub('\:|\;|\*|,|\"|\[|\]|\(|\)|\/', '_', name)
			name = re.sub(r"\\", '_', name)
			seq = re.sub('\*|"|\+', '',seq)
			seq = re.sub('Z|B|J', 'X',seq)
			print >>OUT, '>%s\n%s'%(name, seq)
		except:
			print ('ERROR in', name)
	OUT.close()


if proc_type=='metag':
    	proces_metag(seq_file)
elif proc_type=='check_taxid':
	proces_prote_check_taxid(seq_file)
elif proc_type=='add_taxid':
	proces_prote_add_taxid(seq_file)
else:
	print ('ERROR in processing type')
