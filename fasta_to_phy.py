import sys
import random, string
from ete3 import SeqGroup
from tempfile import NamedTemporaryFile

in_file=sys.argv[1]
transform_fasta=sys.argv[2]
out_file=sys.argv[3]
translate_table=open(sys.argv[4],'w')

alg = SeqGroup(in_file)
translate=open(transform_fasta, 'w')

for num,(name,seq,_) in enumerate(alg):
    taxid=name.split('.')[0]
    code=''.join(random.choices(string.ascii_letters + string.digits, k=5))
    #code=format((num+1), '05')
    #nam_t=taxid+'.'+str(code)
    print  >>translate, '>%s\n%s'%(code, seq)
    print >>translate_table, '%s\t%s'%(name, code)
    
translate_table.close()    
translate.close()
translate_alg=SeqGroup(transform_fasta)

translate_alg.write(format="phylip", outfile=out_file)