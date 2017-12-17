import sys
from Bio import SeqIO
#########################################
##makes a new fasta file with select sequences that are listed in a .txt file
##argv[1] is original fasta file
##argv[2] is text file of sequences you want to keep
#########################################
def prune_seq():
	old_seq=open(sys.argv[1],'r')
	select_list=[line.rstrip() for line in open(sys.argv[2])]
	for record in SeqIO.parse(old_seq,"fasta"):
		for i in select_list:
			if i==record.id:
				j=str(record.id)
				new=j[0:-3]
				record.id=new
				SeqIO.write(record,sys.stdout,"fasta")
	old_seq.close()
########################################
prune_seq()