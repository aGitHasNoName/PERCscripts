from ete3 import PhyloTree
import re
import statistics
import os
import sys
from Bio import SeqIO

###This script removes sequences that are significantly shorter than the median sequence length. 
###sys.argv[1] is path to a tree file called {}.fa.tre and alignment file called {}.fa.aln.
###It's also path to where the output tree will be made called {}.2.tre
###sys.argv[2] is user-supplied cut off, as decimal, for short sequence length (e.g. .25)


def main():
	path = sys.argv[1]
	print("For {}:".format(path))
	t = PhyloTree("{}.fa.tre".format(path))
	t2 = removeShort(t, path)
	t2.write(outfile = "{}.2.fa.tre".format(path))
	
def removeShort(t, path):
	print("Filtering out duplicate sequences that are <{} the median length for all sequences...".format(sys.argv[2]))
	with open("{}.fa.aln".format(path), "r") as f:
		aa_list=[str(record.seq) for record in SeqIO.parse(f, "fasta")]
		aaLen_list=map(lambda x: len(re.sub("\-", "", x)), aa_list)
		aaLen_median=statistics.median(aaLen_list)
	keep_list=[]
	with open("{}.fa.aln".format(path), "r") as f:	
		for record in SeqIO.parse(f, "fasta"):
			aa=re.sub("\-","", str(record.seq))
			if len(aa) > float(sys.argv[2])*aaLen_median:
				keep_list.append(record.id)
	dif=len(aa_list)-len(keep_list)
	print("{} gene copies were removed for length.".format(str(dif)))
	t.prune(keep_list, preserve_branch_length=True)
	return(t)
	
main()