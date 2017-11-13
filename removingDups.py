from ete3 import PhyloTree
import re
import statistics
import os
import sys
from Bio import SeqIO

###This script removes within species duplications. 
###sys.argv[1] is path to a tree file called {}.fa.tre and alignment file called {}.fa.aln.
###It's also path to where the output tree will be made called {}.2.tre
###sys.argv[2] is cut off, as decimal, for short sequence length.

def main():
	path=sys.argv[1]
	t = PhyloTree("{}.fa.tre".format(path))
	species_set=make_species_set(t)
	t2=removeShort(t,path)
	copy_list,t3=remove_species_dups(t2)
	if len(copy_list)==len(species_set):
		print("{} is SINGLE".format(path))
		with open("{}_single.txt".format(path),"w") as f:
			for i in copy_list:
				f.write("{}\n".format(i))
	else:
		diff=len(copy_list)-len(species_set)
		print("For {} there are {} remaining duplicate species".format(path,diff))
		t3.write(outfile="{}.2.fa.tre".format(path))

def removeShort(t,path):
	print("Filtering out duplicate sequences that are <{} the average length for all sequences...".format(sys.argv[2]))
	with open("{}.fa.aln".format(path), "r") as f:
		aa_list=[str(record.seq) for record in SeqIO.parse(f, "fasta")]
		aaLen_list=map(lambda x: len(re.sub("\-","",x)),aa_list)
		aaLen_median=statistics.median(aaLen_list)
	keep_list=[]
	with open("{}.fa.aln".format(path), "r") as f:	
		for record in SeqIO.parse(f, "fasta"):
			aa=re.sub("\-","",str(record.seq))
			if len(aa) > float(sys.argv[2])*aaLen_median:
				keep_list.append(record.id)
	t.prune(keep_list,preserve_branch_length=True)
	return(t)
	
def make_species_set(t):
	leaves=[]
	for leaf in t:
		leaves.append(leaf)
	l=[str(i) for i in leaves]
	l=[i.lstrip("\n--") for i in l]
	l=[re.sub("\d","",i) for i in l]
	set_l=set(l)
	return(set_l)

def remove_species_dups(t):
	t.set_species_naming_function(lambda node: re.sub("\d","",node.name))
	t2 = t.collapse_lineage_specific_expansions() 
	leaves2=[]
	for leaf in t2:
		leaves2.append(leaf)
	l2=[str(i) for i in leaves2]
	l2=[i.lstrip("\n--") for i in l2]
	return(l2,t2)

main()