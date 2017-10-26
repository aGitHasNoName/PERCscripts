from ete3 import PhyloTree
import re
import os

###This script removes within species duplications. 
###sys.argv[1] is path to a tree file called {}.fa.tre.
###It's also path to where the output tree will be made called {}.2.tre

def main():
	with open("{}.fa.tre".format(sys.arg[1]), "r") as f:
	copy_list=remove_species_dups(sys.argv[1])
	species_set=make_species_set(sys.argv[1])
	if len(copy_list)==len(species_set):
		print("{} is SINGLE".format(i))
	else:
		diff=len(copy_list)-len(species_set)
		print("For {} there are {} remaining duplicate species".format(i,diff))
		species_copy_list=[re.sub("\d","",j) for j in copy_list]
		for x in species_set:
			if species_copy_list.count(x)>1:
				print("Extra copy is in {}".format(x))

def make_species_set(path):
	t = PhyloTree("{}.fa.tre".format(path))
	leaves=[]
	for leaf in t:
		leaves.append(leaf)
	l=[str(i) for i in leaves]
	l=[i.lstrip("\n--") for i in l]
	l=[re.sub("\d","",i) for i in l]
	set_l=set(l)
	
	return(set_l)

def remove_species_dups(path):
	t = PhyloTree("{}.fa.tre".format(path))
	t.set_species_naming_function(lambda node: re.sub("\d","",node.name))

	t2 = t.collapse_lineage_specific_expansions() 
	leaves2=[]
	for leaf in t2:
		leaves2.append(leaf)
	l2=[str(i) for i in leaves2]
	l2=[i.lstrip("\n--") for i in l2]

	return(l2)

