from ete3 import PhyloTree
import re
import os
import sys

###This script removes within species duplications by randomly choosing one copy. 
###sys.argv[1] is path to a tree file called {}.2.fa.tre.
###It's also path to where the output tree will be made called {}.3.tre

def main():
	path = sys.argv[1]
	t = PhyloTree("{}.2.fa.tre".format(path))
	species_set = make_species_set(t)
	copy_list,t2 = remove_species_dups(t)
	if len(copy_list) == len(species_set):
		print("{} is SINGLE".format(path))
		t2.write(outfile = "{}.3.fa.tre".format(path))
	else:
		diff = len(copy_list)-len(species_set)
		print("For {} there are {} remaining duplicate species".format(path, diff))
		t2.write(outfile = "{}.3.fa.tre".format(path))
	
def make_species_set(t):
	leaves = []
	for leaf in t:
		leaves.append(leaf)
	l = [str(i) for i in leaves]
	l = [i.lstrip("\n--") for i in l]
	l = [re.sub("\d", "", i) for i in l]
	set_l = set(l)
	return(set_l)

def remove_species_dups(t):
	t.set_species_naming_function(lambda node: re.sub("\d", "", node.name))
	t2 = t.collapse_lineage_specific_expansions() 
	leaves2 = []
	for leaf in t2:
		leaves2.append(leaf)
	l2 = [str(i) for i in leaves2]
	l2 = [i.lstrip("\n--") for i in l2]
	return(l2, t2)

main()