from ete3 import PhyloTree
import re
import os

def main():
	with open(sys.arg[1], "r") as f:
		orthogroup_list=[line.rstrip() for line in f.readlines()]
	for i in orthogroup_list:
		copy_list=remove_species_dups(i)
		species_set=make_species_set(i)
		if len(copy_list)==len(species_set):
			print("{} is SINGLE".format(i))
		else:
			diff=len(copy_list)-len(species_set)
			print("For {} there are {} remaining duplicate species".format(i,diff))


def make_species_set(gene):
	t = PhyloTree("/Users/colbyjones/Documents/singleCopy/genes/{}/{}.fa.tre".format(gene,gene))
	leaves=[]
	for leaf in t:
		leaves.append(leaf)
	l=[str(i) for i in leaves]
	l=[i.lstrip("\n--") for i in l]
	l=[re.sub("\d","",i) for i in l]
	set_l=set(l)
	
	return(set_l)


def remove_species_dups(gene):
	t = PhyloTree("/Users/colbyjones/Documents/singleCopy/genes/{}/{}.fa.tre".format(gene,gene))
	t.set_species_naming_function(lambda node: re.sub("\d","",node.name))

	t2 = t.collapse_lineage_specific_expansions() 
	leaves2=[]
	for leaf in t2:
		leaves2.append(leaf)
	l2=[str(i) for i in leaves2]
	l2=[i.lstrip("\n--") for i in l2]

	return(l2)

