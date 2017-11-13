import re
import sys
import itertools
import os
import ete3
from ete3 import PhyloTree

if (sys.version_info > (3,0)):
	print("Python 2 required!\n")
	sys.exit(1)

###sys.argv[1] is path to gene folders without trailing /
###sys.argv[2] is gene list as .txt file


def main():
	path=sys.argv[1]
	gene_list=[line.rstrip() for line in open(sys.argv[2],"r")]
	for gene in gene_list:
		try:
			newpath="{}/{}/{}".format(path,gene,gene)
			print("Pruning {}...".format(gene))
			copy_list,species_list=make_species_list(newpath)
			new_copy_list=cut_stray_other(newpath,copy_list,species_list)
			with open(newpath+"_single.txt","w") as f:
				for i in new_copy_list:
					f.write(i+"\n")
		except ete3.parser.newick.NewickError:
			print("No tree for "+newpath)

def make_species_list(path):
	t = PhyloTree("{}.2.fa.tre".format(path))
	leaves=[]
	for leaf in t:
		leaves.append(leaf)
	l=[str(i) for i in leaves]
	l=[i.lstrip("\n--") for i in l]
	l2=[re.sub("\d","",i) for i in l]
	return(l,l2)
	
######Viewing rooted tree in png file#####################################
def view_rooted_tree(clade_tree):
	try:
		R=clade_tree.get_midpoint_outgroup()
		clade_tree.set_outgroup(R)
		clade_tree.render("treeimage.png")
		os.system("open treeimage.png")
	except:
		print("\nUnable to root tree. Showing unrooted tree.")
		clade_tree.render("treeimage.png")
		os.system("open treeimage.png")

######Remove kalanchoe and switchgrass duplicates##########################
def remove_some_dups(count_dict,species_keep):
	unwanteds=["tetraploidKalanchoe","switchgrass","diploidKalanchoe"]
	unwanteds=[i for i in unwanteds if i in count_dict.keys()]
	new_keep=[]
	for j in species_keep:
		if re.sub("\d","",j) not in unwanteds:
			new_keep.append(j)
	for i in unwanteds:
		n=(count_dict[i])
		i_keep=[]
		for j in species_keep:
			if re.sub("\d","",j)==i:
				i_keep.append(j)
		for k in range(1,2):
			new_keep.append(i_keep[k])
	return(new_keep)

######Viewing the number of gene copies for each species###################
def view_counts(cut_list, species_list):
	l=[re.sub("\d","",i) for i in cut_list]
	count_dict={species:(l.count(species)) for species in species_list if l.count(species)>1}
	return (count_dict)

######Making a group from one monophyletic clade###########################
def choose_clade(clade_tree):
	c="y"
	while c=="y":
		leaf1=raw_input("\nEnter the gene on one end of the group.")
		leaf2=raw_input("\nEnter the gene on the other end of the group.")
		try:
			mrca = clade_tree.get_common_ancestor(leaf1,leaf2)
			group_list=mrca.get_leaf_names()
			c="n"
		except ValueError:
			print("\nGene name not found. Try again.")
			c="y"
	return(group_list)
	
######Remove stray genes from other tree##############################
def cut_stray_other(gene, species_keep2, species_list):
	######Showing the tree######
	count_dict=view_counts(species_keep2, species_list)
	species_keep=remove_some_dups(count_dict,species_keep2)
	print ("\nNumber of extra gene copies per species:")
	count_dict=view_counts(species_keep, species_list)
	print(count_dict)
	clade_tree=PhyloTree(gene+".2.fa.tre")
	clade_tree.prune(species_keep,preserve_branch_length=True)
	if len(species_keep)>1:
		view_rooted_tree(clade_tree)
		print("\nThis is the clade tree. There are "+str(len(species_keep))+" total gene copies.\n")
	else:
		print("\nSpecies tree only contains 1 species. Tree will not be shown.")
	cut_list=species_keep
	######Removing stray within-clade gene copies from the clade######
	cut_question=raw_input("\nAre there stray genes to cut? (y/n)")
	while cut_question[0]== "y":
		choice4=raw_input("\nIf this group is a monophyletic clade, type c.\nOtherwise, type n.")
		if choice4[0]=="c":
			cut_gene_list=choose_clade(clade_tree)
		else:
			cut_gene_str=raw_input("\nEnter genes to cut, separated by a space: ")
			cut_gene_list=[item for item in cut_gene_str.split()]
		cut_list=[i for i in cut_list if i not in cut_gene_list]
		if set(cut_gene_list).issubset(species_keep):
			try:
				clade_tree.prune(cut_list,preserve_branch_length=True)
				view_rooted_tree(clade_tree)
				print ("\nNumber of extra gene copies per species:")
				count_dict=view_counts(cut_list, species_list)
				print(count_dict)
			except ValueError:
				print ("\nSomething is wrong with the way the genes were entered. You entered:\n"+cut_gene_str+"\nCut abandoned.")
		else:
			print ("\nAt least one gene is not found on the tree. You entered:\n"+cut_gene_str+"\nCut abandoned.")
		cut_question=raw_input("\nAre there more genes to cut? (y/n)")
	return (cut_list)
	
	
main()