import re
import sys
import os
from ete3 import PhyloTree

###sys.argv[1] is path to folder where gene folders are located.
###sys.argv[2] is file with list of genes to be done.
###sys.argv[3] is output file for list of genes already split and ready to be pruned.


def main():
	with open(sys.argv[2],"r") as f:
		gene = f.readline().rstrip()
	file_name = "{}/{}/{}".format(sys.argv[1], gene, gene)
	tree_file_name = "{}.3.fa.tre".format(file_name)
	t = PhyloTree(tree_file_name)
	choice = user_choice(t, gene)
	if choice[0] == "y":
		yes_choice(tree_file_name, gene)
	else:
		no_choice(gene)
	
	
def yes_choice(tree_file_name, gene):
	t=PhyloTree(tree_file_name)
	R = t.get_midpoint_outgroup()
	t.set_outgroup(R)
	gene_names = t.get_leaf_names()
	group_list = clade_to_tree(t)
	###tree1
	cut_list = [i for i in gene_names if i not in group_list]
	tree=PhyloTree(tree_file_name)
	tree.prune(cut_list, preserve_branch_length=True)
	if gene[-3] == "_":
		n = int(gene[-1])+1
		gene1 = gene[:-1]+str(n)
	else:
		gene1 = gene+"_10"
	new_file = "{}/{}/{}.3.fa.tre".format(sys.argv[1], gene1, gene1)
	directory = ("{}/{}".format(sys.argv[1], gene1))
	os.system("mkdir {}".format(directory))
	tree.write(format=1, outfile=new_file)
	with open(sys.argv[2], "a") as todo:
		todo.write("\n{}".format(gene1))
	###tree2
	cut_list2 = [i for i in gene_names if i not in cut_list]
	tree2=PhyloTree(tree_file_name)
	tree2.prune(cut_list2, preserve_branch_length=True)
	if gene1[-3] == "_":
		n = int(gene1[-1])+1
		gene2 = gene1[:-1]+str(n)
	else:
		gene2 = gene1+"_10"
	directory2 = ("{}/{}".format(sys.argv[1], gene2))
	os.system("mkdir {}".format(directory2))
	new_file2 = "{}/{}/{}.3.fa.tre".format(sys.argv[1], gene2, gene2)
	tree2.write(format=1, outfile=new_file2)
	with open(sys.argv[2], "a") as todo:
		todo.write("\n{}".format(gene2))

def no_choice(gene):
	with open(sys.argv[3], "a") as master:
		master.write("\n{}".format(gene))
	todo_list=[line.rstrip() for line in open(sys.argv[2], "r")]
	todo_new=[i for i in todo_list if i != gene]
	with open(sys.argv[2], "w") as todo:
		for i in todo_new:
			todo.write("{}\n".format(i))
	
def user_choice(t, gene):
	print("\nShowing the gene tree for {}.".format(gene))
	view_rooted_tree(t)
	choice = raw_input("\nWould you like to split this gene family into multiple families? (y/n)")
	return(choice)	
	
######Viewing rooted tree in png file#####################################
def view_rooted_tree(t):
	try:
		R=t.get_midpoint_outgroup()
		t.set_outgroup(R)
		t.render("treeimage.png")
		os.system("open treeimage.png")
	except:
		print("\nUnable to root tree. Showing unrooted tree.")
		clade_tree.render("treeimage.png")
		os.system("open treeimage.png")
		
######Making clade list for large trees###################################
def clade_to_tree(tree):
	c="y"
	while c == "y":
		leaf1 = raw_input("\nEnter the gene on one end of the group.")
		leaf2 = raw_input("\nEnter the gene on the other end of the group.")
		try:
			mrca = tree.get_common_ancestor(leaf1, leaf2)
			group_list = mrca.get_leaf_names()
			c = "n"
		except ValueError:
			print("\nGene name not found. Try again.")
			c = "y"
	return(group_list)
	
main()