import re
import sys
import os
from ete3 import PhyloTree

###sys.argv[1] is path to folder where gene folders are located.
###sys.argv[2] is file with list of genes to be done.
###sys.argv[3] is output file for list of genes already split and ready to be pruned.


def main():
	next_choice = "y"
	while next_choice[0] == "y":
		with open(sys.argv[2],"r") as f:
			gene_list = [line.rstrip() for line in f]
		if len(gene_list) == 0:
			next_choice = "n"
			print("List complete")
		else:
			gene = gene_list[0]
			print("Working on gene {}".format(gene))
			file_name = "{}/{}/{}".format(sys.argv[1], gene, gene)
			tree_file_name = "{}.3.fa.tre".format(file_name)
			t = PhyloTree(tree_file_name)
			algae_choice, split_choice = user_choice(t, gene)
			if split_choice[0] == "y":
				yes_choice(tree_file_name, gene, algae_choice)
			else:
				no_choice(gene)
	
	
def yes_choice(tree_file_name, gene, algae_choice):
	t=PhyloTree(tree_file_name)
	R = t.get_midpoint_outgroup()
	t.set_outgroup(R)
	gene_names = t.get_leaf_names()
	if algae_choice[0] == "y":
		print("\nFirst, let's define the algae clade.")
		algae_list = clade_to_tree(t)
	else:
		algae_list = []
	outlier_choice = raw_input("\nIs there another group that is sister to multiple families? (y/n)")
	if outlier_choice[0] == "y":
		print("\nLet's define the outlier clade.")
		outlier_list = clade_to_tree(t)
	else:
		outlier_list=[]
	print("Select one family to define. \nYou will have a later chance to split this family again if needed.")
	group_list = clade_to_tree(t)
	###tree1
	cut_list = [i for i in gene_names if i not in group_list]
	cut_list = cut_list + algae_list + outlier_list
	gene1 = yesMake(cut_list, gene, tree_file_name)
	###tree2
	cut_list1 = [i for i in gene_names if i not in cut_list]
	cut_list1 = cut_list1 + algae_list + outlier_list
	gene2 = yesMake(cut_list1, gene1, tree_file_name)
	with open(sys.argv[2], "r") as f:
		todo_list=[line.rstrip() for line in f]
	todo_list=[i for i in todo_list if i != gene]
	todo_list.append(gene1)
	todo_list.append(gene2)
	with open(sys.argv[2], "w") as todo:
		for i in todo_list:
			todo.write(i+"\n")

def yesMake(cut_list, gene, tree_file_name):
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
	with open(sys.argv[3], "a") as master:
		master.write("{}\n".format(gene1))
	return(gene1)

def no_choice(gene):
	with open(sys.argv[3], "a") as master:
		master.write("{}\n".format(gene))
	todo_list=[line.rstrip() for line in open(sys.argv[2], "r")]
	todo_new=[i for i in todo_list if i != gene]
	with open(sys.argv[2], "w") as todo:
		for i in todo_new:
			todo.write("{}\n".format(i))
	
def user_choice(t, gene):
	print("\nShowing the gene tree for {}.".format(gene))
	view_rooted_tree(t)
	algae_choice = raw_input("\nIs there an algae group that is sister to multiple families? (y/n)")
	choice = raw_input("\nWould you like to split this gene family into multiple families? (y/n)")
	return(algae_choice, choice)	
	
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