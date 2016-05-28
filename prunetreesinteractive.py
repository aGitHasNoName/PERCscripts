import re
import sys
from ete3 import PhyloTree
import itertools
########################################################
#This script takes a fasta file and walks you through the pruning of the gene.

#Products are individual lists of genes for each duplication in grasses, brassicaceae,
#fabideae, and the seedfree plants, as well as full lists containing each possible
#combination of groups.
#Also produces a master list file containing each file name.

#sys.argv[1] is gene name
########################################################

######Checks for Python3################################
if (sys.version_info < (3,0)):
print("Python 3 required!\n")
sys.exit(1)


######MAIN FUNCTIONS###############################################################

######Runs all functions################################################
def prune_gene(gene):
	count_summarize(gene)
	if gene_type="single":
		single_copy(gene)
	elif gene_type="large":
		pre_prune(gene)
	else:
		choice=input("\nContinue with pruning? (y/n)")
	if choice[0]=="y":
		make_grass_groups(gene)
		make_brass_groups(gene)
		make_fab_groups(gene)
		make_seedfree_groups(gene)
		make_other_groups(gene)
		make_all_lists(gene)

######GENE SUMMARY#######################################################
def count_summarize(gene):
	gene=str(gene)
	species_list=[line.rstrip() for line in open("species_list_all.txt")]
	grass_list=["Sbi","Zma","Sit","Svi","Pvi","Pha","Osa","Bdi","Bsta"]
	brass_list=["Ath","Aly","Cru","Cgr","Bst","Bra","Esa"]
	fab_list=["Mtr","Pvu","Gma","Csa","Ppe","Mdo","Fve"]
	seedfree_list=["Smo","Ppa","Sfa","Cre","Vca","Csu","CCM","RCC","Olu"]
	groups_list=seedfree_list+fab_list+brass_list+grass_list
	with open(gene+"/"+gene+".dup.fa.tre") as file:
		f=file.read()
	f=re.sub("\d|,|:|\(|\)|\.|;", " ", f)
	l=[i for i in f.split()]
	count_dict={species:(l.count(species)) for species in species_list}	

	grass_count={key:value for key,value in count_dict.items() if key in grass_list}
	brass_count={key:value for key,value in count_dict.items() if key in brass_list}
	fab_count={key:value for key,value in count_dict.items() if key in fab_list}
	seedfree_count={key:value for key,value in count_dict.items() if key in seedfree_list}
	other_count={key:value for key,value in count_dict.items() if key not in groups_list}
	
	words="\nFor gene "+gene+", the number of copies per species:"
	wordsg="\n\nFor grasses:\n"
	wordsb="\n\nFor brassicaceae:\n"
	wordsf="\n\nFor fabidiae:\n"
	wordssf="\n\nFor seedfree:\n"
	wordso="\n\nFor all other species:\n"
	print (words, wordsg, grass_count, wordsb, brass_count, wordsf, fab_count, wordssf, seedfree_count, wordso, other_count)


	return (gene_type)

########SINGLE COPY GENES#################################################
def single_copy(gene):
	gene=str(gene)

########SPLITTING LARGE GENES FAMILIES####################################
def pre_prune(gene):
	gene=str(gene)

########GRASSES###########################################################
def make_grass_groups(gene):
	gene=str(gene)
	######Getting all grass gene copies######
	with open(gene+"/"+gene+".dup.fa.tre") as file:
		f=file.read()
	gene_names=re.findall("[A-Z][a-z][a-z][0-9][0-9][0-9]|[A-Z][a-z][a-z][a-z][0-9][0-9][0-9]|[A-Z][A-Z][A-Z][0-9][0-9][0-9]", f)	
	species_list=["Sbi","Zma","Sit","Svi","Pvi","Pha","Osa","Bdi","Bsta"]
	species_keep=[i for i in gene_names if i[0:3] in species_list or i[0:4] in species_list]
	######Checking if the list is empty######
	if len(species_keep) == 0:
		print ("\nThere are no grass genes in this gene family. We will continue with the next clade.")
	else:
		######Removing stray within-clade gene copies from the clade######
		cut_stray_genes(gene, species_keep, species_list)
		######Designating whole clade duplications######
		define_groups(gene, cut_list, species_list)
		
########BRASSES###########################################################
def make_brass_groups(gene):
	gene=str(gene)
	######Getting all brass gene copies######
	with open(gene+"/"+gene+".dup.fa.tre") as file:
		f=file.read()
	gene_names=re.findall("[A-Z][a-z][a-z][0-9][0-9][0-9]|[A-Z][a-z][a-z][a-z][0-9][0-9][0-9]|[A-Z][A-Z][A-Z][0-9][0-9][0-9]", f)	
	species_list=["Ath","Aly","Cru","Cgr","Bst","Bra","Esa"]
	species_keep=[i for i in gene_names if i[0:3] in species_list and i[0:4] != "Bsta"]
	######Checking if the list is empty######
	if len(species_keep) == 0:
		print ("\nThere are no brass genes in this gene family. We will continue with the next clade.")
	else:
		######Removing stray within-clade gene copies from the clade######
		cut_stray_genes(gene, species_keep, species_list)
		######Designating whole clade duplications######
		define_groups(gene, cut_list, species_list)

##########FABS############################################################
def make_fab_groups(gene):
	gene=str(gene)
	######Getting all fab gene copies######
	with open(gene+"/"+gene+".dup.fa.tre") as file:
		f=file.read()
	gene_names=re.findall("[A-Z][a-z][a-z][0-9][0-9][0-9]|[A-Z][a-z][a-z][a-z][0-9][0-9][0-9]|[A-Z][A-Z][A-Z][0-9][0-9][0-9]", f)
	species_list=["Mtr","Pvu","Gma","Csa","Ppe","Mdo","Fve"]
	species_keep=[i for i in gene_names if i[0:3] in species_list]
	######Checking if the list is empty######
	if len(species_keep) == 0:
		print ("\nThere are no fab genes in this gene family. We will continue with the next clade.")
	else:
		######Removing stray within-clade gene copies from the clade######
		cut_stray_genes(gene, species_keep, species_list)
		######Designating whole clade duplications######
		define_groups(gene, cut_list, species_list)
		
##########SEEDFREE#########################################################
def make_seedfree_groups(gene):
	gene=str(gene)
	######Getting all seedfree gene copies######
	with open(gene+"/"+gene+".dup.fa.tre") as file:
		f=file.read()
	gene_names=re.findall("[A-Z][a-z][a-z][0-9][0-9][0-9]|[A-Z][a-z][a-z][a-z][0-9][0-9][0-9]|[A-Z][A-Z][A-Z][0-9][0-9][0-9]", f)
	species_list=["Smo","Ppa","Sfa","Cre","Vca","Csu","CCM","RCC","Olu"]
	species_keep=[i for i in gene_names if i[0:3] in species_list]
	######Checking if the list is empty######
	if len(species_keep) == 0:
		print ("\nThere are no seedfree genes in this gene family. We will continue with the rest of the tree.")
	else:
		######Removing stray within-clade gene copies from the clade######
		cut_stray_genes(gene, species_keep, species_list)
		######Designating whole clade duplications######
		define_groups(gene, cut_list, species_list)

########OTHERS##############################################################
def make_other_groups(gene):
	gene=str(gene)
	######Getting all other gene copies######
	with open(gene+"/"+gene+".dup.fa.tre") as file:
		f=file.read()
	gene_names=re.findall("[A-Z][a-z][a-z][0-9][0-9][0-9]|[A-Z][a-z][a-z][a-z][0-9][0-9][0-9]|[A-Z][A-Z][A-Z][0-9][0-9][0-9]", f)
	species_list=['Cpa', 'Gra', 'Tca', 'Csi', 'Ccl', 'Mes', 'Rco', 'Lus', 'Ptr', 'Spu', 'Egr', 'Vvi', 'Kma', 'Stu', 'Sly', 'Mgu', 'Aco', 'Mac', 'Spo', 'Atr']
	species_keep=[i for i in gene_names if i[0:3] in species_list]
	######Checking if the list is empty######
	if len(species_keep) == 0:
		print ("\nThere are no other genes in this gene family.")
	else:
		######Showing the tree######
		clade_tree=PhyloTree(gene+"/"+gene+".dup.fa.tre")
		clade_tree.prune(species_keep,preserve_branch_length=True)
		view_rooted_tree(clade_tree)
		print("\nThis is the OTHER tree. There are "+str(len(species_keep))+" total gene copies.\n")
		t=[item[0:3] for item in species_keep]
		count_dict={species:(t.count(species)) for species in species_list}
		print (count_dict)
		######Removing stray within-clade gene copies from the clade######
		cut_question=input("\nAre there stray genes to cut? (y/n)")
		cut_list=species_keep
		while cut_question[0]== "y":
			cut_gene_str=input("Enter genes to cut, separated by a space: ")
			cut_gene_list=[item for item in cut_gene_str.split()]
			cut_list=[i for i in cut_list if i not in cut_gene_list]
			if set(cut_gene_list).issubset(species_keep):
				try:
					clade_tree.prune(cut_list,preserve_branch_length=True)
					view_rooted_tree(clade_tree)
				except ValueError:
					print ("\nSomething is wrong with the way the genes were entered. You entered:\n"+cut_gene_str+"\nCut abandoned.")
			else:
				print ("\nAt least one gene is not found on the tree. You entered:\n"+cut_gene_str+"\nCut abandoned.")
			cut_question=input("\nAre there more genes to cut? (y/n)")
		######Making it a group######
		group_list=cut_list
		######Checking that there is only one gene per species######
		check_set={str(item[0:3]) for item in group_list}
		while len(group_list) != len(check_set):
			print (group_list)
			group_str=input("\nYou can only have one gene per species. Enter genes to cut, separated by a space: ")
			group_list=[item for item in group_list if item not in group_str]
			check_set={str(item[0:3]) for item in group_list}
		print("\nThere are "+str(len(group_list))+" genes in this group.\nGroup looks like:")
		print (group_list)
		print ("\nMaking the group.")
		######Saving gene group as a file######
		with open(gene+"/"+gene+"_other_prune.txt", "a") as group_file:
				for i in group_list:
					group_file.write(i+"\n")
		######Saving name of group to a master list######
		with open(gene+"/"+gene+"_master_tree_list.txt", "a") as master:
			master.write(gene+"_other\n")

############ALL#########################################################
def make_all_lists(gene):
	gene=str(gene)
	master_list=[line.rstrip() for line in open(gene+"/"+gene+"_master_tree_list.txt")]
	g=[i for i in master_list if "grass" in i]
	b=[i for i in master_list if "brass" in i]
	f=[i for i in master_list if "fab" in i]
	s=[i for i in master_list if "seedfree" in i]
	o=[i for i in master_list if "other" in i]
	clades=[g,b,f,s,o]
	clades=[i for i in clades if len(i)>0]
	n=1
	for i in itertools.product(*clades):
		filenames=[]
		for x in range(len(i)):
			filenames.append(gene+"/"+i[x]+"_prune.txt")
		with open(gene+"/"+gene+"_"+str(n)+"_all.txt", "w") as allfile:
			for fname in filenames:
				with open(fname) as infile:
					allfile.write(infile.read())
		with open(gene+"/"+gene+"_master_tree_list.txt", "a") as master:
			master.write(gene+"_"+str(n)+"_all\n")
		n=n+1	


######SUBFUNCTIONS###################################################################

######Remove stray genes from clade tree##############################
def cut_stray_genes(gene, species_keep, species_list):
	######Showing the tree######
	clade_tree=PhyloTree(gene+"/"+gene+".dup.fa.tre")
	clade_tree.prune(species_keep,preserve_branch_length=True)
	view_rooted_tree(clade_tree)
	print("\nThis is the clade tree. There are "+str(len(species_keep))+" total gene copies.\n")
	cut_list=species_keep
	view_counts(cut_list)
	######Removing stray within-clade gene copies from the clade######
	cut_question=input("\nAre there stray genes to cut? (y/n)")
	while cut_question[0]== "y":
		cut_gene_str=input("\nEnter genes to cut, separated by a space: ")
		cut_gene_list=[item for item in cut_gene_str.split()]
		cut_list=[i for i in cut_list if i not in cut_gene_list]
		if set(cut_gene_list).issubset(species_keep):
			try:
				clade_tree.prune(cut_list,preserve_branch_length=True)
				view_rooted_tree(clade_tree)
				view_counts(cut_list)
			except ValueError:
				print ("\nSomething is wrong with the way the genes were entered. You entered:\n"+cut_gene_str+"\nCut abandoned.")
		else:
			print ("\nAt least one gene is not found on the tree. You entered:\n"+cut_gene_str+"\nCut abandoned.")
		cut_question=input("\nAre there more genes to cut? (y/n)")
	return (cut_list)
	
######Making groups#####################################################
def define_groups(gene, cut_list, species_list):
	clade_tree=PhyloTree(gene+"/"+gene+".dup.fa.tre")
	clade_tree.prune(cut_list,preserve_branch_length=True)
	n=1
	######Designating whole clade duplications######
	choice=input("\nWould you like to make a group? (y/n)")
	while choice[0] == "y":
		choice4=input("\nIf this group includes all the genes left on the tree, type a.\nIf this group is a monophyletic clade, type c.\nOtherwise, type n.")
		if choice4[0]=="a":
			group_list=cut_list
		elif choice4[0]=="c":
			choose_clade(clade_tree)
		else:
			group_str=input("\nEnter genes for the group, separated by a space: ")
			group_list=[item for item in group_str.split()]
		######Checking that there is only one gene per species######
		check_single_group(group_list)
		######Checking for typos######
		if set(group_list).issubset(species_keep):
			######Allow a chance to back out, for example if user forgot to enter spaces######
			print("\nThere are "+str(len(group_list))+" genes in this group.\nGroup looks like:")
			print (group_list)
			choice3=input("\nMake the group? (y/n)")
		else:	
			print ("\nAt least one gene is not found on the tree. You entered:\n")
			print (group_list)
			choice3=input("\nEnter n to abandon this list and start again.")
		if choice3[0]=="y":
			######Saving group as file and add group to master list######
			saving_group(gene, group_list)
			n=n+1
			cut_list=[i for i in cut_list if i not in group_list]
			######Checking to see if the tree is empty######
			if len(cut_list) == 0:
				print("\nThe tree is now empty. We will continue with the next clade.")
				choice="n"
			######Preparing for next group######
			else:
				choice2=input("\nWould you like to view the tree with the group removed? (y/n)")
				if choice2[0] == "y":
					clade_tree.prune(cut_list,preserve_branch_length=True)
					view_rooted_tree(clade_tree)
					view_counts(cut_list)
				choice=input("\nWould you like to make another group for this clade? (y/n)")
		else: 
			print("\nGroup abandoned.")
			choice=input("\nWould you like to make a group for this clade? (y/n)")
			
######Checking that there is only one gene per species####################			
def check_single_group(group_list):		
	check_set={str(item[0:3]) for item in group_list}
	while len(group_list) != len(check_set):
		print (group_list)
		group_str=input("\nYou can only have one gene per species. Enter genes for the group, separated by a space: ")
		group_list=[item for item in group_str.split()]
		check_set={str(item[0:3]) for item in group_list}
	return (group_list)
		
######Saving group as file and add group to master list###################
def saving_group(gene, group_list):
	######Saving gene group as a file######
	with open(gene+"/"+gene+"_"+str(n)+"_grass_prune.txt", "a") as group_file:
		for i in group_list:
			group_file.write(i+"\n")
	######Saving name of group to a master list######
	with open(gene+"/"+gene+"_master_tree_list.txt", "a") as master:
		master.write(gene+"_"+str(n)+"_grass\n")

######Viewing rooted tree in png file#####################################
def view_rooted_tree(clade_tree):
	R=clade_tree.get_midpoint_outgroup()
	clade_tree.set_outgroup(R)
	clade_tree.render("treeimage.png")
	os.system("open treeimage.png")
	
######Making a group from one monophyletic clade###########################
def choose_clade(clade_tree):
	leaf1=input("\nEnter the gene on one end of the group.")
	leaf2=input("\nEnter the gene on the other end of the group.")
	mrca = clade_tree.get_common_ancestor(leaf1,leaf2)
	group_list=mrca.get_leaf_names()
	return(group_list)
	
######Viewing the number of gene copies for each species###################
def view_counts(cut_list):
	l=[i[0:3] for i in cut_list]
	count_dict={species:(l.count(species[0:3])) for species in species_list}
	print ("\nNumber of gene copies per species:")
	print (count_dict)

############RUN THE PROGRAM##########################################################
gene=sys.argv[1]
prune_gene(gene)



