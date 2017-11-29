import re
import sys
import itertools
import os
from ete3 import PhyloTree
########################################################
#This script takes a fasta file and walks you through the pruning of the gene.

#Products are individual lists of genes for each duplication in grasses, brassicaceae,
#fabideae, and the seedfree plants, as well as full lists containing each possible
#combination of groups.
#Also produces a master list file containing each file name.

#sys.argv[1] is a file that includes all genes.
#sys.argv[2] is a .txt file that includes all species in column 0 and clade name
#in column 1, separated by a tab. If species is not in a clade, write noclade in column 1.
########################################################

######Checks for Python3################################
if (sys.version_info > (3,0)):
	print("Python 2 required!\n")
	sys.exit(1)

######MAIN FUNCTIONS###############################################################

def main():
	genes_todo=[line.rstrip() for line in open(sys.argv[1])]
	speciesList,cladeDict=makeSpeciesList()
	for gene in genes_todo:
		c=raw_input("\nPrune gene {}? (y/n)".format(gene))
		if c[0]=="y":
			prune_main(gene,speciesList,cladeDict)
			gtd_copy=[line.rstrip() for line in open(sys.argv[1])]
			with open(sys.argv[1], "w") as f:
				for i in gtd_copy:
					if i != gene:
						f.write(i+"\n")
			with open("genes_done.txt", "a") as g:
				g.write(gene+"\n")
			genes_todo=gtd_copy
		else:
			print("\nExiting. Program will return to this gene next time.")
			break

######Makes species list and dictionary of clades and the species in each clade######
def makeSpeciesList():
	speciesList=[]
	cladeDict={}
	with open(sys.argv[2],"r") as f:
		for line in f.readlines():
			species=line.split("\t")[0].rstrip()
			clade=line.split("\t")[1].rstrip()
			speciesList.append(species)
			if clade not in cladeDict.keys():
				cladeDict[clade]=[]
			cladeDict[clade].append(species)
	return(speciesList,cladeDict)

######Runs all functions################################################
def prune_main(gene,speciesList,cladeDict):
	gene=str(gene)
	erase_previous_files(gene)
	copy_list=copies_in_group(gene)
	gene_type=count_summarize(gene,copy_list,speciesList,cladeDict)
	choice="n"
	if gene_type=="small":
		small_family(gene)
	elif gene_type=="single":
		single_copy(gene,copy_list,cladeDict)
	else:
		print ("\nShowing the gene tree.")
		clade_tree=PhyloTree(gene+"/"+gene+".3.fa.tre")
		view_rooted_tree(clade_tree)
		choice2=raw_input("\nWould you like to split this gene family into multiple families? (y/n)")
		if choice2[0]=="y":
			pre_prune(gene)
		else:
			choice=raw_input("\nContinue with pruning as single gene family? (y/n)")
	if choice[0]=="y":
		make_clade_groups(gene,cladeDict,copy_list,speciesList)
		make_all_lists(gene,cladeDict)

######MAKE LIST OF COPIES IN ORTHOGROUP########################################
def copies_in_group(gene):
	with open("{}/{}.3.fa.tre".format(gene,gene),"r") as f:
		f=f.read()
		list1=[i for i in f.split(",")]
		list2=[i.split(":")[0] for i in list1]
		copy_list=[i.strip("()") for i in list2]
		return(copy_list)

######GENE SUMMARY#######################################################
def count_summarize(gene,copy_list,species_list,cladeDict):
	l=[re.sub("\d", "", i) for i in copy_list]
	count_dict={species:(l.count(species)) for species in species_list}	
	slist=[key for key in count_dict if count_dict[key] > 0]
	n_species=len(slist)
	clist=[value for value in count_dict.values()]
	n_copies=sum(clist)
	species_total=len(species_list)
	print("For gene {}:".format(gene))
	print("number of unique gene copies in the orthogroup: {}".format(n_copies))
	print("number of species represented (out of {} possible): {}".format(species_total,n_species))
	
	for key,value in cladeDict.items():
		clade_count={k:v for k,v in count_dict.items() if k in value}
		print("\nthe number of copies per species in clade {}:".format(key))
		print(clade_count)

	if n_species < (.6*species_total):
		gene_type="small"
	elif n_species == n_copies:
		gene_type="single"
	else:
		gene_type="large"

	return (gene_type)

########SMALL GENE FAMILIES###############################################
def small_family(gene):
	with open("gene_progress_list.txt", "a") as p:
		p.write("{}\tSMALL\tNOTDONE\n".format(gene))
	print ("\nGene family is too small. Gene has been added to gene_progress_list.txt as SMALL NOTDONE")

########SINGLE COPY GENES#################################################
def single_copy(gene,copy_list,cladeDict):
	for key,value in cladeDict.items():
		clade_name=key+"_1"
		group_list=[i for i in copy_list if re.sub("\d","",i) in value]
		if len(group_list)>0:
			saving_group(gene, group_list, clade_name)
	######Saving whole tree######
	clade_list=["all_1","single_1"]
	for clade in clade_list:
		with open(gene+"/"+gene+"_"+clade+"_prune.txt", "a") as allfile:
			for copy in copy_list:
				allfile.write(copy+"\n")
		with open(gene+"/"+gene+"_master_tree_list.txt", "a") as master:
			text="{}_{}\n".format(gene,clade)
			master.write(text)
	with open("gene_progress_list.txt", "a") as p:
		p.write("{}\tSINGLE\tDONE\n".format(gene))
	print ("\nGene family is single copy. No pruning required. All clades are saved.\nGene has been added to gene_progress_list.txt as SINGLE DONE")

########SPLITTING LARGE GENES FAMILIES####################################
def pre_prune(gene):
	full_tree=PhyloTree(gene+"/"+gene+".3.fa.tre")
	gene_names=full_tree.get_leaf_names()
	m=100
	start_gene="{}_all{}".format(gene,str(m))
	os.system("mkdir {}".format(start_gene))
	full_tree.write(format=1, outfile="{}/{}.3.fa.tre".format(start_gene,start_gene))
	m=m+1
	l=[start_gene]
	for item in l:
		full_tree=PhyloTree("{}/{}.3.fa.tre".format(item,item))
		view_rooted_tree(full_tree)
		print("Tree for {}".format(item))
		algae_choice = raw_input("\nIs there an algae group that is sister to all shown families? (y/n)")
		outlier_choice = raw_input("\nIs there another monophyletic or non-monophyletic outlier group that is sister to all families shown? (y/n)")
		c=raw_input("Split off a monophyletic clade? (y/n)")
		while c[0]=="y":
			if algae_choice[0] == "y":
				print("\nFirst, let's define the algae clade.")
				algae_list = clade_to_tree(clade_tree)
			else:
				algae_list = []
			if outlier_choice[0] == "y":
				print("\nLet's define the outlier group. \nFirst you can add any spceices that are in a monophyletic clade")
				outlier_list = clade_to_tree(clade_tree)
				other_copies = raw_input("If there are other genes in the outlier group, enter them here, separated by a space, or else enter n.")
				if other_copies != "n":
					other_list = other_copies.split(" ")
					outlier_list = outlier_list + other_list
			else:
				outlier_list=[]
			b="{}_all{}".format(gene, str(m))
			l.append(b)
			tree1=PhyloTree("{}/{}.3.fa.tre".format(item,item))
			R=tree1.get_midpoint_outgroup()
			tree1.set_outgroup(R)
			group_list=clade_to_tree(tree1)
			group_list=group_list + algae_list + outlier_list
			gene_names=tree1.get_leaf_names()
			if len(group_list)==len(gene_names):
				c1=raw_input("\nList includes all copies on tree.\nMake gene with all copies? (y/n)")
				if c1=="y":
					c="n"
				else:
					print("\nGroup crosses root. Unable to make group.\nChoose new group.")
					c="y"
			else:
				cut_list=[i for i in gene_names if i not in group_list]
				cut_list = cut_list + algae_list + outlier_list
				os.system("mkdir {}".format(b))
				tree2=PhyloTree("{}/{}.3.fa.tre".format(item,item))
				R=tree2.get_midpoint_outgroup()
				tree2.set_outgroup(R)
				tree2.prune(group_list,preserve_branch_length=True)
				tree2.write(format=1, outfile="{}/{}.3.fa.tre".format(b,b))
				tree1.prune(cut_list,preserve_branch_length=True)
				tree1.write(format=1, outfile="{}/{}.3.fa.tre".format(item,item))
				m=m+1
				print ("\nTree now looks like this.")
				view_rooted_tree(tree1)
				c=raw_input("Split off a monophyletic clade? (y/n)")
				if c[0] == "y":
					algae_choice = raw_input("\nIs there an algae group that is sister to all shown families? (y/n)")
					outlier_choice = raw_input("\nIs there another monophyletic or non-monophyletic outlier group that is sister to all families shown? (y/n)")

	with open("genes_todo.txt", "a") as p:
		for i in l:
			p.write(i+"\n")
	
########CLADES###########################################################
def make_clade_groups(gene,cladeDict,copy_list,species_list):
	######Getting all grass gene copies######
	for key,value in cladeDict.items():
		species_keep=[i for i in copy_list if re.sub("\d","",i) in value]
		######Checking if the list is empty######
		if key=="noclade":
			make_other_groups(gene, species_keep, species_list)
		else:
			if len(species_keep) == 0:
				print ("\nThere are no genes in clade {} for {}. We will continue with the next clade.".format(key,gene))
			else:
				######Removing stray within-clade gene copies from the clade######
				cut_list=cut_stray_genes(gene, species_keep, value)
				######Designating whole clade duplications######
				define_groups(gene, cut_list, value, species_keep, key)


########OTHERS##############################################################
def make_other_groups(gene, species_keep, species_list):
	full_tree=PhyloTree(gene+"/"+gene+".3.fa.tre")
	######Checking if the list is empty######
	if len(species_keep) == 0:
		print ("\nThere are no other genes in this gene family.")
	else:
		######Removing stray within-clade gene copies from the clade######
		cut_list=cut_stray_other(gene, species_keep, species_list)
		######Making it a group######
		group_list=cut_list
		######Checking that there is only one gene per species######
		check_set={str(item[0:3]) for item in group_list}
		while len(group_list) != len(check_set):
			view_counts(cut_list, species_list)
			group_str=raw_input("\nYou can only have one gene per species. Enter more genes to cut, separated by a space: ")
			group_list=[item for item in group_list if item not in group_str]
			check_set={str(item[0:3]) for item in group_list}
		print("\nThere are "+str(len(group_list))+" genes in this group.\nGroup looks like:")
		print (group_list)
		print ("\nMaking the group.")
		######Saving gene group as a file######
		with open(gene+"/"+gene+"_noclade_prune.txt", "a") as group_file:
				for i in group_list:
					group_file.write(i+"\n")
		######Saving name of group to a master list######
		with open(gene+"/"+gene+"_master_tree_list.txt", "a") as master:
			master.write(gene+"_noclade\n")

############ALL#########################################################
def make_all_lists(gene, cladeDict):
	master_list=[line.rstrip() for line in open(gene+"/"+gene+"_master_tree_list.txt")]
	clades=[]
	for key in cladeDict:
		clade_list=[i for i in master_list if key in i]
		if len(clade_list)>0:
			clades.append(clade_list)
	n=1
	for i in itertools.product(*clades):
		filenames=[]
		for x in range(len(i)):
			filenames.append(gene+"/"+i[x]+"_prune.txt")
		with open(gene+"/"+gene+"_all_"+str(n)+".txt", "w") as allfile:
			for fname in filenames:
				with open(fname) as infile:
					allfile.write(infile.read())
		with open(gene+"/"+gene+"_master_tree_list.txt", "a") as master:
			master.write(gene+"_all_"+str(n)+"\n")
		n=n+1	
	if n==2:
		os.system("cp {}/{}_all_1.txt {}/{}_single.txt".format(gene,gene,gene,gene))
	with open("gene_progress_list.txt", "a") as p:
		p.write("{}\tNORMAL\tDONE\n".format(gene))
	print("\nComplete.\nGene has been added to gene_progress_list.txt as NORMAL DONE.")

######SUBFUNCTIONS###################################################################

######Remove stray genes from clade tree##############################
def cut_stray_genes(gene, species_keep, species_list):
	######Showing the tree######
	clade_tree=PhyloTree(gene+"/"+gene+".3.fa.tre")
	clade_tree.prune(species_keep,preserve_branch_length=True)
	if len(species_keep)>1:
		view_rooted_tree(clade_tree)
		print("\nThis is the clade tree. There are "+str(len(species_keep))+" total gene copies.\n")
	else:
		print("\nSpecies tree only contains 1 species. Tree will not be shown.")
	cut_list=species_keep
	view_counts(cut_list, species_list)
	######Removing stray within-clade gene copies from the clade######
	cut_question=raw_input("\nAre there stray genes to cut? (y/n)")
	while cut_question[0]== "y":
		cut_gene_str=raw_input("\nEnter genes to cut, separated by a space: ")
		cut_gene_list=[item for item in cut_gene_str.split()]
		cut_list=[i for i in cut_list if i not in cut_gene_list]
		if set(cut_gene_list).issubset(species_keep):
			try:
				clade_tree.prune(cut_list,preserve_branch_length=True)
				view_rooted_tree(clade_tree)
				view_counts(cut_list, species_list)
			except ValueError:
				print ("\nSomething is wrong with the way the genes were entered. You entered:\n"+cut_gene_str+"\nCut abandoned.")
		else:
			print ("\nAt least one gene is not found on the tree. You entered:\n"+cut_gene_str+"\nCut abandoned.")
		cut_question=raw_input("\nAre there more genes to cut? (y/n)")
	return (cut_list)

######Remove stray genes from other tree##############################
def cut_stray_other(gene, species_keep, species_list):
	######Showing the tree######
	clade_tree=PhyloTree(gene+"/"+gene+".3.fa.tre")
	clade_tree.prune(species_keep,preserve_branch_length=True)
	if len(species_keep)>1:
		view_rooted_tree(clade_tree)
		print("\nThis is the clade tree. There are "+str(len(species_keep))+" total gene copies.\n")
	else:
		print("\nSpecies tree only contains 1 species. Tree will not be shown.")
	cut_list=species_keep
	view_counts(cut_list, species_list)
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
				view_counts(cut_list, species_list)
			except ValueError:
				print ("\nSomething is wrong with the way the genes were entered. You entered:\n"+cut_gene_str+"\nCut abandoned.")
		else:
			print ("\nAt least one gene is not found on the tree. You entered:\n"+cut_gene_str+"\nCut abandoned.")
		cut_question=raw_input("\nAre there more genes to cut? (y/n)")
	return (cut_list)

######Making groups#####################################################
def define_groups(gene, cut_list, species_list, species_keep, clade_name):
	clade_tree=PhyloTree(gene+"/"+gene+".3.fa.tre")
	clade_tree.prune(cut_list,preserve_branch_length=True)
	n=1
	######Designating whole clade duplications######
	choice=raw_input("\nWould you like to make a group? (y/n)")
	while choice[0] == "y":
		choice4=raw_input("\nIf this group includes all the genes left on the tree, type a.\nIf this group is a monophyletic clade, type c.\nOtherwise, type n.")
		if choice4[0]=="a":
			group_list=cut_list
		elif choice4[0]=="c":
			group_list=choose_clade(clade_tree)
		else:
			group_str=raw_input("\nEnter genes for the group, separated by a space: ")
			group_list=[item for item in group_str.split()]
		######Checking that there is only one gene per species######
		group_list2,subclade_name=check_single_group(group_list)
		######Checking for typos######
		for i in group_list2:
			if set(i).issubset(species_keep):
				######Allow a chance to back out, for example if user forgot to enter spaces######
				print("\nThere are "+str(len(i))+" genes in this group.\nGroup looks like:")
				print (i)
				choice3=raw_input("\nMake the group? (y/n)")
			else:	
				print ("\nAt least one gene is not found on the tree. You entered:\n")
				print (i)
				choice3=raw_input("\nEnter n to abandon this list and start again.")
			if choice3[0]=="y":
				######Saving group as file and add group to master list######
				if subclade_name=="ynyn":
					clade_filename="{}_{}".format(clade_name, n)
					saving_group(gene, i, clade_filename)
				else:
					clade_filename="{}_{}_{}".format(clade_name, subclade_name, n)
					saving_group(gene, i, clade_filename)
				n=n+1
			else: 
				print("\nGroup abandoned.")
				choice=raw_input("\nWould you like to make a group for this clade? (y/n)")
			cut_list=[j for j in cut_list if j not in i]
		######Checking to see if the tree is empty######
		if len(cut_list) == 0:
			print("\nThe tree is now empty. We will continue with the next clade.")
			choice="n"
			######Preparing for next group######
		else:
			choice2=raw_input("\nWould you like to view the tree with the group removed? (y/n)")
			if choice2[0] == "y":
				clade_tree.prune(cut_list,preserve_branch_length=True)
				view_rooted_tree(clade_tree)
				view_counts(cut_list, species_list)
				choice=raw_input("\nWould you like to make another group for this clade? (y/n)")
			else: 
				print("\nGroup abandoned.")
				choice=raw_input("\nWould you like to make a group for this clade? (y/n)")
			
######Checking that there is only one gene per species####################			
def check_single_group(group_list):		
	group_list2=[]
	subclade_name="ynyn"
	check_set={str(re.sub("\d","",i)) for i in group_list}
	if len(group_list) != len(check_set):
		print (group_list)
		choice=raw_input("\nYou can only have one gene per species. If there is a subclade duplication that you would like to designate, enter y. If you would like to retype the gene names, enter n:" )
		if choice[0]=="y":
			subclade_name=raw_input("\nEnter a name for the subclade:")
			subclade_count=raw_input("\nHow many copies of the subclade are in the clade?")
			if int(subclade_count)==1:
				pass
			else:
				for i in range(int(subclade_count)):
					subclade_str=raw_input("\nDuplication number {}. \nEnter the gene names for this duplication, separated by a space:".format(str(i+1)))
					subclade_list=[j for j in subclade_str.split()]
					subclade_species=[re.sub("\d","",m) for m in subclade_list]
					for g in group_list:
						if re.sub("\d","",g) not in subclade_species:
							subclade_list.append(g)
					group_list2.append(subclade_list)
		else:
			group_str=raw_input("Enter genes for the group, separated by a space: ")
			group_list=[i for i in group_str.split()]
			group_list2.append(group_list)		
	else:
		group_list2.append(group_list)
	return (group_list2,subclade_name)
		
######Saving group as file and add group to master list###################
def saving_group(gene, group_list, clade_name):
	######Saving gene group as a file######
	with open(gene+"/"+gene+"_"+clade_name+"_prune.txt", "a") as group_file:
		for i in group_list:
			group_file.write(i+"\n")
	######Saving name of group to a master list######
	with open(gene+"/"+gene+"_master_tree_list.txt", "a") as master:
		master.write(gene+"_"+clade_name+"\n")

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
	
######Viewing the number of gene copies for each species###################
def view_counts(cut_list, species_list):
	l=[re.sub("\d","",i) for i in cut_list]
	count_dict={species:(l.count(species)) for species in species_list}
	print ("\nNumber of gene copies per species:")
	print (count_dict)
	
######Making clade list for large trees###################################
def clade_to_tree(tree):
	c="y"
	while c=="y":
		leaf1=raw_input("\nEnter the gene on one end of the group.")
		leaf2=raw_input("\nEnter the gene on the other end of the group.")
		try:
			mrca = tree.get_common_ancestor(leaf1,leaf2)
			group_list=mrca.get_leaf_names()
			c="n"
		except ValueError:
			print("\nGene name not found. Try again.")
			c="y"
	return (group_list)
	
######Erasing previous files##############################################
def erase_previous_files(gene):
	os.system("rm {}/{}_*".format(gene,gene))

############RUN THE PROGRAM##########################################################

if __name__ == '__main__': main()



