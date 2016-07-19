
import sys
import plotly.plotly as py
from plotly.tools import FigureFactory as FF
import plotly.graph_objs as go
import os
import json
import pandas as pd


Identify genes/branch lengths of interest
DONEcalculate pearson
calculate pvalue
DONEheatmap in plotly

############################################
##Lets you choose genes to compare ERC for, calculates Pearson’s r, displays heatmap
##argv[1]is 
############################################


def chooseGenes():
	with open("~/Users/colbyjones/PhD/PERCscripts/fam_dict.txt","r") as t:
		fam_dict=json.load(t)
		input_start=input("\nSearch for genes?\nPress y to begin.")
		while input_start=="y":
			input_search=input("\nSearch for genes by which category?\nEnter the correct number.\n\n0 Arabidopsis locus\n1 Gene name\n2 Number of genes in family\n3 Source (where gene was listed as meiosis)\n4 Type (single or large)\n5 Number of pruned genes within family\n6 GO terms\n7 Process\n8 Linked genes\n9 Produce 2N gametes?\n10 Involved in apomixis?\n11 I know the family ID numbers I would like to use.\n12 Single copy genes\n13 All genes\n")
			if input_search==0:
			if input_search==1:
			if input_search==2:
			if input_search==3:
			if input_search==4:
			if input_search==5:
			if input_search==6:
			if input_search==7:
			if input_search==8:
			if input_search==9:
			if input_search==10:
			if input_search==11:
			if input_search==12:
			if input_search==13:
		
		cladeNum=input("\nWhich clade would you like to compare?\nEnter the correct number:\n1 Single copy genes for all species\n2 All species\n3 Grass clade\n4Brassicaceae clade\n5Fabidae clade\n6 Seedfree clade (non-monophyletic)\n")
		
	return(gene_list, cladeNum)




##works
def makeStudyDict(gene_list, cladeNum):
	with open("/Users/colbyjones/PhD/PERCscripts/gene_dict.txt","r") as t:
		gene_dict=json.load(t)
	study_dict={}
	for gene in gene_list:
		v=gene_dict[gene]
		clade=v[cladeNum]
		for i in clade:
			i_list=[t.strip() for t in i.split()]
			study_dict[i_list[0]]=i_list[1:]
	return(study_dict)



#works
def calcPearson(study_dict):
	df=pd.DataFrame.from_dict(study_dict,"columns",float)
	pearson=df.corr()
	pearson=pearson.round(3)
	#Reverse matrix:
	pearson=pearson[::-1]
	x=pearson.columns.values.tolist()
	y=pearson.index.values.tolist()
	return(pearson,x,y)
	
#works
def makeHeatMap(pearson,x,y):
	#Get matrix from pandas style:
	z=pearson.values.tolist()
	#Create heatmap:
	fig = FF.create_annotated_heatmap(z, x=x, y=y, colorscale=([0, 'rgb(249,236,235)'], [1, 'rgb(236,58,38)']), font_colors=(['rgb(0,0,0)', 'rgb(255,255,255)']))
	return(fig)

	
main()
	study_dict=makeStudyDict(gene_list, cladeNum)
	pearson,x,y=calcPearson(study_dict)
	fig=makeHeatMap(pearson,x,y)
	fig['layout'].update(margin=go.Margin(l=200,t=150))
	py.iplot(fig, filename='annotated_heatmap3')




	

def calculatePearson():
	genes=open(sys.argv[1], "r")
	genelist = [x.strip() for x in genes.readlines()]
	length_dict={}
	
	samplegene=genelist[0]
	samplefile=open(samplegene+"_branchlengths.txt", "r")
	samplel=samplefile.readlines()
	s=samplel[0]
	branches=[t.strip() for t in s.split()]
	samplefile.close()
	
	for item in genelist:
		pamlfile=open(item+"_branchlengths.txt", "r")	
		lines=pamlfile.readlines()
		l=lines[1]
		length=[t.strip() for t in l.split()]
		length_dict[item]=length
		pamlfile.close()
	df=pd.DataFrame.from_dict(length_dict,'columns',float)
	df.index=branches
	pearson=df.corr()
	#print (pearson)
	sys.stdout.write(pearson)
############################################
calculatePearson()



