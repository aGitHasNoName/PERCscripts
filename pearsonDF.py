#import plotly.plotly as py
#from plotly.tools import FigureFactory as FF
#import plotly.graph_objs as go
import sys
import os
import json
import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr 


#fix dictionary - fabs getting added to seedfree?
#MSH->MLH linear like yeast (and maybe other well known meiosis genes?)
#quantify differences between All and clades for single copy
#scikit-bio package sklearn bicluster in python???
#functional sorting from Tair and look for patterns
#TALK
#branch length outliers
#dup A vs. dup B in a clade that is otherwise single


############################################
## Calculate ERC: calculates Pearson’s r, calculates p value, makes heatmap, displays on plotly website if using command line OR in notebook if using jupyter notebook.
##Takes output file name, gene list .txt file, clade number, and optional sort list .txt file
##If no sorting required, enter ‘none’ as 4th argument
#(choose clade number from below):
##1 Single copy genes for all species
##2 All species (For single copy genes, this is same as single)
##3 Grass clade
##4 Brassicaceae clade
##5 Fabidae clade
##6 Seedfree clade (non-monophyletic)
############################################

def makeStudyDict(gene_list_file, cladeNum):
	with open("/Users/colbyjones/PhD/PERCscripts/gene_dict.txt","r") as t:
		gene_dict=json.load(t)
	study_dict={}
	gene_list=[line.strip("\n") for line in open(gene_list_file)]
	for gene in gene_list:
		v=gene_dict[gene]
		clade=v[int(cladeNum)]
		for i in clade:
			i_list=[t.strip() for t in i.split()]
			study_dict[i_list[0]]=i_list[1:]
	return(study_dict)

def statsMatrix(study_dict):
	x=[]
	results=[]
	pvalues=[]
	done=[]
	study_dict2=study_dict
	for key in study_dict:
		x.append(key)
		for i in study_dict2:
			check=set()
			check.add(key)
			check.add(i)
			if check not in done:
				v=study_dict[key]
				v=np.array(v).astype(np.float)
				v2=study_dict2[i]
				v2=np.array(v2).astype(np.float)
				try:
					tup=pearsonr(v,v2)
					results.append(tup[0])
					if tup[1]>0:
						pvalues.append(tup[1])
					done.append(check)
				except ValueError:
					print(len(v))
					print(len(v2))
	l=len(x)
	a=np.arange(l*l).reshape(l,l).astype(np.object)
	iu=np.triu_indices(l)
	il=np.tril_indices(l,-1)
	results2=[float(round(i,3)) for i in results]
	pvalues2=[]
	for i in pvalues:
		if i<0.00005:
			pvalues2.append("<0.00005")
		elif i<0.0005:
			pvalues2.append("<0.0005")
		elif i<0.005:
			pvalues2.append("<0.005")
		elif i<0.05:
			pvalues2.append("<0.05")
		else:
			pvalues2.append(">0.05")
	a[iu]=results2
	a[il]=pvalues2
	return(a,x)

def saveDataFrame(a,x,outname,sort_list_file):
	df=pd.DataFrame(data=a,index=x,columns=x)
	if sort_list_file!="none":
		sort_list=[line.strip("\n") for line in open(sort_list_file)]
		cols=sort_list
		df=df[cols]
	df.to_csv("/Users/colbyjones/Desktop/MeiosisGenesTAIR/{}.csv".format(outname))


def main(outname,gene_list_file,cladeNum,sort_list_file):
	study_dict=makeStudyDict(gene_list_file, cladeNum)
	a,x=statsMatrix(study_dict)
	saveDataFrame(a,x,outname,sort_list_file)
	
main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
	

##python /Users/colbyjones/PhD/PERCscripts/pearsonDF.py SingleAll /Users/colbyjones/Desktop/MeiosisGenesTAIR/SortedSingleSingleGenes.txt 1 none



#Code for plotly:	
#	#Get matrix from pandas style:
#	z=df.values.tolist()
#	#Create heatmap:
#	fig = FF.create_annotated_heatmap(z, x=x, y=y, colorscale=([0, 'rgb(249,236,235)'], [1, 'rgb(236,58,38)']), font_colors=(['rgb(0,0,0)', 'rgb(255,255,255)']))
#	return(fig)	
#	
#	fig=makeHeatMap(pearson,x,y)
#	fig['layout'].update(margin=go.Margin(l=200,t=150))
#	py.iplot(fig, filename='annotated_heatmap5')




