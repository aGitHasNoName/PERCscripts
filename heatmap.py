import plotly.plotly as py
from plotly.tools import FigureFactory as FF
import plotly.graph_objs as go
import os
import json
import pandas as pd


fix dictionary - fabs getting added to seedfree?
calculate pvalues
MSH->MLH linear like yeast (and maybe other well known meiosis genes?)
quantify differences between All and clades for single copy
scikit-bio package sklearn bicluster in python???
functional sorting from Tair and look for patterns
TALK
branch length outliers
dup A vs. dup B in a clade that is otherwise single


############################################
## Calculate ERC: calculates Pearson’s r, calculates p value, makes heatmap, displays on plotly website if using command line OR in notebook if using jupyter notebook.
##takes gene list and Clade Number:
##1 Single copy genes for all species
##2 All species (For single copy genes, this is same as single)
##3 Grass clade
##4 Brassicaceae clade
##5 Fabidae clade
##6 Seedfree clade (non-monophyletic)
############################################


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

#df=df[cols]


def calcPvalue():
	

	
#works
def makeHeatMap(pearson,x,y):
	#Get matrix from pandas style:
	z=pearson.values.tolist()
	#Create heatmap:
	fig = FF.create_annotated_heatmap(z, x=x, y=y, colorscale=([0, 'rgb(249,236,235)'], [1, 'rgb(236,58,38)']), font_colors=(['rgb(0,0,0)', 'rgb(255,255,255)']))
	return(fig)

def main(gene_list,cladeNum):
	study_dict=makeStudyDict(gene_list, cladeNum)
	a,x=statsMatrix(study_dict)
	fig=makeHeatMap(pearson,x,y)
	fig['layout'].update(margin=go.Margin(l=200,t=150))
	py.iplot(fig, filename='annotated_heatmap5')




