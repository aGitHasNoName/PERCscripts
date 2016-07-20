import plotly.plotly as py
from plotly.tools import FigureFactory as FF
import plotly.graph_objs as go
import os
import json
import pandas as pd



calculate pvalue


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
def calcPearson(study_dict):
	df=pd.DataFrame.from_dict(study_dict,"columns",float)
	pearson=df.corr()
	pearson=pearson.round(3)
	#Reverse matrix:
	pearson=pearson[::-1]
	x=pearson.columns.values.tolist()
	y=pearson.index.values.tolist()
	return(pearson,x,y)

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
	pearson,x,y=calcPearson(study_dict)
	fig=makeHeatMap(pearson,x,y)
	fig['layout'].update(margin=go.Margin(l=200,t=150))
	py.iplot(fig, filename='annotated_heatmap5')




