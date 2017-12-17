
##1 prune tree and align
python pythonScripts/fastaToPrunedFasta.py $1.2.fa $1_prune.txt | tee $1.pr.fa | mafft --auto - > $1.pr.fa.aln
##2 trim tree
trimal -in $1.pr.fa.aln -out $1.pr.fa.trm -gt 0.8
##3 fill in blanks for missing species and make phylip
python pythonScripts/fastaToPhylipWithBlanks.py $1.pr.fa.trm componentFiles/species_list_all.txt > $1.pr.phy
##4 make .ctl aaml file for each gene
python pythonScripts/makeAamlCtlForGenes.py $1 componentFiles/ERCaaml.ctl > $1.ctl
##5 run paml
codeml $1.ctl
##6 get branch lengths from paml output
python pythonScripts/pamlToBranchLengthList.py $1.paml > $1_branchlengths.txt

##This takes a fasta file and a list of chosen sequences, prepares it for paml, and runs paml.
##Prunes unwanted sequences, aligns, trims, adds blank alignments for missing taxa, and makes phylip file.
##Takes the .phy file, runs it through paml, then parses the .paml file to make a list of branch lengths

##Execute from ERC folder
##So command line should read: parallel -j1 bash bashScripts/bashERC2paml.sh {}/{} :::: componentFiles/gene_list.txt