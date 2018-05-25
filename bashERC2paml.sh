
##1 prune tree and align
python ~/ERCpaper/pythonScripts/fastaToPrunedFasta.py $1.fa $1.txt | tee $1.pr.fa | mafft --auto - > $1.pr.fa.aln
##2 trim tree
trimal -in $1.pr.fa.aln -out $1.pr.fa.trm -gt 0.8
##3 fill in blanks for missing species and make phylip
python ~/ERCpaper/pythonScripts/fastaToPhylipWithBlanks.py $1.pr.fa.trm species_list_all.txt > $1.pr.phy
##4 make .ctl aaml file for each gene
python ~/ERCpaper/pythonScripts/makeAamlCtlForGenes.py $1 ~/ERCpaper/componentFiles/ERCaaml.ctl > $1.ctl
##5 run paml
codeml $1.ctl
##6 get branch lengths from paml output
python ~/ERCpaper/pythonScripts/pamlToBranchLengthList.py $1.paml > $1_branchlengths.txt

##This takes a fasta file and a list of chosen sequences, prepares it for paml, and runs paml.
##Prunes unwanted sequences, aligns, trims, adds blank alignments for missing taxa, and makes phylip file.
##Takes the .phy file, runs it through paml, then parses the .paml file to make a list of branch lengths

##Execute from correct group folder
##So command line should read: parallel -j1 bash bashScripts/bashERC2paml.sh {}/{} :::: ~/ERCpaper/componentFiles/gene_list.txt