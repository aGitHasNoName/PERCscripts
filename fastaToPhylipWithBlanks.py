import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#########################################
##adds blank sequences to alignment for missing taxa, output is relaxed phylip format
##argv[1] is fasta file of aligned sequences
##argv[2] is list of all species
#########################################
def add_blanks():
	trimmed_file=open(sys.argv[1],'r')
	species_list=[line.rstrip() for line in open(sys.argv[2])]
	trimmed_dict=SeqIO.to_dict(SeqIO.parse(trimmed_file,"fasta"))
	sp1=species_list[0]
	l=len(trimmed_dict[sp1].seq)
	blank=str("-"*l)
	l2=len(species_list)
	sys.stdout.write("  "+str(l2)+"  "+str(l)+"\n")
	for species in species_list:
		if species not in trimmed_dict:
			trimmed_dict[species]=SeqRecord(Seq(blank),id=species)
	
	sorted_species_list = sorted(species_list)
	for record in sorted_species_list:
		sys.stdout.write("{}  {}\n".format(trimmed_dict[record].id, trimmed_dict[record].seq))

	
	trimmed_file.close()
#########################################
add_blanks()

