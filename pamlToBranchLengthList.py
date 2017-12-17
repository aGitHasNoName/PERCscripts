import sys
##############################################
##parses paml file, makes a file with two lines, one of branch numbers and one of branch lengths
##argv[1] is paml filename
##############################################
def branchlengthlist():
	file=open(sys.argv[1],"r")
	lines=file.readlines()
	for i, line in enumerate(lines):
		if "lnL(ntime" in line:
			branches= lines[i+1]
			lengths=lines[i+2]
	sys.stdout.write(branches)
	sys.stdout.write(lengths)
	file.close()
##############################################
branchlengthlist()