import sys
#####################################
##Creates aaml .ctl file for each gene with correct input and output file names
##argv[1] is gene name
##argv[2] is original aaml file
#####################################
def makectl():
	example1=open(sys.argv[2],"r")
	example=example1.readlines()
	for line in example:
		new=line.replace("uniquegenename",sys.argv[1])
		sys.stdout.write(new)
	example1.close()
#####################################
makectl()