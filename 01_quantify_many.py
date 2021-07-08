

## Loop over each file



# Get chromosome

# Get position

# get number of REF, HET, ALT, OTHER

from operator import itemgetter 

chromlist = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','X','Y']

outfile = open("data_deepvar/0_all_deepvar_table_ALT.csv", 'w')
outfile2 = open("data_deepvar/0_all_deepvar_table_HET.csv", 'w')
outfile3 = open("data_deepvar/0_all_deepvar_table_ALT_q30.csv", 'w')

outfile.write("chr,pos,ref,new,REF,HET,ALT,NONE,OTHER,lowQ\n")
outfile2.write("chr,pos,ref,new,REF,HET,ALT,NONE,OTHER,lowQ\n")

for chromNumber in chromlist:

	print("Tabilizing ", chromNumber)
	location = 'data_deepvar/deepvariant_chr'+chromNumber+'_48_samples.gvcf.adj' #deepvariant_wmiwli_chrY_6_samples.gvcf

	table_file = open(location, 'r')
	
	first = 1
	for line in table_file:
	
		if first == 1:
			first = 0
			continue
		
		split_line = line.rstrip().split("\t")
		name = split_line[0]		
		chr, pos, ref, new = itemgetter(0, 1, 2, 3,)(name.split("_"))
		#print chr, pos, ref, new

		REF = split_line.count("0/0")
		HET = split_line.count("0/1")
		ALT =  split_line.count("1/1")
		Qlow = split_line.count("x/x")
		NONE = split_line.count("./.")
		OTHER =  len(split_line) - (REF+HET+ALT+NONE)
		
		
		if ALT > 0:
			outfile.write(",".join([str(chr), str(pos), str(ref), str(new), str(REF), str(HET), str(ALT), str(NONE), str(OTHER),str(Qlow)]) + "\n")
		if HET > 0:
			outfile2.write(",".join([str(chr), str(pos), str(ref), str(new), str(REF), str(HET), str(ALT), str(NONE), str(OTHER),str(Qlow)]) + "\n")

outfile.close()
outfile2.close()


print("all done")