



def getImpact(location, label):

	vcf_file = open(location, 'r') 
	for line in vcf_file:
		if not line.startswith("##"):
			#print line
			header = line.split()
			break			
	print line

	for line_nr, line in enumerate(vcf_file, start=0):		
			
		parts = line.split()
		
		chromosome = parts[0]
		pos = parts[1]
		REF = parts[3]
		ALT = parts[4]
		info = parts[7].split("|")
		#print info
		
		
		location = info[1].replace("_"," ").replace("&"," & ")
		impact = info[2]
		gene = info[3]
		gene_ID = info[4]
		impact2 = info[5]
		geneType = info[7]
		
		
		if impact =="HIGH" or impact == "MODERATE":
			outList = [label, chromosome, pos, REF, ALT, gene, gene_ID, impact, location, impact2, geneType]
			print impact
			print info
			#1/0
			outline = "\t".join(outList)+ "\n"
			print outline
			
			outfile.write(outline)
			
			
	vcf_file.close()

# Set op the outFile
outfile = open("writing_draft/table3_highImpact.csv","w")

## Get the high impact SNP from WMI
outfile.write("WMI high impact SNPs\n")
outfile.write("\t".join(["Strain","Chromosome","position","REF","ALT","Gene name","Ensembl ID","Impact","Modification","type","genetype"])+"\n")
location1 = "results_final_Snpeff/ALT_WMI_calls_30_ann.vcf"
getImpact(location1,"WMI")
outfile.write("\n")

## Get the high impact SNP from WLI
outfile.write("WLI high impact SNPs\n")
outfile.write("\t".join(["Strain","Chromosome","position","REF","ALT","Gene name","Ensembl ID","Impact","Modification","type","genetype"])+"\n")
location2 = "results_final_Snpeff/ALT_WLI_calls_30_ann.vcf"
getImpact(location2,"WLI")





outfile.close()
