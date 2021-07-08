
import glob, os, subprocess

chromosomes = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y"]

run_SNPeff = open("07_run_SNPeff.sh", "w")

for f in glob.glob("results_final_selection/*.list"):
	
	print "processing " + f
	
	outSplit = f.split("/")[1].split("_")
	outName = "_".join(outSplit[0:4])
	
	command2 = ["java -jar picard.jar MergeVcfs "]
	for chr in chromosomes:
	
		print "running chromosome " + chr

		#print("_".join(outSplit[0:4]))
		outNameChr = outName +"_"+"chr"+ chr
		
		command1 = "gatk-4.1.7.0/gatk SelectVariants -O results_final_vcf/"+outNameChr+".vcf -V data/deepvariant_wmiwli_chr"+chr+"_6_samples.gvcf -ids " + f
		
		print command1
		runCommand = subprocess.check_output(command1, shell=True)
		command2.append("I=results_final_vcf/"+outNameChr+".vcf ")

	print "merging " + f
	command2.append("O=results_final_vcf_merged/"+outName+".vcf ")
	runCommand2 = subprocess.check_output("".join(command2), shell=True)
	command2 = []
	
	jar="/lustre/haven/proj/UTHSC0013/Tristan_BN_GATK/snpEff/snpEff/snpEff.jar"
	
	command3 = "java -Xmx4g -jar "+jar+" Rnor_6.0.86 "+"results_final_vcf_merged/"+outName+".vcf"+ " > "+ "results_snpEff/" + outName + "_ann.vcf\n" + "mv snpEff_genes.txt results_snpEff/snpEff_genes_"+outName+".txt\n" +  "mv snpEff_summary.html results_snpEff/snpEff_summary_"+outName+".html\n\n" 

	run_SNPeff.write(command3)
	
	
run_SNPeff.close()	
print("All done")


	
	
	
		#
		#testline = 

