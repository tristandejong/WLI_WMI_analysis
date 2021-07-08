
import itertools
import operator
import sys
import os

chromlist = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','X','Y']


testList = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

# method to consistently get the left and right samples without mistakes
def is_good_variant(gt_list, ad_list, gq_list, quality_cutoff=10, get_left = operator.itemgetter(0, 1, 2), get_right = operator.itemgetter(3, 4, 5)):
	# if there are not at least 10 reads of all 6 samples, TRASH
	sample_read_sums = list(map(sum, ad_list))
	total_read_count = sum(sample_read_sums)
	
	
	# if all 3 samples within strain agree, but are called differently between the strains
	unique_left  = set(get_left( gt_list))
	unique_right = set(get_right(gt_list))
	if unique_left != unique_right and len(unique_left) == 1 == len(unique_right):
		testList[4] +=1
		return True
	
    # If less than 10 counts
	if total_read_count < 10:
		testList[0] +=1
		return False
		
    # if higher than 300 then throw away.
	if max(sample_read_sums) > 300: 
		testList[1] +=1
		return False
	
	# if the calls agree within methods across strains
	if get_left(gt_list) == get_right(gt_list):
		testList[2] +=1
		return False
	
	# if all genotypes are the same in all 6 samples, trash
	unique_gt = set(gt_list)
	if len(unique_gt) == 1:
		testList[3] +=1
		return False
	
	# if 5 out of 6 genotypes are all the same with one ./.   , trash #### might have to be re-evaluated later
	gt_ordered_gt_list, gt_ordered_ad_list, gt_ordered_gq_list = zip(*sorted(zip(gt_list, ad_list, gq_list), key=operator.itemgetter(0)))
	
	if gt_ordered_gt_list[0] == "./." and gt_ordered_gt_list[1] == gt_ordered_gt_list[5]:
		testList[5] +=1
		return False
		
	# make a slightly filtered version of the gt_list
	#new_gt_list = []
	#for idx in range(6):
	#	if gq_list[idx] < quality_cutoff:
	#		gt = "./."
	#	else:
	#		gt = gt_list[idx]
	#	new_gt_list.append(gt)
	#	testList[5] +=1
	
	# if 5 out of 6 genotypes are all the same with one ./.   , trash #### might have to be re-evaluated later
	#new_row_sorted = sorted(new_gt_list)
	#if new_row_sorted[0] == "./." and new_row_sorted[1] == new_row_sorted[5]:
	#	testList[6] +=1
	#	return False
	
	#if get_left(new_gt_list) == get_right(new_gt_list):
	#	testList[7] +=1
	#	return False
	
	# find the indexes of the highest quality calls per strain
	best_gq_idx_left  = max(range(0,3), key=lambda idx: gq_list[idx])
	best_gq_idx_right = max(range(3,6), key=lambda idx: gq_list[idx])	
	if gq_list[best_gq_idx_left] >= 30 and gq_list[best_gq_idx_right] >= 30:
		if gt_list[best_gq_idx_left] == gt_list[best_gq_idx_right]:
			testList[6] +=1
			return False
	
	if gt_ordered_gt_list[0] == gt_ordered_gt_list[4]:
		if gt_ordered_gq_list[5] < quality_cutoff:
			testList[7] +=1
			return False
	elif gt_ordered_gt_list[1] == gt_ordered_gt_list[5]:
		if gt_ordered_gq_list[0] < quality_cutoff:
			testList[8] +=1
			return False
	
	summed_depth_left  = [ sum(allelic_depth) for allelic_depth in itertools.izip_longest(*get_left( ad_list), fillvalue=0) ]
	summed_depth_right = [ sum(allelic_depth) for allelic_depth in itertools.izip_longest(*get_right(ad_list), fillvalue=0) ]
	merged_summed_depth = [ sum(allelic_depth) for allelic_depth in itertools.izip_longest(summed_depth_left, summed_depth_right, fillvalue=0) ]
	
	sorted_summed_depth_left  = sorted(summed_depth_left)
	sorted_summed_depth_right = sorted(summed_depth_right)
	sorted_merged_summed_depths = sorted(merged_summed_depth)
	
	
	#ad_ordered_gt_list, ad_ordered_ad_list, ad_ordered_gq_list = zip(*sorted(zip(ad_list, gt_list, gq_list), key=operator.itemgetter(0)))
	fraction_cutoff = sum(sorted_merged_summed_depths[:-1]) * 10                               # Calculate the depth of the lowest 5 combined * 10
	if sorted_merged_summed_depths[-1]  > fraction_cutoff:                                     # If the highest call is larger than the lowest 5
		per_sample_fraction_cutoff = sorted_merged_summed_depths[-1] / 5.0                     # calculate the depth of the highest / 5
		if sorted_summed_depth_left[-1]   > per_sample_fraction_cutoff and sorted_summed_depth_right[-1] > per_sample_fraction_cutoff and summed_depth_left.index(sorted_summed_depth_left[-1]) == summed_depth_right.index(sorted_summed_depth_right[-1]) and summed_depth_left.index(sorted_summed_depth_left[-1]) == merged_summed_depth.index(sorted_merged_summed_depths[-1]):	# If the highest depth of left and right are greater than the cutoff
			testList[9] +=1

			return False
		

	testList[10] +=1	
	return True

# Here it is decided what the final call will be.	
def make_strain_call(gt_list, ad_list, gq_list, quality_cutoff=10, get_left = operator.itemgetter(0, 1, 2), get_right = operator.itemgetter(3, 4, 5)):

	summed_depths_per_sample = list(map(sum, ad_list))
	
	gt_ordered_gt_list, gt_ordered_ad_list, gt_ordered_gq_list = zip(*sorted(zip(gt_list, ad_list, gq_list), key=operator.itemgetter(0)))
	
	#gt_ordered_gt_list, gt_ordered_ad_list, gt_ordered_gq_list = zip(*sorted(zip(gt_list, ad_list, gq_list, summed_depths_per_sample), key= operator.itemgetter(3)))
	gq_ordered_gt_list_left,  gq_ordered_ad_list_left,  gq_ordered_gq_list_left,  summed_depths_per_sample_left  = zip(*sorted(zip(get_left( gt_list), get_left( ad_list), get_left( gq_list), get_left( summed_depths_per_sample)), key= operator.itemgetter(2)))
	gq_ordered_gt_list_right, gq_ordered_ad_list_right, gq_ordered_gq_list_right, summed_depths_per_sample_right = zip(*sorted(zip(get_right(gt_list), get_right(ad_list), get_right(gq_list), get_right(summed_depths_per_sample)), key= operator.itemgetter(2)))

	#gt_ordered_gt_list, gt_ordered_ad_list, gt_ordered_gq_list = zip(*sorted(zip(gt_list, ad_list, gq_list, summed_depths_per_sample), key= operator.itemgetter(3)))
	ad_ordered_gt_list_left,  ad_ordered_ad_list_left,  ad_ordered_gq_list_left,  summed_depths_per_sample_left  = zip(*sorted(zip(get_left( gt_list), get_left( ad_list), get_left( gq_list), get_left( summed_depths_per_sample)), key= operator.itemgetter(3)))
	ad_ordered_gt_list_right, ad_ordered_ad_list_right, ad_ordered_gq_list_right, summed_depths_per_sample_right = zip(*sorted(zip(get_right(gt_list), get_right(ad_list), get_right(gq_list), get_right(summed_depths_per_sample)), key= operator.itemgetter(3)))
	
	# # find the indexes of the highest quality calls per strain
	# best_gq_idx_left  = max(range(0,3), key=lambda idx: gq_list[idx])
	# best_gq_idx_right = max(range(3,6), key=lambda idx: gq_list[idx])

	summed_depth_left   = [ sum(allelic_depth) for allelic_depth in itertools.izip_longest(*get_left( ad_list), fillvalue=0) ]
	summed_depth_right  = [ sum(allelic_depth) for allelic_depth in itertools.izip_longest(*get_right(ad_list), fillvalue=0) ]
	merged_summed_depth = [ sum(allelic_depth) for allelic_depth in itertools.izip_longest(summed_depth_left, summed_depth_right, fillvalue=0) ]
	
	call_left  = gq_ordered_gt_list_left[ -1]
	call_right = gq_ordered_gt_list_right[-1]
	merged_call_gq = call_left + "_" + call_right
	
	call_left  = ad_ordered_gt_list_left[ -1]
	call_right = ad_ordered_gt_list_right[-1]
	merged_call_ad = call_left + "_" + call_right
	
	merged_call = [ merged_call_gq, merged_call_ad ]
	
	call_left  = gq_ordered_gt_list_left[ -1]
	call_right = gq_ordered_gt_list_right[-1]
	
	#if force_call_right:
	#	call_right = force_call_right
	#else:
		#call_right = custom_ordered_gt_list_right[-1]
	
	# Decide if the IonProton is Het, If another is above 30, pick that one instead (for left)
	if gt_list[0] == "0/1" and max(gq_list[1:2]) >= 30:
		index_max_left = max(range(len(gq_list[1:2])), key=gq_list[1:2].__getitem__)
		call_left = gt_list[index_max_left+1]
	
	# Decide if the IonProton is Het, If another is above 30, pick that one instead (for right)
	if gt_list[3] == "0/1" and max(gq_list[4:5]) >= 30:
		index_max_right = max(range(len(gq_list[4:5])), key=gq_list[4:5].__getitem__)
		call_right = gt_list[index_max_right+4]	
		
	merged_call_custom = call_left + "_" + call_right
	merged_call.append(merged_call_custom)
	
	return merged_call
	
for chromNumber in chromlist:
	location = 'data/deepvariant_wmiwli_chr'+chromNumber+'_6_samples.gvcf' #The original gvcf

	print("running", location)

	vcf_file = open(location, 'r')

	for line in vcf_file:
			if not line.startswith("##"):
					##print line
					header = line.split()
					break
					
	names = header[9:]
	out_header = "\t" + "\t".join(names) + "\n"
	out_header2 = "\t" + "\t".join(names) +  "\t" + "call_qual" +  "\t" + "call_depth" +  "\t" + "call_custom" + "\tRNC_" + "\tRNC_".join(names) + "\n"
	#print(out_header2)
	
	#1/0
	gt_out_file = open(location+".gt", "w")
	gt_out_file.write(out_header2)
	ad_out_file = open(location+".ad", "w")
	ad_out_file.write(out_header)
	gq_out_file = open(location+".gq", "w")
	gq_out_file.write(out_header)
	
	depth_ref_file  = open(location+".ref", "w")
	depth_ref_file.write(out_header)
	depth_alt_file  = open(location+".alt", "w")
	depth_alt_file.write(out_header)
	
	complex_out_file = open(location+".complex_out.txt", "w")
	complex_out_file.write("This is the file with complex calls that were removed",)
	
	#rnc_out_file = open(location+".rnc", "w")
	#rnc_out_file.write(out_header)
	
	id_format = "{chromo}_{pos}_{ref}_{alt}".format

	round_dict = {0.0 : "0", 1.0 : "1"}

	for line_nr, line in enumerate(vcf_file, start=0):
			#if line_nr > 300:
			#	break

			cols = line.rstrip().split()

			alt_alleles = cols[4].split(",")
			allele_ids = cols[2]
			
			line_format = cols[8].split(":")
			genotype_index = line_format.index("GT") # index of the genotypes in the individual
			allelic_depth_index = line_format.index("AD") # index of the allelic depth estimation in the individual
			quality_index = line_format.index("GQ") # index of the allelic depth estimation in the individual
			RNC_index = line_format.index("RNC") 
			indiv_data = cols[9:]

			gt_list = []
			ad_list = []
			gq_list = []
			rnc_list = []
			ref_depth_list = []
			alt_depth_list = []
					
			ad_str_list = []
			
			for indiv in indiv_data:
					#print(indiv)
					indiv_info = indiv.split(":")
					
					# Genotype calls
					indiv_gt= indiv_info[genotype_index] # .replace(".", "0") was in between
					gt_list.append(indiv_gt)
					
					# Quality calls
					indiv_gq = indiv_info[quality_index].replace(".", "0")
					#print(indiv_gq)
					gq_list.append(int(indiv_gq))
					
					# Reason Not Called (RNC) 
					indiv_RNC = indiv_info[RNC_index]
					rnc_list.append(indiv_RNC)
					
					# allelic depths
					indiv_ad = list(map(int, indiv_info[allelic_depth_index].replace(".", "0").split(",")))
					
					indiv_ref = indiv_ad[0]
					indiv_alt = sum(indiv_ad[1:len(indiv_ad)])
					
					ad_list.append(indiv_ad)
					ad_str_list.append(str(sum(indiv_ad)))
					
					ref_depth_list.append(str(indiv_ref))
					alt_depth_list.append(str(indiv_alt))
					
					#print(ad_list)
					#print(ad_str_list)
			
			
			
			if is_good_variant(gt_list, ad_list, gq_list) == True:
				strain_calls = make_strain_call(gt_list, ad_list, gq_list)
				gt_out_file.write(allele_ids+"\t")
				ad_out_file.write(allele_ids+"\t")
				gq_out_file.write(allele_ids+"\t")

				depth_ref_file.write(allele_ids+"\t")
				depth_alt_file.write(allele_ids+"\t")
			
				gt_out_file.writelines("\t".join(gt_list))
				ad_out_file.writelines("\t".join(ad_str_list))
				gq_out_file.writelines("\t".join(map(str, gq_list)))
			
				depth_ref_file.writelines("\t".join(ref_depth_list))
				depth_alt_file.writelines("\t".join(alt_depth_list))
			
				gt_out_file.write("\t" + "\t".join(strain_calls))# + "\n")
				gt_out_file.write("\t" + "\t".join(rnc_list) + "\n")
			
			
				gt_out_file.write("\n")
				ad_out_file.write("\n")
				gq_out_file.write("\n")
				
				depth_ref_file.write("\n")
				depth_alt_file.write("\n")
				
			else:
				#sys.stderr.write(str(line_nr)+"\n")
				pass
			
	print("all done")

	print(testList)
	
	gt_out_file.close()
	ad_out_file.close()
	gq_out_file.close()
	depth_ref_file.close()
	depth_alt_file.close()
	
	
	
