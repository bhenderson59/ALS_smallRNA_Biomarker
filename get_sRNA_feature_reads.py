#!/gpfs/gpfs2/software/python-3.7.6/bin/python

import re
import pysam
import sys

#first input is the BAM
inbam_file = sys.argv[1]

#second input is the tsv of feature locations
feature_tsv = sys.argv[2]

#third input is the out tsv of counts
out_counts_tsv = sys.argv[3]

#open BAM
inBAM=pysam.AlignmentFile(inbam_file,"rb")

#open feature tsv 
#feature tsv has a header so save it for later
header_line=1

read_cnts_dict = {}

all_ref_names = inBAM.references

with open(feature_tsv) as features_list:
	for feature in features_list:
		feature = feature.strip()
		if header_line == 1:
			header_split = feature.split("\t")
			header_text = "\t".join(header_split[5:10])
			header_line = 0
		else:
			feature_split = feature.split("\t")
			query_ref_name = eval(feature_split[0])
			query_start = int(feature_split[1]) -1 #pysam wants everything 0 indexed
			query_ends = int(feature_split[2]) - 1
			
			feature_name = str(feature_split[5])
			feature_id = str(feature_split[6])
			rna_type = str(feature_split[7])
			annot_name = str(feature_split[8])
			annot_id = str(feature_split[9])
			
			#concat the feature info with a tab as a key
			feature_key = "\t".join([feature_name,feature_id,rna_type,annot_name,annot_id])
			
			#now count reads that match this query
			
			#check to see if the ref exists in the BAM
			if query_ref_name in all_ref_names:
				match_align_num_reads = inBAM.count(contig=query_ref_name,start=query_start,stop=query_ends)
			
			else:
				match_align_num_reads = 0
			
			#check if this feature has already been added to read_cnts_dict
			if feature_key in read_cnts_dict:
				#add the match_align_num_reads to the value already in the dictionary
				old_num_reads = read_cnts_dict[feature_key]
				read_cnts_dict[feature_key] = old_num_reads + match_align_num_reads
				
			else:
				read_cnts_dict[feature_key] = match_align_num_reads

#close input BAM
inBAM.close()

#now print out hash

with open(out_counts_tsv,"w") as out_file:
	out_file.write(header_text+"\tcounts\n")
	for feature_annots,read_cnts in read_cnts_dict.items():
		out_file.write(feature_annots+"\t"+str(read_cnts)+"\n")


quit()

			
			
			

