#!/usr/bin/python

import re
import pysam
import sys

#first argument is the max number of mismatches to allow for a valid alignment
num_mismatch = int(sys.argv[1])

#second argument is the raw bam file
infile = sys.argv[2]


#this is a dictionary of kept aligns keyed by read id
read_id_dict = {}


inBAM=pysam.AlignmentFile(infile,"rb")

for align_seg in inBAM.fetch():
	#print(str(align_seg))
	#print(align_seg.get_reference_positions())
	mismatches = align_seg.get_tag(tag="XM")
	
	if mismatches <= num_mismatch:
		align_key = align_seg.query_name
		if align_key in read_id_dict:
			#this read id has been found before. check if current read has fewer mismatches then one in dictionary
			prev_align_seg = read_id_dict[align_key]
			prev_align_num_mismatch = prev_align_seg.get_tag(tag="XM")
			
			#replace the entry in dictionary if current align has lower mismatches
			if mismatches < prev_align_num_mismatch:
				read_id_dict[align_key] = align_seg
				
			elif mismatches == prev_align_num_mismatch:
				#current align has same num mismatches; keep the highest align score
				current_align_score = align_seg.get_tag(tag="AS")
				prev_align_score = prev_align_seg.get_tag(tag="AS")
				
				if current_align_score > prev_align_score:
					#if the current align has a better align score than what is in the dictionary, replace
					read_id_dict[align_key] = align_seg
					
				elif current_align_score == prev_align_score:
					#to get here, both current and previous reads have the same mismatches and same align scores
					#to resolve these, take the one that aligns to the reference alphabetically first
					
					prev_ref_name = prev_align_seg.reference_name
					current_ref_name = align_seg.reference_name
					sorted_refs = sorted([prev_ref_name,current_ref_name])
					
					if sorted_refs[0] == current_ref_name:
						read_id_dict[align_key] = align_seg
		
		else:
			read_id_dict[align_key] = align_seg


#read_id_dict now has the picked alignment for every valid read
#print out as a bam

#write the picked aligns to STDOUT
picked_aligns_outBAM = pysam.AlignmentFile("-", "wb", template=inBAM)


for pick_align in read_id_dict.values():
	picked_aligns_outBAM.write(pick_align)

#close open files
inBAM.close()

quit()