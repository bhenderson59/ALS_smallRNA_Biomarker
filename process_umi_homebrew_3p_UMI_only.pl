#!/usr/bin/perl

#with the v3 5' ad and v5 3' ad only the last 4 bases of the read are N's and will be added to the fastq header to be used as a UMI

$line_count=0;

@fq_record_line = ();

while($fastq_line = <STDIN>){
	chomp($fastq_line);
	
	if($line_count<4){
		push(@fq_record_line,$fastq_line);
	}
	
	$line_count++;
	
	if($line_count==4){
		#UMI is first 4 bases of read and last 4
		#read line is second line (array index 1)
		$read_length = length($fq_record_line[1]);
		$umi_right = substr($fq_record_line[1],($read_length-4),4);
		$main_read_trimmed = substr($fq_record_line[1],0,($read_length-4));
		$qualities_trimmed = substr($fq_record_line[3],0,($read_length-4));
		
		@first_line_fq_split = split(" ",$fq_record_line[0]);
		$first_line_fq_UMI_added = $first_line_fq_split[0].":UMI:".$umi_right." ".$first_line_fq_split[1];
		
		#print out new fq record
		print $first_line_fq_UMI_added."\n".$main_read_trimmed."\n".$fq_record_line[2]."\n".$qualities_trimmed."\n";
		
		$line_count = 0;
		@fq_record_line = ();
	}
}

