#!/usr/bin/perl

#SAM coming in STDIN must be sorted be read id!!!!

use POSIX;
use List::Util qw[min max];

%read_id_hash =();

%line_hash = ();

$first_read=1;
$first_write=1;

while ($line = <>){
	chomp($line);
	##header line start with an @ in sam format, print those out
	if($line =~ m/^\@/){
		print $line."\n";
	}else{
		@line_split = split("\t",$line);
		$read_id = $line_split[0];
		if($first_read==1){
			$current_read_id = $read_id;
			$first_read=0;
		}
		##keep adding line to %current_hash until read_id does not equal current_read
		if($read_id eq $current_read_id){
			push(@line_array,$line);
		}else{
			$num_aligns = scalar @line_array;
			for $print_line (@line_array){
				$out_print_line=$print_line."\t"."TA:i:".$num_aligns;
				chomp($out_print_line);
				if($first_write==1){
					print $out_print_line;
					$first_write=0;
				}else{
					print "\n".$out_print_line;
				}
			}

			#reset %current_hash for new read_id
			undef @line_array;
			$current_read_id = $read_id;
			push(@line_array,$line);
		}
	}
}
