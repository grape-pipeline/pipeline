#!/soft/bin/perl
# DGK

#    GRAPE
#    Copyright (C) 2009-2011 Centre for Genomic Regulation (CRG)
#
#    This file is part of GRAPE.
#
#    GRAPE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    GRAPE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with GRAPE.  If not, see <http://www.gnu.org/licenses/>.

#    Author : David Gonzalez, david.gonzalez@crg.eu

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

# This script should parse the output of the GEM splitmapper at the moment it
# returns only the unique maps

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_files_from_table_sub);
use RNAseq_pipeline_settings3 ('read_config_file','read_file_list','get_dbh');
use RNAseq_GEM3 qw(parse_gem_line);

# Get command line options
my $mismatches;
my $closest_hit=1;
my $file_list;
my $multi;

GetOptions('mismatches=i' => \$mismatches,
	   'multi' => \$multi,
	   'closest=i' => \$closest_hit);

my %options=%{read_config_file()};
$mismatches=$options{'MISMATCHES'} - 1;
$file_list=$options{'FILELIST'};

# Get some subroutines
my $dbh=get_dbh();
my $split_map_table=$options{'PREFIX'}.'_split_mapping';
*get_split_file=get_files_from_table_sub($dbh,
					 $split_map_table);
*split_parser=get_split_parser($mismatches);
my $coord_parser=get_split_coord_parser();

my %files=%{read_file_list($file_list)};

# Get for each lane the files corresponding to the junctions mapping
my %lane_files=%{get_lane_files(\%files)};

# Process the files
print STDERR "Processing files as single reads\n";
foreach my $pair (keys %lane_files) {
    foreach my $lane (keys %{$lane_files{$pair}}) {
	print STDERR "Processing\t",$lane_files{$pair}{$lane},"\n";
	split_parser($lane_files{$pair}{$lane},
		     $lane);
    }
}

exit;

sub get_coords {
    my $string=shift;
    my $pos=shift;
    my $read_length=shift;
    my %hit;

    # First we need to parse everything
    # Get the strand of the upstream and of the downstream
    my ($up_string,$down_string)=split(/~/,$string);
    my ($up_coords,$down_coords);
    ($hit{'up_chr'},$up_coords)=split(':',$up_string);
    ($hit{'down_chr'},$down_coords)=split(':',$down_string);
    $hit{'up_strand'}=substr($up_coords,0,1,'');
    $hit{'down_strand'}=substr($down_coords,0,1,'');

    # Depending of the way the strands are combined we must parse in different
    # ways
    my $parser_key=$hit{'up_strand'}.$hit{'down_strand'};
    if (exists $coord_parser->{$parser_key}){
	($hit{'up_start'},
	 $hit{'down_end'})=$coord_parser->{$parser_key}->($up_coords,
							  $down_coords,
							  $pos,
							  $read_length);
    } else {
	warn "Unknown strand combination $parser_key\n";
    }

    return(\%hit);
}

sub check_multi_start {
    my $coord_string=shift;
    my $pos=shift;
    my $read_length=shift;
    my $down=shift;
    my @ranges;

    my %starts;
    $pos=~s/(\[|])//g;
    $coord_string=~s/(\[|])//g;
    
    my @pos=sort {$a <=> $b} (split(/[-;]/,$pos));
    my ($pos_start,$pos_end)=($pos[0],$pos[-1]);
    my @splits=($pos_start .. $pos_end);
    if ($down) {
	foreach my $split (@splits) {
	    $split=$read_length - $split;
	}
	@splits=sort {$b <=> $a} @splits;
    } else {
	@splits=sort {$a <=> $b} @splits;
    }

    my @coords=split(/[-;]/,$coord_string);
    my @coords_sort=sort {$a <=> $b} @coords;
    my ($hit_start,$hit_end)=($coords_sort[0],$coords_sort[-1]);
    # account for negative coordinates. If this happens one of the numbers
    # will be undefined
    unless ($coords_sort[0]) {
	warn "WARNING: Possible problem parsing $pos\n";
	warn "I think there is a negative coordinate...\n";
	$hit_start=$hit_end;
    }
    my @starts=($hit_start .. $hit_end);

    if ($coords[1] &&
	($coords[0] > $coords[1])) {
	# This means we are in the minus strand, so we invert the array
	@starts=reverse(@starts);
    }

    unless (@splits == @starts) {
	warn "Wrong number of starts $pos, $coord_string\n";
	print STDERR join("\t",
			  @splits),"\n";
	print STDERR join("\t",
			  @starts),"\n";
	return();
    }

    for (my $i=0;$i<@splits;$i++) {
	push @ranges, [$starts[$i],$starts[$i] + $splits[$i] - 1];
    }

    return(\@ranges);
}

sub check_single_start {
    my $coord_string=shift;
    my $pos=shift;
    my $read_length=shift;
    my $down=shift;
    my @ranges;

    my %starts;
    $pos=~s/(\[|])//g;
    $coord_string=~s/(\[|])//g;

    my @pos=sort {$a <=> $b} (split(/[-;]/,$pos));
    my ($pos_start,$pos_end)=($pos[0],$pos[-1]);
    my @splits=($pos_start .. $pos_end);
    if ($down) {
	foreach my $split (@splits) {
	    $split= $read_length - $split;
	}
	@splits=sort {$b <=> $a} @splits;
    } else {
	@splits=sort {$a <=> $b} @splits;
    }
    my @starts=split(/[-;]/,$coord_string);
    @starts=($starts[0] .. $starts[-1]);

    unless ((@splits == @starts) ||
	    (@starts == 1)){
	warn "Wrong number of starts $pos, $coord_string\n";
	return();
    }

    for (my $i=0;$i<@splits;$i++) {
	my $start;
	if (@starts > 1) {
	    $start=$starts[$i];
	} else {
	    $start=$starts[0]
	}
	push @ranges, [$start,$start + $splits[$i] - 1];
    }

    return(\@ranges);
}

sub get_split_coord_parser {
    my %parser;

    $parser{'FF'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;
	my $read_length=shift;

	# Check input
	# First the up string
	my $ranges1=check_single_start($up_coord_string,
				       $pos,
				       $read_length);

	# Check the down string
	my $ranges2=check_multi_start($down_coord_string,
				      $pos,
				      $read_length,
				      1);

	return($ranges1,$ranges2);
    };
    
    $parser{'FR'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;
	my $read_length=shift;

	# Check input
	# First the up string
	my $ranges1=check_single_start($up_coord_string,
				       $pos,
				       $read_length);

	# Check the down string
	my $ranges2=check_single_start($down_coord_string,
				       $pos,
				       $read_length,
				       1);

	return($ranges1,$ranges2);
    };

    $parser{'RR'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;
	my $read_length=shift;

	# Check input
	# First the up string
	my $ranges1=check_multi_start($up_coord_string,
				      $pos,
				      $read_length);

	# Check the down string
	my $ranges2=check_single_start($down_coord_string,
				       $pos,
				       $read_length,
				       1);

	return($ranges1,$ranges2);
    };
    
    $parser{'RF'}=sub{
	my $up_coord_string=shift;
	my $down_coord_string=shift;
	my $pos=shift;
	my $read_length=shift;

	# Check input
	# First the up string
	my $ranges1=check_multi_start($up_coord_string,
				      $pos,
				      $read_length);

	# Check the down string
	my $ranges2=check_multi_start($down_coord_string,
				      $pos,
				      $read_length,
				      1);

	return($ranges1,$ranges2);
    };

    return(\%parser);
}

sub get_split_parser {
    my $mismatches=shift;
    my $parser;
    my %stats;

    $parser=sub {
	my $infile=shift;
	my $lane=shift;

	unless($infile) {
	    die "A transcriptome mapping file is required\n";
	}

	my $outfileunique=$lane.'.single.unique.gtf.gz';
	my $outfilemulti=$lane.'.single.multi.gtf.gz';

	print STDERR "Printing results to $lane.single.[unique|multi].gtf.gz\n";

	# Open the required files
	# Output files
	my $outuniquefh=get_fh($outfileunique,1);
	my $outmultifh=get_fh($outfilemulti,1);

	# Inputfiles
	my $infh=get_fh($infile);

	my %strands;
	my $read_length;
	print STDERR 'Parsing...';
	while (my $line=<$infh>) {
	    my %line=%{parse_gem_line($line)};
	    $read_length=$line{'length'};
	    my @map=split(':',$line{'matches'});
	    if (@map < $mismatches - 1) {
		$mismatches=@map - 1;
		warn "File was mapped with $mismatches mismatches, reducing mismatch threshold accordingly \n";
	    }

	    # Get best match position
	    my $best_match=0;
	    for (my $i=0;$i<=$mismatches;$i++) {
		# Cope with both possibilities (! and -) in case the mapper
		# the mapper version is older
		if ($map[$i] eq '-') {
		    $best_match=-1;
		    last;
		} elsif ($map[$i] eq '!') {
		    $best_match=-1;
		    last;
		} elsif ($map[$i] != 0) {
		    $best_match=$i;		    
		    last;
		} 
	    }

	    # Skip hits with too many matches
	    if ($best_match < 0) {
		$stats{'too_many'}++;
		next;
	    }

	    # From the @map array get the number of matches we have to take
	    # First remove the first part of the @map array up to the first
	    # match
	    splice(@map,0,$best_match);

	    # After remove the second part from the position after the first
	    # match plus the closest hit number we want
	    my $offset=$closest_hit + 1;
	    unless ($offset > @map) {
		splice(@map,$closest_hit + 1);
	    }
	    my $matches=0;
	    foreach my $hits (@map) {
		$matches+=$hits;
	    }

	    if ($matches < 1) {
		# If there are no matches skip the rest
		$stats{'no_hits'}++;
		next;
	    } elsif ($matches > 1) {
		# Unless we are going for multi matches, if the best match is
		# not unique discard it.
		$stats{'multi_hits'}++;
		unless ($multi) {
		    next;
		}
	    } else {
		$stats{'unique_hits'}++;
	    }

	    # If the script reaches this point it has found at least one good
	    # hit
	    my @matches=split(',',$line{'hits'});
	    if ($matches > @matches) {
		if (@matches) {
		    $stats{'missing_close'}++;
		    $matches=@matches;
		} else {
		    $stats{'too_many'}++;
		    next;
		}
	    }
	    splice(@matches,$matches);

	    # Parse the hit in this case we are taking only the best one
	    if ($matches == 1) {
		my ($pos,$coord_string)=split('=',$matches[0],2);
		my %hit=%{get_coords($coord_string,
				     $pos,
				     $read_length)};
		
		$strands{$hit{'up_strand'}.$hit{'down_strand'}}++;
		
		# Print two Gtf entries per split read using the same id to be able
		# to group them after
		my $read_id=join('; ',
				 'ID "'.$line{'id'}.'"',
				 'Identity "'.$line{'matches'}.'"',
				 'Split "'.$pos.'"');
		
		# set the strands to + or -
		my $up_strand=$hit{'up_strand'};
		if ($up_strand eq 'F') {
		    $hit{'up_strand'}='+';
		} elsif ($up_strand eq 'R') {
		    $hit{'up_strand'}='-';
		} else {
		    warn "Unknown strand $read_id\n";
		}
		my $down_strand=$hit{'down_strand'};
		if ($down_strand eq 'F') {
		    $hit{'down_strand'}='+';
		} elsif ($down_strand eq 'R') {
		    $hit{'down_strand'}='-';
		} else {
		    warn "Unknown strand $read_id\n";
		}
		
		$infile=~s/.*\///;
		my @range1;
		my @range2;
		for (my $i=0;$i<@{$hit{'up_start'}};$i++) {
		    unless ($hit{'up_start'}->[$i] &&
			    $hit{'down_end'}->[$i]) {
			warn "Problem parsing line: $i $line\n";
			print STDERR join("\t",@{$hit{'up_start'}->[$i]}),"\n";
			next;
		    }
		    if (@range1) {
			if ($range1[0] > $hit{'up_start'}->[$i]->[0]) {
			    $range1[0] = $hit{'up_start'}->[$i]->[0];
			}
			if ($range1[1] < $hit{'up_start'}->[$i]->[1]) {
			    $range1[1] = $hit{'up_start'}->[$i]->[1];
			}
		    } else {
			@range1=@{$hit{'up_start'}->[$i]};
		    }
		    if (@range2) {
			if ($range2[0] > $hit{'down_end'}->[$i]->[0]) {
			    $range2[0] = $hit{'down_end'}->[$i]->[0];
			}
			if ($range2[1] < $hit{'down_end'}->[$i]->[1]) {
			    $range2[1] = $hit{'down_end'}->[$i]->[1];
			}
		    } else {
			@range2=@{$hit{'down_end'}->[$i]};
		    }
		}
		print $outuniquefh join("\t",
					$hit{'up_chr'},
					$infile,
					'splitread',
					@range1,
					'.',
					$hit{'up_strand'},
					'.',
					$read_id),"\n";
		print $outuniquefh join("\t",
					$hit{'down_chr'},
					$infile,
					'splitread',
					@range2,
					'.',
					$hit{'down_strand'},
					'.',
					$read_id),"\n";
	    } elsif ($matches > 1) {
		foreach my $match (@matches) {
		    my ($pos,$coord_string)=split('=',$match,2);
		    my %hit=%{get_coords($coord_string,
					 $pos,
					 $read_length)};
		
		    $strands{$hit{'up_strand'}.$hit{'down_strand'}}++;
		
		    # Print two Gtf entries per split read using the same id to be able
		    # to group them after
		    my $read_id=join('; ',
				     'ID "'.$line{'id'}.'"',
				     'Identity "'.$line{'matches'}.'"',
				     'Split "'.$pos.'"');
		
		    # set the strands to + or -
		    my $up_strand=$hit{'up_strand'};
		    if ($up_strand eq 'F') {
			$hit{'up_strand'}='+';
		    } elsif ($up_strand eq 'R') {
			$hit{'up_strand'}='-';
		    } else {
			warn "Unknown strand $read_id\n";
		    }
		    my $down_strand=$hit{'down_strand'};
		    if ($down_strand eq 'F') {
			$hit{'down_strand'}='+';
		    } elsif ($down_strand eq 'R') {
			$hit{'down_strand'}='-';
		    } else {
			warn "Unknown strand $read_id\n";
		    }
		    
		    $infile=~s/.*\///;
		    my @range1;
		    my @range2;
		    for (my $i=0;$i<@{$hit{'up_start'}};$i++) {
			unless ($hit{'up_start'}->[$i] &&
				$hit{'down_end'}->[$i]) {
			    warn "Problem parsing line: $i $line\n";
			    print STDERR join("\t",@{$hit{'up_start'}->[$i]}),"\n";
			    next;
			}
			if (@range1) {
			    if ($range1[0] > $hit{'up_start'}->[$i]->[0]) {
				$range1[0] = $hit{'up_start'}->[$i]->[0];
			    }
			    if ($range1[1] < $hit{'up_start'}->[$i]->[1]) {
				$range1[1] = $hit{'up_start'}->[$i]->[1];
			    }
			} else {
			    @range1=@{$hit{'up_start'}->[$i]};
			}
			if (@range2) {
			    if ($range2[0] > $hit{'down_end'}->[$i]->[0]) {
				$range2[0] = $hit{'down_end'}->[$i]->[0];
			    }
			    if ($range2[1] < $hit{'down_end'}->[$i]->[1]) {
				$range2[1] = $hit{'down_end'}->[$i]->[1];
			    }
			} else {
			    @range2=@{$hit{'down_end'}->[$i]};
			}
		    }
		    print $outmultifh join("\t",
					   $hit{'up_chr'},
					   $infile,
					   'splitread',
					   @range1,
					   '.',
					   $hit{'up_strand'},
					   '.',
					   $read_id),"\n";
		    print $outmultifh join("\t",
					   $hit{'down_chr'},
					   $infile,
					   'splitread',
					   @range2,
					   '.',
					   $hit{'down_strand'},
					   '.',
					   $read_id),"\n";
		}
	    } else {
		warn "I should not be here\n";
	    }
	}
	close($outuniquefh);
	close($outmultifh);

	print STDERR "done\n";
    };
    
    return($parser);
}

sub get_lane_files {
    my $files=shift;

    my %lane_files;

    foreach my $file (keys %{$files}) {
	my $pair=$files->{$file}->[0];
	my $lane=$files->{$file}->[1];
	my $splitfile=get_split_file($lane);
	my $filepath=$options{'SPLITMAPDIR'}.'/'.$splitfile;
	if (-r $filepath) {
	    $lane_files{$pair}{$lane}=$filepath;
	} elsif (-r $filepath.'.gz') {
	    print STDERR "$splitfile is gzipped\n";
	    $lane_files{$pair}{$lane}=$filepath.'.gz';
	} else {
	    die "Can't find $splitfile or $splitfile.gz\n";
	}
    }

    return(\%lane_files);
}
