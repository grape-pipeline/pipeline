#!/soft/bin/perl

use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

# Objective
# This script will take the output of the flux capacitor and it will create a
# table to load in the database with this information

use Getopt::Long;
use RNAseq_pipeline3 qw(get_fh get_log_fh parse_gff_line);
use RNAseq_pipeline_settings3 qw(read_file_list get_gene_from_trans_sub);

# declare some variables
my $debug=0;

# get subroutines
*trans2gene=get_gene_from_trans_sub();

# Get the log_fh
my $log_fh=get_log_fh('build_transcript_expression_levels.RNAseq.log',
		      $debug);

# First get a list of the bed files we are going to process
my %files=%{read_file_list()};
my %lanes=%{get_lanes(\%files)};

# read the results from the files named after the lane ids
foreach my $lane (keys %lanes) {
    my $filename=$lane.'.flux.gtf';

    if (-r $filename) {
	process_file($filename,
		     $lane,
		     $log_fh);
	my $command="rm $filename";
	print STDERR "Executing:\t$command\n";
	system($command);
    } else {
	die "I can't find $filename\n";
    }
}

close($log_fh);

exit;

sub process_file {
    my $infile=shift;
    my $lane=shift;
    my $log_fh=shift;

    my $infh=get_fh($infile);

    while (my $line=<$infh>) {
	my %line=%{parse_gff_line($line,
				  $log_fh)};

	my $type=$line{'type'};

	if ($type=~/transcript/) {
	    my $flux_id=$line{'feature'}{'gene_id'};
	    my $trans_id=$line{'feature'}{'transcript_id'};
	    my $rpkm=$line{'feature'}{'RPKM'};
	    my $gene_id=trans2gene($trans_id);
	    if ($rpkm > 0) {
		print join("\t",
			   $gene_id,
			   $trans_id,
			   $flux_id,
			   $rpkm,
			   $lane),"\n";
	    }
	}
    }
    close($infh);
}

sub get_lanes {
    my $files=shift;
    my %lanes;

    foreach my $file (keys %{$files}) {
	push @{$lanes{$files->{$file}->[0]}},$files->{$file}->[1];
    }

    return(\%lanes);
}
