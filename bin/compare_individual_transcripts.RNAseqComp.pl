#!/soft/bin/perl
# DGK

#    GRAPE
#    Copyright (C) 2011 Centre for Genomic Regulation (CRG)
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

# Objective
# This script should take a pair of gene ids and it will extract for each of
# Them the expression values of each of the corresponding transcripts
# Also correlation between these two genes across all the datasets in the
# database will be calculated.
# The value used for calulating the correlation will be provided by using the
# suffix of the table containing it

use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file','get_dbh',
			       'get_trans_expression_data',
			       'get_gene_from_trans_sub');
use RNAseq_pipeline_comp3 ('get_tables','check_tables','get_labels_sub',
			   'get_samples');
use Getopt::Long;

# Declare some variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $debug=1;
my $breakdown;
my $transfile;
my $tabsuffix='transcript_expression_levels_pooled';
my @trans_needed;
my $all=0;

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'breakdown|b' => \$breakdown,
	   'transcriptfile|f=s' => \$transfile,
	   'trans|t=s' => \@trans_needed,
	   'all|a' => \$all);

if ($breakdown) {
    $tabsuffix='transcript_expression_levels';
}

# Get subs
*gene2trans=gene2trans_sub();
*trans2gene=get_gene_from_trans_sub();

my %trans2gene;
if ($transfile) {
    my $transfh=get_fh($transfile);
    while (my $line=<$transfh>) {
	chomp($line);
	my ($id,$type)=split("\t",$line);
	if ($type &&
	    $type eq 'gene') {
	    # Get all the transcripts from the gene of interest
	    my @transcripts=@{gene2trans($id)};
	    push @trans_needed,@transcripts;
	} else {
	    push @trans_needed, $id;
	}
    }
    close($transfh);
}

unless (@trans_needed >= 1) {
    die "I need at least 1 transcript (provided with the -t option\n";
}

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_ind_trans.RNAseqComp.log',
		      $debug);
print $log_fh "Extracting expression info for the following transcripts: ",
    join(",",@trans_needed),"\n";

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get subroutines
*get_labels=get_labels_sub($dbhcommon);

# Get all experiment tables from the database
my %tables=%{get_tables($dbhcommon,
			$project,
			$tabsuffix,
			'',
			$all)};

# Remove any tables that do not exist
check_tables($dbh,
	     \%tables);

# For each of tables extract the RPKMs of interest
my %samples=%{get_samples(\%tables,
			  $dbh,
			  $breakdown)};
my @experiments=sort {$a cmp $b} keys %samples;
my @values;
my %all_trans;
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment,2);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_trans_expression_data($dbh,
				       $table,
				       \%all_trans,
				       $sample,
				       $breakdown);
    if ($data) {
	push @values, [$experiment,$data];
    } else {
	print STDERR "Skipping $experiment\n";
    }
}

# Get the human readable lables
foreach my $experiment (@experiments) {
    my $label;
    if ($nolabels) {
	$label=$samples{$experiment}->[1];
    } else {
	$label=get_labels($experiment);
    }
    if ($label) {
	$experiment=$label;
    }
}

# Print the expression values for each gene of interest
my $tmpfn="Trans.Expression.subset.txt";
my $tmpfh=get_fh($tmpfn,1);
foreach my $exp (@values) {
    my @row;
    my $dataset=$exp->[0];
    $dataset=~s/_transcript_expression_levels_pooled_sample//;
    my ($project,$group)=(split('_',$dataset))[0,-1];
    foreach my $trans (@trans_needed) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$trans})) {
	    $value=$exp->[1]->{$trans};
	}
	my $gene_id=trans2gene($trans);

	print $tmpfh join("\t",
			  $dataset,
			  $project,
			  $group,
			  $gene_id,
			  $trans,
			  $value),"\n";
    }
}
close($tmpfh);

exit;

sub gene2trans_sub {
    my %options=%{read_config_file()};
    my $table=$options{'EXONSCLASSTABLE'} || die "No exons table defined\n";

    my $dbh=get_dbh(1);
    my %cache;

    my ($query,$sth);
    $query ='SELECT distinct transcript_id ';
    $query.="FROM $table ";
    $query.='WHERE gene_id = ?';
    $sth=$dbh->prepare($query);
    
    my $gene2trans=sub {
	my $gene=shift;

	unless ($cache{$gene}) {
	    $sth->execute($gene);
	    while (my ($trans_id)=$sth->fetchrow_array()) {
		push @{$cache{$gene}},$trans_id;
	    }
	}
	return($cache{$gene});
    };

    return($gene2trans);
}
