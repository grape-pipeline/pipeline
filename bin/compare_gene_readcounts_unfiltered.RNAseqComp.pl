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
# This script should take a project from the database and for each of the
# experiments belonging to it build a comparison table that can be run
# through R
# We have to be able to decide what table to get the data from (pooled or not)
# And also we have to be able to select genes that are expressed 
# In order to get set the descriptions after using the script get_EnsEMBL_gene_info.pl and a list of genes:
# gawk -F"\t" '{print "UPDATE 1_H_sapiens_EnsEMBL_55_parsed_gt_geneclass set description=\""$4"\" WHERE gene_id=\""$1"\";"}' all.genes.desc.txt > add.description.sql


use RNAseq_pipeline3 qw(get_fh get_log_fh run_system_command get_list);
use RNAseq_pipeline_settings3 ('get_dbh','read_config_file',
			       'get_gene_readcount_data','get_gene_info_sub');
use RNAseq_pipeline_comp3 ('get_tables','check_tables','get_labels_sub',
			   'get_samples','remove_tables');
use Getopt::Long;

# Declare some variables
my $nolabels;
my $dbh;
my $dbhcommon;
my $project;
my $debug=1;
my $tabsuffix='gene_readcount_pooled';
my $fraction='';
my $subset='';

# Get command line options
GetOptions('nolabels|n' => \$nolabels,
	   'debug|d' => \$debug,
	   'limit|l=s' => \$fraction,
	   'subset|s=s' => \$subset);

# read the config file
my %options=%{read_config_file()};
$project=$options{'PROJECTID'};

# get a log file
my $log_fh=get_log_fh('compare_gene_readcounts.RNAseqComp.log',
		      $debug);

# First connect to the database
$dbh=get_dbh();
$dbhcommon=get_dbh(1);

# Get subroutines
*get_labels=get_labels_sub($dbhcommon);
*gene2chr=get_gene_info_sub('chr');
*gene2desc=get_gene_info_sub('description');
*gene2type=get_gene_info_sub('type');

# Get the tables belonging to the project
my %tables=%{get_tables($dbhcommon,
			$project,
			$tabsuffix,
			$fraction)};

# Remove any tables that do not exist
check_tables($dbh,
	     \%tables);

# If a subset has been provided remove any tables that are not included in the
# subset
if ($subset && -r $subset) {
    my %subset=%{get_list($subset)};
    remove_tables(\%tables,
		  \%subset);
}

# For each of tables extract the RPKMs of interest and get for each of the
# tables the different samples present in them
my %samples=%{get_samples(\%tables,
			  $dbh,
			  1)};# currently set to one until we fix the table naming problem
my @experiments=sort {$a cmp $b} keys %samples;
my @values;
my %all_genes;
foreach my $experiment (@experiments) {
    my ($table,$sample)=split('_sample_',$experiment);
    print $log_fh "Extracting $sample, data from $table\n";
    my $data=get_gene_readcount_data($dbh,
				     $table,
				     \%all_genes,
				     $sample,
				     1); # currently set to one until we fix the table naming problem
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

# Print the expression values for each gene in each of the tables into a
# temporary file
my $tmpfn="Gene.ReadCount.unfiltered.$project.txt";
my $tmpfh=get_fh($tmpfn,1);
print $tmpfh join("\t",@experiments),"\n";
foreach my $gene (keys %all_genes) {
    my @row;
    my $no_print=0;
    foreach my $exp (@values) {
	my $value=0;
	if ($exp->[1] &&
	    ($exp->[1]->{$gene})) {
	    $value=$exp->[1]->{$gene};
	}
	push @row,$value;
    }

    # Skip mitochondrial and ribosomal genes
    my $desc=gene2desc($gene) || '';
    print $tmpfh join("\t",
		      $gene,
		      @row),"\n";
}
close($tmpfh);

exit;

