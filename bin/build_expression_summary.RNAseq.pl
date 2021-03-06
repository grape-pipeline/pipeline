#!/soft/bin/perl

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
# This script should take the information from the annotation and determine
# from the annotated features what fractions are detected in each lane and the
# total

use RNAseq_pipeline3 qw(get_fh parse_gff_line check_table_existence);
use RNAseq_pipeline_settings3 ('read_config_file','get_dataset_id',
			       'get_dbh');

# Declare some variables
my @files;
my %expressed;

my $annotation;
my $prefix;
my $gene_rpkm_table;
my $trans_rpkm_table;
my $exon_rpkm_table;
my $exonclasstab;
my $expression_summary_table;
my $datasets_table;
my $db;
my $threshold=1;

my %options=%{read_config_file()};
$annotation=$options{'ANNOTATION'};
$prefix=$options{'PREFIX'};
$gene_rpkm_table=$prefix.'_gene_RPKM_pooled';
$trans_rpkm_table=$prefix.'_transcript_expression_levels_pooled';
$exon_rpkm_table=$prefix.'_exon_RPKM_pooled';
$exonclasstab=$options{'EXONSCLASSTABLE'};
$expression_summary_table=$prefix.'_expression_summary';
$datasets_table=$prefix.'_dataset';
$db=$options{'DB'};

my $dbh=get_dbh();
my $dbhcommon=get_dbh(1);

my %distribution;

my %lanes=%{get_dataset_id($dbh,
			   $datasets_table)};

# Get the number of detected features in each lane
my @lanes=sort keys %lanes;

# Get the detected features for each of the lanes
for (my $i=0;$i<@lanes;$i++) {
    my $gene_no=get_detected_genes($gene_rpkm_table,
				    $dbh,
				    $lanes[$i],
				    $threshold) || 0;
    my $trans_no=get_detected_transcripts($trans_rpkm_table,
					  $dbh,
					  $lanes[$i],
					  $threshold) || 0;
    my $exon_no=get_detected_exons($exon_rpkm_table,
				   $dbh,
				   $lanes[$i],
				   $threshold) || 0;
    $distribution{'genes'}[$i]=$gene_no;
    $distribution{'transcripts'}[$i]=$trans_no;
    $distribution{'exons'}[$i]=$exon_no;
}

# Get the total detected features after pooling the lanes
my $total_gene_no=get_detected_genes($gene_rpkm_table,
				     $dbh,
				     '',
				     $threshold) || 0;
my $total_trans_no=get_detected_transcripts($trans_rpkm_table,
					    $dbh,
					    '',
					    $threshold) || 0;
my $total_exon_no=get_detected_exons($exon_rpkm_table,
				     $dbh,
				     '',
				     $threshold) || 0;
$distribution{'total'}=[$total_gene_no,
			$total_trans_no,
			$total_exon_no];

# get the total number of genes transcripts and exons in the annotation
my %features;
$features{'genes'}=get_gene_number($dbhcommon,
				   $exonclasstab);
$features{'transcripts'}=get_trans_number($dbhcommon,
					  $exonclasstab);
$features{'exons'}=get_exon_number($dbhcommon,
				   $exonclasstab);

# print out the results
my $taboutfh=get_fh($expression_summary_table.'.txt',1);
my $gene_frac=($distribution{'total'}[0] * 100) / $features{'genes'} || '-';
my $trans_frac=($distribution{'total'}[1] * 100) / $features{'transcripts'} || '-';
my $exons_frac=($distribution{'total'}[2] * 100) / $features{'exons'} || '-';
print $taboutfh join("\t",
		     'Genes',
		     $features{'genes'},
		     $distribution{'total'}[0],
		     $gene_frac,
		     @{$distribution{'genes'}}),"\n";
print $taboutfh join("\t",
		     'Transcripts',
		     $features{'transcripts'},
		     $distribution{'total'}[1],
		     $trans_frac,
		     @{$distribution{'transcripts'}}),"\n";
print $taboutfh join("\t",
		     'Exons',
		     $features{'exons'},
		     $distribution{'total'}[2],
		     $exons_frac,
		     @{$distribution{'exons'}}),"\n";
close($taboutfh);

# Build the table in the database
build_db_table($dbh,
	       $expression_summary_table,
	       \@lanes,
	       \%lanes);

### TO DO
# Turn this into a subroutine
# Insert into the database
my $command="mysql $db < $expression_summary_table.sql";
print STDERR "Executing: $command\n";
system($command);
$command="mysqlimport -L $db $expression_summary_table.txt";
print STDERR "Executing: $command\n";
system($command);
$command="rm $expression_summary_table.sql $expression_summary_table.txt";
print STDERR "Executing: $command\n";
system($command);

exit;

sub build_db_table {
    my $dbh=shift;
    my $table=shift;
    my $lanes=shift;
    my $lane2id=shift;
    my $outtable=$table.'.sql';

    my ($query,$sth);

    $query ="DROP TABLE IF EXISTS $table;";
    $query.="CREATE TABLE $table( ";
    $query.='type varchar(50) NOT NULL,';
    $query.='total mediumint unsigned NOT NULL,';
    $query.='detected mediumint unsigned NOT NULL,';
    $query.='fraction float unsigned NOT NULL,';
    foreach my $lane (@{$lanes}) {
	my $lane_id=$lane2id->{$lane};
	$lane_id=$lane;
	if ($lane_id=~/^\d*\.\d*$/) {
	    $lane_id='l'.$lane_id;
	}
	$query.="$lane_id mediumint unsigned NOT NULL,";
	$query.="INDEX idx_${lane_id} ($lane_id),";
    }
    $query.='INDEX idx_type (type)';
    $query.=');';

    my $outfh=get_fh($outtable,1);
    print $outfh $query,"\n";
    close($outtable);
}

sub get_exon_number {
    my $dbh=shift;
    my $table=shift;

    my $exon_number=0;

    if (check_table_existence($dbh,$table)) {
	my ($query,$sth);
	$query ='SELECT count(DISTINCT exon_id) ';
	$query.="FROM $table";
	$sth=$dbh->prepare($query);
	$sth->execute();
	($exon_number)=$sth->fetchrow_array();
    }

    return($exon_number);
}

sub get_trans_number {
    my $dbh=shift;
    my $table=shift;

    my $trans_number=0;

    if (check_table_existence($dbh,$table)) {
	my ($query,$sth);
	$query ='SELECT count(DISTINCT transcript_id) ';
	$query.="FROM $table";
	$sth=$dbh->prepare($query);
	$sth->execute();
	($trans_number)=$sth->fetchrow_array();
    }

    return($trans_number);
}

sub get_gene_number {
    my $dbh=shift;
    my $table=shift;

    my $gene_number=0;

    if (check_table_existence($dbh,$table))  {
	my ($query,$sth);
	$query ='SELECT count(DISTINCT gene_id) ';
	$query.="FROM $table";
	$sth=$dbh->prepare($query);
	$sth->execute();
	($gene_number)=$sth->fetchrow_array();
    }
    
    return($gene_number);
}

sub get_detected_genes {
    my $table=shift;
    my $dbh=shift;
    my $lane=shift;
    my $threshold=shift || 1;
    
    my $detected=0;

    if (check_table_existence($dbh,$table)) {
	my ($query,$sth);
	$query ='SELECT count(DISTINCT gene_id) ';
	$query.="FROM $table ";
	if ($lane) {
	    $query.='WHERE Sample = ?';
	}
	$sth=$dbh->prepare($query);
	if ($lane) {
	    $sth->execute($lane);
	} else {
	    $sth->execute();
	}
	($detected)=$sth->fetchrow_array();
    }
    
    return($detected)
}

sub get_detected_transcripts {
    my $table=shift;
    my $dbh=shift;
    my $lane=shift;
    my $threshold=shift || 1;
    
    my $detected=0;

    if (check_table_existence($dbh,$table)) {
	my ($query,$sth);
	$query ='SELECT count(DISTINCT transcript_id) ';
	$query.="FROM $table ";
	if ($lane) {
	    $query.='WHERE Sample = ?';
	}
	$sth=$dbh->prepare($query);
	if ($lane) {
	    $sth->execute($lane);
	} else {
	    $sth->execute();
	}
	($detected)=$sth->fetchrow_array();
    }
    
    return($detected)
}

sub get_detected_exons {
    my $table=shift;
    my $dbh=shift;
    my $lane=shift;
    my $threshold=shift || 1;
    
    my $detected=0;

    if (check_table_existence($dbh,$table)) {
	my ($query,$sth);
	$query ='SELECT count(DISTINCT exon_id) ';
	$query.="FROM $table ";
	if ($lane) {
	    $query.='WHERE Sample = ?';
	}
	$sth=$dbh->prepare($query);
	if ($lane) {
	    $sth->execute($lane);
	} else {
	    $sth->execute();
	}
	($detected)=$sth->fetchrow_array();
    }
    
    return($detected)
}


