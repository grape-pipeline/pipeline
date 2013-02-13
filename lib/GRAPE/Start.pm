package GRAPE::Start;

#  GRAPE
#  Copyright (C) 2008-2011 Centre for Genomic Regulation (CRG)
#
#  This file is part of GRAPE.
#
#  Author: David Gonzalez
#

# Export subroutines to the caller namespace
use Exporter;
@ISA=('Exporter');
@EXPORT_OK=('base_table_build','get_existing_data_subs');
push @EXPORT_OK,('get_tables_hash','build_table_files');
push @EXPORT_OK,('create_directory_structure');
push @EXPORT_OK,('print_config_file','print_pipeline_file','build_file_list');
push @EXPORT_OK,('add_project','add_experiment','add_proj_info','add_exp_info');
push @EXPORT_OK,('clear_experiment');

# Set strict and warnings
use strict;
use warnings;

# Add the path to the library to be used
BEGIN {
    use Cwd 'abs_path';
    my $libdir=abs_path($0);
    $libdir=~s/bin\/[^\/]*$/lib/;
    unshift @INC, "$libdir";
}

# Load other modules required
## POSIX 
use POSIX qw(uname);
use Cwd;
use Cwd 'abs_path';

## Batabase interface
use DBI;
use DBD::mysql;

## Bio::Perl modules
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::SeqFeature::Gene::GeneStructure;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Exon;

## Pipeline specific modules
use GRAPE::Base ('get_fh','MySQL_DB_Connect','check_file_existence',
		 'run_system_command',
		 'check_table_existence','check_field_value');
use GRAPE::Formats::GFF ('check_gff_file');
use GRAPE::Formats::Fasta ('check_fasta_file');

### Subroutines

### Subroutines for checking the options as well as what is present already in
# the database
# Check the options that are provided to the script in order to decide if there 
# are any unadvisable characters
## OK


## Subroutines used for building the tables for the common database and
# populating it
# Build the tables that will be shared between the projects and experiments


# Actual table syntax
sub base_table_build {
    my %tables=('projects' => 'CREATE TABLE IF NOT EXISTS projects (
       project_id varchar(50) NOT NULL,
       species varchar(50) NOT NULL,
       proj_description mediumtext,
       secret varchar(50) NULL,
       PRIMARY KEY (project_id)
);',
		'experiments' => 'CREATE TABLE IF NOT EXISTS experiments (
       experiment_id varchar(50) NOT NULL,
       project_id varchar(50) NOT NULL REFERENCES projects(project_id),
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       genome_id mediumint unsigned NOT NULL REFERENCES genome_files(genome_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       template_file varchar(200) NOT NULL,
       read_length smallint unsigned NULL,
       mismatches smallint unsigned NOT NULL,
       exp_description mediumtext,
       expDate date NULL,
       CellType varchar(50) NULL,
       RNAType varchar(50) NULL,
       Compartment varchar(50) NULL DEFAULT "CELL",
       Bioreplicate varchar(10) DEFAULT 1,
       partition varchar(50) NULL,
       md5sum varchar(45) NULL,
       Preprocessing varchar(45) NULL,
       annotation_version varchar(45) NULL,
       genome_assembly varchar(45) NULL,
       paired bit(1) DEFAULT 1, 
       lab varchar(45) NULL,
       protocol_id mediumint unsigned DEFAULT 1,
       PRIMARY KEY (project_id,experiment_id)
);',
		'annotation_files' => 'CREATE TABLE IF NOT EXISTS annotation_files (
       annotation_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       annotation varchar(60) NOT NULL,
       location varchar(200) NOT NULL,
       version varchar(50) NULL,
       source mediumtext NULL,
       index idx_annotation (annotation),
       PRIMARY KEY (annotation_id)
) ENGINE=INNODB;',
		'genome_files' => 'CREATE TABLE IF NOT EXISTS genome_files (
       genome_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       genome varchar(60) NOT NULL,
       location varchar(200) NOT NULL,
       assembly varchar(50) NULL,
       source mediumtext NULL,
       gender varchar(20) NULL,
       index idx_genome (genome),
       PRIMARY KEY (genome_id)
);',
		'indices' => 'CREATE TABLE IF NOT EXISTS indices (
       index_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       genome_id mediumint unsigned NOT NULL REFERENCES genome_files(genome_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       type varchar(100) NOT NULL,
       location varchar(200) NOT NULL,
       PRIMARY KEY (index_id)
) ENGINE=INNODB;',	
		'exclusion_files' => 'CREATE TABLE IF NOT EXISTS exclusion_files (
       exclusion_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       species_id mediumint unsigned NOT NULL REFERENCES species_info(species_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       location varchar(200) NOT NULL,
       PRIMARY KEY (exclusion_id)
);',
		'annotation_tables' => 'CREATE TABLE IF NOT EXISTS annotation_tables (
       genome_id mediumint unsigned NOT NULL REFERENCES genome_files(genome_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       type varchar(20) NOT NULL,
       table_name varchar(100) NOT NULL,
       index idx_genome (genome_id),
       index idx_annotation (annotation_id),
       index idx_type (type)
) ENGINE=INNODB;',
		'fasta_files' => 'CREATE TABLE IF NOT EXISTS fasta_files (
       genome_id mediumint unsigned NOT NULL REFERENCES genome_files(genome_id),
       annotation_id mediumint unsigned NOT NULL REFERENCES annotation_files(annotation_id),
       filetype varchar(20) NOT NULL,
       file_name varchar(200) NOT NULL,
       index idx_genome (genome_id),
       index idx_annotation (annotation_id),
       index idx_type (filetype),
       PRIMARY KEY (genome_id,annotation_id,filetype)
);',
		'species_info' => 'CREATE TABLE IF NOT EXISTS species_info (
       species_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       species varchar(200) NOT NULL,
       genus varchar(100) NOT NULL,
       sp_alias varchar(100) NOT NULL,
       abbreviation varchar(20) NOT NULL,
       PRIMARY KEY (species_id),
       index idx_alias (sp_alias)
);',
		'protocol_info' => 'CREATE TABLE IF NOT EXISTS protocol_info (
       protocol_id mediumint unsigned NOT NULL AUTO_INCREMENT,
       abbreviation varchar(20) NULL,
       description varchar(200) NOT NULL,
       PRIMARY KEY (protocol_id)
);'
);
    return(\%tables);
}


# This script should for each of the necessary files
sub guess_values {
    my $guessmaster=shift;
    my $options=shift;

    print STDERR "Some values missing; time for some guesswork...\n";
    foreach my $key (keys %{$guessmaster}) {
	print STDERR "Guessing $key ...";
	my $value=guess($key);
	${$options->{$key}}=$value;
	print STDERR "done\n";
    }


    print STDERR "done\n";
}

## The following subs will build the tables for each of the steps in the
# pipeline
sub get_tables_hash {
    my $prefix=shift;
    # Build a hash with the table root as the key and the name as the value
    my %tables=('_start' => '',
		'_read_stats' => '',
		'_dataset' => '',
		'_junctions' => '',
		'_genes' => '',
		'_detected_genes' => '',
		'_transcripts' => '',
		'_detected_transcripts' => '',
		'_exon_seqs' => '',
		'_junction_seqs' => '',
		'_transcript_seqs' => '',
		'_indices' => '',
		'_genome_mapping' => '',
		'_transcriptome_mapping' => '',
		'_junctions_mapping' => '',
		'_unmapped_reads' => '',
		'_split_mapping' => '',
		'_recursive_mapping' => '',
		'_unique_maps_genome' => '',
		'_unique_maps_transcripts' => '',
		'_unique_maps_junctions' => '',
		'_unique_maps_split' => '',
		'_initial_clusters' => '',
		'_split_mapping_breakdown' => '',
		'_extract_annotation' => '',
		'_gene_coverage' => '',
		'_proj_coverage' => '',
		'_exon_coverage' => '',
		'_junction_coverage' => '',
		'_splicing_summary' => '',
		'_exon_RPKM_pooled' => '',
		'_gene_readcount_pooled' => '',
		'_gene_RPKM_pooled' => '',
		'_gene_RPKM_dist' => '',
		'_EJEI' => '',
		'_store_reads' => '',
		'_exon_classification' => '',
		'_junction_classification' => '',
		'_novel_junctions_summary' => '',
		'_all_junctions_class' => '',
		'_all_junctions_class_pooled' => '',
#		'_transcript_expression_levels' => '',
		'_transcript_expression_levels_pooled' => '',
#		'_exon_inclusion' => '',
		'_exon_inclusion_reads' => '',
		'_exon_inclusion_pooled' => '',
		'_inclusion_correlation' => '',
		'_inclusion_dist' => '',
		'_merged_SAM' => '',
		'_bed_files' => '',
		'_fusion_transcripts' => '',
		'_fusion_transcripts_support' => '',
		'_read_classification' => '',
		'_qualitiespos' => '',
		'_ambiguous' => '',
		'_read_dist_transcripts' => '',
		'_merged_mapping' => '',
		'_junction_maps_class' => '',
		'_completion_status' => '',
		'_summaries' => '',
		'_register_results' => '',
		'_gene_RPKM_pooled_flux' => ''
		);

    # Add the species prefix to the table
    foreach my $table (keys %tables) {
	$tables{$table}=$prefix.$table;
    }
    return(\%tables);
}

# This sub will create the mysql command files necessary for building the tables
sub build_table_files {
    my $options=shift;
    my $mysql_dir=shift;
    my $prefix=$options->get_prefix();

    # Description of the tables to be created
    my %tables=('_start' => "DROP TABLE IF EXISTS ${prefix}_start;
CREATE TABLE ${prefix}_start (
    filename varchar(100),
    presence varchar(20)
    );",
		'_summaries' => "DROP TABLE IF EXISTS ${prefix}_summaries;
CREATE TABLE ${prefix}_summaries (
    summary varchar(120)
    );",
		'_completion_status' => "DROP TABLE IF EXISTS ${prefix}_completion_status;
CREATE TABLE ${prefix}_completion_status (
    step varchar(100),
    status varchar(20)
    );",
		'_inclusion_correlation' => "DROP TABLE IF EXISTS ${prefix}_inclusion_correlation;
CREATE TABLE ${prefix}_inclusion_correlation (
    LaneName1 varchar(50) NOT NULL,
    LaneName2 varchar(50) NOT NULL,
    Lane1Inc smallint unsigned NOT NULL,
    Lane2Inc smallint unsigned NOT NULL
    );",
		'_dataset' => "DROP TABLE IF EXISTS ${prefix}_dataset;
CREATE TABLE ${prefix}_dataset (
    lane_idx varchar(50) NOT NULL,
    lane_id varchar(50) NOT NULL,
    pair_id varchar(50) NOT NULL,
    sample_id varchar(50) NOT NULL,
    INDEX index_lane (lane_idx)
    );",
		'_fusion_transcripts' => "DROP TABLE IF EXISTS ${prefix}_fusion_transcripts;
CREATE TABLE ${prefix}_fusion_transcripts (
    readId varchar(50) not null,
    readlength mediumint unsigned not null,
    type varchar(10) not null,
    trans1mapstart varchar(20) not null,
    trans1 varchar(50) not null,
    chr1 varchar(20) not null,
    genomic1mapstart varchar(20) not null,
    gene1 varchar(50) not null,
    gene1rpkm float unsigned not null,
    trans2mapstart varchar(20) not null,
    trans2 varchar(50) not null,
    chr2 varchar(20) not null,
    genomic2mapstart varchar(20) not null,
    gene2 varchar(50) not null,
    gene2rpkm float unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_trans1 (trans1),
    INDEX idx_trans2 (trans2),
    INDEX idx_gene1 (gene1),
    INDEX idx_gene2 (gene2),
    INDEX idx_lane (LaneName)
    );",
		'_fusion_transcripts_support' => "DROP TABLE IF EXISTS ${prefix}_fusion_transcripts_support;
CREATE TABLE ${prefix}_fusion_transcripts_support (
    type varchar(10),
    location varchar(100),
    number mediumint unsigned not null,
    start_chr varchar(50) NULL,
    start int unsigned NULL,
    end_chr varchar(50) NULL,
    end int unsigned NULL,
    reference varchar(10) NOT NULL,
    LaneName varchar(20) NOT NULL,
    INDEX idx_location (location),
    INDEX idx_lane (LaneName)
    );",
		'_read_classification' => "DROP TABLE IF EXISTS ${prefix}_read_classification;
CREATE TABLE ${prefix}_read_classification (
       lane_id varchar(50) NOT NULL,
       Exonic double unsigned NOT NULL,
       ExFrac double unsigned NOT NULL,
       Intronic double unsigned NOT NULL,
       IntrFrac double unsigned NOT NULL,
       Intergenic double unsigned NOT  NULL,
       InterFrac double unsigned NOT  NULL,
       Total double unsigned NOT NULL,
       TptFrac double unsigned NOT NULL,
       index idx_lane (lane_id)
);",
		'_exon_inclusion_reads' => "DROP TABLE IF EXISTS ${prefix}_exon_inclusion_reads;
CREATE TABLE ${prefix}_exon_inclusion_reads (
       exon_id varchar(50) NOT NULL,
       chr varchar(50) NOT NULL,
       start int unsigned NOT NULL,
       end int unsigned NOT NULL,
       ExIncl int unsigned NOT NULL,
       JuncInc int unsigned NOT NULL,
       JuncExc int unsigned NOT  NULL,
       inc_rate double unsigned NULL,
       sample_id varchar(100) NOT NULL,
       index idx_exon (exon_id),
       index idx_sample (sample_id)
);",
		'_exon_inclusion_pooled' => "DROP TABLE IF EXISTS ${prefix}_exon_inclusion_pooled;
CREATE TABLE ${prefix}_exon_inclusion_pooled (
       gene_id varchar(50) NOT NULL,
       exon_id varchar(50) NOT NULL,
       ExIncl double unsigned NOT NULL,
       JuncInc double unsigned NOT NULL,
       JuncExc double unsigned NOT  NULL,
       inc_rate double unsigned NULL,
       sample varchar(50) NOT NULL,
       index idx_gene (gene_id),
       index idx_exon (exon_id),
       index idx_sample (sample)
);",
		'_inclusion_dist' => "DROP TABLE IF EXISTS ${prefix}_inclusion_dist;
CREATE TABLE ${prefix}_inclusion_dist (
    support mediumint unsigned not null,
    incl_percent smallint unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_lane (LaneName)
);",
		'_transcript_expression_levels_pooled' => "DROP TABLE IF EXISTS ${prefix}_transcript_expression_levels_pooled;
CREATE TABLE ${prefix}_transcript_expression_levels_pooled (
       gene_id varchar(50) NOT NULL,
       transcript_id varchar(50) NOT NULL,
       flux_locus varchar(50) NOT NULL,
       rpkm double unsigned NULL,
       sample varchar(100) NOT NULL,
       index idx_gene (gene_id),
       index idx_transcript (transcript_id),
       index idx_flux_locus (flux_locus),
       index idx_lane (sample)
);",
		'_junction_classification' => "DROP TABLE IF EXISTS ${prefix}_junction_classification;
CREATE TABLE ${prefix}_junction_classification (
       table_id varchar(100) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_exon_classification' => "DROP TABLE IF EXISTS ${prefix}_exon_classification;
CREATE TABLE ${prefix}_exon_classification (
       table_id varchar(100) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_store_reads' => "DROP TABLE IF EXISTS ${prefix}_store_reads;
CREATE TABLE ${prefix}_store_reads (
    filename varchar(100),
    uncompressed varchar(100) not null,
    md5uncomp char(32) not null,
    compressed varchar(100) not null,
    md5comp char(32) not null,
    md5global char(32) not null
    );",
		'_gene_readcount_pooled' => "DROP TABLE IF EXISTS ${prefix}_gene_readcount_pooled;
CREATE TABLE ${prefix}_gene_readcount_pooled (
    gene_id varchar(100) not null,
    readcount int unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_gene (gene_id),
    INDEX idx_lane (LaneName)
    );",
		'_gene_RPKM_pooled' => "DROP TABLE IF EXISTS ${prefix}_gene_RPKM_pooled;
CREATE TABLE ${prefix}_gene_RPKM_pooled (
    gene_id varchar(100) not null,
    RPKM double unsigned not null,
    sample varchar(50) not null,
    INDEX idx_gene (gene_id),
    INDEX idx_sample (sample)
    );",
		'_gene_RPKM_pooled_flux' => "DROP TABLE IF EXISTS ${prefix}_gene_RPKM_pooled_flux;
CREATE TABLE ${prefix}_gene_RPKM_pooled_flux (
    gene_id varchar(100) not null,
    RPKM double unsigned not null,
    sample varchar(50) not null,
    INDEX idx_gene (gene_id),
    INDEX idx_sample (sample)
    );",
		'_gene_RPKM_dist' => "DROP TABLE IF EXISTS ${prefix}_gene_RPKM_dist;
CREATE TABLE ${prefix}_gene_RPKM_dist (
    support varchar(100) not null,
    RPKM mediumint unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_lane (LaneName)
    );",
		'_exon_RPKM_pooled' => "DROP TABLE IF EXISTS ${prefix}_exon_RPKM_pooled;
CREATE TABLE ${prefix}_exon_RPKM_pooled (
    exon_id varchar(100) not null,
    RPKM double unsigned not null,
    sample varchar(50) not null,
    INDEX idx_exon (exon_id),
    INDEX idx_sample (sample)
    );",
		'_EJEI' => "DROP TABLE IF EXISTS ${prefix}_EJEI;
CREATE TABLE ${prefix}_EJEI (
    gene_id varchar(100) not null,
    junction_id varchar(100) not null,
    EJEI double unsigned not null,
    total int not null,
    LaneName varchar(50) not null,
    INDEX idx_gene (gene_id),
    INDEX idx_junction (junction_id),
    INDEX idx_lane (LaneName)
    );",
		'_merged_SAM' => "DROP TABLE IF EXISTS ${prefix}_merged_SAM;
CREATE TABLE ${prefix}_merged_SAM (
    LaneName varchar(100) not null,
    BAM_file varchar(100) not null,
    INDEX idx_lane (LaneName)
    );",

		'_junction_coverage' => "DROP TABLE IF EXISTS ${prefix}_junction_coverage;
CREATE TABLE ${prefix}_junction_coverage (
    filename varchar(100),
    features mediumint unsigned not null
    );",
		'_splicing_summary' => "DROP TABLE IF EXISTS ${prefix}_splicing_summary;
CREATE TABLE ${prefix}_splicing_summary (
    junc_type varchar(50) not null,
    detected mediumint unsigned not null,
    total mediumint unsigned NULL
    );",
		'_novel_junctions_summary' => "DROP TABLE IF EXISTS ${prefix}_novel_junctions_summary;
CREATE TABLE ${prefix}_novel_junctions_summary (
    junction_id varchar(50) not null,
    chr varchar(50) not null,
    start int unsigned not null,
    end int unsigned not null,
    junc_type varchar(50) not null,
    support int unsigned not null,
    sample varchar(50) not null
    );",
		'_all_junctions_class' => "DROP TABLE IF EXISTS ${prefix}_all_junctions_class;
CREATE TABLE ${prefix}_all_junctions_class (
    chr1 varchar(50) not null,
    start int unsigned not null,
    chr2 varchar(50) not null,
    end int unsigned not null,
    junc_type varchar(50) not null,
    support int unsigned not null,
    exons1 text,
    exons2 text,
    LaneName varchar(50) not null
    );",
		'_all_junctions_class_pooled' => "DROP TABLE IF EXISTS ${prefix}_all_junctions_class_pooled;
CREATE TABLE ${prefix}_all_junctions_class_pooled (
    chr1 varchar(50) not null,
    start int unsigned not null,
    chr2 varchar(50) not null,
    end int unsigned not null,
    junc_type varchar(50) not null,
    support int unsigned not null,
    exons1 text,
    exons2 text,
    sample varchar(50) not null
    );",
		'_exon_coverage' => "DROP TABLE IF EXISTS ${prefix}_exon_coverage;
CREATE TABLE ${prefix}_exon_coverage (
    filename varchar(100),
    features mediumint unsigned not null
    );",
		'_gene_coverage' => "DROP TABLE IF EXISTS ${prefix}_gene_coverage;
CREATE TABLE ${prefix}_gene_coverage (
    filename varchar(100),
    features mediumint unsigned not null
    );",
		'_bed_files' => "DROP TABLE IF EXISTS ${prefix}_bed_files;
CREATE TABLE ${prefix}_bed_files (
    filename varchar(100),
    features int unsigned not null
    );",
		'_proj_coverage' => "DROP TABLE IF EXISTS ${prefix}_proj_coverage;
CREATE TABLE ${prefix}_proj_coverage (
    filename varchar(100),
    features mediumint unsigned not null
    );",
		'_extract_annotation' => "DROP TABLE IF EXISTS ${prefix}_extract_annotation;
CREATE TABLE ${prefix}_extract_annotation (
    filename varchar(40),
    features mediumint unsigned not null
    );",
		'_split_mapping_breakdown' => "DROP TABLE IF EXISTS ${prefix}_split_mapping_breakdown;
CREATE TABLE ${prefix}_split_mapping_breakdown (
    filename varchar(100),
    type varchar(10),
    location varchar(100),
    number mediumint unsigned not null,
    start_chr varchar(50) NULL,
    start int unsigned NULL,
    end_chr varchar(50) NULL,
    end int unsigned NULL,
    INDEX idx_location (location)
    );",
		'_initial_clusters' => "DROP TABLE IF EXISTS ${prefix}_initial_clusters;
CREATE TABLE ${prefix}_initial_clusters (
    ClustertType varchar(100),
    chromosome varchar(40),
    ClusterNumber  mediumint unsigned not null
    );",
		'_unique_maps_split' => "DROP TABLE IF EXISTS ${prefix}_unique_maps_split;
CREATE TABLE ${prefix}_unique_maps_split (
    filename varchar(100),
    type varchar(30),
    number mediumint unsigned not null
    );",
		'_unique_maps_junctions' => "DROP TABLE IF EXISTS ${prefix}_unique_maps_junctions;
CREATE TABLE ${prefix}_unique_maps_junctions (
    filename varchar(100),
    type varchar(10),
    number mediumint unsigned not null
    );",
		'_read_dist_transcripts' => "DROP TABLE IF EXISTS ${prefix}_read_dist_transcripts;
CREATE TABLE ${prefix}_read_dist_transcripts (
    length_cat varchar(50) not null,
    start smallint unsigned not null,
    position tinyint unsigned not null,
    hits mediumint unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_lane (LaneName),
    INDEX idx_length (length_cat)
    );",
		'_unique_maps_transcripts' => "DROP TABLE IF EXISTS ${prefix}_unique_maps_transcripts;
CREATE TABLE ${prefix}_unique_maps_transcripts (
    filename varchar(100),
    transcript_id varchar(60),
    uniqueMaps mediumint unsigned NOT NULL,
    INDEX idx_trans_id (transcript_id)
    );",
		'_unique_maps_genome' => "DROP TABLE IF EXISTS ${prefix}_unique_maps_genome;
CREATE TABLE ${prefix}_unique_maps_genome (
    filename varchar(100),
    chromosome varchar(40),
    uniqueMaps mediumint unsigned not null
    );",
		'_junction_maps_class' => "DROP TABLE IF EXISTS ${prefix}_junction_maps_class;
CREATE TABLE ${prefix}_junction_maps_class (
    junction_id varchar(100),
    type varchar(20),
    support mediumint unsigned not null,
    LaneName varchar(50) not null,
    INDEX idx_junc (junction_id),
    INDEX idx_type (type),
    INDEX idx_lane (LaneName)
    );",
		'_split_mapping' => "DROP TABLE IF EXISTS ${prefix}_split_mapping;
CREATE TABLE ${prefix}_split_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_unmapped_reads' => "DROP TABLE IF EXISTS ${prefix}_unmapped_reads;
CREATE TABLE ${prefix}_unmapped_reads (
    filename varchar(100),
    unmappedReads mediumint unsigned not null
    );",
		'_merged_mapping' => "DROP TABLE IF EXISTS ${prefix}_merged_mapping;
CREATE TABLE ${prefix}_merged_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_junctions_mapping' => "DROP TABLE IF EXISTS ${prefix}_junctions_mapping;
CREATE TABLE ${prefix}_junctions_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_transcriptome_mapping' => "DROP TABLE IF EXISTS ${prefix}_transcriptome_mapping;
CREATE TABLE ${prefix}_transcriptome_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_genome_mapping' => "DROP TABLE IF EXISTS ${prefix}_genome_mapping;
CREATE TABLE ${prefix}_genome_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_recursive_mapping' => "DROP TABLE IF EXISTS ${prefix}_recursive_mapping;
CREATE TABLE ${prefix}_recursive_mapping (
    filename varchar(100),
    totalReads int unsigned not null,
    mappedReads int unsigned not null,
    uniqueReads int unsigned not null,
    100uniqueReads int unsigned not null,
    LaneName varchar(50) not null
    );",
		'_indices' => "DROP TABLE IF EXISTS ${prefix}_indices;
CREATE TABLE ${prefix}_indices (
    filename varchar(200) NOT NULL,
    type varchar(20) NOT NULL,
    created varchar(20) NOT NULL
    );",
		'_transcript_seqs' => "DROP TABLE IF EXISTS ${prefix}_transcript_seqs;
CREATE TABLE ${prefix}_transcript_seqs (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_junction_seqs' => "DROP TABLE IF EXISTS ${prefix}_junction_seqs;
CREATE TABLE ${prefix}_junction_seqs (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_genes' => "DROP TABLE IF EXISTS ${prefix}_genes;
CREATE TABLE ${prefix}_genes (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_detected_genes' => "DROP TABLE IF EXISTS ${prefix}_detected_genes;
CREATE TABLE ${prefix}_detected_genes (
       type varchar(200) NOT NULL,
       reliability varchar(10) NOT NULL,
       sample varchar(50) NOT NULL,
       detected int NOT NULL
);",
		'_transcripts' => "DROP TABLE IF EXISTS ${prefix}_transcripts;
CREATE TABLE ${prefix}_transcripts (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_detected_transcripts' => "DROP TABLE IF EXISTS ${prefix}_detected_transcripts;
CREATE TABLE ${prefix}_detected_transcripts (
       type varchar(200) NOT NULL,
       reliability varchar(10) NOT NULL,
       sample varchar(50) NOT NULL,
       detected int NOT NULL
);",
		'_exon_seqs' => "DROP TABLE IF EXISTS ${prefix}_exon_seqs;
CREATE TABLE ${prefix}_exon_seqs (
       table_id varchar(200) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_junctions' => "DROP TABLE IF EXISTS ${prefix}_junctions;
CREATE TABLE ${prefix}_junctions (
       table_id varchar(100) NOT NULL,
       creation varchar(20) NOT NULL
);",
		'_qualitiespos' => "DROP TABLE IF EXISTS ${prefix}_qualitiespos;
CREATE TABLE ${prefix}_qualitiespos (
       LaneName varchar(100) NOT NULL,
       position smallint UNSIGNED NOT NULL,
       mean double UNSIGNED NOT NULL
);",
		'_ambiguous' => "DROP TABLE IF EXISTS ${prefix}_ambiguous;
CREATE TABLE ${prefix}_ambiguous (
       LaneName varchar(100) NOT NULL,
       position smallint UNSIGNED NOT NULL,
       ambiguous int UNSIGNED NOT NULL,
       A int UNSIGNED NOT NULL,
       C int UNSIGNED NOT NULL,
       G int UNSIGNED NOT NULL,
       T int UNSIGNED NOT NULL
);",
		'_read_stats' => "DROP TABLE IF EXISTS ${prefix}_read_stats;
CREATE TABLE ${prefix}_read_stats (
       Filename varchar(100)NOT NULL,
       Species varchar(50) NOT NULL,
       Project varchar(100) NOT NULL,
       ReadLength SMALLINT UNSIGNED NOT NULL,
       TotalReads INT UNSIGNED NOT NULL,
       NoAmbiguousBases INT unsigned NOT NULL,
       AmbiguousBases INT unsigned NOT NULL,
       UniqueReads INT unsigned NOT NULL,
       LaneName varchar(50) NOT NULL,
       pair_id varchar(50) NOT NULL
);",
		'_register_results' => "DROP TABLE IF EXISTS ${prefix}_register_results;
CREATE TABLE ${prefix}_register_results (
    filetype varchar(100),
    filename varchar(100),
    filemd5sum char(32) not null,
    registeredfile varchar(200),
    runmd5sum char(32) not null
    );",
	);

    my $working_dir=$options->get_projectdir();
    unless ($mysql_dir) {
	$mysql_dir=$working_dir.'/mysql/table_build';
    }
    chdir($mysql_dir) ||	
	die "Could not chdir to $mysql_dir from $working_dir\n";
    foreach my $table (sort keys %tables) {	
	my $file_name=$prefix.$table.'.sql';
	my $table_fh=get_fh($file_name,1);
	print $table_fh $tables{$table},"\n";
	close($table_fh);
    }
    chdir($working_dir) ||
	die "Could not chdir from $mysql_dir to $working_dir\n";
    return(keys %tables);
}

# Create the directory structure necessary for the running of the pipeline
sub create_directory_structure {
    my $log_fh=shift;
    my $options=shift;
    my $table_fh=shift;

    my $prefix=$options->get_prefix();

    # First check the directory structure and create any directories that are
    # not there
    $log_fh->printlog("Checking directory structure and attempting creation of missing dirs...");

    my $localdir=$options->get_localdir();
    $log_fh->printlog("Setting localdir as: $localdir");
    my @directories=('mysql',
		     'mysql/table_build',
		     'mysql/table_data',
#		     'bin',
		     'GEMIndices',
		     'genome',
		     'sequence',
		     'transcriptome',
		     'junctions',
		     'splitmapping',
		     'recursivemap',
		     'exons',
		     'readData',
		     'trimmedReads',
		     'graphs',
		     'logs',
		     "$localdir",
		     'clusters',
		     'exclusion',
		     'results',
		     'SAM',
		     'work');
    # The -L option is necessary to print the same path for all nodes
    my $project_dir=$options->get_projectdir();
    $log_fh->printlog("Using dir: ".$project_dir);

    my %vars=('COMMONDB' => $options->get_commondb(),
	      'DB'=> $options->get_database(),
	      'PROJECT' => $project_dir,
	      'TABLES' => 'mysql/table_build',
	      'TAB_DAT' => 'mysql/table_data',
	      'BIN' => 'bin',
	      'LOGS' => 'logs',
	      'PAIRED' => $options->get_paired(),
	      'READDIR' => 'readData',
	      'SEQUENCEDIR' => 'sequence', # (Raw data if not fa or fastq)
	      'GEMINDICES' => 'GEMIndices',
	      'SPLITMAPDIR' => 'splitmapping',
	      'RECMAPDIR' => 'recursivemap',
	      'GRAPHS' => 'graphs',
	      'SPECIES' => $options->get_species(),
	      'PROJECTID' => $options->get_project(),
	      'EXPID' => $options->get_experiment(),
	      'PREFIX' => $prefix,
	      'ANNOTATION' => (abs_path($options->get_annotation()) || die "SORRY, I can't find the specified annotation file, please check the path"),
	      'GENOMESEQ' => (abs_path($options->get_genome()) || die "SORRY, I can't find the specified genome file, please check the path\n"),
	      'TRANSDIR' => 'transcriptome',
	      'EXONDIR' => 'exons',
	      'JUNCTIONSDIR' => 'junctions',
	      'MISMATCHES' => $options->get_mismatches(),
	      'LOCALDIR' => $options->get_localdir(),
	      'GENOMEINDEX' => $options->get_genomeindex(),
	      'TRANSCRIPTOMEINDEX' => $options->get_transcriptomeindex(),
	      'JUNCTIONSINDEX' => $options->get_junctionsindex(),
	      'EXCLUSIONFILE' => $options->get_exclusionfile(),
	      'EXONSFASTA' => $options->get_exonsfasta(),
	      'JUNCTIONSFASTA' => $options->get_junctionsfasta(),
	      'TRANSCRIPTOMEFASTA' => $options->get_transcriptomefasta(),
	      'GENOMEDIR' => 'genome',
	      'STRANDED' => $options->get_stranded(),
	      'READLENGTH' => $options->get_readlength(),
	      'GENECLASSTABLE' => $options->get_geneclass(),
	      'TRANSCLASSTABLE' => $options->get_transclass(),
	      'JUNCTIONSTABLE' => $options->get_junctionstable(),
	      'JUNCTIONSCLASSTABLE' => $options->get_junctionsclass(),
	      'EXONSCLASSTABLE' => $options->get_exonsclass(),
	      'FILELIST' => $options->get_files(),
	      'THREADS' => $options->get_threads(),
	      'MAPPER' => $options->get_mapper(),
	      'CLUSTERDIR' => 'clusters',
	      'CLUSTER' => $options->get_cluster(),
	      'SAMDIR' => 'SAM',
	      'LOCALPARALLEL' => 'work',
	      'QUALITIES' =>$options->get_qualities(),
	      'PREPROCESS' => $options->get_preprocess(),
	      'RESULTS' => 'results',
	      'FLUXMEM' => $options->get_fluxmem(),
	      'TRIMLENGTH' => $options->get_trimlength(),
	      'MAXINTRONLENGTH' => $options->get_maxintronlength()
);

    foreach my $dir (@directories) {
	$log_fh->printlog("Checking $dir...");
	check_dir($dir,
		  $log_fh,
		  $table_fh);
	my $command="chmod g+x $dir";
	run_system_command($command,
			   $log_fh);
    }
    $log_fh->printlog("done");
    return(\%vars);
}

sub check_dir {
    my $dir=shift;
    my $log_fh=shift;
    my $table_fh=shift;
    if (-e $dir) {
	$log_fh->printlog("$dir already exists");
	$log_fh->printlog("Skipping");
    } elsif ($dir=~/^bin/) {
	warn "THIS SHOULD NOT HAPPEN\n";
    } else {
	mkdir($dir) &&
	    $log_fh->printlog("Creating dir $dir");
    }

    print $table_fh join ("\t",
			  $dir,
			  'present'),"\n";
    return();
}

## Print the necessary files for running the pipeline:
# These are the configuration file
sub print_config_file {
    my $vars=shift;
    my $log_fh=shift;
    my $conf_file='project_config.txt';
    my $conf_fh=get_fh($conf_file,1);

    $log_fh->printlog("Creating the configuration file at $conf_file...");

    # Print the keys sorted as it makes it easier to check the file if necessary
    foreach my $var (sort keys %{$vars}) {
	my $value=$vars->{$var};
	if (defined $value) {
	    print $conf_fh join("\t",
				$var,
				$value),"\n";
	} else {
	    die "No value found for $var. Please check\n";
	}
    }
    close($conf_fh);
    $log_fh->printlog("done");
}

# Print the pipeline file
sub print_pipeline_file {
    my $tables=shift;
    my $vars=shift;
    my $template=shift;
    my $log_fh=shift;
    my $pipfile='RNAseq_pipeline.txt';

    $log_fh->printlog("Creating the master file at $pipfile...");

    my $template_fh=get_fh($template);
    my $new_pip_fh=get_fh($pipfile,1);

    while (my $line=<$template_fh>) {
	chomp($line);
	# Set the variable lines
	# Ignore comment lines
	if ($line!~/^\#/) {
	    # Set the variable lines
	    if ($line=~/=/) {
		my($variable,$assignment)=split(/\s*=\s*/,$line);
		my $value=$vars->{$variable};
		if ($assignment) {
		    $value=join('/',
				$vars->{'PROJECT'},
				$vars->{$variable});
		}
		$log_fh->printlog(join("\t",
				       $variable,
				       $value));
		$line=join(' = ',
			   $variable,
			   $value);
	    } elsif ($line=~/(^|\s+)(_\w+)/) {
		if (exists $tables->{$2}) {
		    # Table is present in the module. This has to be removed
		    # moving the tables gradually to the rules
		} elsif ($vars->{'PREFIX'}) {
		    # The table does not exist in the module, and has to be
		    # creted by the rule. When all the tables are migrated this
		    # warning should be removed
#		    warn "WARNING: $2 Does not exist as a table\n";
		    my $prefix=$vars->{'PREFIX'};
		    $tables->{$2}=$prefix.$2;
		} else {
		    die "Something whent very wrong...\n";
		}
		$line=~s/(^|\s+)(_\w+)/$1$tables->{$2}/g;
	    }
	}
	print $new_pip_fh $line,"\n";
    }
    $log_fh->printlog("done");
}

###
# Some subroutines to clean a directory and start from scratch
sub clear_tables {
    my $options=shift;
    my $database=$options->get_database();

    my @tables=@{get_tables_list($options)};

    foreach my $table (@tables) {
	my $command="mysql $database -e 'DROP TABLE IF EXISTS $table'";

	run_system_command($command);
    }
}

# Get a list of all the tables starting with a certain prefix
sub get_tables_list {
    my $options=shift;
    my $dbh=$options->get_dbh();
    my $prefix=$options->get_prefix() || die "No table prefix\n";

    my ($query,$sth);
    $query= "SHOW tables";
    $sth=$dbh->prepare($query);
    $sth->execute();

    my @tables;
    while (my ($table)=$sth->fetchrow_array()) {
	if ($table=~/^$prefix/) {
	    push @tables, $table;
	} else {
#	    print STDERR "Skipping $table\n";
	}
    }

    return(\@tables);
}

# These are subs that will clear the directories in the project directory
### TO DO make an additional level of clean so the genome directory and maybe
# some other ones are not completely wiped
sub clear_dirs {
    my $options=shift;

    # First check the directory structure and create any directories that are
    # not there
    print STDERR "Checking directory structure...\n";

    my $localdir=$options->get_localdir();
    my $project_dir=$options->get_projectdir();
    print STDERR "Setting localdir to: $localdir\n";
    my @directories=('mysql',
		     'genome',
		     'transcriptome',
		     'junctions',
		     'splitmapping',
		     'recursivemap',
		     'exons',
		     'trimmedReads',
		     'graphs',
		     'logs',
		     "$localdir",
		     'clusters',
		     'exclusion',
		     'SAM',
		     'work');
    # The -L option is necessary to print the same path for all nodes
    print STDERR "Using dir: ".$project_dir,"\n";

    foreach my $dir (@directories) {
	my $command="rm -r $dir";
	if (-r $dir) {
	    run_system_command($command);
	}
    }
}

# Delete some files that may be left in the directory
sub clear_files {
    my $opts=shift;
    my $prefix=$opts->get_prefix();

    my $command;

    # Remove any table files that are left over as well as log files the
    # pipeline file and other temporary files
    $command="ls $prefix* *.log RNAseq_pipeline.txt *.tmp.sequences.txt execute_pipeline.err exon.strange.list rep.exons.txt 2>/dev/null";
    my @list=`$command`;
    foreach my $file (@list) {
	chomp($file);
	print STDERR $file,"\n";
	if (-r $file) {
	    $command="rm $file";
	    run_system_command($command);
	}
    }
}

# Delete the information from the projects and experiments tables
sub clear_common_tables {
    my $opts=shift;

    my $proj_id=$opts->get_project();
    my $exp_id=$opts->get_experiment();
    my $dbh=$opts->get_commondbh();
    my $proj_tab='projects';
    my $exp_tab='experiments';

    my ($query,$sth,$count);

    # Delete the entry from the experiments table

    # Check if the table is present and create it if not
    my $present=check_table_existence($dbh,
				      $exp_tab);

    if ($present) {
	$query ="DELETE FROM $exp_tab ";
	$query.='WHERE project_id = ? AND experiment_id = ?';
	$sth=$dbh->prepare($query);
	print STDERR "Executing: $query\n";
	$sth->execute($proj_id,$exp_id);

	# Check if there are any entries left in the experiments table
	# corresponding to the project and if not remove it also from the
	# projects table
	$query ="SELECT experiment_id FROM $exp_tab ";
	$query.='WHERE project_id = ?';
	$sth=$dbh->prepare($query);
	$count=$sth->execute($proj_id);
	
	if ($count > 0) {
	    print STDERR $count,"\tExperiments left in $exp_tab corresponding to $proj_id\n";
	} else {
	    print STDERR "No experiments left in $exp_tab corresponding to $proj_id\n";
	    # Delete the entry from the projects table
	    $query ="DELETE FROM $proj_tab ";
	    $query.='WHERE project_id = ?';
	    $sth=$dbh->prepare($query);
	    $sth->execute($proj_id);
	}
    } else {
	print STDERR "No entry for $exp_id in $exp_tab\n";
    }
}

sub clear_experiment {
    my $options=shift;

    # We need to make sure we have values for the different options (so maybe
    # it would be good to add here something to read the config file or at
    # least to have)
    my $project=$options->get_project();
    my $experiment=$options->get_experiment();
    my $commondb=$options->get_commondb();
    my $database=$options->get_database();

    unless ($project && $experiment &&
	    $commondb && $database) {
	die "I need to know both the project and the experiment you want to remove\n";
    }

    # Print message if in clean mode
    print STDERR "WARNING: Clean mode.\n";
    print STDERR "This will remove the following:\n";
    print STDERR "\tAll tables for this experiment ($experiment) in $project\n";
    print STDERR "\tAll direactories that the script created originally with the files within (this includes all mappings), with the exception of the readData directory.\n";
    print STDERR "\tAll entries in the common tables corresponding to this experiment\n";
    print STDERR "Do you want to continue? (y/N)\n";
    my $reply=<STDIN>;
    chomp($reply);
    unless ($reply=~/y/i) {
	die "Aborting\n";
    }

    # Clear the database
    clear_tables($options);

    # Clear the directories
    clear_dirs($options);

    # Clear the entries from the projects and experiments tables in the commondb
    clear_common_tables($options);

    # Clear some files that may be left in the running dir
    clear_files($options);
}

1;
