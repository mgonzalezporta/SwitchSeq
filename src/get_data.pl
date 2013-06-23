#!/usr/bin/perl

use strict;
use warnings;
use LWP::Simple;
use Getopt::Long;
use Data::Dumper;

##############################################################

# add TableTools.min.js instead of TableTools.js
# add TableTools.css

# OPTIONS
my ($data_dir, $species, $ensembl_v);

GetOptions(
        'data_dir|d:s'			=> \$data_dir,
        'species|s:s'			=> \$species,
        'ensembl_v|e:i'			=> \$ensembl_v
);

##############################################################

# GET FILES
print "# Obtaining the necessary files\n";

# # data
# my $ensembl_url="http://www.ebi.ac.uk/~mar/tools/lorem/data/$species/_ensembl$ensembl_v.annot_coding.txt";
# &get_file($ensembl_url, "$data_dir/$species");

# my $prot_seq_url="http://www.ebi.ac.uk/~mar/tools/lorem/data/$species/_ensembl$ensembl_v.prot_seq.tar.gz";
# &get_file($prot_seq_url, "$data_dir/$species");

# my $appris_url="http://www.ebi.ac.uk/~mar/tools/lorem/data/$species/_appris.results.rel15.2May2013.v1.main.tsv";
# &get_file($appris_url, "$data_dir/$species");

# my $pdb_url="http://www.ebi.ac.uk/~mar/tools/lorem/data/$species/_UniPdbCov.txt";
# &get_file($pdb_url, "$data_dir/$species");

# # css + js
# my $css_url="http://www.ebi.ac.uk/~mar/tools/lorem/data/css.tar.gz";
# &get_file($css_url, $data_dir);

# my $js_url="http://www.ebi.ac.uk/~mar/tools/lorem/data/js.tar.gz";
# &get_file($js_url, $data_dir);

##############################################################

sub get_file {
	my $url=$_[0];
	my $dir=$_[1];
	my $file=$dir."/".( split(/\//, $url) )[-1];

	getstore($url, $file);
	if ($file =~ /.tar.gz$/) {
		my $fname=( split(/\//, $file) )[-1];
		system("cd $dir;
				tar xzf $fname;
				rm $fname");
	}
}