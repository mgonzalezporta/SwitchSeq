#!/usr/bin/perl

package GetData;
use strict;
use warnings;
use Exporter;
use LWP::Simple;
use Getopt::Long;
use Data::Dumper;

our @ISA= qw( Exporter );
our @EXPORT = qw( get_data );

sub get_data {	
	## collect arguments
	my $ref_arguments=$_[0];
	my %arguments=%{$ref_arguments};
	my $data_dir=$arguments{'data_dir'};
	my $species=$arguments{'species'};
	my $ensembl_v=$arguments{'ensembl_v'};

	## define urls
	my %url=(
		ensembl 	=> "http://www.ebi.ac.uk/~mar/tools/lorem/$species/_ensembl$ensembl_v.annot_coding.txt",
		prot_seq    => "http://www.ebi.ac.uk/~mar/tools/lorem/$species/_ensembl$ensembl_v.prot_seq.tar.gz",
		appris      => "http://www.ebi.ac.uk/~mar/tools/lorem/$species/_appris.results.rel15.9Jun2013.v2.main.tsv",
		css         => "http://www.ebi.ac.uk/~mar/tools/lorem/css.tar.gz",
		js          => "http://www.ebi.ac.uk/~mar/tools/lorem/js.tar.gz"
	);

	## create directory structure
	unless ( -e "$data_dir" ) { system("mkdir $data_dir") };
	unless ( -e "$data_dir/$species" ) { system("mkdir $data_dir/$species") };
	
	## download files
	print "# Obtaining the necessary files...\n";
	&_get_file($url{'ensembl'}, "$data_dir/$species");
	&_get_file($url{'prot_seq'}, "$data_dir/$species");
	&_get_file($url{'appris'}, "$data_dir/$species");
	&_get_file($url{'css'}, $data_dir);
	&_get_file($url{'js'}, $data_dir);
	print "# Data files saved under $data_dir\n";
}

sub _get_file {
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

1;