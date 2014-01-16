#!/usr/bin/perl

package GetData;
use strict;
use warnings;
use Exporter;
use LWP::Simple;

our @ISA= qw( Exporter );
our @EXPORT = qw( get_data );

sub get_data {	
	## collect arguments
	my $ref_arguments=$_[0];
	my $data_dir=$ref_arguments->{'data_dir'};
	my $species=$ref_arguments->{'species'};
	my $ensembl_v=$ref_arguments->{'ensembl_v'};

	## expand species name
	my %species_full=(
		hsa	    => "homo_sapiens",
		mmu	    => "mus_musculus",
		rno	    => "rattus_norvegicus",
		dre	    => "danio_rerio"
	);

	## define urls
	## by default, work with the latest version of APPRIS
	my %url=(
		ensembl_1   => "http://www.ebi.ac.uk/~mar/tools/lorem/$species.$ensembl_v/ensembl1.txt",
		ensembl_2   => "http://www.ebi.ac.uk/~mar/tools/lorem/$species.$ensembl_v/ensembl2.txt",
		prot_seq    => "http://www.ebi.ac.uk/~mar/tools/lorem/$species.$ensembl_v/prot_seq.tar.gz",
		appris      => "http://appris.bioinfo.cnio.es/download/data/$species_full{$species}/appris_data.principal.txt",
	);

	## create directory structure
	unless ( -e "$data_dir" ) { system("mkdir $data_dir") };
	unless ( -e "$data_dir/$species" ) { system("mkdir $data_dir/$species") };
	
	## obtain files
	print "# Obtaining the necessary files...\n";
	_get_file($url{'ensembl_1'}, "$data_dir/$species");
	_get_file($url{'ensembl_2'}, "$data_dir/$species");
	_get_file($url{'prot_seq'}, "$data_dir/$species");
	_get_file($url{'appris'}, "$data_dir/$species");
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
