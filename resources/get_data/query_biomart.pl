#!/usr/bin/perl

use strict;
use LWP::UserAgent;

my $in_xml=$ARGV[0];
my $species=$ARGV[1];
my $path=$ARGV[2];

my %species_full=(
	hsa	=> "hsapiens_gene_ensembl",
	mmu	=> "mmusculus_gene_ensembl",
	rno	=> "rnorvegicus_gene_ensembl",
	dre	=> "drerio_gene_ensembl"
);

open (FH,$in_xml) || die $!;
my $xml;
while (<FH>){
	if ($_ =~ /^<Dataset/) {
		$xml .= "<Dataset name = \"$species_full{$species}\" interface = \"default\" >\n";
	} else {
	        $xml .= $_;
	}
}
close(FH);

my $request = HTTP::Request->new("POST",$path,HTTP::Headers->new(),'query='.$xml."\n");
my $ua = LWP::UserAgent->new;

my $response;

$ua->request($request, 
	     sub{   
		 my($data, $response) = @_;
		 if ($response->is_success) {
		     print "$data";
		 }
		 else {
		     warn ("Problems with the web server: ".$response->status_line);
		 }
	     },1000);
