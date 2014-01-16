#!/usr/bin/perl

use strict;
use warnings;

## args
my $species=$ARGV[1];
my $outdir=$ARGV[2];

## ensembl api
my $ensembl_v;
my $modules;
BEGIN { 
$ensembl_v=$ARGV[0];
$modules="/homes/mar/system/ensembl.$ensembl_v/modules"; 
}

use lib "/homes/mar/system/bioperl-live";
use lib "$modules";
use Bio::EnsEMBL::Registry;

## IO
my $fa_dir="$outdir/prot_seq/";
unless ( -e $fa_dir ) { 
	system("mkdir $fa_dir");
}

## set up connection
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', 	# alternatively 'useastdb.ensembl.org'
    -user => 'anonymous',
    -db_version => $ensembl_v
);

my %species_full=(
	hsa => "Human",
	mmu => "Mouse",
	rno => "Rat",		# not tested
	dre => "Zebrafish"	# not tested
);
my $slice_adaptor=$registry->get_adaptor( "$species_full{$species}", 'Core', 'Slice' );

## retrieve protein sequences for all protein coding transcripts
print "\t# Retrieving protein sequence info...\n";
my $slices = $slice_adaptor->fetch_all('toplevel');
my $gCount;
foreach my $slice ( @{ $slices } ) {
	my $genes = $slice->get_all_Genes();
	foreach my $gene ( @{$genes} ) {
		my $gId=$gene->stable_id();
		my $gBiotype=$gene->biotype();
		if ($gBiotype eq "protein_coding") {
			$gCount++;
			my @transcripts=@{ $gene->get_all_Transcripts };
			foreach my $transcript (@transcripts) {
				my $tId=$transcript->stable_id();
				my $tBiotype=$transcript->biotype();

    				if ($tBiotype eq "protein_coding") {
					my ($species_id, $numeric_id) = $tId =~ /([a-zA-Z]+)(\d+)/; 
    					my $fa_subdir=$fa_dir.$species_id.substr($numeric_id, 0, 8);
    					unless (-e $fa_subdir) { system("mkdir $fa_subdir") };
        				my $out_fa="$fa_subdir/$tId.fa";
        				my $seq=$transcript->translate()->seq();
        				
        				open(OUT_FA, ">$out_fa") or die "Cannot open $out_fa: $!";
					my $id=">".$gId."_".$tId;
        				print OUT_FA "$id\n";
					while (my $line = substr($seq, 0, 80, "")) {
						print OUT_FA "$line\n";
					}
        				close(OUT_FA);
    				}
			}
		}
	}
}
print "\t# Data retrieved for $gCount protein coding genes\n";
