#!/usr/bin/perl

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use List::Util qw[ max ];
use Getopt::Long;
use Data::Dumper;
use lib './src/';
use myFunctions;

##############################################################

# OPTIONS

my ($data_dir, $species, $ensembl_v, 
	$input, $out_dir, $cond1, $cond2, $threshold_gexp, $plot);

GetOptions(
        'data_dir|d:s'			=> \$data_dir,
        'species|s:s'			=> \$species,
        'ensembl_v|e:i'			=> \$ensembl_v,
        'input|i=s'				=> \$input,
        'out_dir|o=s'			=> \$out_dir,
        'cond1|c1=s'			=> \$cond1,
        'cond2|c2=s'			=> \$cond2,
        'threshold_gexp|g=f'	=> \$threshold_gexp,
        'plot|p=s'				=> \$plot
);

# adjust columns
my @cond1=&adjust_columns($cond1);
my @cond2=&adjust_columns($cond2);

my $out_file="$out_dir/switch.txt";

##############################################################
# check input file: number of pc genes
##############################################################

print "# Obtaining switch events and generating output\n";

# Obtain the major transcript in each sample
my %major_tx;
	# gId => | major_tx_id  => [@columns]
	#	     | major_tx_exp => [@columns]
	#        | gExp         => [@columns]
open (INPUT, "< $input") or die "Could not open $input: $!";
while( my $row = <INPUT>)  {

	next if $.==1;

	chomp ($row);
	my @row=split(/ /, $row);
	my $gId=$row[0];
	my $tId=$row[1];
	my $end=$#row-2;
	my $length=$end+1;

	# initialise vectors if it's the first time you see the gene
	if (!defined $major_tx{$gId}) {
		@{ $major_tx{$gId}{'major_tx_id'} }=('NA') x $length;
		@{ $major_tx{$gId}{'major_tx_exp'} }=(0) x $length;
		@{ $major_tx{$gId}{'gExp'} }=(0) x $length;
	}
	
	for my $i (0..$end) {
		@{ $major_tx{$gId}{'gExp'} }[$i]+=$row[$i+2];
		if ($major_tx{$gId}{'major_tx_exp'}[$i] < $row[$i+2]) {
			@{ $major_tx{$gId}{'major_tx_exp'} }[$i]=$row[$i+2];
			@{ $major_tx{$gId}{'major_tx_id'} }[$i]=$tId;
		}
	}
}
close (INPUT);
# TO DO: what if they are the same?
# /homes/mar/home_microarray/workspace/cagekid/scripts/plots/mar_v5/fig3.switch_v2.R

#print Dumper %major_tx;

##############################################################

# Obtain the most recurrent major transcripts for each condition
my %recurrent_tx;
	# gId => | cond1 => | recurrent_tx_id
	#                   | recurrent_tx_count
	# 		 | cond2 => | recurrent_tx_id
	# 		            | recurrent_tx_count
foreach my $gId (keys %major_tx) {
	my %subset_major_tx=%{ $major_tx{$gId} };
	# @cond1 and @cond2 are defined by the user as arguments (see above)

	$recurrent_tx{$gId}{'cond1'}=&get_most_recurrent_tx(\%subset_major_tx, \@cond1);
	$recurrent_tx{$gId}{'cond2'}=&get_most_recurrent_tx(\%subset_major_tx, \@cond2);
}


# clean up
# if any of the conditions doesn't have any transcript expressed, discard the whole gene
foreach my $gId (keys %recurrent_tx) {
	if (! defined($recurrent_tx{$gId}{'cond1'}{'recurrent_tx_id'}) || 
		! defined($recurrent_tx{$gId}{'cond2'}{'recurrent_tx_id'})) {
		delete $recurrent_tx{$gId};
	}
}

# my $n=(keys %recurrent_tx);
#print Dumper %recurrent_tx;
##############################################################

# Find and annotate switch events

# Load data from files
my $ensembl_input="$data_dir/$species/_ensembl$ensembl_v.annot_coding.txt";
my %ensembl;
open (INPUT, "< $ensembl_input") or die "Could not open $ensembl_input: $!";
while( my $row = <INPUT>)  {
	chomp ($row);

	if ($row =~ /^ENSG/) {
		my @row=split(/\s+/, $row);
		my $gId=$row[0];
		my $gName=$row[1];
		my $nOfT=$row[2];
		my $tId=$row[3];
		my $tBiotype=$row[4];

		$ensembl{$gId}{'gName'}=$gName;
		$ensembl{$gId}{'nOfT'}=$nOfT;
		$ensembl{$gId}{'transcripts'}{$tId}=$tBiotype;
	}
}
close (INPUT);

my (%appris, %pdb);
if ($species eq 'hsa') {
        my $appris_input="$data_dir/$species/_appris.results.rel15.2May2013.v1.main.tsv";
        open (INPUT, "< $appris_input") or die "Could not open $appris_input: $!";
        while( my $row = <INPUT>)  {
        	chomp ($row);
        	my @row=split(/\t/, $row);

        	if ($row[5] eq 'PRINCIPAL') {
        		$appris{$row[2]}++;
        	}
        }
        close (INPUT);

        my $pdb_input="$data_dir/$species/_UniPdbCov.txt";
        open (INPUT, "<$pdb_input") or die "Could not open $pdb_input: $!";
        while( my $row = <INPUT>)  {
        	chomp ($row);
            $row =~ s/\r//g;
        	my @row=split(/\t/, $row);

        	my $gId=$row[1];
        	my $pdbId=$row[2];
        	my $coverage=sprintf("%.2f", $row[3]);
        	$pdb{$gId}{'pdbId'}=$pdbId;
        	$pdb{$gId}{'coverage'}=$coverage;
        }
        close (INPUT);
}

# Prepare output
open (OUT, ">$out_file") or die "Could not open $out_file: $!";

print OUT "gId:ensembl_gene_id ";
print OUT "gName:gene_name ";
print OUT "nOfT:number_of_annotated_transcripts ";
print OUT "C1.tId:major_transcript_-_condition_1 ";
print OUT "C1.principal:is_the_transcript_classified_as_principal_in_APPRIS?_-_condition_1 "
        if ($species eq 'hsa');
print OUT "C1.biotype:major_transcript_biotype_-_condition_1 ";
print OUT "C1.tExp:in_how_many_samples_is_the_transcript_detected_as_major?_-_condition_1 ";
print OUT "C1.gExp:in_how_many_samples_is_the_gene_expressed?_-_condition_1 ";
print OUT "C1.breadth:major_transcript_expression_breadth_-_condition_1 ";
print OUT "C2.tId:major_transcript_-_condition_2 ";
print OUT "C2.principal:is_the_transcript_classified_as_principal_in_APPRIS?_-_condition_2 "
        if ($species eq 'hsa');
print OUT "C2.biotype:major_transcript_biotype_-_condition_2 ";
print OUT "C2.tExp:in_how_many_samples_is_the_transcript_detected_as_major?_-_condition_2 ";
print OUT "C2.gExp:in_how_many_samples_is_the_gene_expressed?_-_condition_2 ";
print OUT "C2.breadth:major_transcript_expression_breadth_-_condition_2 ";
#"DE:is_the_gene_differentially_expressed?",
print OUT "pSimilarity:similarity_between_the_two_coding_sequences ";
print OUT "pdbEntry:is_there_a_PDB_entry_available? "
        if ($species eq 'hsa');
print OUT "pdbId:PDB_id "
        if ($species eq 'hsa');
print OUT "pdbCoverage:%_of_the_gene_represented_in_the_PDB_structure "
        if ($species eq 'hsa');
print OUT "rank:ranking_to_maximise_expression_breadth\n";

# Check for switch events
foreach my $gId (keys %recurrent_tx) {
	my $tId_cond1=$recurrent_tx{$gId}{'cond1'}{'recurrent_tx_id'};
	my $tId_cond2=$recurrent_tx{$gId}{'cond2'}{'recurrent_tx_id'};

	if ( $tId_cond1 ne $tId_cond2) {
		# Calculate major transcript expression breadth
		my $tExp_count_cond1=$recurrent_tx{$gId}{'cond1'}{'recurrent_tx_count'};
		my $gExp_count_cond1=&get_gExp_count([ @{ $major_tx{$gId}{'gExp'} }[@cond1] ]);
		my $breadth_cond1=&get_exp_breadth($tExp_count_cond1, $gExp_count_cond1);

		my $tExp_count_cond2=$recurrent_tx{$gId}{'cond2'}{'recurrent_tx_count'};
		my $gExp_count_cond2=&get_gExp_count([ @{ $major_tx{$gId}{'gExp'} }[@cond2] ]);
		my $breadth_cond2=&get_exp_breadth($tExp_count_cond2, $gExp_count_cond2);

		# Identify switch events
		if ($breadth_cond1 > 50 and $breadth_cond2 > 50) {

	        # Annotate switch events - gene info: gene name + number of transcripts
	        my $gName=$ensembl{$gId}{'gName'};
	        my $nOfT=$ensembl{$gId}{'nOfT'};

	        # Annotate switch events - transcript info: is principal? + transcript biotype
	        my ($isPrincipal_cond1, $tBiotype_cond1)=&get_transcript_info($gId, $tId_cond1);
	        my ($isPrincipal_cond2, $tBiotype_cond2)=&get_transcript_info($gId, $tId_cond2);

	 	# Annotate switch events - general info: pIdentity + pdb info
                my $pIdentity="NA";
        	my $pdbEntry="NA";
                my $pdbId="NA";
                my $pdbCoverage="NA";	 
	        if ($tBiotype_cond1 eq "protein_coding" and $tBiotype_cond2 eq "protein_coding") {

	            # pIdentity
	            my $fa_cond1="$data_dir/$species/_ensembl$ensembl_v.prot_seq/".substr($tId_cond1, 0, 12)."/$tId_cond1.fa";
	            my $fa_cond2="$data_dir/$species/_ensembl$ensembl_v.prot_seq/".substr($tId_cond2, 0, 12)."/$tId_cond2.fa";
	           	my $outdir_aln="$out_dir/prot_aln/".substr($gId, 0, 12);
	           	unless ( -e  "$out_dir/prot_aln/" ) { system("mkdir $out_dir/prot_aln/") };
	           	unless ( -e  $outdir_aln ) { system("mkdir $outdir_aln") };
	           	my $out_aln="$outdir_aln/$gId.needle.out";
				system("needle $fa_cond1 $fa_cond2 -auto stdout > $out_aln");
	            my $pIdentity=`cat $out_aln | grep Identity | awk -F '(' '{print \$2}' | awk -F '%' '{print \$1}' | sed 's/ //g'`;

	            # pdb structure
	            if ($species eq 'hsa' and defined $pdb{$gId}) {
		            $pdbEntry="YES";
		            $pdbId=$pdb{$gId}{'pdbId'};
		            $pdbCoverage=$pdb{$gId}{'coverage'};
	            } 
	        }

	        # Annotate switch events - general info: rank
	        my $exp=[ $tExp_count_cond1, $gExp_count_cond1, $tExp_count_cond2, $gExp_count_cond2 ];
			my $rank=&calculate_rank($exp);

			# Print output
			print OUT "$gId $gName $nOfT ";
			print OUT "$tId_cond1 ";
                        print OUT "$isPrincipal_cond1 " if ($species eq 'hsa');
                        print OUT "$tBiotype_cond1 $tExp_count_cond1 $gExp_count_cond1 $breadth_cond1 ";
			print OUT "$tId_cond2 ";
                        print OUT "$isPrincipal_cond2 " if ($species eq 'hsa');
                        print OUT "$tBiotype_cond2 $tExp_count_cond2 $gExp_count_cond2 $breadth_cond2 ";
			print OUT "$pIdentity ";
                        print OUT "$pdbEntry $pdbId $pdbCoverage " if ($species eq 'hsa');
                        print OUT "$rank\n";
		}		
	}
}
close (OUT);

# generate html
#require "./src/generate_html.pl";

my %arguments=(
	"input" 	=> 	$out_file,
	"out_dir"	=> $out_dir,
	"data_dir"	=> $data_dir,
	"species"		=> $species,
	"ensembl_v"	=> $ensembl_v,
	"cond1"  		=> $cond1,
	"cond2"  		=> $cond2
);

&generate_html(\%arguments);

# ##############################################################

sub adjust_columns {
	my $columns=$_[0];

	my @tmp=split("-", $columns);
	my @columns=($tmp[0]..$tmp[1]);
	foreach (@columns) { $_ += -3 };

 	return(@columns);
}

sub get_most_recurrent_tx {
	my $ref_subset_major_tx=$_[0];
	my $ref_columns=$_[1];

	my %subset_major_tx=%$ref_subset_major_tx;
	my @columns=@$ref_columns;
	my %count;	# how many times each transcript is detected as major?

	foreach my $i (@columns) {
		if (@{ $subset_major_tx{'gExp'} }[$i] >= $threshold_gexp) {
			my $tId=@{ $subset_major_tx{'major_tx_id'} }[$i];
			$count{$tId}++;
		}
	}

	my %result;
	my $recurrent_tx_id=(sort {$count{$b} <=> $count{$a}} keys %count)[0];
	if (defined $recurrent_tx_id) {
		$result{'recurrent_tx_id'}=$recurrent_tx_id;
		$result{'recurrent_tx_count'}=$count{$recurrent_tx_id};
	}
	return \%result;
}

sub get_gExp_count {
	my $ref_subset_gExp=$_[0];
	my @subset_gExp=@$ref_subset_gExp;

	my $gExp_count=0;
	foreach my $g (@subset_gExp) {
		if ($g >= $threshold_gexp) {
			$gExp_count++;
		}
	}
	return($gExp_count);
}

sub get_exp_breadth {
	my $tExp_count=$_[0];
	my $gExp_count=$_[1];

	my $breadth_tExp=sprintf("%.2f", $tExp_count/$gExp_count*100);
	return($breadth_tExp);
}

sub get_transcript_info {
    my $gId=$_[0];
    my $tId=$_[1];
  
  	# isPrincipal
    my $isPrincipal="NO";
    if (defined $appris{$tId}) { $isPrincipal="YES" };

  	# tBiotype
    my $tBiotype=$ensembl{$gId}{'transcripts'}{$tId};

    my @result=($isPrincipal, $tBiotype);
    return(@result);
}

sub calculate_rank {
	my $ref_exp=$_[0];
	my @exp=@$ref_exp;

	my $a=($exp[0]+$exp[1])*(1-abs($exp[0]-$exp[1])/max($exp[0],$exp[1]));
	my $b=($exp[2]+$exp[3])*(1-abs($exp[2]-$exp[3])/max($exp[2],$exp[3]));
	#my $b=($tExp_count_cond2+$gExp_count_cond2)*(1-abs($tExp_count_cond2-$gExp_count_cond2)/max($tExp_count_cond2,$gExp_count_cond2));
	my $rank=sprintf("%.2f", $a+$b);
	
	return($rank);
}