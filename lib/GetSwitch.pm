#!/usr/bin/perl

package GetSwitch;
use strict;
use warnings;
use Exporter;
use List::Util qw[ max ];
use Text::Template;
use FindBin qw($Bin);
use Data::Dumper;
use JSON;

our @ISA= qw( Exporter );
our @EXPORT = qw( get_switch );

sub get_switch {
	## collect arguments
	my $ref_arguments=$_[0];
	my %arguments=%{$ref_arguments};
	my $input=$arguments{'input'};
	my $out_dir=$arguments{'out_dir'};
	my $data_dir=$arguments{'data_dir'};
	my $species=$arguments{'species'};
	my $ensembl_v=$arguments{'ensembl_v'};
	my $cond1=$arguments{'cond1'};
	my $cond2=$arguments{'cond2'};
	my $threshold_gexp=$arguments{'threshold_gexp'};

	## prepare dir structure
	_prepare_dir_structure($ref_arguments);

	## progress
	print "# Obtaining and annotating alternative splicing switch events...\n";

	## calculations
	my $ref_major_tx=_obtain_major_tx($arguments{'input'});
	#print Dumper %$ref_major_tx;
	my $ref_recurrent_major_tx=_obtain_recurrent_major_tx($ref_major_tx, $ref_arguments);
	#print Dumper %$ref_recurrent_major_tx;
	my $ref_switch=_obtain_switch_events($ref_major_tx, $ref_recurrent_major_tx, $ref_arguments);
	#print Dumper %$ref_switch;

	## print output
	_print_txt($ref_switch, $ref_arguments);
	_print_json($ref_switch, $ref_arguments);
	_print_html($ref_switch, $ref_arguments);

	## progress
	my $count=keys %$ref_switch;
	print "# Switch events obtained for $count protein coding genes.\n";

}

sub _prepare_dir_structure {
	my $ref_arguments=$_[0];
	my %arguments=%$ref_arguments;
	my $out_dir=$arguments{'out_dir'};

	unless ( -e $out_dir ) { system("mkdir $out_dir") };
	unless ( -e "$out_dir/data" ) { system("mkdir $out_dir/data") };
	unless ( -e  "$out_dir/data/prot_aln/" ) { system("mkdir $out_dir/data/prot_aln/") };
	unless ( -e  "$out_dir/data/plots/" ) { system("mkdir $out_dir/data/plots/") };

	system("cp -R $Bin/resources/css $out_dir");
	system("cp -R $Bin/resources/js $out_dir");
}

sub _adjust_columns {
	my $columns=$_[0];

	my @tmp=split("-", $columns);
	my @columns=($tmp[0]..$tmp[1]);
	foreach (@columns) { $_ += -3 };

 	return(\@columns);
}

sub _obtain_major_tx {
	## obtain the major transcript in each sample
	my $input=$_[0];
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

	return \%major_tx;
}

sub _obtain_recurrent_major_tx {
	## obtain the most recurrent major transcripts in each condition
	my $ref_major_tx=$_[0];
	my $ref_arguments=$_[1];

	my %major_tx=%$ref_major_tx;
	my %arguments=%{$ref_arguments};
	my $cond1=$arguments{'cond1'};
	my $ref_cond1=_adjust_columns($cond1);
	my $cond2=$arguments{'cond2'};
	my $ref_cond2=_adjust_columns($cond2);
	my $threshold_gexp=$arguments{'threshold_gexp'};
	my %recurrent_tx;
		# gId => | cond1 => | recurrent_tx_id
		#                   | recurrent_tx_count
		# 		 | cond2 => | recurrent_tx_id
		# 		            | recurrent_tx_count

	foreach my $gId (keys %major_tx) {
		my %subset_major_tx=%{ $major_tx{$gId} };

		$recurrent_tx{$gId}{'cond1'}=_get_most_recurrent_tx(\%subset_major_tx, $ref_cond1, $threshold_gexp);
		$recurrent_tx{$gId}{'cond2'}=_get_most_recurrent_tx(\%subset_major_tx, $ref_cond2, $threshold_gexp);
	}

	# if any of the conditions doesn't have any transcript expressed, discard the whole gene
	foreach my $gId (keys %recurrent_tx) {
		if (! defined($recurrent_tx{$gId}{'cond1'}{'recurrent_tx_id'}) or
			! defined($recurrent_tx{$gId}{'cond2'}{'recurrent_tx_id'})) {
			delete $recurrent_tx{$gId};
		}
	}

	return \%recurrent_tx;
}

sub _get_most_recurrent_tx {
	my $ref_subset_major_tx=$_[0];
	my %subset_major_tx=%$ref_subset_major_tx;
	my $ref_columns=$_[1];	
	my @columns=@$ref_columns;
	my $threshold_gexp=$_[2];
	my %count;
	my %result;

	## count how many times each transcript is detected as major
	foreach my $i (@columns) {
		if (@{ $subset_major_tx{'gExp'} }[$i] >= $threshold_gexp) {
			my $tId=@{ $subset_major_tx{'major_tx_id'} }[$i];
			$count{$tId}++;
		}
	}

	## get the transcript with the highest count
	my $recurrent_tx_id=(sort {$count{$b} <=> $count{$a}} keys %count)[0];
	if (defined $recurrent_tx_id) {
		$result{'recurrent_tx_id'}=$recurrent_tx_id;
		$result{'recurrent_tx_count'}=$count{$recurrent_tx_id};
	}
	return \%result;
}

sub _obtain_switch_events {
	my $ref_major_tx=$_[0];
	my $ref_recurrent_major_tx=$_[1];
	my $ref_arguments=$_[2];

	my %arguments=%$ref_arguments;
	my $data_dir=$arguments{'data_dir'};
	my $species=$arguments{'species'};
	my $ensembl_v=$arguments{'ensembl_v'};
	my $out_dir=$arguments{'out_dir'};
	my $output="$out_dir/switch.txt";
	
	## load data
	my $ensembl_input1="$data_dir/$species.$ensembl_v/ensembl1.txt";
	my $ensembl_input2="$data_dir/$species.$ensembl_v/ensembl2.txt";
	my $ref_ensembl=_load_ensembl($ensembl_input1, $ensembl_input2);
	#print Dumper %$ref_ensembl;

	my $appris_input="$data_dir/$species.$ensembl_v/appris_data.principal.txt";
	my $ref_appris=_load_appris($appris_input);

	## find and annotate
	my $ref_switch=_find_switch($ref_major_tx, $ref_recurrent_major_tx, $ref_arguments);
	
	$ref_switch=_annotate_switch($ref_switch, $ref_ensembl, $ref_appris, $ref_arguments);
	return $ref_switch;
}

sub _load_ensembl {
	my $ensembl_input1=$_[0];
	my $ensembl_input2=$_[1];
	
	my %ensembl;

	open (INPUT, "< $ensembl_input1") or die "Could not open $ensembl_input1: $!";
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

	open (INPUT, "<$ensembl_input2") or die "Could not open $ensembl_input2: $!";
	while( my $row = <INPUT>)  {
		chomp ($row);

		if ($row =~ /^ENSG/) {
			my @row=split(/\s+/, $row);
			my $gId=$row[0];
			my $uniprotId=$row[1];

			$ensembl{$gId}{'uniprotId'}=$uniprotId;
		}
	}
	close (INPUT);

	return \%ensembl;
}

sub _load_appris {
	my $appris_input=$_[0];

	my %appris;

        open (INPUT, "<$appris_input") or die "Could not open $appris_input: $!";
        while( my $row = <INPUT>)  {
        	chomp ($row);
        	my @row=split(/\s+/, $row);
    
		$appris{$row[1]}++;
        }
        close (INPUT);
   
        return \%appris;
}

sub _get_header {
	my $txt_output=$_[0];
	my @header;

	## general gene info
	push(@header, "gId:ensembl_gene_id ");
	push(@header, "gName:gene_name ");
	push(@header, "nOfT:number_of_annotated_transcripts ");

	## condition 1
	push(@header, "C1.tId:major_transcript_-_condition_1 ");
	push(@header, "C1.principal:is_the_transcript_classified_as_principal_in_APPRIS?_-_condition_1 ");
	push(@header, "C1.biotype:major_transcript_biotype_-_condition_1 ");
	push(@header, "C1.tExp:in_how_many_samples_is_the_transcript_detected_as_major?_-_condition_1 ");
	push(@header, "C1.gExp:in_how_many_samples_is_the_gene_expressed?_-_condition_1 ");
	push(@header, "C1.breadth:major_transcript_expression_breadth_-_condition_1 ");

	## condition 2
	push(@header, "C2.tId:major_transcript_-_condition_2 ");
	push(@header, "C2.principal:is_the_transcript_classified_as_principal_in_APPRIS?_-_condition_2 ");
	push(@header, "C2.biotype:major_transcript_biotype_-_condition_2 ");
	push(@header, "C2.tExp:in_how_many_samples_is_the_transcript_detected_as_major?_-_condition_2 ");
	push(@header, "C2.gExp:in_how_many_samples_is_the_gene_expressed?_-_condition_2 ");
	push(@header, "C2.breadth:major_transcript_expression_breadth_-_condition_2 ");

	## general info again
	push(@header, "pIdentity:percentage_identity_between_the_two_coding_sequences ");
	push(@header, "pdbEntry:is_there_any_PDB_entry_available? ");
	if (!$txt_output) {
		push(@header, "distrplot:link_to_distrplot");
		push(@header, "starplot:link_to_starplot");	
	}
	push(@header, "rank:ranking_to_maximise_expression_breadth");

	return(\@header);
}

sub _find_switch {
	my $ref_major_tx=$_[0];
	my $ref_recurrent_major_tx=$_[1];
	my $ref_arguments=$_[2];
	
	my %major_tx=%$ref_major_tx;
	my %recurrent_tx=%$ref_recurrent_major_tx;
	my %switch;
	my %arguments=%$ref_arguments;
	my $cond1=$arguments{'cond1'};
	my $ref_cond1=_adjust_columns($cond1);
	my $cond2=$arguments{'cond2'};
	my $ref_cond2=_adjust_columns($cond2);
	my $threshold_gexp=$arguments{'threshold_gexp'};

	foreach my $gId (keys %recurrent_tx) {
		my $tId_cond1=$recurrent_tx{$gId}{'cond1'}{'recurrent_tx_id'};
		my $tId_cond2=$recurrent_tx{$gId}{'cond2'}{'recurrent_tx_id'};

		if ( $tId_cond1 ne $tId_cond2) {
			## calculate major transcript expression breadth
			my $tExp_count_cond1=$recurrent_tx{$gId}{'cond1'}{'recurrent_tx_count'};
			my $gExp_count_cond1=_get_gexp_count([ @{ $major_tx{$gId}{'gExp'} }[@$ref_cond1] ], $threshold_gexp);
			my $breadth_cond1=_get_exp_breadth($tExp_count_cond1, $gExp_count_cond1);

			my $tExp_count_cond2=$recurrent_tx{$gId}{'cond2'}{'recurrent_tx_count'};
			my $gExp_count_cond2=_get_gexp_count([ @{ $major_tx{$gId}{'gExp'} }[@$ref_cond2] ], $threshold_gexp);
			my $breadth_cond2=_get_exp_breadth($tExp_count_cond2, $gExp_count_cond2);

			## report switch events
			if ($breadth_cond1 > 50 and $breadth_cond2 > 50) {
				$switch{$gId}{'C1.tId'}=$recurrent_tx{$gId}{'cond1'}{'recurrent_tx_id'};
				$switch{$gId}{'C1.tExp'}=$tExp_count_cond1;
				$switch{$gId}{'C1.gExp'}=$gExp_count_cond1;
				$switch{$gId}{'C1.breadth'}=$breadth_cond1;

				$switch{$gId}{'C2.tId'}=$recurrent_tx{$gId}{'cond2'}{'recurrent_tx_id'};
				$switch{$gId}{'C2.tExp'}=$tExp_count_cond2;
				$switch{$gId}{'C2.gExp'}=$gExp_count_cond2;
				$switch{$gId}{'C2.breadth'}=$breadth_cond2;
			}
		}
	}

	return \%switch;
}

sub _get_gexp_count {
	my $ref_subset_gExp=$_[0];
	my $threshold_gexp=$_[1];

	my @subset_gExp=@$ref_subset_gExp;

	my $gExp_count=0;
	foreach my $g (@subset_gExp) {
		if ($g >= $threshold_gexp) {
			$gExp_count++;
		}
	}
	return($gExp_count);
}

sub _get_exp_breadth {
	my $tExp_count=$_[0];
	my $gExp_count=$_[1];

	my $breadth_tExp=sprintf("%.2f", $tExp_count/$gExp_count*100);
	return($breadth_tExp);
}

sub _annotate_switch {
	my $ref_switch=$_[0];
	my $ref_ensembl=$_[1];
	my $ref_appris=$_[2];
	my $ref_arguments=$_[3];

	my %switch=%$ref_switch;
	my %ensembl=%$ref_ensembl;
	# args: data_dir species, ensembl_v, outdir

	foreach my $gId (keys %switch) {
		$switch{$gId}{'gName'}=$ensembl{$gId}{'gName'};
		$switch{$gId}{'nOfT'}=$ensembl{$gId}{'nOfT'};
        $switch{$gId}{'C1.principal'}=_is_principal( $switch{$gId}{'C1.tId'}, $ref_appris );
        $switch{$gId}{'C2.principal'}=_is_principal( $switch{$gId}{'C2.tId'}, $ref_appris );
        $switch{$gId}{'C1.biotype'}=_get_tx_biotype( $gId, $switch{$gId}{'C1.tId'}, $ref_ensembl );
        $switch{$gId}{'C2.biotype'}=_get_tx_biotype( $gId, $switch{$gId}{'C2.tId'}, $ref_ensembl );
 		$switch{$gId}{'rank'}=_calculate_rank( $switch{$gId} );	 		
        $switch{$gId}{'pIdentity'}="NA";
    	$switch{$gId}{'pdbEntry'}="NO";

        if ($switch{$gId}{'C1.biotype'} eq "protein_coding" and 
        	$switch{$gId}{'C2.biotype'} eq "protein_coding") {
	            ## pIdentity
	            $switch{$gId}{'pIdentity'}=_get_prot_identity($ref_arguments, $gId, $switch{$gId});

	            ## pdbEntry
	            if (defined $ensembl{$gId}{'uniprotId'}) {
		            $switch{$gId}{'pdbEntry'}=$ensembl{$gId}{'uniprotId'};
	            } 
        }
	       
	}	
	return \%switch;
}

sub _get_prot_identity {
	my $ref_arguments=$_[0];
	my $gId=$_[1];
	my $ref_subset_switch=$_[2];

	my %arguments=%{$ref_arguments};
	my $out_dir=$arguments{'out_dir'};
	my $data_dir=$arguments{'data_dir'};
	my $species=$arguments{'species'};
	my $ensembl_v=$arguments{'ensembl_v'};

	my %subset_switch=%$ref_subset_switch;
	my $tId_cond1=$subset_switch{'C1.tId'};
	my $tId_cond2=$subset_switch{'C2.tId'};
        my ($species_id_cond1, $numeric_id_cond1) = $tId_cond1 =~ /([a-zA-Z]+)(\d+)/;
	my $fa_cond1="$data_dir/$species.$ensembl_v/prot_seq/".$species_id_cond1.substr($numeric_id_cond1, 0, 8)."/$tId_cond1.fa";
	my ($species_id_cond2, $numeric_id_cond2) = $tId_cond2 =~ /([a-zA-Z]+)(\d+)/;
        my $fa_cond2="$data_dir/$species.$ensembl_v/prot_seq/".$species_id_cond2.substr($numeric_id_cond2, 0, 8)."/$tId_cond2.fa";

	my ($species_id, $numeric_id) = $gId =~ /([a-zA-Z]+)(\d+)/;
   	my $outdir_aln="$out_dir/data/prot_aln/".$species_id.substr($numeric_id, 0, 8);
   	unless ( -e  $outdir_aln ) { system("mkdir $outdir_aln") };
   	my $out_aln="$outdir_aln/$gId.needle_mod.out";
   	my $pIdentity="NA";

   	## run needle to get protein identity
	open(PIPE, "needle $fa_cond1 $fa_cond2 -auto stdout |") or die "Cannot open needle output: $!";
	my @needle_output=<PIPE>;
	close(PIPE);

	my $ref_fa_cond1=_read_fa($fa_cond1);
	my $ref_fa_cond2=_read_fa($fa_cond2);

	## print output
	open(OUT, ">$out_aln") or die "Cannot open $out_aln: $!";
	foreach my $row (@needle_output) {
		print OUT $row;

		if ($row =~ /# Identity.+\((\d+\.\d+)/) {
			$pIdentity=$1;
		}
	}
	
	print OUT "\n# Protein sequences in FASTA format\n";

	foreach my $row (@$ref_fa_cond1) {
		print OUT $row;
	}

	foreach my $row (@$ref_fa_cond2) {
		print OUT $row;
	}

    close(OUT);

    ##
    return($pIdentity);
}

sub _read_fa {
	my $fa=$_[0];

	open(IN_FA, "<$fa") or die "Cannot open $fa: $!";
	my @fa=<IN_FA>;
	close(IN_FA);

	return(\@fa);
}

sub _is_principal {
	my $tId=$_[0];
	my $ref_appris=$_[1];

	my %appris=%$ref_appris;
    my $isPrincipal="NO";

    if (defined $appris{$tId}) { $isPrincipal="YES" };

    return($isPrincipal);
}

sub _get_tx_biotype {
	my $gId=$_[0];
 	my $tId=$_[1];
 	my $ref_ensembl=$_[2];

 	my %ensembl=%$ref_ensembl;

    my $tBiotype=$ensembl{$gId}{'transcripts'}{$tId};
    
    return $tBiotype;
}

sub _calculate_rank {
	my $ref_subset_switch=$_[0];

	my %subset_switch=%$ref_subset_switch;
	my @exp=($subset_switch{'C1.tExp'}, $subset_switch{'C1.gExp'}, $subset_switch{'C2.tExp'}, $subset_switch{'C2.gExp'});
	#my $exp=[ $tExp_count_cond1, $gExp_count_cond1, $tExp_count_cond2, $gExp_count_cond2 ];

	my $a=($exp[0]+$exp[1])*(1-abs($exp[0]-$exp[1])/max($exp[0],$exp[1]));
	my $b=($exp[2]+$exp[3])*(1-abs($exp[2]-$exp[3])/max($exp[2],$exp[3]));

	my $rank=sprintf("%.2f", $a+$b);
	
	return($rank);
}

sub _print_txt {
	my $ref_switch=$_[0];
	my $ref_arguments=$_[1];

	my %switch=%$ref_switch;
	my %arguments=%$ref_arguments;
	my $out_dir=$arguments{'out_dir'};
	my $out_file="$out_dir/data/switch.txt";

	## prepare output
	open(my $fh, ">$out_file") or die "Could not open $out_file: $!";
	my $ref_header=_get_header(1);
	print $fh "@$ref_header \n";

	foreach my $gId (keys %switch) {
		print $fh "$gId $switch{$gId}{'gName'} $switch{$gId}{'nOfT'} ";
		print $fh "$switch{$gId}{'C1.tId'} $switch{$gId}{'C1.principal'} ";
		print $fh "$switch{$gId}{'C1.biotype'} $switch{$gId}{'C1.tExp'} ";
		print $fh "$switch{$gId}{'C1.gExp'} $switch{$gId}{'C1.breadth'} ";
		print $fh "$switch{$gId}{'C2.tId'} $switch{$gId}{'C2.principal'} ";
		print $fh "$switch{$gId}{'C2.biotype'} $switch{$gId}{'C2.tExp'} ";
		print $fh "$switch{$gId}{'C2.gExp'} $switch{$gId}{'C2.breadth'} ";
		print $fh "$switch{$gId}{'pIdentity'} $switch{$gId}{'pdbEntry'} ";
		print $fh "$switch{$gId}{'rank'}\n";
	}
    close($fh);
}

sub _print_json {
	my $ref_switch=$_[0];
	my $ref_arguments=$_[1];

	my %arguments=%$ref_arguments;
	my $out_dir=$arguments{'out_dir'};
	my $out_file="$out_dir/js/data.js";

	## convert initial hash to array of hashes
	my %switch=%$ref_switch;
	my @array_input;
	foreach my $gId (keys %switch) {
		my %tmp_hash;
		
		$tmp_hash{'gId'}=$gId;
		## js doesn't allow . in hash keys
		foreach my $k (keys %{$switch{$gId}}) {
			my $value=$switch{$gId}{$k};
			$k =~ s/\./_/;
			$tmp_hash{$k}=$value;
		}
		
		push(@array_input, \%tmp_hash)
	}

	## print js variable
	open(OUT, ">$out_file") or die "Could not open $out_file: $!";
	print OUT "var data = \n";
	print OUT to_json(\@array_input, { pretty => 1 });
	close(OUT);
}

sub _print_html {
	my $ref_switch=$_[0];	
	my $ref_arguments=$_[1];

	my %switch=%{$ref_switch};
	my %arguments=%{$ref_arguments};
	my $out_dir=$arguments{'out_dir'};

	my $colnames=_get_header(0);
	my %count=(
		pc_to_pc => 0,
		pc_to_nmd => 0,
		pc_to_ri => 0,
		pc_to_pt => 0,
		nmd_to_pc => 0,
		ri_to_pc => 0,
		pt_to_pc => 0,
		other => 0
		);

	foreach my $gId (keys %switch) {
		if ($switch{$gId}{'C1.biotype'} eq 'protein_coding' 
			and $switch{$gId}{'C2.biotype'} eq 'protein_coding') {
			$count{'pc_to_pc'}++;
		} elsif ($switch{$gId}{'C1.biotype'} eq 'protein_coding' 
			and $switch{$gId}{'C2.biotype'} eq 'nonsense_mediated_decay') {
			$count{'pc_to_nmd'}++;
		} elsif ($switch{$gId}{'C1.biotype'} eq 'protein_coding' 
			and $switch{$gId}{'C2.biotype'} eq 'retained_intron') {
			$count{'pc_to_ri'}++;
		} elsif ($switch{$gId}{'C1.biotype'} eq 'protein_coding' 
			and $switch{$gId}{'C2.biotype'} eq 'processed_transcript') {
			$count{'pc_to_pt'}++;
		} elsif ($switch{$gId}{'C1.biotype'} eq 'nonsense_mediated_decay' 
			and $switch{$gId}{'C2.biotype'} eq 'protein_coding') {
			$count{'nmd_to_pc'}++;
		} elsif ($switch{$gId}{'C1.biotype'} eq 'retained_intron' 
			and $switch{$gId}{'C2.biotype'} eq 'protein_coding') {
			$count{'ri_to_pc'}++;
		} elsif ($switch{$gId}{'C1.biotype'} eq 'processed_transcript' 
			and $switch{$gId}{'C2.biotype'} eq 'protein_coding') {
			$count{'pt_to_pc'}++;
		} else {
			$count{'other'}++;
		}

		## generate plots
		my $expdata=$arguments{'input'};
		my $annot=$arguments{'data_dir'}."/".
			"$arguments{'species'}.$arguments{'ensembl_v'}/".
			"ensembl1.txt";
		my $cond1=$arguments{'cond1'};
		my $cond2=$arguments{'cond2'};
		my $command="R CMD BATCH --no-save ".
		  "\"--args bin='$Bin' gId='$gId' expdata='$expdata' annot='$annot' cond1='$cond1' cond2='$cond2' outdir='$out_dir'\" ". 
		  "$Bin/scripts/generate_plots.R /dev/null";
		print $command."\n";
		system($command);
	}

	foreach my $category (keys %count) {
		$count{'total'} += $count{$category};
	}

	## print index.html
	my %to_template = (
		info 	   => \%arguments,
		colnames   => $colnames,
		count 	   => \%count,
	);
	my $outfile="$out_dir/index.html";
	_fill_template(\%to_template, $outfile);
}

sub _fill_template {
	my $ref_to_template=$_[0];
	my $outfile=$_[1];

	my $template = Text::Template->new(SOURCE => "$Bin/templates/index.tmpl")
	 	or die "Couldn't construct template: $Text::Template::ERROR";
	my $result = $template->fill_in(HASH => $ref_to_template);

	open (OUT, ">$outfile");
	if (defined $result) { print OUT $result }
		else { die "Couldn't fill in template: $Text::Template::ERROR" };
	close (OUT);
}

1;
