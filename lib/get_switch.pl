#!/usr/bin/perl



# Check for switch events
# foreach my $gId (keys %recurrent_tx) {
# 	my $tId_cond1=$recurrent_tx{$gId}{'cond1'}{'recurrent_tx_id'};
# 	my $tId_cond2=$recurrent_tx{$gId}{'cond2'}{'recurrent_tx_id'};

# 	if ( $tId_cond1 ne $tId_cond2) {
# 		# Calculate major transcript expression breadth
# 		my $tExp_count_cond1=$recurrent_tx{$gId}{'cond1'}{'recurrent_tx_count'};
# 		my $gExp_count_cond1=&get_gExp_count([ @{ $major_tx{$gId}{'gExp'} }[@cond1] ]);
# 		my $breadth_cond1=&get_exp_breadth($tExp_count_cond1, $gExp_count_cond1);

# 		my $tExp_count_cond2=$recurrent_tx{$gId}{'cond2'}{'recurrent_tx_count'};
# 		my $gExp_count_cond2=&get_gExp_count([ @{ $major_tx{$gId}{'gExp'} }[@cond2] ]);
# 		my $breadth_cond2=&get_exp_breadth($tExp_count_cond2, $gExp_count_cond2);

# 		# Identify switch events
# 		if ($breadth_cond1 > 50 and $breadth_cond2 > 50) {

# 	        # Annotate switch events - gene info: gene name + number of transcripts
# 	        my $gName=$ensembl{$gId}{'gName'};
# 	        my $nOfT=$ensembl{$gId}{'nOfT'};

# 	        # Annotate switch events - transcript info: is principal? + transcript biotype
# 	        my ($isPrincipal_cond1, $tBiotype_cond1)=&get_transcript_info($gId, $tId_cond1);
# 	        my ($isPrincipal_cond2, $tBiotype_cond2)=&get_transcript_info($gId, $tId_cond2);

# 	 		# Annotate switch events - general info: pIdentity + pdb info
#             my $pIdentity="NA";
#         	my $pdbEntry="NA";
#             my $pdbId="NA";
#             my $pdbCoverage="NA";	 
# 	        if ($tBiotype_cond1 eq "protein_coding" and $tBiotype_cond2 eq "protein_coding") {

# 	            # pIdentity
# 	            my $fa_cond1="$data_dir/$species/_ensembl$ensembl_v.prot_seq/".substr($tId_cond1, 0, 12)."/$tId_cond1.fa";
# 	            my $fa_cond2="$data_dir/$species/_ensembl$ensembl_v.prot_seq/".substr($tId_cond2, 0, 12)."/$tId_cond2.fa";
# 	           	my $outdir_aln="$out_dir/prot_aln/".substr($gId, 0, 12);
# 	           	unless ( -e  "$out_dir/prot_aln/" ) { system("mkdir $out_dir/prot_aln/") };
# 	           	unless ( -e  $outdir_aln ) { system("mkdir $outdir_aln") };
# 	           	my $out_aln="$outdir_aln/tmp.$gId.needle.out";
# 				system("needle $fa_cond1 $fa_cond2 -auto stdout > $out_aln");
# 	            $pIdentity=`cat $out_aln | grep Identity | awk -F '(' '{print \$2}' | awk -F '%' '{print \$1}' | sed 's/ //g'`;
# 	            chomp $pIdentity;

# 	            ## append maistas-formatted fa
# 	            system("echo '\n# Input for MAISTAS\n# (http://maistas.bioinformatica.crs4.it/)\n' >> $out_aln");
# 	            my $out="$outdir_aln/$gId.needle_mod.out";
# 	            system("cat $out_aln $fa_cond1 $fa_cond2 > $out");
# 	            system("rm $out_aln");

# 	            # pdb structure
# 	            if ($species eq 'hsa' and defined $pdb{$gId}) {
# 		            $pdbEntry="YES";
# 		            $pdbId=$pdb{$gId}{'pdbId'};
# 		            $pdbCoverage=$pdb{$gId}{'coverage'};
# 	            } 
# 	        }

# 	        # Annotate switch events - general info: rank
# 	        my $exp=[ $tExp_count_cond1, $gExp_count_cond1, $tExp_count_cond2, $gExp_count_cond2 ];
# 			my $rank=&calculate_rank($exp);

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



