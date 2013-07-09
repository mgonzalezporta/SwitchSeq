#!/usr/bin/perl

package GetSwitch;
use strict;
use warnings;
use Exporter;
#use Text::Template;
use Bio::EnsEMBL::Registry;
use Bio::SeqIO;
use List::Util qw[ max ];

our @ISA= qw( Exporter );
our @EXPORT = qw( get_switch );

sub get_switch {

}

sub generate_html {

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

	## define extra variables
	my %data;
	my %colnames_index;
	my $total=0;

	## process switch.txt file
	## create plots
	open INPUT, $input or die "Could not open $input: $!";
	while( my $row = <INPUT>)  {
	    chomp($row);
	    my @row=split(/\s+/, $row);

	    if ($.==1) {
	    	$data{'header'} = \@row;
			push @{ $data{"header"} }, "boxplot:link_to_boxplot", "starplot:link_to_starplot";


			my $n=@{ $data{'header'} };
			for (my $i=0; $i<$n; $i++) {
				my @id=split(/:/, ${ $data{'header'} }[$i]);
				$colnames_index{$id[0]}=$i;
			}

	    } else {
	    	$total++;

	    	## generate plots
	    	my $gId=$row[0];
	    	my $plot="both";	## change this
	    	my $command="Rscript ./src/plots.R $plot $input $out_dir $data_dir $species $ensembl_v $cond1 $cond2 $gId";
	    	#system($command);

	    	## classify switch events based on transcript biotype info
		    my $tBiotype_A=$row[5];
		    my $tBiotype_B=$row[11];
		    
	    	if ($tBiotype_A eq 'protein_coding' and $tBiotype_B eq 'protein_coding') {
	    		push @{$data{'pc_to_pc'}}, [ @row ] ;
	    	} elsif ($tBiotype_A eq 'protein_coding' and $tBiotype_B eq 'nonsense_mediated_decay') {
	    		push @{$data{'pc_to_nmd'}}, [ @row ];
	    	} elsif ($tBiotype_A eq 'protein_coding' and $tBiotype_B eq 'retained_intron') {
	    		push @{$data{'pc_to_ri'}}, [ @row ];
	    	} elsif ($tBiotype_A eq 'protein_coding' and $tBiotype_B eq 'processed_transcript') {
				push @{$data{'pc_to_pt'}}, [ @row ];
			} elsif ($tBiotype_A eq 'nonsense_mediated_decay' and $tBiotype_B eq 'protein_coding') {
				push @{$data{'nmd_to_pc'}}, [ @row ];
			} elsif ($tBiotype_A eq 'retained_intron' and $tBiotype_B eq 'protein_coding') {
				push @{$data{'ri_to_pc'}}, [ @row ];
			} elsif ($tBiotype_A eq 'processed_transcript' and $tBiotype_B eq 'protein_coding') {
				push @{$data{'pt_to_pc'}}, [ @row ];
			} else {
				push @{$data{'other'}}, [ @row ];
			} 
	    }
	}
	close INPUT;

	## save summary statistics
	my %count = (
		pc_to_pc => scalar( @{ $data{'pc_to_pc'} } ),
		pc_to_nmd => scalar( @{ $data{'pc_to_nmd'} } ),
		pc_to_ri => scalar( @{ $data{'pc_to_ri'} } ),
		pc_to_pt => scalar( @{ $data{'pc_to_pt'} } ),
		nmd_to_pc => scalar( @{ $data{'nmd_to_pc'} } ),
		ri_to_pc => scalar( @{ $data{'ri_to_pc'} } ),
		pt_to_pc => scalar( @{ $data{'pt_to_pc'} } ),
		other => scalar( @{ $data{'other'} } ),
		total => $total
	);

	## output

	## print subpages
	## collect all data
	my @all;
	for my $comparison (keys %data) {
		if ($comparison ne 'header') {
			push @all, @{ $data{$comparison} };

			my %to_template = (
				info 			=> \%arguments,		# hash	
				colnames 		=> $data{"header"},		# array
				colnames_index 	=> \%colnames_index,	# hash
				count 			=> \%count,				# hash
				query 			=> $data{$comparison}	# array of arrays
			);
			my $outfile="$out_dir/$comparison.html";
			&fill_template(\%to_template, $outfile);
		}
	}

	## print index.html
	my %to_template = (
		info 			=> \%arguments,
		colnames 		=> $data{"header"},
		colnames_index 	=> \%colnames_index,
		count 			=> \%count,
		query 			=> \@all
	);
	my $outfile="$out_dir/index.html";
	&fill_template(\%to_template, $outfile);

	## copy css + js
	system("cp -R $data_dir/css $out_dir");
	system("cp -R $data_dir/js $out_dir");

}

sub fill_template {
	my $ref_to_template=$_[0];
	my $outfile=$_[1];

	my $template = Text::Template->new(SOURCE => './src/index.tmpl')
	 	or die "Couldn't construct template: $Text::Template::ERROR";
	my $result = $template->fill_in(HASH => $ref_to_template);

	open (OUT, ">$outfile");
	if (defined $result) { print OUT $result }
		else { die "Couldn't fill in template: $Text::Template::ERROR" };
	close (OUT);
}

1;