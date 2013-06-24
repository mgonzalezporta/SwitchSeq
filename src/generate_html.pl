#!/usr/bin/env perl

use strict;
use warnings;
use DBI;
use Text::Template;
use Getopt::Long;
use Data::Dumper;

# TO DO
# add info on the execution (date, options...) at the top of the html
# open prot_aln file
# R args
##############################################################

my ($input, $data_dir, $out_dir, $cond1, $cond2, $species, $ensembl_v, $plot);

GetOptions(
		'input|i=s'				=> \$input,
		'data_dir|d=s'			=> \$data_dir,
        'out_dir|o=s'			=> \$out_dir,
        'cond1|c1=s'			=> \$cond1,
        'cond2|c2=s'			=> \$cond2,
        'species|s=s'			=> \$species,
        'ensembl_v|e:i'			=> \$ensembl_v,
        'plot|p=s'				=> \$plot
);

# Arguments
my $input="$out_dir/switch.txt";

##############################################################

# Input
my %data;
my $total=0;
open INPUT, $input or die "Could not open $input: $!";
while( my $row = <INPUT>)  {
    chomp($row);
    my @row=split(/ /, $row);

    if ($.==1) {
    	# Header
    	$data{'header'} = \@row;
    	if ($plot eq 'starplots') {
			push @{ $data{"header"} }, "starplot:link_to_starplot";
		} elsif ($plot eq 'boxplots') {
			push @{ $data{"header"} }, "boxplot:link_to_boxplot";
		} elsif ($plot eq 'both') {
			push @{ $data{"header"} }, "boxplot:link_to_boxplot", "starplot:link_to_starplot";
		}
    } else {
    	$total++;

    	# Generate plots
    	my $gId=$row[0];
    	my $command="Rscript ./src/plots.R $plot $input $out_dir $data_dir $species $ensembl_v $cond1 $cond2 $gId";
    	#system($command);

    	# Classify switch events based on transcript biotype info
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

##############################################################

# Output
my $template = Text::Template->new(SOURCE => './src/index.tmpl')
 	or die "Couldn't construct template: $Text::Template::ERROR";

# print subpages first
my @all;
for my $comparison (keys %data) {
	if ($comparison ne 'header') {
		push @all, @{ $data{$comparison} };

		my %to_template = (
			colnames => $data{"header"},
			count => \%count,
			query => $data{$comparison},
		);

		my $result = $template->fill_in(HASH => \%to_template);
		my $outfile="$out_dir/$comparison.html";
		open (OUT, ">$outfile");
		if (defined $result) { print OUT $result }
			else { die "Couldn't fill in template: $Text::Template::ERROR" };
		close (OUT);
	}
}

# finally print index.html
my %to_template = (
	colnames => $data{"header"},
	count => \%count,
	query => \@all,
	);

print Dumper $data{"header"};
#print Dumper @all;

my $result = $template->fill_in(HASH => \%to_template);
my $outfile="$out_dir/index.html";
open (OUT, ">$outfile");
if (defined $result) { print OUT $result }
	else { die "Couldn't fill in template: $Text::Template::ERROR" };
close (OUT);

# copy css + js
system("cp -R $data_dir/css $out_dir");
system("cp -R $data_dir/js $out_dir");