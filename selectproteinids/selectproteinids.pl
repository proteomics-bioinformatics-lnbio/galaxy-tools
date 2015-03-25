#!/usr/bin/perl -w

use strict;
use warnings;
use FileHandle;


# open file with the ids and get the lines #
my $in = $ARGV[0];
open INFILE, "<", $in or die $!;
seek(INFILE, 0, 0);

my @lines = <INFILE>;
shift @lines;

# flag to see wheter maintain contaminants or not
my $maintain_contaminants = $ARGV[1];

# output file with table filtered
my $out = $ARGV[2];
open OUTFILE, ">", $out or die $!;

# select first id of first column
foreach(@lines){
    my @vec = split ' ', $_;
    my @id = split ';', $vec[0];
    if($id[0] =~ m/^CON__/){
	if($maintain_contaminants eq "yes"){
	    print OUTFILE $id[0] =~ s/^CON__//r, "\n";
	}
    } else{
	print OUTFILE $id[0], "\n";
    }
}

close INFILE;
close OUTFILE;

