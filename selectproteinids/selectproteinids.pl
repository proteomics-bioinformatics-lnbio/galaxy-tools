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

my $out = $ARGV[1];
open OUTFILE, ">", $out or die $!;

foreach(@lines){
    my @vec = split ' ', $_;
    print OUTFILE $vec[0], "\n";
}

close INFILE;
close OUTFILE;

