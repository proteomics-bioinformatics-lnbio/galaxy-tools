#!/usr/bin/perl -w
# Developed by Flavia Vischi Winck and Mateus Bellomo_2014 (flaviavw[at]gmail.com)
# usage : perl cat_2.pl <FASTA file1> <FASTA file2> <outputfile>

open FILE1, "<", $ARGV[0] or die $!;
my @file1 = <FILE1>;

open FILE2, "<", $ARGV[1] or die $!;
my @file2 = <FILE2>;


open SEL, '>', $ARGV[2] or die $!;
foreach (@file1) {
    print SEL $_;
}

foreach (@file2){
    print SEL $_;
}

close SEL;

