#!/usr/bin/perl -w

open FILE, "<Matrix_dataPSEA_test1_03072015.csv";
open OUT, ">input_fixed.csv";
my @line = <FILE>;

foreach(@line){
    $_ =~ s/,/\t/g;
    print OUT $_;
}
close(FILE);
close(OUT);
