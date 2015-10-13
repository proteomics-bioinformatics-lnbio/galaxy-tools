#!/usr/bin/perl

#teste terminal: perl ontologymap.pl primary_213_intensity_visible_tabular.csv GO out.csv

use DBI;
use strict;
use Data::Dumper;
use Data::Dumper::Simple;
use GO::Parser;

# regex for uniprot ID
my $regex = '(\w*)(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-?[0-9]*)$';

my $db_to_use = $ARGV[1];


my $driver   = "SQLite"; 
my $database = "go_annotation.db";
my $dsn = "DBI:$driver:dbname=$database";
my $userid = "";
my $password = "";
my $dbh = DBI->connect($dsn, $userid, $password, { RaiseError => 1 }) 
                      or die $DBI::errstr;
print $dbh->trace(2);

# $hash{GO_ID}{GO_TERM} = [uniprot ids list]
my %hash = ();

my $table = "UNIPROT_ANNOTATION";

my $sth = $dbh->prepare("select GO_ID from $table where DB_OBJECT_ID = ?");


my $parser = new GO::Parser({handler=>'obj'}); # create parser object
$parser->parse("go.obo"); # parse file -> objects
my $graph = $parser->handler->graph;  # get L<GO::Model::Graph> object

open DATA, $ARGV[0];
my $iline = 0;
foreach(<DATA>){
    if($iline < 2){
	$iline++;
	next;
    }
    my @data = split '\t', $_;
    my @aux = split($regex, $data[1]);
    my $uniprot = $aux[2];
    if($uniprot eq ""){ next; }
 
    my @go_ids = select_db($uniprot);

    my %auxhash = map{$_ => 1} @go_ids;
    my @unique = keys %auxhash;
    for my $i(0..@unique-1){
	my $term = $graph->get_term($unique[$i]);   # fetch a term by ID
	push @{ $hash{$unique[$i]}{$term->name} }, $uniprot;
    }
}

#print Dumper(\%hash);
print_to_file();

sub print_to_file{
    open OUT, ">", $ARGV[2];
    print OUT "GO id", "\t", "GO term name", "\t", "ids associated", "\n";

    foreach my $go_id (keys %hash){
	foreach my $go_term (keys %{ $hash{$go_id} }){
	    my @uniprot_ids = @{$hash{$go_id}{$go_term}};
	    print OUT $go_id,"\t", $go_term, "\t";
	    print OUT join("|", @uniprot_ids);
	    print OUT "\n";
	}
    }
    close OUT;
}


sub select_db{

    my $uniprot = shift(@_);
    print "uniprot = ", $uniprot, "\n";

    my $rv = $sth->execute($uniprot) or die $DBI::errstr;
    if($rv < 0){
	print $DBI::errstr;
    }
    my @ret = ();
    while(my $annotation = $sth->fetchrow_array()) {
	push @ret, $annotation;
    }
    return @ret;
}




$dbh->disconnect();
close DATA;
