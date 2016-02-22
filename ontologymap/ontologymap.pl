#!/usr/bin/perl

#teste terminal: perl ontologymap.pl Description_on_data_406_\(intensity\).csv GO test1_map.csv
#teste terminal: perl ontologymap.pl Description_on_data_406_\(intensity\).csv KEGG test1_map.csv

use DBI;
use strict;
use Data::Dumper;
use Data::Dumper::Simple;
use GO::Parser;

# regex for uniprot ID
my $regex = '(\w*)(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-?[0-9]*)$';

my $db_to_use = $ARGV[1];


#my $dbh = DBI->connect('dbi:mysql:ontology;mysql_socket=/tmp/mysql.sock', 'galaxy', '123456')
my $dbh = DBI->connect('dbi:mysql:ontology', 'galaxy', '123456')
    or die "Connection Error: $DBI::errstr\n";


# $hash{PATHWAY_ID}{PATHWAY_TERM} = [uniprot ids list]
my %hash = ();

my $select_GO = $dbh->prepare("select GO_ID from uniprot_GO_annotation where DB_OBJECT_ID = ?");
my $select_id_KEGG = $dbh->prepare("select kegg_id from KEGG_uniprot_to_kegg where uniprot_id = ?");
my $select_pathwayid_KEGG = $dbh->prepare("select pathway_id from KEGG_protein_pathway where KEGG_id = ?");
my $select_pathwayterm_KEGG = $dbh->prepare("select pathway_term from KEGG_pathway_idToTerm where pathway_id = ?");


my $graph;
if($db_to_use eq "GO"){
    my $parser = new GO::Parser({handler=>'obj'}); # create parser object
    $parser->parse("/home/mateus/galaxy/tools/galaxy_proteomics/ontologymap/go.obo"); # parse file -> objects
    $graph = $parser->handler->graph;  # get L<GO::Model::Graph> object
}

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

    my $pathway_id_to_term_ref;
    if($db_to_use eq "GO"){
	$pathway_id_to_term_ref = getpathways_GO($uniprot);	
    } else{
	$pathway_id_to_term_ref = getpathways_KEGG($uniprot);
    }

    foreach my $id (keys %{$pathway_id_to_term_ref}){
	my $term = $$pathway_id_to_term_ref{$id};
	push @{ $hash{$id}{$term} }, $uniprot;
    }

    # my @go_ids = getpathway_GO($uniprot);
    
    # my %auxhash = map{$_ => 1} @go_ids;
    # my @unique = keys %auxhash;
    # for my $i(0..@unique-1){
    # 	my $term = $graph->get_term($unique[$i]);   # fetch a GO term by GO_ID
    # 	push @{ $hash{$unique[$i]}{$term->name} }, $uniprot;
    # }
}



# print Dumper(\%hash);
print_to_file();

sub print_to_file{
    open OUT, ">", $ARGV[2];
    print OUT "$db_to_use id", "\t", "$db_to_use term name", "\t", "ids associated", "\n";

    foreach my $pathway_id (keys %hash){
	foreach my $pathway_term (keys %{ $hash{$pathway_id} }){
	    my @uniprot_ids = @{$hash{$pathway_id}{$pathway_term}};
	    print OUT $pathway_id,"\t", $pathway_term, "\t";
	    print OUT join("|", @uniprot_ids);
	    print OUT "\n";
	}
    }
    close OUT;
}


sub getpathways_GO{

    my $uniprot = shift(@_);

    my $rv = $select_GO->execute($uniprot) or die $DBI::errstr;
    if($rv < 0){
	print $DBI::errstr;
    }

    my %ret = ();
    while(my $pathway_id = $select_GO->fetchrow_array()) {
	my $term = $graph->get_term($pathway_id);   # fetch a GO term by GO_ID
	my $pathway_term = $term->name;
	$ret{$pathway_id} = $pathway_term;
    }

    return \%ret;
}

sub getpathways_KEGG{

    my $uniprot = shift(@_);

    my $rv = $select_id_KEGG->execute($uniprot) or die $DBI::errstr;
    if($rv < 0){
	print $DBI::errstr;
    }

    my %ret = ();
    while(my $kegg_id = $select_id_KEGG->fetchrow_array()){
	$rv = $select_pathwayid_KEGG->execute($kegg_id) or die $DBI::errstr;
	if($rv < 0){
	    print $DBI::errstr;
	}

	while(my $pathway_id = $select_pathwayid_KEGG->fetchrow_array()){
	    $rv = $select_pathwayterm_KEGG->execute($pathway_id) or die $DBI::errstr;
	    if($rv < 0){
		print $DBI::errstr;
	    }

	    while(my $pathway_term = $select_pathwayterm_KEGG->fetchrow_array()){
		$ret{$pathway_id} = $pathway_term;
	    }
	}

    }

    return \%ret;
}



$dbh->disconnect();
close DATA;
