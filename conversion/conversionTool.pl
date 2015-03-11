#!/usr/bin/perl -w

use DBI;
use strict;
use warnings;

my $connection = DBI->connect('dbi:mysql:conversionDavid;mysql_socket=/tmp/mysql.sock', 'galaxy', '123456') 
    or die "Connection Error: $DBI::errstr\n";

my $input_file = $ARGV[0];
open FILE, "<", $input_file or die $!;
seek(FILE, 0, 0);
my @lines = <FILE>;

my $dbToSearch = $ARGV[1];
my $db_id;

if($dbToSearch eq "DavidGeneNameToUniprot"){
    $db_id = "david_gene_name";
} elsif($dbToSearch eq "EnsemblToUniprot_GeneID"){
    $db_id = "ensembl_geneID";
} elsif($dbToSearch eq "EnsemblToUniprot_TranscriptID"){
    $db_id = "ensembl_transcriptID";
} elsif($dbToSearch eq "IPIToUniprot"){
    $db_id = "ipi_id";
} elsif($dbToSearch eq "TairToUniprot"){
    $db_id = "tair_id";
} elsif($dbToSearch eq "TaxUniprot"){
    $db_id = "tax_id";
}

my $out = $ARGV[2];
open OUT, ">", $out or die $!;

my $sql_select_uniprotid = "select uniprot_id 
                            from $dbToSearch
                            where $db_id = ?";
my $select_uniprotid_sth = $connection->prepare($sql_select_uniprotid);


foreach (@lines) {
    my $row = $_;
    $row =~ s/^\s+|\s+$//g;

    $select_uniprotid_sth->execute($row) or die "SQL Error: $DBI::errstr\n";

    while (my $uniprot_id = $select_uniprotid_sth->fetchrow()) {
	print OUT $uniprot_id."\n";
    }

}
$select_uniprotid_sth->finish();

close FILE;

$connection->disconnect();
