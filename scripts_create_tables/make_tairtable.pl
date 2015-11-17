#!/usr/bin/perl -w

use DBI;
use strict;
use warnings;

my $connection = DBI->connect('dbi:mysql:test', 'root', '100889') 
    or die "Connection Error: $DBI::errstr\n";
my $sql;


my $input_file = "DAVIDKnowledgebase/UNIPROT_ACCESSION2TAIR_ID.txt";
open FILE, "<", $input_file or die $!;
seek(FILE, 0, 0);
my @lines = <FILE>;

my $newtable = "TairToUniprot";
$connection->do("DROP TABLE IF EXISTS $newtable");
$sql = "CREATE TABLE $newtable (
       id int(11) NOT NULL AUTO_INCREMENT PRIMARY KEY,
       uniprot_id varchar(11) NOT NULL,
       tair_id varchar(11) NOT NULL
       )";

my $statement = $connection->prepare($sql);
$statement->execute or die "SQL Error: $DBI::errstr\n";
$statement->finish();

foreach my $row (@lines) {
    my @row_data = split('\\t', $row);
    $sql = "insert into $newtable (uniprot_id, tair_id) values (?, ?)";
    $row_data[0] =~ s/\s+$//;
    chomp($row_data[1]);

    $statement = $connection->prepare($sql);
    $statement->execute("$row_data[0]", "$row_data[1]") 
	or die "SQL Error: $DBI::errstr\n";
}

close FILE;
