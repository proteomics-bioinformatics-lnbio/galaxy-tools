#!/usr/bin/perl -w

use DBI;
use strict;
use warnings;

my $connection = DBI->connect('dbi:mysql:test', 'root', '100889') 
    or die "Connection Error: $DBI::errstr\n";

my $input_file = "phosphat_20130429.csv";
open FILE, "<", $input_file or die $!;
seek(FILE, 0, 0);
my @lines = <FILE>;

my $david_map_uniprotid = "TairToUniprot";
my $david_map_taxid = "TaxUniprot";

my @columns = split(',', $lines[0]);
my $newtable = "Phosphat";
$connection->do("DROP TABLE IF EXISTS $newtable");
my $sql_create = "CREATE TABLE $newtable (
                  id int(11) NOT NULL AUTO_INCREMENT PRIMARY KEY,
                  $columns[0] varchar(11) NOT NULL,
                  uniprot_id varchar(11) NOT NULL,
                  tax_id varchar(11) NOT NULL,
                  $columns[1] varchar(30) NOT NULL,
                  $columns[2] varchar(255) NOT NULL,
                  $columns[3] varchar(255) NOT NULL,
                  $columns[4] varchar(30) NOT NULL,
                  $columns[5] varchar(2) NOT NULL,
                  $columns[6] varchar(2) NOT NULL,
                  $columns[7] varchar(8) NOT NULL,
                  $columns[8] varchar(200) NOT NULL,
                  $columns[9] varchar(100) NOT NULL,
                  $columns[10] varchar(30) NOT NULL,
                  $columns[11] varchar(30) NOT NULL,
                  $columns[12] varchar(30) NOT NULL,
                  $columns[13] varchar(30) NOT NULL
                  )";

my $createTable_sth = $connection->prepare($sql_create);
$createTable_sth->execute or die "SQL Error: $DBI::errstr\n";
$createTable_sth->finish();


my $sql_insert = "insert into $newtable ($columns[0], uniprot_id, tax_id, $columns[1], $columns[2], 
                  $columns[3], $columns[4], $columns[5], $columns[6], $columns[7],
                  $columns[8], $columns[9], $columns[10], $columns[11], $columns[12],
                  $columns[13])
                  values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
my $insertData_sth = $connection->prepare($sql_insert);


my $sql_select_uniprotid = "select uniprot_id 
                            from $david_map_uniprotid
                            where tair_id = ?";
my $select_uniprotid_sth = $connection->prepare($sql_select_uniprotid);


my $sql_select_taxid = "select tax_id 
                        from $david_map_taxid
                        where uniprot_id = ?";
my $select_taxid_sth = $connection->prepare($sql_select_taxid);


foreach (@lines) {
    my @row_data = split(',', $_);
    $row_data[0] =~ s/^\s+|\s+$//g;
    $row_data[0] =~ s/\..*//;


    $select_uniprotid_sth->execute($row_data[0]) or die "SQL Error: $DBI::errstr\n";


    while (my $uniprot_id = $select_uniprotid_sth->fetchrow()) {

	$select_taxid_sth->execute($uniprot_id) or die "SQL Error: $DBI::errstr\n";
	my $tax_id = $select_taxid_sth->fetchrow_array();

    	$insertData_sth->execute("$row_data[0]", "$uniprot_id", "$tax_id", "$row_data[1]", "$row_data[2]",
				 "$row_data[3]", "$row_data[4]", "$row_data[5]", "$row_data[6]",
				 "$row_data[7]", "$row_data[8]", "$row_data[9]", "$row_data[10]",
				 "$row_data[11]", "$row_data[12]", "$row_data[13]") 
    	    or die "SQL Error: $DBI::errstr\n";

    }

}

close FILE;

$connection->disconnect();
