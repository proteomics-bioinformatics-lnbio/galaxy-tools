#!/usr/bin/perl -w

# teste no terminal: perl description.pl /home/ABTLUS/mateus.ruivo/testes_Galaxy/teste_simples.csv saidaINFO 213 someth 9606 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 "" 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 "" id,1,gene_symbol id,2,gene_symbol intensity,3-24

# teste(com mudanca de nome de classe e amostra) no terminal: perl description.pl /home/ABTLUS/mateus.ruivo/testes_Galaxy/teste_simples.csv saidaINFO 213 someth 9606 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 id,1,gene_symbol id,2,gene_symbol intensity,3-24

# teste grande no terminal: perl description.pl testeGrande_description.csv saidaINFO 213 someth 9606 2,3,4,6,7,8,10,11,12,13,14,15,17,18,19 "" 2,3,4,6,7,8,10,11,12,13,14,15,17,18,19 "" id,1,ensembl id,5,gene_symbol id,9,tair id,16,ipi intensity,2-4 intensity,6-8 intensity,10-15 intensity,17-19


use DBI;
use strict;
use warnings;
use Data::Dumper;
use Data::Dumper::Simple;

my $connection = DBI->connect('dbi:mysql:conversionMarcelo;mysql_socket=/tmp/mysql.sock', 'galaxy', '123456') 
    or die "Connection Error: $DBI::errstr\n";

my $sql_select_synonym = "select synonyms 
                          from Synonyms2Uniprot
                          where uniprot = ?";
my $sql_select_uniprot = "select uniprot 
                          from Synonyms2Uniprot
                          where synonyms = ?";
my $sql_select_all = "select * 
                      from Synonyms2Uniprot
                      where synonyms = ?";
my $select_all_sth = $connection->prepare($sql_select_all);
my $select_synonym_sth = $connection->prepare($sql_select_synonym);
my $select_uniprot_sth = $connection->prepare($sql_select_uniprot);


my $input = $ARGV[0];
open IN, "<", $input or die $!;
seek(IN, 0, 0);
my @lines = <IN>;

my $outinfo = $ARGV[1];
open OUTINFO, ">", $outinfo or die $!;

my $tax = $ARGV[4];
my @samples_name_cols = split ',', $ARGV[5];
my @samples_new_name_cols = split ',', $ARGV[6];


my @classes_name_cols = split ',', $ARGV[7];
my @classes_new_name_cols = split ',', $ARGV[8];


my @opts = @ARGV[9..$#ARGV];

my (@intensity, @spectral_count, @fold_change, @log_ratio, @p_value) = ();
my %id_hash = ();
foreach (@opts) {
    my @opt = split ',', $_;

    # put the beginning and end columns of every input data in the respective array
    if($opt[0] eq "id"){
	push @{ $id_hash{$opt[2]} }, $opt[1]-1;
    } elsif($opt[0] eq "intensity"){
	push @intensity, [ split '-', $opt[1] ] ;
    } elsif($opt[0] eq "speccount"){
	push @spectral_count, [ split '-', $opt[1] ];
    } elsif($opt[0] eq "foldchange"){
	push @fold_change, [ split '-', $opt[1] ];
    } elsif($opt[0] eq "logratio"){
	push @log_ratio, [ split '-', $opt[1] ];
    } elsif($opt[0] eq "pvalue"){
	push @p_value, [ split '-', $opt[1] ];
    }   

}


# matrix that will have samples name in first row, classes name in second, data category
# in third, and original index in the fourth. It will help getting the values efficiently
my @input_header = ();

add_set_input_header(\@samples_name_cols, \@samples_new_name_cols, \@input_header, 0);
add_set_input_header(\@classes_name_cols, \@classes_new_name_cols, \@input_header, 1);


# set the sample, class and data category in the input_header for id's columns
# and creates a reversed id_hash use later
my %rid_hash;
while(my ($key, $value) = each %id_hash){
    my @v = @$value;
    foreach(@v){
	$input_header[0][$_] = "no_sample";
	$input_header[1][$_] = "no_class";
	$input_header[2][$_] = "id";
	$rid_hash{$_} = $key;
    }
}


# check if exist a specific data category on the input data
if(@intensity){
    add_dataCategory_input_header(\@intensity, \@input_header, "intensity");
} else{
    print OUTINFO "intensity empty\n";
}
if(@spectral_count){
    add_dataCategory_input_header(\@spectral_count, \@input_header, "spectral_count");
} else{
    print OUTINFO "spectral count empty\n";
}
if(@fold_change){
    add_dataCategory_input_header(\@fold_change, \@input_header, "fold_change");
} else{
    print OUTINFO "fold change empty\n";
}
if(@log_ratio){
    add_dataCategory_input_header(\@log_ratio, \@input_header, "log_ratio");
} else{
    print OUTINFO "log-Ratio empty\n";
}
if(@p_value){
    add_dataCategory_input_header(\@p_value, \@input_header, "p_value");
} else{
    print OUTINFO "no p-value\n";
}

# close the handles and finish connection with database
sub close_connection_and_files{
    my $val = @_;

    close IN;
    close IN;
    close OUTINFO;
    $select_all_sth->finish();
    $select_synonym_sth->finish();
    $select_uniprot_sth->finish();
    $connection->disconnect();
    exit $val;
}


# put an index in the row 3 to know the original position of every data and check the id's
my $numberofCols = @{$input_header[1]}-1;
my $error_in_id = 0;
for(0..$numberofCols){
    $input_header[3][$_] = $_;

    # check if there is any column without sample or class
    if(!defined $input_header[0][$_] || !defined $input_header[1][$_]){
	print OUTINFO "Your table must have class and sample for every data (unless id's)\n";
	close_connection_and_files(1);
    }
    # check if the first row of data has the id's coming from the source that user is telling
    if($input_header[2][$_] eq "id"){
	my @data = split '\t', $lines[2];
	my $key = $rid_hash{$_};

	if($key eq "uniprot"){
	    if($data[$_] !~ /^[OPQ][0-9][A-Z0-9]{3}[0-9] |
                             ^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
                            /x){ $error_in_id = 1; }
	} elsif($key eq "ipi"){
	    if($data[$_] !~ /^IPI[\d]+$/){ $error_in_id = 1; }
	} elsif($key eq "tair"){
	    if($data[$_] !~ /^A[Tt][1-5MC]\w[\d]+$/){ $error_in_id = 1; }
	} elsif($key eq "ensembl"){
	    if($data[$_] !~ /^ENS[\w]+$/){ $error_in_id = 1; }
	} elsif($key eq "refseq"){
	    if($data[$_] !~ /^AC_|^N[CGTWSZMR]_|^X[MR]_|^[ANYXZ]P_/){ $error_in_id = 1; }
	}
    }
    if($error_in_id) {
	print OUTINFO "Your first row has an id that does not match the original source\n";
	close_connection_and_files(1);
    }
}

#transpose a matrix
sub transpose{
    my @in = @_;
    my @out = ();
    foreach my $row(0..@in-1){
	foreach my $col(0..@{$in[$row]}-1){
	    $out[$col][$row] = $in[$row][$col];
	}
    }
    return @out;
}



# sort by data category, class and sample
@input_header = transpose
                sort { $a->[2] cmp $b->[2] or
		       $a->[1] cmp $b->[1] or
		       $a->[0] cmp $b->[0]
                     }
                transpose @input_header;


# for(0..3){
#     my $i = $_;
#     for(0..$numberofCols){
# 	print OUTINFO $input_header[$i][$_]," ";
#     }
#     print OUTINFO "\n";
# }


if(@intensity){
    write_on_file(\%id_hash, \@input_header, "intensity", $numberofCols, $ARGV[2], \@lines);
}
if(@spectral_count){
    write_on_file(\%id_hash, \@input_header, "spectral_count", $numberofCols, $ARGV[2], \@lines);
}
if(@fold_change){
    write_on_file(\%id_hash, \@input_header, "fold_change", $numberofCols, $ARGV[2], \@lines);
}
if(@log_ratio){
    write_on_file(\%id_hash, \@input_header, "log_ratio", $numberofCols, $ARGV[2], \@lines);
}


# close the files and finish connection with database
close_connection_and_files(0);


# opens and writes on file the information about specific data category
sub write_on_file{
    my($id_hash, $header, $name, $maxcol, $arg, $lines) = @_;

    my $filename = $name =~ s/_/./r;
    open my $fh, ">", "primary_".$arg."_".$filename."_visible_tabular" or die $!;

    # print the table header (sample and class)
    foreach(0..1){
	print $fh "\t";
	my @row_data = split '\t', @$lines[$_];
	my $j = 0;
	while($$header[2][$j] ne $name){ $j++; }
	while($j <= $maxcol && $$header[2][$j] eq $name){
	    chomp($$header[$_][$j]);
	    print $fh "\t", $$header[$_][$j].".$name";
	    $j++;
	}
	if($name eq "fold_change" || $name eq "log_ratio"){
	    $j = 0;
	    while($j <= $maxcol && $$header[2][$j] ne "p_value"){ $j++; }
	    while($j <= $maxcol && $$header[2][$j] eq "p_value"){
		chomp($$header[$_][$j]);
		print $fh "\t", $$header[$_][$j].".p_value";
		$j++;
	    }
	}
	print $fh "\n";
    }


    my $id_source = "";
    # search for another id if uniprot accession is not in the table
    if(!defined $$id_hash{"uniprot"}){
	if(defined $$id_hash{"tair"}){ $id_source = "tair"; }
	elsif(defined $$id_hash{"ensembl"}){ $id_source = "ensembl"; }
	elsif(defined $$id_hash{"ipi"}){ $id_source = "ipi"; }
	elsif(defined $$id_hash{"refseq"}){ $id_source = "refseq"; }
	else { $id_source = "gene_symbol"; }
    }
    
    # iterates through the input table data
    foreach(2..$#$lines){
	my @row_data = split '\t', @$lines[$_];

	my $id_found = $row_data[$$id_hash{$id_source}[0]];
	my $gene_symbol = "";
	my $uniprot_id = ($id_source eq "uniprot") ? $id_found : "";

	# if the only id in the table is gene_symbol we have to find the specific uniprot
	if($id_source eq "gene_symbol"){
	    $gene_symbol = $id_found;
	    $select_all_sth->execute($id_found) or die "SQL Error: $DBI::errstr\n";
	  FOUND_UNIPROT: while(my @results = $select_all_sth->fetchrow()){
	      # if the tax_id is the one that the user is looking for and the id is reviewed, then we found it
	      if($results[4] eq $tax){
		  $uniprot_id = $results[2];
		  if($results[5] eq "YES"){ last FOUND_UNIPROT; }
	      }
	  }
	    # in case of uniprot not found
	    # if($gene_symbol eq $id_found){ $id_found = ""; }
	} else{
	    # if the id is not from uniprot, try to find the mapping to uniprot
	    if($id_source ne "uniprot"){
		$select_uniprot_sth->execute($id_found) or die "SQL Error: $DBI::errstr\n";
	      UNIPROT_MAP: while(my $un = $select_uniprot_sth->fetchrow()){
		  if($un =~ /^[OPQ][0-9][A-Z0-9]{3}[0-9] |
                             ^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
                            /x){
		      $uniprot_id = $un;
		      last UNIPROT_MAP;
		  }
	      }
	    }

	    my @possible_gene_symbol = ();
	    # if the got id is from uniprot, search for the gene symbol
	    if($uniprot_id ne ""){
		$uniprot_id =~ s/-[0-9]//;
		$select_synonym_sth->execute($uniprot_id) or die "SQL Error: $DBI::errstr\n";
	      GENESYMBOL_SEARCH: while (my $gs = $select_synonym_sth->fetchrow()) {
		  push @possible_gene_symbol, $gs;
		  if ($gs !~ /^IPI|^ENS|^A[Tt]/){
		      $gene_symbol = $gs;
		      last GENESYMBOL_SEARCH;
		  }
	      }
	    }
	    
	    if($gene_symbol eq ""){
		# in case of gene_symbol not found until here, pick one in the list of possibilities with the 
		# priority of start with AT, ENS, IPI. Thats why the use of sort.
		if(@possible_gene_symbol){
		    @possible_gene_symbol = sort @possible_gene_symbol;
		    $gene_symbol = $possible_gene_symbol[0];
		}
	    }
	}
	# prints the first and second table columns
	print $fh $id_found,"\t",$gene_symbol."_".$uniprot_id;


	# since the header is sorted by data category, class, sample; once we found
	# the first column with the name all the next columns will have same data category
	# until it changes
	my $j = 0;
	while($$header[2][$j] ne $name){ $j++; }
	while($j <= $maxcol && $$header[2][$j] eq $name){
	    chomp($row_data[$$header[3][$j]]);
	    print $fh "\t", $row_data[$$header[3][$j]];
	    $j++;
	}
	# in case of Fold Change or Log_Ratio, we have to add the P-value columns if they exist
	if($name eq "fold_change" || $name eq "log_ratio"){
	    $j = 0;
	    while($j <= $maxcol && $$header[2][$j] ne "p_value"){ $j++; }
	    while($j <= $maxcol && $$header[2][$j] eq "p_value"){
		chomp($row_data[$$header[3][$j]]);
		print $fh "\t", $row_data[$$header[3][$j]];
		$j++;
	    }	    
	}
	print $fh "\n";
    }

    close $fh;
}

# set the columns on the header matrix with sample or class
sub add_set_input_header{
    my ($vec_cols, $vec_newnames, $header, $line) = @_;

    my @l = split '\t', $lines[$line];

    if(@{$vec_newnames}){
	foreach my $v (@{$vec_cols}) {
	    $$header[$line][$v-1] = shift(@{$vec_newnames});
	}
    } else {
	foreach my $v (@{$vec_cols}) {
	    $$header[$line][$v-1] = $l[$v-1];
	}
    }

}


# set the columns on the header matrix to this data category
sub add_dataCategory_input_header{
    my ($data_category, $header, $name) = @_;

    foreach my $a (@$data_category) {
	my $beg = @$a[0]-1;
	my $end = @$a[1]-1;
	if($end < $beg){ my $aux = $beg; $beg = $end; $end = $aux; }
	while($beg <= $end){
	    $$header[2][$beg] = $name;
	    $beg++;
	}
    }

}
