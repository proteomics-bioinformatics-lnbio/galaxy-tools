#!/usr/bin/perl -w

# test flavia in terminal: perl pes.pl primary_213_intensity_visible_tabular.csv mapping.csv 0.5 100 2 out.csv


use List::Util qw(min max);
use Data::Dumper;
use Data::Dumper::Simple;
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use Data::Compare;

use Statistics::Multtest qw(bonferroni holm hommel hochberg BH BY qvalue);
use Statistics::Multtest qw(:all);
use Storable qw(dclone);

# CONSTANTS
my $EPS = 0.00000001;
my $regex = '(\w*)(([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})-?[0-9]*)$';
my $PVALUE = 0.05;
my $QVALUE = 0.005;

open DATA, $ARGV[0]; #file with protein data
open GOTERM, $ARGV[1]; #file with map protein id to GO id

my $f = $ARGV[2]; # tolerance factor
my $numberOfDistributions = $ARGV[3]; # number of null distributions
my $min_proteins_GO = $ARGV[4]; # minimun number of proteins in GO term

if($f > 1 || $f < 0){
    print STDERR "f must be in interval [0,1]\n";
    exit;
}

open OUT, ">", $ARGV[5]; # file with output
open OUTPVALUES, ">", $ARGV[6]; # file with p-values for future heatmap


my %hash_mean; # protein -> mean
my %hash_CV; # protein -> CVp
my %hash_W; # protein -> W
my %hash_ESW; # protein -> enrichment score weight

my %hash_pvalue; # GO_ID -> pvalue
my %hash_qvalue; # GO_ID -> qvalue
my %hash_PES; # GO_ID -> PES
my %hash_PESratio; # GO_ID -> PES_condA/PES_condB
my %hash_GO_Core; # GO_ID -> proteins Core
my %hash_nullInsideRange; # GO_ID -> random null distributions inside of range (1-f)C <= C' <= (1+f)C
my %hash_GO = createGOhash(); # GO_ID -> proteins
my %hash_GOid_GOterm = create_GOid_GOterm(); # GO_ID -> GO_Term


my @ids = get_ids();
my ($ref_hash, $ref_vec) = get_ttest();
my %hash_ttest = %$ref_hash;
my @ttests = @$ref_vec;
separateConditions();


close(DATA);
close(GOTERM);
close(OUT);
close(OUTPVALUES);


# generate a random permutation of @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
	my $j = int rand ($i+1);
	next if $i == $j;
	@$array[$i,$j] = @$array[$j,$i];
    }
}

# get ttest values from table and creates hash and array
sub get_ttest{
    my %hash = ();
    my @vec = ();

    seek DATA, 0, 0;
    my $i = 0;
    while(<DATA>){
	$i++;
	if($i <= 2){ next; }
	my @col = split '\t', $_;
	my $ttest_val = pop @col;
	chomp $ttest_val;
	$hash{$ids[$i-3]} = $ttest_val;
	push @vec, $ttest_val;
    }
    @vec = sort{$b <=> $a} @vec;

    return (\%hash, \@vec);
}

# get ids from table and creates array
sub get_ids{
    my @ids = ();

    seek DATA, 0, 0;
    my $i = 0;
    while(<DATA>){
	if($i < 2){ $i++; next; }
	my @col = split '\t', $_;
	my @aux = split($regex, $col[1]);
	my $id = $aux[2];
	push @ids, $id;
    }
    return (@ids);
}

# create means and CVps sorted (inverse) array of certain condition. Useful to get the max
sub create_mean_and_CVp_array{
    my $condition = shift;
    my @means = values %{$hash_mean{$condition}};
    @means = sort{$b <=> $a} @means;

    my @CVps = values %{$hash_CV{$condition}};
    @CVps = sort{$b <=> $a} @CVps;

    return (\@means, \@CVps);
}

# separates data between conditions and call functions to create hash of each value that must be calculated to generate
# score
sub separateConditions{
    my ($ref_hash, $ref_array);

    seek DATA, 0, 0;
    my @lines = <DATA>;
    shift(@lines);
    my @classes = split '\t', shift @lines;
    chomp @classes;

    my @cond = @classes;
    shift @cond; shift @cond; pop @cond;
    @cond = uniq @cond;

    foreach(@cond){
	my $current_cond = $_;
	my @condition = ();
	for my $i(0..@lines-1){
	    my @data = split '\t', $lines[$i];
	    chomp(@data);
	    $condition[$i][0] = $data[0];
	    $condition[$i][1] = $data[1];

	    my $j = 2;
	    my $k = 2;
	    my $sz = @data;
	    while($j < $sz && ($classes[$j] ne $current_cond)) { $j++; }
	    while(1){
		if($j >= $sz || ($classes[$j] ne $current_cond)){ last; }
		$condition[$i][$k++] = $data[$j++];
	    }
	}

	$ref_hash = calculate_means(\@condition);
	$hash_mean{$current_cond}{$_} = $$ref_hash{$_} foreach keys %$ref_hash;

	$ref_hash = calculate_CVp(\@condition, $current_cond);
	$hash_CV{$current_cond}{$_} = $$ref_hash{$_} foreach keys %$ref_hash;

	my($ref_means, $ref_CVps) = create_mean_and_CVp_array($current_cond);

	$ref_hash = calculate_W($current_cond, $$ref_means[0], $$ref_CVps[0]);
	$hash_W{$current_cond}{$_} = $$ref_hash{$_} foreach keys %$ref_hash;

	$ref_hash = calculate_ESW($current_cond);
	$hash_ESW{$current_cond}{$_} = $$ref_hash{$_} foreach keys %$ref_hash;

    }
    calculate_ratioPES(\@cond);
    print_to_file(\@cond);
}

# I think it has to be the matrix described in the paper
sub calculate_ESW{
    my $condition = shift;

    my %hash = ();
    while(my($k, $v) = each %{$hash_W{$condition}}){
	$hash{$k} = $v;
    }
    return (\%hash);
}

# function that assign a weight to a given protein according to a specific formula
sub calculate_Weight{
    my $condition = shift;
    my $uniprot = shift;
    my $max_mean = shift;
    my $max_CVp = shift;
    my $max_ttest = shift;

    my $CVpt = $hash_CV{$condition}{$uniprot}*$max_mean / $max_CVp;
    my $ttestt = $hash_ttest{$uniprot}*$max_mean / $max_ttest;

    # Score from Yates paper
    # my $weight = ($hash_mean{$condition}{$uniprot} - $CVpt + $max_mean - 1 )/(2* $max_mean - 2);

    # Score modified by us
    my $weight = ($hash_mean{$condition}{$uniprot} - 0.5*$CVpt - 0.5*$ttestt + $max_mean - 1 )/(2* $max_mean - 2);
 
    return $weight;
}

# create hash of weight to each protein in a given condition
sub calculate_W{
    my $condition = shift;
    my $max_mean = shift;
    my $max_CVp = shift;

    my $max_ttest = $ttests[0];
 
    my %hash = ();
    foreach(@ids){
	my $uniprot = $_;
	$hash{$uniprot} = calculate_Weight($condition, $uniprot, $max_mean, $max_CVp, $max_ttest);
    }
    return (\%hash);
}

# create hash of means to a given condition
sub calculate_means{
    my $data = shift;

    my %hash = ();
    my $i = 0;
    foreach(@{$data}){
	my @vec = @$_;
	shift(@vec); shift(@vec);
	my $uniprot = $ids[$i++];

	@vec = map{ $_ =~ s/[Nn][Aa][Nn]/0/g; $_ } @vec;
	my $sum = sum(@vec);
	$sum = ($sum == 0) ? 0.00001 : $sum;

	my $sz = @vec;
	$hash{$uniprot} = $sum/$sz;
    }
    return (\%hash);
}

# create hash of CVp to a given condition
sub calculate_CVp{
    my $data = shift;
    my $condition = shift;

    my %hash = ();
    my $i = 0;
    foreach(@{$data}){
	my @vec = @$_;
	shift(@vec); shift(@vec);
	my $uniprot = $ids[$i++];

	@vec = map{ $_ =~ s/[Nn][Aa][Nn]/0/g; $_ } @vec;

	my $mu = $hash_mean{$condition}{$uniprot};
	$hash{$uniprot} = stddev(@vec)/$mu;
    }
    return (\%hash);
}

sub printdata{
    my $data = shift;

    foreach(@{$data}){
	print Dumper(@$_);
    }

}

sub normalize{
    my $data = shift;

    foreach(@{$data}){
	my @vec = @$_;
	shift(@vec); shift(@vec);

	@vec = map{ $_ =~ s/[Nn][Aa][Nn]/0/g; $_ } @vec;
	my $sum = sum(@vec);
	$sum = ($sum == 0) ? 0.00001 : $sum;

	for my $i(0..@vec-1){
	    $$_[$i+2] /= $sum;
	}
    }
}

# calculate PES null distribution for a specific number of proteins (n)
sub calculateNullDistribution{
    my $n = shift;
    my $originalSetCV = shift;
    my $condition = shift;
    my $GO_term = shift;
    

    my (@acc, @rej, @nullDistribution) = ();
    my $j = 0;
    my $countloop = 0;
    while(1){
	my @list = @ids;
	# create random Set
	fisher_yates_shuffle( \@list );
    	my @set = @list[0..$n-1];
	$countloop++;

	if($j == $numberOfDistributions){
	    push @nullDistribution, $_ foreach @acc;
	    $hash_nullInsideRange{$condition}{$GO_term} = 0;
	    last;
	}

	if($countloop > 3*$numberOfDistributions){
	    push @nullDistribution, $_ foreach @acc;
	    my $sz = @nullDistribution;
	    my $left = $numberOfDistributions - $sz;
	    $hash_nullInsideRange{$condition}{$GO_term} = $left;


	    while($left--){
		push @nullDistribution, shift @rej;
	    }
	    $sz = @nullDistribution;
	    last;
	}

    	# check if this set is acceptable and calculate its PES
    	my $CV = calculate_proteinSetCV(\@set, $condition);
    	if($CV >= (1-$f)*$originalSetCV && $CV <= (1+$f)*$originalSetCV){
	    my $PES = 0;
	    foreach(@set){
		$PES += $hash_ESW{$condition}{$_};
	    }
	    push @acc, $PES;
    	    $j++;
	# rejected ones are stored
    	} else{
	    my $PES = 0;
	    foreach(@set){
		$PES += $hash_ESW{$condition}{$_};
	    }
	    push @rej, $PES;
	}
    }

    return @nullDistribution;
}

# calculate abundance means CV of a given set of proteins
sub calculate_proteinSetCV{
    my $ids = shift;
    my $condition = shift;

    # make the mean array
    my @means = ();
    foreach (@{$ids}){
	push @means, $hash_mean{$condition}{$_};
    }
    my $mean = mean(@means);
    $mean = ($mean == 0) ? 0.00001 : $mean;

    my $CV = stddev(@means) / $mean;

    return $CV;
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

# delete element from array
sub delete_from_array{
    my $vec_ref = shift;
    my $el = shift;

    my $i = 0;
    $i++ until $$vec_ref[$i] eq $el;
    splice(@$vec_ref, $i, 1);
}

# look for Core in a given protein set. Continue to removing the protein with lower weight that leads to a better
# p-value. If new p-value >= last p-value stop
sub find_CoreS{
    my $p_value = shift;
    my $ref_idsSet = shift;
    my $condition = shift;
    my $ref_means = shift;
    my @means = @$ref_means;
    my $ref_CVps = shift;
    my @CVps = @$ref_CVps;
    my $GO_term = shift;

    my ($newPES, $retPES);
    my $new_pvalue = $p_value;
    my @set_core = @$ref_idsSet;


    my $index_max_mean = 0;
    my $index_max_CVp = 0;
    my $index_max_ttest = 0;

    my %hash_Wset;
    foreach(@{$ref_idsSet}){
	$hash_Wset{$_} = $hash_ESW{$condition}{$_};
    }

    my @ret_nullDistribution = ();
    foreach my $uniprot (sort {$hash_Wset{$a} <=> $hash_Wset{$b}} keys %hash_Wset){
	if(@set_core == 1){ last; }
	delete_from_array(\@set_core, $uniprot);

    	if($hash_mean{$condition}{$uniprot} == $means[$index_max_mean]){
    	    $index_max_mean++;
    	}
    	if($hash_CV{$condition}{$uniprot} == $CVps[$index_max_CVp]){
    	    $index_max_CVp++;
    	}
	if($hash_ttest{$uniprot} == $ttests[$index_max_ttest]){
	    $index_max_ttest++;
	}

    	my $max_mean = $means[$index_max_mean];
    	my $max_CVp = $CVps[$index_max_CVp];
	my $max_ttest = $ttests[$index_max_ttest];
    	$newPES = 0;
	foreach(@set_core){
	    $newPES += calculate_Weight($condition, $_, $max_mean, $max_CVp, $max_ttest);
    	}

	my $CV = calculate_proteinSetCV(\@set_core, $condition);
	my @PES_randomSets = calculateNullDistribution($#set_core+1, $CV, $condition, $GO_term);

	my @bigger = grep {$_ > $newPES} @PES_randomSets;
	my $num = @bigger;
	$new_pvalue = $num/$numberOfDistributions;
	if($new_pvalue >= $p_value){
	    push @set_core, $uniprot;
	    last;
	}
	$p_value = $new_pvalue;
	$retPES = $newPES;
	@ret_nullDistribution = @PES_randomSets;

    }
    return ($retPES, $p_value, \@set_core, \@ret_nullDistribution);
}

# calculate PES for a given protein set, get a null distribution and calculate the p-value
sub calculate_PES{
    my $ids_ref = shift;
    my @ids_GOset = @$ids_ref;

    my $GO_term = shift;
    my $condition = shift;
    my $ref_means = shift;
    my $ref_CVps = shift;

    my $PES = 0;
    $PES += $hash_ESW{$condition}{$_} foreach @ids_GOset;
    $hash_PES{$condition}{$GO_term} = $PES;

    my $CV = calculate_proteinSetCV(\@ids_GOset, $condition);
    my @PES_randomSets = calculateNullDistribution($#ids_GOset+1, $CV, $condition, $GO_term);


    # get number of elements in null distribution that are greater than original value -> p-value
    my @bigger = grep {$_ > $PES} @PES_randomSets;
    my $num = @bigger;
    $hash_pvalue{$condition}{$GO_term} = $num/$numberOfDistributions;

    if($hash_pvalue{$condition}{$GO_term} <= $PVALUE){
	# Look for Core(protein set S)
	my ($newPES, $new_pvalue, $ref_core, $ref_nullDistribution) = find_CoreS($hash_pvalue{$condition}{$GO_term}, \@ids_GOset, $condition, $ref_means, $ref_CVps, $GO_term);

	# print "new_pvalue = ", $new_pvalue, " pvalue = ", $hash_pvalue{$condition}{$GO_term}, "\n";
	if($new_pvalue < $hash_pvalue{$condition}{$GO_term}){
	    $hash_PES{$condition}{$GO_term} = $newPES;
	    $hash_GO_Core{$condition}{$GO_term} = \@$ref_core;
	    $hash_pvalue{$condition}{$GO_term} = $new_pvalue;
	    @PES_randomSets = @$ref_nullDistribution;
	}
    }
    return(\@PES_randomSets, $hash_PES{$condition}{$GO_term});    
}

# calculate pvalue of ratio between PES's for a given protein set. Uses the null distributions calculated in 
# each condition
sub calculate_pvalue_ratio{
    my $ratio_condAB = shift;
    my $ref_nullDistributionA = shift;
    my $ref_nullDistributionB = shift;

    my @nullDistribution_ratio = ();
    my $range = $numberOfDistributions;
    for my $j (0..$numberOfDistributions){
	my $random_number1 = int(rand($range));
	my $random_number2 = int(rand($range));

	push @nullDistribution_ratio, $$ref_nullDistributionA[$random_number1] / $$ref_nullDistributionB[$random_number2];
    }
    my @bigger = grep {$_ > $ratio_condAB} @nullDistribution_ratio;
    my $num = @bigger;

    return $num/$numberOfDistributions;
}

sub calculate_qvalue{
    my $condition = shift;

    my %hash = ();
    $hash{$_} = $hash_pvalue{$condition}{$_} foreach keys %{$hash_pvalue{$condition}};

    my @seq = ();
    my $val = 0.0;
    while($val < 0.95){
	push @seq, $val;
	$val+=0.05;
    }

    my %setup = (
	lambda => \@seq,
	robust => 0,
	);
    my $qval = qvalue(\%hash, %setup);

    while(my($k,$v) = each %$qval){
	$hash_qvalue{$condition}{$k} = $v;
    }

}

# call function to calculate PES of a given protein set in two conditions. Then take the ratio between the two PES's
# calculated and call function to calculate the p-value ratio
sub calculate_ratioPES{
    my $ref_conditions = shift; # ref array with two conditions

    my($ref_meansA, $ref_CVpsA) = create_mean_and_CVp_array($$ref_conditions[0]);
    my($ref_meansB, $ref_CVpsB) = create_mean_and_CVp_array($$ref_conditions[1]);

    my $already_calculated = 0;
    my $size = keys %hash_GO;
    while(my($key, $value) = each %hash_GO){
	# print_percentage_done($already_calculated, $size);
	$already_calculated++;

	my @ids_GOset = @{$value};
	if(@ids_GOset < $min_proteins_GO){ next; }

	my($ref_nullDistributionA, $PES_condA) = calculate_PES(\@ids_GOset, $key, $$ref_conditions[0], $ref_meansA, $ref_CVpsA);
	my($ref_nullDistributionB, $PES_condB) = calculate_PES(\@ids_GOset, $key, $$ref_conditions[1], $ref_meansB, $ref_CVpsB);

	$hash_PESratio{$key} = $PES_condA/$PES_condB;
	my $pvalue_ratio = calculate_pvalue_ratio($hash_PESratio{$key}, $ref_nullDistributionA, $ref_nullDistributionB);
	$hash_pvalue_ratio{$key} = $pvalue_ratio;
    }
    calculate_qvalue($$ref_conditions[0]);
    calculate_qvalue($$ref_conditions[1]);
}

sub print_percentage_done(){
    my $already_calc = shift;
    my $size = shift;
    printf "$already_calc $size: %.3f %% GO terms analyzed.\n", $already_calc/$size * 100;
}

sub print_to_file{
    my $ref_conditions = shift;

    print OUT "DATA IN SINGLE CONDITION:\n";
    foreach(@{$ref_conditions}){
	my $cond = $_;
	my $count = keys %{$hash_pvalue{$cond}};
	my @aux = values %{$hash_pvalue{$cond}};
	my $sig = grep {$_ <= $PVALUE} @aux;
	print OUT $cond,"\n";
	print OUT "% p-value <= $PVALUE: ", $sig/$count *100, " %\n";
	print OUT "GO_id", "\t", "GO_term", "\t", "proteins", "\t", "Core", "\t", "p-value", "\t", "q-value", "\t", "PES", "\t", "Number of null distributions inside range", "\n";
	while(my($k,$v) = each %{$hash_pvalue{$cond}}){
	    if($v <= $PVALUE){
		my $protein_list = join("|", @{$hash_GO{$k}});
		my $protein_list_core = (defined $hash_GO_Core{$cond}{$k}) ? join("|", @{$hash_GO_Core{$cond}{$k}}) : "Core not found";

		print OUT $k, "\t", $hash_GOid_GOterm{$k}, "\t", $protein_list, "\t", $protein_list_core,
		"\t", $v, "\t", $hash_qvalue{$cond}{$k}, "\t", $hash_PES{$cond}{$k},"\t",$hash_nullInsideRange{$cond}{$k},"\n";
	    }
	}
	print OUT "\n";
    }

    print OUT "DATA RATIO BETWEEN CONDITIONS:\n";
    my $count = keys %hash_pvalue_ratio;
    my $sig = grep {$_ <= $PVALUE} values %hash_pvalue_ratio;
    print OUT "% p-value <= $PVALUE: ", $sig/$count *100, " %\n";
    print OUT "GO_term", "\t", "p-value", "\t", "PES\n";
    while(my($k,$v) = each %hash_pvalue_ratio){
	if($v <= $PVALUE){
	    print OUT $k, "\t", $v, "\t", $hash_PESratio{$k},"\n";
	}
    }

    my $cond1 = $$ref_conditions[0] =~ s/([-_\w]+).([-_\w]+)/$1/r;
    my $cond2 = $$ref_conditions[1] =~ s/([-_\w]+).([-_\w]+)/$1/r;

    print OUTPVALUES "\t", $cond1, "\t", $cond2,"\n";
    while(my($goid, $pv1) = each %{$hash_pvalue{$$ref_conditions[0]}}){
	my $pv2 = $hash_pvalue{$$ref_conditions[1]}{$goid};
	if($pv1 <= $PVALUE && $pv2 <= $PVALUE){
	    print OUTPVALUES $goid,"\t", $pv1,"\t", $pv2,"\n";
	}
    }


}

# create a hash of GO terms map
sub createGOhash{
    my %hash;

    my $i = 0;
    seek GOTERM, 0, 0;
    foreach(<GOTERM>){
	if($i < 1){
	    $i++;
	    next;
	}
	chomp $_;
	my @data = split '\t', $_;
	push @{ $hash{$data[0]} }, split '\|', $data[2];
    }
    return %hash;
}

# create a hash of GO_id to GO_term
sub create_GOid_GOterm{
    my %hash;

    my $i = 0;
    seek GOTERM, 0, 0;
    foreach(<GOTERM>){
	if($i < 1){
	    $i++;
	    next;
	}
	chomp $_;
	my @data = split '\t', $_;
	$hash{$data[0]} = $data[1];
    }
    return %hash;
}

# calculate mean of array
sub mean{
    return sum(@_)/@_;
}

# calculate standard deviation of array
sub stddev{
    my $mean = mean(@_);
    my @v = ();
    foreach(@_){
	push @v, ($_ - $mean)**2;
    }
    my $sz = @v;
    return sqrt(sum(@v)/$sz);
}

    # My implementation of q-value
    # my $sz = keys %{$hash_pvalue{$condition}};
    # my $i = 1;
    # my $cuttof_found = 0;
    # my %hash;
    # foreach my $GO_term(sort{$hash_pvalue{$condition}{$a} <=> $hash_pvalue{$condition}{$b}} keys %{$hash_pvalue{$condition}}){
    # 	if($i == $sz || ($hash_pvalue{$condition}{$GO_term} > ($i * $ALPHA)/$numberOfDistributions)){ $cuttof_found = 1; }
	
    # 	$hash_qvalue{$condition}{$GO_term} = ($i * $ALPHA)/$numberOfDistributions;

    # 	# if($cuttof_found){
    # 	#     $hash_qvalue{$condition}{$GO_term} = '';
    # 	# } else{
    # 	#     $hash_qvalue{$condition}{$GO_term} = '+';
    # 	# }

    # 	$hash{$GO_term} = $hash_pvalue{$condition}{$GO_term};
    # 	$i++;
    # }
