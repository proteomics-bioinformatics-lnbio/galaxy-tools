#!/usr/bin/perl -w

# pra rodar no terminal: perl aminomodification.pl protein_ids_tester.csv "mod_res,disulf" "None" "S-nitrosocysteine" "None" "None" "None" out_file anything

use strict;
use warnings;
use LWP::UserAgent;
use XML::Twig;
use FileHandle;

# open file with the ids and get the lines #
my $in = $ARGV[0];
open INFILE, "<", $in or die $!;
seek(INFILE, 0, 0);

my $agent = LWP::UserAgent->new;
$agent->env_proxy;

# this file will have the ID's that don't exist in Uniprot
open OUT1, ">", $ARGV[7] or die $!;

# get the features key to be searched and put them into an array 
my @features_key = split ',', $ARGV[1];

# this hash contains the post-translational modifications types as keys and its content as values
my %features;
foreach my $feat (@features_key){
    if($feat eq "nonstd_res"){ $features{"nonstd_res"} = [ split ',', $ARGV[2] ]; }
    if($feat eq "mod_res"){ $features{"mod_res"} = [ split ',', $ARGV[3] ]; }
    if($feat eq "lip"){ $features{"lip"} = [ split ',', $ARGV[4] ]; }
    if($feat eq "glyco"){ $features{"glyco"} = [ split ',', $ARGV[5] ]; }
    if($feat eq "disulf"){ $features{"disulf"} = [ "all" ]; }
    if($feat eq "cross"){ $features{"cross"} = [ split ',', $ARGV[6] ]; }
}

my @lines = <INFILE>;

# array that stores all output information
my @output;

foreach (@lines) {
    # hash for each protein in the input list
    my %protein;

    my $uniprot_id = $_;
    chomp($uniprot_id);
    $uniprot_id =~ s/\r$//;
    my $url = "http://www.uniprot.org/uniprot/" . $uniprot_id . ".xml";
    my $response = $agent->get($url);

    $protein{PROTEIN_ID} = $uniprot_id;
    
    $protein{PROTEIN_NAME} = "";
    $protein{MODIFICATION_LIST} = [];
    if ($response->is_success && $response->content) {

	my $t = XML::Twig->parse($response->content);

	# get the name of protein
	my $node_rec = $t->first_elt('protein')->first_child('recommendedName');
	my $node_sub = $t->first_elt('protein')->first_child('submittedName');
	my $node_alt = $t->first_elt('protein')->first_child('alternativeName');
	my $protein_name;
	if ($node_rec) {
	    $protein_name = $node_rec->first_child_text('fullName');
	} elsif ($node_sub) {
	    $protein_name = $node_sub->first_child_text('fullName');
	} elsif ($node_alt) {
	    $protein_name = $node_alt->first_child_text('fullName');
	}

	searchInUniprot(\%protein, $uniprot_id, $protein_name, $t);
    } else {
	print OUT1 "Cannot find id $uniprot_id in Uniprot\n";
    }
    push @output, \%protein;

}	

printHashToFile();


sub printHashToFile{

    my $header = "uniprot_id\tprotein_name\tmodification_position\tmodification_description\treference\n";

    my ($fh_nonres, $fh_modres, $fh_lip, $fh_glyco, $fh_disulf, $fh_cross) = ();
    foreach my $key (keys %features) {

	if($key eq "nonstd_res"){
	    $fh_nonres = FileHandle->new(">primary_".$ARGV[8]."_"."NonStandardRes_visible_tabular");
	    print $fh_nonres $header;
	} elsif($key eq "mod_res"){
	    $fh_modres = FileHandle->new(">primary_".$ARGV[8]."_"."ModifiedRes_visible_tabular");
	    print $fh_modres $header;
	} elsif($key eq "lip"){
	    $fh_lip = FileHandle->new(">primary_".$ARGV[8]."_"."Lipitadtion_visible_tabular");
	    print $fh_lip $header;
	} elsif($key eq "glyco"){ 
	    $fh_glyco = FileHandle->new(">primary_".$ARGV[8]."_"."Glycosylation_visible_tabular");
	    print $fh_glyco $header;
	} elsif($key eq "disulf"){ 
	    $fh_disulf = FileHandle->new(">primary_".$ARGV[8]."_"."DisulfideBond_visible_tabular");
	    print $fh_disulf $header;
	} elsif($key eq "cross"){
	    $fh_cross = FileHandle->new(">primary_".$ARGV[8]."_"."CrossLink_visible_tabular");
	    print $fh_cross $header;
	}

    }

    foreach my $p(@output){
	my %p_hash = %{$p};

	my @mod_list = @{$p_hash{MODIFICATION_LIST}};
	foreach my $mod(@mod_list){
	    my %h = %{$mod};

	    my $fh;
	    if($h{MODIFICATION} eq "non-standard amino acid"){
		$fh = $fh_nonres;
	    } elsif($h{MODIFICATION} eq "modified residue"){
		$fh = $fh_modres;
	    } elsif($h{MODIFICATION} eq "lipid moiety-binding region"){
		$fh = $fh_lip;
	    } elsif($h{MODIFICATION} eq "glycosylation site"){
		$fh = $fh_glyco;
	    } elsif($h{MODIFICATION} eq "disulfide bond"){
		$fh = $fh_disulf;
	    } elsif($h{MODIFICATION} eq "cross-link"){
		$fh = $fh_cross;
	    }

	    my @feat_list = @{$h{FEATURE_DATA}};
	    foreach my $feat(@feat_list){
		my %f = %{$feat};
		my @ref_list = @{$f{all_refs}};

		if(!@ref_list){
		    print $fh "$p_hash{PROTEIN_ID}\t$p_hash{PROTEIN_NAME}\t$f{modification_position}\t$f{modification_description}\t\n";
		}
		foreach my $ref(@ref_list){
		    print $fh "$p_hash{PROTEIN_ID}\t$p_hash{PROTEIN_NAME}\t$f{modification_position}\t$f{modification_description}\t$ref\n";
		}
	    }

	}

    }


    if(defined $fh_nonres){ $fh_nonres->close; }
    if(defined $fh_modres){ $fh_modres->close; }
    if(defined $fh_lip){ $fh_lip->close; }
    if(defined $fh_glyco){ $fh_glyco->close; }
    if(defined $fh_disulf){ $fh_disulf->close; }
    if(defined $fh_cross){ $fh_cross->close; }
}

sub searchInUniprot {
    my ($modification_position, $description, $ps) = ();

    my ($prot_hashref, $uniprot_id, $protein_name, $t) = @_;

    $$prot_hashref{PROTEIN_NAME} = $protein_name;

    my @modif_data = ();
    foreach my $key (keys %features) {
	my %modif_hash;
	my $modif_type;
	if($key eq "nonstd_res"){
	    $modif_type = "non-standard amino acid";
	} elsif($key eq "mod_res"){
	    $modif_type = "modified residue";
	} elsif($key eq "lip"){
	    $modif_type = "lipid moiety-binding region";
	} elsif($key eq "glyco"){ 
	    $modif_type = "glycosylation site";
	} elsif($key eq "disulf"){ 
	    $modif_type = "disulfide bond";
	} elsif($key eq "cross"){
	    $modif_type = "cross-link";
	}
	$modif_hash{MODIFICATION} = $modif_type;

	my @feat_data = ();
	foreach my $feat_search (@{$features{$key}}){

	    my ($modification_position, $modification_description, $ps) = ();

	    $feat_search =~ s/X__sq__|__sq__/'/g;
	    $feat_search =~ s/X/,/g;
	    my $name = $feat_search =~ s/\s|'/./rg;
	    chomp($feat_search);
	    chomp($name);

	    # For disulfide bond, the search will be by type and not by description
	    my @nodes = ();
	    if($feat_search eq "all") {
		@nodes = $t->findnodes("//entry/feature[\@type =~ /$modif_type/]");
	    } else{
		@nodes = $t->findnodes("//entry/feature[\@description =~ /$feat_search/]");
	    }


	    my @node_data = ();
	    foreach my $feature (@nodes) {
		my %hash_node;

		my $position = $feature->first_child('location')->first_child('position');
		if($position){
		    $modification_position = $position->att('position');
		} else {
		    $modification_position = $feature->first_child('location')->first_child('begin')->att('position')
			."-". $feature->first_child('location')->first_child('end')->att('position');
		}
		$hash_node{modification_position} = $modification_position;

		$modification_description = $feature->att('description') ? $feature->att('description') : '';
		$hash_node{modification_description} = $modification_description;
		

		my @refs = ();
		my @evs = ();
		# Array with all references and evidences, or the status if there isnt't any references
		my @all_refs = ();

		my $status = $feature->att('status');
		if (!$status) {
		    # take the number of evidences if there's any #
		    my $evidence = $feature->att('evidence');
		    if ($evidence) {
			my @numbers = split ' ', $evidence;

			my $pubmed_link = "http://www.ncbi.nlm.nih.gov/pubmed/";
			for (@numbers) {
			    my $e;
			    my $dbReference = $t->first_elt("evidence[\@key=$_]");
			    if($dbReference->first_child('source')){
				$dbReference = $dbReference->first_child('source')->first_child('dbReference');
				my $type = "";
				my $id = "";
				if($dbReference){
				    $type = $dbReference->att('type');
				    $id = $dbReference->att('id');
				}
				if($type eq "PubMed"){
				    $e = $pubmed_link.$id;
				} else{
				    $e = $type.": ".$id;
				}
			    } else{
				$e = $dbReference->att('type');
			    }
			    push (@evs, $e);
			}
			push @all_refs, @evs;
		    }

		    
		    # take the number of references if there's any #
		    my $reference = $feature->att('ref');
		    if($reference){
			my @numbers = split ' ', $reference;

			my $pubmed_link = "http://www.ncbi.nlm.nih.gov/pubmed/";
			for (@numbers) {
			    my $r;
			    foreach my $dbRef ($t->findnodes("//entry/evidence[\@key=$_]/source/dbReference")) {
				my $type = $dbRef->att('type');
				my $id = $dbRef->att('id');
				if($type eq "PubMed"){
				    $r = $pubmed_link.$id;
				} else{
				    $r = $type.": ".$id;
				}
				push(@refs, $r);
			    }
			}
			push @all_refs, @refs;
		    }

		} else {
		    # put by simlarity in the array of references (when the modification doesn't have a reference)
		    push @all_refs, $status;
		}
		$hash_node{all_refs} = \@all_refs;
		push(@node_data, \%hash_node);
	    }
	    push @feat_data, @node_data;

	}
	$modif_hash{FEATURE_DATA} = \@feat_data;
	push @{$$prot_hashref{MODIFICATION_LIST}}, \%modif_hash;
    }

}

close INFILE;
close OUT1;
