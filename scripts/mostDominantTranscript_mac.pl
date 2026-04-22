#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 19.02.2017

###############################################################################
###############################################################################
### Determines most dominant transcript of a gene.                          ###
###############################################################################
###############################################################################

use strict;
use warnings;
use Getopt::Long;
my (
    # variable for parameters which are read in from commandline
    $help,
    $gene2transFile,
    $expressFile,
    $minEnrichment,
    $verbose,
    $expressOut,
   );

##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "help!"                 => \$help,           # print this help
    "expressFile=s"         => \$expressFile,    # e.g. freeze3_v2.kallisto.lib.trans.fpkm.wl.aliquot_id.tsv.gz
    "ensgFile=s"            => \$gene2transFile, # e.g. ensg_ensp_enst_ense_geneName_v75.tsv.gz
    "minEnrichment=f"       => \$minEnrichment,  # e.g. 2
    "out=s"                 => \$expressOut,     # filename of gzipped file to be created for expression values for MDTs. 
    "verbose!"              => \$verbose,        # print out additional information on calculation progress plus warning messages
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

##############################################################################
### SETTINGS #################################################################
##############################################################################
$| = 1;

##############################################################################
### SUBROUTINES ##############################################################
############################################################################## 

###############################################################################
sub printHelp {
###############################################################################
    # prints a help about the using and parameters of this scripts 
    # (execute if user types commandline parameter -h)
    # param:  no paramaters
    # return: no return value

    my (
	$usage,
	$sourceCode,
	@rows,
	$row,
	$option,
	$scriptInfo,
	$example,
	);

    $usage = "$0\n";


    print "\nUsage: " .  $usage . "\n";

    print "Valid options are:\n\n";
    open(MYSELF, "$0") or
	die "Cannot read source code file $0: $!\n";
    $sourceCode .= join "", <MYSELF>;
    close MYSELF;
    $sourceCode =~ s/.*\&GetOptions\(//s;
    $sourceCode =~ s/\).+//s;
    @rows = split /\n/, $sourceCode;
    foreach $row (@rows){
        $option = $row;
	$option =~ s/\s+\"//g;
	$option =~ s/\"\s.+\#/\t\#/g;
	$option =~ s/=./\t<value> [required]/;
	$option =~ s/:./\t<value> [optional]/;
	$option =~ s/!/\t<non value> [optional]/;

	$row =~ s/^.*//;
	print "\t";
	printf("%-1s%-30s%-30s\n", "-",$option,$row);

    } # end of foreach $row (@rows)
    print "\n";
    print "Options may be abreviated, e.g. -h for --help\n\n";

    $example  = "$0";
}

###############################################################################
sub readGene2transFile {
    my %gene2trans = (); 
    open(F1, "gzcat $gene2transFile|") or die "\nERROR: Failed to open $gene2transFile: $!\n";
    while(my $l = <F1>) {
	next if($l !~ /^ENS/);
	chomp($l);
	my @a = split(/\t/, $l);
	$gene2trans{$a[2]} = $a[0];
    }
    close(F1);
    return \%gene2trans;
}
###############################################################################
sub readExpressionFile {
    my $gene2trans = &readGene2transFile();

    my %express = ();
    open(F2, "gzcat $expressFile|") or die "\nERROR: Failed to open $expressFile: $!\n";
    my $n = 0;
    my @samples = ();
    while(my $l = <F2>) {
	chomp($l);
	if($n == 0) {
	    @samples = split(/\t/, $l);
	    $n++;
	    next;
	}

	my @a = split(/\t/, $l);
	my $enst = $a[0];
	$enst =~ s/\..*//;

	if(!exists $gene2trans->{$enst}) {
	    print STDERR "WARNING: $enst does not exist in $gene2transFile. Skipping transcript.\n";
	    next;
	}
	my $ensg = $gene2trans->{$enst};

	for(my $i = 1; $i < @a; $i++) {
	    my $expression = $a[$i];
	    $expression = 0 if($expression eq "NA" or $expression < 0.001);
	    $express{$ensg}->{$samples[$i]}->{$enst} = $expression;
	}
    }
    close(F2);

    return (\%express);
}
##############################################################################
### END OF SUBROUTINES########################################################
############################################################################## 


############
### MAIN ###
############

if(!defined $gene2transFile) {
    print STDERR "\n\tPlease provide an Ensembl ID mapping data file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $expressFile) {
    print STDERR "\n\tPlease provide an expression table file. Try $0 -help to get more information\n\n";
    exit;
}
if(!defined $minEnrichment) {
    print STDERR "\n\tPlease provide a minimum fold difference value. Try $0 -help to get more information\n\n";
    exit;
}
my ($express) = &readExpressionFile();

# print header
foreach my $ensg (sort keys %$express) {
    print "#Gene";
    foreach my $sample (sort keys %{$express->{$ensg}}) {
	print "\t$sample";
    }
    print "\n";
    last;
}
if(defined $expressOut) {
    open(O, "| gzip > $expressOut");
    print O "#Sample\tENSG\tNoOfENST\tENST1\tENST2\tTPM1\tTPM2\tEnrichment\n";
}

foreach my $ensg (sort keys %$express) {
    print $ensg;
    foreach my $sample (sort keys %{$express->{$ensg}}) {
	my $mostDominantTrans = "-";
	my $preTransExpress = -1;
	my $preTrans = "-";
	my $trans = "-";
	my $transExpress = -1;
	my $enrichment = -1;
	my $noTrans = keys %{$express->{$ensg}->{$sample}};

	foreach my $enst (sort {$express->{$ensg}->{$sample}->{$b} <=> $express->{$ensg}->{$sample}->{$a}}
			  keys %{$express->{$ensg}->{$sample}}) {

	    $trans = $enst;
	    $transExpress = $express->{$ensg}->{$sample}->{$trans};

	    if($preTransExpress == -1) {
		$preTransExpress = $transExpress;
		$preTrans = $trans;
		$enrichment = $transExpress;
		last if($transExpress < 0.001); # don't call most dominant transcript for genes that are not expressed.
		next if($noTrans > 1); # if there is no 2nd transcript, finish most dominant transcript determination immediatly.
	    }
	    $enrichment /= $transExpress if($transExpress > 0);
	    
	    # if transcript is the only transcript, or 2nd MDI is no too little expressed than by default it is the most dominant one.
	    if($enrichment >= $minEnrichment or $noTrans == 1 or $transExpress < 0.001) {
		$mostDominantTrans = $preTrans;
	    }
	    last;
	}
	print O "$sample\t$ensg\t$noTrans\t$preTrans\t$trans\t$preTransExpress\t$transExpress\t$enrichment\n" if(defined $expressOut and $preTransExpress > 0);
	print "\t$mostDominantTrans";
    }
    print "\n";
}
close(O) if(defined $expressOut);
