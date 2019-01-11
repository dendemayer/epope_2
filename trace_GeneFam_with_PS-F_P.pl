#!/usr/bin/perl -w
#
# File:   trace_GeneFam.pl    
# Author: jana@vodka.bioinf.uni-leipzig.de
# Date:   2016/05/02 09:05:47
#

use strict;
use Getopt::Std;

my %args;
getopts('hvi:', \%args);

if(exists($args{h})) { usage(); exit(0); }
if(exists($args{v})) { version(); exit(0); }

if(! exists($args{i})) { usage(); exit(1); }

# --- MY BEGIN ----------------------------------------------------------------

# create the tree data structure
# holding only the info about geneFam and geneLoss
# key: label of node
my ($k,%tree);
open(I, "<", $args{i}) or die $!;
while(<I>) {
    chomp();
    my @x = split /\s+/,$_;

 #   print "lab: ",$x[3], " par: ", $x[9], " gainF: ",$x[19], " lossF: ", $x[21], " node: ", $x[24], "\n";
    
    $k = $x[3]; # label of node
    $tree{$k}{par} = $x[9]; # label of parent node
    $tree{$k}{gF}  = $x[17]; # gainFam
    $tree{$k}{lF}  = $x[19]; # lossFam
    $tree{$k}{name}= $x[24]; # name of node
   
}
close(I);

## trace gain and loss of families for each node
my ($Nfam, $j);
foreach my $i (sort {$a<=>$b} keys %tree) {

    $Nfam = $tree{$i}{gF} - $tree{$i}{lF};

    $tree{$i}{Nfam} = traceIt($i, 0, %tree);
    print "label: ", $i,"\t",$tree{$i}{name},"\tN_families: ", $tree{$i}{Nfam},"\n";
}

# --- MY subroutines -----------------------------------------------------------

sub traceIt ($,$,%) {
    my ($p, $Nfam, %T) = @_;


    
    $Nfam += $T{$p}{gF} - $T{$p}{lF};

#    print "\tp: ",$p, "(", $T{$p}{name},") Nfam: ", $Nfam, "\n";

    if($T{$p}{par} >= 0) {
	$Nfam = traceIt($T{$p}{par}, $Nfam, %T);
    }
    
    return $Nfam;
    
}


# --- MY END -------------------------------------------------------------------

sub usage { 

    printf "\n\n"; 

    printf "trace_GeneFam.pl - Reads a summarized ePoPE text file and
    reports the number of gene families for each of the nodes\n";

    version();

    printf "CALL:\n";
    printf "\ttrace_GeneFam.pl [-h|-v] -i FILE\n\n";

    printf "OPTIONS:\n";
    printf "\t-help   Print a brief help message and exits.\n\n";

    # list of options + description
    printf "\t%-10s\tFILE the text file that comes from ePoPE.summariz.pl.\n","-i";
    
    printf "DESCRIPTION:\n";
    
    printf "trace_GeneFam.pl - Reads a summarized ePoPE text file and                                                                                                                                                                        
    reports the number of gene families for each of the nodes\n";
    
    printf "\n\n";

    return;
}

sub version {
    
    printf "\nVERSION: 1.0\n\n";
    return;
}

sub HELP_MESSAGE {
    usage();
    exit(0);
}

sub VERSION_MESSAGE {
    version();
    exit(0)
}

