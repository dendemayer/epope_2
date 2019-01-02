#!/usr/bin/perl 
use strict;
use Getopt::Long;
use Pod::Usage;

my $man = 0;
my $help = 0;
my $both = 0;
my $version = 0;
my ($inDir, $tree,$outfile,$treefile) = ('','GLtreeTemplate','','');



### ------------------------------------------------------------------------- ###
# manpage, help, usage..

GetOptions(
			"help|?"=> \$help,
			"man"   => \$man,
			"d:s"   => \$inDir,
			"o:s"   => \$outfile,
			"v|version" 	=> \$version
          ) or pod2usage(1);

pod2usage(1) if $help;
if($version == 1) { version(); exit(0); }
pod2usage(-verbose => 2) if $man;

#if($tree eq '') {
#  pod2usage(1);
#}

if(! (-d $inDir)) 
{
  printf STDERR "\nERROR: Cannot acces input directory '%s'\n", $inDir;
  pod2usage(1);
}

		print "\nsearching for output files from ePoPE in Parsimony mode in: $inDir \n(Files with an '_PS.out' -extension)\n\n";
	
		# read directory of data
			
		opendir(D, $inDir) or die $!;
		my @entries = grep {/\_PS.out$/} readdir(D);
		closedir(D);
		my $PSlaenge=@entries;
	if ($PSlaenge==0)
	{
		printf("No '_PS.out' Files found\n\n")
	}
		
	if ($PSlaenge!=0)
	{
		# fill tree hash with data
		my %PStree;
		my %Famtree;
		my ($i,$nodes_n) = (0,0);
		foreach my $e (@entries) 
		{
			my $file = $inDir."".$e;
			printf "%s\n",$file; 
			open(E, $file) or die $!;
			
			while(<E>) 
			{
				chomp();
				my $line = $_;
			#   print $line;
				my @x = split /\s+/,$line;
			#    next;
				$i = $x[3]; #/* use i as identifier and lab of tree */ -->lab wise sort
				# write first instance in tree
				
				if(! exists($PStree{$i})) 
				{
					$PStree{$i}{name} = $x[scalar @x -1]; chomp($PStree{$i}{name}); #letzter eintrag ist wohl species name, sollte demnach passen...
					$PStree{$i}{pos} = $x[1];
					$PStree{$i}{lab} = $x[3];
					$PStree{$i}{p} = $x[9];
					$PStree{$i}{c} = $x[5];
					$PStree{$i}{m} = $x[11];
					$PStree{$i}{gain} = $x[13];
					$PStree{$i}{loss} = $x[15];
					$PStree{$i}{gainF} = $x[17];
					$PStree{$i}{lossF} = $x[19];
					$nodes_n++;
				}
				else 
				{ # add values to existing tree
					$PStree{$i}{m} += $x[11];
					$PStree{$i}{gain} += $x[13];
					$PStree{$i}{loss} += $x[15];
					$PStree{$i}{gainF} += $x[17];
					$PStree{$i}{lossF} += $x[19];
				}
				
				if(! exists($Famtree{$i}))
				{
					$Famtree{$i}{p} = $x[9]; # label of parent node
					$Famtree{$i}{gainF}  = $x[17]; # gainFam
					$Famtree{$i}{lossF}  = $x[19]; # lossFam
					$Famtree{$i}{name}= $x[20]; # name of node
					$Famtree{$i}{Nfam}=0;
				}             
				else          
				{             
					$Famtree{$i}{gainF}  += $x[17]; # gainFam
					$Famtree{$i}{lossF}  += $x[19]; # lossFam
				}
			}
			close(E);
		}
	
	 	my ($Nfam);
		sub traceIt;
		
		foreach my $i (sort {$a<=>$b} keys %Famtree)       #The spaceship <=> operator is used for comparing numbers, here sort for lab (reihenfolge eigtl egal, es darf nur kein schlüssel gewählt werden den es nicht gibt...)
		{
		
			$Famtree{$i}{Nfam} = traceIt($i, 0, %Famtree); # traceIt(label,Nfam,%tree) $i ist hier lab number da beim einlesen diese spalte für die cshlüssel verwendet wird
			 #print "label: ", $i,"\t name:",$Famtree{$i}{name}, " \t  lossF: ",$Famtree{$i}{lossF}," gainF: ", $Famtree{$i}{gainF},  "  parent: ", $Famtree{$i}{p} , "\t Nfam:  ", $Famtree{$i}{Nfam} ,  "\n"  ;
		
		}	
		
		my $outfilePS;
		
		if( $outfile =~ /\./)
		{
			$outfilePS = substr($outfile,0, rindex( $outfile, '.') );
			$outfilePS = $outfilePS."_PS.sum";
			
		}
		else
		{	$outfilePS = $outfile."_PS.sum";}
		
		open(O, ">", $outfilePS) or die $!;
		#for($i=0; $i<$nodes_n; $i++) 
		#	{
				#um kompatibel mit ePoPe -c zu sein, muss in pos reihefolge gelistet werden
				#sortierung nach $PStree{$i}{pos} = $x[1];
			
		for my $PSkey ( sort {$PStree{$a}->{pos} <=> $PStree{$b}->{pos}} keys %PStree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
		{			
			#my $value = $PStree{$PSkey};
			printf O "pos: %d lab: %d kids: %s pInP: %d par: %d ",$PStree{$PSkey}{pos},$PStree{$PSkey}{lab},$PStree{$PSkey}{c},$PSkey,$PStree{$PSkey}{p};
			printf O "genes: %d gain: %d loss: %d ",$PStree{$PSkey}{m},$PStree{$PSkey}{gain},$PStree{$PSkey}{loss};
			#print "the NFam of ",$PStree{$PSkey}{name} , " is ", $Famtree{$PSkey}{Nfam} , " (also : " ,$Famtree{$PSkey}{name}, ") \$nodes_n=", $nodes_n, " \n";
			printf O "gainFam: %d lossFam: %d N_Fam: %d %s\n",$PStree{$PSkey}{gainF},$PStree{$PSkey}{lossF},$Famtree{$PSkey}{Nfam},$Famtree{$PSkey}{name};
		}
					# printf O "count: %d\n", $tree{$i}{count}; # test count
					#printf("famtest: species: %s NFAM: %d\n", $PStree{$j}{name} , $PStree{$j}{Nfam} )
		
		
		printf("\nThe Summary File is named: %s\n\n", $outfilePS);
		close(O);
	}	
	
	######################## now do the same for the files with the _PF extension:
		print "\nsearching for output files from ePoPE in Partition Function mode in: $inDir \n(Files with an '_PF.out' -extension)\n";
				
		opendir(D, $inDir) or die $!;
		my @PFentries = grep {/\_PF.out$/} readdir(D);
		closedir(D);
		my $laenge=@PFentries;
	if ($laenge==0)
	{
		printf("No '_PF.out' Files found\n\n")
	}	
		
		
	if ($laenge!=0)
	{
		#printf("länge von entries: %d\n", $laenge);
		
		# fill tree hash with data
		my %PFtree;
		my %FamtreePF;
		my ($ii,$Nnodes_n) = (0,0);
		foreach my $e (@PFentries) 
		{
			my $file = $inDir."".$e;
			printf "%s\n",$file; 
			open(E, $file) or die $!;
			
			while(<E>) 
			{
				chomp();
				my $line = $_;
			#   print $line;
				my @x = split /\s+/,$line;
			#    next;
			
			#pos:  16 lab: 329 kids: -1 pInP:  16 par: 327 mirs:   0 gain:   0 loss   0       seenBW:   1 gainFam:   0        lossFam:   0    egr
			
				$ii = $x[3]; 
				my $count = 0; 
				# write first instance in tree
				if(! exists($PFtree{$ii})) 
				{
					$PFtree{$ii}{name} = $x[scalar @x -1]; chomp($PFtree{$ii}{name}); #letzter eintrag ist wohl species name, sollte demnach passen...
					$PFtree{$ii}{pos} = $x[1];
					$PFtree{$ii}{lab} = $x[3];
					$PFtree{$ii}{p} = $x[9];
					$PFtree{$ii}{c} = $x[5];
					$PFtree{$ii}{m} = $x[11];
					$PFtree{$ii}{gain} = $x[13];
					$PFtree{$ii}{loss} = $x[15];
					$PFtree{$ii}{gainF} = $x[17];
					$PFtree{$ii}{lossF} = $x[19];
					if($PFtree{$ii}{m} == 0) #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit
					{
						$PFtree{$ii}{count} = 0;
						$PFtree{$ii}{pf_P} = 0;
						
					}
					else
					{
						$PFtree{$ii}{count} = 1;		
						$PFtree{$ii}{pf_P} = $x[21];
						#print "m von tree.$i ist $x[11]   zu addierende PF_P: $x[23] aus datei $e\n";
					}
					$Nnodes_n++;
				} 	
				else 
				{ # add values to existing tree
					$PFtree{$ii}{m} += $x[11];
					$PFtree{$ii}{gain} += $x[13];
					$PFtree{$ii}{loss} += $x[15];
					$PFtree{$ii}{gainF} += $x[17];
					$PFtree{$ii}{lossF} += $x[19];
					if($x[11] != 0 && $x[5] ne '-1') #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit//bei kids($x[5]) die mit -1 auslassen ist aber char!, also stringcompare
					{
						#printf( "m von tree.[%d] ist %d   zu addierende PF_P: %f aus datei %s mit namen %s\n",$ii,$x[11],$x[21],$e,$x[scalar @x -1]);
						$PFtree{$ii}{pf_P} += $x[21];  # um relative wahrscheinlichkeit auszurechnen, pf_P mitzählen wenn m über 0 (wenn 0, dann ist knoten außerhalb des alignments)
						$PFtree{$ii}{count} += 1;
					}
				}

				if(! exists($FamtreePF{$ii}))
				{
					$FamtreePF{$ii}{p} = $x[9]; # label of parent node
					$FamtreePF{$ii}{gainF}  = $x[17]; # gainFam
					$FamtreePF{$ii}{lossF}  = $x[19]; # lossFam
					$FamtreePF{$ii}{name}= $x[22]; # name of node
					$FamtreePF{$ii}{Nfam}=0;
					
				}
				else
				{
					$FamtreePF{$ii}{gainF}  += $x[17]; # gainFam
					$FamtreePF{$ii}{lossF}  += $x[19]; # lossFam
				}
				
			}
			close(E);
		}
			
		foreach my $i (sort {$a<=>$b} keys %FamtreePF)       #The spaceship <=> operator is used for comparing 
		{
			
				$FamtreePF{$i}{Nfam} = traceIt($i, 0, %FamtreePF);
				#print "label: ", $i,"\t name:",$FamtreePF{$i}{name}, " \t  lossF: ",$FamtreePF{$i}{lossF}," gainF: ", $FamtreePF{$i}{gainF},  "  parent: ", $FamtreePF{$i}{p} ,"\t Nfam: ", $FamtreePF{$i}{Nfam} ,  "\n"  ;
			
		}	
			
		my $outfilePF ;
		if( $outfile =~ /\./)
		{
			$outfilePF = substr($outfile,0, rindex( $outfile, '.') );
			$outfilePF = $outfilePF."_PF.sum";
			
		}
		else
		{$outfilePF = $outfile."_PF.sum";}
		
		open(O, ">", $outfilePF) or die $!;
		
		for my $PFkey ( sort {$PFtree{$a}->{pos} <=> $PFtree{$b}->{pos}} keys %PFtree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
		{			
			#my $value = $PStree{$PSkey};
			printf O "pos: %d lab: %d kids: %s pInP: %d par: %d ",$PFtree{$PFkey}{pos},$PFtree{$PFkey}{lab},$PFtree{$PFkey}{c},$PFkey,$PFtree{$PFkey}{p};
			printf O "genes: %d gain: %d loss: %d ",$PFtree{$PFkey}{m},$PFtree{$PFkey}{gain},$PFtree{$PFkey}{loss};
			#print "the NFam of ",$PFtree{$PFkey}{name} , " is ", $Famtree{$PFkey}{Nfam} , " (also : " ,$Famtree{$PFkey}{name}, ") \$nodes_n=", $nodes_n, " \n";
			if ($PFtree{$PFkey}{count} ==0){$PFtree{$PFkey}{count} =1;}
			printf O "gainFam: %d lossFam: %d N_Fam: %d PF_P: %f %s\n",$PFtree{$PFkey}{gainF},$PFtree{$PFkey}{lossF},$FamtreePF{$PFkey}{Nfam},$PFtree{$PFkey}{pf_P}/$PFtree{$PFkey}{count},$FamtreePF{$PFkey}{name};
		}
		
				
		printf("\nThe Summary File is named: %s\n\n", $outfilePF);
		close(O);
		
	}
	
	######################## now do the same for the files with the _PSb extension:
		print "\nsearching for output files from ePoPE in both mode in: $inDir \n(Files with an '_PSb.out' -extension)\n\n";
				
		opendir(D, $inDir) or die $!;
		my @PSbentries = grep {/\_PSb.out$/} readdir(D);
		closedir(D);
		$laenge=@PSbentries;
	if ($laenge==0)
	{
		printf("No '_PSb.out' Files found\n\n")
	}	
		
		
	if ($laenge!=0)
	{
		#printf("länge von entries: %d\n", $laenge);
		
		# fill tree hash with data
		my %PSbtree;
		my %FamtreePSb;
		my ($ii,$Nnodes_n) = (0,0);
		foreach my $e (@PSbentries) 
		{
			my $file = $inDir."".$e;
			printf "%s\n",$file; 
			open(E, $file) or die $!;
			
			while(<E>) 
			{
				chomp();
				my $line = $_;
			#   print $line;
				my @x = split /\s+/,$line;
			#    next;
			
			#pos:  16 lab: 329 kids: -1 pInP:  16 par: 327 mirs:   0 gain:   0 loss   0       seenBW:   1 gainFam:   0        lossFam:   0    egr
			
				$ii = $x[3]; 
				my $count = 0; 
				# write first instance in tree
				if(! exists($PSbtree{$ii})) 
				{
					$PSbtree{$ii}{name} = $x[scalar @x -1]; chomp($PSbtree{$ii}{name}); #letzter eintrag ist wohl species name, sollte demnach passen...
					$PSbtree{$ii}{pos} = $x[1];
					$PSbtree{$ii}{lab} = $x[3];
					$PSbtree{$ii}{p} = $x[9];
					$PSbtree{$ii}{c} = $x[5];
					$PSbtree{$ii}{m} = $x[11];
					$PSbtree{$ii}{gain} = $x[13];
					$PSbtree{$ii}{loss} = $x[15];
					$PSbtree{$ii}{gainF} = $x[17];
					$PSbtree{$ii}{lossF} = $x[19];
					if($PSbtree{$ii}{m} == 0) #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit
					{
						$PSbtree{$ii}{count} = 0;
						$PSbtree{$ii}{pf_P} = 0;
						
					}
					else
					{
						$PSbtree{$ii}{count} = 1;		
						$PSbtree{$ii}{pf_P} = $x[21];
						#print "m von tree.$i ist $x[11]   zu addierende PF_P: $x[23] aus datei $e\n";
					}
					$Nnodes_n++;
				} 	
				else 
				{ # add values to existing tree
					$PSbtree{$ii}{m} += $x[11];
					$PSbtree{$ii}{gain} += $x[13];
					$PSbtree{$ii}{loss} += $x[15];
					$PSbtree{$ii}{gainF} += $x[17];
					$PSbtree{$ii}{lossF} += $x[19];
					if($x[11] != 0 && $x[5] ne '-1') #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit//bei kids($x[5]) die mit -1 auslassen ist aber char!, also stringcompare
					{
						#printf( "m von tree.[%d] ist %d   zu addierende PF_P: %f aus datei %s mit namen %s\n",$ii,$x[11],$x[21],$e,$x[scalar @x -1]);
						$PSbtree{$ii}{pf_P} += $x[21];  # um relative wahrscheinlichkeit auszurechnen, pf_P mitzählen wenn m über 0 (wenn 0, dann ist knoten außerhalb des alignments)
						$PSbtree{$ii}{count} += 1;
					}
				}

				if(! exists($FamtreePSb{$ii}))
				{
					$FamtreePSb{$ii}{p} = $x[9]; # label of parent node
					$FamtreePSb{$ii}{gainF}  = $x[17]; # gainFam
					$FamtreePSb{$ii}{lossF}  = $x[19]; # lossFam
					$FamtreePSb{$ii}{name}= $x[22]; # name of node
					$FamtreePSb{$ii}{Nfam}=0;
					
				}
				else
				{
					$FamtreePSb{$ii}{gainF}  += $x[17]; # gainFam
					$FamtreePSb{$ii}{lossF}  += $x[19]; # lossFam
				}
				
			}
			close(E);
		}
			
		foreach my $i (sort {$a<=>$b} keys %FamtreePSb)       #The spaceship <=> operator is used for comparing 
		{
			
				$FamtreePSb{$i}{Nfam} = traceIt($i, 0, %FamtreePSb);
				#print "label: ", $i,"\t name:",$FamtreePSb{$i}{name}, " \t  lossF: ",$FamtreePSb{$i}{lossF}," gainF: ", $FamtreePSb{$i}{gainF},  "  parent: ", $FamtreePSb{$i}{p} ,"\t Nfam: ", $FamtreePSb{$i}{Nfam} ,  "\n"  ;
			
		}	
		
		my $outfilePSb;
		if( $outfile =~ /\./)
		{
			$outfilePSb = substr($outfile,0, rindex( $outfile, '.') );
			$outfilePSb = $outfilePSb."_PSb.sum";
			
		}
		else
		{ $outfilePSb = $outfile."_PSb.sum";}
			
		open(O, ">", $outfilePSb) or die $!;
		
		for my $PSbkey ( sort {$PSbtree{$a}->{pos} <=> $PSbtree{$b}->{pos}} keys %PSbtree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
		{			
			#my $value = $PStree{$PSkey};
			printf O "pos: %d lab: %d kids: %s pInP: %d par: %d ",$PSbtree{$PSbkey}{pos},$PSbtree{$PSbkey}{lab},$PSbtree{$PSbkey}{c},$PSbkey,$PSbtree{$PSbkey}{p};
			printf O "genes: %d gain: %d loss: %d ",$PSbtree{$PSbkey}{m},$PSbtree{$PSbkey}{gain},$PSbtree{$PSbkey}{loss};
			#print "the NFam of ",$PSbtree{$PSbkey}{name} , " is ", $Famtree{$PSbkey}{Nfam} , " (also : " ,$Famtree{$PSbkey}{name}, ") \$nodes_n=", $nodes_n, " \n";
			if ($PSbtree{$PSbkey}{count} ==0){$PSbtree{$PSbkey}{count} =1;}
			printf O "gainFam: %d lossFam: %d N_Fam: %d PF_P: %f %s\n",$PSbtree{$PSbkey}{gainF},$PSbtree{$PSbkey}{lossF},$FamtreePSb{$PSbkey}{Nfam},$PSbtree{$PSbkey}{pf_P}/$PSbtree{$PSbkey}{count},$FamtreePSb{$PSbkey}{name};
		}
		
				
		printf("\nThe Summary File is named: %s\n\n", $outfilePSb);
		close(O);
		
	}
		
	######################## now do the same for the files with the _PFb extension:
		print "\nsearching for output files from ePoPE in both mode in: $inDir \n(Files with an '_PFb.out' -extension)\n\n";
				
		opendir(D, $inDir) or die $!;
		my @PFbentries = grep {/\_PFb.out$/} readdir(D);
		closedir(D);
		$laenge=@PFbentries;
	if ($laenge==0)
	{
		printf("No '_PFb.out' Files found\n\n")
	}	
			
	if ($laenge!=0)
	{
		#printf("länge von entries: %d\n", $laenge);
		
		# fill tree hash with data
		my %PFbtree;
		my %FamtreePFb;
		my ($ii,$Nnodes_n) = (0,0);
		foreach my $e (@PFbentries) 
		{
			my $file = $inDir."".$e;
			printf "%s\n",$file; 
			open(E, $file) or die $!;
			
			while(<E>) 
			{
				chomp();
				my $line = $_;
			#   print $line;
				my @x = split /\s+/,$line;
			#    next;
			
			#pos:  16 lab: 329 kids: -1 pInP:  16 par: 327 mirs:   0 gain:   0 loss   0       seenBW:   1 gainFam:   0        lossFam:   0    egr
			
				$ii = $x[3]; 
				my $count = 0; 
				# write first instance in tree
				if(! exists($PFbtree{$ii})) 
				{
					$PFbtree{$ii}{name} = $x[scalar @x -1]; chomp($PFbtree{$ii}{name}); #letzter eintrag ist wohl species name, sollte demnach passen...
					$PFbtree{$ii}{pos} = $x[1];
					$PFbtree{$ii}{lab} = $x[3];
					$PFbtree{$ii}{p} = $x[9];
					$PFbtree{$ii}{c} = $x[5];
					$PFbtree{$ii}{m} = $x[11];
					$PFbtree{$ii}{gain} = $x[13];
					$PFbtree{$ii}{loss} = $x[15];
					$PFbtree{$ii}{gainF} = $x[17];
					$PFbtree{$ii}{lossF} = $x[19];
					if($PFbtree{$ii}{m} == 0) #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit
					{
						$PFbtree{$ii}{count} = 0;
						$PFbtree{$ii}{pf_P} = 0;
						
					}
					else
					{
						$PFbtree{$ii}{count} = 1;		
						$PFbtree{$ii}{pf_P} = $x[21];
						#print "m von tree.$i ist $x[11]   zu addierende PF_P: $x[23] aus datei $e\n";
					}
					$Nnodes_n++;
				} 	
				else 
				{ # add values to existing tree
					$PFbtree{$ii}{m} += $x[11];
					$PFbtree{$ii}{gain} += $x[13];
					$PFbtree{$ii}{loss} += $x[15];
					$PFbtree{$ii}{gainF} += $x[17];
					$PFbtree{$ii}{lossF} += $x[19];
					if($x[11] != 0 && $x[5] ne '-1') #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit//bei kids($x[5]) die mit -1 auslassen ist aber char!, also stringcompare
					{
						#printf( "m von tree.[%d] ist %d   zu addierende PF_P: %f aus datei %s mit namen %s\n",$ii,$x[11],$x[21],$e,$x[scalar @x -1]);
						$PFbtree{$ii}{pf_P} += $x[21];  # um relative wahrscheinlichkeit auszurechnen, pf_P mitzählen wenn m über 0 (wenn 0, dann ist knoten außerhalb des alignments)
						$PFbtree{$ii}{count} += 1;
					}
				}

				if(! exists($FamtreePFb{$ii}))
				{
					$FamtreePFb{$ii}{p} = $x[9]; # label of parent node
					$FamtreePFb{$ii}{gainF}  = $x[17]; # gainFam
					$FamtreePFb{$ii}{lossF}  = $x[19]; # lossFam
					$FamtreePFb{$ii}{name}= $x[22]; # name of node
					$FamtreePFb{$ii}{Nfam}=0;
					
				}
				else
				{
					$FamtreePFb{$ii}{gainF}  += $x[17]; # gainFam
					$FamtreePFb{$ii}{lossF}  += $x[19]; # lossFam
				}
				
			}
			close(E);
		}
			
		foreach my $i (sort {$a<=>$b} keys %FamtreePFb)       #The spaceship <=> operator is used for comparing 
		{
			
				$FamtreePFb{$i}{Nfam} = traceIt($i, 0, %FamtreePFb);
				#print "label: ", $i,"\t name:",$FamtreePFb{$i}{name}, " \t  lossF: ",$FamtreePFb{$i}{lossF}," gainF: ", $FamtreePFb{$i}{gainF},  "  parent: ", $FamtreePFb{$i}{p} ,"\t Nfam: ", $FamtreePFb{$i}{Nfam} ,  "\n"  ;
			
		}	
	
		my $outfilePFb;
		if( $outfile =~ /\./)
		{
			$outfilePFb = substr($outfile,0, rindex( $outfile, '.') );
			$outfilePFb = $outfilePFb."_PFb.sum";
			
		}
		else
		{ $outfilePFb = $outfile."_PFb.sum";}
	
	
		open(O, ">", $outfilePFb) or die $!;
		
		for my $PFbkey ( sort {$PFbtree{$a}->{pos} <=> $PFbtree{$b}->{pos}} keys %PFbtree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
		{			
			#my $value = $PStree{$PSkey};
			printf O "pos: %d lab: %d kids: %s pInP: %d par: %d ",$PFbtree{$PFbkey}{pos},$PFbtree{$PFbkey}{lab},$PFbtree{$PFbkey}{c},$PFbkey,$PFbtree{$PFbkey}{p};
			printf O "genes: %d gain: %d loss: %d ",$PFbtree{$PFbkey}{m},$PFbtree{$PFbkey}{gain},$PFbtree{$PFbkey}{loss};
			#print "the NFam of ",$PFbtree{$PFbkey}{name} , " is ", $Famtree{$PFbkey}{Nfam} , " (also : " ,$Famtree{$PFbkey}{name}, ") \$nodes_n=", $nodes_n, " \n";
			if ($PFbtree{$PFbkey}{count} ==0){$PFbtree{$PFbkey}{count} =1;}
			printf O "gainFam: %d lossFam: %d N_Fam: %d PF_P: %f %s\n",$PFbtree{$PFbkey}{gainF},$PFbtree{$PFbkey}{lossF},$FamtreePFb{$PFbkey}{Nfam},$PFbtree{$PFbkey}{pf_P}/$PFbtree{$PFbkey}{count},$FamtreePFb{$PFbkey}{name};
		}
						
		printf("\nThe Summary File is named: %s\n\n", $outfilePFb);
		close(O);
		
	}
	

# --- MY subroutines -----------------------------------------------------------
		
sub traceIt($,$,%) 
		{
			#kommt hier auch alles an?
			my ($i, $Nfam, %T) = @_;
			
			$Nfam += $T{$i}{gainF} - $T{$i}{lossF};
			#print "\ti: ",$i, "(", $T{$i}{name},") Nfam: ", $Nfam;
			#print " genes: ", "parent: " ,$T{$i}{p}; #genes is not initialized
					   
			if($T{$i}{p} >= 0) 
			{
				#print "traceIt(" ,  $T{$i}{p},",",$Nfam,",", $T{$i}{name},")\n";
			    $Nfam = traceIt($T{$i}{p}, $Nfam, %T);
				#print "\ti: ",$i, "(", $T{$i}{name},") Nfam: ", $Nfam;
				#print " m: ",$T{$i}{m}, "p: " ,$T{$i}{p};
				
			}
			return $Nfam;
		}	
		
	

sub version 
{
	printf "\nVERSION: 2.0\n\n";
    return;
}	
	
	

### ------------------------------------------------------------------------- ###

__END__

=head1 NAME

ePoPE.summarize.pl

summarizing the data output files of ePoPE.

=head1 SYNOPSIS

ePoPE.summarize.pl [-d DIR/ -o FILE ]
 


=head1 OPTIONS

=over 8

=item B<-d DIR/>

Directory that contains the outputfiles from a ePoPE run

=item B<-o FILE>

the output file for the final summarized tree data.

=item B<--help | -h>

Print a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<ePoPE.summarize.pl> summarizes the data output files of
ePoPE, returns the summarized tree in space seperated table
format which can easily be read by the ePoPE -c option tool that
creates a PostScript tree, or with other tools that analyse the data
(e.g. R).

=cut
