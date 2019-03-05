#!/usr/bin/perl 
use strict;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;


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


######################## now do the same for the files with the _PF extension:
my @ePoPE_options=qw(PF.out PFb.out PSb.out); #PF -> PF -P, PFb default, PSb -> PS
my $ref_opt=\@ePoPE_options;
#print Dumper(@ePoPE_options);
#print Dumper($ref_opt);
#print Dumper(@{$ref_opt});

#~ Summarize("PF.out");
#~ Summarize("PFb.out");
#~ Summarize("PSb.out");
	
Summarize(@ { $ref_opt});	

my %Ges_Tree; #alles zwischenspeichern für Baum mit allen optionen
	
sub Summarize
{
	foreach my $Suffix(@_)
	{	
		my $sum_PF;
			
		print "the suffix: ", @_ ,"\n";	
		print "the suffix: ", Dumper(@_), "\n";	
		print "\nsearching for output files from ePoPE in Partition Function mode in: $inDir \n(Files with an $Suffix -extension)\n";
				
		opendir(D, $inDir) or die $!;
		my @PFentries = grep {  /$Suffix$/g } readdir(D);
		print "entries: \n ";
		print  Dumper(@PFentries),"\n";
		
		print "#################################\n";
		closedir(D);
		my $laenge=@PFentries;
		if ($laenge==0)
		{
			printf("No '$Suffix' Files found\n\n")
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
					my @x = split /\s+/,$line;		#in zeile @x sind spalten $x[] 
				#    next;
				
				#pos:  16 lab: 329 kids: -1 pInP:  16 par: 327 mirs:   0 gain:   0 loss   0       seenBW:   1 gainFam:   0        lossFam:   0    egr
				
					$ii = $x[3]; #lab
					my $count = 0; 
					# write first instance in tree
					if(! exists($PFtree{$ii})) 
					{
						$PFtree{$ii}{name} = $x[scalar @x -1]; chomp($PFtree{$ii}{name}); #letzter eintrag ist wohl species name, sollte demnach passen...
						$Ges_Tree{$Suffix}{$ii}{name} = $x[scalar @x -1]; chomp($PFtree{$ii}{name});
						
						$PFtree{$ii}{pos} = $x[1];
						$Ges_Tree{$Suffix}{$ii}{pos} = $x[1];
						
						$PFtree{$ii}{lab} = $x[3];
						$Ges_Tree{$Suffix}{$ii}{lab} = $x[3];
						
						$PFtree{$ii}{p} = $x[9];
						
						
						$PFtree{$ii}{c} = $x[5];	# child array
						$Ges_Tree{$Suffix}{$ii}{c} = $x[5];
						
						#print $PFtree{$ii}{c};
						
						$PFtree{$ii}{m} = $x[11];
						$Ges_Tree{$Suffix}{$ii}{m} = $x[11];
						
						$PFtree{$ii}{gain} = $x[13];
						$Ges_Tree{$Suffix}{$ii}{gain} = $x[13];
						
						$PFtree{$ii}{loss} = $x[15];
						$Ges_Tree{$Suffix}{$ii}{loss} = $x[15];
						
						$PFtree{$ii}{gainF} = $x[17];
						$PFtree{$ii}{lossF} = $x[19];
						#print Dumper($PFtree{$ii}{name},$ii);
						
						if($PFtree{$ii}{m} == 0) #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit
						{
							$PFtree{$ii}{count} = 0;
							$Ges_Tree{$Suffix}{$ii}{count} = 0;
							
							$PFtree{$ii}{pf_P} = 0;
							
						}
						else
						{
							$PFtree{$ii}{count} = 1;
							$Ges_Tree{$Suffix}{$ii}{count} = 1;
									
							$PFtree{$ii}{pf_P} = $x[21];
							#print "m von tree.$i ist $x[11]   zu addierende PF_P: $x[23] aus datei $e\n";
						}
						$Nnodes_n++;
					} 	
					else 
					{ # add values to existing tree
						$PFtree{$ii}{m} += $x[11];
						$Ges_Tree{$Suffix}{$ii}{m} += $x[11];
					
						$PFtree{$ii}{gain} += $x[13];
						$Ges_Tree{$Suffix}{$ii}{gain} += $x[13];
					
						$PFtree{$ii}{loss} += $x[15];
						$Ges_Tree{$Suffix}{$ii}{loss} += $x[15];
						
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
			
			#print ("\n\n Dumper(\%Pftree)\n\n" );
			#print Dumper(%PFtree);
				
			foreach my $i (sort {$a<=>$b} keys %FamtreePF)       #The spaceship <=> operator is used for comparing ob numbers
			{
				
					$FamtreePF{$i}{Nfam} = traceIt($i, 0, %FamtreePF);
					#print "label: ", $i,"\t name:",$FamtreePF{$i}{name}, " \t  lossF: ",$FamtreePF{$i}{lossF}," gainF: ", $FamtreePF{$i}{gainF},  "  parent: ", $FamtreePF{$i}{p} ,"\t Nfam: ", $FamtreePF{$i}{Nfam} ,  "\n"  ;
				
			}	
				
			my $outfilePF ;
			if( $outfile =~ /\./)
			{
				$outfilePF = substr($outfile,0, rindex( $outfile, '.') );			#rindex STR,SUBSTR
				$outfilePF = $outfilePF.$Suffix.".sum";				#Works just like index except that it returns the position of the last occurrence of SUBSTR in STR. If POSITION is specified,returns the last occurrence beginning at or before that position.
	
				
			}
			else
			{$outfilePF = $outfile.$Suffix.".sum";}
			
			print "outfilePF: ", Dumper($outfilePF), "\n\n";
			
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
		
			$sum_PF= $outfilePF;
			
			
			
			###################################################	
			#new .dot output: with graphviz a svg file is created using this .dot file
			
			
			
			$outfilePF = $Suffix."_dot_tree.dot";
			open(O, ">", $outfilePF) or die $!;
			
			printf O ("graph graphname{\n\nnode [ fontsize=25,shape=record,nodesep=0,];edge [fontsize=25, decorate=true, penwidth=3];graph [ranksep=0];rankdir=LR;ratio=\"auto\";\n\n");
				
			for my $PFkey ( sort {$PFtree{$a}->{pos} <=> $PFtree{$b}->{pos}} keys %PFtree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
			{			
				if ($PFtree{$PFkey}{count} ==0){$PFtree{$PFkey}{count} =1;}
				if ($PFtree{$PFkey}{c} ne '-1')
				{				
					my @kidsarray = split /\,/,$PFtree{$PFkey}{c};
					#print join(", ",@kidsarray, "\n");					#für jedes element in @kidsarray wird gesamte baum durchsucht wenn $kidsarray[$p] eq lab print $PFtree{$PFkey}{c} -> $PFtree{$ii}{name} 
					print O $PFtree{$PFkey}{name}, "[ label=\" {{",$PFtree{$PFkey}{name},"| k = ",$PFtree{$PFkey}{m} ,"}}\"];\n";  #  T [shape=record, label="{T|000}"];
					foreach my $line (@kidsarray)
					{ 	
						for my $newkey ( sort {$PFtree{$a}->{pos} <=> $PFtree{$b}->{pos}} keys %PFtree) 
						{
								if ($line eq $PFtree{$newkey}{lab}) # $PFtree{$PFkey} is father   $PFtree{newkey}  ist kind
								{
									print O $PFtree{$newkey}{name}, "[ label=\"",$PFtree{$newkey}{name}, " k= " ,$PFtree{$newkey}{m}, "\"];\n";
									#print O $PFtree{$newkey}{name}, "[ label=\"{{",$PFtree{$newkey}{name},"| k = ",$PFtree{$newkey}{m},"}}\"];\n";
																	
									#tad[label=<<font color="darkgreen">middle</font>       <font color="green">middle</font>>]; 
									print O $PFtree{$PFkey}{name},  " -- ", $PFtree{$newkey}{name}  , "[label=<<font color=\"darkgreen\">gain:", $PFtree{$newkey}{gain},"</font>","<br/>", "<font color=\"red\">loss:",$PFtree{$newkey}{loss} ,"</font>>];\n"; 
									
									#print O $PFtree{$PFkey}{name},  " -- ", $PFtree{$newkey}{name}, "[label=\"gain: ", $PFtree{$newkey}{gain},"  loss: ", $PFtree{$newkey}{loss}, " \"]; \n"; 
									
								} # [shape=record, label="{{species b}|{gain: 1}}"];
								#"[shape=record, label=\"{{", $PFtree{$PFkey}{name},"}|{k=",$PFtree{$PFkey}{m},"}}\"]
						}
					
					}	
					
				#printf O "%s\n",$FamtreePF{$PFkey}{name};
				#print(" Name:  ", $FamtreePF{$PFkey}{name}, "  labs --> ", $PFtree{$PFkey}{lab} ,"\n" );
				
				}
			} 
			# jetzt müssen noch alle blätte rauf ein level gesetzt werden:
			
			print O "{rank = same;";
			for my $PFkey ( sort {$PFtree{$a}->{pos} <=> $PFtree{$b}->{pos}} keys %PFtree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
			{			
				if ($PFtree{$PFkey}{c} eq '-1')
				{
					print O $PFtree{$PFkey}{name},";";
				}
			}
			print O "}\n";	
			
				printf O ("}");
				printf("\nThe DOT - tree-file is named: %s\n\n", $outfilePF);
				close(O);
		
		
			####################################################
			
			my $file;
			
				
				#print "innerhalb perl\n";
				#`echo Außerhalb von perl\n`;
			#system("echo ergebnis von pwd");
			if ($Suffix eq "PF.out")
			{
				system("dot -Tsvg PF.out_dot_tree.dot > PF_output_mit_rank.svg");
				$file = "PF_output_mit_rank.svg";
			
			}
			elsif ($Suffix eq "PFb.out")
			{
				system("dot -Tsvg PFb.out_dot_tree.dot > PFb_output_mit_rank.svg");
				$file = "PFb_output_mit_rank.svg";
			}	
			elsif ($Suffix eq "PSb.out")
			{
				system("dot -Tsvg PSb.out_dot_tree.dot > PSb_output_mit_rank.svg");
				$file = "PSb_output_mit_rank.svg";
			}	
		
			print "\n\n";
			my @firstkeys;
			my @secondkeys;
			my @thirdkeys;
	
			for my $key_all_1 ( keys %Ges_Tree) 
			{
				#print "\$key_all_1: " ,$key_all_1, "\n\n" ;
				push (@firstkeys,$key_all_1);
				for my $key_all_2(keys %{$Ges_Tree{$key_all_1}})
				{
					#print "\$key_all_2: " ,$key_all_2, "\n\n" ;
					push (@secondkeys,$key_all_2);
				}
				
			}	
			
			#~ my $firstlength=@firstkeys;
			#~ my $secondlength=@secondkeys;
			#~ my $methodlength= $secondlength / $firstlength  ;
			#~ print "\@firstkeys:  ",@firstkeys, "länge: ", $firstlength,   "\n";
			#~ print "\@secondkeys:  ",@secondkeys, "länge: ", $secondlength, "\n";
			#~ print "länge je methode:  ", $methodlength;
			
			print "\n\n";
			
		#	print " PF tree: \n", Dumper(\%PFtree);
			
			$outfilePF = "all_dot_tree.dot";
			open(O, ">", $outfilePF) or die $!;
			
			printf O ("graph graphname{\n\nnode [ fontsize=25,shape=record,nodesep=0,];edge [fontsize=25, decorate=true, penwidth=3];graph [ranksep=0];rankdir=LR;ratio=\"auto\";\n\n");
				
			for my $PFkey ( sort {$PFtree{$a}->{pos} <=> $PFtree{$b}->{pos}} keys %PFtree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
			{			
				if ($PFtree{$PFkey}{count} ==0){$PFtree{$PFkey}{count} =1;}
				if ($PFtree{$PFkey}{c} ne '-1')
				{				
					my @kidsarray = split /\,/,$PFtree{$PFkey}{c};
					#print join(", ",@kidsarray, "\n");					#für jedes element in @kidsarray wird gesamte baum durchsucht wenn $kidsarray[$p] eq lab print $PFtree{$PFkey}{c} -> $PFtree{$ii}{name} 
					# }}"];
					#Obscuragr [ label=" {{ Obscuragr |   3   3    2 }}"];
					
					print O $PFtree{$PFkey}{name},"[label=<<font color=\"black\">",$PFtree{$PFkey}{name},"</font>","<br/>",
					"<font color=\"red\">",$Ges_Tree{'PSb.out'}{$PFkey}{m},"</font>","<br/>",
					"<font color=\"blue\">",$Ges_Tree{'PFb.out'}{$PFkey}{m},"</font>","<br/>",
					"<font color=\"darkgreen\">",$Ges_Tree{'PF.out'}{$PFkey}{m},"</font>>];\n";  #  T [shape=record, label="{T|000}"];
					foreach my $line (@kidsarray)
					{ 	
						for my $newkey ( sort {$PFtree{$a}->{pos} <=> $PFtree{$b}->{pos}} keys %PFtree) 
						{
							if ($line eq $PFtree{$newkey}{lab}) # $PFtree{$PFkey} is father   $PFtree{newkey}  ist kind
							{
								if ($PFtree{$newkey}{c} eq "-1") #at leaves print k
								{
									if ($PFtree{$newkey}{m} eq "0" )
									{
										print O $PFtree{$newkey}{name}, "[fontcolor=\"grey\", color=\"grey\" ,label=\"",$PFtree{$newkey}{name}, " k= " ,$PFtree{$newkey}{m} , "\"];\n";
									}
									else
									{	
									print O $PFtree{$newkey}{name}, "[label=\"",$PFtree{$newkey}{name}, " k= " ,$PFtree{$newkey}{m} , "\"];\n";
									}
								}
								#print O $PFtree{$newkey}{name}, "[ label=\"{{",$PFtree{$newkey}{name},"| k = ",$PFtree{$newkey}{m},"}}\"];\n";
																
								#tad[label=<<font color="darkgreen">middle</font>       <font color="green">middle</font>>]; 
								print O $PFtree{$PFkey}{name},  " -- ", $PFtree{$newkey}{name}  ,
								 "[label=<<font color=\"black\">gain:</font>", 
								 "<font color=\"red\">", $Ges_Tree{'PSb.out'}{$newkey}{gain},"</font>",
								 "<font color=\"blue\">" ,$Ges_Tree{'PFb.out'}{$newkey}{gain},"</font>",
								 "<font color=\"darkgreen\">",$Ges_Tree{'PF.out'}{$newkey}{gain},"</font>","<br/>",
								 "<font color=\"black\">loss:</font>",
								 "<font color=\"red\">",$Ges_Tree{'PSb.out'}{$newkey}{loss},"</font>",
								 "<font color=\"blue\">",$Ges_Tree{'PFb.out'}{$newkey}{loss},"</font>",
								 "<font color=\"darkgreen\">",$Ges_Tree{'PF.out'}{$newkey}{loss},"</font>>];\n"; 
								
								#print O $PFtree{$PFkey}{name},  " -- ", $PFtree{$newkey}{name}, "[label=\"gain: ", $PFtree{$newkey}{gain},"  loss: ", $PFtree{$newkey}{loss}, " \"]; \n"; 
								
							} # [shape=record, label="{{species b}|{gain: 1}}"];
						}
					}	
				}
			} 
			
			#~ sub fontcolor
			#~ {	
				#~ if ($Ges_Tree{'PFb.out'}{$newkey}{gain} eq "0"){ return "<font color=\"grey\">" }else{ return "<font color=\"blue\">" }
				#~ print $_;
				
			#~ }
			# jetzt müssen noch alle blätte rauf ein level gesetzt werden:
			
			print O "{rank = same;";
			for my $PFkey ( sort {$PFtree{$a}->{pos} <=> $PFtree{$b}->{pos}} keys %PFtree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
			{			
				if ($PFtree{$PFkey}{c} eq '-1')
				{
					print O $PFtree{$PFkey}{name},";";
				}
			}
					
			print O "}\n";	
			print O "Legend[ label=<<font color=\"red\">Sankoff</font><br/><font color=\"blue\">Partition Function</font><br/><font color=\"darkgreen\">Partition Function -P</font>>] \n ";
				printf O ("}");
				print("\nThe DOT - tree-file is named: all_dot_tree.dot\n\n");
				close(O);
		
	
	
	#~ ####################################################
	
			
				#~ #print "innerhalb perl\n";
				#~ #`echo Außerhalb von perl\n`;
				#~ #system("echo ergebnis von pwd");
			 system("dot -Tsvg all_dot_tree.dot > all_output_mit_rank.svg");
		
	
	
	
	
	
	
	
	
	
			####################################################### öffnen der erstellten svg um koordinaten zu erhalten
			#print "öffnen der output_mit_rank.svg \n\n";
			#~ print "Suffix:", Dumper($Suffix), "\n\n";
			#~ print Dumper($file), "\n\n";
			my @svgzeilen;
			#my $count=0;		
			#open(E, $file) or die $!;
			open E,'<',$file or die $!;
			while(<E>) 
			{
				chomp();
				#my $line = $_;
				#print $_;
				push(@svgzeilen,$_);
				#my @svglines[$count]=$_
				#grep 
				#my @x = split /\s+/,$line;
				#$count=$count +1;
			}
			close(E);	#    next;
			
			#my $svgzeilen = @svgzeilen;
			#print $svgzeilen, "\n";
			my $counter =0;
			my %nodeorder=();
			my %nodeorderx=();
			my $height;
			my $width	;
			foreach my $p (@svgzeilen)
			{
				if (grep  {/class="node"/} $p)		#knoten und deren y koordinate werden festgehalten
				{
					#print $p,"\n";
					
					(my $nodename)=$svgzeilen[1+$counter] =~/<title>(.*)<\/title>/;
					#my $nodename= $&;
					(my $y)=$svgzeilen[3+$counter] =~/y="(.*)/;
					$y =~/" font.*/;
					$y= $`;
					
					(my $x)=$svgzeilen[3+$counter] =~/x="(.*)/;
					$x =~/" font.*/;
					$x= $`;
					
					
					#print "nodename= ",$nodename, "  y= ", $y, "\n";
					#print "nodename= ",$nodename, "  x= ", $x, "\n";
					$nodeorder{$nodename}=$y;
					$nodeorderx{$nodename}=$x;
					
					
					#print "nodename:",$nodename,"\n";
					#print "sind in zeile ", $p,"\n";
					#print "nächste zeile: ", $svgzeilen[$counter+1],"\n"; 
					#print "nächste zeile: ", $svgzeilen[$counter+1] =~,"\n"; 
				} 
				$counter++;
				
				if (grep  {/height/} $p)		#höhe des baus wird festgehalten
				{
					#print $p,"\n";
					($height)=$p =~/height="(.*)pt/;
					print "höhe: ", $height,"\n";
					
					($width)=$p =~/width="(.*)pt" /;
					print "breite: ", $width,"\n";
					#my $height=
					
					
				}
				
				
			}
			
			#print Dumper(%nodeorder);
			
							
			my @keys=sort{$nodeorder{$a} <=> $nodeorder{$b} or $nodeorderx{$a} <=> $nodeorderx{$b}} keys (%nodeorder);  #sortiert nach value
			my $count=@keys;
			for (@keys)					#koordianten werden mit int werten überschrieben Reihenfolge muss umgekehrt werden, da in R wiederum sonst falsche ausgabe:
			{
				$nodeorder{$_}=$count;
				$count--;
				#print("key: ", $_ ,"\tvalue: ",$nodeorder{$_},  "\n" );
			}
	
			#~ print("\nvalues sortiert:\n\n");
			#~ for (@keys)
			#~ {
				#~ #$nodeorder{$_}=$count;
				#~ print("key: ", $_ ,"\tvalue: ",$nodeorder{$_},  "\n" );
			#~ }
			
			
			#########################################  bestehende _PFSUM.out wird geöffnet und in fertige 
			
			############## 
			my $outfilePF ;
				if( $outfile =~ /\./)
				{
					$outfilePF = substr($outfile,0, rindex( $outfile, '.') );			#rindex STR,SUBSTR
					$outfilePF = $outfilePF.$Suffix;				#Works just like index except that it returns the position of the last occurrence of SUBSTR in STR. If POSITION is specified,returns the last occurrence beginning at or before that position.
		
					
				}
				else
				{$outfilePF = $outfile.$Suffix;}		
					
			print "geöffnet wird: ", $sum_PF, "  aber outfile= ", $outfilePF, "\n";		

	
	
			print "einlesen der sum", $sum_PF, "\n";
	
			my @sum_zeilen;
			#my $count=0;		
			#open(E, $file) or die $!;				#erneutes einlesen der .sum, mit alter .sum muss svg gebaut werden um koordinaten zu bekommen, 
													#dann wieder einlesen des baumes um koordinaten nachzutragen
			open E,'<',$sum_PF or die $!;
			while(<E>) 
			{
				chomp();
				#my $line = $_;
				#print $_;
				push(@sum_zeilen,$_);
				#my @svglines[$count]=$_
				#grep 
				#my @x = split /\s+/,$line;
				#$count=$count +1;
			}
			close(E);
			
	#print Dumper(@sum_zeilen);	#datei eingelesen
		
			my %final_sum;
			my $counter=0;
			for my $i (@sum_zeilen)				#baum wird in hash gelegt
			{
				my @sum_spalten = split /\s+/,$i;		#in zeile @x sind spalten $x[] 
				$final_sum{$counter}{lab} 			= $sum_spalten[3]; 
				$final_sum{$counter}{kids} 			= $sum_spalten[5];
				$final_sum{$counter}{parents} 		= $sum_spalten[9];
				$final_sum{$counter}{genes} 		= $sum_spalten[11];
				$final_sum{$counter}{gain} 			= $sum_spalten[13];
				$final_sum{$counter}{loss} 			= $sum_spalten[15];
				$final_sum{$counter}{Nr_Fams} 		= $sum_spalten[21];
				$final_sum{$counter}{gainFam} 		= $sum_spalten[17];
				$final_sum{$counter}{lossFam} 		= $sum_spalten[19];
				$final_sum{$counter}{PF_P} 			= $sum_spalten[23];
				$final_sum{$counter}{species} 		= $sum_spalten[-1]; 
				$final_sum{$counter}{Printorder}	= 0;
				$counter++;
			}	
		
		
			#print Dumper(%nodeorder);
			
			#print Dumper(\%final_sum);
			$counter=0;
		
			for ($counter;$counter< scalar @keys ;$counter++)		#jeweiligen species wird koordinate zugeordnet
			{
				for (@keys)
				{
					if ($_ eq $final_sum{$counter}{species})
					{	#print "keys= ",$_," \$final_sum{\$counter}{species} = ", $final_sum{$counter}{species}  ,"\n";
						
						$final_sum{$counter}{Printorder}=$nodeorder{$_}; 
						
					}
					
				}
			}
			
			#print Dumper(%final_sum);
			
			
			
			#~ $counter=0;
			#~ for my $i (@sum_zeilen)	
			#~ {
				#~ print Dumper(
				
				#~ $final_sum{$counter}{lab} 			,
				#~ $final_sum{$counter}{kids} 			,
				#~ $final_sum{$counter}{parents} 		,
				#~ $final_sum{$counter}{genes} 		,
				#~ $final_sum{$counter}{gain} 			,
				#~ $final_sum{$counter}{loss} 			,
				#~ $final_sum{$counter}{Nr_Fams} 		,
				#~ $final_sum{$counter}{gainFam} 		,
				#~ $final_sum{$counter}{lossFam} 		,
				#~ $final_sum{$counter}{PF_P} 			,
				#~ $final_sum{$counter}{species} 		,
				#~ $final_sum{$counter}{Printorder}	
				#~ );         
				#~ $counter++;
			#~ }
				
			################################# fertige lst_4R_ schreiben
					
						
			$outfilePF = $Suffix."_lst_4R.dat";
			
			print "lst_4R wird gespeichert in:", $outfilePF,"\n";
			
			my @sorted_entries 	= sort {$a->{Printorder} <=> $b->{Printorder}} values %final_sum; 			#.out liste wird nach koordianten sortiert
			#print Dumper(@sorted_entries);
					#print Dumper(\%final_sum);
			#~ $final_sum{$counter}{species} 
			
			
			
			open(O, ">", $outfilePF) or die $!;
			print O "lab\tkids\tparents\tgenes\tgain\tloss\tNr_Fams\tgainFam\tlossFam\tPF_P\tspecies\tPrintorder","\n";
			my $counter=0;
			
			
			for (@sorted_entries)
			{			 
				if (defined $sorted_entries[$counter]{lab})
				{
				print O (	$sorted_entries[$counter]{lab},"\t",$sorted_entries[$counter]{kids},"\t",$sorted_entries[$counter]{parents},"\t",
							$sorted_entries[$counter]{genes},"\t",$sorted_entries[$counter]{gain},"\t",$sorted_entries[$counter]{loss},"\t",
							$sorted_entries[$counter]{Nr_Fams},"\t",$sorted_entries[$counter]{gainFam},"\t",$sorted_entries[$counter]{lossFam},"\t",
							$sorted_entries[$counter]{PF_P},"\t",$sorted_entries[$counter]{species},"\t",$sorted_entries[$counter]{Printorder},"\n") ;
				}			
				$counter++;
			}	 
			close O;
			#~ open(O, ">", $outfilePF) or die $!;
			
		#~ printf O ("lab\tkids\tparent\tgenes\tgain\tloss\tNr_Fams\tgainFam\tlossFam\tPF_P\tspecies\tPrintorder\n");
		#~ $counter=0;
		#~ for (@sorted_entries)
		#~ {
			#~ print $sorted_entries[$counter]{species}, "\n";
			#~ print O $sorted_entries[$counter]{kids}, "\n";
				  
			#~ print O $sorted_entries[$counter]{lab} ;
			
			#~ ,"\t",	$sorted_entries[$counter]{kids} ,"\t", $sorted_entries{$counter}{parents} ,"\t", $sorted_entries[$counter]{genes} ,"\t",
						#~ $sorted_entries[$counter]{gain} ,"\t",$sorted_entries[$counter]{loss} ,"\t",$sorted_entries[$counter]{Nr_Fams} ,"\t", $sorted_entries[$counter]{gainFam} ,"\t",
						#~ $sorted_entries[$counter]{lossFam} , "\t",$sorted_entries[$counter]{PF_P} ,"\t", $sorted_entries[$counter]{species} ,"\t",
						#~ $sorted_entries[$counter]{Printorder} ,"\n";
			
			#~ $counter++;
		#~ }
		
		#~ close O;
		
		
		
		
		 
			
		#~ my $set_print=0;
		#~ my $counter=0;
		#~ my $len=@sum_zeilen;
		#~ my $counter=@sum_zeilen;
		#~ for ($len; $len>=0 ;$len-- )
		#~ {	for ($counter; $counter>=0 ;$counter-- )
			#~ {	print "len: ", $len, "counter: ", $counter,"\n";
				#~ if ($final_sum{$counter}{Printorder} == $set_print)
				#~ {	
					#~ print 	"final_sum{counter}{Printorder}", $final_sum{$counter}{Printorder}, " set_print",$set_print,"\n";
						#~ print O ($final_sum{$counter}{lab} ,"\t",$final_sum{$counter}{kids} ,"\t", $final_sum{$counter}{parents} ,"\t", $final_sum{$counter}{genes} ,"\t",
						#~ $final_sum{$counter}{gain} ,"\t",$final_sum{$counter}{loss} ,"\t",$final_sum{$counter}{Nr_Fams} ,"\t", $final_sum{$counter}{gainFam} ,"\t",
						#~ $final_sum{$counter}{lossFam} , "\t",$final_sum{$counter}{PF_P} ,"\t", $final_sum{$counter}{species} ,"\t",
						#~ $final_sum{$counter}{Printorder} ,"\n");
						#~ $counter=0;
						#~ next;
				#~ }	
				
			#~ }
			#~ $set_print++;
			#~ $counter=@sum_zeilen;
		#~ }
		
		#close(O);
		}	
				
	}

}














#   $Ges_Tree{$Suffix}{$ii}{gain} = $x[13];




#~ print "\n\n\n";
#print Dumper(\%Ges_Tree);

#foreach 
#


#~ foreach my $Suffix(@_)
#~ {	
	#~ for my $PFkey ( sort {$Ges_Tree{$Suffix}{$a}->{pos} <=> $Ges_Tree{$Suffix}{$b}->{pos}} keys %Ges_Tree{$Suffix})
	#~ {
		#~ print "sortierter Key:  ", $PFkey, "\n";
	#~ } 
	
#~ }


#print Dumper(%Ges_Tree);
#print "\n\n\n";

#print Dumper(\%Ges_Tree);

#~ for my $PFkey ( sort {$Ges_Tree{"PF.out"}{$a}->{pos} <=> $Ges_Tree{"PF.out"}{$b}->{pos}} keys %Ges_Tree)
#~ {
	#~ print $Ges_Tree{"PF.out"}{$PFkey}{pos}, "\n";
	
#~ }


#~ my @sorted_entries 					#.out liste wird nach koordianten sortiert
	#~ = sort {$a->{Printorder} <=> $b->{Printorder}} 
		   #~ values %final_sum;
#~ #print Dumper(@sorted_entries);
		   #~ #print Dumper(\%final_sum);
#$final_sum{$counter}{species} 
		


#~ my $first_key="PF.out";

#~ my $outfilePF = "Gesamt_dot_tree.dot";
#~ open(O, ">",  $outfilePF) or die $!;

 #~ printf O ("graph graphname{\n\nnode [ fontsize=25,shape=record,nodesep=0,];edge [fontsize=25, decorate=true, penwidth=3];graph [ranksep=0];rankdir=LR;ratio=\"auto\";\n\n");
	
#~ for my $PFkey ( sort {$Ges_Tree{"PF.out"}{$a}->{pos} <=> $Ges_Tree{"PF.out"}{$b}->{pos}} keys %Ges_Tree{"PF.out"}) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
#~ {			
	#~ if ($Ges_Tree{$first_key}{$PFkey}{count} ==0){$Ges_Tree{$first_key}{$PFkey}{count} =1;}
	#~ if ($Ges_Tree{$first_key}{$PFkey}{c} ne '-1')
	#~ {				
		#~ my @kidsarray = split /\,/,$Ges_Tree{$first_key}{$PFkey}{c};
		#~ #print join(", ",@kidsarray, "\n");					#für jedes element in @kidsarray wird gesamte baum durchsucht wenn $kidsarray[$p] eq lab print $Ges_Tree{$first_key}{$PFkey}{c} -> $Ges_Tree{$first_key}{$ii}{name} 
		#~ print O $Ges_Tree{$first_key}{$PFkey}{name}, "[ label=\" {{",$Ges_Tree{$first_key}{$PFkey}{name},"| k = ",$Ges_Tree{$first_key}{$PFkey}{m} ,"}}\"];\n";  #  T [shape=record, label="{T|000}"];
		#~ foreach my $line (@kidsarray)
		#~ { 	
			#~ for my $newkey ( sort {$Ges_Tree{$first_key}{$a}->{pos} <=> $Ges_Tree{$first_key}{$b}->{pos}} keys %Ges_Tree) 
			#~ {
					#~ if ($line eq $Ges_Tree{$first_key}{$newkey}{lab}) # $Ges_Tree{$first_key}{$PFkey} is father   $Ges_Tree{$first_key}{newkey}  ist kind
					#~ {
						#~ print O $Ges_Tree{$first_key}{$newkey}{name}, "[ label=\"",$Ges_Tree{$first_key}{$newkey}{name}, " k= " ,$Ges_Tree{$first_key}{$newkey}{m}, "\"];\n";
						#~ #print O $Ges_Tree{$first_key}{$newkey}{name}, "[ label=\"{{",$Ges_Tree{$first_key}{$newkey}{name},"| k = ",$Ges_Tree{$first_key}{$newkey}{m},"}}\"];\n";
														
						#~ #tad[label=<<font color="darkgreen">middle</font>       <font color="green">middle</font>>]; 
						#~ print O $Ges_Tree{$first_key}{$PFkey}{name},  " -- ", $Ges_Tree{$first_key}{$newkey}{name}  , "[label=<<font color=\"darkgreen\">gain:", $Ges_Tree{$first_key}{$newkey}{gain},"</font>","<br/>", "<font color=\"red\">loss:",$Ges_Tree{$first_key}{$newkey}{loss} ,"</font>>];\n"; 
						
						#~ #print O $Ges_Tree{$first_key}{$PFkey}{name},  " -- ", $Ges_Tree{$first_key}{$newkey}{name}, "[label=\"gain: ", $Ges_Tree{$first_key}{$newkey}{gain},"  loss: ", $Ges_Tree{$first_key}{$newkey}{loss}, " \"]; \n"; 
						
					#~ } # [shape=record, label="{{species b}|{gain: 1}}"];
					#~ #"[shape=record, label=\"{{", $Ges_Tree{$first_key}{$PFkey}{name},"}|{k=",$Ges_Tree{$first_key}{$PFkey}{m},"}}\"]
			#~ }
		
		#~ }	
		
	#~ #printf O "%s\n",$FamtreePF{$PFkey}{name};
	#~ #print(" Name:  ", $FamtreePF{$PFkey}{name}, "  labs --> ", $Ges_Tree{$first_key}{$PFkey}{lab} ,"\n" );
	
	#~ }
#~ } 
#~ # jetzt müssen noch alle blätte rauf ein level gesetzt werden:

#~ print O "{rank = same;";
#~ for my $PFkey ( sort {$Ges_Tree{$first_key}{$a}->{pos} <=> $Ges_Tree{$first_key}{$b}->{pos}} keys %Ges_Tree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
#~ {			
	#~ if ($Ges_Tree{$first_key}{$PFkey}{c} eq '-1')
	#~ {
		#~ print O $Ges_Tree{$first_key}{$PFkey}{name},";";
	#~ }
#~ }
#~ print O "}\n";	

	#~ printf O ("}");
	#~ printf("\nThe DOT - tree-file is named: %s\n\n", $outfilePF);
	#~ close(O);















############################################################## searching for _PS.out files


		#~ print "\nsearching for output files from ePoPE in Parsimony mode in: $inDir \n(Files with an '_PS.out' -extension)\n\n";
	
		#~ # read directory of data
			
		#~ opendir(D, $inDir) or die $!;
		#~ my @entries = grep {/\_PS.out$/} readdir(D);
		#~ closedir(D);
		#~ my $PSlaenge=@entries;
	#~ if ($PSlaenge==0)
	#~ {
		#~ printf("No '_PS.out' Files found\n\n")
	#~ }
		
	#~ if ($PSlaenge!=0)
	#~ {
		#~ # fill tree hash with data
		#~ my %PStree;
		#~ my %Famtree;
		#~ my ($i,$nodes_n) = (0,0);
		#~ foreach my $e (@entries) 
		#~ {
			#~ my $file = $inDir."".$e;
			#~ printf "%s\n",$file; 
			#~ open(E, $file) or die $!;
			
			#~ while(<E>) 
			#~ {
				#~ chomp();
				#~ my $line = $_;
			#~ #   print $line;
				#~ my @x = split /\s+/,$line;   # $line ganze eingelesene zeile,  @x array bestehend aus elementen der zeile, welche an leerzeichen getrennt wurden
			#~ #    next;
				#~ $i = $x[3]; #/* use i as identifier and lab of tree */ -->lab wise sort
				#~ # write first instance in tree
				
				#~ if(! exists($PStree{$i})) 
				#~ {
					#~ $PStree{$i}{name} = $x[scalar @x -1]; chomp($PStree{$i}{name}); #letzter eintrag ist wohl species name, sollte demnach passen...
					#~ $PStree{$i}{pos} = $x[1];
					#~ $PStree{$i}{lab} = $x[3];
					#~ $PStree{$i}{p} = $x[9];
					#~ $PStree{$i}{c} = $x[5];
					#~ $PStree{$i}{m} = $x[11];
					#~ $PStree{$i}{gain} = $x[13];
					#~ $PStree{$i}{loss} = $x[15];
					#~ $PStree{$i}{gainF} = $x[17];
					#~ $PStree{$i}{lossF} = $x[19];
					#~ $nodes_n++;
				#~ }
				#~ else 
				#~ { # add values to existing tree
					#~ $PStree{$i}{m} += $x[11];
					#~ $PStree{$i}{gain} += $x[13];
					#~ $PStree{$i}{loss} += $x[15];
					#~ $PStree{$i}{gainF} += $x[17];
					#~ $PStree{$i}{lossF} += $x[19];
				#~ }
				
				#~ if(! exists($Famtree{$i}))
				#~ {
					#~ $Famtree{$i}{p} = $x[9]; # label of parent node
					#~ $Famtree{$i}{gainF}  = $x[17]; # gainFam
					#~ $Famtree{$i}{lossF}  = $x[19]; # lossFam
					#~ $Famtree{$i}{name}= $x[20]; # name of node
					#~ $Famtree{$i}{Nfam}=0;
				#~ }             
				#~ else          
				#~ {             
					#~ $Famtree{$i}{gainF}  += $x[17]; # gainFam
					#~ $Famtree{$i}{lossF}  += $x[19]; # lossFam
				#~ }
			#~ }
			#~ close(E);
		#~ }
	
	 	#~ my ($Nfam);
		#~ sub traceIt;
		
		#~ foreach my $i (sort {$a<=>$b} keys %Famtree)       #The spaceship <=> operator is used for comparing numbers, here sort for lab (reihenfolge eigtl egal, es darf nur kein schlüssel gewählt werden den es nicht gibt...) ( use cmp for literals, <=> fpr numerals )
		#~ {
		
			#~ $Famtree{$i}{Nfam} = traceIt($i, 0, %Famtree); # traceIt(label,Nfam,%tree) $i ist hier lab number da beim einlesen diese spalte für die cshlüssel verwendet wird
			 #~ #print "label: ", $i,"\t name:",$Famtree{$i}{name}, " \t  lossF: ",$Famtree{$i}{lossF}," gainF: ", $Famtree{$i}{gainF},  "  parent: ", $Famtree{$i}{p} , "\t Nfam:  ", $Famtree{$i}{Nfam} ,  "\n"  ;
		
		#~ }	
		
		#~ my $outfilePS;
		
		#~ if( $outfile =~ /\./)
		#~ {
			#~ $outfilePS = substr($outfile,0, rindex( $outfile, '.') );
			#~ $outfilePS = $outfilePS."_PS.sum";
			
		#~ }
		#~ else
		#~ {	$outfilePS = $outfile."_PS.sum";}
		
		#~ open(O, ">", $outfilePS) or die $!;
		#~ #for($i=0; $i<$nodes_n; $i++) 
		#~ #	{
				#~ #um kompatibel mit ePoPe -c zu sein, muss in pos reihefolge gelistet werden
				#~ #sortierung nach $PStree{$i}{pos} = $x[1];
			
		#~ for my $PSkey ( sort {$PStree{$a}->{pos} <=> $PStree{$b}->{pos}} keys %PStree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
		#~ {			
			#~ #my $value = $PStree{$PSkey};
			#~ printf O "pos: %d lab: %d kids: %s pInP: %d par: %d ",$PStree{$PSkey}{pos},$PStree{$PSkey}{lab},$PStree{$PSkey}{c},$PSkey,$PStree{$PSkey}{p};
			#~ printf O "genes: %d gain: %d loss: %d ",$PStree{$PSkey}{m},$PStree{$PSkey}{gain},$PStree{$PSkey}{loss};
			#~ #print "the NFam of ",$PStree{$PSkey}{name} , " is ", $Famtree{$PSkey}{Nfam} , " (also : " ,$Famtree{$PSkey}{name}, ") \$nodes_n=", $nodes_n, " \n";
			#~ printf O "gainFam: %d lossFam: %d N_Fam: %d %s\n",$PStree{$PSkey}{gainF},$PStree{$PSkey}{lossF},$Famtree{$PSkey}{Nfam},$Famtree{$PSkey}{name};
		#~ }
					#~ # printf O "count: %d\n", $tree{$i}{count}; # test count
					#~ #printf("famtest: species: %s NFAM: %d\n", $PStree{$j}{name} , $PStree{$j}{Nfam} )
		
		
		#~ printf("\nThe Summary File is named: %s\n\n", $outfilePS);
		#~ close(O);
	#~ }	
	
	
	
	
	
	
	
	
	
			
	
	#~ ######################## now do the same for the files with the _PSb extension:
		#~ print "\nsearching for output files from ePoPE in both mode in: $inDir \n(Files with an '_PSb.out' -extension)\n\n";
				
		#~ opendir(D, $inDir) or die $!;
		#~ my @PSbentries = grep {/\_PSb.out$/} readdir(D);
		#~ closedir(D);
		#~ $laenge=@PSbentries;
	#~ if ($laenge==0)
	#~ {
		#~ printf("No '_PSb.out' Files found\n\n")
	#~ }	
		
		
	#~ if ($laenge!=0)
	#~ {
		#~ #printf("länge von entries: %d\n", $laenge);
		
		#~ # fill tree hash with data
		#~ my %PSbtree;
		#~ my %FamtreePSb;
		#~ my ($ii,$Nnodes_n) = (0,0);
		#~ foreach my $e (@PSbentries) 
		#~ {
			#~ my $file = $inDir."".$e;
			#~ printf "%s\n",$file; 
			#~ open(E, $file) or die $!;
			
			#~ while(<E>) 
			#~ {
				#~ chomp();
				#~ my $line = $_;
			#~ #   print $line;
				#~ my @x = split /\s+/,$line;
			#~ #    next;
			
			#~ #pos:  16 lab: 329 kids: -1 pInP:  16 par: 327 mirs:   0 gain:   0 loss   0       seenBW:   1 gainFam:   0        lossFam:   0    egr
			
				#~ $ii = $x[3]; 
				#~ my $count = 0; 
				#~ # write first instance in tree
				#~ if(! exists($PSbtree{$ii})) 
				#~ {
					#~ $PSbtree{$ii}{name} = $x[scalar @x -1]; chomp($PSbtree{$ii}{name}); #letzter eintrag ist wohl species name, sollte demnach passen...
					#~ $PSbtree{$ii}{pos} = $x[1];
					#~ $PSbtree{$ii}{lab} = $x[3];
					#~ $PSbtree{$ii}{p} = $x[9];
					#~ $PSbtree{$ii}{c} = $x[5];
					#~ $PSbtree{$ii}{m} = $x[11];
					#~ $PSbtree{$ii}{gain} = $x[13];
					#~ $PSbtree{$ii}{loss} = $x[15];
					#~ $PSbtree{$ii}{gainF} = $x[17];
					#~ $PSbtree{$ii}{lossF} = $x[19];
					#~ if($PSbtree{$ii}{m} == 0) #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit
					#~ {
						#~ $PSbtree{$ii}{count} = 0;
						#~ $PSbtree{$ii}{pf_P} = 0;
						
					#~ }
					#~ else
					#~ {
						#~ $PSbtree{$ii}{count} = 1;		
						#~ $PSbtree{$ii}{pf_P} = $x[21];
						#~ #print "m von tree.$i ist $x[11]   zu addierende PF_P: $x[23] aus datei $e\n";
					#~ }
					#~ $Nnodes_n++;
				#~ } 	
				#~ else 
				#~ { # add values to existing tree
					#~ $PSbtree{$ii}{m} += $x[11];
					#~ $PSbtree{$ii}{gain} += $x[13];
					#~ $PSbtree{$ii}{loss} += $x[15];
					#~ $PSbtree{$ii}{gainF} += $x[17];
					#~ $PSbtree{$ii}{lossF} += $x[19];
					#~ if($x[11] != 0 && $x[5] ne '-1') #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit//bei kids($x[5]) die mit -1 auslassen ist aber char!, also stringcompare
					#~ {
						#~ #printf( "m von tree.[%d] ist %d   zu addierende PF_P: %f aus datei %s mit namen %s\n",$ii,$x[11],$x[21],$e,$x[scalar @x -1]);
						#~ $PSbtree{$ii}{pf_P} += $x[21];  # um relative wahrscheinlichkeit auszurechnen, pf_P mitzählen wenn m über 0 (wenn 0, dann ist knoten außerhalb des alignments)
						#~ $PSbtree{$ii}{count} += 1;
					#~ }
				#~ }

				#~ if(! exists($FamtreePSb{$ii}))
				#~ {
					#~ $FamtreePSb{$ii}{p} = $x[9]; # label of parent node
					#~ $FamtreePSb{$ii}{gainF}  = $x[17]; # gainFam
					#~ $FamtreePSb{$ii}{lossF}  = $x[19]; # lossFam
					#~ $FamtreePSb{$ii}{name}= $x[22]; # name of node
					#~ $FamtreePSb{$ii}{Nfam}=0;
					
				#~ }
				#~ else
				#~ {
					#~ $FamtreePSb{$ii}{gainF}  += $x[17]; # gainFam
					#~ $FamtreePSb{$ii}{lossF}  += $x[19]; # lossFam
				#~ }
				
			#~ }
			#~ close(E);
		#~ }
			
		#~ foreach my $i (sort {$a<=>$b} keys %FamtreePSb)       #The spaceship <=> operator is used for comparing 
		#~ {
			
				#~ $FamtreePSb{$i}{Nfam} = traceIt($i, 0, %FamtreePSb);
				#~ #print "label: ", $i,"\t name:",$FamtreePSb{$i}{name}, " \t  lossF: ",$FamtreePSb{$i}{lossF}," gainF: ", $FamtreePSb{$i}{gainF},  "  parent: ", $FamtreePSb{$i}{p} ,"\t Nfam: ", $FamtreePSb{$i}{Nfam} ,  "\n"  ;
			
		#~ }	
		
		#~ my $outfilePSb;
		#~ if( $outfile =~ /\./)
		#~ {
			#~ $outfilePSb = substr($outfile,0, rindex( $outfile, '.') );
			#~ $outfilePSb = $outfilePSb."_PSb.sum";
			
		#~ }
		#~ else
		#~ { $outfilePSb = $outfile."_PSb.sum";}
			
		#~ open(O, ">", $outfilePSb) or die $!;
		
		#~ for my $PSbkey ( sort {$PSbtree{$a}->{pos} <=> $PSbtree{$b}->{pos}} keys %PSbtree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
		#~ {			
			#~ #my $value = $PStree{$PSkey};
			#~ printf O "pos: %d lab: %d kids: %s pInP: %d par: %d ",$PSbtree{$PSbkey}{pos},$PSbtree{$PSbkey}{lab},$PSbtree{$PSbkey}{c},$PSbkey,$PSbtree{$PSbkey}{p};
			#~ printf O "genes: %d gain: %d loss: %d ",$PSbtree{$PSbkey}{m},$PSbtree{$PSbkey}{gain},$PSbtree{$PSbkey}{loss};
			#~ #print "the NFam of ",$PSbtree{$PSbkey}{name} , " is ", $Famtree{$PSbkey}{Nfam} , " (also : " ,$Famtree{$PSbkey}{name}, ") \$nodes_n=", $nodes_n, " \n";
			#~ if ($PSbtree{$PSbkey}{count} ==0){$PSbtree{$PSbkey}{count} =1;}
			#~ printf O "gainFam: %d lossFam: %d N_Fam: %d PF_P: %f %s\n",$PSbtree{$PSbkey}{gainF},$PSbtree{$PSbkey}{lossF},$FamtreePSb{$PSbkey}{Nfam},$PSbtree{$PSbkey}{pf_P}/$PSbtree{$PSbkey}{count},$FamtreePSb{$PSbkey}{name};
		#~ }
		
				
		#~ printf("\nThe Summary File is named: %s\n\n", $outfilePSb);
		#~ close(O);
		
	#~ }
		
	#~ ######################## now do the same for the files with the _PFb extension:
		#~ print "\nsearching for output files from ePoPE in both mode in: $inDir \n(Files with an '_PFb.out' -extension)\n\n";
				
		#~ opendir(D, $inDir) or die $!;
		#~ my @PFbentries = grep {/\_PFb.out$/} readdir(D);
		#~ closedir(D);
		#~ $laenge=@PFbentries;
	#~ if ($laenge==0)
	#~ {
		#~ printf("No '_PFb.out' Files found\n\n")
	#~ }	
			
	#~ if ($laenge!=0)
	#~ {
		#~ #printf("länge von entries: %d\n", $laenge);
		
		#~ # fill tree hash with data
		#~ my %PFbtree;
		#~ my %FamtreePFb;
		#~ my ($ii,$Nnodes_n) = (0,0);
		#~ foreach my $e (@PFbentries) 
		#~ {
			#~ my $file = $inDir."".$e;
			#~ printf "%s\n",$file; 
			#~ open(E, $file) or die $!;
			
			#~ while(<E>) 
			#~ {
				#~ chomp();
				#~ my $line = $_;
			#~ #   print $line;
				#~ my @x = split /\s+/,$line;
			#~ #    next;
			
			#~ #pos:  16 lab: 329 kids: -1 pInP:  16 par: 327 mirs:   0 gain:   0 loss   0       seenBW:   1 gainFam:   0        lossFam:   0    egr
			
				#~ $ii = $x[3]; 
				#~ my $count = 0; 
				#~ # write first instance in tree
				#~ if(! exists($PFbtree{$ii})) 
				#~ {
					#~ $PFbtree{$ii}{name} = $x[scalar @x -1]; chomp($PFbtree{$ii}{name}); #letzter eintrag ist wohl species name, sollte demnach passen...
					#~ $PFbtree{$ii}{pos} = $x[1];
					#~ $PFbtree{$ii}{lab} = $x[3];
					#~ $PFbtree{$ii}{p} = $x[9];
					#~ $PFbtree{$ii}{c} = $x[5];
					#~ $PFbtree{$ii}{m} = $x[11];
					#~ $PFbtree{$ii}{gain} = $x[13];
					#~ $PFbtree{$ii}{loss} = $x[15];
					#~ $PFbtree{$ii}{gainF} = $x[17];
					#~ $PFbtree{$ii}{lossF} = $x[19];
					#~ if($PFbtree{$ii}{m} == 0) #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit
					#~ {
						#~ $PFbtree{$ii}{count} = 0;
						#~ $PFbtree{$ii}{pf_P} = 0;
						
					#~ }
					#~ else
					#~ {
						#~ $PFbtree{$ii}{count} = 1;		
						#~ $PFbtree{$ii}{pf_P} = $x[21];
						#~ #print "m von tree.$i ist $x[11]   zu addierende PF_P: $x[23] aus datei $e\n";
					#~ }
					#~ $Nnodes_n++;
				#~ } 	
				#~ else 
				#~ { # add values to existing tree
					#~ $PFbtree{$ii}{m} += $x[11];
					#~ $PFbtree{$ii}{gain} += $x[13];
					#~ $PFbtree{$ii}{loss} += $x[15];
					#~ $PFbtree{$ii}{gainF} += $x[17];
					#~ $PFbtree{$ii}{lossF} += $x[19];
					#~ if($x[11] != 0 && $x[5] ne '-1') #sollte alignment außerhalb des inneren knoten liegen hier 0 --> knoten hat keinen einfluss auf wahrscheinlichkeit//bei kids($x[5]) die mit -1 auslassen ist aber char!, also stringcompare
					#~ {
						#~ #printf( "m von tree.[%d] ist %d   zu addierende PF_P: %f aus datei %s mit namen %s\n",$ii,$x[11],$x[21],$e,$x[scalar @x -1]);
						#~ $PFbtree{$ii}{pf_P} += $x[21];  # um relative wahrscheinlichkeit auszurechnen, pf_P mitzählen wenn m über 0 (wenn 0, dann ist knoten außerhalb des alignments)
						#~ $PFbtree{$ii}{count} += 1;
					#~ }
				#~ }

				#~ if(! exists($FamtreePFb{$ii}))
				#~ {
					#~ $FamtreePFb{$ii}{p} = $x[9]; # label of parent node
					#~ $FamtreePFb{$ii}{gainF}  = $x[17]; # gainFam
					#~ $FamtreePFb{$ii}{lossF}  = $x[19]; # lossFam
					#~ $FamtreePFb{$ii}{name}= $x[22]; # name of node
					#~ $FamtreePFb{$ii}{Nfam}=0;
					
				#~ }
				#~ else
				#~ {
					#~ $FamtreePFb{$ii}{gainF}  += $x[17]; # gainFam
					#~ $FamtreePFb{$ii}{lossF}  += $x[19]; # lossFam
				#~ }
				
			#~ }
			#~ close(E);
		#~ }
			
		#~ foreach my $i (sort {$a<=>$b} keys %FamtreePFb)       #The spaceship <=> operator is used for comparing 
		#~ {
			
				#~ $FamtreePFb{$i}{Nfam} = traceIt($i, 0, %FamtreePFb);
				#~ #print "label: ", $i,"\t name:",$FamtreePFb{$i}{name}, " \t  lossF: ",$FamtreePFb{$i}{lossF}," gainF: ", $FamtreePFb{$i}{gainF},  "  parent: ", $FamtreePFb{$i}{p} ,"\t Nfam: ", $FamtreePFb{$i}{Nfam} ,  "\n"  ;
			
		#~ }	
	
		#~ my $outfilePFb;
		#~ if( $outfile =~ /\./)
		#~ {
			#~ $outfilePFb = substr($outfile,0, rindex( $outfile, '.') );
			#~ $outfilePFb = $outfilePFb."_PFb.sum";
			
		#~ }
		#~ else
		#~ { $outfilePFb = $outfile."_PFb.sum";}
	
	
		#~ open(O, ">", $outfilePFb) or die $!;
		
		#~ for my $PFbkey ( sort {$PFbtree{$a}->{pos} <=> $PFbtree{$b}->{pos}} keys %PFbtree) #ausgabe nach pos in tree (lückenlos, pInP kann lücken haben) für späteren aufruf von ePoPE -c
		#~ {			
			#~ #my $value = $PStree{$PSkey};
			#~ printf O "pos: %d lab: %d kids: %s pInP: %d par: %d ",$PFbtree{$PFbkey}{pos},$PFbtree{$PFbkey}{lab},$PFbtree{$PFbkey}{c},$PFbkey,$PFbtree{$PFbkey}{p};
			#~ printf O "genes: %d gain: %d loss: %d ",$PFbtree{$PFbkey}{m},$PFbtree{$PFbkey}{gain},$PFbtree{$PFbkey}{loss};
			#~ #print "the NFam of ",$PFbtree{$PFbkey}{name} , " is ", $Famtree{$PFbkey}{Nfam} , " (also : " ,$Famtree{$PFbkey}{name}, ") \$nodes_n=", $nodes_n, " \n";
			#~ if ($PFbtree{$PFbkey}{count} ==0){$PFbtree{$PFbkey}{count} =1;}
			#~ printf O "gainFam: %d lossFam: %d N_Fam: %d PF_P: %f %s\n",$PFbtree{$PFbkey}{gainF},$PFbtree{$PFbkey}{lossF},$FamtreePFb{$PFbkey}{Nfam},$PFbtree{$PFbkey}{pf_P}/$PFbtree{$PFbkey}{count},$FamtreePFb{$PFbkey}{name};
		#~ }
						
		#~ printf("\nThe Summary File is named: %s\n\n", $outfilePFb);
		#~ close(O);
		
	#~ }
	

# --- MY subroutines -----------------------------------------------------------
		
sub traceIt($,$,%) 
		{
			
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
