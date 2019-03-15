#!/bin/bash

clear;

PATHEPOPE="/homes/biertruck/gabor/""epope_2_paper/epope_2/"		#DIR
ALIGN="/homes/biertruck/gabor/epope_2_paper/epope_bash_run/mipf_temp/"		#DIR
RESULTS="/homes/biertruck/gabor/epope_2_paper/ergebnisse/script_test_2/"		#DIR
TREE="/homes/biertruck/gabor/epope_2_paper/epope_bash_run/collapsed_tree.dat"	#FILE

echo -e "path for epope binary is:\t $PATHEPOPE "
echo -e "path for alignments is:\t\t $ALIGN "
echo -e "path for results is:\t\t $RESULTS"
echo -e "path to the tree:\t\t $TREE \n\n"

datei="$PATHEPOPE""ePoPE"
#echo $datei

if [ -f "$datei" ]; then
        echo -e "found ePoPE:\t\t\t" "$datei"
        else 	printf "\nePoPE not found, exiting now\n"
				exit
				
fi


if [ -d "$ALIGN" ] 
	then 
		echo -e "found alignmentpath:\t\t $ALIGN "
	else  printf "\nalignmentpath not found, exiting now\n"
			exit
		
fi



if [ -d "$RESULTS" ] 
	then 
		echo -e "found results path:\t\t $RESULTS  "
	else	printf "\nresults path not found, exiting now\n"
			exit
			
fi


if [ -f "$TREE" ] 
	then 
		echo -e "found tree file:\t\t $TREE  "
	else 
		printf "\ntree file not found, exiting now\n"
		exit
		
fi

read -p "Press [Enter] key to start calculations..."



#for i in "$ALIGN"* ; do echo "$datei" "$i" "$TREE" "$RESULTS""${i##*/}" "$RESULTS""${i##*/}" ; done


# doing the main epope calculation, parsimony and pf: 

#for i in "$ALIGN"*; do "$datei" -i "$i" -t "$TREE" -o "$RESULTS""${i##*/}" -p "$RESULTS""${i##*/}" -b; done


for i in "$ALIGN"* ; do "$datei" -i "$i" -t "$TREE" -o "$RESULTS""${i##*/}" -p "$RESULTS""${i##*/}" -z -P  ; done

for i in "$ALIGN"* ; do "$datei" -i "$i" -t "$TREE" -o "$RESULTS""${i##*/}" -p "$RESULTS""${i##*/}" -b ; done

#for i in "$ALIGN"* ; do "$datei" -i "$i" -t "$TREE" -o "$RESULTS""${i##*/}" -p "$RESULTS""${i##*/}" ; done

#checking whether the summarize script is located in the ePoPE folder as it should be:

if [ -f "$PATHEPOPE""ePoPE.summarize.pl" ]; then
       echo -e "\nfound the summarize perl script: " "$PATHEPOPE""ePoPE.summarize.pl"
      else printf "\nno summarize script found, it's supposed to be in the ePoPE source folder, exiting now\n"
      
fi

# ordne ranlegen, in die ePoPE.summarize.pl die ergebnisse legen soll
if [ !  -d "$RESULTS""PF" ]; then 
		echo -e "\ncreating a $RESULTS""PF directory\n"
		mkdir "$RESULTS""PF"
fi



if [ !  -d "$RESULTS""PFb" ]; then 

		 echo -e "\ncreating a $RESULTS""PFb directory\n"
		mkdir "$RESULTS""PFb"
fi


if [ !  -d "$RESULTS""PSb" ]; then 
		echo -e "\ncreating a $RESULTS""PSb directory\n"
		mkdir "$RESULTS""PSb"
fi






"$PATHEPOPE""ePoPE.summarize.pl" -d "$RESULTS" -o "$RESULTS"

#alle SUM_<...> dateien sind jetzt fertig berechnet, abfragen welche da 
#sind, diese dann in die jeweiligen ordner schieben:




#trace gene fam ist schon on ePoPE.summarize.pl mit eingearbeitet

Rscript ~/epope_2_paper/epope_2/Rscript.R "$RESULTS""PF.out_lst_4R.dat" "$RESULTS""PFb.out_lst_4R.dat" "$RESULTS""PSb.out_lst_4R.dat"

mv "$RESULTS"PFb*dat "$RESULTS"PFb/
mv "$RESULTS"PF*dat "$RESULTS"PF/
mv "$RESULTS"PSb*dat "$RESULTS"PSb/
rm $RESULTS*.out
rm $RESULTS*.ps
