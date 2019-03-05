#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "read.h"
#include "calc.h"

#include <sys/stat.h>
#include <unistd.h>
//#include <dirent.h>
#include <error.h>


extern struct gl_nodeN treeN[];
extern int nodesN;

int main(int argc, char *argv[]) 
{
	struct gl_arguments arg;
	
		
	/*optionen und input bzw output namen werden gecheckt und gesetzt*/
	arg = gl_readArguments(argc, argv);  
	arg = getFilenameExtension(arg);
	
	
	
	if(arg.collectfile != NULL)
		{
			
	
			gl_buildTree(arg.treefile);		
			
			
			
			gl_readSummaryTreeN(arg.collectfile); 
			int lca, kMin, kMax;	
			lca = gl_getLCA(&kMin, &kMax); // label of root of subtree 
			//printf("least common ancestor: kmin: %d    kmax: %d\n\n", kMin,kMax);
			//float **pku=NULL;
			gl_printTreePSN(kMax, lca, arg);
			//gl_printTreeN(arg); eigentlich überflüssig, all das steh ja schon in der summary!!!
			gl_freeArguments(arg);
			return(0);
			
		}
	
	//
			
	
			
	gl_printArguments(arg);
	
	gl_readDATA(arg);  
	
	// baum ist hier soweit fertig gebaut, jetzt noch die weights einfügen:
	//~ if(arg.scores !=NULL)
	//~ {	
		//~ set_weights(arg);	
	//~ }
	
	
	
	
	
		if(arg.collectfile==NULL)
		{
			if (arg.b==1 || arg.z==0) //das ist ok, wenn z==0 und b==0 scores, dann beenden in calc vor PF
			{	
				gl_calc(arg);
			}	
		
				
			if (arg.z ==1) //
			{
				pf_with_z(arg);
			}	
		}
		
		gl_freeArguments(arg);
		
		return 0;
	
}
