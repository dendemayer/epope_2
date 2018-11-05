#include "read.h"
#include "newicktolist2testbaumok.h"

#define NEWICKD

/* weight to gain a miRNA for each node */
float treeGainWeight[N_NODES];			//N_NODES in tree.h definiert (400)

/* weight to loose a miRNA for each node */
float treeLossWeight[N_NODES];

/* tree with max. N_CHILDREN number of child nodes */
struct gl_nodeN treeN[N_NODES];

struct collapse treeC[N_NODES];
/* preorder of treeN */
int treeOrderN[N_NODES];

/* position of nodes in treeN[] according to preorder */
int treePosN[N_NODES];

/* number of used nodes in built tree */
int nodesN;


/* -------------------------------------------------------------------------- */
/* gl_readArguments fragt epope aufruf ab und belegt mit funktion
 * gl_initArguments die inhalte der struct gl_Arguments*/
struct gl_arguments gl_readArguments(int argc, char *argv[]) {

  int i;
  gl_arguments args;    // args -> variabĺe vom typ gl_arguments wird initialisiert 

  if(argc < 2) usage(); //wenn nur ePoPE aufgerufen wird, usage()
  /* printf("argc: %d, argv[0]: %s\n",argc,argv[0]); */
  args = gl_initArguments();   //args werden initialisiert

  for(i=0; i<argc; i++) 
{
	if(argv[i][0]=='-') 
	{
		switch(argv[i][1]) 
		{
			case 'h':
			usage();
			break;
			case 'v':
			version();
			break;
			case 'i':
	/* get alignment file */
			if(argc < i+2) usage();
			args.infile = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));         // für die länge des strings argv[i+1] wird platz reserviert in args.infile
			strcpy(args.infile,argv[i+1]);
			break;
			case 't':
	/* get tree */
			if(argc < i+2) usage();
			args.treefile = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));
			strcpy(args.treefile, argv[i+1]);
			break;
			
			case 'd':
	/* get directory of alignment files */
			if(argc < i+2) usage();
			args.directory = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));
			strcpy(args.directory, argv[i+1]);
			break;
			
			
			case 'n':			//-n newickfile wird eingelesen
			if(argc < i+2) usage();
			args.treefileNew = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));
			strcpy(args.treefileNew, argv[i+1]);
			break;
						
			
					
			case 'z':			//-z Partition function is on, NO scores are calculated anymore
			if(argc < i+1) usage();
			args.z = 1;
			break;
		
			case 'b':			//-b Partition function AND scores are calculated
			if(argc < i+1) usage();
			args.b = 1;
			break;
						
			case 'P':			//-P k is set to the maximum calculatet p in PF modus
			if(argc < i+1) usage();
			args.P = 1;
			break;
			
			case 'C':			//-C do not collapse the tree
			if(argc < i+1) usage();
			args.C = 1;
			break;
			
			
			case 'w':
	/* get weight array */
			if(argc < i+2) usage();
			args.scores = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));
			strcpy(args.scores, argv[i+1]);
			break;
			case 'o':
	/* get dat-outfile name */
			if(argc < i+2) usage();
			args.outfile = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char)); 
			strcpy(args.outfile, argv[i+1]);
			break;
			case 'p':
	/* get ps-outfile name */
			if(argc < i+2) usage();
			args.psfile = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));
			strcpy(args.psfile, argv[i+1]);
			break;
			case 'c':
	/* get a collection of ePoPE output created with ePoPE.summarize.pl and draw tree */
			if(argc < i+2) usage();
			args.collectfile = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));
			strcpy(args.collectfile, argv[i+1]);
			break;
			case 'l':
	/* activate loss-backtracing */
			args.lossBT = 1;
			break;
			case '-':
			if(strcmp(argv[i],      "--help"     ) ==0) usage();
			else if(strcmp(argv[i], "--version"  ) ==0) version();
			else if(strcmp(argv[i], "--clw"      ) ==0) args.clw = 1;
			else if(strcmp(argv[i], "--stk"      ) ==0) args.stk = 1;
			else if(strcmp(argv[i], "--type"     ) ==0) 
			{
				if(argc < i+2) usage();
				args.type = (char *) calloc(strlen(argv[i+1]) + 1, sizeof(char));
				strcpy(args.type, argv[i+1]);
			}
			else usage();
			break;
			default:
			usage();
		}
    }
    
}

  if(args.infile == NULL) {
    if(args.collectfile == NULL) {
      gl_freeArguments(args);
      usage();
    }
  }

  /* if infile is set, collectfile is ignored ! */
  if(args.infile != NULL && args.collectfile !=NULL) {
    free(args.collectfile);                                 
    args.collectfile = NULL;
  }
  if(args.treefile == NULL && args.treefileNew == NULL) {
    gl_freeArguments(args);                                
    usage();
  }
 

  if(args.clw == 1 && args.stk == 1) { /* ?? -> check type of aln */
    args.clw = 0; args.stk = 0;
  }

  return args;              
}

/* init functions ----------------------------------------------------------- */

struct gl_arguments gl_initArguments() {

  gl_arguments a;

  a.infile = NULL;
  a.treefile = NULL;
  a.treefileNew = NULL;
  a.outfile = NULL;
  a.outfileFlag = NULL;
  a.outfileFlag_b = NULL;
  a.psfile = NULL;
  a.psfileFlag = NULL;
  a.scores = NULL;
  a.collectfile = NULL;
  a.type = NULL;
  a.directory=NULL;
  a.clw = 0;
  a.stk = 0;
  a.z = 0;
  a.b = 0;
  a.P = 0;
  a.C = 0;
  a.lossBT = 0;
  a.collaps = 0;
  a.seenBWp =0;

  return a;

}

/* print functions ---------------------------------------------------------- */

struct gl_arguments gl_printArguments(gl_arguments ga) {
	
	
	if( ga.z == 1 && ga.b==1)
	{
		printf("please choose o ne of the options -b OR -z\nexiting programm\n");
		exit(0);
	}
	
	if( ga.directory != NULL && ga.outfile!= NULL)
	{
		printf("when -d is set, you don't have to set -o \nexiting programm");
		exit(0);
	}
	
	if( ga.directory != NULL && ga.psfile!= NULL)
	{
		printf("when -d is set, you don't have to set -p \nexiting programm");
		exit(0);
	}
	
	//ga =getFilenameExtension(ga);
	
	//printf("\n ga.psfileFlag  %s\n",ga.psfileFlag );
	//printf("\n ga.outfileFlag %s\n",ga.outfileFlag );
	//printf("\n ga.outfileFlagb %s\n\n",ga.outfileFlag_b );
	
	
	

	
  fprintf(stdout, "\nArguments:\n-----------\n");                    //
  fprintf(stdout, "%-35s %s\n","Input aln file:" ,ga.infile);          //exampl.stk
 
  if (ga.treefileNew!= NULL)
  fprintf(stdout, "%-35s %s\n","Input NEWICK tree file:", ga.treefileNew);		//
  else
  fprintf(stdout, "%-35s %s\n","Input tree file:", ga.treefile);        //example.tree.dat
  
	if(ga.outfile != NULL) 
	{
		if(ga.b==0)
		{
			fprintf(stdout, "%-35s %s\n","Output file:", ga.outfileFlag);         //example.ePoPE.out
		}
		else if(ga.b==1)
		{
			fprintf(stdout, "%-35s %s\n","Output file:", ga.outfileFlag);
			fprintf(stdout, "%-35s %s\n","2. Output file:", ga.outfileFlag_b);
		}	
	}	
	else fprintf(stdout,  "%-35s no output file applied, writing to standard out\n", "Output file:");
  
  if(ga.psfile != NULL) fprintf(stdout, "%-35s %s\n","PS-output file:", ga.psfileFlag);           //example.ePoPE.out.ps
  else fprintf(stdout, "%-35s %s\n","PS-output file:", ga.psfileFlag);
  
  fprintf(stdout, "%-35s %s\n","Score array file:", ga.scores);     //(null)
  
  if(ga.type != NULL) fprintf(stdout, "%-35s %s\n","Type:", ga.type); //all
  else fprintf(stdout, "%-35s all\n", "Type:");
  
  fprintf(stdout, "%-35s %d\n","-P is set to:", ga.P);
  fprintf(stdout, "%-35s %d\n","-C is set to:", ga.C);
  fprintf(stdout, "%-35s %d\n","-z is set to:", ga.z);
  fprintf(stdout, "%-35s %d\n","-b is set to:", ga.b);
  fprintf(stdout, "%-35s %s\n","-c is set to:", ga.collectfile);
 // fprintf(stdout, "%-35s %s\n","-d is set to:", ga.directory);
  
    
  if(ga.b ==1)
  fprintf(stdout, "\nCalculating Partition Funktion and Parsimony Scores\n\n");
  
  if(ga.z ==1 && ga.b==0)
  fprintf(stdout, "\nCalculating Partition Function only\n\n");
  
  if(ga.b==0 && ga.z ==0)
  fprintf(stdout, "\nCalculating Parsimony Scores only\n\n");
  
  return ga;
 //~ free(ga.outfileFlag);
 //~ free(ga.outfileFlag_b);
 //~ free(ga.psfileFlag);
  
  
}

/* -------------------------------------------------------------------------- */

void usage() {

  fprintf(stdout, "\n+--------------------------------------------------+");
  fprintf(stdout, "\n| ePoPE 2.0                                        |");
  fprintf(stdout, "\n|                                                  |");
  fprintf(stdout, "\n| ePoPE - efficent Prediction of Paralog Evolution |");
  fprintf(stdout, "\n+==================================================+\n");
  fprintf(stdout, "\nePoPE predicts a maximal parsimony solution of gain and loss events of a gene family with paralogs.\n");

  fprintf(stdout,"\nUsage: ePoPE [ arguments ] -i ALNFILE -t TREEFILE (or -n NEWICK-TREEFILE)\n\n");
  fprintf(stdout, "arguments: [-o OUTFILE] [-p PS-OUTFILE]\n");
  fprintf(stdout, "           [-c COLLECTFILE]\n");
  fprintf(stdout, "           [-h,--help] [-v,--version]\n");
  fprintf(stdout, "           [--type TYPE]\n");
  fprintf(stdout, "           [-z] [-b] [-C] [-P]\n\n");
  
  

  fprintf(stdout, "%-20s Input alignment FILE in CLUSTALW/STOCKHOLM format. [REQUIRED]\n", "-i FILE");
  fprintf(stdout, "%-20s Input tree FILE see example.tree.dat format. [REQUIRED]\n", "-t FILE");
  fprintf(stdout, "%-20s Input tree FILE in newick format see example.newick.tree.dat format. [REQUIRED]\n", "-n FILE");
  fprintf(stdout, "%-20s Input weight array FILE. [OPTIONAL]\n", "-w FILE");
  fprintf(stdout, "%-20s Output FILE for tree data. Default is writing to stdout.[OPTIONAL]\n", "-o FILE");
  fprintf(stdout, "%-20s Output FILE for PS-tree data. Default is 'INFILE.ps'. [OPTIONAL]\n", "-p FILE");
  fprintf(stdout, "%-20s FILE is a collection of calls to ePoPE with the same tree on a set of gene families created via 'ePoPE.summarize.pl'. This option forces ePoPE to draw this summarized tree. You must provide the tree file you used for the single ePoPE calls with -t option. Example call: ePoPE -c COLLECTFILE -t TREEFILE -p PS-OUTFILE [OPTIONAL]\n", "-c FILE");
  fprintf(stdout, "%-20s TYPE is one of {genes, gainFam, lossFam, gain, loss, all}. Is the type of values that are plotted in the tree. Default: 'all'. [OPTIONAL]\n", "--type TYPE");
  fprintf(stdout, "%-20s Calculating partition function.\n", "-z");
  fprintf(stdout, "%-20s Calculating both, partition function and parsimony scores.\n", "-b");
  fprintf(stdout, "%-20s suppress the collapse of inner nodes with a degree of two.\n", "-C");
  fprintf(stdout, "%-20s setting the number of genes to the highest probability value in back-recursion.\n", "-P");
  fprintf(stdout, "%-20s Show this help message.\n", "-h,--help");
  fprintf(stdout, "%-20s Show version information.\n", "-v,--version");
  fprintf(stdout, "\nExample call:\n\n");
  fprintf(stdout, "\t./ePoPE -n example_newick -i example_ali1 -o example_ali1 -p example_ali1 -b\n\n");
  fprintf(stdout, "\nPlease feel free to contact me for comments, bug-reports, etc.\n\n");

  fprintf(stdout, "--\n");
  version();

}

/* -------------------------------------------------------------------------- */

void version() {

  fprintf(stdout, "ePoPE 2.0\n\n");
  fprintf(stdout, "Auhthor: Jana Hertel and Gabor Balogh:\n\n");
  fprintf(stdout, "         jana@bioinf.uni-leipzig.de\n");
  fprintf(stdout, "         gabor@bioinf.uni-leipzig.de\n\n");
  fprintf(stdout, "Date:    December, 2016\n\n");

  exit(EXIT_SUCCESS);

}

/* free functions ----------------------------------------------------------- */

void gl_freeArguments(gl_arguments ga) {

  if(ga.infile != NULL)   free(ga.infile);
  if(ga.treefile != NULL) free(ga.treefile);
  if(ga.treefileNew != NULL) free(ga.treefileNew);
  if(ga.scores != NULL)   free(ga.scores);
  if(ga.outfile != NULL)  free(ga.outfile);
  if(ga.psfile != NULL)   free(ga.psfile);
  if(ga.collectfile != NULL)   free(ga.collectfile);
  if(ga.type != NULL)   free(ga.type);
  if(ga.psfileFlag != NULL)   free(ga.psfileFlag);
  if(ga.outfileFlag != NULL)   free(ga.outfileFlag);
  if(ga.outfileFlag_b != NULL)   free(ga.outfileFlag_b);
  if(ga.directory != NULL)   free(ga.directory);
  
}

void collapstest( char *filename, struct gl_arguments *ga)	//die funktion muss entweder neues struct gl_arguments zurückgeben, mit dem weitergearbeitet wird, oder ihr muss der pointer auf die struct übergeben werden, damit änderungen auch in der aufrufenden funktion wirksam werden (call by reference, value))
{
	//printf("vor änderung: in collapstest angekommen, der flag  ist : %d\n", ga->collaps);
	//ga->collaps=1;
	//printf("nach änderung: in collapstest angekommen, der flag  ist : %d\n", ga->collaps);
	//alles was man wissen muss, ist ob es einen knoten gibt, der nur ein kind hat, dann kann test abgebrochen werden, gl_argumentsflag auf 1 setzen und zurück
	//return ga;
	
	
	FILE *fp=NULL;
	char line[512];
	int i,n,n2, par=-1, read=0;
	int nodes=0,lines=0;
	char ch, **splitLine, **kids;
	int checkiftocollapse =0;




	fp = fopen(filename, "r");
	while(!feof(fp)) 
	{				//Zeilen werden gezählt
		ch = fgetc(fp);	//mit fgetc(FILE *fp) wird pointer automatisch eine position vorgerückt
		if(ch == '\n') lines++;
	}
	
	if(lines > N_NODES) 
	{
		fprintf(stderr, "ERROR: you have more than '%d' nodes in your treefile. I cannot handle this!\nSorry.\n",N_NODES);
		fclose(fp);
		exit(0);
	}
	
	//fclose(fp);	
	rewind(fp);   //The C library function void rewind(FILE *stream) 
				//	sets the file position to the beginning of the 
				//	 file of the given stream.
					
					//--> nach dem zeilenzählen, wieder zum anfang der datei
	
	// while(fgets(line, 512, fp)) { 
	
	int labtoc[lines];		//labtocollapse
	int postoc =0;
	
	while(fgets(line, 512, fp)) 
	{
		splitLine = gl_splitString(line, ' ',&n); //zwischen den leerzeichen wird string geiteilt, und z.B. splitline[0] entspricht species,
		strcpy(treeC[nodes].o,splitLine[0]); 
		read = sscanf(splitLine[2], "%d", &treeC[nodes].n);    //label wird gesetzt
		assert(read);	//
		read = sscanf(splitLine[4], "%d", &par);    //par wird gesetzt
		assert(read);
		if(par == -1) 
		{
			treeC[nodes].parentLab = -1;
			treeC[nodes].p = NULL;
		} 
		else 
		{ 
				treeC[nodes].parentLab = par;
				//printf("nodes: %d,  par: %d\n", nodes,par);
		}
       	read = sscanf(splitLine[8], "%d", &treeC[nodes].level); assert(read);
		read = sscanf(splitLine[10], "%f", &treeC[nodes].gainW);  assert(read);
		read = sscanf(splitLine[12], "%f", &treeC[nodes].lossW);  assert(read);
	
		kids = gl_splitString(splitLine[6], ',', &n2); //wenn in diesem string kein komma, gibts nur ein kind (außer bei leaves, da -1 bei kind)
		
		
		
		if (n2==1 && strcmp("-1",kids[0])!=0) //wenn anzahl kinder ==1 und kind ist kein leave, dann knoten gefunden, welcher collapsed werden muss
		{
			checkiftocollapse =1;
			//printf(" zu collapsen lab: %d\n", treeC[nodes].n);
			labtoc[postoc]= treeC[nodes].n;
			treeC[nodes].tocallapse=1;
			//printf("postoc:%d\nund lab:%d\n", postoc,labtoc[postoc]);
			postoc++;
		}
		else
		treeC[nodes].tocallapse=0;
					
		for(i=0; i<n2; i++) 
		{
		read = sscanf(kids[i], "%d", &treeC[nodes].c[i]); assert(read);
				//printf("kids[%d]=%s\n",i,kids[i]); 
		}
	
		//treeN[nodes].m = 0;
		//treeN[nodes].P = 0;
		if(n2 == 1 && treeC[nodes].c[0] == -1) 
		{
		treeC[nodes].nc = 0;
		// treeN[nodes].seenBW = 1; 
		} else 
		{
		treeC[nodes].nc = n2;
		// treeN[nodes].seenBW = -1; 
		}
		// liste von childnodes etc
		//printf("tree[%d].n=%3d\tKIDS: %d tree[%d].o=%s\n",nodes,treeN[nodes].n,treeN[nodes].nc,nodes,treeN[nodes].o); 
		//	for(i=0; i<n2; i++) { printf("treeN[%d].c[%d] = %d\n",nodes,i,treeN[nodes].c[i]);    } 
		//	fflush(stdout); 
		nodes++;
		for(i=0; i<n; i++) { free(splitLine[i]); }
		free(splitLine);
		for(i=0; i<n2; i++) { free(kids[i]); }
		free(kids);
	}
	fclose(fp);	
	
	
	
	//wenn kein zu kollapsender knoten gefunden wird, kann hier return, sonst wird ga.treefile geändert in output "collapsed_tree.dat"
	if (checkiftocollapse==1)
	{
		ga->treefile="collapsed_tree.dat";
		ga->collaps=1;
		//fprintf(stdout, "The tree got collapsed, new nodelist in 'collapsed_tree.dat' and the deleted nodes in 'collapsed_nodes.log'\nPlease note, if ePoPE in collect-mode shall be used, apply the 'collapsed_tree.dat' file with the '-t' option \n\n");
	}
	else
	return;
	
//	for (int i=0; i<postoc; i++)
//	printf("labtoc[%d]: %d\n", i, labtoc[i]);
	
/*	for (int i=0; i<lines; i++)
	{
		printf("\n\nnode in line: %d\n\n",i);
		printf("species names: %s\n", treeC[i].o);
		printf("label n: %d\n", treeC[i].n);
		printf("parentlab : %d\n", treeC[i].parentLab);
		printf("1. kind: c[0]: %d\n", treeC[i].c[0]);
		printf("kinder insgesamt nc: %d\n", treeC[i].nc);
		printf("level: %d\n", treeC[i].level);
		printf("gainW: %f\n", treeC[i].gainW);
		printf("lossW: %f\n", treeC[i].lossW);
	}
*/


		//Knoten werden "umgehangen" 
	for (int l=0; l<postoc; l++)
	{
		for (int i=0; i<lines; i++)
		{
			if (labtoc[l] == treeC[i].n)
			{	//printf("zu collapsen: %s\n", treeC[i].o);
				//ab hier müssen noch mal alle child arrays aller knoten nach dem lab des zu löschenden knoten abgesucht werden, dieser muss dann ausgetauscht werden mit dem child lab das labtoc
				for (int j=0; j<lines; j++)
				{
					for (int c=0; c<treeC[j].nc; c++)
					{	
						if (treeC[j].c[c]== treeC[i].n) //suche labtoc in allen kind listen
						{
							treeC[j].c[c]= treeC[i].c[0];	//wenn in einem kind array labtoc auftaucht, wird dieser ersetzt mit kind des labtoc 
						}
					}
					if (treeC[j].parentLab==treeC[i].n)	//suche labtoc in allen parentLab
					{
						treeC[j].parentLab=treeC[i].parentLab;	
					}
				}
			}
		} 
	}
	//level setzen:
	
		
	//struct collapse treeC[] ist fertig gefüllt bis lines, und man hat alle zu collapsende labs herausgeschrieben in labtoc[]
	//level setzen, zählen, solange ein knoten einen vater hat
	
	//int node = 0;	//wird zeile für zeile erhöht für suche des levels des knoten der jeweiligen zeile
	//int foundroot =0;
	int templabparent = 0;	//setzen der vaäter, wenn -1, wird level gesetzt
	int level=0;	//hochzählen des levels
	
	
	for (int i=0; i<lines; i++)
	//printf("%s lab %d par %d  d2root %d  tocollapse  %d\n", treeC[i].o, treeC[i].n, treeC[i].parentLab,  treeC[i].level, treeC[i].tocallapse);
	
	
	for (int i= 0; i<lines; i++)// jede zeile um level zu setzen
	{	
		if (treeC[i].tocallapse==1)			//wenn knoten collapsed wird, überspringen
		continue;
		templabparent = treeC[i].parentLab;
		
		//newsearch:
		for (int j=0; j<lines;j++)// suche des vaters eines labs in allen zeilen
		{
			
			if (treeC[j].tocallapse==1)		//wenn knoten collapsed wird, einfach überspringen
			continue;
			
			else if (templabparent == -1)	//bei wurzel angekommen, level wird gesetzt
			{
				treeC[i].level=level;
				//printf("\nlevel von %s ist %d\n",treeC[i].o,treeC[i].level);
				//printf("lines:%d\n",lines);
			//	printf("i=%d\n",i);
			//	printf("j=%d\n",j);
				level=0;
				break; //schleife kann abgevrochen werden, zum nächsten i
			}
			else if (treeC[j].n==templabparent) //vater ist nicht wurzel, level wird erhöht, templabparent wird neu gesetzt und erneute suche nach neuen vater
			{	//for (int k=0; k<postoc;k++)
				//if()
			//	printf("vater von lab %d ist %d\n",templabparent, treeC[j].parentLab);
				templabparent=treeC[j].parentLab;
				level++;
				j=-1; 		//schleife darf hier nicht abgebrochen werden, da sich i nicht erhöhen darf!!!!sie soll sich wiederholen!!!
				continue;
				//goto newsearch; //schleife darf hier nicht abgebrochen werden, da sich i nicht erhöhen darf!!!!sie soll sich wiederholen!!!
			
			}
			
		}
	}
/*	
	for (int i=0; i<postoc; i++)
	{
		printf("labtoc[%d]=%d\n",i,labtoc[i] );
	}
*/

/*	for (int i=0; i<lines; i++)
	{
		printf("\n\nnode in line: %d\n\n",i);
		printf("species names: %s\n", treeC[i].o);
		printf("label n: %d\n", treeC[i].n);
		printf("parentlab : %d\n", treeC[i].parentLab);
		printf("1. kind: c[0]: %d\n", treeC[i].c[0]);
		//printf("kinder insgesamt nc: %d\n", treeC[i].nc);
		printf("level: %d\n", treeC[i].level);
		printf("tocpllapse: %d\n", treeC[i].tocallapse);
		//printf("gainW: %f\n", treeC[i].gainW);
		//printf("lossW: %f\n", treeC[i].lossW);
	}
*/		
	FILE *fp2=NULL;
	fp2 = fopen("collapsed_tree.dat", "w");
						
	for (int i=0; i<lines; i++)
	{
		if (treeC[i].tocallapse==0)
		{		
			fprintf(fp2,"%s lab %d par %d children ", treeC[i].o, treeC[i].n, treeC[i].parentLab);
			if (treeC[i].c[0]==-1)
				{
					//printf("das ist der fall");
					fprintf(fp2,"-1 ");
				}
			for (int j=0; j<treeC[i].nc;j++)
			{
				if(j==treeC[i].nc-1)
				fprintf(fp2,"%d ",treeC[i].c[j]);
				else 
				fprintf(fp2,"%d,",treeC[i].c[j]);
				
			}
			fprintf(fp2,"d2root %d gainW %.0f lossW %.0f\n", treeC[i].level, treeC[i].gainW, treeC[i].lossW); 
			
			//fprintf(fp2,"1. kind: c[0]: %d\n", treeC[i].c[0]);
			//printf("kinder insgesamt nc: %d\n", treeC[i].nc);
			
			//printf("gainW: %f\n", treeC[i].gainW);
			//printf("lossW: %f\n", treeC[i].lossW);
		
		}
	}
	fclose(fp2);
	
	FILE *fp3=NULL;
	fp3 = fopen("collapsed_nodes.log", "w");
	
	for (int i=0; i<lines; i++)
	{	
		if (treeC[i].tocallapse==1)
		{	
			fprintf(fp3,"%s lab %d par %d children ", treeC[i].o, treeC[i].n, treeC[i].parentLab);
			if (treeC[i].c[0]==-1)
				{
					//printf("das ist der fall");
					fprintf(fp3,"-1 ");
				}
			for (int j=0; j<treeC[i].nc;j++)
			{
				if(j==treeC[i].nc-1)
				fprintf(fp3,"%d ",treeC[i].c[j]);
				else 
				fprintf(fp3,"%d,",treeC[i].c[j]);
			}
			fprintf(fp3,"d2root %d gainW %.0f lossW %.0f\n", treeC[i].level, treeC[i].gainW, treeC[i].lossW); 
			
		}
	}
	fclose(fp3);

	
}

void gl_readDATA(struct gl_arguments ga)
{

	
		
	//wenn newick gesetzt, wird hier übersetzte knotenliste übergeben und der char pointer auf ga.treefile gesetzt
	

	if (ga.treefileNew != NULL)
		{
			newick(ga.treefileNew);
			
			//printf("kommen wir hier an?\n\n\n");
			
			ga.treefile="newick_to_nodelist.dat";
		}
	
// nachdem newick eingelesen wurde, oder auch nicht, wird knotenliste überprüft, ob sie callapst werden muss:
	//printf("in r eadaDATA angekommen, der flag  ist : %d\n", ga.collaps);
	//printf("wie groß ist die struct gl_argument s?:%lu\n", sizeof(ga)); //96 byte
	
	//kann  bei -C option übersprungen werden
	
	
	if (ga.C == 0)
	{
		collapstest(ga.treefile, &ga);
		
	}
		printf("\nThe nodelist tree used is '%s', please use this file, if you want to print a summarized tree with the ePoPE collect option '-c'\n\n",ga.treefile);
		if (ga.z==0 && ga.b==0 && ga.P==1)
	{
		printf("Please be aware, that using the '-P' option in default parsimony mode does not have any effect. The calculation will start nevertheless.\n\n");
	}
	//printf("2. mal in readaDATA angekommen, der flag  ist : %d\n", ga.collaps);
		
		gl_buildTree(ga.treefile);
	
	int i=0, read;
	FILE * fp;        //Stream welcher mit fgets ausgelesen wird fp ist der FILE Pointer auf den Stream 
	char line[2000],* name, abbr[4];
	float g=-1, l=-1;		//
 /* read tree, build datastructure gl_node tree[] */
 /* read alignment */
	fp = fopen(ga.infile, "r");       //  "r"	Opens a file for reading. The file must exist.

	if(fp == NULL) 
	{
		fprintf(stderr, "Cannot access input file: '%s'.\n", ga.infile);
		gl_freeArguments(ga);
		exit(EXIT_FAILURE);
	}



  /* read file linewise */
	char *pch;

	while(fgets(line, 2000, fp))
	{        //char *fgets(char *str, int n, FILE *stream)
                                        // str -- This is the pointer to an array of chars where the string read is stored.
                                        //n -- This is the maximum number of characters to be read (including the final null-character). Usually, the length of the array passed as str is used.
                                        //stream -- This is the pointer to a FILE object that identifies the stream where characters are read from.

    /* extract first word */

		pch = strchr(line, ' ');        //char *strchr(char *string, int c);
		if(pch == NULL) continue;
		name = (char *) calloc(pch-line+2, sizeof(char));
		strncpy(name, line, pch-line+1);
    //    read = sscanf(line, "%s %s", name, seq);
    //if(read == -1) continue;

		if(strcmp(name, "STOCKHOLM")     ==0) { free(name); continue;}
		else if(strcmp(name, "CLUSTAL")  ==0) { free(name); continue;}
		else if(strcmp(name, "CLUSTALW") ==0) { free(name); continue;}
		else if(strcmp(name, "//")       ==0) { free(name); continue;}
		else if(strlen(name) == 0) { free(name); continue;}
    /* exclude lines starting with '#' (comments) */
		if(name[0] == 35)
		{   
			free(name);
			continue;
		}
  //printf("line: %s\nname1: %s\n",line,name);

    /* extract species abbreviation */
		int namesplit_n = 0;
		
		
		char **namesplit = gl_splitString(name, ':', &namesplit_n);
		
		
    /* printf("namesplit_n = %d -> %s\n",namesplit_n, name); */
		if(namesplit_n == 1)
		{
      //      if(strncpy(abbr, name, 3)) {
			if(strncpy(abbr, name, 3))
			{
				abbr[3] = '\0';
	//printf("abbr: %s\n",abbr);
				for(i=0; i<nodesN; i++)
				{
	  /* printf ("%d %s %s\n", i,treeN[i].o,abbr); */
					if(strncmp(treeN[i].o, abbr, 3) == 0)
					{
						treeN[i].m++;
						treeN[i].P++;
						break;
					}	
	  /* if(strcmp(tree[i].o, abbr) == 0) { */
	  /*   tree[i].m++; */
	  /*   break; */
	  /* } */
				}
			}
		} 
		else
		{
			if(strncpy(abbr,namesplit[1],3))
			{
				abbr[3] = '\0';
				for(i=0; i<N_NODES; i++)
				{
	  /* printf ("%d %s %s\n", i,tree[i].o,abbr); */
					if(strcmp(treeN[i].o, abbr) == 0)
					{
						treeN[i].m++;
						treeN[i].P++;
						break;
					}
	  /* if(strcmp(tree[i].o, abbr) == 0) { */
	  /*   tree[i].m++; */
	  /*   break; */
	  /* } */
				}
			}

		}
    int z;
    for(z=0; z<namesplit_n; z++) 
    { 
		free(namesplit[z]); 
	}
    free(namesplit);
    free(name);

    abbr[0] = '\0';

	}

	fclose(fp);

	/* calculate P_v (max No of genes in subtree rooted at v) */
	
	
	gl_setPatInnerNodes();


  /* for(i=0; i<nodesN; i++) { */
  /*    printf("treeN[%d].m = %d\t",i,treeN[i].m); */
  /*    printf("treeN[%d].P = %d\t",i,treeN[i].P); */
  /*    printf("treeN[%d].o = %s\n",i,treeN[i].o); */
  /* } */
  /* exit(0); */

  /* read weights for gain and loss */
	if(ga.scores !=NULL)
	{
		fp = fopen(ga.scores, "r");
		while(fgets(line, 512, fp))
		{
			read = sscanf(line, "%d %f %f\n", &i, &g, &l);
			if(read == 0) continue;

			treeGainWeight[i] = g;
			treeLossWeight[i] = l;
		}

		fclose(fp);
	} 
	else
    {
        for(i=0;i<nodesN;i++)
        {
            treeGainWeight[i] = 1.;
            treeLossWeight[i] = 1.;
        }
    }

  return;

}
/*---------------------------ENDE void gl_readDATA(struct gl_arguments ga)--------------*/
/* -------------------------------------------------------------------------- */


void gl_buildTree(char *filename)
{		

	FILE *fp=NULL;
	char line[512];
	int i,n,n2, par=-1, read=0;
	int nodes=0,lines=0;
	char ch, **splitLine, **kids;





	fp = fopen(filename, "r");
	while(!feof(fp)) 
	{				//Zeilen werden gezählt
		ch = fgetc(fp);	//mit fgetc(FILE *fp) wird pointer automatisch eine position vorgerückt
		if(ch == '\n') lines++;
	}
	
	if(lines > N_NODES) 
	{
		fprintf(stderr, "ERROR: you have more than '%d' nodes in your treefile. I cannot handle this!\nSorry.\n",N_NODES);
		fclose(fp);
		exit(0);
	}
	
	rewind(fp);   /*The C library function void rewind(FILE *stream) 
					sets the file position to the beginning of the 
					* file of the given stream.*/
					
					//--> nach dem zeilenzählen, wieder zum anfang der datei
	
	/* while(fgets(line, 512, fp)) { */
	
	//printf("number of lines: %d \n",lines);
	
	while(fgets(line, 512, fp)) 
	{
		splitLine = gl_splitString(line, ' ',&n);
		strcpy(treeN[nodes].o,splitLine[0]);
		//printf("treeN[%d].o=%s\n",nodes,treeN[nodes].o);
		read = sscanf(splitLine[2], "%d", &treeN[nodes].n);    assert(read);	//
		read = sscanf(splitLine[4], "%d", &par);    assert(read);
		if(par == -1) 
		{
			treeN[nodes].parentLab = -1;
			treeN[nodes].p = NULL;
		} 
		else 
		{ 
				treeN[nodes].parentLab = par;
				//printf("nodes: %d,  par: %d\n", nodes,par);
		}
    
    	read = sscanf(splitLine[8], "%d", &treeN[nodes].level); assert(read);
		read = sscanf(splitLine[10], "%f", &treeN[nodes].gainW);  assert(read);
		read = sscanf(splitLine[12], "%f", &treeN[nodes].lossW);  assert(read);
	
		kids = gl_splitString(splitLine[6], ',', &n2);
		for(i=0; i<n2; i++) 
		{
		read = sscanf(kids[i], "%d", &treeN[nodes].c[i]); assert(read);
				//printf("kids[%d]=%s\n",i,kids[i]); 
		}
	
		treeN[nodes].m = 0;
		treeN[nodes].P = 0;
		if(n2 == 1 && treeN[nodes].c[0] == -1) 
		{
		treeN[nodes].nc = 0;
		/* treeN[nodes].seenBW = 1; */
		} else 
		{
		treeN[nodes].nc = n2;
		/* treeN[nodes].seenBW = -1; */
		}
		treeN[nodes].score = 0.;
		treeN[nodes].gain = 0;
		treeN[nodes].loss = 0;
		treeN[nodes].pInT = nodes;
		treeN[nodes].seenFW = -1;
		treeN[nodes].seenBW = -1;
				
		treeN[nodes].pfGain = 0;
		treeN[nodes].pfLoss = 0;
		treeN[nodes].pfm = 0;
		
		treeN[nodes].pfGainFam = 0;
		treeN[nodes].pfLossFam = 0 ;
		treeN[nodes].pfP = 0. ;
		treeN[nodes].psP = 0. ;
		
		/* liste von childnodes etc*/
		//printf("tree[%d].n=%3d\tKIDS: %d tree[%d].o=%s\n",nodes,treeN[nodes].n,treeN[nodes].nc,nodes,treeN[nodes].o); 
		//	for(i=0; i<n2; i++) { printf("treeN[%d].c[%d] = %d\n",nodes,i,treeN[nodes].c[i]);    } 
		//	fflush(stdout); 
		nodes++;
		for(i=0; i<n; i++) { free(splitLine[i]); }
		free(splitLine);
		for(i=0; i<n2; i++) { free(kids[i]); }
		free(kids);
	}
	fclose(fp);	
	nodesN = nodes;
  
  
	gl_label2pos();
  
  
	//  for(i=0;i<nodesN;i++) { printf("treePosN[%d]=%d\n",i,treePosN[i]); }
	/* printf("nodesN: %d\n",nodesN); */
	/* initialize treeOrderN */
	for(i=0; i<nodesN; i++) 
	{
		treeOrderN[i] = -1;
	}

	int root=0, p=0;
	for(i=0; i<nodesN; i++) 
	{
		if(treeN[i].parentLab == -1) 
		{
			root = i;
			break;
		}
	}

	p = gl_preorder(treeN[root], p);
	//  printf("nodesN: %d\nTree in pre-order:\n",nodesN);
	
	for(i=0;i<nodesN;i++) 
	{
		p = treePosN[treeOrderN[i]];
		if(treeN[p].parentLab >=0 ) 
		{
		int plab = treeN[p].parentLab;
		int pposintree = treePosN[plab];
		treeN[p].p = &treeN[pposintree];
		} else { treeN[p].p = NULL; }
		 //~ printf("treeOrderN[%3d] = %d -> pos in tree: %4d (%s) ",i,treeOrderN[i], p, treeN[p].o); 
		 //~ if(treeN[p].p) { printf(" parent: %d (%s)\n",(treeN[p].p)->n,(treeN[p].p)->o);} 
		 //~ else { printf("at root\n"); } 
	}
	
	

  //  exit(0);

}

/*----------------------------------------ENDE void gl_buildTree(char *filename)--------*/
/* -------------------------------------------------------------------------- */

/* fills array of treePosN with positions of labels */
/* index: label, value: position */
void gl_label2pos() {
  int i, j;

  for(i=0; i<nodesN; i++) {
    j = treeN[i].n;
    //    printf("T: treeN[%d].n=%-5d\t",i,treeN[i].n);
    treePosN[j] = i;
    //    printf("treePosN[%d]=%d\n",j,i);
  }

  return;
}

/* -------------------------------------------------------------------------- */

int gl_preorder(struct gl_nodeN node, int o) {

  int k, i;
  /* o is index in pre-order, value: label of node in treeN */
  treeOrderN[o] = node.n;
  treeN[treePosN[node.n]].pInP = o; /* stores position in preorder in corresp. node */

    //treePosN[node.n] = node.pInT;

  // printf("order[%d] = %d\t",o,treeOrderN[o]); 
   //printf("  pos[%d] = %d\t",treeOrderN[o],treePosN[treeOrderN[o]]); 
   //printf("%s -- %s\n", node.o, treeN[treePosN[node.n]].o); 
	
  o = o + 1;

  if(node.c[0] == -1) return o;/* leaf reached */

  for(k=0; k<node.nc; k++) {
    for(i=0; i<nodesN; i++) {
      if(treeN[i].n == node.c[k]) {
	o = gl_preorder(treeN[i], o);
	break;
      }
    }
  }

  return o;
}


char *strdup(const char *src) 
{
    size_t len = strlen(src) + 1;
    char *s = malloc(len);
    if (s == NULL)
        return NULL;
    return (char *)memcpy(s, src, len);
}


/* -------------------------------------------------------------------------- */

char** gl_splitString(char * str2split, const char a_delim, int *counter) {

  char ** result = 0;
  size_t count = 0;
  char * tmp = str2split;
  char * last_comma = 0;
  char delim[2];
  delim[0] = a_delim; delim[1] = 0;

      /* Count how many elements will be extracted. */
    while (*tmp)
    {
        if (a_delim == *tmp)
        {
            count++;
            last_comma = tmp;
        }
        tmp++;
    }

    /* Add space for trailing token. */
    count += last_comma < (str2split + strlen(str2split) - 1);

    /* Add space for terminating null string so caller
       knows where the list of returned strings ends. */
    count++;

    result = malloc(sizeof(char*) * count);

    if (result)
    {
        size_t idx  = 0;
        char* token = strtok(str2split, delim);

        while (token)
        {
            assert(idx < count);
            *(result + idx++) =  strdup(token);
            token = strtok(0, delim);
        }
        assert(idx == count - 1);
        *(result + idx) = 0;
	if(token != NULL) free(token);
    }

    * counter = count-1;

    return result;

}

/* -------------------------------------------------------------------------- */
   /* calculate P_v (max No of genes in subtree rooted at v) */
void gl_setPatInnerNodes() {

  int i, lab;

  for(i=0; i<nodesN; i++) {
    lab = treePosN[treeOrderN[i]];

    /* printf("%d: treeOrderN[%d]=",i,i); */
    /* printf("%d (pos %d)",treeOrderN[i],lab); */
    /* printf("P=%d\n",treeN[lab].P); */

    treeN[lab].P = gl_getPv(lab);
  }

  return;
}

int gl_getPv(int label_v) {

  int u, maxP, curP;
  /* printf("treeN[%d].nc = %d\n",label_v, treeN[label_v].nc); */
  if(treeN[label_v].nc == 0) {
    /* printf("\t\tLEAVE: P=%d (m=%d)\n",treeN[label_v].P,treeN[label_v].m); */
    return treeN[label_v].P;
  } else {
    maxP = 0;
    for(u=0; u<treeN[label_v].nc; u++) {
      /* printf("\tu: %d\n",u); */
      curP = gl_getPv(treePosN[treeN[label_v].c[u]]);
      if(curP > maxP) {
	maxP = curP;
      }
    }
    return maxP;
  }

}

/* -------------------------------------------------------------------------- */
/* reads the summarized tree and stores values in treeN */
void gl_readSummaryTreeN(char * file) {
  FILE *fp=NULL;
  char line[512];
  int j,i,n,n2, par=-1, read=0;
  int lines=0;
  char ch, **splitLine, **kids;


  fp = fopen(file, "r");
  while(!feof(fp)) {
    ch = fgetc(fp);
    if(ch == '\n') lines++;
  }

  if(lines > N_NODES) {
    fprintf(stderr, "ERROR: you have more than '%d' nodes in your treefile. I cannot handle this!\nSorry.\n",N_NODES);
    fclose(fp);
    exit(0);
  }

  rewind(fp);

  /* while(fgets(line, 512, fp)) { */
  j = 0;
  while(fgets(line, 512, fp)) 
{
    splitLine = gl_splitString(line, ' ',&n);
    read = sscanf(splitLine[1], "%d", &j); assert(read);
    read = sscanf(splitLine[3], "%d", &treeN[j].n); assert(read);
    read = sscanf(splitLine[7], "%d", &treeN[j].pInP); assert(read);
    read = sscanf(splitLine[11], "%d", &treeN[j].m); assert(read);

    read = sscanf(splitLine[13], "%d", &treeN[j].gain); assert(read);
    read = sscanf(splitLine[15], "%d", &treeN[j].loss); assert(read);
    read = sscanf(splitLine[17], "%d", &treeN[j].gainFam); assert(read);
    read = sscanf(splitLine[19], "%d", &treeN[j].lossFam); assert(read);

    /* parent */
    read = sscanf(splitLine[9], "%d", &par);    assert(read);
    if(par == -1) {
      treeN[j].parentLab = -1;
      treeN[j].p = NULL;
    } else {treeN[j].parentLab = par;
						//printf("gelesener vater: %d\n", treeN[j].parentLab);
						}
    /* name */
    //treeN[j].o = (char *) calloc(strlen(splitLine[22]) + 1, sizeof(char));
    //~ strcpy(treeN[j].o, splitLine[22]); //problem, \n wird mit eingefügt...
    //~ treeN[j].o[strlen(treeN[j].o)-1]='\0'; //austausch von \n mit \0
    
    //--->der spiciesname muss nicht nochmal gesetzt werden, der wurde schon in buildtree gesetzt....

    /* kids */
    kids = gl_splitString(splitLine[5], ',', &n2);
    for(i=0; i<n2; i++) {
      read = sscanf(kids[i], "%d", &treeN[j].c[i]); assert(read);
      /* printf("kids[%d]=%s\n",i,kids[i]); */
    }
    
    for(i=0; i<n2; i++) { free(kids[i]); }
	free(kids);
	for(i=0; i<n; i++) { free(splitLine[i]); }
	free(splitLine);
    j++;
    
  }

  fclose(fp);

  nodesN = j;
  gl_label2pos();
   // for(i=0;i<nodesN;i++) { printf("treePosN[%d]=%d\n",i,treePosN[i]); }
   //printf("nodesN: %d\n",nodesN); 
  /* initialize treeOrderN */
  for(i=0; i<nodesN; i++) {
    treeOrderN[i] = -1;
  }

  int root=0, p=0;
  for(i=0; i<nodesN; i++) {
    if(treeN[i].parentLab == -1) {
      root = i;
      break;
    }
  }

  p = gl_preorder(treeN[root], p);
   // printf("nodesN: %d\nTree in pre-order:\n",nodesN);

  for(i=0;i<nodesN;i++) 
  {
	p = treePosN[treeOrderN[i]];
	if(treeN[p].parentLab >=0 ) 
	{
		int plab = treeN[p].parentLab;
		int pposintree = treePosN[plab];
		treeN[p].p = &treeN[pposintree];
	} 
	else 
	{ 
		treeN[p].p = NULL; 
	}
		//printf("treeOrderN[%3d] = %d -> pos in tree: %4d (%s) the strlen is: %lu \n ",i,treeOrderN[i], p, treeN[p].o, strlen(treeN[p].o)); 
		if(treeN[p].p) 
		{ 
			//printf(" parent: %d (%s)\n",(treeN[p].p)->n,(treeN[p].p)->o);
		} 
		else 
		{ //printf("at root\n"); 
		} 
  }		
}


struct gl_arguments getFilenameExtension(gl_arguments ga)
 {
	//printf("z: %d\t b: %d\t filename: %s\n",ga.z, ga.b, ga.outfile);
	
	if(ga.z==1)   //
	{
		if(ga.outfile!=NULL)
		{
			int point =0;
			for(int i=0; i<strlen(ga.outfile); i++)
			{
				if(ga.outfile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			if(point==0 || ga.outfile[point-1]=='.' )
			{
				point=strlen(ga.outfile);
			}	
			//char newflag [100] ;  
			char *newflag  = (char*) calloc(100,sizeof(char));
			strncpy(newflag,ga.outfile,point);		//outfile bis letzten punkt
			newflag[point]='\0';
		//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_PF.out",8);
		//	printf("flag with flag: %s\n",newflag);
			//~ for (int j=point; j< strlen(ga.outfile) ;j++)
			//~ {
				//~ newflag[3+j]=ga.outfile[j];
				//~ //printf("newflag[%d]=%c \t ga.outfile[%d]=%c\n",4+j,newflag[4+j],j,ga.outfile[j]);
				
			//~ }
					
			ga.outfileFlag = (char *) calloc(strlen(newflag) + 1, sizeof(char));         // für die länge des strings argv[i+1] wird platz reserviert in args.infile
			strncpy(ga.outfileFlag,newflag, strlen(newflag)+1);
			free(newflag);
		}
				
		if (ga.psfile != NULL) //psout_file vorhanden
		{
			int point =0;
			for(int i=0; i<strlen(ga.psfile); i++)
			{
				if(ga.psfile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			if(point==0 || ga.psfile[point-1]=='.' )
			{
				point=strlen(ga.psfile);
			}	
			
			char newflag [100] ;  
			
			
			strncpy(newflag,ga.psfile,point);		//outfile bis letzten punkt
			newflag[point]='\0';
		//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_PF.ps",7);
		//	printf("flag with flag, complete: %s\n",newflag);
						
			ga.psfileFlag = (char *) calloc(strlen(newflag) + 1, sizeof(char)); 
			strncpy(ga.psfileFlag,newflag, strlen(newflag)+1);
			
			
			
		}
		if (ga.psfile == NULL) //psout_file nicht vorhanden
		{
			int point =0;
			for(int i=0; i<strlen(ga.infile); i++)
			{
				if(ga.infile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			if(point==0 || ga.psfile[point-1]=='.' )
			{
				point=strlen(ga.infile);
			}	
			
			char newflag [100] ;  
			strncpy(newflag,ga.infile,point);		//outfile bis letzten punkt
			newflag[point]='\0';
		//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_PF.ps",7);
		//	printf("flag with flag, complete: %s\n",newflag);
			
			
			ga.psfileFlag = (char *) calloc(strlen(newflag) + 1, sizeof(char)); 
			strncpy(ga.psfileFlag,newflag, strlen(newflag)+1);
		}
		
		
		/////////////////////////wenn -d gesetzt, darf auch kein -o outfile und kein -p psfile gesetzt sein, outfileFlag wird zu: dir/infile_PS.out und psfileFlag wird zu dir/infile_PS.ps
		//~ if (ga.directory != NULL)
		//~ {
			//~ //printf("the directory is: %s\n", ga.directory);
			//~ char * outfile = (char *) calloc(100, sizeof(char));
			//~ char * outpsfile = (char *) calloc(100, sizeof(char));
			//~ strncpy(outfile,ga.directory,strlen(ga.directory)+1);
			//~ strncpy(outpsfile,ga.directory,strlen(ga.directory)+1);
			//~ //printf("the new flag:%s\n", outfile);
			
			//~ int point =0;
			//~ for(int i=0; i<strlen(ga.infile); i++)
			//~ {
				//~ if(ga.infile[i]=='.')
				//~ {
					//~ //printf("point found at position %d\n",i);
					//~ point=i;
				//~ }
			//~ }
			//~ if(point==0)
			//~ {
				//~ point=strlen(ga.infile);
			//~ }
			
			//~ strncat(outfile,ga.infile,point);
			//~ strncat(outpsfile,ga.infile,point);
			//~ strncat(outfile,"_PF.out",7);
			//~ strncat(outpsfile,"_PF.ps",6);
			//~ printf("the new flag until point:%s\n", outfile);
			//~ printf("the new psflag until point:%s\n", outpsfile);
					
			//~ free(outfile);
			//~ free(outpsfile);
				
	}
	   
		
	if(ga.b==1)
	{
		
		if(ga.outfile!=NULL)
		{
			int point =0;
			for(int i=0; i<strlen(ga.outfile); i++)
			{
				if(ga.outfile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			if(point==0 || ga.outfile[point-1]=='.' )
			{
				point=strlen(ga.outfile);
			}	
			
			//char newflag [100] ;  //	printf("strlen(newflag): %lu\n", strlen(newflag));
			//char newflagb [100] ; //printf("strlen(newflag): %lu\n", strlen(newflagb));
			
			char *newflag = (char*) calloc(100,sizeof(char));
			char *newflagb = (char*) calloc(100,sizeof(char));
			
			
			strncpy(newflag,ga.outfile,point);
			newflag[point]='\0';
			strncpy(newflagb,ga.outfile,point);		//outfile bis letzten punkt
			newflagb[point]='\0';
			printf("flag until last point: %s\n",newflag);
			printf("flag until last point: %s\n",newflagb);
			strncat(newflag,"_PSb.out",9);
			strncat(newflagb,"_PFb.out",9);
			
			
		//	printf("flag with flag: %s\n",newflag);
		//	printf("flag with flag: %s\n",newflagb);
		//	printf("strlen ga.outfile: %lu \t %s\n",strlen(ga.outfile), ga.outfile);
			//int j=0;
			//~ for ( j=point; j< strlen(ga.outfile) ;j++)
			//~ {
				//~ newflag[4+j]=ga.outfile[j];
				//~ newflagb[4+j]=ga.outfile[j];
				//~ //printf("newflag[%d]=%c \t ga.outfile[%d]=%c\n",4+j,newflag[4+j],j,ga.outfile[j]);
			//~ }
			
			//printf("complete flag: %s, strlen: %lu\n", newflag, strlen(newflag));
		
			
			ga.outfileFlag = (char *) calloc(strlen(newflag) + 1, sizeof(char)); 
			strncpy(ga.outfileFlag,newflag, strlen(newflag)+1);
			
			ga.outfileFlag_b = (char *) calloc(strlen(newflag) + 1, sizeof(char));         // für die länge des strings argv[i+1] wird platz reserviert in args.infile
			strncpy(ga.outfileFlag_b,newflagb, strlen(newflag)+1);
			
			
			free(newflag);
			free(newflagb);
			//ga.outfileFlag_b = (char *) &newflagb;
			//ga.outfileFlag = (char *) &newflag;
			
						
			
	//	printf("new outfileFlag: %s\n", ga.outfileFlag);
			
		//~ //	printf("sizeof ga.outfileFlag: %lu\n", sizeof(ga.outfileFlag));
		//~ //	printf("sizeof newflag: %lu\n", sizeof(newflag));
							
		}
		if(ga.psfile!=NULL)
		{
			int point =0;
			for(int i=0; i<strlen(ga.psfile); i++)
			{
				if(ga.psfile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			if(point==0 || ga.psfile[point-1]=='.' )
			{
				point=strlen(ga.psfile);
			}	
			
			char newflag [strlen(ga.psfile)+9] ;  
			strncpy(newflag,ga.psfile,point);		//outfile bis letzten punkt
			newflag[point]='\0';
		//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_PS_PF.ps",10);
		//	printf("flag with flag, complete: %s\n",newflag);
			
			ga.psfileFlag = (char *) calloc(strlen(newflag) + 2, sizeof(char)); 
			strncpy(ga.psfileFlag,newflag, strlen(newflag)+1);
		

			
		}
		
		if (ga.psfile == NULL) //psout_file nicht vorhanden, nimm infile
		{
			int point =0;
			for(int i=0; i<strlen(ga.infile); i++)
			{
				if(ga.infile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			if(point==0 || ga.infile[point-1]=='.' )
			{
				point=strlen(ga.infile);
			}	
			
			char newflag [100] ;  
			strncpy(newflag,ga.infile,point);		//outfile bis letzten punkt
			newflag[point]='\0';
			
		//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_PS_PF.ps",10);
		//	printf("flag with flag, complete: %s\n",newflag);
			
			ga.psfileFlag = (char *) calloc(strlen(newflag) + 20, sizeof(char)); 
			strncpy(ga.psfileFlag,newflag, strlen(newflag)+1);
			
		/////////////////////////wenn -d gesetzt, darf auch kein -o outfile und kein -p psfile gesetzt sein, outfileFlag wird zu: dir/infile_PS.out und psfileFlag wird zu dir/infile_PS.ps
		//~ if (ga.directory != NULL)
		//~ {
			//~ //printf("the directory is: %s\n", ga.directory);
			//~ char * outfile = (char *) calloc(100, sizeof(char));
			//~ char * outfile_b = (char *) calloc(100, sizeof(char));
			//~ char * outpsfile = (char *) calloc(100, sizeof(char));
			//~ strncpy(outfile,ga.directory,strlen(ga.directory)+1);
			//~ strncpy(outfile_b,ga.directory,strlen(ga.directory)+1);
			//~ strncpy(outpsfile,ga.directory,strlen(ga.directory)+1);
			//~ //printf("the new flag:%s\n", outfile);
			
			//~ int point =0;
			//~ for(int i=0; i<strlen(ga.infile); i++)
			//~ {
				//~ if(ga.infile[i]=='.')
				//~ {
					//~ //printf("point found at position %d\n",i);
					//~ point=i;
				//~ }
			//~ }
			//~ if(point==0)
			//~ {
				//~ point=strlen(ga.infile);
			//~ }
			
			//~ strncat(outfile,ga.infile,point);
			//~ strncat(outfile_b,ga.infile,point);
			//~ strncat(outpsfile,ga.infile,point);
			//~ strncat(outfile,"_PSb.out",8);
			//~ strncat(outfile_b,"_PFb.out",8);
			//~ strncat(outpsfile,"_PS_PF.ps",9);
			//~ printf("the new flag until point:%s\n", outfile);
			//~ printf("the new flagb until point:%s\n", outfile_b);
			//~ printf("the new psflag until point:%s\n", outpsfile);
					
			//~ free(outfile);
			//~ free(outfile_b);
			//~ free(outpsfile);
				
		//~ }				
			
		}
				
	}
	
	if(ga.b==0 && ga.z==0)
	{
		if (ga.collectfile != NULL)
		{
			int point =0;
			for(int i=0; i<strlen(ga.psfile); i++)
			{
				if(ga.psfile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			
			if(point==0 || ga.psfile[point-1]=='.' )
			{
				point=strlen(ga.psfile);
			}	
			//printf("last point at %d\n", point);
			
			char *newflag=(char*)calloc(100,sizeof(char*)); 
			strncpy(newflag,ga.psfile,point);		//outfile bis letzten punkt
			
			
			//char newflag [100] ;  
			//strncpy(newflag,ga.infile,point);		//outfile bis letzten punkt
			//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_collect.ps",12);
		//	printf("flag with flag, complete: %s\n",newflag);
		
			ga.psfileFlag = (char *) calloc(strlen(newflag) + 1, sizeof(char)); 
			strncpy(ga.psfileFlag,newflag, strlen(newflag)+1);
			
			free(newflag);
			return ga;
			
		}
		
		
		if(ga.outfile != NULL)
		{
			int point =0;
			for(int i=0; i<strlen(ga.outfile); i++)
			{
				if(ga.outfile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
		    if(point==0 || ga.outfile[point-1]=='.' )
			{
				point=strlen(ga.outfile);
			}	
			
			//printf("strlen ga.outfile in getfileextension: %lu\n",  );
			
			char *newflag=(char*)calloc(100,sizeof(char*)); 
			strncpy(newflag,ga.outfile,point);		//outfile bis letzten punkt
			newflag[point]='\0';
		//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_PS.out",8);
		//	printf("flag with flag: %s\n",newflag);
			//~ for (int j=point; j< strlen(ga.outfile) ;j++)
			//~ {
				//~ newflag[3+j]=ga.outfile[j];				
			//~ }
		//~ //	printf("completedfdfdfdf flag: %s\n", newflag);
			//newflag[strlen(newflag)]='\0';
			
			ga.outfileFlag = (char *) calloc(strlen(newflag) + 1, sizeof(char)); 
			strncpy(ga.outfileFlag,newflag, strlen(newflag)+1);
			
			free(newflag);
						
		}
		
		if(ga.psfile!=NULL)
		{
			int point =0;
			for(int i=0; i<strlen(ga.psfile); i++)
			{
				if(ga.psfile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			
			if(point==0 || ga.psfile[point-1]=='.' ) 
			{
				point=strlen(ga.psfile);
			}	
		
		//	printf("last point at %d\n", point);
			
			char *newflag=(char*)calloc(100,sizeof(char*)); 
			strncpy(newflag,ga.psfile,point);		//outfile bis letzten punkt
			
			
			//char newflag [100] ;  
			//strncpy(newflag,ga.psfile,point);		//outfile bis letzten punkt
			
		//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_PS.ps",7);
	
		//	printf("flag with flag, complete: %s\n",newflag);
			
			ga.psfileFlag = (char *) calloc(strlen(newflag) + 1, sizeof(char)); 
			strncpy(ga.psfileFlag,newflag, strlen(newflag)+1);
			free(newflag);
		}
		
		if (ga.psfile == NULL) //psout_file nicht vorhanden, nimm infile
		{
			int point =0;
			for(int i=0; i<strlen(ga.infile); i++)
			{
				if(ga.infile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			
			if(point==0 || ga.infile[point-1]=='.' )
			{
				point=strlen(ga.infile);
			}	
		//	printf("last point at %d\n", point);
			
			char *newflag=(char*)calloc(100,sizeof(char*)); 
			strncpy(newflag,ga.infile,point);		//outfile bis letzten punkt
			
			
			//char newflag [100] ;  
			//strncpy(newflag,ga.infile,point);		//outfile bis letzten punkt
			newflag[point]='\0';
		//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_PS.ps",7);
		//	printf("flag with flag, complete: %s\n",newflag);
		
			ga.psfileFlag = (char *) calloc(strlen(newflag) + 1, sizeof(char)); 
			strncpy(ga.psfileFlag,newflag, strlen(newflag)+1);
			
			free(newflag);
		}
		
		
		//////////////////////////setzen der flags wenn -d gegeben:
		
		//~ if (ga.directory != NULL)
		//~ {
			//~ //printf("the directory is: %s\n", ga.directory);
			//~ char * outfile = (char *) calloc(100, sizeof(char));
			//~ char * outpsfile = (char *) calloc(100, sizeof(char));
			//~ strncpy(outfile,ga.directory,strlen(ga.directory)+1);
			//~ strncpy(outpsfile,ga.directory,strlen(ga.directory)+1);
			//~ //printf("the new flag:%s\n", outfile);
			
			//~ int point =0;
			//~ for(int i=0; i<strlen(ga.infile); i++)
			//~ {
				//~ if(ga.infile[i]=='.')
				//~ {
					//~ //printf("point found at position %d\n",i);
					//~ point=i;
				//~ }
			//~ }
			//~ if(point==0)
			//~ {
				//~ point=strlen(ga.infile);
			//~ }
			
			//~ strncat(outfile,ga.infile,point);
			//~ strncat(outpsfile,ga.infile,point);
			//~ strncat(outfile,"_PS.out",7);
			//~ strncat(outpsfile,"_PS.ps",6);
			//~ printf("the new flag until point:%s\n", outfile);
			//~ printf("the new psflag until point:%s\n", outpsfile);
			
			//~ strncpy(ga.outfileFlag,outfile,strlen(outfile)+1);
			//~ strncpy(ga.psfileFlag,outpsfile,strlen(outpsfile)+1);
			
			//~ free(outfile);
			//~ free(outpsfile);		
		
		
		//~ }
		
		if (ga.collectfile != NULL)
		{
			int point =0;
			for(int i=0; i<strlen(ga.psfile); i++)
			{
				if(ga.psfile[i]=='.')
				{
					//printf("point found at position %d\n",i);
					point=i;
				}
			}
			
			if(point==0 || ga.psfile[point-1]=='.' )
			{
				point=strlen(ga.psfile);
			}	
			 //printf("last point at %d\n", point);
			
			char *newflag=(char*)calloc(100,sizeof(char*)); 
			strncpy(newflag,ga.psfile,point);		//outfile bis letzten punkt
			
			
			//char newflag [100] ;  
			//strncpy(newflag,ga.infile,point);		//outfile bis letzten punkt
			//	printf("flag until last point: %s\n",newflag);
			strncat(newflag,"_collect.ps",12);
		//	printf("flag with flag, complete: %s\n",newflag);
		
			ga.psfileFlag = (char *) calloc(strlen(newflag) + 1, sizeof(char)); 
			strncpy(ga.psfileFlag,newflag, strlen(newflag)+1);
			
			free(newflag);
			
			
		}
			
	}	
	
	
			
	return ga;
}			
		

	
	
	

