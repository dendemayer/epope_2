#include "newicktolist2testbaumok.h"

#undef NEWICKD

static int lab=0;
static int d2root=0;
static int return_on =0;

void newick(char * filein)
{	
	
	//char filein[]="newick.tre";
	test(filein);					//zaehlt öffnende und schließende klammern
	int l = count(filein);		 //zählt länge des arrays aus newickformat, ohne eoa \0 zeichen
	int *pl = &l;
	//printf("count:\t%d\n", l);
	//char sp[l];
	//char *p_sp=sp;
	
	
				//newickarray soll in die main:
	//printf("l:%d\n", l);
	//char newickarray[2000];
	char* newickarray = (char*) calloc(2000, sizeof(char));
	//char *pnewickarray = newickarray;
	newicktoarray(filein, pl, newickarray);   //l mus noch angepasst werden, um evt distanzen kürzen wenn fälle wie: (,   ,,   ),  auftreten, muss dazwischen 
		
	arrayplus(pl,newickarray); //dummy lab gesetzt werden, pl wird geändert
//	printf("newickhääääää: %s\n",newickarray); 	
//	printf("newick vor setroot: %s\n",newickarray); 	
//	printf("*pl:%d\n", *pl);															//hier noch ok
	//=strlen(newickarray); 			//diese fkt scheint array zu ändern -> l wird innerhalb von arrayplus angepasst...
	
	
	
	//printf("laenge: %d\n",l);
	//exit (0);
	//printf("um distanzen gekürzte länge: %d \n",l);
/*	printf("Newickarray aus erstelltem array:\n");
	printf("%s", newickarray);					
	printf("\n");
	printf(" newickarray[%d]:%c\n",l,newickarray[l]);
*/	
	//hier neue funktion, die an newickarray vor ';' ROOT einfügt und l um 4 verlängert

	setroot(newickarray, *pl);									//hier falsch
//	printf("newick nach setroot: %s\n",newickarray); 	
	
	l=strlen(newickarray);
	//printf("neues l aus strlen: %d\n",l);
	
//	printf(" neues l:%d\n",l);
//	printf(" newickarray[%d]:%c\n",l,newickarray[l-1]);
	
	int spcount = countsp(newickarray, l);
//	printf("spcount %d\n", spcount);
	
	
	
	///////////bis hier ist nur newickstring in array geschrieben worden, struct für knoten anlegen://///
	
	/*
	nodeNEW test;
	strcpy(test.parent , "parentname");
	*/
	
	///array aus nodeNEW anlegen so lang wie spcount
	
	//nodeNEW list[spcount+1];						//array der länge |species| aus sutruct nodeNEW
	
	
	nodeNEW *list= (nodeNEW*) calloc(spcount+1, sizeof(nodeNEW));
	//initialisieren;
	
	

	//printf("Adresse der list in main:%p\n", &list);
	
	//nodeNEW *plist = list;
/*	
	strcpy(spstructar[0].parent, "parentone");
	printf("parentone: %s\n", spstructar[0].parent);
*/
	//printf("pos pnewickarray:%c\n", *pnewickarray);
	
//	printf("array vor find leaves: %s\n", newickarray);
	
	findleaves(list, newickarray, l, filein, spcount );
	
//	printf("nach findleaves\n");
	free(list);
	return;
}

char* setroot(char *newickarray, int  pl)
{
	//pl kommt falsch an!?!?! setlen bestimmt letzte pos!
	
	
//	printf("newickarray in setroot: %s\n", newickarray);
//	printf(" last pos in setroot: %d\n", pl);
//	printf("last with strlen; %lu\n  ",strlen(newickarray));
	
	pl=strlen(newickarray);
	
//	printf("newickarray[%d-1]: %c \n newickarray[%d]: %c \n",pl,newickarray[pl-2],pl,newickarray[pl-1] );
	if (newickarray[pl-2]==')' && newickarray[pl-1]==';')
	{
		//printf("root muss eingefügt werden");
		char root[]="Root;";
		newickarray[pl-1]='\0';
		strncat(newickarray, root, 5);
		//printf(" neues array: %s\n", newickarray);
	}	
	return newickarray;
}




	////start find leaves und belegung der ersten nodeNEW struct elemente
	
void findleaves(nodeNEW list[], char newickarray[], int l, char* filein, int spcount)
{	
	//printf("Adresse der list in findleaves:%p\n", &list);
//nodeNEW *plist=list;	
	//printf("newickarray:%s\n", newickarray);
	//printf("strlennewickarray: %lu\n", strlen(newickarray));
	
	l=strlen(newickarray);
	//printf("l in findleaves: %d\n",l);
	
	//printf("newickarray in findleaves: %s\n",newickarray);
	
	for (int i=0; i<l; i++)
	{
		if (newickarray[i]=='(')			//solange Klammern gefunden werden, gibts noch spezies im tree
		{
			break;
		}
	
		else
		{	if (strlen(newickarray)>1)		//letzten knoten(root) einfügen
			{	
		//		printf("\n\n\noder hier\n\n\n\n\n");
			//	printf("\n\n\%s\n\n\n\n\n", newickarray);
				//char name[strlen(newickarray)];
				char *name = (char*) calloc(strlen(newickarray), sizeof(char));
				char dad[]="-1";
				int j;
				for (j=0; j<strlen(newickarray)-1;j++)
					name[j]=newickarray[j];
				name[j]='\0';
		//		printf("name:%s\n",name);
				list[lab] = writeleave(name,dad,lab,d2root,list[lab]);
				//printf("\nd2root: %d, filein: %s, spcount: %d\n\n",d2root ,filein ,spcount );
				free(newickarray);
				free(name);
			}
			
			else //if (newickarray[0]==';')
				{	//printf("sind wir hier???\n");
					char name[]="ROOT";
					char dad[]="-1";
				list[lab] = writeleave(name,dad,lab,d2root,list[lab]);
				}
			//printf("\nd2root: %d, filein: %s, spcount: %d\n\n",d2root ,filein ,spcount );	
			print(list, d2root, filein, spcount);
			//free(newickarray);
			int return_on =1;
		//	printf("newickarray:%s\n",newickarray);
			if (return_on==1)			//sonst unused value	
			return;
						
			
			
		}
	}	
	
	for (int i=0; i<l; i++)
	{	//printf("i in erster for schleife: %d\n", i);
		if (newickarray[i]=='(')	/*wenn öffnende klammer gefunden wird, 
										ist ihre dazugehörige schließende
									nicht unterbrochen von einer weiteren öffnenden
									-> alles dazwischen sind leaves*/
		{	
			for (int j=i+1; j<l; j++)
			{	//printf("j in zweiter for schleife: %d\n",j);
				if (newickarray[j]=='(')
					{
						i=j-1;
						 
						break;
					}
				
				else if (  newickarray[j]==')')
				{	//printf("\nleave von %d bis %d \n", i, j);
				/*	
					printf(" speziestest mit komma:\n");
					for (int t = i+1;t<j;t++)
					printf("%c", newickarray[t]);
					printf("\n");
					* 
				*/
					/*bis hier gefunden, öffnende und schließende klammern
					jetzt noch bereiche zwischen kommas aufdröseln*/
					int kposalt=i+1;
					for (int k=i+1; k<=j;k++)
					{	if (newickarray[k+1]==',')		//wenn nach leave löschung zwei kommas hintereinander...
							continue;
						if (newickarray[k]==',' || newickarray[k]==')')
						{	
							
							int namepos =0;    
							int leavenamelength= k-kposalt;
							char leavename[leavenamelength];
							leavename[leavenamelength]='\0';				//laenge string festlegen
							//printf("stringlenth: %d\n", k-kposalt);
							for (int l = kposalt; l<k; l++)  	//create leavename string
								{
								 	leavename[namepos]=newickarray[l];
								 	namepos++;
								}
							
						
						/*	printf("leavenametest:\n\n");
							for (int n=0; n<leavenamelength; n++)
							printf("%c",leavename[n]);
							printf("\n\n");
						*/	
								
							///////////////////////////////////////////////////////////ab hier werden leaves in strings geschrieben, was ist mit vätern?
							//neue schleife, es wird davon ausgegangen, dass väter immer rechts nach schließender klammer stehen, (pos j schon bekannt)
							// z.B. (222,333)Virilisgr , also lesen ab ')' bis (c=='(' || c==')' ||c==',' || c==';')  falls c==leer -> ROOT
							// 
							//erstmal laenge dadname testen
						
						/*	
							printf("j: %d \n",j);
							printf("count(): %d\n", count());
							printf("count()-j=%d\n", count()-j);
							printf("newickarray[j+1]=%c\n",newickarray[j+1]);
						*/	
						//	printf("strlen(leavename):%lu\n", strlen(leavename));
							int dadlengthcount;
							int k_dad=j+1;
							for (dadlengthcount = 0; dadlengthcount < l-j; dadlengthcount++ )
							{	
								if (newickarray[k_dad]=='(' || newickarray[k_dad]==')' ||newickarray[k_dad]==',' || newickarray[k_dad]==';')
									break;
								k_dad++;
							}		
							
						
							//printf("dadlengthcount: %d\n", dadlengthcount);  //passt
							
								
							//dadname in string schreiben:
														
							//char dadname[dadlengthcount+1];				//laenge string festlegen
							char *dadname=(char*)calloc(dadlengthcount+1,sizeof(char));
							dadname[dadlengthcount]='\0';
						/*		dadname[0]=newickarray[j+1];
								dadname[0+dadlengthcount]=newickarray[j+1+dadlengthcount];
								printf("dadnamefirst: %c\n", dadname[0]);
								printf("dadnamelast: %c\n", dadname[0+dadlengthcount]);
						*/	
							int nameposdad =0;    
							for (int k_dadtostringstart=j+1; k_dadtostringstart<j+1+dadlengthcount; k_dadtostringstart++)  	//create dadname string
								{
								 	dadname[nameposdad]=newickarray[k_dadtostringstart];
								 	nameposdad++;
								}
								 
													//~ printf("dadnametest:\n\n");
													//~ for (int n=0; n<dadlengthcount; n++)
													//~ printf("%c",dadname[n]);
													//~ printf("\nstrleng: %lu\n",strlen(dadname));
													//~ printf("\n\n");
						 
							////// leaves und deren väter bekannt -> erzeugen der untersten structs
							//char *pleavename=leavename;
							//char *pdadname=dadname;
							if (dadlengthcount!=0)
							list[lab] = writeleave(leavename,dadname,lab,d2root, list[lab]);
							else
							{
								char tempdadname[]="ROOT";
								list[lab] = writeleave(leavename,tempdadname,lab,d2root, list[lab]);
							}
							
							lab++;
						
							kposalt=k+1;
							
							free(dadname);
												
						}
					
					}break;
				}
			}
		}
	}	
	
//neue funktion cutarray, diese übergibt neues array an find leaves	
d2root++;
//char *pnewickarray=newickarray;

//hier noch mal list überprüfen:
/*
printf("lab in findleaves:\t%d\n", lab);

for (int i=0;i<lab;i++)
	printf("lab(findleaves:\t\t%d\n",list[i].lab);
*/

//printf("newick vor cut: %s\n", newickarray); //hier schon falsch

cutarray(newickarray,list, filein, spcount );
}
	
	////print////////
	
void print(nodeNEW list[],int d2root, char *filein, int l)
{
	//int l=countsp(filein);
	
	/*
	char empty[]="-1";
	for (int i=0;i<l;i++)				//leaves get -1 in child string array
		if (list[i].d2root==0)
			strcpy(list[i].children, empty);
		else
		list[i].children[0]='\0';		//warum auch immer, sonst fehler bei Dro, Repag und Sp9
	*/
	for (int i=0;i<l;i++)			
	list[i].children[0]='\0';     //warum auch immer, sonst undefinierte zeichen in child array, wenn es nicht -1 gesetzt wird
	
	for (int i=0; i<l; i++)			//nodeNEWs getting their children numbers
	{	int foundchild=0;
		
		for (int j=0; j<l; j++)
		{	
			
			if (strcmp(list[i].name, list[j].parent) == 0)	//wenn der name eines nodeNEWs(i) woanders als parent(j) geführt wird
			{	
				foundchild=1;											//wird lab(j) als children(i) eingetragen
				char temp[4];
				
				sprintf(temp, "%d,", list[j].lab); 			//lab(j) wird in temp string geschrieben
				//printf("temp: %s \t j:%d \ti:%d\n", temp, j,i );
				strcat(list[i].children,temp);				//an bestehendes children(i) string wird lab(j) angehängt		
				//printf("list[i.children]:\t%s\n", list[i].children);
				//printf("Name: %s\tlab:%d\tparent:%s\tlab:%d\n", list[i].name, list[i].lab, list[j].parent, list[j].lab);
			}
			
		}
		if (foundchild==0)
			{
				//printf("-1 ist hier zu setzen: %s\n ", list[i].name);
				char empty[]="-1";
				strcpy(list[i].children, empty);
				
			}
		//letztes komma muss gelöscht werden:
		//printf("%s\n", list[i].children);
		char komma =','; 
		char *find;
		find = strrchr(list[i].children, komma);
		if( find != NULL)
		{
			*find = '\0';
			/*
			printf("%c\n", *find);
			find = find -1;
			printf("%c\n", *find);
			*/
		}
		//printf("%s\n", list[i].children);
	}	
		
	for (int i=0; i<l; i++)		
		for (int j=0; j<l; j++)									//nodeNEWs getting their par numbers
		{	if (i==j)
			continue;
			if (strcmp(list[i].parent, list[j].name) == 0)		//parent[](i) ist schon geschrieben, einfach lab von parent in par eintragen
				list[i].par=list[j].lab;
				//printf("list[%d].par:\t%d\n",i, list[i].par);
		}	
		
		
	//für root noch ma extra, dort ist char dad[] = "-1" gesetzt:
	for (int i=0; i<l; i++)		
		if (strcmp(list[i].parent, "-1") == 0)	
			list[i].par= -1;					
		
	//root anpassen:
	
	int root =0;
	int i;		//i braucht man außerhalb for schleife, ist root nodeNEW
	for (i=0; i<l-1; i++)		
		{
			if (strcmp(list[i].parent, "-1") == 0)	
				list[i].d2root= root;	
				//printf("list[i].d2root=%d\n", list[i].d2root);
				
		}	
		
		//bis hier nodeNEW array pos root gefunden und d2root auf 0 gesetzt
		//bei allen anderen müssen d2root aktualisiert werden
	//printf("rootlab for findroot aufruf: %d\n",list[i-1].lab );
	//printf("i: %d\n",i-1 );	


    
	findroot(list[i].lab, root, list, filein, l);
	
FILE  *pfile;
		pfile  =  fopen  ("newick_to_nodelist.dat",  "w");
		if  (pfile  ==  NULL)
		{
			printf("Fehler  beim  oeffnen  der  Datei.");
			
		}		
/*	
	for (int i=0; i<l; i++)	
	{
		printf("%s ", list[i].name);
		printf("lab %d ", list[i].lab);
		printf("par %d ", list[i].par);
		printf("children %s ", list[i].children);
		printf("d2root %d \n", list[i].d2root);
	}
*/			
	for (int i=0; i<l; i++)
	{	
		if (strcmp(list[i].parent, "-1") == 0)	//FEHLER finden! d2root setzt sich sonst auf 3
				list[i].d2root= 0;				//hier wieder auf 0 gesetzt
			
/*		//printf("%s ", list[i].name);
		for (int j =0; j<strlen(list[i].children);j++)
			printf("%c",list[i].children[j]);
			 
		
		printf("%s ", list[i].name);
		printf("lab %d ", list[i].lab);
		printf("par %d ", list[i].par);
		printf("children %s ", list[i].children);
		printf("d2root %d ", list[i].d2root);
		printf("gainW 0 lossW 0\n");
*/		
		fprintf(pfile,"%s ", list[i].name);
		fprintf(pfile,"lab %d ", list[i].lab);
		fprintf(pfile,"par %d ", list[i].par);
		fprintf(pfile,"children %s ", list[i].children);
		fprintf(pfile,"d2root %d ", list[i].d2root);
		fprintf(pfile,"gainW 0 lossW 0\n");
		
	}
	fclose  (pfile);
	return;

}		
	
	////////////////findroot
void findroot(int root_lab, int root, nodeNEW list[], char *filein, int l)
{
	//printf("rootlab: %d\nroot: %d\n\n", root_lab,root);
	//int l = countsp(filein);
	
	
	root++;
	//~ for (int i=0;i<l-1;i++)
		//~ printf("list[%d].par:\t%d\n",i,list[i].par);
	
	
	for (int i=0; i<l-1; i++)
	{
		if (list[i].par==root_lab)
			{	
				list[i].d2root=root;
				
				if (strcmp(list[i].children, "-1") == 0)
				continue;
				
				else
				findroot(list[i].lab, root, list, filein, l);
			}
			
		
	}
	
	
	
}

	////cut array
	
void cutarray(char newickarray[], nodeNEW list[], char *filein, int spcount)
{	//printf("l als count:\t%d\n",count());
	//printf("l als strlen:\t%lu\n", strlen(newickarray));
	if (return_on == 1)
	return;
	
	#ifdef NEWICKD
	printf("newick in cut: %s \n", newickarray);
	#endif
	int l = strlen(newickarray);
	
	for (int i=0; i<l; i++)
	{	
		if (newickarray[i]=='(')	
		{	
			for (int j=i+1; j<l; j++)
			{	//printf("j in zweiter for schleife: %d\n",j);
				if (newickarray[j]=='(')
				{
					i=j-1;
					break;
				}
				
				else if (  newickarray[j]==')')
				{
					newickarray[i]='!';
					newickarray[j]='!';
					i=j+1;
					break;					
				}
			}
		}
	}
/*	
	printf("\ncutarray:\n");	
	for(int i=0; i<l; i++)
		printf("%c",newickarray[i]);
		printf("\n");
*/		
	char *bnewickarray = (char*) calloc(l+1,sizeof(char));	
	
	
	int m =0;
	for (int i=0; i<l; i++)
		{	
			if (newickarray[i]=='!')
			{	//printf("i:%d\n",i);
				for (int j=i+1;j<l;j++)
					if (newickarray[j]=='!')
					{
						//printf("j:%d\n",j);
						i=j;
						break;
					}	
			}
			else
			{	
				bnewickarray[m]=newickarray[i];	
				m++;
			}	
		}
		
	for (int i=0; i<strlen(bnewickarray);i++)
		if (bnewickarray[i]==';')
		{
			bnewickarray[i+1]='\0';
			break;
		}
	
	#ifdef NEWICKD
	printf("\ncutarray:\n");	
	for(int i=0; i<strlen(bnewickarray); i++)
		printf("%c",bnewickarray[i]);	
	printf("\n");
	#endif
	
	int len = strlen(bnewickarray);
	//char *pnewickarray=bnewickarray;
	//nodeNEW *plist=list;
	
	//~ printf("newick before find leaves; %s\n",bnewickarray);
	//~ printf("len: %d\n filein: %s\n spcount: %d\n",len,filein, spcount);
	
	
	free(newickarray);
	
	findleaves(list,bnewickarray,len, filein,spcount );
}	



	
//////////writeleaves////////

nodeNEW writeleave(char leavename[],char dadname[], int lab,int d2root, nodeNEW leave)
{	
	if (return_on  ==1)
	return leave;
	//test ob alles da ist:
	//~ printf("leavenametest:\t");
							//~ for (int n=0; n<(strlen(leavename)); n++)
							//~ printf("%c",leavename[n]);
							//~ printf("\n");
	//~ printf("dadnametest:\t");
							//~ for (int n=0; n<strlen(dadname); n++)
							//~ printf("%c",dadname[n]);
							//~ printf("\n");
	//~ printf("labnr.:\t %d\n\n",lab);
	
	
	strncpy(leave.name, leavename, strlen(leavename)+1);
//	printf("strlen(leave.name):%lu\n",strlen(leave.name));
//	printf("leave.name:\t%s\n", leave.name);
	
/*	printf("strcpytest:\t");
	for (int i=0; i<strlen(leave.child);i++)
		printf("%c", leave.child[i]);
	printf("\n");
*/	

	strncpy(leave.parent, dadname, strlen(dadname)+1);
//	printf("strlen(leave.parent):%lu\n",strlen(leave.parent));
//	printf("leave.parent:\t%s\n", leave.parent);
	
	leave.lab=lab;
//	printf("leave.lab:\t%d\n", leave.lab);
	
	
	leave.d2root=d2root;
//	printf("leave.d2root:\t%d\n\n\n", leave.d2root);
	
	return leave;
}
	

	///////////////////////////////////test auf grobe konsistenz
void test(char *filein)
{	
	FILE * pfile;
	pfile = fopen(filein, "r");
	char c;
	int auf =0;
	int zu =0;
	while(1)
   {
      c = fgetc(pfile);
      if( feof(pfile) )
      { 
		break ;
      }
		//printf("%c", c);		
		
		if (c=='(')
		auf++;
		if (c==')')
		zu++;
   }
  //printf("\n%d\n", auf);
  // printf("%d\n", zu);
   
   if (auf != zu)
   {
		printf("something is wrong with the newick file, you got %d opening and %d closing braces\nexiting now\n", auf,zu );
		exit(0);
		
	}
 
	fclose(pfile);
}		
		
	/////// laenge zählen: 
	// um groesse des char[count] arrays festzulegen
	//weitergehen, bis 1. spezies

int count(char *filein)
{
//	printf("filein:%s\n",filein);
	
	FILE * pfile;
	pfile = fopen(filein, "r");
	int count=0;
	char c;
	while(1)
	{ 
		c=fgetc(pfile);
		if(c==':')
		{
			count++;
			//printf("%c bei %d\n",c ,count );
			
		}
		
		else if( feof(pfile) )
		{ 
			fclose(pfile);
				//printf("count %d\n", count);
			break ;
		}
		else if(c==';')
		{	count++;
			break;
		}
		else 
		{
			count++;
			//printf("c ist: %c count ist :%d\n",c,count);
		}
	}
	fclose(pfile);
	return count-1;
	
	
}
	/////////////ende laenge zaehlen
	
	
	
	////spezies zaehlen
	
int countsp(char *newickarray, int l)
{
	//FILE * pfile;
	//pfile = fopen(filein, "r");
	
//	printf("countsp array: %s\n", newickarray);
	
	int count = 0;
	char c;
	int i;
	int a=0;
	int j;
	startlabel:
	for ( i = a; i<l;i++)	//weitergehen, bis 1. spezies
	{
		c = newickarray[i];
	/*	printf("erster loop c:\t%c\n", c);
		printf("i: %d\n ", i); 
		printf("j: %d\n ", j);
		printf("count: %d\n ", count);
	*/
		if( i == l-1)
				return count;
				//printf("countsp %d\n", count);
			//exit(0);
			//break ;
		
		if(c=='(' || c==')' ||c==',' || c==';')
		continue;
		else 
		{	count++;
			break;
		}	
		
	}
	
	for ( j=i+1; j<l;j++)	//weitergehen bis spezies aufhört
	{	
		
		c = newickarray[j];
		if( j==l-1)
			return count;
			
		if(c=='(' || c==')' ||c==',' || c==';')
		{	
			break;
		}
		else 
		{	
			
	/*	printf(" zweiter loop sp %d \tc: %c\n",count, c);
		printf("i: %d\n ", i); 
		printf("j: %d\n ", j);
		printf("count: %d\n ", count);
	*/
			a=j+1;
			continue;
		}
		
	}
	goto startlabel;
}
	///////////ende anzahl spezies
void newicktoarray(char *filein, int *laenge, char array[])
{	FILE * pfile;
	pfile = fopen(filein, "r");
	int i;
	//distanzen haben die form: 
	/*(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);  distances and leaf names (popular)
	 * (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;   distances and all names
	 * (:0.1,:0.2,(:0.3,:0.4):0.5);           all but root node have a distance to parent
	 * (:0.1,:0.2,(:0.3,:0.4):0.5):0.0;       all have a distance to parent
	 *   also muss alles raus, von einschl. doppelpunkt bis ausschließlich komma
	 * 
	 * ((FFF:0.1,GGG:0.2)CCC:0.4,(EEE:0.5,DDD:0.6)BBB:0.7)AAA:0.8;
	 * :....    , |  ) |  ; 
	 */
	 
	 int cutcount = 0;		//nur zur info, eigtl nicht wichtig
	 int cutcounttemp = 0;	//nur zur info, eigtl nicht wichtig
	 char c;
	// printf("lange: %d\n",*laenge);
	for (i=0;i<*laenge+1;i++)
	{	//in das newick array sollen keine distanzen, länge der zeichen die fehlen werden dann gekürzt
		c=fgetc(pfile);
		if (c!=':')
		{
			array[i]=c;
			//printf("array ohne kkkdistanz: %c bei i= %d\n",c,i);
			//wenn keine distanzen vorkommen, muss hier 
		}		
		else //hier ist : gefunden, eintragen in array überspringen bis , |  ) |  ;  //nach austritt muss i um eins verringert werden, damit : gelöscht wird
		{
			while(1)
			{
				c=fgetc(pfile);
				if((c==',') | (c==')') | (c==';'))
				{
					array[i]=c;
					break; //while schleife wird abgebrochen
				}
				else 
				cutcount++;
				cutcounttemp++;
			}
			
			//printf(" c nach distanz: %c bei cutcount %d cutcount temp: %d\n",c,cutcount, cutcounttemp);
			//noch mal check nach ';', da dann auch for schleife abgevrochen werden muss!, zeichen nach ; auf \setzen
			//
			
		}
		if(c==';')
			{
				//printf("newickarray beendet bei c= %c , bei i= %d \n ",c,i);
				array[i]=c; //also ';'
				array[i+1]='\0';
				break;
				
			}
			cutcounttemp=0;
			//i--;
		
		
	}	
	
	
/*	printf("Arraytest:\n");
	int j=0;
	while(array[j]!='\0')
	{
		printf("%c",array[j]);
		j++;
	}

	laenge=j;
	 
*/	
	//printf("\n\nj ist %d und damit die neue laenge %d!\n\n", j,*laenge);
	
		
	fclose(pfile);
	
}

void arrayplus(int *laenge, char newickarray[])
{
	int dummycount=0;
	//int templen=*laenge+3; 
	for (int i=0; i<*laenge;i++)
	{
		if ((newickarray[i]=='(' && newickarray[i+1]== ',') || (newickarray[i]==',' && newickarray[i+1]==',') || (newickarray[i]== ')' && newickarray[i+1]== ',') ||  (newickarray[i]== ',' && newickarray[i+1]== ')')  ||  (newickarray[i]== ')' && newickarray[i+1]== ')'))
		{	
			dummycount++;
		}	
	
	}	
	if (dummycount==0)
	{
		//printf("sdjkfnsdjkbf\n");
		return ;
	}
	
	//	printf("\ndummycount=%d\n", dummycount);
		
	//printf("\ntime for a dummy after %c and %c\n",newickarray[i],newickarray[i+1]);
	char newickarray2[*laenge+dummycount*3];
	char lab[]="lab";
	//printf("old length:%d\tnew length:%d\n", *laenge, *laenge+ dummycount*3);
	//printf("lab = %s,  strlen(newickarray2) = %lu, sizeof(newickarray2)= %lu  strlen(lab)=%lu\n",lab,strlen(newickarray2), sizeof(newickarray2),strlen(lab));
	
	//erste schleife liest altes array, bis dummy pos i, i muss festgehalten werden, bei zweiten dummy von i alt bis i array übernehmen, dann u.s.w.
	int itemp= 0;
	int dummyvorhanden=0;
	int j=0;
	for (int i=0; i<*laenge;i++)//newickarray wird durchsucht
	{
		if ((newickarray[i]=='(' && newickarray[i+1]== ',') || (newickarray[i]==',' && newickarray[i+1]==',') || (newickarray[i]== ')' && newickarray[i+1]== ',') ||  (newickarray[i]== ',' && newickarray[i+1]== ')') ||  (newickarray[i]== ')' && newickarray[i+1]== ')') )
		{	
		//	printf("itemp=%d\n",itemp);
		//	printf("i=    %d\n",i);		
			//von 0 bis i kann array übernommen werden, dann lab, weitersuchen...
			if(itemp!=0)
			itemp+=1;												//FEHLER
			for (j=itemp ; j<=i ;j++)
			{
		//		printf("newickarray[%d]=%c\n",j,newickarray[j]);
				newickarray2[j+dummyvorhanden*strlen(lab)]=newickarray[j];
				newickarray2[j+dummyvorhanden*strlen(lab)+1]='\0';
				//printf("newickarray2  %s\n",newickarray2);
				//printf("newickarray%s\n",newickarray);
							
				
			}
			
			itemp=i; 	//i wird festgehalten
		/*	int y=0;
			for (int x=0;x<100;x++)
			{
				y=x%10;
				printf(" modulo von %d ist %d\n ", x, y);
				printf(" 10 er  von %d ist %d\n ", x, x/10);
			}
		*/	
			newickarray2[j+dummyvorhanden*strlen(lab)+0]='L';
			
			if (dummyvorhanden < 10)
			{
				newickarray2[j+dummyvorhanden*strlen(lab)+1]='0' ;
				newickarray2[j+dummyvorhanden*strlen(lab)+2]=dummyvorhanden + '0'; //mit +'0' wird aus int ein char 
			}
			else
			{
				newickarray2[j+dummyvorhanden*strlen(lab)+1]=dummyvorhanden/10 +'0'; 	//10er stelle wird gesetzt
				newickarray2[j+dummyvorhanden*strlen(lab)+2]=dummyvorhanden%10 +'0';	//einer stelle wird gesetzt
			}
			
					//zugefügte labs werden durch nummeriert
			
			
		/*			
			newickarray2[j+dummyvorhanden*strlen(lab)+1]='0';
			newickarray2[j+dummyvorhanden*strlen(lab)+2]=dummyvorhanden +'0';		
			newickarray2[j+dummyvorhanden*strlen(lab)+3]='\0';
			*/ 
			//printf(" newickarray 2 bei i=%d ist %s\n",i,newickarray2);
			dummyvorhanden++;
			//printf("dummyvorhanden:%d\tdummycount:%d\n",dummyvorhanden,dummycount);
			if(dummyvorhanden==dummycount) //das letzte auftreten wurde bereits gefunden, der rest muss jetzt noch in array 2 kopiert werden:
			{//von aktuelem i bis zu schluß
				int k=0;
				for (k=i;k<=strlen(newickarray)-1;k++)
				{
					//printf("k=%d\n",k);
					//printf("strlen(newickarray)=%lu\n",strlen(newickarray));
					
					newickarray2[k+dummyvorhanden*strlen(lab)+1]=newickarray[k+1];
					newickarray2[k+dummyvorhanden*strlen(lab)+2]='\0';
				//	printf("newickarray2; %s\n", newickarray2);
				}
				//newickarray2[k+dummyvorhanden*strlen(lab)+1]='\0';
				int lastpos=strlen(newickarray2);
			//	printf("vorletzte pos: %s\n\n",newickarray2);
			//	printf("strlen: %lu\n",strlen(newickarray2));
			//	printf("char[strlen-1]: %c\n\n",newickarray2[strlen(newickarray2)-1]);
				int *pl=&lastpos;
			//	printf("laenge vor: %d\n\n", *laenge);
				laenge=pl;
			//	printf("laenge nach: %d\n\n", *laenge);
				
			}
			
			
		}
		
	}	
//	printf("newickarray2; %s\n", newickarray2);
	//exit(0);	
	
	
//	printf("laenge nach2: %d\n", *laenge);
	memmove(newickarray, newickarray2, strlen(newickarray2));
	//printf("newickarray:%s\n",newickarray);
	//jetzt müssen noch 
	
	
	
	//
		
}

