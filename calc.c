



#include "calc.h"

extern struct gl_nodeN treeN[];
extern int nodesN;
extern int treeOrderN[];
extern int treePosN[]; /* position of  nodes in treeN[] according to preorder */

extern float treeGainWeight[];
extern float treeLossWeight[];

int FIRST;
int LAST; 
#define BETA -1

#undef PRINTTABLE

//#define zkMax 3
/* -------------------------------------------------------------------------- */

/* do forward and backward recursion */
void gl_calc(struct gl_arguments ga) 
{	


	
  int lca, lca_pos, kMin, kMax, m, i, j;	
  
  lca = gl_getLCA(&kMin, &kMax);
  	

  if(ga.collectfile != NULL) 
  {						
     // label of root of subtree 
    //float **pku=NULL;
    gl_printTreePSN(kMax, lca, ga);
    gl_printTreeN(ga);
    return;						
  }




  //  int done[N_NODES];
  float **S=NULL, inf=9999., SS; /*, p, q, x;*/
  int k;/*kk, left, right, l, leaves[N_NODES];*/
														//char * psfname=NULL;
  if(lca == -1) { printf("Your species are not in the tree\nExiting program.\n"); return; }
  
  lca_pos = treePosN[lca]; //lca_pos ist die einlesepos. (zeilennummer) des LCA
  //printf("LCA:   root of subtree: pInT: %d label: %d = %s\n", lca_pos ,lca , treeN[lca_pos].o);
  //printf("kMin: %d kMax: %d\n", kMin, kMax);
  
  /* get position in treeN of right-most last node */
  //r = treePosN[treeOrderN[nodesN-1]];
  
  /* allocate Score array k-> rows (# mirs), j -> cols (nodes, in preorder)*/
  S = (float **) calloc(kMax+1, sizeof(float *));				// calloc für arrays, erstes argument gibt array länge an
  for(i=0; i<=kMax; i++) {
    S[i] = (float *) calloc(nodesN+1, sizeof(float));
    for(j=0; j<nodesN; j++) {
      S[i][j] = -1.; /* preorder does not matter jet */
    }
  }
  
  /* initialize leaves */
  //  l = 0;
 // printf("kMax = %d\n",kMax);
	for(j=0; j<nodesN; j++) 
	{
		if(treeN[j].nc == 0) 
		{ 	/* this is a leaf node for it has no children */
			//      leaves[l] = j; l++;    /* remember order of leaves - do i need that?*/
			for(i=0; i<=kMax; i++) 
			{ /* check the number of miRNAs at this leaf */
				if(treeN[j].m == i) S[i][j] = 0.;
				else S[i][j] = inf;
			}
		} 
	}
			/*
							printf("\nScore array after initialization (before forward recursion):\n"); 
													printf("\nS_kl:\n");
													gl_printS(S, kMax+1, nodesN);
			*/
	int c,pos, cpos, mink;
	float minscore,score,sum;
	
	/* start forward recursion in lca to get optimal score */
	/* work on positions in tree, parse labels to positions before */
	lca_pos = treePosN[lca];
	j = lca_pos; /* j is position in treeN */
	while(S[kMax][lca_pos] == -1) 
	{
		/* printf("FW: j=%d nc: %d treeN[%d].o=%s\n",j,treeN[j].nc, j,treeN[j].o); fflush(stdout); */
		if(treeN[j].nc == 0) 
		{ /* this is a leaf, that has already been scored */
			treeN[j].seenFW = 1;
			j = treePosN[(treeN[j].p)->n];
			goto nextWhileFW;
		}

		/* go through all child nodes, dig into those that have not been scored so far */
		for(i=0; i<treeN[j].nc; i++) 
		{
			pos = treePosN[treeN[j].c[i]];
			/* if(S[1][pos] == -1) {  *//* child that has not been scored so far */
			if(treeN[pos].seenFW == -1) 
			{
				j = pos;
				/* printf("go down to child to pos: %d\n",j); */
				goto nextWhileFW;
			}
		}
		/* all childnodes have been scored - check inner node */
		if(treeN[j].P == 0) 
		{ 
			S[0][j] = 0; 
		}
		else 
		{ 
			S[0][j] = inf; 
		}
		for(i=1; i<=kMax; i++) 
		{ /* for each row */
			sum=0;
			for(c=0; c<treeN[j].nc; c++) 
			{
				/* child i */  //neee das ist eher child j
				pos = treePosN[treeN[j].c[c]]; /* position of child in tree and in j of S[k][j]*/
				/* printf("position: pos = treePosN[%d] = %d\n",treeN[j].c[c],treePosN[treeN[j].c[c]]); */
				minscore = inf;
				score = inf;
				for(k=0; k<=kMax; k++) 
				{
					score = S[k][pos] + gl_delta(i,k,pos);
					//printf("score=%f + %f = %f (i: %d k: %d pos: %d)\n",S[k][pos],gl_delta(i,k,pos),score,i,k,pos); 
					if(score < minscore) 
					{ 
						minscore = score; 
					}
				}
				sum += minscore;
			}
			S[i][j] = sum; 
			
			//gl_printS(S, kMax+1, nodesN); 
			
		}
		treeN[j].seenFW = 1;
		
		/* END of scoring, go to parent node */
		if(treeN[j].p == NULL) 
		{ 
			break; 
		} /* stop if root of T is reached */
		j = treePosN[(treeN[j].p)->n]; 
		goto nextWhileFW;
		
		nextWhileFW:
		continue;
  
	}
			#ifdef PRINTTABLE
						printf("Score array after forward recursion:\n");
						printf("\nS_kv:\n");
						gl_printS(S, kMax+1, nodesN);
			#endif			
			//~ printf("Score array after forward recursion:\n");
						//~ printf("\nS_kv:\n");
						//~ gl_printS(S, kMax+1, nodesN);			
		
  /* adjust gain and loss of complete gene family */
	int jP;
	for(j=0; j<nodesN; j++) 
	{
		if((treeN[j].p) == NULL) 
		continue;
		jP = treePosN[(treeN[j].p)->n];
		if(S[0][j] == 0 && S[0][jP] >= inf) 
		{
			treeN[j].lossFam = 1;
		} 	
		else 
		{
			treeN[j].lossFam = 0;
		}
		treeN[j].gainFam = 0;
	}

  /* start backtracing to label inner nodes */
  //  for(j=0; j<N_NODES; j++) done[j] = -1;
  //  printf("hello\n"); exit(0);  
  /* find optimal score for lca */
  i = 1;
  SS = inf;
  /* this variant results in first minimum! -> favours gain after lca */
  /* start from kMax to k=1 to favour loss after lca instead */
  //~ for(k=kMax; k>=0; k--) 
  //~ {
	for(k=0; k<=kMax; k++) 
	{ 
		if(S[k][lca_pos] <= SS) 
		{ 
			SS = S[k][lca_pos];
			//printf(" SS= %f lca_pos: %d\n",SS,lca_pos);
			
			i = k;
		}
	}
	//fprintf(stdout,"Optimal score at LCA=S[%d][%d]: %.2f (# of miRNAs: %d)\n\n",i,lca_pos,SS,i);
	/* i carries score of inner node */
	/* add gain of whole family */
	if(treeN[lca_pos].p == NULL && i == 0) { // members of this family are not in tree
		treeN[lca_pos].gainFam = 0;
	} else {
		treeN[lca_pos].gainFam = 1;
	}
	
	j = lca_pos;
	k = kMax*inf; //kk = kMax*inf;
	
	/**** start bw recursion */
	
	/* initalize lca and done list for lca */
	//done[j] = 1;
	treeN[j].m = i;
	treeN[j].seenBW = 1;
	treeN[j].gain = i;
	treeN[j].loss = 0;
	int seenNodes = 1;
  /* initialize all nodes before lca and their seenBW with 0 and 1 */
  /* for(c=0; c<nodesN; c++) { */
  /*   if(treeOrderN[c] == treeN[lca_pos].n) { break; } */
  /*   if(treeOrderN[c]  */
  /*   lab = treeOrderN[c]; */
  /*   pos = treePosN[lab]; */
  /*   treeN[pos].seenBW = 1; */
  /*   treeN[pos].m = 0; */
  /* } */


					//printf("mir at lca: %d\n",treeN[j].m); exit(0);
	for(c=0; c<nodesN; c++) 
	{
		//werte vor bwrecursion:
			    	//printf("p: %d l: %d m: %d seenBW: %d (%s)\n",c,treeN[c].n,treeN[c].m,treeN[c].seenBW,treeN[c].o);
		if(treeN[c].seenBW == 1) { seenNodes++; } /* every leaf and all nodes outside of tree rooted at lca need not to be visited */
	}

					// printf("dre: %d (%s) (seenBW: %d)\n",treeN[192].m,treeN[192].o,treeN[192].seenBW);     
	while(seenNodes <= nodesN) 
	{
		/* printf("j: %d i: %d (seenNodes: %d nodesN: %d) seenBW: %d\n",j,i, seenNodes, nodesN,treeN[j].seenBW); */
		/* gl_printTreeN(ga.outfile); */
	
		if(treeN[j].seenBW == 1) 
		{ 	/* already checked node */
			
			/* check if all children have been seen */
			for(c=0; c<treeN[j].nc; c++) 
			{
				cpos = treePosN[treeN[j].c[c]];
				if(treeN[cpos].seenBW == -1) 
				{
					j = cpos;
					i = treeN[j].m;
					goto nextWhile;
				}
			}
	
			/* all children scored, go to parent */
			if(j == lca_pos) { i = 0; }
			else { i = (treeN[j].p)->m; }
			if(treeN[j].p) 
			{
				//	j = (treeN[j].p)->pInT;
				j = treePosN[(treeN[j].p)->n];
				goto nextWhile;
			} else 
			{
				//printf("breaking\n");
				break;
			}
		
		} else 
		{ /* node to be checked */
	
			/* if not checked leave */
			if(treeN[j].nc == 0) 
			{
				i = treeN[j].m;
			} else 
			{	 
				/* since I select only the minimum (over k) of S[k][j] in fw recursion,
				I can do this here as well */
				/* yes, but include delta function!! */ // warum denn hier?
				minscore = inf;
				mink = kMax;
				for(k=0; k<=kMax; k++) 
				{
					if(S[k][j]  < minscore) 
					{
						minscore = S[k][j];
						
						mink = k;
					}
					
				}
				i = mink;
				
			}
			/* m is difference of parent.m (NULL at lca -> -i) and current i (number of miRNAs) */
			if(j == lca_pos || treeN[j].p == NULL) { m = -i; }
			else { m = (treeN[j].p)->m - i; }
		
			//printf("j:%d ist %s  k:%d m: %d\n",j,treeN[j].o,i,m); 
			/* adjust gain/loss */
			if(m < 0) { treeN[j].gain = abs(m); }
			else { treeN[j].loss = m; }
			
			//if(treeN[j].nc == 0) {
				treeN[j].m = i;
				
			
					/* gabor:
					* Problem: Dollo wird nicht immer beachtet, in manchen situationen, wo alle kinder leaves sind und keine gene haben wird
					* trotzdem am vaterknoten k!=0 gesetzt...
					* an der stelle ist der fehler schon passiert, man kann noch mal childs checken, 
					* wenn alle childs k=0 haben, muss auch aktueller (vater)knoten auf 0 gesetzt werden
					*/
					/*								lca	(label)
			
											treePosN[lca]  (Zeile)   
									treeN[			Zeile				].pInP  (pInP)
					treeOrderN	[				pInP++								] (label)		//hier kann man pInP erhöhen um in preorder abzulaufen
			treePosN[									label										] (Zeile)
			//~ */	
														
							//~ int testleave=0;
							//~ if (treeN[j].nc != 0 )	//knoten hat Kinder
							//~ {	//test ob alle kinder leaves sind UND diese leaves alle 0 gene haben
								
								//~ for (int i=0;i<treeN[j].nc;i++)
								//~ {	
									//~ if (treeN[treePosN[treeN[j].c[i]]].c[0]==(-1) && treeN[treePosN[treeN[j].c[i]]].m==0)		
									//~ {	//kind ist leave								&& leave hat keine mirs
										//~ testleave=testleave+1;
										//~ //printf("%s hat %d mirs\n",treeN[treePosN[treeN[j].c[i]]].o, treeN[treePosN[treeN[j].c[i]]].m );
										
									//~ }
								//~ }
								//~ if (testleave==treeN[j].nc)	//alle kinder sind leaves und all diese leaves haben 0 gene--> k kann auf 0 gesetzt werden
								//~ {
									//~ //printf("node der nur leaves als child hat, die KEINE mirs haben: %s\n",treeN[j].o);
									//~ treeN[j].m = 0;
									//~ treeN[j].gain = 0;
									
								//~ }   
								//~ testleave=0;
							//~ }	
							
			
				
								
			
		
		
																
			treeN[j].seenBW = 1;
			seenNodes++;
			goto nextWhile;
		}
    	nextWhile:
		continue;
	}
  
		/*ab hier sind alle k´ s gesetzt
		 *  
		 * */
		 printf("anzahl knoten: %d\n",nodesN);
			 
			//gl_printTreePS(lca, psfname);
			
	//abfrage notwendig, wenn -b==1 wird erst später ausgedruckt und verglichen
	if (ga.b==0)
	{	//float **pku=NULL;
		gl_printTreePSN(kMax, lca, ga);
		//fprintf(stdout, "The Parsimony Score tree is VISUALIZED in '%s'.\n", psfname);
		
	}
	//  gl_printTree(ga.outfile);
	if(ga.b==0 && ga.z==0)
	{
		gl_printTreeN(ga);
	}
	for(i=0; i<=kMax; i++) free(S[i]);
	free(S);
	
	if (ga.b==1)
	pf(kMax, ga, lca);
		
	return;

}

void pf_with_z(struct gl_arguments ga)
{
	int kMin, kMax, lca;
	
	lca = gl_getLCA(&kMin, &kMax);
	
	if(lca == -1) { printf("Your species are not in the tree\nExiting program.\n"); return; }
	
	//printf("LCA: %d\n", lca);
	
	
	pf(kMax, ga, lca);
	
}


void pf(int kMax, struct gl_arguments ga, int lca) //pf(int kMax) wird aufgerufen aus gl_calc bzw pf_with_z(int kMax)
{	
	//printf("KMAX: %d\n", kMax);
	
	/*  Z-Matrix reservieren */

	float **Z=NULL; 
	  Z = (float **) calloc(kMax+1, sizeof(float *));				// calloc für arrays, erstes argument gibt array länge an (kMax+1 ist länge i in matrix))
	for(int i=0; i<=kMax; i++) 
	{		//nodesN+1 ensprechen spalten, also spezies
			//in jeder Zeile Z[i] muss spalte der länge nodesN reserviert werden
		Z[i] = (float *) calloc(nodesN+1, sizeof(float));  
		for(int j=0; j<nodesN; j++) 
		{
			Z[i][j] = -1.; /* preorder does not matter jet */
		}
	}
  
  /* initialize leaves */
  
  
	for(int j=0; j<nodesN; j++) 
	{
		if(treeN[j].nc == 0) 
		{ /* this is a leaf node for it has no children */
			
			for(int i=0; i<=kMax; i++) 
			{ /* check the number of miRNAs at this leaf */
				if(treeN[j].m == i) Z[i][j] = 1.;
				else Z[i][j] = 0;
			}
		} 
	}
		
		#ifdef PRINTTABLE	
			printf("\nZ array after initialization (before forward recursion): \n"); 
			printf("\nZ_kl:\n");
			gl_printS(Z, kMax+1, nodesN);
		#endif		
	
	////////////////forward recursion///////////////////////////
	
	
	int bet=BETA;
	//int zkM = zkMax;
	int found;
	labeltwo:
	found = 0; //bleibt 0, wenn alle kinder von v gescored sind
	//begin:
	for (int i=0; i<nodesN; i++ ) //durchsuche Z nach ersten eintrag mit -1 (die zu scorenden vater knoten)
	{	labelone:
		//printf("i:%d\n",i);
		//printf("1. kind von TreeN[%d]: %d\n", i, treePosN[treeN[i].c[0]]);
		if (Z[0][i]==-1.) 
		{	if (treeN[i].pInP < treeN[treePosN[lca]].pInP || treeN[i].pInP > LAST) //wenn sich ein vaterknoten außerhalb des LCA befindet wird er nicht berechnet, 1. splte auf 0, rest auf 1
			{
				//printf("der vater %s and preordpos %d wird ausgelassen\n lca:%d last:%d\n", treeN[i].o, treeN[i].pInP, treeN[treePosN[lca]].pInP,LAST);
				for (int n=0;n<=kMax;n++)
				{
					if (n==0)
					{
						Z[n][i]=1;
					}
					else
					Z[n][i]=0;	//da sich vaterknoten außerhalb das LCA befindet, keine berechnung, alles null setzen
				}
			continue;	
			}
			int childn=treeN[i].nc;  //überprüfen ob alle kinder gescored sind
			int y=0;
			do
			{	//printf("kinder von TreeN[%d]:%d. %d\n", i, y+1, treePosN[treeN[i].c[y]]);
				if(Z[0][treePosN[treeN[i].c[y]]] == -1.  ) //wenn kind noch nicht gescored wurde (-1), zum nächsten vater gehen
				{	found = 1;	
					//printf("überspringen bei i=%d\n", i);							//schleife wird verlassen, wenn ein kind -1 hat
					i++;
					goto labelone;
				}
				y++;	
			}while (y<childn);	
			
			
			////an der stelle sollte skaliert werden denn, v noch nicht berechnet, aber alle seine kinder, diese scores werden jetzt erst angepasst und dann wird v ausgerechnet
			//wenn u von v leave nix machen, sonst anpassen
			//sobald ein u von v ein leave ist, werden die scores von u nicht skaliert



			y=0;
			int leave=0;
			do
			{
				if (treeN[treePosN[treeN[i].c[y]]].c[0]== -1.)
				{
					leave=1;
					break;
				}
				//printf("vater ist: %s und dessen leaves sind: %s\n", treeN[i].o, treeN[treePosN[treeN[i].c[y]]].o);
				y++;
			}while (y<childn);	
			
			if (leave == 0) //alle knoten u von v sind innere knoten und müssen skaliert werden		
			{	float tempscoresum =0.;	//summe der scores eines u's
				float tempscoretotal =0.; //summe der scores aller u's eines v's
				for (y=0; y<childn; y++) //summe der scores aller u von v
				{
					for (int zeile = 0; zeile <=kMax; zeile++)
					{
						//printf("alter score bei Z[%d][%d]= %f\n", zeile, treePosN[treeN[i].c[y]], Z[zeile][treePosN[treeN[i].c[y]]]);
						tempscoresum= tempscoresum + 	Z[zeile][treePosN[treeN[i].c[y]]];
					}
					tempscoretotal = tempscoretotal + tempscoresum;	
				}
				//hier sind alle scores aller u's eines v's ausummiert, setzen der skalierten scores in Z:
				for (y=0; y<childn; y++) 
				{
					for (int zeile = 0; zeile <=kMax; zeile++)
					{
						Z[zeile][treePosN[treeN[i].c[y]]] = Z[zeile][treePosN[treeN[i].c[y]]] / tempscoretotal;
						//printf("neuer score bei Z[%d][%d]= %f\n", zeile, treePosN[treeN[i].c[y]], Z[zeile][treePosN[treeN[i].c[y]]]);
					}	
				}
				//printf("wir sind hier bei: %s\n", treeN[i].o);
			}	
				
			//printf(" kind treePosN von %s ist %d\n ", treeN[i].o, treeN[treePosN[treeN[i].c[0]]].n);
			//printf("was haben wir bis jetzt bei i = %d:\n\n",i);
			//gl_printS(Z, kMax+1, nodesN); 
			//für den vater v bei i wurden alle kinder gescored, der score der kinder muss jetzt aufsummiert werden, dann ist gescalte score = score alt/ tempsumscale
			//
			
		//
			float tempZ=1.;
			float tempSum=0.;
			int k;
			for (k=0; k<=treeN[i].P; k++) //für ein k werden kombinationen aller kinder berechnet 
			{
				for (int j=0; j<childn ;j++)		 
				{	int kchild=0;
					while (kchild<=treeN[i].P)
					{					//hier Z[kchild][j(kindspalte)]
						
						tempSum= tempSum+Z[kchild][treePosN[treeN[i].c[j]]]*exp((bet*abs(k-kchild)));   //version ohne delta funktion
						
						//tempSum= tempSum+Z[kchild][treePosN[treeN[i].c[j]]]*exp(bet*gl_delta(k,kchild,treePosN[treeN[i].c[j]]));
						
						
						//printf("Z[kchild][treePosN[treeN[%d].c[%d]]] = Z[%d][%d]%f\n",kchild,treePosN[treeN[i].c[j]] ,kchild,treePosN[treeN[i].c[j]], Z[kchild][treePosN[treeN[i].c[j]]]);
						//printf("i:%d\n", i);
						//printf("treeN[%d].P:%d\n",i,treeN[i].P);
						kchild++;
					}
					//hier ist eine kindspalte fertig aufsummiert
					tempZ=tempZ*tempSum;
					//printf("tempZ:%f\ttempSum:%f\n",tempZ,tempSum);
					tempSum=0.;
				//	printf("child in treeN pos: %d\n", treePosN[treeN[i].c[j]]);
				
				}
				//printf("tempZ in Z[%d][%d]:%f\t\n",k,i,tempZ);
				//double tempZlog=log(1/tempZ-1);
				Z[k][i]=tempZ;		
				
		//hier test auf nans	
		//printf("Z[%d][%d]=%f\n",k,i,Z[k][i]);	
		if (isnan(Z[k][i]))
		{
			printf(" this is a nan:%f\n",Z[k][i]);
			Z[k][i] = 0.;
			printf(" now it's:%f\n", Z[k][i]);
		}
				
				//ACHTUNG FALSCH FALSCH FALSCH !!!! SO NICHT:
				//um dollo constraint zu beachten: müssen alle werte für k=0 an v auf null gesetzt werden, wenn tempZ nicht 1 ist, das bedeutet, wenn nicht alle kinder k=0 belegt haben, darf bei v k nicht 0 gewählt werden können ( if (P_v != 1)  Z_0v=0 
				//nee moment das ist doof, da scores bei k=0 mit einbezogen werden müssen für berechnung aller scores bei v...
				
	/*			if(k==0)
				{
					//test: finde alle werte, die nicht 1 sind in k=0 zeile
					if (Z[k][i]!= 1 )
					{
					 //printf("0 setzen bei v=%s\n ",treeN[i].o); //passt
					 Z[k][i]=0.;
					}
				}
				//so würde man zu wenig loss berechnen, dollo constraint wird erst beim setzen des k's berücksichtig, ohne die scores zu beeinflussen
	*/			
				
				
										//wenn Pm kleiner als kmax bleiben ab Pm -1 stehen:
				if(treeN[i].P<kMax)	//werden mit 0 gefüllt
					for(int m =treeN[i].P+1;m<=kMax;m++)
						Z[m][i]=0.;
				
				tempZ=1.;
				//printf("wdwdwdtreeN[%d].o: %s\n", i,treeN[i].o);
			}
		}	
	}
	if (found==1)	//sollte vaterknoten gefunden sein, welcher noch nicht bearbeitet wurde, suche von vorn
	goto labeltwo;	//nochmals durchlauf aller inneren knoten
	
	 
	#ifdef PRINTTABLE
		 printf("\nZ array after FR :\n"); 
		 printf("\nZ_kv:\n");
		 gl_printS(Z, kMax+1, nodesN);
	#endif
	
	
	
		pfb(kMax,nodesN, Z, ga, lca);
		

	
	
	
}


void pfb(int kMax,int nodesN, float **Z, struct gl_arguments ga, int lca )  
{
						/*	for (int i=1; i<nodesN;i++)
							{   int postreeN= treePosN[treeOrderN[i]];
								printf("ausgabe nach preorder %s davon der Vater: %s\n", treeN[postreeN].o, treeN[treePosN[treeN[postreeN].parentLab]].o);
							}
						*/
	//neue matrix für Z*_ku 
		//gl_printS(Z,kMax+1, nodesN);	
	
	float **Zaussen=NULL; 
	Zaussen = (float **) calloc(kMax+1, sizeof(float *));			
	for(int i=0; i<=kMax; i++) 
	{		
		Zaussen[i] = (float *) calloc(nodesN+1, sizeof(float));  
		for(int j=0; j<nodesN; j++) 
		{
			Zaussen[i][j] = -1.;
		}
	}
	
	
	//Außenmatrix einmalig bei LCA!!!! belegen: , wenn man schon durch schleife geht, kann gleich summe der spalte ausgerechnet werden:
	float tempsumroot=0;
	for (int k=0; k<=kMax;k++)
	{	//printf("Z[k][treePosN[lca]]  %f  \n",Z[k][treePosN[lca]]);
		Zaussen[k][treePosN[lca]]=Z[k][treePosN[lca]];			//belegen der LCA-Spalte
		tempsumroot=tempsumroot + Zaussen[k][treePosN[lca]];	//summe der LCA festhalten
	}
	//printf("tempsumroot für LCA:%f\n",tempsumroot);
	//P an LCA kann direkt ausgerechnet werden:
	//P- Matrix wird nicht mehr benötigt, skalierte Z^{*S} entspricht den p's
	for (int i=0; i<=kMax; i++)
	{
		//printf("Zaussen lca bei k vor skalierung: %f\n",Zaussen[i][treePosN[lca]]);
		Zaussen[i][treePosN[lca]]=Zaussen[i][treePosN[lca]]/tempsumroot;	//p für jedes k in lca-spalte eintragen
		//printf("Zaussen lca bei k nach skalierung: %f",Zaussen[i][treePosN[lca]]);																						//ist gleichzeitig skalierung des lca node
	}	
	
	
		
	/*		
	///////////////////////////////////////////////////////////////////////		
					Änderung gabor: 
			  		Z^*_w=\prod_{s}^{}\sum_{k}^{}Z_{ks}\cdot\sum_{i}^{} Z^*_{iv}\cdot e^{-\delta(i,k)}
	*/
	/*										lca	(label)
	    
										treePosN[lca]  (Zeile)   
								 treeN[			Zeile				].pInP  (pInP)
				 treeOrderN	[				pInP++								] (label)		//hier kann man pInP erhöhen um in preorder abzulaufen
		treePosN[									label										] (Zeile)
	*/	
	
	
/*	
	printf("nodes in preorder:\n");
	for (int u= treeN[treePosN[lca]].pInP; u<=LAST;u++)		
	{	
		printf("%s :pInP:%d\n",treeN[treePosN[treeOrderN[u]]].o,treeN[treePosN[treeOrderN[u]]].pInP );
	}
*/	
	int bet=BETA;
			
	for (int u= treeN[treePosN[lca]].pInP+1; u<=LAST;u++)			//lca+1 starten um lca auszulassen, ablaufen in perorder bis node LAST in perorder erreicht ist 
	{	//printf("pInP:%d , ist: %s\n",u, treeN[treePosN[treeOrderN[u]]].o);
			//wichtig, bei den sib loops immer auf childarray des vaters zugreifen das u bleibt ab hier fixiert
		struct gl_nodeN *parent_p;
		parent_p=&treeN[treePosN[treeN[treePosN[treeOrderN[u]]].parentLab]]; //treeN struct des vaters des aktuellen u´ s
		int parent_zeile=treePosN[treeN[treePosN[treeOrderN[u]]].parentLab]; //Zeile vater, wird in matrix angegeben um seine außen zu erhalten
		for (int sib_u1=0;sib_u1<=parent_p->nc; sib_u1++)	//knoten wählen aus childs des vaters
		{	
						
			if (treeN[treePosN[parent_p->c[sib_u1]]].c[0]==-1 || treeN[treePosN[parent_p->c[sib_u1]]].seenBWp == 1)
			{		//wenn knoten leave oder bereits gescored wurde, weiter
				//printf("weiter bei: %s\n",treeN[treePosN[parent_p->c[sib_u1]]].o);
				continue;
			}
			
			else	
			{	//printf("berechnet wird score für: %s\n",treeN[treePosN[parent_p->c[sib_u1]]].o);
				
				float tempsum =0.;		//lt. recursion: summe über i
				float tempsumk =0.;		//lt. recursion innen mal summe für ein k
				float tempsumk_hold=0.;	//lt. recursion: summe über alle k´s
				float temp_prod_sib=1.; //lt. recursion: hält prod über s
					
				for (int sib_u=0;sib_u<parent_p->nc; sib_u++) //ablaufen der kindknoten des vaters parent_p
				{
					if (treeN[treePosN[parent_p->c[sib_u]]].n != treeN[treePosN[parent_p->c[sib_u1]]].n) //auslassen, wenn zu scorender knoten = ein geschwister
					{	//printf("geschwister von %s die beachtet werden: %s\n",treeN[treePosN[parent_p->c[sib_u1]]].o,treeN[treePosN[parent_p->c[sib_u]]].o);		//passt
						for (int k=0; k<=kMax; k++)
						{	
							for(int s=0; s<=kMax;s++)								
							{	//an stelle k bei u werden mit s alle wege ausprobiert(vom vater zu u) für alle geschister sib_u von sib_u1
								
							//printf("sib_u1=%d: %s aufsummieren der geschw.:sib_u=%d %s\n",sib_u1,treeN[treePosN[parent_p->c[sib_u1]]].o,sib_u,treeN[treePosN[parent_p->c[sib_u]]].o);
								
							//printf("tempsum: %f + Zaussen[%d][%d] (%f)=\n",tempsum,s,*parent_p, Zaussen[s][*parent_p]);	
														
							tempsum= tempsum + Zaussen[s][parent_zeile]*exp((bet*abs(k-s))); //wege von v nach u(k)      //version ohne delta fkt
							//printf("Zaussen: %f * exp((%d * abs(%d - %d)  ((%d)) =  %f\n",Zaussen[s][*parent_p], bet, k,s,abs(k-s), exp((bet*abs(k-s))));
							//printf(" tempsum = %f\n",tempsum);
							//printf("dad species: %s\n",parent_p->o);
							
							}
							tempsumk=Z[k][treePosN[parent_p->c[sib_u]]]*tempsum;	//innen mal summe der wege eines k´s vom vater zu allen k´s von u
							tempsumk_hold=tempsumk_hold+tempsumk;			
							
							//printf("innen %f * tempsum %f = tempsumk %f\n",Z[k][treePosN[parent_p->c[sib_u]]],tempsum, tempsumk);
							//printf("tempsumk: %f\n",tempsumk);
							//printf("tempsumk_hold: %f\n",tempsumk_hold);
						}
						
						//hier ist ein geschwisterknoten komplett aufgerechnet -> multiplikation aller geschwister 
						//printf(" temp_prod_sib %f * tempsumk_hold  %f", temp_prod_sib,tempsumk_hold);
						temp_prod_sib = temp_prod_sib*tempsumk_hold;
						//printf("=temp_prod_sib %f\n",temp_prod_sib);
					
					}
							
				}	//prod aller siblings fertig	--> gesp. in temp_prod_sib		das entsprich Z^*_w in recursion
				tempsum=0;
				//weiter mit "normaler recursion um u auszurechnen", hier kann prod aller siblings einbezogen werden (temp_prod_sib)
				for (int k=0; k<=kMax; k++)
					{	
						for(int s=0; s<=kMax;s++)								
						{	//an stelle k bei u werden mit s alle wege ausprobiert
						
						tempsum= tempsum + Zaussen[s][parent_zeile]*exp((bet*abs(k-s))); //wege von v nach u(k)      //version ohne delta fkt
						//printf("tempsum %f, bei s=%d  und k=%d für %s\n", tempsum,s,k,treeN[treePosN[parent_p->c[sib_u1]]].o);
						}	
							//hier sib_u1 beachten, das ist der knoten um den es gerade geht
					
						Zaussen[k][treePosN[parent_p->c[sib_u1]]] = Z[k][treePosN[parent_p->c[sib_u1]]]*tempsum*temp_prod_sib; 
						//printf(" innen    %f  mal tempsum   %f   mal tempprod  %f  =  außen %f  bei k=%d   für  %s \n",Z[k][treePosN[parent_p->c[sib_u1]]], tempsum, temp_prod_sib,Zaussen[k][treePosN[parent_p->c[sib_u1]]],k,treeN[treePosN[parent_p->c[sib_u1]]].o ); 
						tempsum=0;
					}
					//printf("sind gerade bei: %s\n",treeN[treePosN[parent_p->c[sib_u1]]].o);	
					treeN[treePosN[parent_p->c[sib_u1]]].seenBWp = 1;
				
			}
			
		}
		for (int sib_skalieren=0;sib_skalieren<=parent_p->nc; sib_skalieren++)	//knoten wählen aus childs des vaters
		{
			//skalieren für alle kinder von parent_p sib_u1
			float tempsump=0.;
			for (int pk=0; pk<=kMax;pk++)
			{
				tempsump=tempsump + Zaussen[pk][treePosN[parent_p->c[sib_skalieren]]];	
				
			}
			//jetzt steht erst mal summe je spalte, für jedes k wahrscheinlichkeit ausrechnen:
			for (int pk=0; pk<=kMax;pk++)
			{
				Zaussen[pk][treePosN[parent_p->c[sib_skalieren]]]= Zaussen[pk][treePosN[parent_p->c[sib_skalieren]]]/tempsump;	
			}
		}
	}
		
	#ifdef PRINTTABLE
	printf("\nZ_aussen array after calculation: \n"); 
	gl_printS(Zaussen, kMax+1, nodesN);
	#endif
	
	printf("\nZ_aussen array after calculation: \n"); 
	gl_printS(Zaussen, kMax+1, nodesN);
	
	//funktionsaufruf für gain/loss nach fertiger PF:
	pfgl(Zaussen, kMax, nodesN, lca, ga);
	
	
	for(int i=0; i<=kMax; i++) free(Z[i]);
  free(Z);
  
  for(int i=0; i<=kMax; i++) free(Zaussen[i]);
  free(Zaussen);
  
//	for (int i=0; i<nodesN;i++)
//	printf("mir von %s ist %d\n",treeN[i].o ,treeN[i].pfm );
  
  
	
}

void pfgl(float **Pku,int kMax,int nodesN, int lca, struct gl_arguments ga )	// partitionfunctionGainLoss
{														//was ist mit gainFam, lossFam???
	//für jedes k einer spalte (aus matrix Pku))muss diff(k) ausgerechnet werden 
	// für jede spalte wird array angelegt, länge k+1, einmaliges initialisieren, kleinstes diff ermitteln, k in TreeN[i].pfm  eintragen
	// wenn alle k eingetragen, diff der k zw u und v jeweils in TreeN[i].pfGain und TreeN[i].pfLoss eintragen
	//ALSO pfGain(u) und pfLoss(u) immer in bezug auf v --> v.m-u.m>0, dann ergebnis in pfLoss(u) (von v nach u verlust von k)
	//wenn v.m-u.m<0, dann betrag des ergebnisses in pfGain(u) (von v nach u gewinn von k)
	//in buildtree() initialisierung von pfGain bzw pfLoss mit 0, ebenso pfm und pfGainFam und pfLossFam
	
	//ACHTUNG um dollo constraint zu beachten, bei beiden methoden treeM[xy].P berücksichtigen!!! wenn .P != 0 darf .pfm nicht 0 gesetzt werden!!!(jedesmal abfragen!)
	
	
	if (ga.P==0)	//k is set to the smallest Diff value
	{
		float *arrayDiff=NULL;
		arrayDiff= (float *) calloc(kMax+1, sizeof(float*));	//kMax+1,z.B. wenn k=1, dann zeile für k=0 und k=1
	
		for (int u=0; u<nodesN;u++)		//ablaufen der spalten aus Pku
		{	if (treeN[u].c[0] == -1 || treeN[u].pInP< treeN[treePosN[lca]].pInP || treeN[u].pInP> treeN[treePosN[treeOrderN[LAST]]].pInP )	//leaves überspringen, diff für alle inneren knoten, reihenfolge egal || und überspringen aller nodes außerhalb der LCA und LAST!!!
			continue;
			else
			{
				for (int i=0;i<=kMax;i++)	//ablaufen jeder zeile einer spalte
				{	
					if (i==0)				//fall 1: i==0 -> diff ist summe aller pku einträge ab i+1 bis i=kMax
					{	
						float difftemp=0;
						for (int n=0; n<kMax;n++)
						{					
							difftemp=difftemp+Pku[n+1][u];	//geht los bei 1, endet bei kMax, da[n+1] und n<kMax
							
						}	
						arrayDiff[i]=difftemp;
												
						//printf("u=%d\n",u );
						//printf("arrayDiff[%d]= %f\n",i ,arrayDiff[i] );
					}
					
					else if (i==kMax)		//fall 2: i==kMax -> diff ist summe von i=0 bis i=kMax-1
					{
						float difftemp=0;
						for (int n=0; n<kMax;n++)
						{
							difftemp=difftemp+Pku[n][u];
						}
						arrayDiff[i]=difftemp;
						
						//printf("arrayDiff[%d]= %f\n",i ,arrayDiff[i] );
					}
					
					else 					//fall 3: alles dazwischen
					{
						float diffneg=0;	// summe der einträge für n<i
						float diffpos=0;	// summe der einträge für n>i
						for (int n=0; n<=kMax;n++)
						{	
							if (n==i)
							continue;
							else if (n<i)
							{
								diffneg=diffneg+Pku[n][u];
							}
							else 
							{
								diffpos=diffpos+Pku[n][u];
							}
						}
						//printf("diffneg bei u=%d uns i=%d = %f\n", u,i,diffneg);
						//printf("diffpos bei u=%d uns i=%d = %f\n", u,i,diffpos);
						arrayDiff[i]=fabs(diffneg-diffpos);	
						
						//printf("arrayDiff[%d]=%f\n",i,arrayDiff[i]);	
					}
				}
				//hier steht jetzt arraydiff[] von u, jetzt muss i von min diff gesetzt werden an treeN[u].pfm			
				float tempmindiff=1.;
				int tempi=0;
				for (int i=0;i<=kMax;i++)
				{
					if (tempmindiff>=arrayDiff[i])
					{
						tempmindiff=arrayDiff[i];
						tempi=i;
					}
					else
					continue;
				}
				//printf("minDiff ist %f bei k=%d\n",tempmindiff ,tempi );
				//hier dollo constraint beachten: sobald ein k auf 0 gesetzt ist können alle kinder darunter und deren kinder auch gleich 0 gesetzt werden --> Pku an k=0 auf 1 an allen kindern...
				
				treeN[u].pfm=tempi;	
				//DOLLO CONSTRAINT:
				if (tempi==0 && treeN[u].P!=0 )
				{
					treeN[u].pfm=1;
					
					
					
				}
				if (ga.b == 1) //ACHTUNG: hier setzen des zusätzlichen spalteneintrags für .out Datei wenn -b , auch für PS:
					//nodes außerhalb LCA und LAST dürften hier gar nicht auftreten:
					{
						treeN[u].pfP=Pku[treeN[u].pfm][u];
						treeN[u].psP=Pku[treeN[u].m][u];
						//printf("fuer %s ist pfP=%f bei genes=%d\n",treeN[u].o ,treeN[u].pfP, treeN[u].pfm);
						//printf("fuer %s ist psP=%f bei genes=%d\n",treeN[u].o ,treeN[u].psP, treeN[u].m);
					}
					
				if (ga.z == 1) //ACHTUNG: hier setzen des zusätzlichen spalteneintrags für .out Datei wenn -z 
					//nodes außerhalb LCA und LAST dürften hier gar nicht auftreten:
					{
						treeN[u].pfP=Pku[treeN[u].pfm][u];
						//printf("fuer %s ist pfP=%f bei genes=%d\n",treeN[u].o ,treeN[u].pfP, treeN[u].pfm);
					}
			}
		}
		free(arrayDiff);
	}
	
	if (ga.P==1)	//k is set to the maximum p
	{
		for (int u=0; u<nodesN;u++)		//ablaufen der spalten aus Pku
		{	if (treeN[u].c[0] == -1 || treeN[u].pInP< treeN[treePosN[lca]].pInP || treeN[u].pInP> treeN[treePosN[treeOrderN[LAST]]].pInP )	//leaves überspringen, diff für alle inneren knoten, reihenfolge egal || und überspringen aller nodes außerhalb der LCA und LAST!!!
			continue;
			else
			{
				int tempi =0;
				float tempP=0.;
				if (treeN[u].c[0] == -1)	//leaves überspringen, diff für alle inneren knoten, reihenfolge egal
				continue;
				else
				{
					for (int i=0;i<=kMax;i++)	//ablaufen jeder zeile einer spalte
					{
						if (Pku[i][u]>=tempP)
						{
							tempP=Pku[i][u];
							tempi=i;
						}	
					}
					treeN[u].pfm=tempi;	//größter P in tempP gespeichert, setzen des k's an stelle des gefundenen tempP
					//DOLLO CONSTRAINT:
					if (tempi==0 && treeN[u].P!=0)
					{
						treeN[u].pfm=1;
					}
					if (ga.b == 1)	//ACHTUNG: hier setzen des zusätzlichen spalteneintrags für .out Datei wenn -b , auch für PS:
					{
						treeN[u].pfP=Pku[treeN[u].pfm][u];
						treeN[u].psP=Pku[treeN[u].m][u];
						//printf("fuer %s ist pfP=%f bei genes=%d\n",treeN[u].o ,treeN[u].pfP, treeN[u].pfm);
						//printf("fuer %s ist psP=%f bei genes=%d\n",treeN[u].o ,treeN[u].psP, treeN[u].m);
					}
					if (ga.z == 1) //ACHTUNG: hier setzen des zusätzlichen spalteneintrags für .out Datei wenn -z 
					//nodes außerhalb LCA und LAST dürften hier gar nicht auftreten:
					{
						treeN[u].pfP=Pku[treeN[u].pfm][u];
						//printf("fuer %s ist pfP=%f bei genes=%d\n",treeN[u].o ,treeN[u].pfP, treeN[u].pfm);
					}
				}
			}		
		}		
	}

					/*	for (int u=0; u<nodesN;u++)	
						{
							printf("pfm: %d, bei %s und m:%d\n",  treeN[u].pfm, treeN[u].o, treeN[u].m);
						}
					*/
	//da hier alle k's gesetzt wurden, können pfLossFam gesetzt werden, wenn k=0 bei einem u, dessen v k!=0 ist, pfLossFam an u=1, aber ACHTUNG, bei leaves dürfen nur m betrachtet werden, nicht pfm, da diese auf 0 gesetzt bleibt (initialer wert, der nicht geändert wird, da leaves bei pfgl übersprungen werden!!!) 
	for (int u=0; u<nodesN;u++)	
	{	
		if (treeN[u].n==lca) //HIER pfGain setzen!!! an LCA!!! ist gain= mirs (bzw. genes) erstes auftreten der gene
		{
			treeN[u].pfGain=treeN[u].pfm;
			continue;
		}
		//printf("der Vater von %s ist %s an einleseposition: %d und par: %d \n", treeN[u].o, treeN[treePosN[treeN[u].parentLab]].o,treePosN[treeN[u].parentLab], treeN[u].parentLab );
		if (treeN[u].pfm==0 && treeN[treePosN[treeN[u].parentLab]].pfm!=0 && treeN[u].c[0] != -1) //vgl d. k's zw v und u wenn u kein leave(an u wird pfm verwedent)
		{
			//printf("pfFamLoss bei %s  k bei Vater:%d k bei u:%d \n", treeN[u].o, treeN[treePosN[treeN[u].parentLab]].pfm ,treeN[u].pfm );
			treeN[u].pfLossFam=1;
		}
		else if (treeN[u].m==0 && treeN[treePosN[treeN[u].parentLab]].pfm!=0 && treeN[u].c[0] == -1) //vgl d. k's zw v und u wenn u leave ist (an u wird m verwedent)
		{	
			//printf("pfFamLoss bei %s  k bei Vater:%d k bei u:%d \n", treeN[u].o, treeN[treePosN[treeN[u].parentLab]].pfm ,treeN[u].m );
			treeN[u].pfLossFam=1;
		}
	}
	//test ob alles richtig eingetragen wurde:
	//for (int u=0; u<nodesN; u++)
	//printf("treeN[%d].pfm=%d\n",u , treeN[u].pfm );	
	
	/* belegen von pfGain, pfLoss: */
	//pfm steht jetzt, für alle u von v müssen jetzt pfGain und pfLoss gesetzt werden (s.o.:
	// ALSO pfGain(u) und pfLoss(u) immer in bezug auf v --> v.m-u.m>0, dann ergebnis in pfLoss(u) (von v nach u verlust von k)
	//wenn v.m-u.m<0, dann betrag des ergebnisses in pfGain(u) (von v nach u gewinn von k) )
	
	//reihenfolge ist egal, an jedem knoten steht lab des vaters, ABER: an LCA muss auf ast zu LCA gain=pfm gesetzt werden! dort ist gainFAMpf=1 

	
	int lca_pos = treePosN[lca];		//lca_pos ist die einlesepos. (zeilennummer) des LCA	
	//printf("LCA:%d\n", lca);
	//printf("LCApos:%d\n", lca_pos);
	
	//ACHTUNG!!! vorm nächsten schritt müssen treeN[u].m aller leaves übernommen werden (gegebene daten, die sich nicht ändern) für treeN[u].pfm (nur bei leaves!) 
	for (int u=0; u<nodesN;u++)
	{
		if (treeN[u].c[0]==-1) //leave gefunden
		{
			treeN[u].pfm=treeN[u].m;
			//printf("das ist ein leave : %s\n" , treeN[u].o);
			//printf("%s hat %d mirs\n",treeN[u].o,treeN[u].pfm); 
		}
	} 

					
	for (int u=0; u<nodesN;u++)
	{
		if (u==lca_pos)
		{	
			//ACHTUNG: beim LCA einer Familie wird gainFAM auf 1 gesetzt, sollte sich nicht von parsimony Score unterscheiden
			//printf("LCA wird ausgelassen: %s\n", treeN[u].o);
			treeN[u].pfGainFam=1;		
			if (treeN[lca_pos].c[0]== -1)		//wenn LCA ein leave ist, muss pfGain auf pfm gesetzt werden!!!:
			{
				treeN[lca_pos].pfGain=treeN[lca_pos].pfm;
			}	
			continue;
		}
		//else if (treeN[u].pInP< treeN[treePosN[lca]].pInP || treeN[u].pInP> treeN[treePosN[treeOrderN[LAST]]].pInP)  //gain und loss muss auch für knoten außerhalb der gefundenen paralogen berechnent werden!!!
		//continue;
		else 	//ALSO pfGain(u) und pfLoss(u) immer in bezug auf v --> v.m-u.m>0, dann ergebnis in pfLoss(u) (von v nach u verlust von k)
		{		//wenn v.pfm-u.pfm<0, dann betrag des ergebnisses in pfGain(u) (von v nach u gewinn von k)
			if (treeN[treePosN[treeN[u].parentLab]].pfm-treeN[u].pfm > 0)	 //hier muss man noch mal checken ob treeN[u] root ist, dann darf kein loss gesetzt werden
			{	//printf("die root ist %s\n", treeN[treePosN[treeOrderN[0]]].o);
				if(treeN[u].pInP != 0)
				{
					//printf("loss bei %s:\n gerechnet wird m von vater %d,  (%s)  - m child: %d (%s) \n\n ", treeN[u].o, treeN[treePosN[treeN[u].parentLab]].pfm , treeN[treePosN[treeN[u].parentLab]].o , treeN[u].pfm, treeN[u].o);
					treeN[u].pfLoss=treeN[treePosN[treeN[u].parentLab]].pfm-treeN[u].pfm; //abs nicht nötig, ist ja eh positiv
					//sonst nichts zu tun, pfGain bleibt auf 0 gesetzt
				}
			}	
			else if (treeN[treePosN[treeN[u].parentLab]].pfm-treeN[u].pfm < 0)  //in if bedingung kommt man nicht rein(?), so wird pfGain gesetzt auch wenn v.pfm-u.pfm=0
			{	//AUFPASSEN an root wurde pfGain bereits gesetzt. darf nicht überschrieben werden!
				if (treeN[u].parentLab==-1)
				continue;
				//printf("gain bei %s:\n", treeN[u].o);
				treeN[u].pfGain=abs(treeN[treePosN[treeN[u].parentLab]].pfm-treeN[u].pfm);
			}
			//else
			//printf("keine änderung bei %s:\n", treeN[u].o);
		}
	}
/*	
	for (int u=0; u<nodesN;u++)
	{
		printf("%s hat P: %d\n", treeN[u].o, treeN[u].P);
	}
*/

//~ for (int u=0; u<nodesN;u++)
//~ {
	//~ //printf("%s hat %d gains\n", treeN[u].o, treeN[u].pfGain);
	//~ //printf("%s hat %d loss\n", treeN[u].o, treeN[u].pfLoss);
	//~ //printf("parent von %s (lab:%d) (mirs:%d) ist %s\n",treeN[u].o, treeN[u].n,treeN[u].pfm ,treeN[treePosN[treeN[u].parentLab]].o);
	//~ printf("%s hat %d PSmirs und \t%f .psP\n", treeN[u].o, treeN[u].m, treeN[u].psP );
	
	//~ printf("%s hat %d PFmirs und \t%f .pfP\n", treeN[u].o, treeN[u].pfm, treeN[u].pfP);
//~ }
 
		if(ga.z==1 || ga.b==1)
		{	
			gl_printTreePSN(kMax, lca, ga);
			gl_printTreeN(ga);
		}
		
	//printf("P_ku array after calculation: \n"); 
	//	gl_printS(Pku, kMax+1, nodesN);
	
}

/* ---------------------------------------------t----------------------------- */

/* delta function that returns corresponding weight */
/*
 * 
 * name: gl_delta
 * @param int i, int k, int node
 * @return float
 * 
 */

float gl_delta(int i, int k, int node) 
{
	float inf=9999.; 
	float x, delta=0.;

  //if(i==0 && k == 0) { delta = inf; }
  //else {
    //~ x = fabs((float)(k - i));			
    //~ x = (float)(k - i);
    //~ if(x == 0) { delta = 0.; }
    //~ else if(x < 0) 
    //~ {
      //~ delta = fabs(x)*treeGainWeight[node];
    //~ } 
    //~ else 
    //~ {
      //~ delta = x*(treeLossWeight[node]);		//wozu diese 1??? ACHTUNG!!: delta = x*(1+treeLossWeight[node]);	im original!?!?!
    //~ }
  
				//~ //printf("delta = %f (%s) gainW: %f lossW: %f\n",delta,treeN[node].o,treeGainWeight[node],treeLossWeight[node]);
  //~ return delta;
  
		if(i==0 && k == 0) { delta = inf; }
		else 
		{
			x = (k - i);	
		}			
			
		if(x == 0) { delta = 0.; }
		else if(x < 0) 
		{
			delta = abs(x)*treeGainWeight[node];
		} 
		else 
		{
			delta = x*(treeLossWeight[node]);		//wozu diese 1??? ACHTUNG!!: delta = x*(1+treeLossWeight[node]);	im original!?!?!
		}
  
				//printf("delta = %f (%s) gainW: %f lossW: %f\n",delta,treeN[node].o,treeGainWeight[node],treeLossWeight[node]);
	return delta;
	
  
}

/* -------------------------------------------------------------------------- */

/* find index of node that is last common ancestor of current miRNA family
   assign min and max number of miRNAs in current alignment 
   returns label (!)
 */
int gl_getLCA(int *kMin, int *kMax) 
{

	int i, j, lab=0, pos=0, labFirst=-1, labLast=-1, lca=0;
	int posFirst=-1,posLast=-1; /* position of first and last sp. with gene members in preorder */
	/* these are the possible values for min and max number of miRNAs in the leafs */
	*kMin = 9999;
	*kMax = 0;

	/* find first (acc. to preorder) leaf with miRNAs and last one */
	for(i=0; i<nodesN; i++) 
	{
		lab = treeOrderN[i]; /* label of node */
		pos = treePosN[lab]; /* position of node in treeN */
		if(treeN[pos].m > 0) 
		{
			labFirst = lab; /* hold label of node */
			posFirst = i;   /* in preorde */
			break;
		}
	}
	FIRST=posFirst;
	printf("FIRST: preorder: %d pInT: %d label: %d (%s)\n",posFirst,pos,lab,treeN[pos].o);
  
	//int path2root1[nodesN-1];
	int * path2root1=  (int*) calloc(nodesN-1, sizeof(int));
  
	int steps1=0;
	j = treePosN[labFirst];
	while(treeN[j].p) 
	{
			//printf("treeN[%3d].n = %d\n",j,treeN[j].n); 
		path2root1[steps1] = treeN[j].n;
		j = treePosN[(treeN[j].p)->n];
		steps1++;
	}
	path2root1[steps1] = treeN[j].n;
				//for(i=0; i<=steps1; i++) { printf("%d,",path2root1[i]); } printf("\n");
	for(i=nodesN-1; i>=posFirst; i--) 
	{
		lab = treeOrderN[i]; /* label of node */
		pos = treePosN[lab]; /* position of node in treeN */
		if(treeN[pos].m > 0) 
		{
			labLast = lab; 
			posLast = i; /* in preorder */
			break;
		}
	}
	LAST=posLast;
	printf("LAST:  preorder: %d pInT: %d label: %d (%s)\n",posLast,pos,lab,treeN[pos].o);
	/* no mirs in tree ! species probably not in tree?! and this is the only species in the input file!!!*/

	if(posFirst > posLast) 
	{
		return -1;
	}
	for(i=posFirst; i<=posLast; i++) 
	{
		if(treeN[treePosN[treeOrderN[i]]].m > *kMax) *kMax = treeN[treePosN[treeOrderN[i]]].m;	
		if(treeN[treePosN[treeOrderN[i]]].m < *kMin) *kMin = treeN[treePosN[treeOrderN[i]]].m;
	}
	
	/*  not because of: ==9394== Conditional jump or move depends on uninitialised value(s)
												==9394==    at 0x409340: gl_getLCA (calc.c:1293)
												==9394==    by 0x401086: main (main.c:35)
								*/
	int * path2root2=  (int*) calloc(nodesN-1, sizeof(int));
	//int path2root2[nodesN-1];
								
	int steps2=0;
	j = treePosN[labLast];
	while(treeN[j].p) 
	{
		path2root2[steps2] = treeN[j].n;
		j = treePosN[(treeN[j].p)->n];
		steps2++;
	}
	path2root2[steps2] = treeN[j].n;
// for(i=0; i<=steps2; i++) { printf("%d,",path2root2[i]); } printf("\n");

  /* int maxsteps; */
  /* if(steps1 > steps2) { maxsteps = steps1; }  */
  /* else { maxsteps = steps2; } */

	if(posLast == posFirst) { lca = labFirst; }
	else 
	{
		printf("kMax: %d kMin: %d\n", *kMax, *kMin); 
    
		i = steps1; j = steps2;
		while(path2root1[i] == path2root2[j]) 
		{ 
			i--; 
			j--; 
		}
		i++; 
		j++;
    //printf("i: %d j: %d lca: %d (%s)\n",i,j,path2root1[i],treeN[treePosN[path2root1[i]]].o);
    
    
		lca = path2root1[i]; /*  label of LCA node! */
	}
  
	int lca_pos = treePosN[lca];
	int lca_ord = treeN[lca_pos].pInP;
	int posX;
		
	for(j=0; j<lca_ord; j++) 
	{
		lab = treeOrderN[j];
		posX = treePosN[lab];
		treeN[posX].seenBW = 1;
		treeN[posX].m = 0;
    /* printf("B: treeN[%d].o = %s .. initialized\n",posX,treeN[posX].o); */
	}
	//printf("--\n");

  /* find rightmost first pos in preorder, not element of subtree! */
  /* start in lca */

  //  int lastleafpos = lca_pos; /* pos in treeN */
	int max,maxI,curpos;

	curpos = lca_pos;
    
	while(treeN[curpos].nc > 0) 
	{
		lab = treeN[curpos].c[0];
		pos = treePosN[lab];
		max = treeN[pos].pInP;
		maxI = pos;
		for(i=1;i<treeN[curpos].nc;i++) 
		{
			lab = treeN[curpos].c[i];
			pos = treePosN[lab];
			if(treeN[pos].pInP > max) 
			{
				max = treeN[pos].pInP;
				maxI = pos;
			}
		}
		curpos = maxI;
	}

  /* printf("maxI: %d (%s)\n", curpos, treeN[curpos].o); */

	for(j=treeN[curpos].pInP+1; j<nodesN; j++) 
	{
  // for(j=lastleafpos+1; j<nodesN; j++) {
		lab = treeOrderN[j];
		posX = treePosN[lab];
		treeN[posX].seenBW = 1;
		treeN[posX].m = 0;
    /* printf("A: treeN[%d].o = %s .. initialized\n",posX,treeN[posX].o); */
	}
	
	free(path2root2);
	free(path2root1);
    return lca;
}


/* -------------------------------------------------------------------------- */

/* prints the treeN in all its elements as POSTSCRIPT */
void gl_printTreePSN(int kMax, int lca, struct gl_arguments ga) {

  int i, t, maxL=21, N=0, n_leaves=0, c, lca_pos;
  float x, y, w, h, x1, y1, x2, y2;
  
  //  float xn[N_NODES], yn[N_NODES];
  float *xn, *yn;

  //  printf("treeN[%d].m = %d\n",331,treeN[331].m);
  //printf("type: %s\n",ga.type);
  maxL = 0;
	for(i=0; i<nodesN; i++) 
	{
	    if(treeN[i].nc == 0) 
	    { 
			n_leaves++;
			//printf("number of leaves: %d", n_leaves);
		}
		if(treeN[i].level > maxL) 
		{ 
			maxL = treeN[i].level; 
		}
  }
  xn = (float *) calloc(nodesN+1, sizeof(float));
  yn = (float *) calloc(nodesN+1, sizeof(float));

  //  for(i=0; i<N_NODES; i++) { xn[i] = -1; yn[i] = -1; }
  for(i=0; i<nodesN; i++) { xn[i] = -1; yn[i] = -1; }

  lca_pos = treePosN[lca];
  /* width of level (w = 20) */
  w = 490 / (maxL+1);
  /* height of level (h = 9) */
  //h = 780 / N_LEAVES;  
  h = 820/ n_leaves;

  /* adjust fontsize of inner nodes to fit between the edges */
  int fontsize = h*0.8;
  if(fontsize > 10) fontsize = 10;

	
	char * psfname=NULL;
	psfname = (char *) calloc(strlen(ga.psfileFlag)+1, sizeof(char));
			strcpy(psfname, ga.psfileFlag);
	FILE *T;
	T = fopen(psfname, "w");
	 
	
  /* printing only the parsimony scores*/
  if (ga.b==0 && ga.z==0 )
  {		if (ga.collectfile == NULL)
		{
			printf("The Parsimony Scores only are visualized in:\t  %s\n", psfname);
		}
		else if(ga.collectfile != NULL)
		{
			printf("\nThe collected Tree is visualized in:\t  %s\n\n", psfname);
			
		}
	fprintf(T, "%%!PS-Adobe-2.1 EPSF-2.0\n");			
	fprintf(T, "%%"); fprintf(T, "%%Title: Efficient Prediction of Paralog Evolution - %s\n",psfname);
	fprintf(T, "%%"); fprintf(T, "%%Creator: ePoPE 1.0 -- by Jana Hertel\n");
	fprintf(T, "%%"); fprintf(T, "%%DocumentFonts: Helvetica\n");
	fprintf(T, "%%"); fprintf(T, "%%BoundingBox: 0 0 595 842\n");
	fprintf(T, "%%"); fprintf(T, "%%EndComments\n\n");
	
	/* start after Metazoa */
	x = 45; 
	y = 820;
	
	fprintf(T, "1 setlinecap\n");
	
	/* get nodes by tree preorder */
	//  for(i=0; i<329; i++) {
	//  for(i=0; i<nodesN; i++) {
	int drawnNodes = 0;
	i = 0;
  
	  
	  
		while(drawnNodes < nodesN) 
		{
			t = treePosN[treeOrderN[i]];
			//    printf("treeN[%d].o = %s\n",t,treeN[t].o);
			if(xn[t] != -1 && treeN[t].p !=NULL) 
			{
			i = (treeN[t].p)->pInP;
			goto nextWhile0;
			}
			/* draw the horizontal lines of the leafs */
			//    if(tree[t].l == NULL && tree[t].r == NULL) {
			if(treeN[t].nc == 0) 
			{ /* this is a leaf */
			/* x1 = x + ((tree[t].level) * w); */ /* all leaves end at their level */
			x1 = 505; /* all leaves end at 505 */
			y1 = y - (N * h) + 1;
			fprintf(T, "gsave newpath\n");
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto\n",x1,y1,x1+2,y1);
			fprintf(T, "%.2f %.2f lineto\n",x1+2,y1-2);
			fprintf(T, "%.2f %.2f lineto\n",x1,y1-2);
			fprintf(T, "%.2f %.2f lineto\n",x1,y1);
			fprintf(T, "closepath\n fill stroke\n");
			fprintf(T, "grestore\n");
			
			fprintf(T, "/Helvetica findfont 7 scalefont setfont\n");
			fprintf(T, "%.2f %.2f moveto (%s) show\n",565.,y1-3,treeN[t].o);
			/* node label */
			fprintf(T, "/Helvetica findfont 6 scalefont setfont\n");
			fprintf(T, "%.2f %.2f moveto (%-4d) show\n\n",550.,y1-3,treeN[t].n);
			/* node score (No. of genes) */
		
			fprintf(T, "gsave /Helvetica findfont 6 scalefont setfont\n");
			/*if(strcmp(type,"none") == 0) {
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%4d) show\n",507.,y1-3,treeN[t].m);
			} else {*/
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%4d) show\n",505.,y1-3,treeN[t].m);
			/*}*/
			fprintf(T, "grestore\n");
			
			/* leaf gain/loss */
			fprintf(T, "gsave /Helvetica findfont 6 scalefont setfont\n");
			/*if(strcmp(type,"none") != 0) {*/
			fprintf(T, "0 0 1 setrgbcolor %.2f %.2f moveto (%2d) show\n",520.,y1-3,treeN[t].gain);
			fprintf(T, "1 0 0 setrgbcolor %.2f %.2f moveto (%2d) show\n",535.,y1-3,treeN[t].loss);
			fprintf(T, "grestore\n");
			/*      }*/
			
			drawnNodes++;
			xn[t] = x1; yn[t] = y1;
			N++;
			i = (treeN[t].p)->pInP;
			goto nextWhile0;
		
			} else { /* inner node */
			/* check if all children have coords */
			for(c=0; c<treeN[t].nc; c++) { /* if not step to uncoordinated children recursively */
			if(xn[treePosN[treeN[t].c[c]]] == -1) { 
			i = treeN[treePosN[treeN[t].c[c]]].pInP;
			goto nextWhile0;
			}
			}
			/* all children coordinated */
			/* yn[t] = yn[treeN[treePosN[treeN[t].c[0]]].pInP]; */
			/* float yChildren=yn[treePosN[treeN[t].c[treeN[t].nc-1]]]; */
			float yChildren = 0;
			for(c=0; c<treeN[t].nc; c++) {
			yChildren += yn[treePosN[treeN[t].c[c]]];
			}
			yChildren /= treeN[t].nc;
			yn[t] = yChildren;
			xn[t] = x + (treeN[t].level*w);
			x1 = xn[t]; y1 = yn[t];
			if(t == lca_pos) {
			fprintf(T, "gsave\n");
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto\n",x1,y1-1,x1,y1+3);
			fprintf(T, "%.2f %.2f lineto %.2f %.2f lineto\n",x1-5,y1-1,x1,y1-5);
			fprintf(T, "%.2f %.2f lineto\n",x1,y1-1);
			fprintf(T, "fill closepath\nstroke\n");
			fprintf(T, "grestore\n");
			}
			
			/* node label */
			fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "%.2f %.2f moveto (%s) show\n",x1-20,y1+1,treeN[t].o);
			
			/* node score */
			if(ga.type == NULL) {
			fprintf(T, "gsave /Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%d) show\n",x1+2,y1-7,treeN[t].m);
			fprintf(T, "grestore\n");
			} else {
			if(strcmp(ga.type, "all") == 0 || strcmp(ga.type, "genes") == 0) {
			/* if(treeN[t].m > 0) { */
				fprintf(T, "gsave /Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%d) show\n",x1+2,y1-7,treeN[t].m);
				fprintf(T, "grestore\n");
			/* } */
			} 
			}
			
			/* gain */
			if(ga.type == NULL)
			{
			if(treeN[t].gain > 0) {
			fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "gsave 0 0 1 setrgbcolor\n");
			fprintf(T, "%.2f %.2f moveto (\\(%d\\)) show\n",x1+8,y1-7,treeN[t].gain);
			fprintf(T, "grestore\n");
			}
			} else {
			if(strcmp(ga.type, "all") == 0) {
			if(treeN[t].gain > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 0 0 1 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(%d\\)) show\n",x1+8,y1-7,treeN[t].gain);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "gain") == 0) {
			if(treeN[t].gain > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 0 0 1 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(%d\\)) show\n",x1+8,y1-7,treeN[t].gain);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "gainFam") == 0) {
			if(treeN[t].gain > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 0 0 1 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-7,treeN[t].gainFam);
				fprintf(T, "grestore\n");
			}	
			}
			}
		
			/* loss */
			if(ga.type == NULL) {
			if(treeN[t].loss > 0) {
			fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "gsave 1 0 0 setrgbcolor\n");
			fprintf(T, "%.2f %.2f moveto (\\(%d\\)) show\n",x1+16,y1-7,treeN[t].loss);
			fprintf(T, "grestore\n");
			}
			} else { 
			if(strcmp(ga.type, "all") == 0) {
			if(treeN[t].loss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(%d\\)) show\n",x1+16,y1-7,treeN[t].loss);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "loss") == 0) {
			if(treeN[t].loss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(%d\\)) show\n",x1+16,y1-7,treeN[t].loss);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "lossFam") == 0) {
			if(treeN[t].loss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-7,treeN[t].lossFam);
				fprintf(T, "grestore\n");
			}
			}
			}
			
			/* stroke paths to all children */
			for(c=0; c<treeN[t].nc; c++) {
			x2 = xn[treePosN[treeN[t].c[c]]];
			y2 = yn[treePosN[treeN[t].c[c]]]-1;
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto\n",x1,y1,x1,y2);
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto stroke\n",x1,y2,x2,y2);
			}
			drawnNodes++;
		
			goto nextWhile0;
		
			nextWhile0:
			continue;
			}
		} 
	}
	
	/*printing the PF scores without any comparison*/
	
	if (ga.z==1)
	{	
	printf("The Partition function Scores only are visualized in the Tree: %s\n", psfname);	

	fprintf(T, "%%!PS-Adobe-2.1 EPSF-2.0\n");			
	fprintf(T, "%%"); fprintf(T, "%%Title: Efficient Prediction of Paralog Evolution - %s\n",psfname);
	fprintf(T, "%%"); fprintf(T, "%%Creator: ePoPE 1.0 -- by Jana Hertel\n");
	fprintf(T, "%%"); fprintf(T, "%%DocumentFonts: Helvetica\n");
	fprintf(T, "%%"); fprintf(T, "%%BoundingBox: 0 0 595 842\n");
	fprintf(T, "%%"); fprintf(T, "%%EndComments\n\n");
	
	/* start after Metazoa */
	x = 45; 
	y = 820;
	
	fprintf(T, "1 setlinecap\n");
	
	/* get nodes by tree preorder */
	//  for(i=0; i<329; i++) {
	//  for(i=0; i<nodesN; i++) {
	int drawnNodes = 0;
	i = 0;
		
		
		
		
		while(drawnNodes < nodesN) 
		{
			t = treePosN[treeOrderN[i]];
			//    printf("treeN[%d].o = %s\n",t,treeN[t].o);
			if(xn[t] != -1 && treeN[t].p !=NULL) 
			{
			i = (treeN[t].p)->pInP;
			goto nextWhilez;
			}
			/* draw the horizontal lines of the leafs */
			//    if(tree[t].l == NULL && tree[t].r == NULL) {
			if(treeN[t].nc == 0) 
			{ /* this is a leaf */
			/* x1 = x + ((tree[t].level) * w); */ /* all leaves end at their level */
			x1 = 505; /* all leaves end at 505 */
			y1 = y - (N * h) + 1;
			fprintf(T, "gsave newpath\n");
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto\n",x1,y1,x1+2,y1);
			fprintf(T, "%.2f %.2f lineto\n",x1+2,y1-2);
			fprintf(T, "%.2f %.2f lineto\n",x1,y1-2);
			fprintf(T, "%.2f %.2f lineto\n",x1,y1);
			fprintf(T, "closepath\n fill stroke\n");
			fprintf(T, "grestore\n");
			
			fprintf(T, "/Helvetica findfont 7 scalefont setfont\n");
			fprintf(T, "%.2f %.2f moveto (%s) show\n",565.,y1-3,treeN[t].o);
			/* node label */
			fprintf(T, "/Helvetica findfont 6 scalefont setfont\n");
			fprintf(T, "%.2f %.2f moveto (%-4d) show\n\n",550.,y1-3,treeN[t].n);
			/* node score (No. of genes) */
		
			fprintf(T, "gsave /Helvetica findfont 6 scalefont setfont\n");
			/*if(strcmp(type,"none") == 0) {
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%4d) show\n",507.,y1-3,treeN[t].m);
			} else {*/
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%4d) show\n",505.,y1-3,treeN[t].m);
			/*}*/
			fprintf(T, "grestore\n");
			
			/* leaf gain/loss */
			fprintf(T, "gsave /Helvetica findfont 6 scalefont setfont\n");
			/*if(strcmp(type,"none") != 0) {*/
			fprintf(T, "0 0 1 setrgbcolor %.2f %.2f moveto (%2d) show\n",520.,y1-3,treeN[t].pfGain);
			fprintf(T, "1 0 0 setrgbcolor %.2f %.2f moveto (%2d) show\n",535.,y1-3,treeN[t].pfLoss);
			fprintf(T, "grestore\n");
			/*      }*/
			
			drawnNodes++;
			xn[t] = x1; yn[t] = y1;
			N++;
			i = (treeN[t].p)->pInP;
			goto nextWhilez;
		
			} else { /* inner node */
			/* check if all children have coords */
			for(c=0; c<treeN[t].nc; c++) { /* if not step to uncoordinated children recursively */
			if(xn[treePosN[treeN[t].c[c]]] == -1) { 
			i = treeN[treePosN[treeN[t].c[c]]].pInP;
			goto nextWhilez;
			}
			}
			/* all children coordinated */
			/* yn[t] = yn[treeN[treePosN[treeN[t].c[0]]].pInP]; */
			/* float yChildren=yn[treePosN[treeN[t].c[treeN[t].nc-1]]]; */
			float yChildren = 0;
			for(c=0; c<treeN[t].nc; c++) {
			yChildren += yn[treePosN[treeN[t].c[c]]];
			}
			yChildren /= treeN[t].nc;
			yn[t] = yChildren;
			xn[t] = x + (treeN[t].level*w);
			x1 = xn[t]; y1 = yn[t];
			if(t == lca_pos) {
			fprintf(T, "gsave\n");
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto\n",x1,y1-1,x1,y1+3);
			fprintf(T, "%.2f %.2f lineto %.2f %.2f lineto\n",x1-5,y1-1,x1,y1-5);
			fprintf(T, "%.2f %.2f lineto\n",x1,y1-1);
			fprintf(T, "fill closepath\nstroke\n");
			fprintf(T, "grestore\n");
			}
			
			/* node label */
			fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "%.2f %.2f moveto (%s) show\n",x1-20,y1+1,treeN[t].o);
			
			/* node score */
			if(ga.type == NULL) {
			fprintf(T, "gsave /Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%d) show\n",x1+2,y1-7,treeN[t].pfm);
			fprintf(T, "grestore\n");
			} else {
			if(strcmp(ga.type, "all") == 0 || strcmp(ga.type, "genes") == 0) {
			/* if(treeN[t].m > 0) { */
				fprintf(T, "gsave /Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%d) show\n",x1+2,y1-7,treeN[t].pfm);
				fprintf(T, "grestore\n");
			/* } */
			} 
			}
		
				
			/* gain */
			if(ga.type == NULL)
			{
			if(treeN[t].pfGain > 0) {
			fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "gsave 0 0 1 setrgbcolor\n");
			fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-7,treeN[t].pfGain);
			fprintf(T, "grestore\n");
			}
			} else {
			if(strcmp(ga.type, "all") == 0) {
			if(treeN[t].pfGain > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 0 0 1 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-7,treeN[t].pfGain);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "gain") == 0) {
			if(treeN[t].pfGain > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 0 0 1 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-7,treeN[t].pfGain);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "gainFam") == 0) {
			if(treeN[t].pfGain > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 0 0 1 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-7,treeN[t].pfGainFam);
				fprintf(T, "grestore\n");
			}	
			}
			}
		
			/* loss */
			if(ga.type == NULL) {
			if(treeN[t].pfLoss > 0) {
			fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "gsave 1 0 0 setrgbcolor\n");
			fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-7,treeN[t].pfLoss);
			fprintf(T, "grestore\n");
			}
			} else { 
			if(strcmp(ga.type, "all") == 0) {
			if(treeN[t].pfLoss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-7,treeN[t].pfLoss);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "loss") == 0) {
			if(treeN[t].pfLoss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-7,treeN[t].pfLoss);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "lossFam") == 0) {
			if(treeN[t].pfLossFam > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-7,treeN[t].pfLossFam);
				fprintf(T, "grestore\n");
			}
			}
			}
			
			/* stroke paths to all children */
			for(c=0; c<treeN[t].nc; c++) {
			x2 = xn[treePosN[treeN[t].c[c]]];
			y2 = yn[treePosN[treeN[t].c[c]]]-1;
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto\n",x1,y1,x1,y2);
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto stroke\n",x1,y2,x2,y2);
			}
			drawnNodes++;
		
			goto nextWhilez;
		
			nextWhilez:
			continue;
			}
		} 
	}
	
	//print the comparison between parsimony and PF:
	if (ga.b==1)
	{	
		printf("A comparison of Parsimony Scores and Partition function is visualized in the Tree: %s\n", psfname);

	fprintf(T, "%%!PS-Adobe-2.1 EPSF-2.0\n");			
	fprintf(T, "%%"); fprintf(T, "%%Title: Efficient Prediction of Paralog Evolution - %s\n",ga.psfile);
	fprintf(T, "%%"); fprintf(T, "%%Creator: ePoPE 1.0 -- by Jana Hertel\n");
	fprintf(T, "%%"); fprintf(T, "%%DocumentFonts: Helvetica\n");
	fprintf(T, "%%"); fprintf(T, "%%BoundingBox: 0 0 595 842\n");
	fprintf(T, "%%"); fprintf(T, "%%EndComments\n\n");
	
	/* start after Metazoa */
	x = 45; 
	y = 820;
	
	fprintf(T, "1 setlinecap\n");
	
	/* get nodes by tree preorder */
	//  for(i=0; i<329; i++) {
	//  for(i=0; i<nodesN; i++) {
	int drawnNodes = 0;
	i = 0;
		
		
		
		while(drawnNodes < nodesN) 
		{
			t = treePosN[treeOrderN[i]];
			//    printf("treeN[%d].o = %s\n",t,treeN[t].o);
			if(xn[t] != -1 && treeN[t].p !=NULL) 
			{
			i = (treeN[t].p)->pInP;
			goto nextWhileb;
			}
			/* draw the horizontal lines of the leafs */
			//    if(tree[t].l == NULL && tree[t].r == NULL) {
			if(treeN[t].nc == 0) 
			{ /* this is a leaf */
			/* x1 = x + ((tree[t].level) * w); */ /* all leaves end at their level */
			x1 = 505; /* all leaves end at 505 */
			y1 = y - (N * h) + 1;
			fprintf(T, "gsave newpath\n");
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto\n",x1,y1,x1+2,y1);
			fprintf(T, "%.2f %.2f lineto\n",x1+2,y1-2);
			fprintf(T, "%.2f %.2f lineto\n",x1,y1-2);
			fprintf(T, "%.2f %.2f lineto\n",x1,y1);
			fprintf(T, "closepath\n fill stroke\n");
			fprintf(T, "grestore\n");
			


/*leaves for PS, those are printed before the PF*/
			fprintf(T, "/Helvetica findfont 7 scalefont setfont\n");
			fprintf(T, "%.2f %.2f moveto (%s) show\n",575.,y1-3,treeN[t].o);
			/* node label */
			fprintf(T, "/Helvetica findfont 6 scalefont setfont\n");
			fprintf(T, "%.2f %.2f moveto (%-4d) show\n\n",560.,y1-3,treeN[t].n);
			/* node score (No. of genes) */
		
			fprintf(T, "gsave /Helvetica findfont 6 scalefont setfont\n");
			/*if(strcmp(type,"none") == 0) {
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%4d) show\n",507.,y1-3,treeN[t].m);
			} else {*/
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%4d) show\n",503.,y1-3,treeN[t].m);
			/*}*/
			fprintf(T, "grestore\n");
			
			/* leaf gain/loss */
			fprintf(T, "gsave /Helvetica findfont 6 scalefont setfont\n");
			/*if(strcmp(type,"none") != 0) {*/
			fprintf(T, "0 0 1 setrgbcolor %.2f %.2f moveto (%2d) show\n",524.,y1-3,treeN[t].gain);
			fprintf(T, "1 0 0 setrgbcolor %.2f %.2f moveto (%2d) show\n",543.,y1-3,treeN[t].loss);
/*leaves for PF, those are printed behind the PS*/
			//fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%4d) show\n",511.,y1-3,treeN[t].pfm);  /braucht man nicht, k an leaves bei beiden gleich
			/*}*/
			fprintf(T, "grestore\n");
			
			/* leaf gain/loss */
			fprintf(T, "gsave /Helvetica findfont 6 scalefont setfont\n");
			/*if(strcmp(type,"none") != 0) {*/
			fprintf(T, "0 0 1 setrgbcolor %.2f %.2f moveto (%2d) show\n",532.,y1-3,treeN[t].pfGain);
			fprintf(T, "1 0 0 setrgbcolor %.2f %.2f moveto (%2d) show\n",551.,y1-3,treeN[t].pfLoss);
			fprintf(T, "grestore\n");
			/*      }*/


			
			
			
			
			drawnNodes++;
			xn[t] = x1; yn[t] = y1;
			N++;
			i = (treeN[t].p)->pInP;
			goto nextWhileb;
		
			} else { /* inner node */
			/* check if all children have coords */
			for(c=0; c<treeN[t].nc; c++) { /* if not step to uncoordinated children recursively */
			if(xn[treePosN[treeN[t].c[c]]] == -1) { 
			i = treeN[treePosN[treeN[t].c[c]]].pInP;
			goto nextWhileb;
			}
			}
			/* all children coordinated */
			/* yn[t] = yn[treeN[treePosN[treeN[t].c[0]]].pInP]; */
			/* float yChildren=yn[treePosN[treeN[t].c[treeN[t].nc-1]]]; */
			float yChildren = 0;
			for(c=0; c<treeN[t].nc; c++) {
			yChildren += yn[treePosN[treeN[t].c[c]]];
			}
			yChildren /= treeN[t].nc;
			yn[t] = yChildren;
			xn[t] = x + (treeN[t].level*w);
			x1 = xn[t]; y1 = yn[t];
			if(t == lca_pos) {
			fprintf(T, "gsave\n");
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto\n",x1,y1-1,x1,y1+3);
			fprintf(T, "%.2f %.2f lineto %.2f %.2f lineto\n",x1-5,y1-1,x1,y1-5);
			fprintf(T, "%.2f %.2f lineto\n",x1,y1-1);
			fprintf(T, "fill closepath\nstroke\n");
			fprintf(T, "grestore\n");
			}
			
			/* node label */
			fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "%.2f %.2f moveto (%s) show\n",x1-20,y1+1,treeN[t].o);
			
			/* node score */
			if(ga.type == NULL) {
			fprintf(T, "gsave /Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%d) show\n",x1+2,y1-7,treeN[t].m);
			fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%d) show\n",x1+2,y1-14,treeN[t].pfm);
			fprintf(T, "grestore\n");
			} else {
			if(strcmp(ga.type, "all") == 0 || strcmp(ga.type, "genes") == 0) {
			/* if(treeN[t].m > 0) { */
				fprintf(T, "gsave /Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%d) show\n",x1+2,y1-7,treeN[t].m);
				fprintf(T, "0 0.4 0 setrgbcolor %.2f %.2f moveto (%d) show\n",x1+2,y1-14,treeN[t].pfm);
				fprintf(T, "grestore\n");
			/* } */
			} 
			}
		
		
		
		
			/* gain PS */
			if(ga.type == NULL)
			{
				if(treeN[t].gain > 0) 
				{
					fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
					fprintf(T, "gsave 0 0 1 setrgbcolor\n");
					fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-7,treeN[t].gain);
					fprintf(T, "grestore\n");
				}
			} else 
			{
				if(strcmp(ga.type, "all") == 0) 
				{
					if(treeN[t].gain > 0) 
					{
						fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
						fprintf(T, "gsave 0 0 1 setrgbcolor\n");
						fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-7,treeN[t].gain);
						fprintf(T, "grestore\n");
					}
				} 	
				else if(strcmp(ga.type, "gain") == 0) 
				{
					if(treeN[t].gain > 0) 	
					{
						fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
						fprintf(T, "gsave 0 0 1 setrgbcolor\n");
						fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-7,treeN[t].gain);
						fprintf(T, "grestore\n");
					}
				} 	
					else if(strcmp(ga.type, "gainFam") == 0) 
					{
						if(treeN[t].gain > 0) 	
						{
							fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
							fprintf(T, "gsave 0 0 1 setrgbcolor\n");
							fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-7,treeN[t].gainFam);
							fprintf(T, "grestore\n");
						}	
					}
			}

		/* gain PF */
			if(ga.type == NULL)
			{
				if(treeN[t].pfGain > 0) 
				{
					fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
					fprintf(T, "gsave 0 0 1 setrgbcolor\n");
					fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-14,treeN[t].pfGain);
					fprintf(T, "grestore\n");
				}
			} else 
			{
				if(strcmp(ga.type, "all") == 0) 
				{
					if(treeN[t].pfGain > 0) 
					{
						fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
						fprintf(T, "gsave 0 0 1 setrgbcolor\n");
						fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-14,treeN[t].pfGain);
						fprintf(T, "grestore\n");
					}
				} 	
				else if(strcmp(ga.type, "gain") == 0) 
				{
					if(treeN[t].pfGain > 0) 	
					{
						fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
						fprintf(T, "gsave 0 0 1 setrgbcolor\n");
						fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-14,treeN[t].pfGain);
						fprintf(T, "grestore\n");
					}
				} 	
					else if(strcmp(ga.type, "gainFam") == 0) 
					{
						if(treeN[t].pfGainFam > 0) 	
						{
							fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
							fprintf(T, "gsave 0 0 1 setrgbcolor\n");
							fprintf(T, "%.2f %.2f moveto (\\(+%d\\)) show\n",x1+10,y1-14,treeN[t].pfGainFam);//kann man lassen, gainFam zw PF und OS sind gleich
							fprintf(T, "grestore\n");
						}	
					}
			}
		
			/* loss PS*/
			if(ga.type == NULL) {
			if(treeN[t].loss > 0) {
			fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "gsave 1 0 0 setrgbcolor\n");
			fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-7,treeN[t].loss);
			fprintf(T, "grestore\n");
			}
			} else { 
			if(strcmp(ga.type, "all") == 0) {
			if(treeN[t].loss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-7,treeN[t].loss);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "loss") == 0) {
			if(treeN[t].loss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-7,treeN[t].loss);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "lossFam") == 0) {
			if(treeN[t].loss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-7,treeN[t].lossFam);
				fprintf(T, "grestore\n");
			}
			}
			}
			
			/* loss PF*/
			if(ga.type == NULL) {
			if(treeN[t].pfLoss > 0) {
			fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
			fprintf(T, "gsave 1 0 0 setrgbcolor\n");
			fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-14,treeN[t].pfLoss);
			fprintf(T, "grestore\n");
			}
			} else { 
			if(strcmp(ga.type, "all") == 0) {
			if(treeN[t].pfLoss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-14,treeN[t].pfLoss);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "loss") == 0) {
			if(treeN[t].pfLoss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-14,treeN[t].pfLoss);
				fprintf(T, "grestore\n");
			}
			} else if(strcmp(ga.type, "lossFam") == 0) {
			if(treeN[t].pfLoss > 0) {
				fprintf(T, "/Helvetica findfont %d scalefont setfont\n", fontsize);
				fprintf(T, "gsave 1 0 0 setrgbcolor\n");
				fprintf(T, "%.2f %.2f moveto (\\(-%d\\)) show\n",x1+18,y1-14,treeN[t].pfLossFam);
				fprintf(T, "grestore\n");
			}
			}
			}
				
			
			/* stroke paths to all children */
			for(c=0; c<treeN[t].nc; c++) {
			x2 = xn[treePosN[treeN[t].c[c]]];
			y2 = yn[treePosN[treeN[t].c[c]]]-1;
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto\n",x1,y1,x1,y2);
			fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto stroke\n",x1,y2,x2,y2);
			}
			drawnNodes++;
		
			goto nextWhileb;
		
			nextWhileb:
			continue;
			}
		} 
	}
	
	
	
  fprintf(T, "%.2f %.2f moveto %.2f %.2f lineto stroke\n\n",xn[0],yn[0],xn[0]-w,yn[0]);
  fprintf(T, "showpage\n");
  fflush(T);
  fclose(T);
 // printf("filename: %s\n", fname);
  free(xn);
  free(yn);
  
  free(psfname);
  
}

/* -------------------------------------------------------------------------- */

void gl_printTreeN(struct gl_arguments ga) 
{
	

	int j, c;
	//char* fname= ga.outfile;
	
	//~ printf("b: %d\n", ga.b);
	//~ printf("z: %d\n", ga.z);
	
	if (ga.z==0 && ga.b==0)	//ausdrucken nur des Parsimony scores
	{	
		if(ga.outfile == NULL) 
		{ 	printf("\nParsimony table:\n");
			for(j=0; j<nodesN; j++) 
			{
				printf("pos: %3d\t lab: %3d\t kids: ",treeN[j].pInT,treeN[j].n);
				if(treeN[j].nc == 0) { printf("%d ", -1); }
				else 
				{
					for(c=0; c<treeN[j].nc-1; c++) 
					{ 
						printf("%d,", treeN[j].c[c]);
					}
					
					printf("%d ", treeN[j].c[(treeN[j].nc)-1]);
				}
				
				printf("\tpInP: %3d\t ",treeN[j].pInP);
				if(treeN[j].p == NULL) printf("par:  -1\t "); 
				else printf("par: %3d\t ",(treeN[j].p)->n);
				printf("genes: %3d\t ",treeN[j].m);
				printf("gain: %3d\t loss %3d\t ",treeN[j].gain,treeN[j].loss);
				printf("gainFam: %3d\t lossFam: %3d\t %s\n",treeN[j].gainFam,treeN[j].lossFam,treeN[j].o);
			}
		
		} 
		else /* write to file */
		{ 	//char praefix[]="PS.";
			//strcat(praefix,fname);
			
			//~ char suffix2[strlen(fname)+4];
			//~ strcpy(suffix2, fname);
			//~ strcat(suffix2, "_PS");
			char* fname= ga.outfileFlag;
			printf("The Parsimony Scores only are stored in table:\t  %s\n",  fname);
			FILE *fp=NULL;
			fp = fopen(fname, "w");
			
			for(j=0; j<nodesN; j++) 
			{
				fprintf(fp,"pos: %3d\t lab: %3d\t kids: ",treeN[j].pInT,treeN[j].n);
				if(treeN[j].nc == 0) 
				{ 
					fprintf(fp, "%d ", -1); 
				}
				else 
				{
					for(c=0; c<treeN[j].nc-1; c++) 
					{ 
						fprintf(fp, "%d,", treeN[j].c[c]);
					}
					fprintf(fp, "%d ", treeN[j].c[(treeN[j].nc)-1]);
				}
				fprintf(fp,"\tpInP: %3d\t ",treeN[j].pInP);
				if(treeN[j].p == NULL) fprintf(fp, "par:  -1\t "); 
				else fprintf(fp,"par: %3d\t ",(treeN[j].p)->n);
				fprintf(fp,"genes: %3d\t ",treeN[j].m);
				fprintf(fp,"gain: %3d\t loss: %3d\t ",treeN[j].gain,treeN[j].loss);
				fprintf(fp,"gainFam: %3d\t lossFam: %3d\t %s\n",treeN[j].gainFam,treeN[j].lossFam,treeN[j].o);
		
			}
			fclose(fp);
		}
	}
	
	if (ga.z==1)	//ausdrucken nur des PF scores
	{	
		if(ga.outfile == NULL) 
		{ /* write to stdout */
			printf("\nPartition function table:\n");
			for(j=0; j<nodesN; j++) 
			{
				printf("pos: %3d\t lab: %3d\t kids: ",treeN[j].pInT,treeN[j].n);
				if(treeN[j].nc == 0) { printf("%d ", -1); }
				else 
				{
					for(c=0; c<treeN[j].nc-1; c++) 
					{ 
						printf("%d,", treeN[j].c[c]);
					}
					
					printf("%d ", treeN[j].c[(treeN[j].nc)-1]);
				}
				
				printf("\tpInP: %3d\t ",treeN[j].pInP);
				if(treeN[j].p == NULL) printf("par:  -1\t "); 
				else printf("par: %3d\t ",(treeN[j].p)->n);
				printf("genes: %3d\t ",treeN[j].pfm);
				printf("gain: %3d\t loss %3d\t ",treeN[j].pfGain,treeN[j].pfLoss);
				//printf("gain: %3d loss %3d\t seenBW: %3d ",treeN[j].pfGain,treeN[j].pfLoss, 1);
				printf("gainFam: %3d\t lossFam: %3d\t PF_P: %3f\t %s \n",treeN[j].pfGainFam,treeN[j].pfLossFam,treeN[j].pfP,treeN[j].o  );
			}
		
		} 
		else /* write to file */
		{ 	
			char* fname= ga.outfileFlag;
			
			printf("Partition function score is stored in %s\n", fname);
			FILE *fp=NULL;
			fp = fopen(fname, "w");
			for(j=0; j<nodesN; j++) 
			{
				fprintf(fp,"pos: %3d\t lab: %3d\t kids: ",treeN[j].pInT,treeN[j].n);
				if(treeN[j].nc == 0) 
				{ 
					fprintf(fp, "%d ", -1); 
				}
				else 
				{
					for(c=0; c<treeN[j].nc-1; c++) 
					{ 
						fprintf(fp, "%d,", treeN[j].c[c]);
					}
					fprintf(fp, "%d ", treeN[j].c[(treeN[j].nc)-1]);
				}
				fprintf(fp,"\tpInP: %3d\t ",treeN[j].pInP);
				if(treeN[j].p == NULL) fprintf(fp, "par:  -1\t "); 
				else fprintf(fp,"par: %3d\t ",(treeN[j].p)->n);
				fprintf(fp,"genes: %3d\t ",treeN[j].pfm);
				fprintf(fp,"gain: %3d\t loss %3d\t ",treeN[j].pfGain,treeN[j].pfLoss);
				//fprintf(fp,"gain: %3d loss %3d\t seenBW: %3d ",treeN[j].pfGain,treeN[j].pfLoss,1);
				fprintf(fp,"gainFam: %3d\t lossFam: %3d\t PF_P: %3f\t %s\n",treeN[j].pfGainFam,treeN[j].pfLossFam,treeN[j].pfP,treeN[j].o );
		
			}
			fclose(fp);
		}
	}
	
	
	if (ga.b==1)	//ausdrucken von 2 tabellen mit suffix zum unterscheiden
	{						//erst PS table:
		if(ga.outfileFlag == NULL) 
		{ 	printf("\nParsimony table:\n");
			for(j=0; j<nodesN; j++) 
			{
				printf("pos: %3d\t lab: %3d\t kids: ",treeN[j].pInT,treeN[j].n);
				if(treeN[j].nc == 0) { printf("%d ", -1); }
				else 
				{
					for(c=0; c<treeN[j].nc-1; c++) 
					{ 
						printf("%d,", treeN[j].c[c]);
					}
					
					printf("%d ", treeN[j].c[(treeN[j].nc)-1]);
				}
				
				printf("\tpInP: %3d\t ",treeN[j].pInP);
				if(treeN[j].p == NULL) printf("par:  -1\t "); 
				else printf("par: %3d\t ",(treeN[j].p)->n);
				printf("genes: %3d\t ",treeN[j].m);
				printf("gain: %3d\t loss %3d\t",treeN[j].gain,treeN[j].loss);
				printf("gainFam: %3d\t lossFam: %3d\t PS_P: %3f\t %s \n",treeN[j].gainFam,treeN[j].lossFam,treeN[j].psP,treeN[j].o  );
			}
		
		} 
		else /* write to file */
		{ 	//char praefix[]="_PS";  //praefix erst mal als suffix, muss später noch nach letztem punkt in filename eingefuegt werden!
			
			
			//~ char suffix[strlen(fname)+4];
			//~ strcpy(suffix, fname);
			//~ strcat(suffix, "_PSb");
									
			printf("Parsimony Score is stored in %s\n", ga.outfileFlag );
			char *fname=ga.outfileFlag;
			FILE *fp=NULL;
			fp = fopen(fname, "w");
			for(j=0; j<nodesN; j++) 
			{
				fprintf(fp,"pos: %3d\t lab: %3d\t kids: ",treeN[j].pInT,treeN[j].n);
				if(treeN[j].nc == 0) 
				{ 
					fprintf(fp, "%d ", -1); 
				}
				else 
				{
					for(c=0; c<treeN[j].nc-1; c++) 
					{ 
						fprintf(fp, "%d,", treeN[j].c[c]);
					}
					fprintf(fp, "%d ", treeN[j].c[(treeN[j].nc)-1]);
				}
				fprintf(fp,"\tpInP: %3d\t ",treeN[j].pInP);
				if(treeN[j].p == NULL) fprintf(fp, "par:  -1\t "); 
				else fprintf(fp,"par: %3d\t ",(treeN[j].p)->n);
				fprintf(fp,"genes: %3d\t ",treeN[j].m);
				fprintf(fp,"gain: %3d\t loss %3d\t",treeN[j].gain,treeN[j].loss);
				fprintf(fp,"gainFam: %3d\t lossFam: %3d\t PS_P: %3f\t %s \n",treeN[j].gainFam,treeN[j].lossFam,treeN[j].psP,treeN[j].o );
		
			}
			fclose(fp);
		}
		
		//dann Pf table:
		
		if(ga.outfileFlag_b == NULL) 
		{ /* write to stdout */
			printf("\nPartition function table:\n");
			for(j=0; j<nodesN; j++) 
			{
				printf("pos: %3d\t lab: %3d\t kids: ",treeN[j].pInT,treeN[j].n);
				if(treeN[j].nc == 0) { printf("%d ", -1); }
				else 
				{
					for(c=0; c<treeN[j].nc-1; c++) 
					{ 
						printf("%d,", treeN[j].c[c]);
					}
					
					printf("%d ", treeN[j].c[(treeN[j].nc)-1]);
				}
				
				printf("\tpInP: %3d\t ",treeN[j].pInP);
				if(treeN[j].p == NULL) printf("par:  -1\t "); 
				else printf("par: %3d\t ",(treeN[j].p)->n);
				printf("genes: %3d\t ",treeN[j].pfm);
				printf("gain: %3d\t loss %3d\t",treeN[j].pfGain,treeN[j].pfLoss);
				//printf("gain: %3d loss %3d\t seenBW: %3d ",treeN[j].pfGain,treeN[j].pfLoss, 1);
				printf("gainFam: %3d\t lossFam: %3d\t PF_P: %3f\t %s\n",treeN[j].pfGainFam,treeN[j].pfLossFam,treeN[j].pfP,treeN[j].o );
			}
		
		} 
		else /* write to file */
		{ 	
			//~ char suffix2[strlen(fname)+4];
			//~ strcpy(suffix2, fname);
			//~ strcat(suffix2, "_PFb");
			char* fname=ga.outfileFlag_b;
			
			printf("Partition function score is stored in %s\n",fname);
			FILE *fp=NULL;
			fp = fopen(fname, "w");
			for(j=0; j<nodesN; j++) 
			{
				fprintf(fp,"pos: %3d\t lab: %3d\t kids: ",treeN[j].pInT,treeN[j].n);
				if(treeN[j].nc == 0) 
				{ 
					fprintf(fp, "%d ", -1); 
				}
				else 
				{
					for(c=0; c<treeN[j].nc-1; c++) 
					{ 
						fprintf(fp, "%d,", treeN[j].c[c]);
					}
					fprintf(fp, "%d ", treeN[j].c[(treeN[j].nc)-1]);
				}
				fprintf(fp,"\tpInP: %3d\t ",treeN[j].pInP);
				if(treeN[j].p == NULL) fprintf(fp, "par:  -1\t "); 
				else fprintf(fp,"par: %3d\t ",(treeN[j].p)->n);
				fprintf(fp,"genes: %3d\t ",treeN[j].pfm);
				fprintf(fp,"gain: %3d\t loss %3d\t",treeN[j].pfGain,treeN[j].pfLoss);
				//fprintf(fp,"gain: %3d loss %3d\t seenBW: %3d ",treeN[j].pfGain,treeN[j].pfLoss,1);
				fprintf(fp,"gainFam: %3d\t lossFam: %3d\t PF_P: %3f\t %s\n",treeN[j].pfGainFam,treeN[j].pfLossFam,treeN[j].pfP,treeN[j].o );
		
			}
			fclose(fp);
		}
		
		
	}
	
	
  return;
}

/* -------------------------------------------------------------------------- */

void gl_printS(float **S, int m, int n) {

  int i,j;

  printf("-   \t");
  for(j=0; j<n; j++){
    printf("  %c%c%c \t",treeN[j].o[0],treeN[j].o[1],treeN[j].o[2]);
  }
  printf("\n");

  printf("i\\j \t");
  for(j=0; j<n; j++)
    printf("%5d\t",j);
  printf("\n");

  //  printf("m: %d n: %d\n", m, n); exit(0)
  for(i=0; i<m; i++) {
    printf("%3d:\t",i);
    for(j=0; j<n; j++) {
      printf("%5.2f\t", S[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  return;
  
}

/* -------------------------------------------------------------------------- */

