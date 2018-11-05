#ifndef TREE_H			//wenn TREE_H NICHT definiert ist, definiere TREE_H 
#define TREE_H			// hier wird TREE_H definiert, der folgende text bis endif wird compiliert
//#if defined TREE_H
//#error "TREE_H defined"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define N_NODES 400   /* number of nodes in the phylo tree */
#define N_LEAVES 400   /* number of leaves in the phylo tree */
#define MAX_SPEC 400   /* max. number of species in tree */
#define N_CHILDREN 10 /* max. number of childnodes of a node in the tree
		        NOTE: increasing, slows down bw recursion!!! */

/* -------------------------------------------------------------------------- */

typedef struct gl_nodeN {
	
	
	
	int parentLab;    /* lab of parent, as read from input tree file */   
	struct gl_nodeN * p; /* parent */ 	//eher root ((NUTZLOS?..))   .. übrigens, innerhalb einer struct sind nur pointer auf gleichen struct typ erlaubt, aber keine recursive deklaration der gleichen struct								
	int c[N_CHILDREN]; /* list of labels of childnodes */
						/* -1 if no child (at leaves) */
	
	//  int c_list[N_CHILDREN]; /* list of node-labels of child nodes */
	int m;        /* number of mirs at node */
	int n;        /* number of node == label */   
	int nc;       /* number of actual childnodes */
	char o[200];   /* species/lineage name */
	int pInT;     /* position in treeN array */
	int pInP;     /* position of node in preorder */
	int seenBW;   /* already seen in bw recursion? {-1, 1} */
	int seenFW;   /* already seen in fw recursion? {-1, 1} */
	int seenO;    /* mark seen nodes when going through preorder */
	/* for the traversion */
	/* int label; */
	int P;       /* max No of miRNAs in subtree rooted at this node */
	int level;   /* distance to root */
	float score; /* score in recursions */						//ACHTUNG: nach bw recursion wird score nicht mehr gebraucht, wurde mit p(Zv) überschrieben!
	int gain;    /* number of gained miRNAs */
	int loss;    /* number of lost miRNAs */
	int gainFam; /* 1 if gene family has been invented at this node, 0 otherwise */
	int lossFam; /* 1 if gene family has been lost at this node, 0 otherwise */
		
	/* weights for gain and loss at this node */
	float gainW;
	float lossW; 

	//calculations for pf:
	
	int pfGain;
	int pfLoss;
	int pfm;
	
	int pfGainFam;
	int pfLossFam;
	
	float pfP;
	float psP;
	
	int seenBWp;   /* already seen in partition function  bw recursion? {-1, 1} (that means the top down recursion, root to leaves)*/
	
  
} gl_nodeN;

	
typedef struct collapse
{
	char o[200];			//species
	struct collapse * p; 
	
	int n;				//label
	int parentLab;		//
	int c[N_CHILDREN];	
	int nc;				//anzahl kinder
	int level;			//
	int tocallapse;
	float gainW;
	float lossW;	
} collapse;


/* -------------------------------------------------------------------------- */

#endif
