#ifndef CALC_H
#define CALC_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "read.h"
#include "tree.h"



/* do forward and backward recursion */
void gl_calc(struct gl_arguments ga);

/* delta function that returns corresponding weight */
float gl_delta(int i, int k, int node);

/* find index of node that is last common ancestor of current miRNA family 
   assign min and max number of miRNAs in current alignment */
int gl_getLCA(int *kMin, int *kMax);

/* prints the tree in all its elements as POSTSCRIPT */

void gl_printTreePSN(int kMax, int lca, struct gl_arguments ga);
/* prints the treeN in all its elements */
void gl_printTreeN(struct gl_arguments ga);

/* prints score array */
void gl_printS(float **S, int m, int n);

/* old get lca that uses preorder of tree
   NOTE: labels of the tree were not in preorder !
 */
 
void pf(int kMax, struct gl_arguments ga, int lca);
void pfb(int kMax,int nodesN, float **Z,struct gl_arguments ga, int lca );
void pf_with_z(struct gl_arguments ga);
void pfgl(float **Pku, int kMax,int nodesN, int lca, struct gl_arguments ga);

void set_weights(struct gl_arguments ga) ;

int gl_OLDgetLCA();
#endif
