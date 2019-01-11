#ifndef READ_H
#define READ_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "tree.h"
#include <setjmp.h>


typedef struct gl_arguments {

  char * infile;  /* alignment */
  char * treefile; /* tree (see example.tree.dat) */
  char * type;    /* type of output values */ //kann hier schon auf default all gesetzt werden:
  char * treefileNew;/* NEWICK tree*/
  char * outfile; /* outfile name [OPTIONAL] */
  char * outfileFlag;	//the additional names for the outfiles
  char * outfileFlag_b;
  char * psfile;  /* postscript outfile name [OPTIONAL] */
  char * psfileFlag; //the additional name for the .ps file, no second file needed for -b because both are in one file
  char * collectfile; /* file to collect results over several calls */
  char * scores;  /* array of scores for inner nodes */
  char * directory;
  int clw; //
  int stk; //
  int z; //calculate just PF
  int b; // calculate both, PF and PS
  int P; //k bei höchsten P wählen (PF backward)
  int C; //collapsen der nodes unterdrücken?
  int lossBT;     /* 0|1 that activates loss-backtracing at the end */
  int collaps;
  int seenBWp;

} gl_arguments;


struct gl_arguments gl_readArguments(int argc, char *argv[]);

/* init functions */
struct gl_arguments gl_initArguments();

/* print functions */
struct gl_arguments gl_printArguments(gl_arguments ga);

void usage();

void version();

/* free functions */
void gl_freeArguments(gl_arguments ga);

/* data reading functions */
void gl_readDATA(struct gl_arguments ga);

/* reads the tree file and builds treeorder and tree[] itself */
void gl_buildTree(char *filename);

/* fills array of treePosN with positions of labels */
void gl_label2pos();

/* reads the tree in preorder and stores the labels in treeOrderN[] */
int gl_preorder(struct gl_nodeN node, int o);

char *strdup(const char *src); 

/* splits a string str2split by a character delim, returns an array of char * */
char** gl_splitString(char * str2split, const char delim, int * counter);

/* sets the P at each node in treeN */
void gl_setPatInnerNodes();

/* recursively get the Pmax (max number of genes in subtree rooted at v) */
int gl_getPv(int label_v);

/* reads the summarized tree and stores values in treeN */
void gl_readSummaryTreeN(char * file) ;

void collapstest(char *filename, struct gl_arguments *ga);

struct gl_arguments getFilenameExtension(struct gl_arguments ga);

void set_weights(struct gl_arguments ga) ;

#endif
