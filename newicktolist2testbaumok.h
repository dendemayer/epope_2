#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <setjmp.h>

typedef struct
	{	char name[100];			//nodeNEW name of spezies, set in findleaves
		char parent[100];		//parent name of nodeNEW, set in findleaves
		int par;				//lab number of the parent
		char children[100];		//string array of lab numbers of all children of a nodeNEW
		int d2root;				//d2root is set in print fkt, with recursive fkt findroot
		int lab;				//counting up lab numbers while setting leavenames //Problem, diese müssen einer preorder entsprechen, damit sie dann korrekt in epope eingelesen werden
	}nodeNEW;
	
	
void newick(char * filein);
void test(char *filein);	//counting opening and closing brackets
int count(char *filein);	//number of characters in newickarray
int countsp(char *array, int l);	//number of nodeNEWs
void newicktoarray(char *filein, int *laenge, char array[]);	//reading from file, creating the newickarray, distances are cut away
void arrayplus(int *laenge, char array[]);
void findleaves(nodeNEW list[], char array[], int newarrlen, char *filein, int spcount);	//setting up lab numbers, name, parent in nodeNEW structs
//nodeNEW list[] -> array aus struct nodeNEW der länge spcount()
//char array[] -> newickarray
//int newarrlen -> laenge des newickarrays
nodeNEW writeleave(char leavename[],char dadname[], int lab,int d2root, nodeNEW leave);		//fkt to write the nodeNEW structs
//char leavename[] -> string für speziesnamen
//char dadname[] -> string für vater
//int lab -> lab nr. soll global je aufruf von writeleave erhöht werden
//int d2root -> wird vor jeden cut array aufruf erhöht, muss dann aber genau umgedreht werden...
//nodeNEW leave -> aus findleaves an writeleaves übergebener nodeNEW aus nodeNEW list[]
void cutarray(char newickarray[], nodeNEW list[], char *filein, int spcount);	//gleiches vorgehen wie bei findleaves, inhalte werden jedoch gelöscht
												//nodeNEW list aus findleaves mit übergeben ums wieder zurückzugeben
void print(nodeNEW list[], int d2root, char *filein, int l); //bearbeiten von 												

void findroot(int rootlab, int root, nodeNEW list[],char *filein, int spcount);

char* setroot(char *newickarray, int  pl);
