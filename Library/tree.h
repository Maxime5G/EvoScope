#ifndef _LIBTREE_
#define _LIBTREE_

#include <stdio.h>
#include <stdlib.h>


/**
\file 	tree.h
\brief 	Headers for reading and outputting generic trees in newick format
\author g achaz and madamesophie
\date 	June 28 2012

**/


/**
\struct multinode
\brief  same as node with multi descendants  representation
**/

typedef struct mymultinode {

	/* OLD JOMPO char  id; changed to int on 28 juin 2017 */
	int id;					/**< number id for the node **/
	double time;                            /**< evaluation of time **/

	struct mymultinode *anc;                /**<  only one ancestor **/
	struct mymultinode **descs; 	        /**<  aray of mymultinodes **/
	int nbdesc; 				/**< nbr of descendants (size of descs array)**/

	char *name; 				/**<  nickname (only for leaves) **/

	int *SeqInLeaves;                       /**<   usefull array **/
	int nbMuts;                             /**< nbr of possible mutations **/

	int which_specie;                       /**< nbr to find easily the specie in infoseq X **/
	char *possMut;                          /**< description of possibles allels **/
	double misc;                            /**<misc value to store anything**/

	int nleaf;                              /* nombre de feuilles */
	int nbranches;

	// epics / epocs
	char  evtype;                           /** type of event 0,1,2,3,4,5 */
	int   nevt;                             /* nombre d'evenements dans un noeud */
	int   *evt;                             /* array with all events */

} Node;


/*
	Memory functions
*/
Node * InitMultiTree(  int nleaves );
void FreeMultiTree(struct mymultinode *myNodes, int nleaves);
void FreeTree( Node  * Tree, int n_nodes );
void FreeNode( Node  * n );

/*
	Read From File
*/
Node *ReadTreeFile( char *infile, int *nleaves );                                  /* Embedding function memory included from a file to a Tree structure */
char *readFileNwck(FILE *f,int *nbleaves);                                         /* from a FILE *f to a string */
char **readFileMultiNewick(char *treefile, int opt_read_stdin, int **nbleaves, int *ntrees);	/* from a FILE *f to an array of strings	*/
int lit_newickmyMulti(char *chaine,int la_pos,int *which_node, Node *myNodes);     /* from a string to the structure */
void clean_end(char *str);														   /* removes last value if exists*/
Node *read_tree( char * treefile, int opt_read_stdin, int *nbleaves, char *outgroup ); /* read a tree from the stdin or from a file	*/
Node *read_tree_forest( char * newick_string, int *nbleaves, char *outgroup, int verbose ); /* read a tree from a newick string	*/

/*
	Ouput Tree
*/
void PrintTree( Node *Tree, char *outfile, short opt_FromRoot );                   /* if outfile is NULL --> print to stdout ; if opt_FromRoot=1, track up to root */
void multitreeoutNck( struct mymultinode *pptree, FILE *f );
char *countchar(char *ch);
void print_arbre(Node *n, FILE *f);
void PrintNode( Node *n, char *mssg );


/* utils for tree */
int count_leaf(Node *n);
int count_nodes(Node *n);

int count_evt(Node *n);
int verif_evt(Node *root, int nevt, int verbose);
int count_max_event(Node *n, int nevt, int verbose);
int verif_nevt(Node *n, int nevt, int verbose);
int set_evt_type_count_fork(Node *n, int nevt, int evti, int evtj, int verbose, int output_state_vector);
int verif_polymorphism(Node *n, int *v_evt, int nevt);
char **genere_all_vector_types(Node *n, int last_branch, int nbfork, char **ts, int *nvect);
char **genere_dual_vector_types(Node *n, int last_branch, int nbfork , char **ts, int *nvect);
void add_random_vector_types(char **ts, int last_branch, int nbfork, int maximumTvectors);
void add_random_vector_types_ok(char **ts, int last_branch, int maximumTvectors);
void print_vector_types(char **ts, int len, FILE *fout);
double ln_vraisemblance(Node *n, char *vector_evt_type, double tau[2][4], int *IS );
void change_names_with_events(Node *n, int *v_evt, int nevt);
int count_branches(Node *n);
void set_e1_vectors(Node *n, int **e1);
void set_e1_vectors_mask(Node *n, int **e1, char **mask);

double *ComputeBranchLengths( Node *root, int nbranches, int verbose, double *total_length_all_tree, int t, int runIsEpics);				// compute branch lengths
void SetMinimumBranchLength( Node *n, double min_length );   // set all <= min_length to min_length
void NormalizeTree( Node *n, double Ltot );                  // normalize tree to 1
void NormalizeTreeAndLongestBranch( Node *n, double Ltot, double *longest);	// normalize to 1 and extract the longest branch

double TotalLength( Node * mynode );   /* Sum all branch lengths from Root */
double SumLength( Node * mynode );     /* Sum all branch lengths from mynode included */
double SubTreeLength( Node * mynode );  /* Sum all branch lengths from mynode excluded */


/*
Utils for sofiz maybe not here
*/
void BuildSeqInLeaves(Node *pptree,int N);
void recopie(int *dest, int* src, int N);
int *joinLeaves(int *s1, int* s2, int N);
void initLeaves(Node *pptree,int N);


int isAnc(Node * Test, Node * Ancestor);
int HowManyAnc(Node  *Test);

int binarise (Node *myNodes,int nbnodes);                        /* make a binary tree from a tree with multiple desc */
void SplitLeaves( Node *pptree, int *ind, int ind_size );        /* track each leaves in the subtree and split leaves with several sequences into 1 leaf for each seq */

void connect( Node * WholeTree, int anc, int desc1, int desc2);  /* need the whole Tree Pointer --global version-- */
void Connect( Node * AncNode , Node * DescNode , float time );   /* new version --only use pointer --local version-- */


void ReRoot( Node * TargetNode );                                    /* reset root to the branch going up from TargetNode */
void ReRootOutgroup( char *outgroup, Node *init_tree, int nnodes );  /* reroot using an outgroup string */

Node *SetRootEpics( Node *tree, char *outgroup, int output_tree);
void ReLabelNodeId( Node *node, int *id );

/*
	In findnode.c
*/

Node * FindRoot(   Node *mynode );                               /* reccursive faster version of findHead() */
Node * findHead( Node *myNode, int nnodes );                      /* deprecated */
Node * FindLeaf(  char *name, Node *NodeArray, int nnodes );
Node * FindMRCA(  Node *n1, Node *n2  );
Node * FindMRCA_char( char *request, Node *NodeArray, int nnodes );
Node *FindMRCA_intarray( int *array, int array_size, Node *NodeArray, int nnodes );

/*
	In developpement
*/

int *GenerateLeavesArray( char *leaves_str, Node *node_array, int nnodes, char opt_invert, int *size_array);

#endif
