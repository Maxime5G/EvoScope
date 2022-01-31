#ifndef _COEVOL_H_
#define _COEVOL_H_

#include "tree.h"
#include <time.h>

#define MAX_MU 1000.0
#define MIN_MU 1e-7
#define MIN_GRAD 1e-7
#define MIN_C 1e-4

extern int verbose;

/* PART1: Objects initializations and input manipulations (input_output.c) */

/* Integer matrix */

typedef struct integer_matrix {
	int ncol, nrow;
	int *val;
} IntegerMatrix;

/* Sparse matrix */
typedef struct sparse_matrix {
	int *col;
	int *row;
} sparseMatrix;

typedef struct tvector_struct{
	char **tvector;
	int nbfork;
	int nvect;
	char **current_vector;
} TVector;

typedef struct events_array{
	int **inpe2;
	int *e1;
	int length;
	int *e2L;
}EventsArray;

struct EpocsData {

	Node *tree;
	int nbleaves_tot;

	Node *root;
	double *branch_lengths;
	int nbranches;
	int nbleaves;

	// Epocs
	int changed;
	int nevt;

};

struct CoevolData {

	Node *tree;
	int nbleaves_tot;

	Node *root;
	double *branch_lengths;
	int nbranches;
	int nbleaves;
	double longestbranch;
	char **mask;

	// Epics

	IntegerMatrix *chronomat;
	sparseMatrix *sparseChronomat;
	int *v_evts;
	int **e1;
	int *e1M;

	// Epocs
	int changed;
	int nevt;

    TVector theTVector;

};

/* The next structure is basically an array of trees		*/
/* If the user inputs a file with multiple trees inside		*/
/* One line = one tree and the whole structure is a forest	*/
/* If multiple files: 1 entry in the array but multiple 	*/
/* CoevolData because I have independent trees				*/

struct Trees {
	struct CoevolData *MyCoevolData; /* Forest of trees contain 1,2,...n trees in the forest -> 1,2,...,n CoevolData	*/
	double branch_lengths;	/* sum of branch lengths for each tree in the forest	*/
	double lik_ntrees;	/* ML of each tree in the forest (e.g., as output from BEAST)	*/

	// epocs
	double lik_epocs;	/* ML of epocs	*/
	double final_rates_epocs[2][4]; /* final_rate is [2][4], we have ntrees * final_rate	*/
	int IS_epocs[2];
	int imax_epocs;
	int jmax_epocs;
	int *o2;
	int nbfork;
	int *forkVector;				/* number of forks per tree - for recursion	*/
	int *maxForkVector;

	//epics
	double pvalue_of_pair; /* pvalue of pair in epics	*/

	// mode calculation
	// double weight;	/* likelihood_tree*likelihood_epocs for each tree - DONOTKEEP	*/
	// double *weights_norm;	/* normalized weights such that they sum to 1 - DONOTKEEP	*/

};

/* object initializations */

void InitData(  struct CoevolData *t );

/* Parse and read inputs */

void Read_Coevol_Data( struct CoevolData *MyCoevolData, int t, char *infile, char *outgroup, int opt_read_stdin, int *nevt, double *total_length_all_tree, int IsItEpics, int mat_type, int verbose );
void Read_Coevol_InForest_Data( struct CoevolData *MyCoevolData, int t, char *newick_string, char *outgroup, int opt_read_stdin, int *nevt, double *total_length_all_tree, int IsItEpics, int mat_type, int verbose );

struct CoevolData *ParseInput( int argc, int optind, char **argv, char *infile, int verbose, int *nevt, int *maxtree, int mat_type, int IsItEpics  );
struct Trees *ParseInputForest( int argc, int optind, char **argv, char *infile, int verbose, int *nevt, int *maxtree, int **maxtreeforest, int mat_type, int IsItEpics  );

double * ComputeLikelihoodDistribution(struct Trees *MyForest, int ntrees, int isItEpics, int evti, int evtj);

// init
char * initializeNewFile(char * pref, char scen, char *matID, int whichOutput);
char ** createEventsArray(char * myFname);
int ** createArrayEventsToRun(char *MyInputName, int *fileLength, int verbose);
EventsArray *createArrayEventsToRunEpics(char *MyInputName, int verbose);

int isDuplicate(int number1, int number2, int **storedArray, int lastPosition);

// writing in file
void writeInFile(int numOfTree, int forestLength, char *eventIds, FILE *filePointer, char **inputIds, char *matId, int e1, int e2, int treeId, int treeInForest, double pvalue, int obsPairs, int maxObs);

// free memory
void freeArray(int **a, int m);
void FreeCoevolData( struct CoevolData *t, int n, int nevts, int mat_type );

/* PART2: EPOCS - coevol_FindML and coevol_FindML_multi	*/
/* function headers */

double ML_by_NR( Node * root, char **ts, double tau[2][4], int *IS, char scenario );
double ML_by_Grad_n_NR( Node * root, char **ts, double rate[2][4], int *IS, char scenario );
double ML_by_Grad_n_NR_multi( struct CoevolData *MyEpocsData, int ntree, double rate[2][4], int *IS, char scenario_multi );
void ML_multi(struct CoevolData *MyEpocsData, int ntree, double *BestMaxlnL, double rate[2][4], double final_rate[2][4], int *IS, char scenario_multi, int *maxVector, int level, int *resVector, int *o, int *imax, int *jmax, int opt_IS, int round);

void final_state( int ev, int *IS, int *FS, int nstates);
double lnLbranch( double l, int ev, int *IS, double tau[2][4], int *FS, int nstates );
double lnLbranch_ori( double l, int ev, int IS1, int IS2, double *tau1, double *tau2, int *FS1, int *FS2  );

double ln_vraisemblance(Node *n, char *vector_evt_type, double tau[2][4], int *IS );

double minimization(Node *root, char *vect_events, double tau[2][4], int *IS, int verbose);

/* PART3: EPICS - coevol_EpicsMat, coevol_sparseMat	*/
/* fonctions de mat.c appelees en dehors */

IntegerMatrix *InitIntegerMatrix(int nrow, int ncol);
void FreeIntegerMatrix(IntegerMatrix *m);
double fill_evt_matrix(Node *n, IntegerMatrix *matrice_a_completer, double *longueurs, int nbranches);
void chrono(Node *n, int *t, IntegerMatrix *chronomat);
void chrono_old(Node *n, int *t, IntegerMatrix *chronomat);
void chrono_adj(Node *n, IntegerMatrix *chronomat);
void transpose_a_matrix(IntegerMatrix *mat);
void transpose_a_matrix_antidiagonal(IntegerMatrix *mat);
IntegerMatrix * chrono_excl(Node *n, int nbranches);
void print_matrice(IntegerMatrix *m, FILE *f);
int sumv(int *v, int lv);
int sumvP(int *v, int lv);
int nnzv(int *v, int lv);
int iprodsca(int *u, int *v, int l);
void vectmat(int *v, IntegerMatrix *m, int *pv);
void init_sTableau(int t);
void init_sTableau_t(int n, int t);
void free_sTableau(int T);
void init_sFacto(int n);
void free_sFacto(void);
void print_ivect(FILE *f, int *v, int l);

void init_slnFacto(int n);
void free_slnFacto(void);

/* partie main	*/

int *SetEvents( Node *root, int *Nevts, int verbose );
IntegerMatrix *SetMatrix( int mat_type, int nbranches, Node *root, int output_mat );
int **InitEvents( Node *root, int nbranches, int nevt, int mat_type, int **e1M, int output_vectors );
int **InitEventsAndMaskFlag( Node *root, int nbranches, int nevt, int mat_type, int **e1M, char ***mask, int output_vectors);

/* partie multinomiale */
void set_vect_classes(int *e1m, int *e2, int nbranches, int *nclasses, int **ne1m, int **ne2, double *proba, double **nproba, int t);
void set_vect_classes_mask(int *e1m, int *e2, char *maske1, char *maske2, int nbranches, int *nclasses, int **ne1m, int **ne2, double *proba, double **nproba, int t);
void generate_multinomiale(int *veck, double *vproba, int lveck, int current_niveau, int prevn, int reste,
		double fact_ratio, double *distrib, int *e1m, int sum, int t);

void generate_multinomial(int *veck, double *vproba, int lveck, int sum_e2, int rest, int level, double *distrib, int* e1m, double proba, double denom);
// int generate_multinomial(int *veck, double *vproba, int lveck, int sum_e2, int rest, int level, double *distrib, int* e1m, double proba, double denom, clock_t *start_time, int res);
long long binomial(int n, int r);
void normalize_a_vector (double *v, int lv);
void generate_approached_multinomial(int *veck, double *vproba, int lveck, int sum_e2, int rest, int level, double *distrib, int* e1m, double proba, double denom);
/* coevol_sparseMat.c */

sparseMatrix *convertToSparseMatrix(int nrow, IntegerMatrix *input_mat);
void vectSparseMat (int *v, sparseMatrix *m, int lv, int *pv);
void vectSparseMatAndMask (int *v, sparseMatrix *m, int lv, int *pv, char *maske1, char *maske2);
void print_sparse_matrix(sparseMatrix *m, FILE *f);

#endif
