#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tree.h"
#include "coevol.h"

#include <sys/stat.h> 		/* for the creation of directories and files */
#include <sys/types.h> 		/* for the creation of directories and files */

/**
\file 	input_output.c
\brief 	initialization, input and output options for epics and epocs
\author Maxime Godfroid
\date 	December 12 2020
**/

/* Objects initialization */

/*	Initializing the EpicsData object	*/
void InitData(  struct CoevolData *t ){

	t->tree=NULL;
	t->nbleaves_tot=-1;

	t->root=NULL;
	t->branch_lengths=NULL;
	t->nbranches=-1;
	t->nbleaves=-1;
	t->longestbranch=-1;
	t->mask=NULL;

	t->chronomat=NULL;
	t->sparseChronomat=NULL;
	t->v_evts=NULL;
	t->e1=NULL;
	t->e1M=NULL;

    t->changed=-1;
    t->nevt=-1;
    t->theTVector.tvector=NULL;
    t->theTVector.nbfork=-1;
    t->theTVector.nvect=-1;
    t->theTVector.current_vector=NULL;

}

void InitForestData(struct Trees *t){
	t->branch_lengths=0;
	t->lik_ntrees=0;
	t->lik_epocs=0;
	t->imax_epocs=0;
	t->jmax_epocs=0;
	t->pvalue_of_pair=0;
	t->o2=NULL;
	t->nbfork=0;
	t->forkVector=NULL;
	t->maxForkVector=NULL;
}

/* This function determines the chronology matrix to generate and store in chronomat. */
/* INPUTS: the mat_type (determined from the command line argument), the number of    */
/* branches in the tree, the root node and the verbose option.						  */
/* The idea is to generate the A, S, X, P and B matrices depending on the option and  */
/* store it in chronomat, from the MyEpicsData structure. Special case is 0, the 	  */
/* identity matrix where nothing has to be done. */
IntegerMatrix *SetMatrix( int mat_type, int nbranches, Node *root, int output_mat ){

	IntegerMatrix *chronomat=NULL;
	int *chrono_vect = NULL;
	int i;

  /* ---------------------------- Chronology matrix -------------------------------- */
  /* Initialization of the chronology matrix - determines which matrix will be used. */

  switch(mat_type) {

  	case 0:
		/* Matrix is identity, so			*/
		/* no need of matrix, only scalar product e1.e2	*/
		fprintf(stderr, "Matrix Id will be used (co-occurences)\n");
		break;

	case 1:
		fprintf(stderr, "Matrix S will be used (chronologies)\n");
		chronomat = InitIntegerMatrix(nbranches,nbranches);

		/* Allocation of the chrono_vect vector.						   */
		/* This vector is important in the chrono function to determine    */
		/* whether a node is a descendant, of any order, of any given node */
		if ((chrono_vect = calloc(nbranches,sizeof(int))) == NULL)
			fprintf(stderr, "memory error for time vectors, exit\n"),	exit(3);

		/* Filling the chronology matrix 								   */
		chrono(root,chrono_vect, chronomat);

		/* No further need of chrono_vect - it is freed 				   */
		free(chrono_vect);

		if (output_mat) {
			fprintf(stderr,"S matrix :\n");
			print_matrice(chronomat, stderr);
		}
		break;

  	case 2:
		fprintf(stderr, "Matrix Id+S will be used (co-occurences plus chronology)\n");
		chronomat = InitIntegerMatrix(nbranches,nbranches);

		/* Allocation of the chrono_vect vector.						   */
		/* This vector is important in the chrono function to determine    */
		/* whether a node is a descendant, of any order, of any given node */
		if ((chrono_vect = calloc(nbranches,sizeof(int))) == NULL) {
			fprintf(stderr, "emory error for time vectors, exit\n");
			exit(1);
		}

		/* Filling the chronology matrix 								   */
		chrono(root,chrono_vect, chronomat);

		/* No further need of chrono_vect - it is freed 				   */
		free(chrono_vect);

		/* add co-occurences (i.e. add matrix Identity) */
		for (i = 0; i < nbranches; i++)
			chronomat->val[i+i*chronomat->nrow] += 1;

		if (output_mat) {
			fprintf(stderr,"S+Id matrix :\n");
			print_matrice(chronomat, stderr);
		}
		break;

	case 3:
        fprintf(stderr, "Matrix A+Id will be used (co-occurences plus adjacent pairs)\n");
		chronomat = InitIntegerMatrix(nbranches,nbranches);

		/* Filling the chronology matrix 								   */
		chrono_adj(root, chronomat);

        /* add co-occurences (i.e. add matrix Identity) 				   */
		for (i = 0; i < nbranches; i++)
			chronomat->val[i+i*chronomat->nrow] += 1;

        if (output_mat) {
			fprintf(stderr,"A+Id matrix :\n");
			print_matrice(chronomat, stderr);
		}

		break;

    case 4:
		fprintf(stderr, "Matrix P+Id will be used (co-occurences plus precedence)\n");
		chronomat = InitIntegerMatrix(nbranches,nbranches);

		/* Allocation of the chrono_vect vector.						   */
		/* This vector is important in the chrono function to determine    */
		/* whether a node is a descendant, of any order, of any given node */
		if ((chrono_vect = calloc(nbranches,sizeof(int))) == NULL) {
			fprintf(stderr, "Memory error for time vectors, exit\n");
			exit(1);
		}

		/* Filling the chronology matrix 								   */
		chrono(root,chrono_vect, chronomat);

		/* No further need of chrono_vect - it is freed 				   */
		free(chrono_vect);

		/* add co-occurences (i.e. add matrix Identity) 				   */
		for (i = 0; i < nbranches; i++)
			chronomat->val[i+i*chronomat->nrow] += 1;

        /* now transposing the matrix 						   			   */

        transpose_a_matrix(chronomat);

		if (output_mat) {
			fprintf(stderr,"P+Id matrix :\n");
			print_matrice(chronomat, stderr);
		}
		break;

      case 5:
        fprintf(stderr, "Matrix of exclusion will be used\n");

		/* Filling the chronology matrix 								   */
		/* We don't need to initialize an empty matrix, it is performed in */
		/* the function and is returned by it.							   */

        chronomat = chrono_excl(root, nbranches);

        if (output_mat) {
			fprintf(stderr,"X matrix :\n");
			print_matrice(chronomat, stderr);
		}
		break;

  	default:
		fprintf(stderr, "Should not happen: matrix type is %d (not 0,1,2,3,4 or 5)\n", mat_type);
		break;
  }

	return chronomat;
}

/* Function to read the events and check a couple of informations: 		*/
/* INPUTS: the root node, the number of events and the verbose option   */
/* - nevts is counting the number of events on each branch. 			*/
/* Exits the program if numbers are different. 							*/
/* If there are unspecified events -> verif_evt will check the tree 	*/
/* recursively and where events are missing will add 0's. 				*/
/* finally, verif_polymorphism will check if events are polymorphic.	*/
/* if yes, will store a -2 in the array v_event, 0 otherwise.			*/
int *SetEvents( Node *root, int *Nevts, int verbose ){

	int *v_evt, nevts;
	int n_polymorph=0;
	int changed = 0;


	/* count all events on the tree */
	nevts = count_evt(root);

	if ( *Nevts != -1 && *Nevts != nevts)
		fprintf(stderr, "All trees must have the same # of event types, exciting\n"), exit(3);
	else
		*Nevts = nevts;


	/* Check and eventually add missing zeros (unspecified events) */
	if (nevts > 0)
	{
		changed = verif_evt(root, nevts, verbose);
		if ( changed && verbose )
			fprintf(stderr,"-> %d nodes have been initialised with list of %d events (filled with 0)\n",changed, nevts);
	}


	/* check if events are polymorphic		*/
	/* -1: polymorphism, other: polymorphism	*/

	if ((v_evt=(int *)calloc((size_t)nevts, sizeof(int)))==NULL)
		fprintf(stderr,"Not enough memory for vector of polymorphism of events\n"),exit(2);

	n_polymorph = verif_polymorphism(root, v_evt, nevts);

	if (verbose)
		fprintf(stderr,"Number of event sites: %d from which %d are polymorph\n",nevts,n_polymorph);

	return v_evt;
}

/* This function intializes and fills the e1 vectors (i.e. the mutation vectors) per branch 				   */
/* e1 is a pointer of pointer of ints, where rows are the number of events and the columns are the branches.   */
/* The function also initializes the e1M array, which will be used to store the multiplication of e1 by the    */
/* chronology matrix chosen.										 										   */
/* INPUTS: the root node, the number of branches, the number of events, the chronology matrix, the pointer     */
/* to e1M and the verbose option to print the vectors. 														   */
/* As said, the function will initialize the e1 array and then will fill it with the function set_e1_vectors.  */
/* The function is straightforward and just reads the tree recursively and stores the state of each event at   */
/* branch. 																									   */
int **InitEvents( Node *root, int nbranches, int nevt, int mat_type, int **e1M, int output_vectors ){

	int **e1, i, j;

	/* We draw the e1/e2/... vectors from the tree as an array of vector e[i]. */
	/* The number of events is nevt.										   */

	if ((e1=(int **)malloc((size_t)nevt*sizeof(int *)))== NULL) {
		fprintf(stderr,"Not enough memory for e1 vectors\n");
	exit(1);
	}
	for (i = 0; i < nevt; i++) {
	  if ((e1[i]=(int *)malloc((size_t)nbranches*sizeof(int)))== NULL) {
		  fprintf(stderr,"Not enough memory for e1 vectors\n");
		  exit(1);
	  }
	}

	/* if Matrice is Id, no e1M is needed	(e1M is Matrix.e1)	*/
	/* as Id.e1 = e1 											*/
	if (mat_type != 10) {
	  if (((*e1M)=(int *)malloc((size_t)nbranches*sizeof(int)))== NULL) {
		  fprintf(stderr,"Not enough memory for e1M vectors\n");
		  exit(1);
	  }
	}

	/* fill the event vectors	*/
	set_e1_vectors(root, e1);

	/* print the vectors		*/
	if (output_vectors) {
	  for (i = 0; i < nevt; i++) {
		  fprintf(stderr,"vector %2d (sum=%2d): ",i,sumv(e1[i],nbranches));
		  for (j = 0; j < nbranches; j++) {
			  fprintf(stderr,"%1d ",e1[i][j]);
		  }
		  fprintf(stderr,"\n");
	  }
	}

	return e1;
}

/* WIP: same as above but with the mask vector	*/

int **InitEventsAndMaskFlag( Node *root, int nbranches, int nevt, int mat_type, int **e1M, char **mask, int output_vectors){

	int **e1, i, j;

	/* We draw the e1/e2/... vectors from the tree as an array of vector e[i]. */
	/* The number of events is nevt.										   */

	if ((e1=(int **)malloc((size_t)nevt*sizeof(int *)))== NULL) {
		fprintf(stderr,"Not enough memory for e1 vectors\n");
	exit(1);
	}
	for (i = 0; i < nevt; i++) {
	  if ((e1[i]=(int *)malloc((size_t)nbranches*sizeof(int)))== NULL) {
		  fprintf(stderr,"Not enough memory for e1 vectors\n");
		  exit(1);
	  }
	}

	/* if Matrice is Id, no e1M is needed	(e1M is Matrix.e1)	*/
	/* as Id.e1 = e1 											*/
	if (mat_type != 10) {
	  if (((*e1M)=(int *)malloc((size_t)nbranches*sizeof(int)))== NULL) {
		  fprintf(stderr,"Not enough memory for e1M vectors\n");
		  exit(1);
	  }
	}

	if ((mask = (char **)malloc((size_t)nevt*sizeof(char *))) == NULL){
		fprintf(stderr,"Not enough memory for mask vectors\n");
		exit(1);
	}
	for (i = 0; i < nevt; i++) {
	  if ((mask[i]=(char *)malloc((size_t)nbranches*sizeof(char)))== NULL) {
		  fprintf(stderr,"Not enough memory for e1 vectors\n");
		  exit(1);
	  }
	}

	/* fill the event vectors	*/
	set_e1_vectors_mask(root, e1, mask);

	/* print the vectors		*/
	if (output_vectors) {
	  for (i = 0; i < nevt; i++) {
		  fprintf(stderr,"vector %2d (sum=%2d): ",i,sumv(e1[i],nbranches));
		  for (j = 0; j < nbranches; j++) {
			  fprintf(stderr,"%1d ",e1[i][j]);
		  }
		  fprintf(stderr,"\n");
	  }
	  for (i = 0; i < nevt; i++) {
		  fprintf(stderr,"mask vector %d: ",i);
		  for (j = 0; j < nbranches; j++){
			  fprintf(stderr,"%1d ",mask[i][j]);
		  }
		  fprintf(stderr,"\n");
	  }
	}

	return e1;
}

/* Read Coevol Data */

/*	Function to read the coevol data from a given set of parameters, and store the informations in the MyCoevolData object	*/
/*	Make use of: read_tree_forest, SetRootEpics, count_leaf, count_branches, ComputeBranchLengths, count_max_event, verif_nevt	*/
void Read_Coevol_InForest_Data( struct CoevolData *MyCoevolData, int t, char *newick_string, char *outgroup, int opt_read_stdin, int *nevt, double *total_length_all_tree, int IsItEpics, int mat_type, int verbose ){

	MyCoevolData[t].nbleaves=0;
	MyCoevolData[t].nbranches=0;
	MyCoevolData[t].longestbranch=0;

	MyCoevolData[t].tree = read_tree_forest( newick_string, &(MyCoevolData[t].nbleaves_tot), outgroup );
	MyCoevolData[t].root = SetRootEpics( MyCoevolData[t].tree, outgroup, verbose );

	MyCoevolData[t].root->anc=NULL;
	MyCoevolData[t].nbleaves = count_leaf( MyCoevolData[t].root );
	MyCoevolData[t].nbranches = count_branches( MyCoevolData[t].root );

	MyCoevolData[t].branch_lengths = ComputeBranchLengths( MyCoevolData[t].root, MyCoevolData[t].nbranches, verbose, total_length_all_tree, t, IsItEpics);

    if (IsItEpics){
		MyCoevolData[t].e1M = NULL;
		MyCoevolData[t].v_evts = SetEvents( MyCoevolData[t].root, nevt, verbose);
        MyCoevolData[t].chronomat = SetMatrix( mat_type, MyCoevolData[t].nbranches, MyCoevolData[t].root, verbose ); // mat_type 0:Id, 1:S, 2:Both
		if (mat_type != 0)
			MyCoevolData[t].sparseChronomat = convertToSparseMatrix(MyCoevolData[t].nbranches, MyCoevolData[t].chronomat); // convert chronomat to sparse format
		MyCoevolData[t].e1 = InitEvents( MyCoevolData[t].root, MyCoevolData[t].nbranches, *nevt, mat_type, &MyCoevolData[t].e1M, verbose );
		// MyCoevolData[t].e1 = InitEventsAndMaskFlag( MyCoevolData[t].root, MyCoevolData[t].nbranches, *nevt, mat_type, &MyCoevolData[t].e1M, MyCoevolData[t].mask, verbose );
    }else{

		MyCoevolData[t].nevt = count_max_event(MyCoevolData[t].root, 0, verbose);
		MyCoevolData[t].changed = verif_nevt(MyCoevolData[t].root, MyCoevolData[t].nevt, verbose);
    }

}

/* Parse the input and store information in Trees object	*/
/* The Forest object contains n CoevolData depending on the number of trees in each multi newick	*/
/* Multiple possibilities	*/

struct Trees *ParseInputForest( int argc, int optind, char **argv, char *infile, int verbose, int *nevt, int *maxtree, int **maxtreeforest, int mat_type, int IsItEpics  ){

	struct Trees *MyTrees=NULL;

	int usage=argc-optind;
	char *outgroup;
	FILE *f; 				// input file (line = treefile \t outgroup)
	FILE *f2; 				// tree file
	FILE *flikelihoods;		// likelihoods file - if a forest is inputted
	int c; 					// character
	int t=0; 				// number of trees
	int lmax=0,
		l=0;
	double total_branch_length_all_tree=0.0;
	double normalized_total_branch_length=0.0;

	char *filename=NULL,	// input filename
	     *foutgroup=NULL,	// outgroup string
		 *lkfilename=NULL;	// filename for the likelihoods

	*maxtree=0;

	char **tree_strings=NULL; 	// storing the newick strings
	int *nbleaves=NULL;
	int trees_in_forest=0; 		// counting the number of trees in each forest
	int max_trees_in_forest=0; 	// saving the maximum number of trees in all forests to allocate memory
	int t2;						// iterating over all trees in the forest

	switch( usage ){

		/* if the user input a tree file with -f all_trees_file: open the files, count the trees and store all trees	*/
		case 0:
			if ( ! infile ) fprintf(stderr, "Use either 'epocs -f all_trees_file' or 'epocs outgroup [tree_file]', exiting\n"),exit(5);
			f = fopen(infile, "r");
			if(!f)fprintf(stderr, "ParseInput: cannot open all_trees_file, exiting\n" ), exit(4);

			while ( (c=fgetc(f)) != EOF )
				if(c!='\n'){
					l=1;
					while ( (c=fgetc(f)) != '\n' && c != EOF && l++);
					(*maxtree)++;
					lmax=(l>lmax)?l:lmax;
				}
			// if(verbose)fprintf(stderr, "ParseInput: There are %d trees in all_trees_file\n", *maxtree);
			rewind(f);


			filename=(char *)malloc( (size_t) sizeof(char)*lmax);
			foutgroup=(char *)malloc((size_t) sizeof(char)*lmax);
			(*maxtreeforest)=calloc((*maxtree), sizeof(int));

			/* if maxtree == 1: I have a forest (otherwise shouldn't use the -f option)	*/
			/* therefore, I should have a 1-column file containing the tree likelihoods	*/
			if (*maxtree == 1){
				lkfilename=(char *)malloc((size_t) sizeof(char)*lmax);
			}


			/* opening the treefiles listed in the input file and then counting characters, number of trees etc	*/
			for (t=0; t<*maxtree; t++)
			{
				trees_in_forest=0;
				if (*maxtree == 1){
					fscanf(f,"%s%s%s",filename,foutgroup,lkfilename);
					if(verbose){printf("/***\n\ttreefile: %s outgroup: %s likelihoodfile: %s\n***/\n", filename,foutgroup, lkfilename);}
				}
				else{
					fscanf(f,"%s%s",filename,foutgroup);
					if(verbose)printf("/***\n\ttreefile: %s outgroup: %s\n***/\n", filename,foutgroup);
				}

				f2=fopen(filename, "r");
				if(!f2)fprintf(stderr, "ParseInput: cannot open newick file, exiting\n" ), exit(4);
				while ( (c=fgetc(f2)) != EOF )
					if(c!='\n'){
						l=1;
						while ( (c=fgetc(f2)) != '\n' && c != EOF && l++);
						trees_in_forest++;
						lmax=(l>lmax)?l:lmax;
					}
				max_trees_in_forest=(trees_in_forest>max_trees_in_forest)?trees_in_forest:max_trees_in_forest;
			}
			fclose(f2);

			if(verbose)fprintf(stderr, "ParseInput: There are %d %s\n", *maxtree>1?*maxtree:max_trees_in_forest, *maxtree>1?"independent trees":"trees in forest");

			/* Now that I know everything I need to know, I can allocate memory for each of the struct	*/

			MyTrees = (struct Trees *) malloc (max_trees_in_forest*sizeof(struct Trees));
			if(!MyTrees)fprintf(stderr, "ParseInput: cannot allocate MyForest, exiting\n"),exit(3);

			init_sTableau(*maxtree);

			rewind(f);

			for (t2=0; t2<max_trees_in_forest; t2++){
				InitForestData( MyTrees+t2);
				MyTrees[t2].MyCoevolData = (struct CoevolData *) malloc((*maxtree) * sizeof(struct CoevolData));
			}

			// HERE: read lnLik file
			if (*maxtree==1){
				flikelihoods=fopen(lkfilename, "r");
				if (!flikelihoods)fprintf(stderr, "ParseInput: cannot open likelihoods file, exiting\n"), exit(4);

				// For each tree in the forest, I'm storing the likelihood of the tree in lik_ntrees
				for (t2=0; t2<max_trees_in_forest; t2++){
					if(fscanf(flikelihoods, "%lf", &MyTrees[t2].lik_ntrees)!=1){fprintf(stderr, "Error reading entries in the likelihood file - probably not a float!\n"), exit(3);}
				}
				fclose(flikelihoods);
				free(lkfilename);
			}

			for(t=0;t<*maxtree; t++)
			{

				fscanf(f,"%s%s",filename,foutgroup);
				tree_strings = readFileMultiNewick(filename, 0, &nbleaves, &(*maxtreeforest)[t]);
				for (t2=0; t2<(*maxtreeforest)[t]; t2++){
					InitData(MyTrees[t2].MyCoevolData+t);

					MyTrees[t2].MyCoevolData[t].nbleaves_tot=nbleaves[0];
					if(foutgroup[0]=='"' && foutgroup[1]=='"')
						Read_Coevol_InForest_Data( MyTrees[t2].MyCoevolData, t, tree_strings[t2], NULL, 0, nevt, &MyTrees[t2].branch_lengths, IsItEpics, mat_type, verbose );
					else
						Read_Coevol_InForest_Data( MyTrees[t2].MyCoevolData, t, tree_strings[t2], foutgroup, 0, nevt, &MyTrees[t2].branch_lengths, IsItEpics, mat_type, verbose );
				}
				nbleaves=NULL;
			}
			fclose(f);

			free(filename);
			free(foutgroup);

			if (verbose)printf("total length branch tree = %f\n", total_branch_length_all_tree);

			if (!IsItEpics){
				if (max_trees_in_forest > 1){
					for (t=0; t<*maxtree; t++){
						for (t2=0; t2<(*maxtreeforest)[t]; t2++){
							normalized_total_branch_length = SubTreeLength(MyTrees[t2].MyCoevolData[t].root);
							NormalizeTree( MyTrees[t2].MyCoevolData[t].root, normalized_total_branch_length );
							SetMinimumBranchLength( MyTrees[t2].MyCoevolData[t].root, 1e-4 );
							normalized_total_branch_length = SubTreeLength(MyTrees[t2].MyCoevolData[t].root);
							NormalizeTreeAndLongestBranch( MyTrees[t2].MyCoevolData[t].root, normalized_total_branch_length, &MyTrees[t2].MyCoevolData[t].longestbranch );
						}
					}
				}
				if (max_trees_in_forest == 1){
					for(t=0;t<*maxtree; t++){
	    				NormalizeTree( MyTrees[0].MyCoevolData[t].root, MyTrees[0].branch_lengths );   // rescale branches to sum to 1 /!\ remove negative values
	    				SetMinimumBranchLength( MyTrees[0].MyCoevolData[t].root, 1e-4 );         // set all null branches to 1e-4
	    				normalized_total_branch_length+=SubTreeLength(MyTrees[0].MyCoevolData[t].root);	// recompute the total branch lengths after normalization
	    			}
	    			for(t=0;t<*maxtree; t++){
						NormalizeTreeAndLongestBranch( MyTrees[0].MyCoevolData[t].root, normalized_total_branch_length, &MyTrees[0].MyCoevolData[t].longestbranch );   // and re-normalize to 1
	    			}
				}
            }

			break;

		case 1:
		case 2:

			if ( infile ) fprintf(stderr, "Use either 'epocs -f all_trees_file' or 'epocs outgroup [tree_file]', exiting\n"),exit(5);

			outgroup = argv[optind];
			infile = (usage>1)?argv[optind+1]:NULL;

			*maxtree=1;
			init_sTableau(1);

			(*maxtreeforest)=calloc((*maxtree), sizeof(int));

			tree_strings = readFileMultiNewick(infile, 0, &nbleaves, maxtreeforest[0]);

			MyTrees = (struct Trees *) malloc ((*maxtreeforest)[0]*sizeof(struct Trees));
			if(!MyTrees)fprintf(stderr, "ParseInput: cannot allocate MyTrees, exiting\n"),exit(3);

			for (t2=0; t2<(*maxtreeforest)[0]; t2++){
				InitForestData( MyTrees+t2);
				MyTrees[t2].MyCoevolData = (struct CoevolData *) malloc(1 * sizeof(struct CoevolData));

				InitData(MyTrees[t2].MyCoevolData);
				MyTrees[t2].MyCoevolData[0].nbleaves_tot=nbleaves[0];
				Read_Coevol_InForest_Data( MyTrees[t2].MyCoevolData, 0, tree_strings[t2], outgroup, 0, nevt, &MyTrees[t2].branch_lengths, IsItEpics, mat_type, verbose );
				if (!IsItEpics){
					normalized_total_branch_length = SubTreeLength(MyTrees[t2].MyCoevolData[0].root);
					NormalizeTree( MyTrees[t2].MyCoevolData[0].root, normalized_total_branch_length );
					SetMinimumBranchLength( MyTrees[t2].MyCoevolData[0].root, 1e-4 );
					normalized_total_branch_length = SubTreeLength(MyTrees[t2].MyCoevolData[0].root);
					NormalizeTreeAndLongestBranch( MyTrees[t2].MyCoevolData[0].root, normalized_total_branch_length, &MyTrees[t2].MyCoevolData[0].longestbranch );
				}
			}

			break;

		default:
			fprintf(stderr, "Cannot process the input, sorry\n");
			exit(1);

	}

	free(tree_strings);
	return MyTrees;

}

/* Structure to store the pvalues (for epics) or the rates (for epocs)	*/
/* And their index in the vector (sorted in trees order at the beginning	*/
/* More informative than just a vector because I can retrieve the indices at all time	*/
typedef struct {
	double value;
	int index;
} R;

/* The compare_2 function allows to sort on a struct													*/
/* i.e., I have a struct w/ values and indices, I retrieve the sorted values and associated indices!	*/

int compare_2 (const void * a, const void * b)
{
  	R *fa = (R *) a;
  	R *fb = (R *) b;
	if ( fa->value < fb->value)
        return -1;

     else if (fa->value > fb->value)
        return 1;

     else
        return 0;
}

/* Display the results of epics/epocs with the 95% confidence intervals	*/
void DisplayEpoicsWithConfidenceIntervals(struct Trees *MyForest, int ntrees, double *lkvectors, int mode_index, int evti, int evtj, int isItEpics){
	int i,j,z,flag;
	double mode, firstCI, secondCI, sum;

	R *values_vector=NULL;
	values_vector = malloc(sizeof(R)*ntrees);

	if(isItEpics){

		flag=0;
		sum=0.0L;
		mode=0.0L;
		firstCI=0.0L;
		secondCI=0.0L;

		for (z=0; z<ntrees; z++){
			values_vector[z].index=z;
			values_vector[z].value=MyForest[z].pvalue_of_pair;
		}
		// Sorting the values vector in ascending order
		qsort(values_vector, ntrees, sizeof(R), compare_2);

		// And now, I'm retrieving the mode of the distribution
		// I have a flag to mark whether the mode is within the 95%CI (flag=1) or not
		for (z=0; z<ntrees; z++){
			sum+=lkvectors[values_vector[z].index];
			if (values_vector[z].index == mode_index){
				if (sum<0.025 || sum>0.975){
					fprintf(stderr, "mode not within confidence intervals for evts %d vs %d!\n", evti+1, evtj+1);
					mode=values_vector[z].value;
					break;
				}
				flag=1;
				break;
			}
		}
		sum=0.0L;
		if (flag){
			for (z=0; z<ntrees; z++){
				sum+=lkvectors[values_vector[z].index];
				if ((sum>0.025) && (firstCI<=0)){
					firstCI=values_vector[z].value;
				}
				if ((sum>0.975) && (secondCI<=0)){
					secondCI=values_vector[z].value;
				}
				if (values_vector[z].index == mode_index){
					mode=values_vector[z].value;
				}
			}
			fprintf(stdout, "e%-4d vs e%-4d: pval=%.5g[%.5g, %.5g]\n", evti, evtj, mode, firstCI, secondCI);
		}
		else{
			fprintf(stdout, "e%-4d vs e%-4d: pval=%.5g[%c, %c]\n", evti, evtj, mode, '?', '?');
		}
	}
	else{
		fprintf(stdout,"%10c", ' ');
		fprintf(stdout,"%27s %27s %27s %27s\n","mu", "mu*", "nu", "nu*");
		for (i=0; i<2; i++){
			fprintf(stdout,"event[%d]:",(i==0)?evti+1:evtj+1);
			for (j=0; j<4; j++){
				flag=0;
				sum=0.0L;
				mode=0.0L;
				firstCI=0.0L;
				secondCI=0.0L;
				for (z=0; z<ntrees; z++){
					values_vector[z].index = z;
					values_vector[z].value = MyForest[z].final_rates_epocs[i][j];
				}
				qsort(values_vector, ntrees, sizeof(R), compare_2);
				for (z=0; z<ntrees; z++){
					sum+=lkvectors[values_vector[z].index];
					if (values_vector[z].index == mode_index){
						if (sum<0.025 || sum>0.975){
							fprintf(stderr, "mode not within confidence intervals for evts %d vs %d!\n", evti+1, evtj+1);
							mode=values_vector[z].value;
							break;
						}
						flag=1;
						break;
					}
				}
				sum=0.0L;
				if (flag){
					for (z=0; z<ntrees; z++){
						sum+=lkvectors[values_vector[z].index];
						if ((sum>0.025) && (firstCI<=0)){
							firstCI=values_vector[z].value;
						}
						if ((sum>0.975) && (secondCI<=0)){
							secondCI=values_vector[z].value;
						}
						if (values_vector[z].index == mode_index){
							mode=values_vector[z].value;
						}
					}
					if (mode >= 0){
						if (mode >= MAX_MU){
							fprintf(stdout, "%17.5g+[%.5g-%.5g]", mode>MAX_MU?MAX_MU:mode, firstCI, secondCI);
						}
						else{
							fprintf(stdout, "%18.5g[%.5g-%.5g]", mode, firstCI, secondCI);
						}
					}
					else{
						fprintf(stdout, "%28c",'?');
					}

				}
				else{
					fprintf(stdout, "%15.5g[%3c-%3c]", mode, '?', '?');
				}
			}
			fprintf(stdout, "\n");
		}
		fprintf(stdout, "\n");
	}

	free(values_vector);
	return;
}

/*
Here, I want to calculate the weights for the forest analysis
The idea is to either take the likelihoods of the tree and compute the distribution (epics)
Or to compute the product of lk_epocs and lk_tree (epocs)
INPUT: the struct Trees, which contains the likelihood values
*/

double * ComputeLikelihoodDistribution(struct Trees *MyForest, int ntrees, int isItEpics, int evti, int evtj){
	int f,maxf;
	double *res=NULL;
	double sum=0.0L;
	double maxlk=0.0L;

	R *rates_vector=NULL;
	rates_vector = malloc(sizeof(R)*ntrees);

	res = calloc(sizeof(double), ntrees);

	maxf=0;

	// First: retrieving the likelihoods, normalize them and store in array
	// Two different ways possible: epics (only the lk of the tree)
	// epocs: product of epocs and tree likelihoods (here: sum of the logs!)

	if (isItEpics){
		for (f=0; f<ntrees; f++){
			sum+=MyForest[f].lik_ntrees;
			if (MyForest[f].lik_ntrees>maxlk){
				maxlk=MyForest[f].lik_ntrees;
				maxf=f;	// The mode of the distribution at the end
			}
		}
		for (f=0; f<ntrees; f++){
			res[f]=MyForest[f].lik_ntrees/sum;
		}
		DisplayEpoicsWithConfidenceIntervals(MyForest, ntrees, res, maxf, evti, evtj, isItEpics);
	}
	else{
		for (f=0; f<ntrees; f++){
			sum+=MyForest[f].lik_ntrees + MyForest[f].lik_epocs;
			if ((MyForest[f].lik_ntrees+MyForest[f].lik_epocs)>maxlk){
				maxlk=MyForest[f].lik_ntrees+MyForest[f].lik_epocs;
				maxf=f;													// The mode of the distribution at the end
			}
		}
		for (f=0; f<ntrees; f++){
			res[f]=(MyForest[f].lik_ntrees + MyForest[f].lik_epocs)/sum;
		}
		DisplayEpoicsWithConfidenceIntervals(MyForest, ntrees, res, maxf, evti, evtj, isItEpics);
	}

	// Now, I retrieve the rates one by one and sort them in ascending order.
	// Then, I need to retrieve the mode (i.e. the max lk value) and check if it is within the 95% CI
	// If yes, I can proceed -> I display the rate corresponding of the mode and the associated CI [2.5% - 97.5%]

	return res;
}

/* I/O options. 										   */
/* To initialize a new file [prefix]/[prefix]_epocs_[scenario].tab that */
/* will contain the results of epocs.					   */
/* To parse the input later on, the scenario is in the header */

char * initializeNewFile(char * pref, char scen, char *matID, int whichOutput){

    static char MyOutputName[1000]="";
	char *basename;
	char cond=1;

    struct stat st = {0};

	if (strrchr(pref, '/')){
		cond=0;
		basename=strrchr(pref, '/')+1;
	}

	/* I'm checking whether the output folder exists 			 */
	/* if no, I'm creating it so that I can create a file inside */
    if (stat(pref, &st) == -1) {
        mkdir(pref, 0700);
    }

	// switch on whichOutput
	// 0 = epocs output file
	// 1 = epics output file
	// 2 = epics output file - only significant pairs
	// 3 = epics output file - multi tree results
	// 4 = epics output file - multi tree results - significant pairs
	// 5 = epics output file - summary of forest output
	// 6 = epocs output file - summary of forest output

	switch(whichOutput){
		case 0:
			if (scen >= 'a' && scen <= 'z')
		    	sprintf(MyOutputName, "%s/%s_mat_epocs_%c%c.tab", pref, cond?pref:basename, scen, scen);
			else
				sprintf(MyOutputName, "%s/%s_mat_epocs_%c.tab", pref, cond?pref:basename, scen);
			break;

		case 1:
			sprintf(MyOutputName, "%s/%s_mat_epics_%s.tab", pref, cond?pref:basename, matID);
			break;

		case 2:
			sprintf(MyOutputName, "%s/%s_mat_signif_epics_%s.tab", pref, cond?pref:basename, matID);
			break;

		case 3:
			sprintf(MyOutputName, "%s/%s_mat_epics_%s_multitree.tab", pref, cond?pref:basename, matID);
			break;

		case 4:
			sprintf(MyOutputName, "%s/%s_mat_signif_epics_%s_multitree.tab", pref, cond?pref:basename, matID);
			break;

		case 5:
			sprintf(MyOutputName, "%s/%s_mat_epics_%s_forest.tab", pref, cond?pref:basename, matID);
			break;

		case 6:
			if (scen >= 'a' && scen <= 'z')
				sprintf(MyOutputName, "%s/%s_mat_epocs_%c%c_forest.tab", pref, cond?pref:basename, scen, scen);
			else
				sprintf(MyOutputName, "%s/%s_mat_epocs_%c_forest.tab", pref, cond?pref:basename, scen);
			break;

		default:
			fprintf(stderr, "initializeNewFile -> option:<%d> unknown... exiting", whichOutput);
			exit(1);
	}
    return MyOutputName;

}

/* Create an array (pointer of pointers) that will contain the 	   */
/* event/variant ids. If option is chosen, these ids will replace  */
/* the automatic notation of e0, e1, etc.                          */
/* the input file is a 1-column file, with the event names sorted  */
/* as in the epics tree. Maximum string length is 128.             */
char ** createEventsArray(char * myFname){

    char **outArray;
    FILE *fPointer;
    char line[128];
    int count=0;

    fPointer = fopen(myFname, "r");
    while (fgets(line,128,fPointer))
        count++;
    fclose(fPointer);

    outArray = (char**) malloc(count * sizeof(char*));
	for (int i=0; i<count; i++){
		outArray[i] = (char*)malloc(sizeof(char)*128);
	}
    if (!outArray)fprintf(stderr, "createEventsArray: cannot allocate outArray, exiting\n "),exit(3);
    count=0;

    fPointer = fopen(myFname, "r");
    while (fgets(line,128,fPointer)){
        line[strcspn(line, "\n")]=0;
		strcpy(outArray[count], line);
		count++;
	}

	fclose(fPointer);

    return outArray;

}

/* Creates a struct of event pairs to run in epics only				*/
/* The input is a string of a file with double columns e1 e2		*/
/* a pointer to an int for the length of the file and a verbose int	*/
/* The output is an 'EventsArray' with a vector of e1 to run		*/
/* the length of the array, the length of e2 for each e1 and an		*/
/* asymmetrical e1xe2 array for each event to run against each e1	*/
/* I also have a vector e2L containing the lengths of each e2		*/

EventsArray *createArrayEventsToRunEpics(char *MyInputName, int verbose){

	EventsArray *outArray;
	outArray = (EventsArray*) malloc(sizeof(EventsArray));
	outArray->length=0;
	outArray->e1=NULL;
	outArray->e2L=NULL;
	outArray->inpe2=NULL;

	fpos_t position;
	int tmp1,tmp2, n1, n2, cnt1, cnt2;
	FILE *f;
    char line[128];
	char *header=NULL;
	int e1[10000];

	f=fopen(MyInputName, "r");

	while (fgets(line, 128, f)){
		if (!(header=strstr(line, "IDe1"))){
			sscanf(line, "%d%d", &tmp1, &tmp2);
			n1=tmp1;
			e1[outArray->length]=0;
			while (n1==tmp1){
				e1[outArray->length]+=1;
				fgetpos(f, &position);
				if(fgets(line, 128, f) == NULL){
					break;
				}
				sscanf(line, "%d%d", &n1, &n2);
			}
			fsetpos(f, &position);
			outArray->length++;
		}
	}

	rewind(f);

	outArray->e1=calloc(outArray->length, sizeof(int));
	outArray->e2L=calloc(outArray->length, sizeof(int));
	outArray->inpe2=(int**) malloc(sizeof(int **)*outArray->length);
	for (int i=0; i<outArray->length;i++){
		outArray->inpe2[i] = (int*)malloc(e1[i]*sizeof(int));
	}

	cnt1=0;
	while (fgets(line, 128, f)){
		cnt2=0;
		if (!(header=strstr(line, "IDe1"))){
			sscanf(line, "%d%d", &tmp1, &n2);
			n1=tmp1;
			outArray->e1[cnt1]=n1;
			while (n1==tmp1){
				outArray->inpe2[cnt1][cnt2] = n2;
				outArray->e2L[cnt1] += 1;
				cnt2++;
				fgetpos(f, &position);
				if(fgets(line, 128, f) == NULL){
					break;
				}
				sscanf(line, "%d%d", &n1, &n2);
			}
			fsetpos(f, &position);
		}

		cnt1++;
	}

	fclose(f);

	return outArray;

}

/*	When reading an input file, create the array of events to run	*/
int ** createArrayEventsToRun(char *MyInputName, int *fileLength, int verbose)
{
	int **outArray;
    FILE *fPointer;
    char line[128];
	int count=0;
	char *header=NULL;

    fPointer = fopen(MyInputName, "r");
    while (fgets(line,128,fPointer))
        *fileLength+=1;
    // fclose(fPointer);
	rewind(fPointer);
	count=*fileLength;

    outArray = (int**) malloc(*fileLength * sizeof(int*));
	for (int i=0; i<*fileLength; i++){
		outArray[i] = (int*)malloc(2 * sizeof(int));
	}
    if (!outArray)fprintf(stderr, "createArrayEventsToRun: cannot allocate outArray, exiting\n "),exit(3);
    *fileLength=0;

    // fPointer = fopen(MyInputName, "r");
    while (fgets(line,128,fPointer)){

		if (!(header=strstr(line, "IDe1"))) /* Checking that I'm not reading the header. */
		{
			sscanf(line, "%d%d", &outArray[*fileLength][0], &outArray[*fileLength][1]);
			if (outArray[*fileLength][1] < outArray[*fileLength][0])
			{
				int temp=outArray[*fileLength][0];
				outArray[*fileLength][0]=outArray[*fileLength][1];
				outArray[*fileLength][1]=temp;
			}
			if (isDuplicate(outArray[*fileLength][0], outArray[*fileLength][1], outArray, *fileLength))continue;

			*fileLength+=1;
		}
	}

	fclose(fPointer);

	for (int i=*fileLength; i<count; i++)
		free(outArray[i]);

	if(verbose)fprintf(stderr, "there are %d pair of events to run\n", *fileLength);

	return outArray;
}

/*	Function important when inputting a file: need to not have duplicate events	*/
/*	e.g. 1 vs 2 and 2 vs 1. This can happen with epics output.	*/
int isDuplicate(int number1, int number2, int **storedArray, int lastPosition){

	for (int i=0; i<lastPosition; i++){
		if ((storedArray[i][0] == number1) &&  (storedArray[i][1] == number2))return 1;
	}
	return 0;
}

/* Write standard metrics from epics in output file. */
void writeInFile(int numOfTree, int forestLength, char *eventIds, FILE *filePointer, char **inputIds, char *matId, int e1, int e2, int treeId, int treeInForest, double pvalue, int obsPairs, int maxObs){

	if (numOfTree > 1){
		if (eventIds){
			fprintf(filePointer, "%-10d\t%-10d\t%-10s\t%-10s\t%-5s\t%-5d\t%-8.3e\t%-4d\t%-4d\n", e1+1, e2+1, inputIds[e1], inputIds[e2], matId, treeId+1, pvalue, obsPairs, maxObs);
		}else{
			fprintf(filePointer, "%-10d\t%-10d\t%-10s\t%-10s\t%-5s\t%-5d\t%-8.3e\t%-4d\t%-4d\n", e1+1,  e2+1, "NA", "NA", matId, treeId+1, pvalue, obsPairs, maxObs);
		}
	}
	else{
		if (forestLength == 1){
			if (eventIds){
				fprintf(filePointer, "%-10d\t%-10d\t%-10s\t%-10s\t%-5s\t%-8.3e\t%-4d\t%-4d\n", e1+1, e2+1, inputIds[e1], inputIds[e2], matId, pvalue, obsPairs, maxObs);
			}else{
				fprintf(filePointer, "%-10d\t%-10d\t%-10s\t%-10s\t%-5s\t%-8.3e\t%-4d\t%-4d\n", e1+1,  e2+1, "NA", "NA", matId, pvalue, obsPairs, maxObs);
			}
		}
		else{
			if (eventIds){
				fprintf(filePointer, "%-10d\t%-10d\t%-10s\t%-10s\t%-5s\t%-4d\t%-8.3e\t%-4d\t%-4d\n", e1+1, e2+1, inputIds[e1], inputIds[e2], matId, treeInForest, pvalue, obsPairs, maxObs);
			}else{
				fprintf(filePointer, "%-10d\t%-10d\t%-10s\t%-10s\t%-5s\t%-4d\t%-8.3e\t%-4d\t%-4d\n", e1+1,  e2+1, "NA", "NA", matId, treeInForest, pvalue, obsPairs, maxObs);
			}
		}
	}

	return;
}

/*	Free a 2D array	*/
void freeArray(int **a, int m) {
    int i;
    for (i = 0; i < m; ++i) {
        free(a[i]);
    }
    free(a);
}

/*	Freeing the CoevolData object	*/
void FreeCoevolData( struct CoevolData *t, int n, int nevts, int mat_type ){

	int i,j;

	for( i=0; i<n ; i++)
	{

		if( t[i].tree ) free( t[i].tree );
		if( t[i].branch_lengths ) free( t[i].branch_lengths );
		if( t[i].v_evts ) free( t[i].v_evts );
		if( t[i].chronomat ) free( t[i].chronomat );
		if( t[i].sparseChronomat ) free( t[i].sparseChronomat );

		if (t[i].e1){

			for (j = 0; j < nevts; j++)
				if (t[i].e1[j]) free(t[i].e1[j]);

			free(t[i].e1);
		}

		if( mat_type && t[i].e1M ) free( t[i].e1M ); /// CHECK CHECK CHECK

	}

	if( t )free( t );

}
