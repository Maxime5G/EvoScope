/*
	Copyright (C) 2002-2018 K Bedhenna G Achaz & S Brouillet Joel Pothier William Belaid

	This program is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public License
	as published by the Free Software Foundation; either version 2.1
	of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 	for more information, please contact guillaume achaz <guillaume.achaz@mnhn.fr>/<gachaz@gmail.com>
*/
/**
\file 	epics.c
\author K Bedhenna, J Pothier, G Achaz, S Brouillet,  W Bellaid, M Godfroid
\date 	Sept 2015 - May 2021
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"
#include "coevol.h"

#include <sys/stat.h> /* for the creation of directories and files */
#include <sys/types.h> /* for the creation of directories and files */

#include <unistd.h>     /* pour getopt  */

#include <time.h>	/* for terminating the multinomial if it takes too long */

/* This function initializes empty vectors for the calculations of probabilities.  				 */
/* INPUTS: the number of classes (i.e. the possible results when throwing e2 on the tree), 	 	 */
/* the (static) number max of different classes, the (static) nveck vector of the categories, 	 */
/* the length of the probability distribution vector (equal to (n1*n2)+1), the maximum length 	 */
/* of the distrib and the vector of the probabilities.	 										 */
/* This function serves only the purpose of initializing all vectors with zero's everywhere.  	 */
/* For the veck vector, we first find the number of classes and initialize it with this length 	 */
/* For the probabilities vector, we find the maximum length of this vector and initializes it 	 */
/* following the results, with double 0.0L.														 */
void Init_veck_distrib( int nclasses, int *nclassmax, int **nveck, int ldistrib,  int *max_ldistrib, double **distrib )
{

	int k;

	/* allocate if needed vector of work veck */
	if (nclasses > *nclassmax) {
		*nclassmax = nclasses;
		if (*nveck)
			free(*nveck);
		if (( (*nveck) = (int *) calloc((size_t)nclasses,sizeof(int))) == NULL)
			fprintf(stderr, "Init_veck_distrib: Not enough memory for vecteur veck\n"),  exit(3);
		}

	/* allocate the  distribution array only if preceeding array was smaller	*/
	/* else array will be greater than needed, but we do not mind				*/
	if (ldistrib > *max_ldistrib)
	{
		if(*distrib) free(*distrib);

		if ( ( (*distrib) = (double *)malloc((size_t) sizeof(double)*ldistrib ) ) == NULL )
			fprintf(stderr, "Init_veck_distrib: Not enough memory for vecteur distrib\n"),exit(3);

		*max_ldistrib = ldistrib;
	}

	/* init the probability distrib of each pairs */
	for (k = 0; k < ldistrib; k++)
		(*distrib)[k] = 0.0L;

}

/* This function initializes the the first factorials of a given input. 					 */
/* INPUTS: the number of events, the number of trees and a pointer to the struc MyEpicsData. */
/* This function will go through each tree and each event and calculate the cumulative  	 */
/* sum of all mutations to retrieve the maximum number of factorials to compute later on. 	 */
/* Next - check comment below.																 */
void InitFacto( int nevt, int ntree, struct CoevolData *MyEpicsData, int verbose){

	/* To calculate the first n factorials. */
	int nmax = 0;
	int i,n1,t;


	for(t=0;t<ntree;t++){

		for (i = 0; i < nevt; i++)
		{
			n1 = sumvP(MyEpicsData[t].e1[i], MyEpicsData[t].nbranches);
			if (n1 > nmax)
				nmax = n1;
		}
	}

	if (verbose) fprintf(stderr, "Compute the %d first ln factorials\n",nmax);
	  init_slnFacto( nmax );

	return;

}

/* Verbose function to print the probability distribution vector, the number of hits and the probabilities associated.	*/
/* INPUTS: the probability distribution vector, the length of the vector and the type of the chronology matrix.			*/
void PrintDist(double *distrib, int ldistrib, int mat_type ){

	int k=0;
	int last_distvalue=0;
	double pval=0;

	k = ldistrib-1;
	while (k >= 1 && distrib[k] == 0L) {
		k--;
	}
	last_distvalue = k;
	pval = 0.0L;
	for (k = 0; k <= last_distvalue; k++) {

		fprintf(stdout, "Probability of %d %s = %g\n",k, (mat_type==0)?"co-occurences":(mat_type==1)?"chronologies":
                        "co-occurences or chronologies",distrib[k]);

			pval += distrib[k];
	}

	fprintf(stdout, "Sum of probabilites of distribution = %g\n",pval);

}

/* Function to calculate the index of the maximum value in a vector of probabilities	*/
int MaxDistrib( double *distrib, int ldistrib  ){

	int k;

	k = ldistrib-1;
	while (k >= 1 && distrib[k] == 0L) {
		k--;
	}
	return k;
}

/* Function to combine probability distributions when calculating the proba for multiple trees	*/
double *CombineDistrib( double *cum, int *size_cum,  double *new, int size_new  ){

	int i,j;
	int last_distvalue=0;

	double *Comb;
	double Max;


	last_distvalue =  MaxDistrib( new, size_new  );

	Max = (*size_cum==0)?last_distvalue+1:last_distvalue+*size_cum;

	Comb = (double *)calloc( Max, sizeof(double)  );
	if ( ! Comb )fprintf(stderr, "CombineDistrib: cannot allocate Comb, exiting\n "),exit(3);

	if( *size_cum == 0 )
	{
		for(i=0;i<=last_distvalue;i++)
			Comb[i] += new[i];
	}
	else
	{
		for(i=0 ; i<=last_distvalue ; i++)
			for(j=0 ; j<*size_cum ; j++){
				Comb[i+j] += new[i]*cum[j];
//				printf("Comb[%d] [%d] + [%d] = %lf * %lf = %lf\n",i+j,i,j,new[i],cum[j],Comb[i+j] );
			}

		free(cum);
	}

	*size_cum=Max;
	return Comb;

}

/* Function to compute the pvalue from a probability distribution. i.e. sum the proba until you	*/
/* reach the observed value	*/
double ComputePval( double *distrib, int ldistrib, int obs){

	int i;
	double pval=0.0L;

	for(i=obs;i<ldistrib; i++){
		pval += distrib[i];
	}

	return pval;
}

/* -----------------------------------------------------------------	*/
void PrintHelp(char *s)
{
	fprintf(stderr,"Usage: %s [options] outgroup [tree_file]\n", s);
	fprintf(stderr,"or:    %s [options] -f all_trees_file\n", s);

	fprintf(stderr," [first usage: a single tree]\n");
	fprintf(stderr,"     outgroup: specify the outgroup using leaves\n");
	fprintf(stderr,"               format is \"name1[:name2:name3:...]\" \n" );
	fprintf(stderr, "              provide an empty string (\"\") if you have no outgroup\n");
	fprintf(stderr, "              When counting chronologies (S or B matrices), make sure your tree is correctly rooted !!\n");
	fprintf(stderr,"     tree_file: tree file needs to be provided in a modified Newick format\n");
	fprintf(stderr,"                with events provided at the locations usually dedicated to bootstrap values (see README).\n");
	fprintf(stderr,"                if not provided, read the tree from stdin.\n");

	fprintf(stderr," [second usage: many trees with the same events]\n");
	fprintf(stderr,"     all_trees_file: a file with all info on each tree each line as 2 columns\n");
	fprintf(stderr,"                     1: a tree_file (see above) \n");
	fprintf(stderr,"                     2: an outgroup (or \"\") for no outgroup) \n");
	fprintf(stderr,"     all trees must have all their events in the same order !\n");


	fprintf(stderr,"\n [matrix choice]\n");
	fprintf(stderr,"    -I: inseparable pairs only will be computed (Identity matrix)\n");
	fprintf(stderr,"    -S: genealogically ordered pairs only will be computed (S matrix) [default]\n");
	fprintf(stderr,"    -B: (both) inseparable and genealogically ordered pairs will be computed (S+Id matrix)\n");
	fprintf(stderr,"    -A: (both) inseparable and adjacent pairs will be computed (A+Id matrix)\n");
    fprintf(stderr,"    -P: (both) inseparable and precedent pairs will be computed (P+Id matrix)\n");
    fprintf(stderr,"    -X: Exclusion matrix will be computed (X matrix)\n");

	fprintf(stderr,"\n [filters on output pairs]\n");
	fprintf(stderr,"    -s <real value>: set p-value threshold to output on stdout a Newick tree with events in name (default: 0.05)\n");
	fprintf(stderr,"    -M: output also monomorphic event sites (i.e. couple of events that do not occur) [default not]\n");
	fprintf(stderr,"    -0: output also p-value of 0 inseparable/genealogically ordered pairs (p-value=1)\n");
	fprintf(stderr,"    -T <real value>: set p-value threshold to output for pairs (default: 1.0, i.e. show all)\n");

	fprintf(stderr,"\n [output tree with events]\n");
	fprintf(stderr,"   -N: output on stdout a Newick tree with events (only pairs with pval < 0.05) [default not]\n");

	fprintf(stderr, "\n [input and output operations]\n");
    fprintf(stderr, "   -E <events>: reads a 1-column file for displaying event IDs instead of default [default: e1,e2...]\n");
	fprintf(stderr, "   -R <filename>: input filename containing significant events to run \n");
    fprintf(stderr, "   -O <output>: set a value for the output folder and files prefixes [default: out]\n");
	fprintf(stderr, "   -1 <event 1>: set the first event id (integer) to be tested \n");
	fprintf(stderr, "   -2 <event 2>: set the second event id (integer) to be tested \n");

	fprintf(stderr,"\n [verbose options]\n");
	fprintf(stderr,"    -v: (verbose) outputs running progression of the program\n");
	fprintf(stderr,"    -V: (very verbose) outputs running progression of the program and more...\n");
	fprintf(stderr,"    -d: (distribution) output on stderr the inseparable and/or genealogically ordered pairs probabilities distribution \n");
}


/* -----------------------------------------------------------------	*/
int main( int argc , char **argv)
{
   int nevt=-1,
       i, j, k,
       verbose = 0,
       carg,
	   *veck=NULL,
	   *nveck=NULL,
	   n1,
	   n2,
	   ldistrib,
	   mat_type,       // 0: Id, 1: S, 2: Both, 3: A+Id, 4: P+Id, 5: X
	   n_obs_pairs,
	   max_ldistrib,
	   *ne1M = NULL,
	   *ne2 = NULL,
	   nclasses=0,
	   nclassmax,
	   nn2,
	   output_dist = 0,
	   output_monomorph = 0,
	   output_00 = 0,
	   t=0,
	   ntree=1,
	   output_newick=0,
       choice=0,					/* which event you want to compare 		 */
	   do_e1M,						/* to decide if we need to calculate e1M */
	   choice_e1=0,			/* if the user wants to input a particular event to test - e1 */
	   choice_e2=0,			/* if the user wants to input a particular event to test - e2 */
	   mymax=0;				/* max possible observed co-occurrences in the tree	*/


	double *distrib = NULL,			/* proba distribution	*/
	       *nbranch_lengths=NULL,	/* proba vector for each category in the multinomial	*/
           pval,					/* combined pvalue for multiple trees 	*/
           seuil_evt = 0.05L,		/* pvalue threshold to determine whether to output a pair on the newick tree	*/
   		   seuil_pval = 0.05L;		/* pvalue threshold to consider a pair significant	*/

    double *lkDist=NULL;			/* likelihood distribution in case of forest implementation	*/

	char * infile=NULL;
	char *input_events_file=NULL;		/* pointer to open input event to run */
	EventsArray *input_events=NULL;

	char output_prefix[1000]="out"; /* default output prefix           */

	char *fname=NULL; 				/* output file name 			   */
	char *fnames=NULL;				/* output file name for signif pairs */
    FILE *fptr; 					/* pointer to open the output file */
	FILE *fptrs;					/* pointer to open output file for signif pairs */

	char *fname_multi=NULL;			/* same name and file pointers for multiple trees	*/
	FILE *fptr_m;
	char *fname_multis=NULL;
	FILE *fptr_ms;

	char *fnamef=NULL;				/* same name and file pointers for forest of trees	*/
	FILE *fptrf;

    char * input_names=NULL;		/* input file for the event ids    */
    char **input_ids;				/* array of the event ids          */
	char *mat_name=NULL;			/* matrix name for the output file */

	struct Trees *MyTrees=NULL;		/* initializing the forest	*/
	int *forestntree=NULL;
	int f;							/* the forest counter */

	/* default parameters	*/
	mat_type = 1;	/* matrice is chronology	*/
	mat_name = "S";
	verbose = 0;

	/* -- recursive gathering of command line arguments -- */
	while ((carg = getopt(argc, argv, "1:2:0IBAPXMNSVT:E:C:O:R:abdehms:tvf:")) != -1) {

		switch ((char)carg) {

			case '0':
				output_00 = 1; /* output 0-0 chronologies */
				break;

			case 'V':
				verbose++;
			case 'v':
				verbose++;
				break;

			case 'I':
				mat_type = 0;	/* matrice is identity	*/
				mat_name = "I";
				break;

		   case 'S':
				mat_type = 1;	/* matrice is chronology	*/
				mat_name = "S";
				break;

			case 'B':
				mat_type = 2;	/* matrice is Id+chronology	*/
				mat_name = "B";
				break;

			case 'A':
                mat_type = 3;	/* matrice is Id+Adj	*/
				mat_name = "A+Id";
                break;

            case 'P':
                mat_type = 4;	/* matrice is Precedence+Id	*/
				mat_name = "P+Id";
                break;

            case 'X':
                mat_type = 5;	/* matrice is Exclusion	*/
				mat_name = "X";
                break;

			case 'M':
				output_monomorph = 1;	/* output monomorphic sites */
				break;

			case '1':
				if (sscanf(optarg, "%d", &(choice_e1)) != 1){
					fprintf(stderr, "cannot read the choice id. Please enter an integer, exiting... \n"), exit(1);
                    fprintf(stderr, "    setting it to default (all events e1)\n");
				}
				break;

			case '2':
				if (sscanf(optarg, "%d", &(choice_e2)) != 1){
					fprintf(stderr, "cannot read the choice id. Please enter an integer, exiting... \n"), exit(1);
					fprintf(stderr, "    setting it to default (all events e2)\n");
				}
				break;

			case 'R':
                input_events_file=(char *)malloc((strlen(optarg)+1)*sizeof(char));
				if (!(strcpy(input_events_file, optarg))){
                    fprintf(stderr, "Cannot read input file or too long name (reading %s)\n", optarg);
                    fprintf(stderr, "    setting it to none\n");
                }
                break;

			case 's':	/* threshold of p-value to output event pairs in the output tree	*/
				if (sscanf(optarg, "%lf", &(seuil_evt)) != 1) {
					fprintf(stderr,"%s: cannot read threshold of p-value to output event pairs in the output tree (option -s: reading %s (?))\n",argv[0],optarg);
					fprintf(stderr,"    setting it to 0.05\n");
					seuil_evt = 0.05L;
				}
				if (seuil_evt < 0.0L || seuil_evt > 1.0L) {
					fprintf(stderr,"%s: given threshold of p-value to output event pairs is < 0.0 or > 1.0\n",argv[0]);
					fprintf(stderr,"    setting it to 0.05\n");
					seuil_evt = 0.05L;
				}
				/* continue with setting output newick tree as if N option was selected */

			case 'N':
				output_newick = 1;	/* outputs a newick tree with events which pval > 0.05,  readable by NJPlot e.g.	*/
				break;

			case 'T':	/* threshold of p-value to discriminate significant event pairs 	*/
				if (sscanf(optarg, "%lf", &(seuil_pval)) != 1) {
					fprintf(stderr,"Cannot read threshold of p-value to determine significant pairs (reading %s (?))\n",optarg);
					fprintf(stderr,"    setting it to 0.05 (output all significant pairs with p-value <= 0.05)\n");
					seuil_pval = 0.05L;
				}
				break;

			case 'E':
                input_names=(char *)malloc((strlen(optarg)+1)*sizeof(char));
				if (!(strcpy(input_names, optarg))){
                    fprintf(stderr, "Cannot read output prefix or too long name (reading %s)\n", optarg);
                    fprintf(stderr, "    setting it to none\n");
                }
                break;

            case 'O':
                if (!(strcpy(output_prefix, optarg))){
                    fprintf(stderr, "Cannot read output prefix or too long name (reading %s)\n", optarg);
                    fprintf(stderr, "    setting it to out (default)\n");
                }
                break;

			case 'd':
				output_dist = 1;	/* outputs the distribution	*/
				break;

			case 'h':
				PrintHelp(argv[0]);
				exit(0);

			case 'f':
				infile=optarg;
				break;

			default:
				fprintf(stderr,"option:<%c> unknown... exiting\n",(char)carg);
				exit(1);
		}

	}

	if(argc == 1)
		fprintf(stderr, "Usage: ./epics [options] outgroup [tree_file]\nor:    ./epics [options] -f all_trees_file\nfor more help ./epics -h\n"), exit(2);

	/*
		Retrieving the tree informations. The ParseInputForest reads the tree or the treefile and determines whether we have
		a forest, multiple independent trees or a single newick tree. All informations are stored in a Trees structure
	*/
	MyTrees = ParseInputForest(argc, optind, argv, infile, verbose, &nevt, &ntree, &forestntree, mat_type, 1);

	int maxforest=0;
	int forestlength=0;

	for (t=0; t<ntree; t++){
		forestlength=forestntree[t];
		if (forestlength>maxforest) maxforest=forestlength;
	}
	/* ---------------------------------------------------------------------------------------------- */
    /* Check that, if the user chose events to run, e1 < e2 (given the loop) 	  	 		  		  */
    /* This also verifies that if the user mistakenly input -2 only, choice_e2 is copied to choice_e1 */

	if ((input_events_file) && (choice_e1 || choice_e2))fprintf(stderr,"ERROR: Choices not permitted when inputting an events file\n"), exit(1);

    if (!choice_e1 && choice_e2){
 	   choice=choice_e2;
 	   choice_e2=choice_e1=0;
    }

	if (!choice_e2 && choice_e1){
 	   choice=choice_e1;
 	   choice_e2=choice_e1=0;
    }

    if ((choice_e1 && choice_e2) && (choice_e1 == choice_e2))fprintf(stderr,"ERROR: Please input different numbers for the event choice\n"), exit(1);

    /* ---------------------------------------------------------------------------------------------- */

	/* Initializing the factorial distribution - stored in slnFacto	*/
	InitFacto( nevt, ntree, MyTrees[0].MyCoevolData, verbose);
	/* print the tree on stdout (without outgroup) */
	if ( verbose ) {
		for (f=0; f<maxforest; f++){
			for(t=0;t<ntree;t++)
			{
				printf("f=%d, t=%d\n", f, t);
				fprintf(stderr,"\nAscii Tree (%d):\n",t);
				print_arbre(MyTrees[f].MyCoevolData[t].root, stderr);
			}
		}
	}

	/* --- FILE INITIALIZATIONS ---	*/

	fname=initializeNewFile(output_prefix, 'i', mat_name, 1);				/* create the string for the output file and */
    fptr=fopen(fname, "w");													/* opening it */

	fnames=initializeNewFile(output_prefix, 'i', mat_name, 2);
	fptrs=fopen(fnames, "w");

	/* the output format is e0, ev0, e1, ev1, matrix, pval, nobs, maxobs                 */
	/* if multiple trees in input: e0, ev0, e1, ev1, tree id, pval, nobs, maxobs */

    if (ntree>1){
		fname_multi = initializeNewFile(output_prefix, 'i', mat_name, 3);
		fptr_m = fopen(fname_multi, "w");
		fname_multis = initializeNewFile(output_prefix, 'i', mat_name, 4);
		fptr_ms = fopen(fname_multis, "w");
		fprintf(fptr, "%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-5s\t%-8s\t%-4s\t%-4s\n", "IDe1", "IDe2", "Event1", "Event2", "Matrix", "Tree", "Pval", "nobs", "max");
		fprintf(fptrs, "%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-5s\t%-8s\t%-4s\t%-4s\n", "IDe1", "IDe2", "Event1", "Event2", "Matrix", "Tree", "Pval", "nobs", "max");
		fprintf(fptr_m, "%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-8s\t%-4s\t%-4s\n", "IDe1", "IDe2", "Event1", "Event2", "Matrix", "Pval", "nobs", "max");
		fprintf(fptr_ms, "%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-8s\t%-4s\t%-4s\n", "IDe1", "IDe2", "Event1", "Event2", "Matrix", "Pval", "nobs", "max");
	}
    if (ntree == 1){
		if (maxforest>1){
			fnamef = initializeNewFile(output_prefix, 'i', mat_name, 5);
			fptrf = fopen(fnamef, "w");
			fprintf(fptrf, "%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-5s\t%-8s\t%-4s\t%-4s\n", "IDe1", "IDe2", "Event1", "Event2", "Matrix", "TOTHINK", "Pval", "Nobs", "Max");
			fprintf(fptr, "%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-5s\t%-8s\t%-4s\t%-4s\n", "IDe1", "IDe2", "Event1", "Event2", "Matrix", "Tree", "Pval", "Nobs", "Max");
			fprintf(fptrs, "%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-5s\t%-8s\t%-4s\t%-4s\n", "IDe1", "IDe2", "Event1", "Event2", "Matrix", "Tree", "Pval", "Nobs", "Max");
		}else{
			fprintf(fptr, "%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-8s\t%-4s\t%-4s\n", "IDe1", "IDe2", "Event1", "Event2", "Matrix", "Pval", "Nobs", "Max");
			fprintf(fptrs, "%-10s\t%-10s\t%-10s\t%-10s\t%-5s\t%-8s\t%-4s\t%-4s\n", "IDe1", "IDe2", "Event1", "Event2", "Matrix", "Pval", "Nobs", "Max");
		}
	}

    if (input_names){
        input_ids = createEventsArray(input_names); 		/* if -E chosen, creates the array of event ids. */
    }

	/* Verifying if the user gave an input file of significant events to run. If yes, store the events */
	/* in an EventsArray struct. */

	if (input_events_file){
		input_events = createArrayEventsToRunEpics(input_events_file, verbose); 		/* if -R chosen, creates the array of event ids to run. */
	}

	/* for each pair of the nevt*(nevt-1)/2 event vectors */
	ldistrib = 0;
	max_ldistrib = 0;
	nclassmax = 0;

	double weight_cum = 0;
	double *distrib_cum=NULL;
	int max_cum=0;
	int total_obs=0;
	int ntree_skipped=0;

	int cnti=0,
	 	cntj=0;

	/* --- MAIN LOOP --- */

	for (i = 0; i < nevt; i++) {

		if (input_events_file){
			if ((i+1) != input_events->e1[cnti]) continue;
		}

		/* condition to check whether to do e1M in the loop. Modified after the t loop */
		/* such that it is performed only once per tree and per e1 event.              */
		do_e1M = 0;

		if ((choice_e1 && choice_e2) && (choice_e1!=i+1 && choice_e2!=i+1)) continue;

  		for (j = 0; j < nevt; j++) {

			if (input_events_file){
				if ((j+1) != input_events->inpe2[cnti][cntj]) continue;
			}

  			if(i==j) continue;

			/* if one event is chosen, continue the loop until condition is met */
			if ((choice>0) && (choice_e1==0 && choice_e2==0) && ((i != choice-1) && (j != choice-1))) continue;

			if ((choice_e1 && choice_e2) && (choice_e1!=j+1 && choice_e2!=j+1)) continue;

			printf("\n== e%d vs e%d ==\n",i+1,j+1);

			for (f=0; f<maxforest; f++){

				if(distrib_cum){free(distrib_cum); distrib_cum=NULL;}
				max_cum=0;
				weight_cum=0;
				total_obs=0;
				ntree_skipped=0;

				for(t=0;t<ntree;t++)
				{

					/*
						Output monomorph is off and one of the site is monomorph - skip the tree
					*/

					if (output_monomorph==0 && (MyTrees[f].MyCoevolData[t].v_evts[i]==0 || MyTrees[f].MyCoevolData[t].v_evts[j]==0) ){
						ntree_skipped++;
						continue;
					}

					if(verbose>1) fprintf(stderr, "events %d vs %d, tree %d\n", i+1,j+1,maxforest>1?f:t+1);

					/*
						n1 and n2 - cumulative sums of e1 and e2 respectively
					*/

					n1 = sumvP( MyTrees[f].MyCoevolData[t].e1[i],  MyTrees[f].MyCoevolData[t].nbranches); 	/* cumulative sum of e1	*/
					n2 = sumvP( MyTrees[f].MyCoevolData[t].e1[j],  MyTrees[f].MyCoevolData[t].nbranches); 	/* cumulative sum of e2	*/
					if(verbose>1) fprintf(stderr, "n1 %d ; n2 %d\n", n1,n2);

					/*
						Compute e1*M - matrix product between e1 vector and the scenario matrix
					*/

					if (!do_e1M)	// We don't need to compute e1.M at every round
					{
						if (mat_type != 0)	vectSparseMatAndMask( MyTrees[f].MyCoevolData[t].e1[i],  MyTrees[f].MyCoevolData[t].sparseChronomat, MyTrees[f].MyCoevolData[t].nbranches, MyTrees[f].MyCoevolData[t].e1M, MyTrees[f].MyCoevolData[t].mask[i], MyTrees[f].MyCoevolData[t].mask[j]);
						else{
							MyTrees[f].MyCoevolData[t].e1M =  MyTrees[f].MyCoevolData[t].e1[i];
							for (int z=0; z<MyTrees[f].MyCoevolData[t].nbranches; z++){
								MyTrees[f].MyCoevolData[t].e1M[z] *= (int)MyTrees[f].MyCoevolData[t].mask[i][z] * (int)MyTrees[f].MyCoevolData[t].mask[j][z];
							}

						}
						if (ntree==1 && maxforest==1)
							do_e1M=1;
					}

			  		/*
			  			Computation of a "compressed" vector of possible pairs (given M and e1)
			  		*/

					set_vect_classes_mask(MyTrees[f].MyCoevolData[t].e1M, MyTrees[f].MyCoevolData[t].e1[j], MyTrees[f].MyCoevolData[t].mask[i], MyTrees[f].MyCoevolData[t].mask[j], MyTrees[f].MyCoevolData[t].nbranches, &nclasses, &ne1M, &ne2, MyTrees[f].MyCoevolData[t].branch_lengths, &nbranch_lengths, t);

			  		nn2 = sumv(ne2, nclasses); /* new cumulative sum of vecteur e2, identical to n2 a priori */
			  		if(nn2 != n2) fprintf(stderr, "main: this should not happen, exiting\n"),exit(4);

			  		/*
			  			Computation of the max distrib and (re)allocate memory for the distributions
			  		*/
					ldistrib = n1*n2+1;
					Init_veck_distrib( nclasses, &nclassmax, &nveck, ldistrib,  &max_ldistrib, &distrib );

					/*
						Generation of the distributions
					*/

					/* Calculating the multinomial dimension: weak composition of n into k parts	*/
					long long dimension=binomial(nn2-nclasses+1, nclasses-1);

					/* If the dimension is too big, we can't finish the multinomial -> no convergence!	*/
					/* Therefore, we use an approached multinomial, kind of like a MCMC approach	*/
					/* Afterwards we need to normalize the vector to sum to 1	*/
					if (dimension>50000000){
						generate_approached_multinomial(nveck, nbranch_lengths, nclasses, nn2, nn2, 0, distrib, ne1M, 0.0L, 0.0L);
						normalize_a_vector(distrib,ldistrib);
					}
					if (dimension<=50000000){
						generate_multinomial(nveck, nbranch_lengths, nclasses, nn2, nn2, 0, distrib, ne1M, 0.0L, 0.0L);
					}

			  		if (output_dist)  PrintDist(distrib, ldistrib, mat_type );
			  		/*
			  			# of observed pairs (depends on M)
			  		*/
			  		n_obs_pairs = iprodsca(ne1M, ne2, nclasses);
			  		mymax=MaxDistrib( distrib, ldistrib  );

					MyTrees[f].pvalue_of_pair=ComputePval( distrib, ldistrib, n_obs_pairs);

					if (MyTrees[f].pvalue_of_pair==0.0){

						double sum_distrib = 0.0L;

						for (int z=0; z<ldistrib; z++){
							sum_distrib += distrib[z];
						}

						if (sum_distrib != 1){
							fprintf(stderr, "WARNING: for pair %d vs %d, there seems to be no observable/possible %s\n", i+1, j+1, (mat_type==0)?"co-occurences":(mat_type==1)?"chronologies":"co-occurences or chronologies");
							fprintf(stderr, "Setting pvalue to 1.0\n");
							MyTrees[f].pvalue_of_pair=1.0;
						}
					}

					/*
						Writing the results to stdout. If there is an event array, replaces the e0, e1, etc. by the ids
					*/

					if (input_names){
						printf("%s vs %s t%-4d: nobs= %d pval= %.5g [max= %d]\n", input_ids[i],input_ids[j],maxforest>1?f:t+1,n_obs_pairs, MyTrees[f].pvalue_of_pair, mymax );
					}else printf("e%-4d vs e%-4d t%-4d: nobs= %d pval= %.5g [max= %d]\n", i+1,j+1,maxforest>1?f:t+1,n_obs_pairs, MyTrees[f].pvalue_of_pair, mymax );

					/*
						Writing the results in the output file
					*/

					writeInFile(ntree, maxforest, input_names, fptr, input_ids, mat_name, i, j, t, f, MyTrees[f].pvalue_of_pair, n_obs_pairs, mymax);
					if (MyTrees[f].pvalue_of_pair <= seuil_pval)
						writeInFile(ntree, maxforest, input_names, fptrs, input_ids, mat_name, i, j, t, f, MyTrees[f].pvalue_of_pair, n_obs_pairs, mymax);

					/* if the correlation is significant, tagging the event as -2 */
					/* will be used if needed to output the tree with only the    */
					/* significant pairs. Was previously outside of loop		  */

					if (output_newick == 1 && MyTrees[f].pvalue_of_pair <= seuil_evt)
			  		{
			  			MyTrees[f].MyCoevolData[t].v_evts[i] = -2;
			  			MyTrees[f].MyCoevolData[t].v_evts[j] = -2;	/* will output name of pair in tree */
					}

			  		/*
			  			Combine with all previous trees
			  		*/
			  		if(mymax){
				  		total_obs += n_obs_pairs;
				  		distrib_cum = CombineDistrib( distrib_cum, &max_cum,  distrib, ldistrib  );
				  	}


			  		/*
			  			Free tmp vectors
			  		*/

			  		if (ne1M) {free(ne1M); ne1M=NULL;}
					if (ne2)  {free(ne2); ne2=NULL;}
					if (nbranch_lengths){free(nbranch_lengths);nbranch_lengths=NULL;}

				} /* end of for (t=0; t<=ntree; t++) */

				/* We don't need to compute e1.M at each round - when multiple trees, we need to generate all e1M before flagging	*/
				if (ntree>1 && maxforest==1)
					do_e1M=1;

			  	/*
			  		Compare obs with Exp and compute a pval
			  	*/
				if (ntree>1){
				  	if (total_obs >= 0 || output_00==1) {

						double tot=0;
					  	for(k=0;k<max_cum; k++){
							printf("%f, ", distrib_cum[k]);
						  	tot += distrib_cum[k];
						  }
						  printf("\n");
			  			if(total_obs>0 && tot<0.99){
				  			printf("Bug: tot is %lf\n", tot);
				  			exit(3);
				  		}

				  		if(total_obs>0)
				  			pval = ComputePval( distrib_cum, max_cum, total_obs);                              /* the p-value is the sum of all proba of the distribution >= n_obs_pairs */
				  		else
				  		{
				  			pval=1;
				  			max_cum=1;
				  		}

				  		if (ntree>1 && pval <= seuil_pval){
				  			fprintf(stdout, "-------\ne%-4d vs e%-4d all_t: nobs= %d pval= %.5g [max= %d]\n", i+1,j+1,total_obs,pval, max_cum-1 );
						}
						if (ntree>1){
							writeInFile(1, maxforest, input_names, fptr_m, input_ids, mat_name, i, j, 0, 0, pval, total_obs, max_cum-1);
							if (pval <= seuil_pval){
								writeInFile(1, maxforest, input_names, fptr_ms, input_ids, mat_name, i, j, 0, 0, pval, total_obs, max_cum-1);
							}
						}

					} /* end of n_obs_pairs > 0 || output_00==1) */
				}
			} /* end of for f=0; f<maxforest	*/

			/*
				When we have a forest - computation of the 95% confidence intervals around the p-values
			*/
			if (maxforest>1)
			{
				fprintf(stdout, "Estimations and confidence intervals for the whole forest\n");
				lkDist = ComputeLikelihoodDistribution(MyTrees, maxforest, 1, i, j);
			}
			cntj++;
		}	/* end for j */
		cntj=0;
		cnti++;
	} /* end for i */

	fclose(fptr);
	fclose(fptrs);
	if(ntree>1){
		fclose(fptr_m);
		fclose(fptr_ms);
	}

	if (output_newick == 1) {

		for (f=0; f<maxforest; f++){

			for(t=0;t<ntree;t++)
			{
				change_names_with_events(MyTrees[f].MyCoevolData[t].root, MyTrees[f].MyCoevolData[t].v_evts, nevt);
				PrintTree(MyTrees[f].MyCoevolData[t].root, NULL , 1);
			}
		}

	}

	/*
		Freeing memory
	*/

	if (veck) free(veck);
	if (nveck) free(nveck);
	if (distrib) free(distrib);
	free_sFacto();

	if (input_names){
        free(input_ids);
        free(input_names);
    }

	if (input_events_file){

		for (i=0; i<input_events->length; i++){
			free(input_events->inpe2[i]);
		}
		free(input_events->e1);
		free(input_events->e2L);
		free(input_events->inpe2);
		free(input_events);
	}

	for (f=0; f<maxforest; f++){
		FreeCoevolData(MyTrees[f].MyCoevolData, ntree, 0, 0);
		if (MyTrees[f].o2)free(MyTrees[f].o2);
		if(MyTrees[f].nbfork<=20){
			free(MyTrees[f].forkVector);
			free(MyTrees[f].maxForkVector);
		}
	}
	free(MyTrees);

	return 0;
}
