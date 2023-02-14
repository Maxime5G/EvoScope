/*
	Copyright (C) 2015-2020 K Bedhenna G Achaz M Godforid & S Brouillet Joel Pothier William Belaid

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

 	for more information, please contact guillaume achaz <achaz@abi.snv.jussieu.fr>/<gachaz@gmail.com>
*/
/**
\file 	epocs_mcmc.c
\author     G Achaz, M Godfroid, J Pothier, K. Bedhenna and madamesophie and William and Patrice and M Godfroid
\date 	Sept 2015 -> June 2022
**/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"
#include "coevol.h"
#include "MCMC.h"

#include <unistd.h>     /* pour getopt  */

#include <sys/stat.h> 		/* for the creation of directories and files */
#include <sys/types.h> 		/* for the creation of directories and files */

/*	Function to count the number of event chains (number of tvectors)	*/
void count_vector_types(char **ts, int *count){
	while (*ts){
		*count+=1;
		ts++;
	}
	return;
}

/* Function to calculate the proba that >1 event occurs on a given branch	*/
/* k = branch length, m = rate	*/
double probaMoreThanOneEvent(double k, double m){
	return (1-exp(-(m*k))-(m*k)*exp(-(m*k)));
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
	fprintf(stderr,"     tree_file: tree file needs to be provided in a modified Newick format\n");
	fprintf(stderr,"                with events provided at the locations usually dedicated to bootstrap values (see README).\n");
	fprintf(stderr,"                if not provided, read the tree from stdin.\n");

	fprintf(stderr," [second usage: many trees with the same events]\n");
	fprintf(stderr,"     all_trees_file: a file with all info on each tree each line as 2 columns\n");
	fprintf(stderr,"                     1: a tree_file (see above) \n");
	fprintf(stderr,"                     2: an outgroup (or \"\") for no outgroup) \n");
	fprintf(stderr,"     all trees must have all their events in the same order !\n");

	fprintf(stderr,"   options are:\n\n");

	fprintf(stderr," [input options] \n");
	fprintf(stderr,"   -1 <event 1>: set the first event id (integer) to be tested \n");
	fprintf(stderr,"   -2 <event 2>: set the second event id (integer) to be tested \n");
    fprintf(stderr,"   -N <integer>: set the number of rounds for the MCMC chain. Default is 100,000\n");
	fprintf(stderr,"   -S <sampling>: set the intervals for sampling the MCMC chain. Default is sampling at every iteration.\n\n");

	fprintf(stderr," [output options] \n");
	fprintf(stderr,"   -O: <output>: set a value for the output folder and files prefixes [default: out]\n");
	fprintf(stderr,"   -v: (verbose) outputs running progression of the program\n");
	fprintf(stderr,"   -V: (very verbose) outputs running progression of the program and more...\n\n");

	fprintf(stderr,"\n");

}

int main( int argc , char **argv)
{
	char *s=NULL,
	     *outgroup=NULL;

	Node *init_tree=NULL;

	char output_prefix[1000]="out"; 	/* default output prefix           */
    char *fname=NULL; 					/* output file name 			   */
    FILE *fptr; 						/* pointer to open the output file */
	char * infile=NULL;

	int verbose=0;

	int nbleaves=0,
	    nevt=-1,
	    carg,
	    evti,
	    evtj,
		veryverbose=0,
	    output_tree = 0,    /*1= in newick, 2=in ascii*/
	    output_state_vector = 0,
		choice_e1=0,			/* if the user wants to input a particular event to test - e1 */
		choice_e2=0,			/* if the user wants to input a particular event to test - e2 */
		ntree=1,				/* the number of trees in input	*/
		t,						/* the tree counter				*/
		f,						/* the forest counter	*/
		sampling=1;				/* the sampling frequency	*/

	struct Trees *MyTrees=NULL;	/* initializing the trees structure	*/
	int *forestntree=NULL;

	double rate[2][4];			/* initializing the rate vector	*/
    int Nrounds=100000;			/* number of MCMC rounds	*/

	/*
		Recursive catching of command line arguments
	*/
	while ((carg = getopt(argc, argv, "1:2:VO:N:ehtTvf:S:")) != -1)
	{

		switch ((char)carg) {

			case 'v':
				verbose = 1;
				break;

			case 'V':
		   		verbose = 2;
		   		veryverbose = 1;
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

			case 'O':
				if (!(strcpy(output_prefix, optarg))){
					fprintf(stderr, "Cannot read output prefix or too long name (reading %s)\n", optarg);
					fprintf(stderr, "    setting it to out (default)\n");
				}
				break;

            case 'N':
            if (sscanf(optarg, "%d", &(Nrounds)) != 1){
                fprintf(stderr, "cannot read the number of rounds. Please enter an integer, exiting... \n"), exit(1);
                fprintf(stderr, "    setting it to default (100,000 rounds)\n");
            }
            break;

			case 'S':
			if (sscanf(optarg, "%d", &(sampling)) != 1){
                fprintf(stderr, "cannot read the sampling number. Please enter an integer, exiting... \n"), exit(1);
                fprintf(stderr, "    setting it to default (sampling every iteration)\n");
            }
            break;

			case 'e':
		   		output_state_vector = 1;
		   		break;

			case 't':
				output_tree = 1;	         /* outputs the tree */
				break;

			case 'T':
				output_tree = 2;	         /* outputs the tree */
				break;

			case 'h':
				PrintHelp(argv[0]), exit(0);

			case 'f':
				infile=optarg;
				break;

			default:
				fprintf(stderr,"option:<%c> unknown... bye\n",(char)carg);
				exit(1);
	   }

   }

   if(argc == 1)
	   fprintf(stderr, "Usage: ./epocs_mcmc [options] outgroup [tree_file]\nfor more help ./epocs -h\n"), exit(2);

   srand(time(NULL));


   /* ---------------------------------------------------------------------------------------------- */
   /* Check that, if the user chose events to run, e1 < e2 (given the loop) 	  	 		  		  */
   /* This also verifies that if the user mistakenly input -2 only, choice_e2 is copied to choice_e1 */

	if ((!choice_e1 && choice_e2) || (!choice_e2 && choice_e1) || (!choice_e1 && !choice_e2)){
        fprintf(stderr, "You did not select two events, setting by default 1 vs 2\n");
        choice_e1=1;
        choice_e2=2;
	}

	if ((choice_e1 && choice_e2) && (choice_e2 <= choice_e1)){
		int tmp_choice=choice_e1;
		choice_e1=choice_e2;
		choice_e2=tmp_choice;
	}

	if ((choice_e1 && choice_e2) && (choice_e1 == choice_e2))fprintf(stderr,"ERROR: Please input different numbers for the event choice\n"), exit(1);

   /* ---------------------------------------------------------------------------------------------- */

	if(!(infile)){
		if( strlen( argv[optind+0]) > 0 )
			outgroup = argv[optind+0];
	}
	/*
		Retrieving the tree informations. The ParseInputForest reads the tree or the treefile and determines whether we have
		a forest, multiple independent trees or a single newick tree. All informations are stored in a Trees structure
	*/

	MyTrees = ParseInputForest( argc, optind, argv, infile, verbose, &nevt, &ntree, &forestntree, 0, 0 );

	/*
		Setting the total number of events in nevt. Make sure that the number of events is consistent across multiple trees accordingly.
		Taking the opportunity to count the maximum number of trees in forest
	*/

	int maxforest=0;
	int forestlength=0;
	for (t=0; t<ntree; t++){
		forestlength=forestntree[t];
		if (forestlength>maxforest) maxforest=forestlength;
		for (f=0; f<forestntree[t]; f++){
			if (nevt<0){nevt=MyTrees[f].MyCoevolData[t].nevt;}
			else{if (MyTrees[f].MyCoevolData[t].nevt != nevt){fprintf(stderr, "not the same number of events for tree %d! Exiting.", t); exit(3);}}
		}
	}

	/*
		VERBOSE - Displays some information about the input dataset on stderr
	*/

	if(verbose){
		fprintf(stderr, "ntrees in forest=%d\n", maxforest);
		fprintf(stderr, "ntree in input=%d\n", ntree);
		fprintf(stderr, "number of events = %d\n", nevt);
		for (t=0; t<ntree; t++){
			for (f=0; f<forestntree[t]; f++){
				fprintf(stderr, "longest branch length of tree [%d] in forest [%d] = %f\n", t, f, MyTrees[t].MyCoevolData[f].longestbranch);
				fprintf(stderr, "nbranches of tree [%d] in forest [%d] = %d\n", t, f, MyTrees[t].MyCoevolData[f].nbranches);
				fprintf(stderr, "nbleaves of tree [%d] in forest [%d] = %d\n", t, f, MyTrees[t].MyCoevolData[f].nbleaves);
				print_arbre(MyTrees[t].MyCoevolData[f].root, stderr);
			}
		}
		fprintf(stderr,"\n");
	}
    if (maxforest>1){fprintf(stderr, "Please input only a single tree. Exiting...\n"), exit(1);}
    if (ntree>1){fprintf(stderr, "Please input only a single tree. Exiting...\n"), exit(1);}

	/* --- FILE INITIALIZATIONS ---	*/

	/* Opening the output file and writing the header */
	fname = initializeNewFile(output_prefix, '\0', "", 7);				/* create the string for the output file and */
	fptr=fopen(fname, "w");													/* opening it */

	/*
		VERBOSE - printing the changes made to some nodes (e.g., missing events are now 0/0/...)
	*/

	if(verbose){
		fprintf(stderr,"->Verifying nodes events:\n");
		for (t=0; t<ntree; t++){
			for (f=0; f<forestntree[t]; f++){
				fprintf(stderr,"Numbers of species in the tree %d : %d (excluding outgroup)\n",t, MyTrees[f].MyCoevolData[t].nbleaves);
				if (MyTrees[f].MyCoevolData[t].changed != 0)
					fprintf(stderr,"->Verification: %d nodes have been initialised with list of %d events (filled with 0)\n",MyTrees[f].MyCoevolData[t].changed, nevt);
			}
		}
	}

	/*
		VERBOSE - printing the ASCII tree on stderr
	*/

	if (output_tree==2){
		for (t=0; t<ntree; t++){
			for (f=0; f<forestntree[t]; f++){
				fprintf(stderr,"->Tree %d from forest %d resume:\n", t, f);
				print_arbre(MyTrees[f].MyCoevolData[t].root, stderr);
			}
		}
	}

	/* --- MAIN LOOP --- */

	fprintf(stderr, "/*\n\t** epocs_mcmc header **\n\n");
	fprintf(stderr, "\toutgroup is '%s'\n", outgroup);
	fprintf(stderr, "\ttreefile is '%s'\n", (infile)?infile:((argc == optind + 2)?argv[optind+1]:NULL));
	fprintf(stderr, "*/\n");

    evti = choice_e1-1<0?0:choice_e1-1;
    evtj = choice_e2-1<1?1:choice_e2-1;
	fprintf(stderr, "Comparing e_%d vs. e_%d\n",evti+1,evtj+1);

	MyTrees[0].MyCoevolData[0].theTVector.nbfork = set_evt_type_count_fork_gaps(MyTrees[0].MyCoevolData[0].root, nevt, evti, evtj, verbose, output_state_vector);
	MyTrees[0].nbfork+=MyTrees[0].MyCoevolData[0].theTVector.nbfork;
	if(verbose) fprintf(stderr, "nbfork: %d for tree: %d\n", MyTrees[0].MyCoevolData[0].theTVector.nbfork, maxforest>1?f:t);

	MyTrees[0].forkVector = calloc(ntree, sizeof(int));
	MyTrees[0].maxForkVector = calloc(ntree, sizeof(int));

	MyTrees[0].MyCoevolData[0].theTVector.current_vector = calloc((size_t)MyTrees[0].MyCoevolData[0].nbranches, sizeof(size_t)); // allocating the initial tvector

	MyTrees[0].MyCoevolData[0].theTVector.tvector = genere_single_vector(MyTrees[0].MyCoevolData[0].root, MyTrees[0].MyCoevolData[0].nbranches, MyTrees[0].MyCoevolData[0].theTVector.nbfork, (char **)NULL, &MyTrees[0].MyCoevolData[0].theTVector.nvect); 	// filling all tvectors

	MyTrees[0].maxForkVector[0] = MyTrees[0].MyCoevolData[0].theTVector.nvect; 										 // storing the number of tvectors

	if(verbose){
		fprintf(stderr, "tree: %d ; #branches: %d ; #cooccurences: %d #trees+events+root: %d\n", maxforest>1?f:t,  MyTrees[0].MyCoevolData[0].nbranches, MyTrees[0].MyCoevolData[0].theTVector.nbfork, (int)pow(2,MyTrees[0].MyCoevolData[0].theTVector.nbfork)*3);
		print_vector_types(MyTrees[0].MyCoevolData[0].theTVector.tvector, MyTrees[0].MyCoevolData[0].nbranches, stderr);
	}

	MyTrees[0].IS_epocs[0] = MyTrees[0].IS_epocs[1] = 0;									// Initial states
	memcpy((void *)MyTrees[0].final_rates_epocs, (void *)rate, (size_t)8*sizeof(double));	// initializing the rate matrix
	MyTrees[0].imax_epocs=0;																// Initializing the initial state value of e1 (updated if -I chosen)
	MyTrees[0].jmax_epocs=0;																// Initializing the initial state value of e2 (updated if -I chosen)
	MyTrees[0].o2=malloc(sizeof(int)*ntree);												// Allocating memory for the ids of the tvectors
	MyTrees[0].lik_epocs=0;																	// Initializing the maximum likelihood value

	int total_acceptance=0;	/* Not really used now...	*/

	int w_value=5;			/* When MCMC exploring, size of the window	*/

	fprintf(stderr, "Processing...\n\n");
	MCMC2(MyTrees[0].MyCoevolData, MyTrees[0].IS_epocs, Nrounds, w_value, &total_acceptance, sampling, fptr);

	fclose(fptr);

	fprintf(stderr, "Finished! Results written in file: %s!\n\n", fname);

	for (f=0; f<maxforest; f++){
		FreeCoevolData(MyTrees[f].MyCoevolData, ntree, 0, 0);
		if (MyTrees[f].o2)free(MyTrees[f].o2);
	}
	free(MyTrees);

	if (s) free(s);
	if (init_tree) FreeMultiTree(init_tree, nbleaves);
	return 0;
}
