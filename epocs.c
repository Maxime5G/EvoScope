/*
	Copyright (C) 2015-2020 K Bedhenna G Achaz & S Brouillet Joel Pothier William Belaid

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
\file 	epocs.c
\author     G Achaz, J Pothier, K. Bedhenna and madamesophie and William and Patrice and M Godfroid
\date 	Sept 2015 -> May 2021
**/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tree.h"
#include "coevol.h"

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
	fprintf(stderr,"   -R <filename>: input filename containing significant events to run (epics format) \n\n");

	fprintf(stderr," [output options] \n");
	fprintf(stderr,"   -O: <output>: set a value for the output folder and files prefixes [default: out]\n");
	fprintf(stderr,"   -t: (tree) output the input tree on stderr\n");
	fprintf(stderr,"   -e: output on stderr the vector of states)\n");
	fprintf(stderr,"   -v: (verbose) outputs running progression of the program\n");
	fprintf(stderr,"   -V: (very verbose) outputs running progression of the program and more...\n\n");

	fprintf(stderr," [ML options] \n");
	fprintf(stderr,"   -C   : Maximum number of co-occurrences to consider in the analysis. Default is 10\n");
	fprintf(stderr,"   -I   : evaluate the Initial State (otherwise, assume <0,0>)\n");

	fprintf(stderr," [scenario selection] \n");
	fprintf(stderr,"   -S # : depending on the value of #, a different scenario is evaluated\n");
	fprintf(stderr,"      1 : 1 parameter: mu\n");
	fprintf(stderr,"      i : 2 parameters: mu1, mu2 (H0)\n");
	fprintf(stderr,"      x : 2 parameters: mu, lambda -- modified rates are mu*lambda --\n");
	fprintf(stderr,"      u : 3 parameters: mu1, nu1, mu2\n");
	fprintf(stderr,"      U : 3 parameters: mu1, mu2, nu2\n");
	fprintf(stderr,"      X : 3 parameters: mu, lambda, kappa --mu*lambda*kappa--\n");
	fprintf(stderr,"      r : 3 parameters: mu1, mu2, lambda (reciprocal induction: E1<->E2) -- modified rates are mu[12]*lambda --\n");
	fprintf(stderr,"      a : 3 parameters: mu1, mu2, mu2* (one-way induction: E1->E2)\n");
	fprintf(stderr,"      b : 3 parameters: mu1, mu1*, mu2 (one-way induction: E2->E1)\n");
	fprintf(stderr,"      l : 4 parameters: mu1, mu1*, mu2, mu2* (induction both ways)\n");
	fprintf(stderr,"      I : 4 parameters: mu1, mu2, nu1, nu2 (state dependance, no induction)\n");
	fprintf(stderr,"      R : 4 parameters: mu1, mu2, lambda, kappa (reciprocal induction: E1<->E2, and state dependance) --mu[12]*lambda*kappa--\n");
	fprintf(stderr,"      A : 5 parameters: mu1, mu2, mu2*, kappa1, kappa2 (one-way induction: E1->E2)\n");
	fprintf(stderr,"      B : 5 parameters: mu1, mu1*, mu2, kappa1, kappa2 (one-way induction: E2->E1)\n");
	fprintf(stderr,"      L : 6 parameters: mu1, mu1*, kappa1, mu2, mu2*, kappa2 (induction + state dependance)\n");
	fprintf(stderr," default: 8 parameters: mu1, mu1*, nu1, nu1*, mu2, mu2*, nu2, nu2*\n");
	fprintf(stderr,"\n");

}

/* Print a phylogenetic tree in ASCII on stdout	*/
void print_arbre2(Node *n, FILE *f)
{
	static int branchnum = -2, recurs = 0;
	int i;

	branchnum++;
	recurs++;

	if (branchnum == -1) {
		fprintf(f,"[%3d]\n",branchnum);
	}
	else {
		for (i = 0; i < recurs; i++)
			fprintf(f,"\t");
		fprintf(f,"|->[%3d]\n",branchnum);
	}

	if (n->nbdesc > 0)
		for (i = 0; i < n->nbdesc; i++)
			print_arbre (n->descs[i], f);

	recurs--;

	return;
}

/* -----------------------------------------------------------------	*/
int main( int argc , char **argv)
{
	char *s=NULL,
	     opt_IS=0,
	     *outgroup=NULL;

	char scenario='8';               // the selected scenario. By default, it is the full-model (8 params)

	Node *init_tree=NULL;

	char output_prefix[1000]="out"; 	/* default output prefix           */
    char *fname=NULL; 					/* output file name 			   */
    FILE *fptr; 						/* pointer to open the output file */
	char *input_events_file=NULL;		/* pointer to open input event to run */
	int **input_events=NULL;			/* array of events to run */
	int filelength=0;					/* length of the array of events to run */
	char * infile=NULL;

	double maxrate=0.0;

	int verbose=0;

	int nbleaves=0,
		i, j,
	    nevt=-1,
	    carg,
	    evti,
	    evtj,
		veryverbose=0,
	    output_tree = 0,    /*1= in newick, 2=in ascii*/
	    output_newick = 0,
	    output_state_vector = 0,
		choice_e1=0,			/* if the user wants to input a particular event to test - e1 */
		choice_e2=0,			/* if the user wants to input a particular event to test - e2 */
		ntree=1,				/* the number of trees in input	*/
		t,						/* the tree counter				*/
		f,						/* the forest counter	*/
		forkmax=10;				/* maximum number of forks allowed	*/

	struct Trees *MyTrees=NULL;	/* initializing the forest	*/
	int *forestntree=NULL;

	double rate[2][4];

	/*
		Recursive catching of command line arguments
	*/
	while ((carg = getopt(argc, argv, "1:2:VO:R:ehmtTvS:If:C:")) != -1)
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

			case 'R':
                input_events_file=(char *)malloc((strlen(optarg)+1)*sizeof(char));
				if (!(strcpy(input_events_file, optarg))){
                    fprintf(stderr, "Cannot read input file or too long name (reading %s)\n", optarg);
                    fprintf(stderr, "    setting it to none\n");
                }
                break;

			case 'C':
				if (sscanf(optarg, "%d", &(forkmax)) != 1){
					fprintf(stderr, "cannot read your co-occurrence level. Please enter an integer. \n");
					fprintf(stderr, "setting it to default (10)\n");
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

			case 'S':
				scenario = *optarg;

				break;

			case 'I':
				opt_IS = 1;
				break;

			case 'f':
				infile=optarg;
				break;

			default:
				fprintf(stderr,"option:<%c> unknown... bye\n",(char)carg);
				exit(1);
	   }

   }

   if(argc == 1)
	   fprintf(stderr, "Usage: ./epocs [options] outgroup [tree_file]\nfor more help ./epocs -h\n"), exit(2);

   /* ---------------------------------------------------------------------------------------------- */
   /* Check that, if the user chose events to run, e1 < e2 (given the loop) 	  	 		  		  */
   /* This also verifies that if the user mistakenly input -2 only, choice_e2 is copied to choice_e1 */

	if ((input_events_file) && (choice_e1 || choice_e2))fprintf(stderr,"ERROR: Choices not permitted when inputting an events file\n"), exit(1);

	if (!choice_e1 && choice_e2){
		choice_e1=choice_e2;
		choice_e2=0;
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
		printf("\n");
	}

	/* --- FILE INITIALIZATIONS ---	*/

	/* Opening the output file and writing the header */
	fname = initializeNewFile(output_prefix, scenario, "", 0);				/* create the string for the output file and */
	fptr=fopen(fname, "w");													/* opening it */
	fprintf(fptr, "%-10s\t%-10s\t%-10s\t%-5s\t%-14s\t%-14s\t%-4s\t%-4s\t%-4s\t%-4s\t%-4s\t%-4s\t%-4s\t%-4s\n", "Event1", "Event2", "Tree", "Model", "lnML", "ML", "mu1", "mu1star", "nu1", "nu1star", "mu2", "mu2star", "nu2", "nu2star");

	/* Verifying if the user gave an input file of significant events to run. If yes, store the events */
	/* in an nx2 array. */
	if (input_events_file){
		input_events = createArrayEventsToRun(input_events_file, &filelength, verbose); 		/* if -R chosen, creates the array of event ids to run. */
	}

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

	fprintf(stderr, "/*\n\t** epocs header **\n\n");
	fprintf(stderr, "\toutgroup is '%s'\n", outgroup);
	fprintf(stderr, "\ttreefile is '%s'\n", (infile)?infile:((argc == optind + 2)?argv[optind+1]:NULL));
	fprintf(stderr, "\tchosen scenario is '%c'\n", scenario);
	fprintf(stderr, "*/\n");

	/*
		The idea of the loop is to allow the user to input either a file with signif pairs (stored in array).
		If there is no file, then goes along the loop. For this reason, the for loop has been
		rewritten from the inside to a while loop.
	*/

	int count_file=0;
	int count_i=0;
	int count_j=1;

	while (1)
	{
		if (input_events)
		{
			if (count_file == filelength)
				break;
			evti = input_events[count_file][0]-1;
			evtj = input_events[count_file][1]-1;
			count_file+=1;
		}else
		{
			if (count_i == nevt-1)
				break;
			evti = count_i;
			evtj = count_j++;
			if (count_j == nevt)
			{
				count_i += 1;
				count_j = count_i+1;
			}
		}
		if (choice_e1 && ((choice_e1 != evti+1) && (choice_e1 != evtj+1))) continue;
		if (choice_e2 && ((choice_e2 != evtj+1))) continue;
		for (f=0; f<maxforest; f++)
			MyTrees[f].nbfork=0;

		/*
			For the user, it is better to display the events starting at 1 and not 0
		*/
		printf("Comparing e_%d vs. e_%d\n",evti+1,evtj+1);

		/*
			COUNTING the number of forks (i.e., co-occurrences), the main limiting factor of the method so far
			Storing them in nbfork for each tree
		*/

		for (f=0; f<maxforest; f++){
			for (t=0; t<ntree; t++){
				// MyTrees[f].MyCoevolData[t].theTVector.nbfork = set_evt_type_count_fork(MyTrees[f].MyCoevolData[t].root, nevt, evti, evtj, verbose, output_state_vector);
				MyTrees[f].MyCoevolData[t].theTVector.nbfork = set_evt_type_count_fork_gaps(MyTrees[f].MyCoevolData[t].root, nevt, evti, evtj, verbose, output_state_vector);
				MyTrees[f].nbfork+=MyTrees[f].MyCoevolData[t].theTVector.nbfork;
				if(verbose) fprintf(stderr, "nbfork: %d for tree: %d\n", MyTrees[f].MyCoevolData[t].theTVector.nbfork, maxforest>1?f:t);
			}
		}

		if (maxforest==1){
			if (verbose) fprintf(stderr,"number of forks %d for couple events %d and %d \n",MyTrees[0].nbfork,evti+1,evtj+1);
			/* Need to count the number of forks for each of t in maxtree, i.e., there are one tree in each t files */

			/* forkVector and maxForkVector necessary to explore all combinations of tvectors efficiently	*/
			MyTrees[0].forkVector = calloc(ntree, sizeof(int));
			MyTrees[0].maxForkVector = calloc(ntree, sizeof(int));

			/*
				TVECTORS - generate the transition vectors for each tree. In other words, I'm saving the order of mutations on the branches of the tree in arrays of int
				Essentially: 1 on a branch = e1, 2 = e2, 3 = co-occ e1->e2, 4 = co-occ e2->e1.
				When there are co-occurrences -> duplicate the vectors, such that I have one representative w/ a 3 and another with a 4. (see genere_all_vector_types)
			*/

			if (MyTrees[0].nbfork <= forkmax){	/* If I have less than the maximum number of forks allowed - storing every events	*/
				for (t=0; t<ntree; t++){
					MyTrees[0].MyCoevolData[t].theTVector.current_vector = calloc((size_t)MyTrees[0].MyCoevolData[t].nbranches, sizeof(size_t)); // allocating the initial tvector

					MyTrees[0].MyCoevolData[t].theTVector.tvector = genere_all_vector_types(MyTrees[0].MyCoevolData[t].root, MyTrees[0].MyCoevolData[t].nbranches, MyTrees[0].MyCoevolData[t].theTVector.nbfork, (char **)NULL, &MyTrees[0].MyCoevolData[t].theTVector.nvect); 	// filling all tvectors

					MyTrees[0].MyCoevolData[t].theTVector.current_vector = MyTrees[0].MyCoevolData[t].theTVector.tvector;

					MyTrees[0].MyCoevolData[t].theTVector.nvect=0; 																	 // recounting the number of tvectors
					count_vector_types(MyTrees[0].MyCoevolData[t].theTVector.tvector, &MyTrees[0].MyCoevolData[t].theTVector.nvect); // actual counting
					MyTrees[0].maxForkVector[t] = MyTrees[0].MyCoevolData[t].theTVector.nvect; 										 // storing the number of tvectors
				}
			}

			/* If I have more than the allowed number of co-occurrences - generate random transition vectors, at most 1000	*/
			else{
				if(verbose){
					fprintf(stderr, "for pair %d %d -> nbfork = %d\n", evti, evtj, MyTrees[0].nbfork);
					fprintf(stderr, "the number of forks is higher than the maximum allowed (%d)\n", forkmax);
					fprintf(stderr, "Epocs will generate random tvectors to explore + mandatory orders 1->2 and 2->1 at each co-occurrence\n");
				}
				for (t=0; t<ntree; t++){
					MyTrees[0].MyCoevolData[t].theTVector.current_vector = malloc(sizeof(size_t)*MyTrees[0].MyCoevolData[t].nbranches);

					MyTrees[0].MyCoevolData[t].theTVector.tvector = genere_dual_vector_types(MyTrees[0].MyCoevolData[t].root, MyTrees[0].MyCoevolData[t].nbranches, MyTrees[0].MyCoevolData[t].theTVector.nbfork, (char **)NULL, &MyTrees[0].MyCoevolData[t].theTVector.nvect);

					add_random_vector_types_ok(MyTrees[0].MyCoevolData[t].theTVector.tvector, MyTrees[0].MyCoevolData[t].nbranches, 1000);

					MyTrees[0].MyCoevolData[t].theTVector.current_vector = MyTrees[0].MyCoevolData[t].theTVector.tvector;
					MyTrees[0].MyCoevolData[t].theTVector.nvect=1000;
					MyTrees[0].maxForkVector[t] = MyTrees[0].MyCoevolData[t].theTVector.nvect;
				}
			}
		}

		if (maxforest>1){
			if(ntree>1){
				fprintf(stderr, "Please input either 1 forest, 1 tree or multiple independent trees. exiting\n"); exit(3);
			}

			// I need to count the number of forks independently for each tree in each forest

			/*
				SAME PROCEDURE AS ABOVE - counting forks, generating tvectors etc. One vector of each for each tree in the forest
			*/

			for (t=0; t<ntree; t++){
				for (f=0; f<forestntree[t]; f++){

					if (!(MyTrees[f].forkVector)){
						MyTrees[f].forkVector = calloc(ntree, sizeof(int));
						MyTrees[f].maxForkVector = calloc(ntree, sizeof(int));
					}

					if (MyTrees[f].MyCoevolData[t].theTVector.nbfork <= forkmax){
						MyTrees[f].MyCoevolData[t].theTVector.current_vector = calloc((size_t)MyTrees[f].MyCoevolData[t].nbranches, sizeof(size_t)); // allocating the initial tvector

						MyTrees[f].MyCoevolData[t].theTVector.tvector = genere_all_vector_types(MyTrees[f].MyCoevolData[t].root, MyTrees[f].MyCoevolData[t].nbranches, MyTrees[f].MyCoevolData[t].theTVector.nbfork, (char **)NULL, &MyTrees[f].MyCoevolData[t].theTVector.nvect); // filling all tvectors

						MyTrees[f].MyCoevolData[t].theTVector.current_vector = MyTrees[f].MyCoevolData[t].theTVector.tvector;

						MyTrees[f].MyCoevolData[t].theTVector.nvect=0; // recounting the number of tvectors
						count_vector_types(MyTrees[f].MyCoevolData[t].theTVector.tvector, &MyTrees[f].MyCoevolData[t].theTVector.nvect); // actual counting
						MyTrees[f].maxForkVector[t] = MyTrees[f].MyCoevolData[t].theTVector.nvect; // storing the number of tvectors
					}

					else{
						if(verbose){
							fprintf(stderr, "for pair %d %d -> nbfork = %d\n", evti, evtj, MyTrees[0].nbfork);
							fprintf(stderr, "the number of forks is higher than the maximum allowed (%d)\n", forkmax);
							fprintf(stderr, "Epocs will generate random tvectors to explore + mandatory orders 1->2 and 2->1 at each co-occurrence\n");
						}
						MyTrees[f].MyCoevolData[t].theTVector.current_vector = malloc(sizeof(size_t)*MyTrees[f].MyCoevolData[t].nbranches);

						MyTrees[f].MyCoevolData[t].theTVector.tvector = genere_dual_vector_types(MyTrees[f].MyCoevolData[t].root, MyTrees[f].MyCoevolData[t].nbranches, MyTrees[f].MyCoevolData[t].theTVector.nbfork, (char **)NULL, &MyTrees[f].MyCoevolData[t].theTVector.nvect);

						add_random_vector_types_ok(MyTrees[f].MyCoevolData[t].theTVector.tvector, MyTrees[f].MyCoevolData[t].nbranches, 1000);

						MyTrees[f].MyCoevolData[t].theTVector.current_vector = MyTrees[0].MyCoevolData[t].theTVector.tvector;

						MyTrees[f].MyCoevolData[t].theTVector.nvect=1000;
						MyTrees[f].maxForkVector[t] = MyTrees[f].MyCoevolData[t].theTVector.nvect;
					}

				} // end of for (f=0; f<forestntree[t]; f++)
			} // end of for (t=0; t<ntree; t++)
		} // end of if (maxforest>1)

		/*
			VERBOSE - printing the transition vectors generated above
		*/

		if(verbose){
			for (f=0; f<maxforest; f++){
				for (t=0; t<ntree; t++) {
					printf("tree: %d ; #branches: %d ; #cooccurences: %d #trees+events+root: %d\n", maxforest>1?f:t,  MyTrees[f].MyCoevolData[t].nbranches, MyTrees[f].MyCoevolData[t].theTVector.nbfork, MyTrees[f].MyCoevolData[t].theTVector.nbfork>forkmax?(int)pow(2,forkmax)*3: (int)pow(2,MyTrees[f].MyCoevolData[t].theTVector.nbfork)*3);
					print_vector_types(MyTrees[f].MyCoevolData[t].theTVector.tvector, MyTrees[f].MyCoevolData[t].nbranches, stdout);
				}
			}
		}

		/*
			MAIN CALCULATION
			Iteration over the whole forest.
			If one entry -> calculations of the tree or over all independent trees
			If >1 entry (forest) -> calculations for each tree of the forest
			Every values are stored in an array of Trees (Hereafter named MyTrees)
		*/

		for (f=0; f<maxforest; f++){
			MyTrees[f].IS_epocs[0] = MyTrees[f].IS_epocs[1] = 0;									// Initial states
			memcpy((void *)MyTrees[f].final_rates_epocs, (void *)rate, (size_t)8*sizeof(double));	// initializing the rate matrix
			MyTrees[f].imax_epocs=0;																// Initializing the initial state value of e1 (updated if -I chosen)
			MyTrees[f].jmax_epocs=0;																// Initializing the initial state value of e2 (updated if -I chosen)
			MyTrees[f].o2=malloc(sizeof(int)*ntree);												// Allocating memory for the ids of the tvectors
			MyTrees[f].lik_epocs=0;																	// Initializing the maximum likelihood value

			/* Computing the maximum likelihood	*/
			ML_multi(MyTrees[f].MyCoevolData, ntree, &MyTrees[f].lik_epocs, rate, MyTrees[f].final_rates_epocs, MyTrees[f].IS_epocs, scenario, MyTrees[f].maxForkVector, 0, MyTrees[f].forkVector, MyTrees[f].o2, &MyTrees[f].imax_epocs, &MyTrees[f].jmax_epocs, opt_IS, 0);

			/*
				OUTPUT - writes the ML results to stdout and in the output file
			*/

			if (maxforest>1){
				fprintf(stdout, "TREE: %d\n", f);
			}

			printf("Mutation order: vector ids = ");
			for (t=0; t<ntree; t++)
				printf("%d, ", MyTrees[f].o2[t]);
			printf(";");

			for (t=0; t<ntree; t++){
			  for (i = 0; i<MyTrees[f].MyCoevolData[t].nbranches; i++){
				  if (MyTrees[f].MyCoevolData[t].theTVector.tvector[MyTrees[f].o2[t]][i] == 3) fprintf(stdout,"%d->%d ; ", evti+1, evtj+1);
				  if (MyTrees[f].MyCoevolData[t].theTVector.tvector[MyTrees[f].o2[t]][i] == 4) fprintf(stdout,"%d->%d ; ", evtj+1, evti+1);
			  }
			}
			printf("\n");

			for (i=0; i<2; i++){
				for (j=0; j<4; j++){
					if ((!(j%2)) && (MyTrees[f].final_rates_epocs[i][j] > maxrate)) {
						maxrate=MyTrees[f].final_rates_epocs[i][j];
					}
				}
			}

			fprintf(stdout,"IS={%d,%d} : ",MyTrees[f].imax_epocs,MyTrees[f].jmax_epocs);
			fprintf(stdout,"lnML = %.10f ML=%.5g\n",MyTrees[f].lik_epocs,exp(MyTrees[f].lik_epocs));

			fprintf(stdout,"%10c", ' ');
			fprintf(stdout,"%14s %14s %14s %14s\n","mu", "mu*", "nu", "nu*");

			/* Writing in the output file fptr	*/
			fprintf(fptr, "%-10d\t%-10d\t%-10d\t%-5c\t%-14.10f\t%-14.5g\t", evti+1, evtj+1, ntree>1?0:f+1, scenario, MyTrees[f].lik_epocs,exp(MyTrees[f].lik_epocs));

			for (i = 0; i < 2; i++){

			  fprintf(stdout,"event[%d]:",(i==0)?evti+1:evtj+1);

			  for (j = 0; j < 4; j++){

			  	if(MyTrees[f].final_rates_epocs[i][j]>=0)
			  	{
					if(MyTrees[f].final_rates_epocs[i][j] >= MAX_MU){
  					  fprintf(stdout, "%14.5g+",MyTrees[f].final_rates_epocs[i][j]>MAX_MU?MAX_MU:MyTrees[f].final_rates_epocs[i][j]);
  					  fprintf(fptr, "%-4.5g", MyTrees[f].final_rates_epocs[i][j]>MAX_MU?MAX_MU:MyTrees[f].final_rates_epocs[i][j]);
  					  if ((i+j) < 4)fprintf(fptr, "\t");
  					}
  					else{
  					  fprintf(stdout, "%15.5g", MyTrees[f].final_rates_epocs[i][j]);
  					  fprintf(fptr, "%-4.5g", MyTrees[f].final_rates_epocs[i][j]);
  					  if ((i+j) < 4)fprintf(fptr, "\t");
  				  	}
				  // if(MyTrees[f].final_rates_epocs[i][j] == MAX_MU){
					//   fprintf(stdout, "%14.5g+",MyTrees[f].final_rates_epocs[i][j]);
					//   fprintf(fptr, "%-4.5g", MyTrees[f].final_rates_epocs[i][j]);
					//   if ((i+j) < 4)fprintf(fptr, "\t");
					// }
					// else{
					//   fprintf(stdout, "%15.5g", MyTrees[f].final_rates_epocs[i][j]);
					//   fprintf(fptr, "%-4.5g", MyTrees[f].final_rates_epocs[i][j]);
					//   if ((i+j) < 4)fprintf(fptr, "\t");
				  // 	}
				}
				else{
				  fprintf(stdout, "%15c",'?');
				  fprintf(fptr, "%-4c", '-');
				  if ((i+j) < 4)fprintf(fptr, "\t");
			  	}
			  }
			  putchar('\n');
		  } // end of for (i = 0; i < 2; i++)

			printf("\n");
			fputs("\n", fptr);

			/*
				WARNING FUNCTION: Calculating with the final rates whether the probability of having >1 event on the longest branch is higher than 0.01
			*/
			if ((probaMoreThanOneEvent(MyTrees[f].MyCoevolData[0].longestbranch, maxrate)) > 0.01){fprintf(stderr, "WARNING: e1=%d e2=%d tree=%d proba of having >1 mutation in longest branch = %f\n\n", evti+1, evtj+1, f, probaMoreThanOneEvent(MyTrees[f].MyCoevolData[0].longestbranch, maxrate));}

			/* freeing tvectors for next run	*/
			maxrate=0;
			if (verbose)
				fprintf(stderr, "freeing vectors\n");

			for (t=0; t<ntree; t++){
				for (i=0; i<MyTrees[f].MyCoevolData[t].theTVector.nvect; i++){
					if (MyTrees[f].MyCoevolData[t].theTVector.tvector[i]) free(MyTrees[f].MyCoevolData[t].theTVector.tvector[i]);
				}
				if (MyTrees[f].MyCoevolData[t].theTVector.tvector) free(MyTrees[f].MyCoevolData[t].theTVector.tvector);
			}

		} // end of for (f=0; f<maxforest; f++)


		/*
			FOREST - if forest given in input, calculating the mode of the rates and the associated confidence intervals
		*/
		if (maxforest>1){
			fprintf(stdout, "Estimations and confidence intervals for the whole forest\n");
			double *lkDist=NULL;
			lkDist = ComputeLikelihoodDistribution(MyTrees, maxforest, 0, evti, evtj);
		}

	} /* end for while 1 */
	fclose(fptr);


	/*
		END OF THE MAIN LOOP
	*/

	/* Freeing memory	*/

	if (output_newick == 1){
		for (f=0; f<maxforest; f++){
			for (t=0; t<ntree; t++)
				PrintTree(MyTrees[f].MyCoevolData[t].root, NULL , 1);
		}
	}

	for (f=0; f<maxforest; f++){
		FreeCoevolData(MyTrees[f].MyCoevolData, ntree, 0, 0);
		if (MyTrees[f].o2)free(MyTrees[f].o2);
		if(MyTrees[f].nbfork<=forkmax){
			free(MyTrees[f].forkVector);
			free(MyTrees[f].maxForkVector);
		}
	}
	free(MyTrees);

	if (s) free(s);
	if (init_tree) FreeMultiTree(init_tree, nbleaves);
	if (input_events){
		freeArray(input_events, filelength);
		free(input_events_file);
	}

	return 0;
}
