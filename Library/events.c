/*
	Copyright (C) 2002-2018 G Achaz, J Pothier & S Brouillet

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
\file 	events.c
\brief 	basic event tree metrics and manipulations
\author g achaz, j pothier, madamesophie
\date 	May 23 2018
**/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tree.h"
#include "coevol.h"

// extern int verbose;

/* ----------------------------------------------------------------------------	*/
/* compte les evenements a partir du noeud n					*/
/* si un ou plusieurs noeuds n'ont pas de liste d'evenements, ils ne sont pas	*/
/* comptes, mais si deux noeuds different en nombre d'evenement, la fonction	*/
/* signale et sort avec erreur 2						*/
/* a n'appeler qu'(une fois (nevt est static !!!)				*/
/* ----------------------------------------------------------------------------	*/
int count_evt(Node *n)
{
	static int nevt=0;
	int i;

	if(n->anc == NULL)
		 nevt = 0;


	if (n->nevt != 0) {
		if (nevt != 0 && n->nevt != nevt) {
			fprintf(stderr,"No the same event number on node %s, exiting\n",n->name);
			exit(2);
		}
		nevt = n->nevt;
	}
	for (i = 0; i < n->nbdesc; i++) {
			count_evt( n->descs[i]);
	}

	return nevt;
}

/* ----------------------------------------------------------------------------	*/
/* verifie que chaque noeud a partir du noeud n a bien				*/
/* une liste de nevt evenements, sinon						*/
/* alloue cette liste, et initialise a 0					*/
/* renvoie le nombre de noeud ou une liste a ete initialisee			*/
/* ----------------------------------------------------------------------------	*/
int verif_evt(Node *n, int nevt, int verbose)
{
	int changed=0, i;
	static int branchnum = -1;

	if(n->anc == NULL)
		 branchnum = -1;

	branchnum++;

	if (n->nevt == 0) {
	        n->nevt=nevt;
		n->evt=malloc((n->nevt)*sizeof(int));
		if (n->evt==NULL) {
			fprintf(stderr,"Malloc error for event list (node %s),exiting\n",n->name);
			exit(11);
		}
		memset((void *)n->evt, 0, (size_t)n->nevt*sizeof(int));
		if (verbose) fprintf(stderr,"Node %s [%d] has been given a list of %d events (filled with 0)\n",n->name, branchnum, n->nevt);
		changed = 1;
	}
	else if (n->nevt != nevt) {
		fprintf(stderr,"Should not happen: %s [%d] node has not %d events in its list, but %d events ! exiting\n",n->name, branchnum, nevt, n->nevt);
		exit(12);
	}

	/* va voir chez les descendants	*/
	for (i = 0; i < n->nbdesc; i++) {
			changed += verif_evt(n->descs[i], nevt, verbose);
	}


	return changed;
}

/* ----------------------------------------------------------------------------	*/
/* compte les evenements a chaque noeud et renvoie le plus grand nombre		*/
/* d'element trouve a un noeud							*/
/* ----------------------------------------------------------------------------	*/
int count_max_event(Node *n, int nevt, int verbose)
{
	int i, maxnevt;
	static int branchnum = -1;

	if(n->anc == NULL)
		 branchnum = -1;


	branchnum++;

	if (n->nevt > nevt)
		maxnevt = n->nevt;
	else
		maxnevt = nevt;

	if (verbose)
		fprintf(stderr, "max evt is %d at node %d (%d events)\n",maxnevt,branchnum,n->nevt);

	/* va voir chez les descendants	*/
	for (i = 0; i < n->nbdesc; i++) {
			maxnevt = count_max_event(n->descs[i], maxnevt, verbose);
	}

	return maxnevt;
}

/* ----------------------------------------------------------------------------	*/
/* verifie que chaque noeud a partir du noeud n a bien				*/
/* une liste de nevt evenements composee de 0 ou 1				*/
/* si pas de liste, alloue cette liste, et initialise a 0 les nevt evenements	*/
/* si une liste avec un nombre different d'evts, EXIT				*/
/* renvoie le nombre de noeud ou une liste a ete initialisee			*/
/* ----------------------------------------------------------------------------	*/
int verif_nevt(Node *n, int nevt, int verbose)
{
	int changed=0, i;
	static int branchnum = -1;

	if(n->anc == NULL)
		 branchnum = -1;


	branchnum++;


	if (n->nevt == 0) {
	        n->nevt=nevt;
		n->evt=malloc((n->nevt)*sizeof(int));
		if (n->evt==NULL) {
			fprintf(stderr,"Malloc error for event list (node %s),exiting\n",n->name);
			exit(11);
		}
		memset((void *)n->evt, 0, (size_t)n->nevt*sizeof(int));
		if (verbose) fprintf(stderr,"Node %s [%d] has been given a list of %d events (filled with 0)\n",n->name, branchnum, n->nevt);
		changed = 1;
	}
	else if (n->nevt != nevt) {
		fprintf(stderr,"Should not happen: %s [%d] node has not %d events in its list, but %d events:\n",n->name, branchnum, nevt, n->nevt);
		for (i = 0; i < n->nevt; i++) {
		   fprintf(stderr,"event %d -> %d\n",i,n->evt[i]);
		}
		fprintf(stderr,"cannot continue, exiting...\n");
		exit(12);
	}

	/* verifie les valeurs 0 ou 1 */
	for (i = 0; i < n->nevt; i++) {
	   if (n->evt[i] != 0 && n->evt[i] != 1) {
	      fprintf(stderr,"Error: Node %d name %s event %d -> %d is not 0 nor 1\n",branchnum, n->name, i,n->evt[i]);
	      exit(15);
	   }
	}

	/* va voir chez les descendants	*/
	for (i = 0; i < n->nbdesc; i++) {
		changed += verif_nevt(n->descs[i], nevt, verbose);
	}

	return changed;

}

/* -------------------------------------------------------------------------------- */
/* verify that each node from node n consists of a list of nevt events				*/
/* exclusively composed of 0, 1 or -1 (unknown trait)								*/
/* if no list is present, allocatre this list and initialize with 0 the nevt events	*/
/* If a list exists with different number of events, EXIT							*/
/* return the number of nodes needed to be initialized								*/
/* --------------------------------------------------------------------------------	*/
int verif_nevt_gaps(Node *n, int nevt, int verbose)
{
	int changed=0, i;
	static int branchnum = -1;

	if(n->anc == NULL)
		 branchnum = -1;


	branchnum++;


	if (n->nevt == 0) {
	        n->nevt=nevt;
		n->evt=malloc((n->nevt)*sizeof(int));
		if (n->evt==NULL) {
			fprintf(stderr,"Malloc error for event list (node %s),exiting\n",n->name);
			exit(11);
		}
		memset((void *)n->evt, 0, (size_t)n->nevt*sizeof(int));
		if (verbose) fprintf(stderr,"Node %s [%d] has been given a list of %d events (filled with 0)\n",n->name, branchnum, n->nevt);
		changed = 1;
	}
	else if (n->nevt != nevt) {
		fprintf(stderr,"Should not happen: %s [%d] node has not %d events in its list, but %d events:\n",n->name, branchnum, nevt, n->nevt);
		for (i = 0; i < n->nevt; i++) {
		   fprintf(stderr,"event %d -> %d\n",i,n->evt[i]);
		}
		fprintf(stderr,"cannot continue, exiting...\n");
		exit(12);
	}

	/* verifie les valeurs 0 ou 1 */
	for (i = 0; i < n->nevt; i++) {
	   if (n->evt[i] != 0 && n->evt[i] != 1 && n->evt[i] != -1) {
	      fprintf(stderr,"Error: Node %d name %s event %d -> %d is not 0 nor 1\n",branchnum, n->name, i,n->evt[i]);
	      exit(15);
	   }
	}

	/* va voir chez les descendants	*/
	for (i = 0; i < n->nbdesc; i++) {
		changed += verif_nevt_gaps(n->descs[i], nevt, verbose);
	}

	return changed;

}

/* ----------------------------------------------------------------------------	*/
/* met le type d'evt a chaque noeud selon les evt i et j    			*/
/* retourne le nombre de fork							*/
/* ----------------------------------------------------------------------------	*/
int set_evt_type_count_fork(Node *n, int nevt, int evti, int evtj, int verbose, int output_state_vector)
{
	int nbfork=0, i;
	static int branchnum = -1;

	if(n->anc == NULL)
		 branchnum = -1;


	branchnum++;

	/* DEBUG fprintf(stderr,"Noeud %d evt[%d] = %d, evt[%d] = %d\n",n->id,evti,n->evt[evti],evtj,n->evt[evtj]); */
	if (n->evt[evti] == n->evt[evtj]) {
		if (n->evt[evti]==0)
			n->evtype = (char)0; /* pas d'evenement */
		else {
			n->evtype = (char)3; /* 2 evenements (sans ordre, donnera 3: ordre 12 ou 4: ordre 21)  */
			nbfork = 1;
		}
	}
	else {
		if (n->evt[evti]==1)
			n->evtype = (char)1; /* evenement 1 */
		else
			n->evtype = (char)2; /* evenement 2 */
	}

	if (output_state_vector) {
		fprintf(stderr, "Node %d: type event: %d\n",branchnum,n->evtype);
	}

	/* va voir chez les descendants	*/
	for (i = 0; i < n->nbdesc; i++) {
		nbfork += set_evt_type_count_fork(n->descs[i], nevt, evti, evtj, verbose, output_state_vector);
	}

	return nbfork;

}

/* ----------------------------------------------------------------------------	*/
/* met le type d'evt a chaque noeud selon les evt i et j    			*/
/* retourne le nombre de fork							*/
/* ----------------------------------------------------------------------------	*/
int set_evt_type_count_fork_gaps(Node *n, int nevt, int evti, int evtj, int verbose, int output_state_vector)
{
	int nbfork=0, i;
	static int branchnum = -1;

	if(n->anc == NULL)
		 branchnum = -1;


	branchnum++;

	/* DEBUG fprintf(stderr,"Noeud %d evt[%d] = %d, evt[%d] = %d\n",n->id,evti,n->evt[evti],evtj,n->evt[evtj]); */
	if (n->evt[evti] == -1 || n->evt[evtj] == -1){
		n->evtype=(char)5;
	}else{
		if (n->evt[evti] == n->evt[evtj]) {
			if (n->evt[evti]==0)
				n->evtype = (char)0; /* pas d'evenement */
			else {
				n->evtype = (char)3; /* 2 evenements (sans ordre, donnera 3: ordre 12 ou 4: ordre 21)  */
				nbfork = 1;
			}
		}
		else {
			if (n->evt[evti]==1)
				n->evtype = (char)1; /* evenement 1 */
			else
				n->evtype = (char)2; /* evenement 2 */
		}
	}

	if (output_state_vector) {
		fprintf(stderr, "Node %d: type event: %d\n",branchnum,n->evtype);
	}

	/* va voir chez les descendants	*/
	for (i = 0; i < n->nbdesc; i++) {
		nbfork += set_evt_type_count_fork_gaps(n->descs[i], nevt, evti, evtj, verbose, output_state_vector);
	}

	return nbfork;

}

/* --------------------------------------------------------------------------- */
/* check at all nodes below n have a different state than n                    */
/* n is root at the beginning and root has [0,0]                               */
/* thus counts when sites have # events > 0                                    */
/* --------------------------------------------------------------------------- */
int verif_polymorphism(Node *n, int *v_evt, int nevt)
{
	int polymorph=0, i;
	static int branchnum = -1;

	if(n->anc == NULL)
		 branchnum = -1;


	/* At the root, set v_ets to be the type of events from the root */
	if (branchnum == -1) {
	   for (i = 0; i < nevt; i++)
	      v_evt[i] = n->evt[i];
	   branchnum++;
	}
	else {
		branchnum++;
		for (i = 0; i < nevt; i++) {
		  if (v_evt[i] != -1 && n->evt[i] != v_evt[i]) {
		     v_evt[i] = -1;
		     polymorph++;
		  }
		}
	}

	/* then check the descendants */
	for (i = 0; i < n->nbdesc; i++) {
			polymorph += verif_polymorphism(n->descs[i], v_evt, nevt);
	}

	return polymorph;
}

/* ----------------------------------------------------------------------------	*/
/* imprime un vecteur de type d'evenement					*/
/* ----------------------------------------------------------------------------	*/
void print_vector_types(char **ts, int len, FILE *fout) {
	int i;
	while (*ts) {
		for (i = 0; i < len; i++)
			fprintf(fout,"%1d",(int)(*ts)[i]);
		fprintf(fout,"\n");
		ts++;
	}
	return;
}

/* ------------------------------------------------------------------ */
/* calcule les vecteurs d'evtype                                      */
/* actual_evt_type: vecteur des types d'evenement jusque ici          */
/* ------------------------------------------------------------------ */
char **genere_all_vector_types(Node *n, int last_branch, int nbfork , char **ts, int *nvect)
{
	int i;
	char *ns,    /* pointer to an event from the event sequence */
	     **ps,   /* pointer to a particular order */
	     **ps2,  /* */
	     **fin;

	if(nbfork > 20) {
		fprintf(stderr, "cannot handle more than 20 branches with double events, sorry. bye\n"),exit(3);
	}
	/* at the root, nothing to compute, just call childs  */
	/* DEBUG: fprintf(stderr, "node id=%d ev=%d last branch=%d\n",n->id,n->evtype, last_branch); */

	if (n->id == 0) {
		/* alloue un tableau de pointeur de chaine de 2^last_branch + 1 pour la deerniere case a zero (pour signal stop)	*/
		*nvect = (2L<<nbfork)+1;
		if ((ts = (char **)calloc((size_t) *nvect, sizeof(char *))) == NULL) {
			fprintf(stderr,"genere_all_vector_types: Not enough memory\n");exit(3);}
		/* alloue la premiere chaine	*/
		if ((ns = (char *)calloc((size_t)last_branch, sizeof(char))) == NULL) {
			fprintf(stderr,"genere_all_vector_types: Not enough memory\n");exit(3);}
		ts[0]=ns;
	}
	else {
		/* for ev=0,1,2,3 i.e. 00, 10, 02, 12 */
		/* met a jour les chaines d'evenements iavec l'evenement courant */
		ps = ts;
		while (*ps) {
			(*ps)[n->id-1]=n->evtype;
			ps++;
		}
		/* si l'evenement courant est 3 (1->2) il faut dupliquer toutes les chaines 	*/
		/* jsuqu'ici, et mettre l'evenement 4 (2->1) */
		if (n->evtype == 3) {
			/* duplique toutes les chaines d'eveneemnts 	*/
			fin = ps-1;		/* dernier pointeur de la liste courante	*/
			ps2 = ps;		/* sera le pointeur vers la premiere case libre	*/
			ps = ts;		/* ps reprend la premiere case			*/
			while (ps <= fin) {

				//printf("copy from site %d\n", n->id);

				/* duplique */
				//fprintf(stderr, "genere_all_vector_types::allocate %d for n->id=%d\n",last_branch,n->id);

				if ((ns = (char *)calloc((size_t)last_branch, sizeof(char))) == NULL) {fprintf(stderr,"genere_all_vector_types: Not enough memory\n");exit(3);}

				memcpy((void *)ns, (const void *)(*ps), (size_t) (n->id-1) * sizeof(char));
				ns[n->id-1]=4;	/* modifie le dernier type d'evenement (i.e. met 4)	*/
				(*ps2)=ns;	/* ajoute la chaine a la premiere position libre	*/
				ps++;
				ps2++;
			}
		}
	}

	/* appelle tous les fils */
	for ( i = 0; i < n->nbdesc; i++) {
		genere_all_vector_types(n->descs[i], last_branch, nbfork, ts, nvect);
	}

	return ts;
}


/* ------------------------------------------------------------------ */
/* generation of tvectors 1->2 and 2->1 possibilities at forks        */
/* ------------------------------------------------------------------ */

char **genere_dual_vector_types(Node *n, int last_branch, int nbfork , char **ts, int *nvect)
{
	int i,j;
	char *ns, *ns2;

	/* at the root, nothing to compute, just call childs  */
	/* DEBUG: fprintf(stderr, "node id=%d ev=%d last branch=%d\n",n->id,n->evtype, last_branch); */

	if (n->id == 0) {
		/* allocate an array of 1000 (forward, backward, 998 random) + 1 (stop signal)	*/
		*nvect = 1001;

		if ((ts = (char **)calloc((size_t) *nvect, sizeof(char *))) == NULL) {
			fprintf(stderr,"genere_dual_vector_types: Not enough memory\n");exit(3);}

		/* allocating the two arrays directly	*/
		if ((ns = (char *)calloc((size_t)last_branch, sizeof(char))) == NULL) {
			fprintf(stderr,"genere_dual_vector_types: Not enough memory\n");exit(3);}

		if ((ns2 = (char *)calloc((size_t)last_branch, sizeof(char))) == NULL) {
			fprintf(stderr,"genere_dual_vector_types: Not enough memory\n");exit(3);}
		ts[0]=ns;
		ts[1]=ns2;
	}
	else {

		for (j=0; j<2; j++){
			ts[j][n->id-1]=n->evtype;
		}
		if (n->evtype == 3){
			ts[1][n->id-1]=4;
		}
	}

	/* call the descendants */
	for ( i = 0; i < n->nbdesc; i++) {
		genere_dual_vector_types(n->descs[i], last_branch, nbfork, ts, nvect);
	}

	return ts;
}

void add_random_vector_types_ok(char **ts, int last_branch, int maximumTvectors){
	int i,k;
	for (k=2; k<maximumTvectors; k++){														// I want to allocated maximum 1k vectors, starting from 2 ( 0 is 1->2 everywhere, 1 is 2->1 everywhere)
		if ((ts[k] = (char *)calloc((size_t)last_branch, sizeof(char))) == NULL) {			// I need of course to allocate memory (**ts already has the space for 1k sequences)
			fprintf(stderr,"add_random_vector_types: Not enough memory\n");exit(3);}

		memcpy(ts[k], ts[0], last_branch*sizeof(char));
		for (i=0; i<last_branch; i++){
			if (ts[k][i]==3){
				if ((rand()%100)>=50)
					ts[k][i]=4;
			}
		}
	}
}


/* ----------------------------------------------------------------------------	*/
/* change les noms de feuilles et de noeud de maniere a avoir les evenements	*/
/* polymorphes inscrits dans les noms de feuilles et de noeud internes		*/
/* ----------------------------------------------------------------------------	*/
void change_names_with_events(Node *n, int *v_evt, int nevt)
{
	int l=0, i;
#define LSTRTMP 1024
	static int branchnum = -1;
	static char tmpstr[LSTRTMP], s[LSTRTMP], ss[15], *t;

	if(n->anc == NULL)
		 branchnum = -1;


	branchnum++;

	/* copy the old name	*/
        tmpstr[0]='\0';
	l = LSTRTMP;
	t = n->name;
	/* suppress the events from old name */
	while(*t) {
	  if (*t=='[')
	     *t = '\0';
	  else
	    t++;
	}
	strncat(tmpstr, n->name, l);
	l -= strlen(n->name);
	strncat(tmpstr, " ", l);
	l-=1;

	/* construct the events string 	*/
        s[0]='\0';
        ss[0]='\0';
	strncat(ss,"[",2);
	l-=1;
	for (i = 0; i < nevt; i++) {
	   if (v_evt[i] == -2) {
	     strncat(s,ss,l);
	     l-=sprintf(ss,"%d,",n->evt[i]);
	   }
	}
	strncat(s,ss,l);
	strncat(s,"]",l);
	l-=1;

	/* append the events string 	*/
	if (l < 1) {
	   fprintf(stderr,"Name too long with event (function: change_names_with_events, node: %s)\n",n->name);
	   exit(0);
	}
	strcat(tmpstr, s);

	/* copy the new name 	*/
        l = LSTRTMP - l + 1; 	/* longueur du nom maintenant avec le \0 final */

	if ((n->name=realloc(n->name, sizeof(char)*(size_t)l))==NULL) {
		fprintf(stderr,"Malloc error for node name (branch %d)\n",branchnum);
		exit(15);
	}
        n->name[0]='\0';
	strncat(n->name, tmpstr,l+1);

	/* va voir chez les descendants	*/
	for (i = 0; i < n->nbdesc; i++) {
		change_names_with_events(n->descs[i], v_evt, nevt);
	}

	return;
}

/* --------------------------------------------------	*/
/* set_e1_vectors: recupere les vecteurs d'evenements	*/
/* dans l'arbre						*/
/* --------------------------------------------------	*/

void set_e1_vectors(Node *n, int **e1)
{
   static int numBranche=-2 ;
   int i;

	if(n->anc == NULL)
		 numBranche = -2;


   numBranche++;
   if (numBranche != -1) { 	/* sur le root, on ne fait rien (pas de branche)	*/

#if 0
	   fprintf(stderr,"set_e1_vectors:: Branch: %d (node %s)\n",numBranche,n->name);
#endif

	   for (i = 0; i < n->nevt; i++) {
		   e1[i][numBranche]=n->evt[i];
	   }
   }

   if (n->nbdesc > 0) {
	   for (i = 0; i < n->nbdesc; i++) {
		   set_e1_vectors(n->descs[i], e1);
	   }
   }

 return;

}

/* --------------------------------------------------			*/
/* set_e1_vectors_mask: recupere les vecteurs d'evenements		*/
/* dans l'arbre and sets a flag if it encounters an unknown char*/
/* --------------------------------------------------			*/

void set_e1_vectors_mask(Node *n, int **e1, char **mask)
{
   static int numBranche=-2 ;
   int i;

	if(n->anc == NULL)
		 numBranche = -2;


   numBranche++;
   if (numBranche != -1) { 	/* sur le root, on ne fait rien (pas de branche)	*/

#if 0
	   fprintf(stderr,"set_e1_vectors_mask:: Branch: %d (node %s)\n",numBranche,n->name);
#endif

	   for (i = 0; i < n->nevt; i++) {
		   e1[i][numBranche]=n->evt[i];

		   if ((int)n->evt[i] >= 0){
		   	   mask[i][numBranche]=1;
			}
	   }
   }

   if (n->nbdesc > 0) {
	   for (i = 0; i < n->nbdesc; i++) {
		   set_e1_vectors_mask(n->descs[i], e1, mask);
	   }
   }

 return;

}
