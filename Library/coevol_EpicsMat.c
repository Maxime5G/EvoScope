/* ------------------------------------------------------------------------	*/
/* Calcul des matrices de chronologies/cooccurences et des multinomiales	*/
/* ------------------------------------------------------------------------	*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>	/* to terminate the multinomial if it takes too long */

#ifdef LINUX
#include <strings.h>
#endif
#include <string.h>

#include <math.h>
#include "tree.h"
#include "coevol.h"

#define MAXTIME 180 /* to terminate the multinomial if it takes more than half an hour (1800) */

// extern int verbose;

/*-------------------------------------------------------*/
/* statiques						*/
/*-------------------------------------------------------*/
static double *sFacto;
static double *slnFacto;

static int **sTableau=NULL,
           **sTableau2=NULL;			/* de longueur nombre de branches */

/*-------------------------------------------------------*/
/* PART MATRICES / CHRONOLOGIES			*/
/*-------------------------------------------------------*/

/*-------------------------------------------------------	*/
/* -- routine called by qsort -----------------------------	*/
/* returns 1 if i>j, -1 if i<j and 0 if i==j			*/
/*-------------------------------------------------------	*/
int estplusgrand(const void *i, const void *j)
{
	if (*((int *)i) == *((int *)j))
		return 0;
	else if (*((int *)i) > *((int *)j))
		return 1;
	else
		return -1;
}

/*-------------------------------------------------------	*/
/* -- routine called by qsort -----------------------------	*/
/* retourne 1 if i<j, -1 if i>j and 0 if i==j			*/
/*-------------------------------------------------------	*/
int estpluspetit(const void *i, const void *j)
{
	if (*((int *)i) == *((int *)j))
		return 0;
	else if (*((int *)i) < *((int *)j))
		return 1;
	else
		return -1;
}

/*-------------------------------------------------------	*/
/* returns the sum of the integer vector v of length lv	    */
/*-------------------------------------------------------	*/
int sumv(int *v, int lv)
{
   int i, sum = 0;

   for ( i = 0; i < lv; i++)
     sum+=v[i];

   return sum;

}

/*-------------------------------------------------------	*/
/* returns the number of nnz elements of the integer        */
/* vector v of length lv                                    */
/*-------------------------------------------------------	*/

int nnzv(int *v, int lv)
{
    int i, sum = 0;

   for ( i = 0; i < lv; i++){
       if (v[i]!=0)
            sum+=v[i];
    }

    return sum;

}

/*-------------------------------------------------------	*/
/* returns the number of nnz elements of the integer        */
/* vector v of length lv                                    */
/*-------------------------------------------------------	*/

int sumvP(int *v, int lv)
{
    int i, sum = 0;

   for ( i = 0; i < lv; i++){
       if (v[i]>0)
            sum+=v[i];
    }

    return sum;

}

/* ------------------------------------------------------------ */
/* prints an integer vector v of length l in file f			    */
/* ------------------------------------------------------------ */
void print_ivect(FILE *f, int *v, int l)
{
   int i;

   for (i = 0; i < l; i++)
      fprintf(f,"%d,",v[i]);

   fprintf(f,"\n");

   return;

}


/* ------------------------------------------------------------ */
/* product of the integer vector v by the integer matrix m	    */
/* the resulting product is in pv							    */
/* The vector length should be equal to the number of rows of   */
/* the matrix m. 												*/
/* ------------------------------------------------------------ */
void vectmat(int *v, IntegerMatrix *m, int *pv)
{
   int	i,j,k;

   for (j = 0; j < m->nrow; j++) {
   	pv[j] = 0;
   	for (i = 0, k=j; i < m->nrow; i++,k+=m->nrow) {
	  pv[j] += v[i]*m->val[k];
	}
   }

   return;
}

/* ------------------------------------------------------------ */
/* scalar product of u by v, vectors of length l			    */
/* ------------------------------------------------------------ */
int iprodsca(int *u, int *v, int l)
{
   int ps = 0, i;

   for (i = 0; i < l; i++)
     ps += u[i]*v[i];

   return ps;
}

/* -------------------------------------------------------------- */
/* Initialize an empty integer matrix of nrow rows and ncol cols. */
/* The matrix object contains the nrows, ncols and the values in  */
/* vector val.													  */
/* -------------------------------------------------------------- */

IntegerMatrix *InitIntegerMatrix(int nrow, int ncol)
{
  IntegerMatrix *m;

  m =malloc(sizeof(IntegerMatrix));
  if (m==NULL) {
    fprintf(stderr,"Not enough memory for matrix allocation, exiting\n");
    exit(1);
  }
  m->ncol=ncol;
  m->nrow=nrow;
  m->val =malloc(sizeof(int)*m->nrow*m->ncol);
  if (m->val==NULL) {
    fprintf(stderr,"Not enough memory for matrix values allocation, exiting\n");
    exit(1);
  }

  return m;
}

/* -------------------------------------------------------------- */
/* Free up the memory of any integer matrix object				  */
/* -------------------------------------------------------------- */

void FreeIntegerMatrix(IntegerMatrix *m)
{
   if (m->val)
   	free(m->val);
   if (m)
   	free(m);
   return;
}

/* -------------------------------------------------------------- */
/* Fills the event matrix and the branch length vector.			  */
/* Returns the total length of branches.						  */
/* -------------------------------------------------------------- */

double fill_evt_matrix( Node *n, IntegerMatrix *matrix_to_fill, double *lengths, int nbranches)
{
  int i,k;
  static int branche_num=-2;
  static double total_length = 0.0L;

  if( n->anc == NULL ){
	branche_num=-2;
	total_length = 0.0L;
  }


#if 0
  printf("fill_evt_matrix:: Branch: %d (Node %s)\n",branche_num,n->name);
#endif

  /* if (! n) {fprintf(stderr,">>>>>> fill_evt_matrix:: node is NULL NULL NULL\n"); exit(0);} else {fprintf(stderr,"fill_evt_matrix:: node is not null\n");} */

  branche_num++;
  if (branche_num != -1) {	/* root node - doesn't count as a branch. */
      /* We set the length of the branch in the lengths array   */
	  if (branche_num >= nbranches) {
	     fprintf(stderr,"fill_evt_matrix::branch number %d is greater than nbranches (%d)\n",branche_num+1,nbranches);
	     exit(1);
	  }
	  lengths[branche_num]= n->time;
	  /* fprintf(stderr,"longueur[%d]=%g\n",branche_num,n->time); */
	  total_length += n->time;

	  /* line: nbranches, column: j from 0 to nb of events */

	  for (i=0;i<(matrix_to_fill->nrow);i++) {
		  k = (matrix_to_fill->ncol*i)+branche_num ;
		  matrix_to_fill->val[k]=n->evt[i];
	  }

  }
  if (n->nbdesc != 0) {
	  for (i = 0; i < n->nbdesc; i++) {
		  fill_evt_matrix ( n->descs[i],matrix_to_fill,lengths,nbranches);
	  }
  }

  return total_length;

}


/* -------------------------------------------------------------- */
/* Fills the chronological matrix (or S matrix). 				  */
/* I.e. i,j=1 in the matrix if j follows i, 0 else				  */
/* The idea is to set the current node at 1 in vector t.		  */
/* When going down the recursion following the descendants, their */
/* ascendants have the value 1 in the vector. For each node, we   */
/* fill the line of 0's until the node position, and afterwards   */
/* we fill the value of the ascendants (which is 1).              */
/* We reset the current vector position to zero after the         */
/* Since he won't be an ascendant for the following nodes.        */
/* -------------------------------------------------------------- */


void chrono(Node *n, int *t, IntegerMatrix *chronomat)
{
  int i,j,k, mabranche = 0;
  static int branche_num = -2, recurs = -1;

	if(n->anc == NULL)
	{
		branche_num = -2;
	}

  branche_num++;
  if (branche_num != -1) { /* sur la racine, on ne fait rien	*/
	  recurs++;
	  mabranche = branche_num;

	  /* remplissage de la ligne branche_num avec des 0 de	*/
	  /* 0 a la colonne branche_num,  dans la matrice	*/
	  /* note: pour une position i,j k = ncol*i+j		*/
	  for (j = 0, k=chronomat->ncol*branche_num; j <= branche_num; j++,k++) {
		  chronomat->val[k]=0;
	  }


	  /* remplissage de la colonne branch_num dans la matrice	*/
	  /* de la ligne 0 a la ligne branche_num (ligne du noeud	*/
	  /* courant							*/
	  for (i=0,k=branche_num; i<branche_num; i++,k+=chronomat->ncol) {
		  chronomat->val[k]=t[i];
	  }

	  /* on met t a un pour les descendants */
	  t[mabranche]=1;
	  /*
	     for (i=0;i<recurs;i++)
	     printf("     ");
	     printf("->[%3d] %s\n",mabranche,n->name);
	   */
  }
  if (n->nbdesc!= 0) {
	  for (i = 0; i < n->nbdesc; i++) {
		  chrono(n->descs[i],t,chronomat);
	  }
  }

  /* quand on sort, on remet t a 0 pour les noeuds qui seront	*/
  /* a cote (donc pas descendants)				*/
  if (branche_num != -1) {
	  t[mabranche]=0;
	  recurs--;
  }

  /*
  for (i=0;i<mabranche;i++)
    printf(" ");
  printf("<-[%3d] %s\n",mabranche,n->name);
  */
  return;

}

/*----------------------------------------------------	*/
/* Fills in the Adjacency matrix (initialized with 0's  */
/* Everywhere).											*/
/* The idea is to go through the tree recursively,      */
/* similarly than chrono, and fill the matrix of 		*/
/* adjacency by going through the number of columns     */
/* of the chronomat matrix (i.e. #branches) and adding  */
/* zeroes. When the index crosses a descendant, we add  */
/* a one.  												*/
/*----------------------------------------------------	*/

void chrono_adj(Node *n, IntegerMatrix *chronomat){

    int i,j,k,l=0,ind,mybranch=0;
    static int branch_num = -2;

    if (n->anc == NULL){
        branch_num = -2;
    }

    branch_num++;
    if (branch_num != -1){
		mybranch=branch_num;
		if (n->nbdesc != 0){
            // If there are indeed descendants
			for (j=0, k=n->descs[l]->id; j < chronomat->ncol; j++){

				// I'm going through each branch id in the chronomat
				// the k is the index of the first descendant

				ind=chronomat->ncol*mybranch+j;

				if (j+1==k){ 							// If the position in the matrix corresponds to the descendant
					chronomat->val[ind] = 1; 			// The value in the matrix is 1
					l++;
					if (l < n->nbdesc) 					// If I still have a descendant
						k=n->descs[l]->id; 				// I'm incrementing the value of k to be the next descendant.
				}
				else
					chronomat->val[ind] = 0; 			// If the position in the matrix is not the descendant, the value in the matrix is 0
			}
		}
		else{ 											// There are no descendants (i.e. we have a leaf)
			for (j=0; j < chronomat->ncol; j++){
				ind=chronomat->ncol*mybranch+j;
				chronomat->val[ind] = 0;
			}
		}
    }
    if (n->nbdesc != 0){ 								// And going through the tree recursively
        for (i=0; i<n->nbdesc; i++){
            chrono_adj(n->descs[i], chronomat);
        }
    }
    return;
}


/*----------------------------------------------------	*/
/* Function to transpose a matrix. Required to generate */
/* the precedence matrix, since the precedence matrix   */
/* is the transpose of the S matrix.   					*/
/*----------------------------------------------------	*/

void transpose_a_matrix(IntegerMatrix *mat){

    int i,j,k,l;
    int * tmp = malloc(sizeof(int)*mat->nrow*mat->ncol);

    for (i=0; i<mat->nrow; i++){
        for (j=0; j<mat->ncol; j++){
                k=(i*mat->ncol)+j;
                l=(j*mat->nrow)+i;

                tmp[k]=mat->val[l];

        }
    }

    for (i=0; i<mat->nrow*mat->ncol; i++)
        mat->val[i] = tmp[i];

    free(tmp);
    return;

}

/*----------------------------------------------------	*/
/* Function required to generate the exclusion matrix   */
/* Basically, the exclusion matrix is the sum of the S  */
/* matrix and the P matrix, then flipped on the anti-   */
/* diagonal.											*/
/*----------------------------------------------------	*/

void transpose_a_matrix_antidiagonal(IntegerMatrix *mat){

    int i,j,k,l;
    int * tmp = malloc(sizeof(int)*mat->nrow*mat->ncol);

    for (i=0; i<mat->nrow; i++){
        for (j=0; j<mat->ncol; j++){
                k=(i*mat->ncol)+j;
                l = (mat->ncol*mat->ncol) - (j*mat->ncol) - i - 1;

                tmp[k]=mat->val[l];
        }
    }

    for (i=0; i<mat->nrow*mat->ncol; i++)
        mat->val[i] = tmp[i];

    free(tmp);
    return;

}

/*----------------------------------------------------	*/
/* Function to generate the exclusion matrix   			*/
/* The idea is to generate both S and P matrices, sum   */
/* them together, and flip the resulting matrix along   */
/* the anti-diagonal. That way, we have a matrix where  */
/* only branches that are on different clades of the    */
/* tree can have the value 1. 							*/
/*----------------------------------------------------	*/

IntegerMatrix * chrono_excl(Node *n, int nbranches){

        int i,j,k;
        int *chrono_vect1 = NULL;
        int *chrono_vect2 = NULL;
        IntegerMatrix *tmp1=NULL;
        IntegerMatrix *tmp2=NULL;
        IntegerMatrix *tmp3=NULL;

        if ((chrono_vect1 = calloc(nbranches,sizeof(int))) == NULL)
                fprintf(stderr, "memory error for time vectors, exit\n"),	exit(1);

        if ((chrono_vect2 = calloc(nbranches,sizeof(int))) == NULL)
                fprintf(stderr, "memory error for time vectors, exit\n"),	exit(1);

        tmp1 = InitIntegerMatrix(nbranches,nbranches);
        tmp2 = InitIntegerMatrix(nbranches,nbranches);
        tmp3 = InitIntegerMatrix(nbranches,nbranches);

        chrono(n,chrono_vect1, tmp1);
        chrono(n,chrono_vect2, tmp2);
        transpose_a_matrix(tmp2);

        for (i=0; i<tmp3->ncol; i++){
            for (j=0; j<tmp3->nrow; j++){
                k = i+j*tmp3->nrow;
                tmp3->val[k] = tmp1->val[k]+tmp2->val[k];
            }
        }

        free(chrono_vect1);
        free(chrono_vect2);
        FreeIntegerMatrix(tmp1);
        FreeIntegerMatrix(tmp2);
        transpose_a_matrix_antidiagonal(tmp3);

        return tmp3;
}


/*-------------------------------------------------------*/
/* VERBOSE FUNCTION - print the chronology matrix chosen */

void print_matrice(IntegerMatrix *m, FILE *f)
{  int i,j;

  for (i=0;i<m->nrow;i++){
    for (j=0;j<m->ncol;j++){
      fprintf(f,"%1d ",m->val[(i*m->ncol)+j]);
    }
   fprintf(f,"\n");
  }
  return ;
}

/*-------------------------------------------------------*/
/* PART 2: MULTINOMIAL calculations            			*/
/*-------------------------------------------------------*/

/* ------------------------------------------------------------ */
/* gives the next permutation of vector v, such that the value  */
/* is ascending, e.g., 123 -> 132 -> 213 -> 231 -> 312 -> 321   */
/* allows to obtain permutations without repetition             */
/* e.g. 113 -> 131 -> 311                                       */
/* we start with an array sorted in ascending order and we end  */
/* with a vector sorted in descending order                     */
/* ------------------------------------------------------------ */

int get_next_permut(int *v, int l)
{
        int tail, s, i, j, t, encore = 0;

        for (tail = l - 1; tail > 0; tail--) {
                if (v[tail - 1] < v[tail]) {/*still increasing */
			encore = 1;
                        /*find last element which does not exceed ind[tail-1] */
                        s = l - 1;
                        while(v[tail-1] >= v[s])
                                s--;

                        /* swap(v, tail-1, s); */
                        t=v[tail-1];v[tail-1]=v[s];v[s]=t;

                        /*reverse order of elements in the tail */
                        for(i = tail, j = l - 1; i < j; i++, j--){
                                /* swap(v, i, j); */
                                t=v[i];v[i]=v[j];v[j]=t;
                        }
                        break;
                }
        }

        return encore;
}

/*-------------------------------------------------------       */
/* calculates the n+1 first factorials                          */
/* since from o! to n! there are n+1 factorials                 */
/* Minimum 3.                                                   */
/*-------------------------------------------------------       */
void init_sFacto(int n)
{
        int i;

        slnFacto=NULL;   // set the ln factorials to null

        n = n+1;
        if (n < 3)
                n = 3;

        if ((sFacto = malloc(n*sizeof(double)))==NULL) {
                fprintf(stderr, "init_sFacto: Not enough memory to compute factorials\nexiting...\n");
                exit(4);
        }

        sFacto[0]=1.0L;
        sFacto[1]=1.0L;
        for (i = 2 ; i < n ; i++) {
                sFacto[i] = i*sFacto[i-1];
        }
        return;
}

void init_slnFacto(int n)
{
        int i;

        sFacto=NULL;  // set the std factorials to null

        n = n+1;
        if (n < 3)
                n = 3;

        if ((slnFacto = malloc(n*sizeof(double)))==NULL) {
                fprintf(stderr, "init_slnFacto: Not enough memory to compute factorials\nexiting...\n");
                exit(4);
        }

        slnFacto[0]=0;
        slnFacto[1]=0;
        for (i = 2 ; i < n ; i++) {
                slnFacto[i] = log(i)+slnFacto[i-1];
        }
        return;
}

/*-------------------------------------------------------       */
/* free array of factorials (static array)                      */
/*-------------------------------------------------------       */
void free_sFacto()
{
   if (sFacto)
        free(sFacto);

	if(slnFacto)
		free(slnFacto);

   return;
}

/*-------------------------------------------------------	*/
/* initialize the static array sTableau and fill with 0's   */
/*-------------------------------------------------------	*/
void init_sTableau(int t)
{
	if (((sTableau = calloc((size_t)t, sizeof(double *)))==NULL) || ((sTableau2 = calloc((size_t)t, sizeof(double *)))==NULL)) {
		fprintf(stderr, "Not enough memory to create an array of %d integers\nexiting...\n", t);
		exit(4);
	}

	return;
}

void init_sTableau_t(int t, int n)
{
	if (((sTableau[t] = calloc((size_t)n, sizeof(double)))==NULL) || ((sTableau2[t] = calloc((size_t)n, sizeof(double)))==NULL)) {
		fprintf(stderr, "Not enough memory to create an array of %d integers\nexiting...\n", n);
		exit(4);
	}

	return;
}

/*-------------------------------------------------------	*/
/* free array of factorials (static array)			*/
/*-------------------------------------------------------	*/
void free_sTableau(int T)
{
	int t;

	if(sTableau){
		for(t=0;t<T;t++){
			if (sTableau[t])
				free(sTableau[t]);
		}
		free(sTableau[t]);
	}

	if(sTableau2){
		for(t=0;t<T;t++){
			if (sTableau2[t])
				free(sTableau2[t]);
		}
		free(sTableau2[t]);
	}

   return;
}

/* -------------------------------------------------------	*/
/* takes a vector in input, outputs the number of different values,             */
/* generates a vector of classes with a vector of lengths (i.e., probabilities) */
/* of the classes. same for ne2. Length of new vectors = nclasses               */
/* -------------------------------------------------------	*/

void set_vect_classes(int *e1m, int *e2, int nbranches, int *nclasses, int **ne1m, int **ne2, double *proba, double **nproba, int t)
{
   int i, j, last, *e1 = sTableau[t], *cl = sTableau2[t];

   /* first copy the e1m vector	*/
   memcpy((void *)e1, (void *)e1m, (size_t)nbranches*sizeof(int));
   /* trie en ordre decroissant */
   qsort(e1,nbranches,sizeof(int), (int (*)(const void * , const void * ))estpluspetit);

   /* find the number of classes and their value in cl */
   last = e1[0]; /* n will be the last value seen   */


   j = 0;
   for ( i = 1; i < nbranches; i++) {
        if (e1[i] != last) {        /* different from the last seen value   */
           cl[j]=last;         /* keep the number seen in identical values     */
           last = e1[i];            /* new value */
           j++;
        }
   }
   cl[j] = last; /* last ni[j] */

   *nclasses = j+1;        /* j est superieur au dernier indice utilise */


   /* first free old vectors	*/
   if (*ne1m)
	free(*ne1m);
   if (*ne2)
	free(*ne2);
   if (*nproba)
	free(*nproba);

   if ((*ne1m=(int*)calloc((size_t)(*nclasses),sizeof(int)))==NULL || (*ne2=(int*)calloc((size_t)(*nclasses),sizeof(int)))==NULL
	|| (*nproba=(double*)calloc((size_t)(*nclasses),sizeof(double)))==NULL) {
	fprintf(stderr, "Not enough memory to create an array (ne1m,ne2 or nproba)\nexiting...\n");
	exit(1);
   }

   /* parcours de e1m et compression dans ne1m (idem pour ne2 et addition des proba) */
   for (i = 0; i < nbranches; i++) {
	/* find the class of the i-th element	*/
        j=0;
	while (e1m[i]!=cl[j])
		j++;
	(*ne1m)[j]=cl[j];	/* in fact, ne1m is same as cl ! */
	/* (*ne1m)[j] += e1m[i]; */
	/* append the probability */
	(*nproba)[j] += proba[i];
	/* append the e2 value to the ne2 class */
	(*ne2)[j] += e2[i];
   }

#if 0
   fprintf(stderr,"e1m=");
   print_ivect(stderr,e1m,nbranches);
   fprintf(stderr,"ne1m=");
   print_ivect(stderr,*ne1m,*nclasses);
   fprintf(stderr,"ne2=");
   print_ivect(stderr,*ne2,*nclasses);
   fprintf(stderr,"nproba=");
   {
   	double s=0.0L;
   	for (i=0; i < *nclasses; i++) {
     		fprintf(stderr,"%f ",(*nproba)[i]);
     		s+=(*nproba)[i];
   	}
   	fprintf(stderr,"sumprob= %f\n",s);
   }
#endif

   return;

}

/* -------------------------------------------------------	*/
/* takes a vector in input, outputs the number of different values,             */
/* generates a vector of classes with a vector of lengths (i.e., probabilities) */
/* of the classes. same for ne2. Length of new vectors = nclasses               */
/* -------------------------------------------------------	*/

void set_vect_classes_mask(int *e1m, int *e2, char *maske1, char *maske2, int nbranches, int *nclasses, int **ne1m, int **ne2, double *proba, double **nproba, int t)
{
    int i, j, last, *e1 = sTableau[t], *cl = sTableau2[t];

    /* first copy the e1m vector	*/
    memcpy((void *)e1, (void *)e1m, (size_t)nbranches*sizeof(int));
    /* trie en ordre decroissant */
    qsort(e1,nbranches,sizeof(int), (int (*)(const void * , const void * ))estpluspetit);

    /* find the number of classes and their value in cl */
    last = e1[0]; /* n will be the last value seen   */


    j = 0;
    for ( i = 1; i < nbranches; i++) {
        if (e1[i] != last) {        /* different from the last seen value   */
           cl[j]=last;         /* keep the number seen in identical values     */
           last = e1[i];            /* new value */
           j++;
        }
    }
    cl[j] = last; /* last ni[j] */

    *nclasses = j+1;        /* j est superieur au dernier indice utilise */


    /* first free old vectors	*/
    if (*ne1m)
    free(*ne1m);
    if (*ne2)
    free(*ne2);
    if (*nproba)
    free(*nproba);

    if ((*ne1m=(int*)calloc((size_t)(*nclasses),sizeof(int)))==NULL || (*ne2=(int*)calloc((size_t)(*nclasses),sizeof(int)))==NULL
    || (*nproba=(double*)calloc((size_t)(*nclasses),sizeof(double)))==NULL) {
    fprintf(stderr, "Not enough memory to create an array (ne1m,ne2 or nproba)\nexiting...\n");
    exit(1);
    }

    /* parcours de e1m et compression dans ne1m (idem pour ne2 et addition des proba) */
    for (i = 0; i < nbranches; i++) {
    /* find the class of the i-th element	*/
        j=0;
    while (e1m[i]!=cl[j])
    	j++;
    (*ne1m)[j]=cl[j];	/* in fact, ne1m is same as cl ! */
    /* (*ne1m)[j] += e1m[i]; */
    /* append the probability */
    (*nproba)[j] += proba[i]*(int)maske1[i]*(int)maske2[i];
    /* append the e2 value to the ne2 class */
    (*ne2)[j] += e2[i]*(int)maske2[i];
    // printf("ne2[%d] = %d, e2[%d] = %d, mask[%d] = %d\n", j, (*ne2)[j], i, e2[i], i, maske2[i]);
    }
    double s=0.0L;
    for (i=0; i < *nclasses; i++) {
           if ((*nproba)[i]==0) (*nproba)[i]=1e-15;
           s+=(*nproba)[i];
    }
    for (i=0; i < *nclasses; i++) {
        (*nproba)[i]/=s;
    }

#if 0
   fprintf(stderr,"e1m=");
   print_ivect(stderr,e1m,nbranches);
   fprintf(stderr,"ne1m=");
   print_ivect(stderr,*ne1m,*nclasses);
   fprintf(stderr,"ne2=");
   print_ivect(stderr,*ne2,*nclasses);
   fprintf(stderr,"nproba=");
   {
   	double s=0.0L;
   	for (i=0; i < *nclasses; i++) {
     		fprintf(stderr,"%f ",(*nproba)[i]);
     		s+=(*nproba)[i];
   	}
   	fprintf(stderr,"sumprob= %f\n",s);
   }
#endif

   return;

}

/*-----------------------------------MULTINOMIAL--------------------------------------*/
/* the form of the multinomial is entirely via logarithms. 							  */
/* the formula is therefore: lnFacto[#e2] + SUM (log(ni)+log(Pi)) - SUM (lnFacto(ni)) */
/* the idea is to generate the multinomial recursively. At each step of the recursion */
/* we add one more log to the "denom" factor, that is the SUM(lnFacto(ni)),			  */
/* we generate all the permutations possible recursively by going through each        */
/* combination available from the number of e2 events and the number of categories    */
/* given by e1. The probabilities (SUM (log(ni)+log(Pi))) are generated by adding,    */
/* at each recursion step, the log of the branch lengths of the category considered   */
/* Everytime we go down the recursion, we remove what we previously added and we      */
/* continue the recursion.				   											  */
/*------------------------------------------------------------------------------------*/

void generate_multinomial(int *veck, double *vproba, int lveck, int sum_e2, int rest, int level, double *distrib, int* e1m, double proba, double denom)
{

	int i,k; 											/* for each category 				  	   		 	  */
	int index; 											/* the index in the distrib vector 	  	   			  */

	if (rest==0){										/* if we placed all e2 events						  */

		proba = proba-denom+slnFacto[sum_e2];			/* the probability calculation in log      			  */
					                                    /* we test whether it is <0        	                  */
        if (proba < 0.0L){
			// index = iprodsca(e1m, veck, lveck);		/* we find the hits							    	  */

			/* I'm using the version of Joel to compute the index (it is faster apparently) */
			k=0;
			index = 0;
			while(k < lveck && veck[k]==0) k++;
			while(k < lveck) {
				index += e1m[k]*veck[k];
				k++;
				while(k < lveck && veck[k]==0) k++;
			}
			distrib[index] += exp(proba);			/* we add the probability of the hits in the category */
		}
		return;
	}


	for (i=level; i<lveck; i++){					/* for all categories 									  */
		veck[i]+=1;									/* we add one in the category 							  */
		rest-=1;									/* and remove one event from the rest list (= #e2 events) */

		level=i;
		/* we need to save this index in order to go down the recursion correctly. This level index will	   */
		/* be propagated along the recursion and allows us to generate all combinations for the    			   */
		/* multinomial. 																					   */

 		if (vproba[i]>0.0L){
			proba += log(vproba[i]);					/* we add the probability of the category 				   */
        }

		if (veck[i]>0)
		{
			denom += log(veck[i]);
			/* we add the log of the count to the denom variable (the SUM(lnFacto(ni)))	from the header */
		}

		/* and we go down the recursion */
		generate_multinomial(veck, vproba, lveck, sum_e2, rest, level, distrib, e1m, proba, denom);

		/* now we are going up, we need to remove everything we added in the previous instance of the recursion */
		/* and we add one to rest since we have one event to replace now.								    */

		if (veck[i]>0)
			denom-=log(veck[i]);
		veck[i]-=1;
		rest+=1;
		if (vproba[i]>0.0L)
			proba -= log(vproba[i]);

	}
	return;
}

/* normalize a vector of length lv */
void normalize_a_vector (double *v, int lv){
    int i;
    double sum=0;
    for (i=0; i<lv; i++)
        sum += v[i];

    for (i=0; i<lv; i++)
        v[i]/=sum;

    return;
}

/* function to calculate the binomial coefficient */
long long binomial(int n, int r) {
    if(r > n - r) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++) {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}

/*	Function to check if an int vector is present in a table of vectors	*/
int isAVectorDuplicated(int *v, int lv, int **storedArray, int lastPosition){

    if (lastPosition==0)return 0;
	for (int i=0; i<lastPosition; i++){
        int sum=0;
        for (int j=0; j<lv;j++){
            if (storedArray[i][j] == v[j]) sum+=1;
        }
        if (sum==lv) return 1;
	}
	return 0;
}

/* I generate here a list of random numbers of length length in vector v*/
void listOfRandom(int *v, int length, int maxval, int cnt)
{
    int i;
    int sum=0;                      // sum of the random vector - should be maxval at end
    int tmp_maxval=maxval+1-cnt;
    // the maximum number of the range in which I pick random numbers. it changes every 10% iterations. +1 because I also need to include the latest number.

    int safemaxval=maxval;          // the maximum sum that the vector should be
    int myval;                      // the random number

    for (i=0;i<length;i++) v[i]=0;  // re-initializing the vector

    for (i=0; i<length-1; i++){

        if(maxval==0 || sum==safemaxval)
            return;
        if (tmp_maxval==0) return;                  // if I have nothing else to draw, it's over
        if (tmp_maxval==1) myval=rand()%2;          // if there is still a 1 i can draw, i draw either 1 or 0
        if (tmp_maxval>1) myval=rand()%tmp_maxval;  // I draw from 0-tmp_maxval


        v[i]=myval;                                 // position in the vector is my random number
        maxval -= v[i];                             // I subtract the value to the maximum
        if(tmp_maxval>v[i])                         // I need this to continue decreasing tmp_maxval until it is 1
            tmp_maxval -= v[i];
        sum+=v[i];                                  // And I sum the elements (as a double check)

    }

    if (maxval){                                    // If I have something left, I just add it to the last position of the vector
        v[i]+=maxval;
    }

    return;

}

/* Here I generate an approached multinomial ~ monte carlo procedure */
void generate_approached_multinomial(int *veck, double *vproba, int lveck, int sum_e2, int rest, int level, double *distrib, int* e1m, double proba, double denom)
{

	int i,k; 											/* for each category 				  	   		 	  */
	int index; 											/* the index in the distrib vector 	  	   			  */
    int count=0;
    int tablelevel=0;                                        /* the length of the table that includes the random vectors */
    int **storedVectors=NULL;
    int nrounds=100000;


    storedVectors = (int**) malloc(nrounds * sizeof(int*));
	for (int i=0; i<nrounds; i++){
		storedVectors[i] = (int*)malloc(lveck * sizeof(int));
	}
    if (!storedVectors)fprintf(stderr, "approached multinomial: cannot allocate storedVectors, exiting\n "),exit(3);

    for (int z=0; z<nrounds; z++){                                              // I'm generating vectors for z rounds (100.000)
        if (!(z%(nrounds/100)) && (count<sum_e2) && (z>1)){
            count+=1;
        }

        listOfRandom(veck, lveck, sum_e2, count);
        qsort(veck, lveck, sizeof(int), (int (*)(const void * , const void * ))estplusgrand);

        if(isAVectorDuplicated(veck,lveck,storedVectors, tablelevel)) continue; // If I already have the series of numbers stored in the array - generate the next one
        for (int j=0; j<lveck; j++) storedVectors[tablelevel][j] = veck[j]; // Otherwise, add the (sorted) vector to the table
        tablelevel+=1;

        /*
        Here, I'm performing the calculations on the generated vector. To explore as much as possible,
        I iterate over all permutations of the vector (hence the do...while loop)
        */

        do{

            proba=0.0L;
            denom=0.0L;
            for (i=0; i<lveck; i++){
                if (vproba[i]>0.0L)
    				proba += (log(vproba[i]))*veck[i];					/* we add the probability of the category 				   */

    			if (veck[i]>0){
                    for (int w=1; w<(veck[i]+1);w++) denom+=log(w);
                }
            }
            proba=proba-denom+slnFacto[sum_e2];


            // if (exp(proba) != 0.0L){						/* we exponentiate it and test whether it is >0 	  */
            if (proba < 0.0L){                              /* we test whether it is <0 	  */
                /* I'm using the version of Joel to compute the index (it is more efficient) */
                k=0;
                index = 0;
                while(k < lveck && veck[k]==0) k++;
                while(k < lveck) {
                    index += e1m[k]*veck[k];
                    k++;
                    while(k < lveck && veck[k]==0) k++;
                }
                distrib[index] += exp(proba);			/* we add the probability of the hits in the category */
            }

        } while(get_next_permut(veck, lveck) == 1);
    }
    for (int j=0; j<nrounds; j++){

        free(storedVectors[j]);
    }
    free(storedVectors);
    return;
}
