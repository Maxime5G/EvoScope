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
\file 	tree.c
\brief 	basic tree metrics and manipulations
\author g achaz, j pothier madamesophie
\date 	June 28 2012
**/


#include <stdio.h>
#include <stdlib.h>

#ifdef LINUX
#include <strings.h>
#endif
#include <string.h>

#include <math.h>
#include "tree.h"
#include "coevol.h"

/**
\fn *read_tree_forest
\brief read a newick string and return the tree
\param newick_string: a newick string
\param opt_read_stdin: if need to open treefile or stdin
\param nbleaves: Forest.nbleaves, #leaves of the trees in the forest
\param outgroup: outgroup of tree(s)
**/

Node *read_tree_forest( char * newick_string, int *nbleaves, char *outgroup )
{

	Node *init_tree,
	     *root;

	int n, nnodes;

	/*
		Set tree in memory
	*/

	/* init */
	init_tree = InitMultiTree(*nbleaves);

	if (verbose) fprintf(stderr, "Total number of leaves : %d\n",*nbleaves);

	/* reads a newick string and stores it in the tree */
	n = 0;
	// printf("newick_string = %s\n", newick_string);

	lit_newickmyMulti( newick_string, 0, &n, init_tree );
	nnodes = n+1;   // the last slot corresponds to #nodes -1

	/* set a root name if there isn't one */
	root=FindRoot( init_tree );

	if (! root->name){
		if ((root->name=(char *)malloc(sizeof(char)*5))==NULL)
			fprintf(stderr,"Not enough memory for name of root node\n"), exit(1);
		strcpy(root->name,"root");
	}

	return init_tree;
}

/*---------------------------------------------------*/
/**
\fn	struct  mymultinode *InitMultiTree( int nleaves ,struct mymultinode *myNodes)
\brief init and alloc all the struct needed. return it
**/
struct mymultinode *InitMultiTree(  int nleaves ){

	struct mymultinode *myNodes;
	int i;
	int nbt=1+(2 * nleaves-1);
	myNodes=(struct  mymultinode *)malloc( sizeof(struct mymultinode )*nbt );
	if(!myNodes)fprintf(stderr, "InitTree: cannot allocate myNodes, bye\n"), exit(1);

	for(i=0; i< nbt; i++){                                         /* initialise les noeuds en mettant les parametres par defaut pour pouvoir ensuite les remplir */
		myNodes[i].nevt= 0;
		myNodes[i].evt= NULL;
		myNodes[i].anc = NULL;
		myNodes[i].nbdesc=0;
		myNodes[i].descs = NULL;
		myNodes[i].time = 0.0L;
		myNodes[i].id = -1;
		myNodes[i].nbMuts=-1;
		myNodes[i].which_specie=-1;
		myNodes[i].possMut=NULL;
		myNodes[i].misc=0;
		if ((myNodes[i].SeqInLeaves = (int *)calloc( nleaves, sizeof(int))) == NULL) {
			fprintf(stderr,"Not enough memory for SeqInLeaves, node %d\n",i);
			exit(1);
		}
		if ((myNodes[i].name=(char *)malloc(sizeof(char)*16))==NULL) {
			fprintf(stderr,"Not enough memory for name of node %d\n",i);
			exit(1);
		}
		strcpy(myNodes[i].name,"internal-node");
	}


	return myNodes;
}





void PrintNode( Node *n, char *name )
{
	int i;

	printf("\nNode %s\n", name);

	printf("     id:        %d\n", n->id);
	printf("     time:      %f\n", n->time);
	if(!n->anc)
	printf("     anc:       (null)\n");
	else
	printf("     anc:       %d\n", n->anc->id);
	printf("     nbdesc:    %d\n", n->nbdesc);
	for(i=0;i<n->nbdesc;i++)
	printf("       desc[%d]: %d\n", i, n->descs[i]->id);
	printf("     name:      %s\n", n->name);

}


/* ----------------------------------------------------------------------------	*/
/* count leaves from node n					*/
/* ----------------------------------------------------------------------------	*/
int count_leaf(Node *n)
{
	int nleaf=0, i;


	if (n->nbdesc == 0) {
		nleaf = 1;
	}
	else {
		for (i = 0; i < n->nbdesc; i++) {
			nleaf += count_leaf( n->descs[i] );
		}
	}
	return nleaf;
}

/* ----------------------------------------------------------------------------	*/
/* count nodes from node n					*/
/* ----------------------------------------------------------------------------	*/
int count_nodes(Node *n)
{
	int nnodes=1, i;

	if (n->nbdesc > 0)
		for (i = 0; i < n->nbdesc; i++)
			nnodes += count_nodes( n->descs[i] );

	return nnodes;
}



/* ----------------------------------------------------------------------------	*/
/* count branches from node n					*/
/* but not the root one !						*/
/* ----------------------------------------------------------------------------	*/
int count_branches(Node *n)
{

  int nbranches=0;
  int i;

  for (i = 0; i < n->nbdesc; i++)
  {
  		nbranches += ( count_branches (n->descs[i]) + 1);
  }

 return nbranches;
}



/* --------------------------------------------------	*/
int *joinLeaves(int *s1, int* s2, int N)
{
	int *res,i;

	res=malloc(sizeof(int )*N);                 /* met dans un tableau de pointeur le nombre de noeuds et de feuilles ?*/

	if (res==NULL) {
		fprintf(stderr, " Erreur de malloc, bye\n");
		exit(2);
	}
	for (i=0;i<N;i++)
		if (s1[i]==0 && s2[i]==0)
			res[i]=0;
		else
			res[i]=1;

	return res;
}

void recopie(int *dest, int* src, int N)
{
	int i;
	for (i=0;i<N;i++)
		dest[i]=src[i];                       /* crée deux tableaux de pointeurs de meme taille ? */
}

/*
	Must be called after initseqinleaves which initiates seqinleaves for leaves
*/
void BuildSeqInLeaves(Node *pptree,int N)                 /* ????????????? */
{

	if (pptree->descs[0]->nbdesc!=0)
		BuildSeqInLeaves(pptree->descs[0],N);

	if (pptree->descs[1]->nbdesc!=0)
		BuildSeqInLeaves(pptree->descs[1],N);

	recopie(pptree->SeqInLeaves ,joinLeaves(pptree->descs[0]->SeqInLeaves,pptree->descs[1]->SeqInLeaves,N),N);


}



/*---------------------------------------------------*/
/**
\fn 	void clean_end(str)
\brief	remove the value at the end of newick string if present because unexpected error otherwise
\brief	*Warning*: this function changes the value of str
**/
void clean_end(char *str)
{
	int l=strlen(str);
	int i;
	for (i=l-2;i>0;i--) /*str[l-1]=';' so test previous */
	{

	if (str[i]==')')
		break;
	}
	if (i==0) printf("pb cleaning newick string....\n"),exit (11) ;
	if (i!=l-2)
	{

	str[i+1]=';';
	str[i+2]='\0';
	}

}




/**
\fn    double TotalLength( Node * mynode )
\brief Sum from mynode TopDown all mynode->time, including the current branch
**/
double SumLength( Node * mynode ){

	double L=(mynode->anc==NULL)?0:mynode->time;
	int i;

	if(mynode->nbdesc > 0)
	{
		for(i=0;i<mynode->nbdesc; i++)
			L += SumLength( mynode->descs[i] );
	}

	return L;
}

/**
\fn    double TotalLength( Node * mynode )
\brief Sum from Root TopDown all mynode->time
**/
double TotalLength( Node * mynode ) {
	return SumLength( FindRoot( mynode ) );
}

double SubTreeLength( Node * mynode ) {
	return SumLength( mynode )-mynode->time;
}



/**
\fn 	countchar(char *ch)
\brief 	returns a string containing a newick name
\param  ch is the newick string begining a name,
\return the name
**/
char *countchar(char *ch)                                                                 /* Permet de creer un tableau de pointeur de chaine de caractere contenant un nom */
{
	int nb=0;
	char *nom;
	while ((* (ch+nb)!=':')&& (* (ch+nb)!=')')&& (* (ch+nb)!=',') && (* (ch+nb)!='[') && (* (ch+nb)!=']')) nb++;

	nom = malloc(sizeof(char) * (nb +1) );
	if ( !nom )
		fprintf(stderr, "countchar: cannot allocate nom, bye\n"), exit(1);

	nb=0;

	while ((* (ch+nb)!=':') && (* (ch+nb)!=')')&& (* (ch+nb)!=',') && (* (ch+nb)!='[') && (* (ch+nb)!=']')) {
		nom[nb]=ch[nb];
		nb++;
	}

	nom[nb]='\0';

	return nom;
}




/*---------------------------------------------------*/
/**
\fn 	binarise (struct mymultinode *myNodes,int nbnodes)
\brief	if struct multinodes have more than 2 nodes change it for a only '2 nodes by struct'
		by creating more nodes
\param	myNodes array of struct nodes
\param 	nbnodes number of nodes before binarisation
\returns the new number of nodes
**/
int binarise (struct mymultinode *myNodes,int nbnodes)
{
	int i,j,o=0, lastnode=nbnodes,ledesc=-1;

	struct mymultinode 	*temp; /* new node */

	for (i=0;i<nbnodes;i++)
	{
/*		printf("node %d myNodes[0].descs[1]=%d\n",i,myNodes[0].descs[1]); */
		if (myNodes[i].nbdesc>2)
			{
/*		printf("$$$$node %d myNodes[0].descs[1]=%d\n",i,myNodes[0].descs[1]); */
			for (j=2;j<myNodes[i].nbdesc;j++)
				{
				if (myNodes[i].anc!=NULL)
					{
					if (myNodes[i].anc->descs[0] == &myNodes[i])
						ledesc=0; else ledesc=1;
					}

				temp=&myNodes[lastnode]; /* new node */
				/* always only 2 children */
				temp->descs=realloc(temp->descs, sizeof (struct mymultinode) *2);
				temp->nbdesc=2;
				temp->descs[1]=myNodes[i].descs[j];
				temp->descs[0]=&myNodes[i];
				temp->anc=myNodes[i].anc;
				temp->id=lastnode;

				temp->time=myNodes[i].time;

/*			printf("----\ni=%d j=%d-> &myNodes[i]=%d temp=%d temp->descs[0]=%d temp->descs[1]=%d myNodes[i].descs[0]=%d myNodes[i].descs[1]=%d (myNodes[i].anc=%d)\n",i,j,&myNodes[i],temp,temp->descs[0],temp->descs[1],myNodes[i].descs[0],myNodes[i].descs[1],myNodes[i].anc); */
				if (myNodes[i].anc!=NULL)
					myNodes[i].anc->descs[ledesc]=temp;

				myNodes[i].anc=temp;
				myNodes[i].descs[j]->anc=temp;
				myNodes[i].time=0;
/*			    printf("i=%d j=%d-> &myNodes[i]=%d temp=%d temp->descs[0]=%d temp->descs[1]=%d myNodes[i].descs[0]=%d myNodes[i].descs[1]=%d (myNodes[i].anc=%d) myNodes[0].descs[1]=%d\n\n",i,j,&myNodes[i],temp,temp->descs[0],temp->descs[1],myNodes[i].descs[0],myNodes[i].descs[1],myNodes[i].anc,myNodes[0].descs[1]); */

				lastnode++;
				}



			myNodes[i].nbdesc=2;

			}

	}
/*	printf("check binarise %d %d\n",nbnodes,lastnode);		*/
	for (i=0;i<lastnode;i++)
		if (myNodes[i].nbdesc>2)
			o++;
	if (o>0) printf("Tree was not binarizzzzed\n");

return(lastnode);


}




/**
\fn	void SplitLeaves( struct mymultinode *pptree, int *ind, int ind_size );
\brief	in case of unresolved tree do the job

**/

void SplitLeaves( struct mymultinode *pptree, int *ind, int ind_size )  {

	int i,
	    t=0,           /* #of the sequences with id in ind equal to this node->id */
	    max=ind[0],    /* the last used node */
	    first=-1;      /* first sequence with id in ind equal to this node->id */

	if (pptree->descs[0]==NULL && pptree->descs[1]==NULL ){   /* thats a leaf*/

		for (i=0;i<ind_size;i++){

			if (ind[i]==pptree->id){
				t++;

				if(first==-1)
					first=i;

			}

			max=(ind[i]>max)?ind[i]:max;

		}
/*	printf("--->%d %d (%d)\n",max,t,pptree->id);*/
		if(t>1){

			pptree->nbdesc=2;  /* now this node has 2 desc(s) */

			/* mettre la première sequence dans le noeud à l'adresse max-id+1
			   connect noeud à l'adresse max-id+1 avec pptree (descs[0], temps=0)
			   ind[ sequence 1] = max+1
			   SeqInLeaves de noeud à l'adresse max-id+1 vaut 1 pour sequences 1
			   Fill id for noeud à l'adresse max-id+1
			*/

			pptree->descs[0] = pptree + max-pptree->id+1;
			pptree->descs[0]->id = max+1;
			pptree->descs[0]->time = 0.0L;
			pptree->descs[0]->anc = pptree;
			pptree->descs[0]->SeqInLeaves[first] = 1;
			ind[first] = max+1;

			/*
			   mettre tous les autres dans le noeud à l'adresse max-id+2
			   connect noeud à l'adresse max-id+2 avec pptree (descs[1], temps=0)
			   pour tous les ind[ autre_seq ] = max+2
			   SeqInLeaves de noeud à l'adresse max-id+2 vaut 1 pour ttes les sequences
			   Fill id for noeud à l'adresse max-id+2
			*/
/*				printf("max+2 (%d)\n",max+2);      */
			pptree->descs[1] = pptree + max-pptree->id+2;
			pptree->descs[1]->id = max+2;
			pptree->descs[1]->time = 0.0L;
			pptree->descs[1]->anc = pptree;

			for(i=first; i<ind_size; i++)
			{
				if( ind[i]==pptree->id ){
					pptree->descs[1]->SeqInLeaves[i] = 1;
					ind[i] = max+2;
				}

			}
/*		printf("--->(%d)\n",pptree->descs[1]->id); */
			SplitLeaves( pptree->descs[1], ind, ind_size );

		}

		return;


	}

	SplitLeaves( pptree->descs[0], ind, ind_size );
	SplitLeaves( pptree->descs[1], ind, ind_size );

}



/**
\fn	void connect(struct mymultinode  *monNoeud, int anc, int desc1, int desc2)
\brief	Connecting Nodes
\param 	monNoeud the  node to connect
\param 	anc index of ancestor to connect to monNoeud
\param	descs[0] index of 1st descendant to connect to connect to monNoeud
\param descs[1] index of 2nd descendant to connect to connect to monNoeud
**/
void connect(struct mymultinode  *monNoeud, int anc, int desc1, int desc2)
{

	monNoeud[anc].id = anc;
	monNoeud[anc].descs[0] = monNoeud+desc1;
	monNoeud[anc].descs[1] = monNoeud+desc2;

	monNoeud[desc1].id = desc1;
	monNoeud[desc1].anc = monNoeud+anc;

	monNoeud[desc2].id = desc2;
	monNoeud[desc2].anc = monNoeud+anc;
	monNoeud[anc].nbdesc=2;

}



/**
\fn	void Connect( Node * AncNode , Node * DescNode , float time )
\brief	Connecting Desc to Anc, while decreasing the ex ancestors by 1
\param 	monNoeud the  node to connect
\param 	AncNode is the new ancestor
\param	DescNode is the new descendant
\param  time sets to the DescNode.time value
**/
void Connect( Node * AncNode , Node * DescNode , float time )
{

	/*
		reallocate the memory and set pointers
	*/
	if(AncNode->nbdesc == 0)
	 	AncNode->descs = ( Node ** )malloc( (size_t) 1*sizeof( Node * ) );
	else
	 	AncNode->descs = ( Node ** )realloc( (void *)AncNode->descs, (size_t) (AncNode->nbdesc+1)*sizeof( Node * )  );

	AncNode->descs[  AncNode->nbdesc ++  ] =  DescNode;


	/*
		Remove the link from the previous anc to the DescNode [if any]
	*/
	if(DescNode -> anc != NULL){

		Node * ex_anc = DescNode -> anc;
		int i=0;

		while( ex_anc->descs[i] != DescNode )   /* find the link that need to be removed */
			i++;

		while( i < ex_anc->nbdesc-1 ){
			ex_anc->descs[i] = ex_anc->descs[i+1];  /* move the link 1 down above the chosen link */
			i++;
		}

		/*
			free memory and dec the nbdesc
		*/
		ex_anc->descs = ( Node ** )realloc( (void *)ex_anc->descs, (size_t) (ex_anc->nbdesc-1)*sizeof( Node * )  );
		ex_anc->nbdesc--;

	}

	DescNode -> anc = AncNode;
	DescNode -> time = time;


}



/**
\fn int isAnc(struct mymultinode  *Test,struct mymultinode  *Ancestor)
\brief	test if Ancestor is a ancestor of Test
\param Test the node we want to know the ancestor
\param Ancestor the putative ancestor
\return 1 if test is descendant 1 of Ancestor, return 2 if descendant 0
		0  if Ancestor not the ancestor of Test
**/

int isAnc(struct mymultinode  *Test,struct mymultinode  *Ancestor)
{
	while (Test != NULL && Test->anc!=NULL)
	{
		/*printf("test is %d; Anc is %d\n", Test->id, Test->anc->id);*/

		if (Test->anc==Ancestor)
		{
			if( Test->anc->descs[0]== Test)
				return 1;
			else
				return 2;
		}
		else
			Test=Test->anc;
	}

	return 0;
}


/**
\fn int HowManyAnc(struct mymultinode  *Test)
\brief	test how many ancestor have Test (how many generation)
\param	Test the node to test
\return	the number of ancestors
**/
int HowManyAnc(struct mymultinode  *Test)
{
	int x=0;

	while (Test->anc != NULL)
	{
		Test=Test->anc;
		x++;
	}

	return x;
}

/*

*/
void FreeNode( Node *n ){


		if ( n->SeqInLeaves) {free( n->SeqInLeaves ); n->SeqInLeaves=NULL;}
		if ( n->name)        {free( n->name ); n->name=NULL;}
		if ( n->descs )      {free( n->descs ); n->descs=NULL;}
		if ( n->possMut )    {free( n->possMut ); n->possMut=NULL;}
		if ( n->evt )        {free( n->evt ); n->evt=NULL;}

		n->descs=NULL;
		n->nbdesc=0;
		n->time=0.0L;


}

/*
	Free the whole tree
*/
void FreeTree( struct mymultinode  * myNodes, int n_nodes ){

	int i;
	for (i=0; i<n_nodes ;i++)
		FreeNode( myNodes+i );

	free( myNodes );
}

void FreeMultiTree(struct mymultinode *myNodes, int nleaves)
{
   int i;
   int n_nodes=1+(2 * nleaves-1);

   if (myNodes) {
      for (i=0; i<n_nodes ;i++)
		FreeNode( myNodes+i );
      free( myNodes );
   }
   return;
}

/* This function reads the tree and stores the branch lengths in the vector branch_lengths. */
/* INPUTS: the root node, the number of branches of the tree, the output vector and the     */
/* tree id.																				    */
/* This function also stores the total branch length of the tree. 							*/
/* If runIsEpics, the branch lengths are normalized to sum to 1 (i.e. every branch is 		*/
/* divided by the total branch length).    													*/

double *ComputeBranchLengths( Node *root, int nbranches, int output_branch_probas, double *total_length_all_tree, int t, int runIsEpics){

	double total_length;
 	IntegerMatrix *mateve=NULL;
 	double *branch_lengths=NULL;
 	double pval=0;
 	int i;


	/* allocate a static vector in mat.c that will serve to everything */
	init_sTableau_t(t, nbranches);


	/* allocate the matrix num_branches*num_event (number of events taken on the first righthand son) */
	mateve = InitIntegerMatrix(root->nevt, nbranches );

	/* fills the matrix and the branch_lengths vector */
	if ((branch_lengths = malloc(sizeof(double)*(nbranches))) == 0)
		fprintf(stderr, "Memory error for time vector allocation, exit\n"),exit(3);


	if (verbose)fprintf(stderr, "->Filling event matrix\n");

	total_length = fill_evt_matrix(root, mateve, branch_lengths, nbranches);
	if (!(runIsEpics))
		*total_length_all_tree += total_length;

	if (total_length == 0.0)
		fprintf(stderr, "Problem: total length of branches is 0, exiting\n"), exit(2);

	if (runIsEpics){
		/* Normalization of branch probabilities, stored in the branch_lengths vector */
		pval = 0.0L;
		for (i = 0 ; i < nbranches; i++) {

	  		if (output_branch_probas) fprintf(stderr, "branch %3d length=%g ",i,branch_lengths[i]);

	  		branch_lengths[i] /= total_length;

	  		if (output_branch_probas) fprintf(stderr, "proba=%g\n",branch_lengths[i]);

			pval+=branch_lengths[i];

		}

		if (output_branch_probas) fprintf(stderr, "Sum of probabilities of branches: %g\n",pval);
	}

	FreeIntegerMatrix(mateve);

	return branch_lengths;
}


/*
	Reset any branch smaller than min_length to min_length
*/
void SetMinimumBranchLength( Node *n, double min_length ){

	int i;

	if( n->anc != NULL && n->time <= min_length )
		n->time = min_length;

	for (i = 0; i < n->nbdesc; i++)
		SetMinimumBranchLength( n->descs[i], min_length );

}

/*
	Rescale the whole tree to Ltot=1
*/
void NormalizeTree( Node *n, double Ltot ){

	int i;

	if( n->anc != NULL )
		n->time /= Ltot;

	for (i = 0; i < n->nbdesc; i++)
		NormalizeTree( n->descs[i], Ltot );

}

/*
	Rescale the whole tree to Ltot=1 and extract longest branch
*/
void NormalizeTreeAndLongestBranch( Node *n, double Ltot, double *longest){

	int i;

	if( n->anc != NULL ){

		n->time /= Ltot;

		if(n->time > *longest){
			*longest = n->time;

		}
	}


	for (i = 0; i < n->nbdesc; i++)
		NormalizeTreeAndLongestBranch( n->descs[i], Ltot, longest );

}

/*
	leaves_str must be a string with either a single leaf or leaf1:leaf2:leaf3:...

	size_array is the number of leaves in the returned array (evaluated within the fonction)

*/
int *GenerateLeavesArray( char *leaves_str, Node *node_array, int nnodes, char opt_invert, int *size_array){

	int *leaves_array;

	int nleaves,
	    n_leaves_str=1,
	    p=0,
	    len_leaves_str = strlen(leaves_str);

	Node *root;

	int i, n, s;

	root    = FindRoot( node_array );
	nleaves = count_leaf( root ) ;

//	printf("#leaves: %d\n", nleaves);


	for(i=0; i<len_leaves_str; i++)
		if ( leaves_str[i] == ':' )
		{
			n_leaves_str++;
			leaves_str[i] = 0;
		}


	leaves_array = (int *)malloc( (size_t) nleaves * sizeof(int) );
	if (!leaves_array)fprintf(stderr, "GenerateLeavesArray: cannot allocate leaves_array, bye\n"), exit(3);

//	int cc=0;

	for( s=0, n=0 ; n<nnodes ; n++ )
	{
		// foreach leaf
		if( node_array[n].nbdesc == 0 )
		{
//			cc++;
//			printf(">> %s\n", node_array[n].name);

			// scan all names in leaves_char
			for( i=0, p=0 ; i<n_leaves_str ; i++ )
			{
//				printf("name %s vs node %s\n", node_array[n].name , leaves_str+p);


				if ( strcmp( node_array[n].name , leaves_str+p ) == 0 )
					break;
				else
				{
					do{ p++; } while ( *(leaves_str+p) != 0 );
					p++;
				}
			}

//			if( i<n_leaves_str )printf(">> Found !! %d\n", s);

			if( i==n_leaves_str && opt_invert ) 		leaves_array[ s++ ]=n;
			else if( i<n_leaves_str && opt_invert==0 )	leaves_array[ s++ ]=n;

		}


	}

	/*
		For outgroups that did not match any node name
	*/
	*size_array=s;
	leaves_array = (int *)realloc( (void *) leaves_array, (size_t) (*size_array) * sizeof(int) );
	if (!leaves_array)fprintf(stderr, "GenerateLeavesArray: cannot reallocate leaves_array, bye\n"), exit(3);

//	printf("#leaves: %d ; #nodes: %d\n", cc, n);


	/*
		restore the original str
	*/
	for(i=0; i<len_leaves_str; i++)
		if ( leaves_str[i] == 0 )
			leaves_str[i] = ':';


	return leaves_array;
}
