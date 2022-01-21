/**
\file 	readwrite.c
\brief 	everything for reading and outputting generic trees in newick format
\author g achaz, j pothier and madamesophie
\date 	2018
**/


#include <stdio.h>
#include <stdlib.h>

#ifdef LINUX
#include <strings.h>
#endif

#include <string.h>

#include <math.h>
#include "tree.h"


/*
	return how many char were read
*/
static int readBEAST( char *s , Node *node ){

	int ncharlu=0,
		pos=0;

	float rate=0.0;

	printf("I am in BEAST %c%c%c\n",s[0], s[1], s[2]);

 /*
 	check whether the opening '[' ends with a ']'
 */
  while (s[pos]!=']' && s[pos]!='\0' ) {
     pos++, ncharlu++;
  }

  if (s[pos]=='\0') {
   fprintf(stderr, "readBEAST:: error: not opening \"[\" in name \n");
   exit(1);
  }

	pos=0;
	while( s[pos] != ']' && strncmp(s+pos,"rate=",5) != 0 )
		pos++;

	if( s[pos] != ']' )
		sscanf(s+pos+5, "%f", &rate);

	//printf("%f\n", rate );
	node->misc=rate;

	while( s[pos] != ']' )
		pos++;

	return ncharlu;

}


/*-------------------------------------------------------------------	*/
int set_evenement (char *s, Node *noeud, int *fl)
{
  char temps[10];
  int i=0 , pos=0, j=0, slash=0, ncharlu=0;
  static int nbTotalEvenement=0;

 /* Verification de la chaine de caractere */

  while (s[pos]!='['  && s[pos]!='\0' ) {
     pos++, ncharlu++;
  }

  if (s[pos]=='\0') {
   fprintf(stderr, "set_evenement:: error: not opening \"[\" in name \n");
   exit(1);
  }


/* compter les points slash*/

  while (s[pos]!=']'  && s[pos]!='\0') {
     if (s[pos]== '/') {
       slash++;
     }

     pos++, ncharlu++;
  }

  if (s[pos]=='\0') {
   fprintf(stderr, "set_evenement:: error: no closing \"]\" in name\n");
   exit(1);
  }


  /*printf("nbTotalEvenement =  %d \nNombre de slash = %d\n",nbTotalEvenement,slash);*/

  if (slash!=nbTotalEvenement && (nbTotalEvenement!= 0 || slash != 0)) {
	  //fprintf(stderr, "Warning: event numbers differ (at tree chain=%s)\n",s);
	  // exit(1);
  }
  nbTotalEvenement=slash;


  slash=slash  +1 ;
  noeud->nevt=slash;

#if 0
  fprintf(stderr,"Nombre d'evt=%d\n",noeud->nevt);
  fprintf(stderr,"La liste contient %d caracteres (apres le '[')\n",ncharlu);
  fprintf(stderr,"La liste contient %d nombres\n",slash);
#endif

/*conversion*/

  noeud->evt=malloc((noeud->nevt)*sizeof(int));
  if (noeud->evt==NULL) {
    fprintf(stderr,"Malloc error for event list (node %s),exiting\n",noeud->name);
    exit(11);
  }

  pos=0;
  while (s[pos]!='['  && s[pos]!='\0' ) {
     pos++;
  }
  do {

       pos++;
       if (s[pos]!='/' && s[pos]!=']') {
          temps[i]=s[pos];
          i++;
       }
       else {
           temps[i]='\0';
           sscanf( temps,"%d",&(noeud->evt[j]));
	   if (noeud->evt[j] < -1) {
	      fprintf(stderr, "node %s, event %d is negative (%d) !\nexiting...\n",noeud->name,j+1,noeud->evt[j]);
	      exit(5);
	   }
	   if ((noeud->evt[j] < 0) && (*fl)==0) {
		   fprintf(stderr, "you have negative events in your tree.\n In this situation, epics/epocs considers this as unknown!\n");
		   fprintf(stderr, "first occurrence = node %s, event %d (%d)\n", noeud->name,j+1,noeud->evt[j]);
		   (*fl)=1;
	   }
           i=0;
           j++;
#if 0
           printf("Liste====>%s\n",temps);
#endif
       }

   } 	while (s[pos]!='\0' && s[pos]!=']');

   return ncharlu;

}


/* --------------------------------------------------	*/
/* print_arbre: imprime l'arbre a partir du noeud n */
/* --------------------------------------------------	*/

void print_arbre(Node *n, FILE *f)
{
   static int branchnum = -2, recurs = 0;
   int i,j;

   branchnum++;
   recurs++;

   if (branchnum == -1) {
	fprintf(f,"[%3d;id:%3d] %s\n",branchnum, n->id,n->name);
   }
   else {
	   for (i = 0; i < recurs; i++)
		   fprintf(f,"\t");
	   fprintf(f,"|->[%3d;id:%3d] %s (len %f)",branchnum, n->id,n->name,n->time);
	   if (n->evt) {fprintf(f," events=[");for (j=0;j<n->nevt;j++) {fprintf(f,"%d,",n->evt[j]);};fprintf(f,"]");}
	   fprintf(f,"\n");
   }

   if (n->nbdesc > 0)
	   for (i = 0; i < n->nbdesc; i++)
		   print_arbre (n->descs[i], f);

   recurs--;

   return;
}




/**
\fn 	static void multitreeoutNck( struct multinode *pptree )
\brief 	print the tree but not the ";\n" at end. so it should be done outside
\param  pptree array containing all described in struct multinode
**/

void printevents( FILE *f, int *evt, int n ){

	int i=0;

	if (n > 0) {
		fprintf(f,"[");
		for (i=0;i<n-1;i++) {
			fprintf(f, "%d/",evt[i] );
		}
		printf("%d]",evt[i]);
	}

}




void multitreeoutNck( struct mymultinode *pptree, FILE *f )
{
	int i;
#if 0
	fprintf(stderr,"multitreeouNck-->%s (desc:%d)\n", pptree->name, pptree->nbdesc);
#endif
	if (pptree->nbdesc==0)  {                                          /* Permet d'imprimer un arbre de la forme (A:5.0[1/2/3],(C:6.0[3/1/6],D:2.0[8/-1-4])B:3.0[26/5/98]); */

		fprintf(f,"%s", pptree->name);
		printevents( f, pptree->evt, pptree->nevt );
		if (pptree->time!=-1.0L)
			fprintf(f,":%f", pptree->time);

	}
	else {
		fprintf(f,"(");
		for (i=0;i<pptree->nbdesc;i++) {
			multitreeoutNck( pptree->descs[i], f );
			if (i!=pptree->nbdesc -1)
				fprintf(f,",");
		}

			fprintf(f,")");
			printevents( f, pptree->evt, pptree->nevt );
			if (pptree->time!=-1.0L && pptree->anc != NULL )
				fprintf(f,":%f",pptree->time);
	}

}

/*---------------------------------------------------*/
/**
\fn 	readFileMultiNwck(FILE *f,int *nbleaves)
\brief 	read a multi-newick file , count leaves nbr and return the strings
\param 	f a pointer on a file
\param 	nbleaves how many leaves will be in each tree of this file (modified inside fn) call with 0
\returns array with pointers to strings containing the newick description (initialised in this fn)
**/

char **readFileMultiNewick(char *treefile, int opt_read_stdin, int **nbleaves, int *ntrees)
{
	char c=0;
	char *str=NULL;
	int count=0;
	int *l;
	int nbp=0;
	int j=0;
	int maxlength=0;
	char **strout;
	int i;
	FILE *f;
	int flag;

	f = fopen(treefile, "r");
	while (!feof(f)){
		c=fgetc(f);
		count++;

		if (c=='\n'){
			(*ntrees)+=1;
			if (count>maxlength) maxlength=count;
			count=0;
		}
	}

	(*nbleaves) = calloc( (*ntrees), sizeof(int));
	l = calloc( (*ntrees), sizeof(int));
	strout = (char**) malloc( (*ntrees) * sizeof(char*) );
	str = malloc( sizeof(char) * maxlength+1);

	rewind(f);
	count=0;

	while (!(feof(f))){
		fgets(str, maxlength+1, f);
		nbp=0;
		j=0;
		c=0;
		while (c!=';' && c!='\n'){
			c=str[j++];

			l[count]++;

			if (c==',')
				(*nbleaves)[count]+=1;
			if (c=='(')
				nbp++;
			else
				if (c==')')
					nbp--;
		}
		(*nbleaves)[count]+=1;
		count++;

		if (nbp!=0) printf("unbalanced parenthesis - check newick syntax\n"),exit(1);
	}
	for (i=0; i<(*ntrees); i++){
		strout[i] = malloc(sizeof(char) *l[i] +1);
	}
	if (*strout==NULL){
		printf("malloc error in readFileMultiNwck\n"), exit(1);
	}

	rewind(f);
	count=0;

	while (count<*ntrees){
		flag=0;
		fgets(str, maxlength+1, f);
		c=0;
		j=0;
		while (c!=';' && c!='\n'){
			c=str[j];
			if (c!='\n' && c!=EOF){
				strout[count][j++] = c;
			}
			if (c==';'){flag=1;}
		}
		if (flag==0){
			fprintf(stderr, "Warning: No ; at the end of newick string. Adding it.");
			strout[count][j++]=';';
		}
		strout[count][j]='\0';	// End of string
		if (strout[count][j-2]!=')')
			clean_end(strout[count]);
		count++;
	}

	free(str);
	free(l);

	return(strout);
}


/*---------------------------------------------------*/
/**
\fn 	readFileNwck(FILE *f,int *nbleaves)
\brief 	read a newick file , count leaves nbr and return the string
\param 	f a pointer on a file
\param 	nbleaves how many leaves will be in this file (modified inside fn) call with 0
\returns  pointer to a string containing the newick description (initialised in this fn)
**/

char *readFileNwck(FILE *f, int *nbleaves)
{
	int l=0;   /* nbleaves is nbr of comma always a leaf even if no comma */
	char c=0;
	int nbp=0; /* Nbp = Number of parenthesis*/
	int fl=0;
	char *str;

	while (c!=';' && !feof(f))   /* Feof= file End of File */
	{
		c=fgetc(f);
		l++;                /* compte le nombre de caracteres du fichier*/

		if (c==',')
			(*nbleaves)++;  /*compte le nombre de feuille en comptant le nombre de virgules*/

		if (c=='(')            /* compte le nombre de parentheses ouvrantes */
			nbp++;
		else
			if (c==')')        /* qu'on soustrait avec le nombre de parentheses fermantes */
				nbp--;
	}

	if (nbp!=0) printf("unbalanced parenthesis - check newick syntax (nb leaves %d, nb parenthesis %d)\n",*nbleaves,nbp),exit(1);    /* Test d'équilibrage correct des parenthèses*/
	str=malloc(sizeof (char) *(l+1));
	if (str==NULL) {
		fprintf(stderr,"Erreur de malloc, bye\n");
		exit(3);
	}              /* initialise la taille du pointeur a prévoir */

	l=0;
	rewind(f);                                     /* retour au début du fichier */
	c=0;
	while (!feof(f) && c!=';')
	{
		c=fgetc(f);
		if (c!='\n'&& c!=EOF)                                  /* ???*/
			str[l++]=c;    /* skip any CR char */
		if (c==';')
			fl=1;
	}
	if (fl==0) {fprintf(stderr,"*warning*:\";\" added at end of newick string\n");str[l++]=';';}
	str[l]='\0';


	if (str[l-2]!=')')
		clean_end(str);

	(*nbleaves)++;

	return str;
}

/*---------------------------------------------------*/
/**
\fn int lit_newickmyMulti(char *chaine,int la_pos,int *which_node,struct mymultinode *myNodes)
\brief	read a newick string and put it in an array of struct
\param	chaine the newick string
\param	la_pos pointer on string chaine, Must be 0 at 1st call
\param 	which_node pointer with the number of the struct multinode to fill
\param 	myNodes array of struct multinodes
\return	position of string during reading
**/

/*
	annotations are typically provided in [] after the node
	for epics/epocs , evolutionnary events are given as [0/1/1...]
	for BEAST, some annotations are given as [&height=..., rate=..., length=... ]
	thus, depending on '[' vs '[&', it reads it as epoics style or BEAST style
*/

int lit_newickmyMulti( char *chaine, int la_pos, int *which_node, struct mymultinode *myNodes )
{
	static int ok_for_name = 0;
	int leaf=*which_node,
		which_desc=0;

	float t=-1;
	char *nom;
	static int fl=0;

	// printf("node=%d\n", *which_node);
	// printf("here is the tree: %s\n", chaine);

	myNodes[*which_node].id=(*which_node);

	while (chaine[la_pos]!=')' && chaine[la_pos]!=';')
	{
		// fprintf(stderr,"%c",chaine[la_pos]);
		if(chaine[la_pos]==',')  /* beetwen 2 nodes */
		{
			ok_for_name = 1;

			which_desc = which_desc + 1;

			t= -1;
			myNodes[leaf].anc->descs = realloc(myNodes[leaf].anc->descs, sizeof (struct mymultinode) *(which_desc+1));
			myNodes[leaf].anc->nbdesc++;
			myNodes[leaf].anc->descs[which_desc] = &myNodes[(*which_node) +1];
			myNodes[(*which_node)+1].anc = myNodes[leaf].anc;

			myNodes[leaf].anc->descs[which_desc]->id= (*which_node)+1 ;
			(*which_node)++;
			/*printf("==>%d,%d\n",*which_node,which_desc);*/
			la_pos++;

		} /* end of -- if(chaine[la_pos]==',')  --*/
		else
		{

			if (chaine[la_pos]=='(') /* newnode */
			{
				ok_for_name = 1;
				if (myNodes[*which_node].descs == NULL) {
					if ((myNodes[*which_node].descs = malloc (sizeof(struct mymultinode)))==NULL) {
						fprintf(stderr, "Malloc error in reading nodes");
						exit(1);
					}
				}

 				myNodes[*which_node].descs[0] = &myNodes[(*which_node)+1];
				myNodes[*which_node].nbdesc++;
				myNodes[*which_node].evtype = (char)255;
				myNodes[(*which_node)+1].anc = &myNodes[*which_node];
				*which_node = *which_node + 1;
				/* printf("appel avec ==>%d,%d (%s)\n",*which_node,which_desc,chaine+la_pos+1); */
				la_pos=lit_newickmyMulti(chaine,la_pos+1,which_node,myNodes);	/* see recursively what's next */

			}
			else             /* all other chars which are not "," or "(" --> must be distance value or specie name */
			{

				if (chaine[la_pos]==':'){

					sscanf(chaine+la_pos,":%f",&t);
					if (t < 0)
						myNodes[leaf].anc->descs[which_desc]->time = 0;
					else
						myNodes[leaf].anc->descs[which_desc]->time = t;

				}
				else if (chaine[la_pos]=='['){
					/*
					{int k=0; char *ss=chaine,c;  while (*ss != '\0' && *ss != ']') {k++;}
					c=*(chaine+la_pos+k+1);*(chaine+la_pos+k+1)='\0';
					fprintf(stderr,"%s\n",chaine+la_pos); *(chaine+la_pos+k+1)=c;}
					*/
					if(	chaine[la_pos+1] != '&' )
						la_pos += set_evenement(chaine+la_pos,myNodes[leaf].anc->descs[which_desc], &fl); /* fonction lecture des evenements */
					else
						la_pos += readBEAST(chaine+la_pos,myNodes[leaf].anc->descs[which_desc]);     /* get info from the BEAST annotation */


					//printf("fin de chaine crochet %c",*(chaine+la_pos));
				}

				else if (t == -1) {                                          /* we must have a specie name somewhere */


					nom = countchar(chaine+la_pos);                    // count and retrieve the string


					if (! ok_for_name) {

						/*   Could well be a bootstrap or post. prob -- that is for now ignored -- */

						int x=0;
						while( x<strlen(nom) && (nom[x]>=48 && nom[x]<=57) || nom[x]==46 )x++;

						if(x<strlen(nom)){
							fprintf(stderr,"name found after a parenthesis. String=\n%s\nexiting...",chaine+la_pos);
							exit(1);
						}
						else
							free(nom);
					}
					else {

						if (myNodes[leaf].anc->descs[which_desc]->name != NULL) {
							free(myNodes[leaf].anc->descs[which_desc]->name);
						}
						myNodes[leaf].anc->descs[which_desc]->name =nom;
						/*printf("%s ;;;;;%s;%d,%d\n",nom,myNodes[leaf].anc->descs[which_desc]->name,leaf,which_desc);*/
						t = -2;                                          /* no more specie name */

						ok_for_name = 0;
						/* fprintf(stderr,"ok_for_name=:%d\n",ok_for_name); */
					}


				}

				la_pos++;
				// printf("char=%c\n", chaine[la_pos]);
				// printf("next char = %c\n", chaine[la_pos+1]);
				// printf("second next char = %c\n", chaine[la_pos+2]);
			}
		}
	}

	return(la_pos+1);
}
Node *ReadTreeFile( char *infile, int *nleaves ){


	FILE *f;
	Node *InTree;
	char *newick_string;
	int bb=0;

	f=fopen( infile, "r" );

	if(!f)fprintf(stderr, "main: cannot open file %s, please check the file, bye\n", infile);

	newick_string = readFileNwck(f, nleaves);     /* Lit un fichier Newick, compte le nombre de feuilles et le retourne en chaine de caractere */

	fclose(f);

	InTree = InitMultiTree(  *nleaves );           /* Initialise un arbre de la taille trouvée précédement et en mettant tous les parametres par defaut pour ensuite les remplir*/

	lit_newickmyMulti( newick_string, 0, &bb, InTree );         /* ?? Lit une chaine de caractere qui est le nom des especes et le met dans un tableau de structure ?? */


	return InTree;
}

/*
	If outfile is NULL, print out in stdout
*/
void PrintTree( Node *Tree, char *outfile, short opt_FromRoot ){

	FILE *f=stdout;


	if(outfile){
		f= fopen(outfile, "w");
		if(!f)fprintf(stderr, "cannot open file %s for output, please check the directory, bye\n", outfile), exit(2);
	}

	if(opt_FromRoot)
		multitreeoutNck( FindRoot( Tree ), f );
	else
		multitreeoutNck( Tree, f );

	fprintf(f,";\n");

	if(outfile)
		fclose(f);

}
