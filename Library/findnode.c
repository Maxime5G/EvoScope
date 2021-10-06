/*
	Copyright (C) 2002-2013 G Achaz & S Brouillet

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
\brief 	find specific node in the tree
\author g achaz and madamesophie
\date 	June 18 2013
**/

#ifdef LINUX
#include <strings.h>
#endif

#include <string.h>

#include "tree.h"

/**
\fn    Node * FindRoot(   Node *mynode )
\brief Recurrsive Search BottomUp --faster than findHead() --
**/
Node * FindRoot( Node * mynode ) 
{
	return ( mynode->anc==NULL) ? mynode : FindRoot( mynode->anc );
}

/**
\fn	struct  mymultinode *findHead(struct mymultinode *myNode,  int nnodes )
\brief find the root of the tree in case not the 1st elt
**/
struct  mymultinode *findHead(struct mymultinode *myNode,  int nnodes )            /* Permet de trouver la racine de l'arbre*/
{
	int i;
	for (i=0;i<nnodes;i++)
		if (myNode[i].anc==NULL)
			return myNode+i;
			
	fprintf(stderr,"Error no root was found\n"),exit(1);		
}

/*---------------------------------------------------*/


/*
	From a char name of a node, finds it into the Node array (ie the tree)
*/
Node * FindLeaf(  char *name, Node *NodeArray, int nnodes )
{

	int i=0;

	while( i < nnodes ){
	
		if( NodeArray[i].name != NULL && strcmp(NodeArray[i].name, name) == 0 ){
			return NodeArray+i;
		}
		i++;
	}
	
	return NULL;
}



/*
	From two nodes, find their MRCA
*/
Node * FindMRCA(  Node *n1, Node *n2  )
{

	if(n1 == NULL || n2 == NULL)
		return NULL;

	if( isAnc( n1, n2 ) )
		return n2;
	
	while( isAnc( n2, n1 ) == 0 && n1->anc != NULL)
		n1 = n1->anc;

	return n1;
}


/*
	From 1 to x nodes given as char *name, find their MRCA
	request is either "node_name" or "node_name1:node_name2:node_name3:..."
*/
Node * FindMRCA_char( char *request, Node *NodeArray, int nnodes ){


	char *p=request,
	    *mycopy,
	     *q;
	Node *pn;
	
	
	mycopy=(char *)malloc(strlen(request)*sizeof(char));
	strcpy(mycopy,request);
	p=mycopy;
			
	/*
		Scan first leaf
	*/
	p=index(mycopy, ':' );
	if(p)*p=0;

	pn = FindLeaf(mycopy,NodeArray, nnodes);
	if(p)*p=':';
		
	/*
		And all subsequent ones and iteratively the MRCA's
	*/
	while (p)
	{
		q=p+1;
		p = index(q, ':' );
		if(p)*p=0;
		pn = FindMRCA( pn  ,  FindLeaf(q, NodeArray, nnodes)  );
		if(p)*p=':';
	}
	
	free(mycopy);
	
	return pn;

}

/*
	From 1 to x nodes given as int *array, indexed by their node number
	the array has size 'array_size'
*/
Node *FindMRCA_intarray( int *array, int array_size, Node *NodeArray, int nnodes ){

	Node *pn;
	int i=0;
	
	pn = NodeArray+array[i];
	
//	printf("pn = %d (name=%s)\n", pn->id, pn->name);
		
	/*
		Find iteratively the MRCA's
	*/
	while (i+1<array_size)
	{
		pn = FindMRCA( pn  ,  NodeArray+array[i+1]  );
//		printf("pn = %d (name=%s)\n", pn->id, pn->name);
		
		i++;
	}
	
//	printf("<> pn = %d (name=%s) %d leaves\n", pn->id, pn->name, count_leaf(pn));
	
	return pn;

}

