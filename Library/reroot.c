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
\brief 	everything for rerooting a tree
\author g achaz and madamesophie
\date 	June 28 2012
**/


#include "tree.h"



/*
	Reverse ancestrality and descendance from a given node up to the root
*/
void ReverseAncLink( Node * Target ){

	Node *tmp;
	int i;

	if( Target==NULL || Target->anc == NULL )
		return;

//	printf("now swap node %d with its ancestor\n", Target->id);

	ReverseAncLink( Target->anc );

	tmp=Target->anc;

	Connect( Target , Target->anc , Target->time );   // exchange anc and desc status

	Target->anc=NULL;
	for(i=0;i<tmp->nbdesc;i++)
		if( tmp->descs[i]==tmp->anc )
		{
			tmp->descs[i] = tmp->descs[tmp->nbdesc-1];
			tmp->nbdesc--;

		}

//	PrintNode( Target, "Target" );
//	PrintNode( tmp, "Tmp" );

}


Node *FindFirstFreeNode( Node *p  ){

	while ( p->anc != NULL || p->descs != NULL )
		p++;

	return p;

}

int FindMaxIdSubTree( Node *p ){

	int id=p->id, i;

	for(i=0;i<p->nbdesc; i++)
		if( FindMaxIdSubTree( p->descs[i] ) > id )
			id = FindMaxIdSubTree( p->descs[i] );

	return id;
}

int FindMaxId( Node *p ){
	return FindMaxIdSubTree( FindRoot(p) );
}

void ReRootMulti( Node * TargetNode ){

	Node *pN, *pNFree;

	pNFree = FindFirstFreeNode( TargetNode  );
	pNFree->id = FindMaxId( TargetNode )+1;

	pN     = TargetNode->anc;

	Connect(  pNFree, TargetNode, TargetNode->time/2.0 );
	Connect( pN, pNFree, TargetNode->time/2.0 );

	ReverseAncLink( pNFree );

}


void ReRoot( Node * TargetNode ){


	Node *RootNode, *pN;

	int *tmp=NULL;

	/*
		First thing, remove the rootand rewire the nodes under the root adequatly
	*/
	RootNode = FindRoot( TargetNode );     /* Find the root of the tree */


	if( RootNode->nbdesc > 2 ){

		ReRootMulti( TargetNode );
		return;
	}


	if( TargetNode == RootNode || (TargetNode->anc == RootNode) )   /* in that case, do nothing */
		return;


	/*
		Find a node under root that is not a leaf and set it as 'pN'
	*/
	pN=RootNode->descs[0];

	while( pN->nbdesc < 2 && pN-RootNode->descs[0] < RootNode->nbdesc )
		pN++;

	if(pN-RootNode->descs[0] >= RootNode->nbdesc)
		fprintf(stderr, "ReRoot: The tree has only leaves, cannot reroot it, sorry\n"), exit(5);


	/*
		Deal with the Root Node
		connect the desc != pN to node pN as an ancestor
		desc[0] is a temporary 'root' TempRoot
	*/
	pN->anc = NULL;

	int x=( RootNode->descs[0] == pN )?1:0;
	Connect( pN, RootNode->descs[x], pN->time+ RootNode->descs[x]->time );


	/*
		Then empty Root Node
	*/
	tmp = RootNode->SeqInLeaves;
	RootNode->SeqInLeaves = NULL;
	FreeNode( RootNode );
	RootNode->SeqInLeaves = tmp;


	/*
		Rewire all nodes from TargetNode->anc to the temporary 'root' TempRoot
		in opposite anc/desc connections
	*/
	ReverseAncLink( TargetNode->anc );      /* at this point, there is no root ! */

	/*
		Insert Root above TargetNode
	*/

	Connect( RootNode, TargetNode->anc, TargetNode->time/2.0 );
	Connect( RootNode, TargetNode, TargetNode->time/2.0 );

}

/*
	Reroot with outgroup
*/
void ReRootOutgroup( char *outgroup, Node *init_tree, int nnodes ){

	int *list_node,
	     size_list;

	Node *pNode=NULL, *pMRCA;


//	print_arbre(FindRoot(init_tree), stdout);


	/*
		Check whether outgroup is monophyletic
		use opt_inv = 0 [ 4th argument ]
	*/
 	list_node = GenerateLeavesArray( outgroup, init_tree, nnodes, 0, &size_list);
	pMRCA = FindMRCA_intarray( list_node, size_list, init_tree, nnodes );
	free(list_node);

	if( count_leaf(pMRCA) == size_list )
		pNode = pMRCA;
	else
	{

		/*
			Otherwise check whether the in group is monophyletic
			use opt_inv = 1 [ 4th argument ]
		*/
		list_node = GenerateLeavesArray( outgroup, init_tree, nnodes, 1, &size_list);
		pMRCA = FindMRCA_intarray( list_node, size_list, init_tree, nnodes );
		free(list_node);

		if( count_leaf(pMRCA) == size_list )
			pNode = pMRCA;
	}

	if( pNode )
		ReRoot( pNode );
	else
		fprintf(stderr, "either the ingroup or the outgroup is monophyletic. Cannot reroot. Bye\n"),exit(1);

//	print_arbre(FindRoot(init_tree), stdout);

}


void ReLabelNodeId( Node *node, int *id ){

	int i;
	node->id = *id;
	for(i=0;i<node->nbdesc;i++)
	{
		(*id)++;
		ReLabelNodeId( node->descs[i], id );
	}

}
Node *SetRootEpics( Node *tree, char *outgroup, int output_tree){

	Node *pMRCA, *Root=NULL;

	int ingroup_size,
	    *ingroup_lst;

//	printf("look for outgroup %s\n",outgroup);

	int nnodes = count_nodes(FindRoot(tree));

	if( outgroup )
	{

		ingroup_lst = GenerateLeavesArray( outgroup, tree, nnodes, 1, &ingroup_size);       // retrieve all nodes, except the outgroup(s)
		pMRCA = FindMRCA_intarray( ingroup_lst, ingroup_size, tree, nnodes );

		free(ingroup_lst);

		if( count_leaf(pMRCA) != ingroup_size )
			fprintf(stderr, "SetRootEpics: the rooted tree is paraphyletic (outgroup excluded).\
			    ingroup_size: %d, #leaves: %d\nCannot proceed further, bye\n", ingroup_size, count_leaf(pMRCA)), exit(1);

		Root = pMRCA;

	}
	else
		Root=FindRoot(tree);

	if (output_tree )
		PrintTree(Root, NULL , 0 );


//	fprintf(stderr,"Numbers of species in the tree : %d (outgroup excluded)\n", count_leaf(Root));


	/*
		Reset the node as the root
		no length, no anc, no evts, with id=0
	*/
	int xx=0;
	ReLabelNodeId( Root, &xx );
	Root->anc=NULL;
	Root->nevt=0;
	free(Root->evt);
	Root->evt=NULL;
	Root->evtype=0;
	Root->time=0;

	return Root;

}
