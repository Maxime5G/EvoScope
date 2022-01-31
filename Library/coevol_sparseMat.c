/* ------------------------------------------------------------------------	*/
/* Generation and operations on sparse matrices							    */
/* ------------------------------------------------------------------------	*/

#include <stdio.h>
#include <stdlib.h>

#ifdef LINUX
#include <strings.h>
#endif
#include <string.h>

#include <math.h>
#include "tree.h"
#include "coevol.h"

// extern int verbose;

/* ------------------------------------------------------------------------------------------ */
/* Convert an existing matrix (struct IntegerMatrix) into sparse format.    				  */
/* The sparse format follows the Yale convention and is implemented in the  				  */
/* Compressed Sparse Column (CSC) format.								    				  */
/* The choice for the CSC format has been made in order to be able to 	    				  */
/* perform easily the e1T.M operations.									    				  */
/* how does it work (for a CSC format)?														  */
/* for a given matrix, a sparse format contains three vectors:			    				  */
/* the vector V, containing the non-zero (nnz) values of a matrix.		    				  */
/* the vector ROW (length=nnz), containing the rows indices of the nnz values.		    	  */
/* the vector COL (length=ncol+1), has one element per col in the matrix, and encodes   	  */
/* the indices in V where the given column starts (Definition from Wikipedia).		    	  */
/* The last element in COL is nnz, by convention.									    	  */
/* In the current version, V is not implemented, as we work only with binary matrices. 		  */
/* ------------------------------------------------------------------------------------------ */

sparseMatrix *convertToSparseMatrix(int nrow, IntegerMatrix *input_mat){
	sparseMatrix *m;
	int i,j,k, cnt_col, cnt_row, check_col, check_row_val, nnz;
	/* i = rows, j = cols, k = values in input_mat.															 */
	/* cnt_row = indices of non-zero elements to store in the row vector of m.								 */
	/* cnt_col = indices corresponding to the index of the non-zero elements where the given column starts.  */

	if ((m=malloc(sizeof(sparseMatrix))) == NULL){
		fprintf(stderr,"Not enough memory for matrix allocation, exiting\n");
		exit(1);
	}

	/* non-zeroes elements of the matrix. Since elements are only 1 or 0, it is the cumulative sum of the elements. */
	// nnz=sumv(input_mat->val, input_mat->nrow*input_mat->ncol);

	/* in the case of masking, negative values have to be recorder as well: nnzv does that */
	nnz=nnzv(input_mat->val, input_mat->nrow*input_mat->ncol);

	m->row=calloc(nnz, sizeof(int)*nnz);				/* the ROW vector 						   */
	m->col=calloc(nrow+1, sizeof(int)*nrow+1);	/* the COL vector, initialized with zeroes */

	cnt_row = cnt_col = 0;						/* setting the indices to fill COL and ROW */

	for (j = 0; j < input_mat->ncol; j++) {									/* for each column in the matrix 							*/
		check_row_val=0;
		for (i = 0, k=j; i < input_mat->nrow; i++,k+=input_mat->ncol) { 	/* and for each row (remember, for a given i,j, k=ncol*i+j) */
			if (input_mat->val[k]){											/* if I have a non-zero value 								*/
				check_row_val=1;
				m->row[cnt_row++] = i;										/* I store the index of the row in ROW 						*/
				if (!(check_col)){												/* if I don't yet have a value for the COL 					*/
					check_col=1;												/* I mark the flag as done 									*/
					m->col[j]=cnt_col;										/* I add to the COL vector the current column considered	*/
				}
				cnt_col++;													/* and I move forward of one for each nnz value encountered */
			}
		}
		if (!check_row_val)
			m->col[j]=cnt_col;
		check_col=0;															/* after having done the j column - reinit the check 		*/
	}
	m->col[j]=cnt_row;														/* the last element of COL is nnz. 							*/
	return m;
}

/* ------------------------------------------------------------------------------------------ */
/* Performs the multiplication of a vector with a sparse matrix.							  */
/* Essentially, I'm retrieving the indices of the nnz values in the sparse and adding 	      */
/* together the corresponding values in the vector (since we have only 1's as nnz).			  */
/* ------------------------------------------------------------------------------------------ */

void vectSparseMat (int *v, sparseMatrix *m, int lv, int *pv){
	int i, j; /* i=columns, j=rows */

	/* lv is the length of the vector, which is also the	*/
	/* number of rows and columns of the matrix 			*/

	for (i=0; i<lv; i++){						/* For each column/row 											 */
		pv[i]=0;
		for (j=m->col[i]; j<m->col[i+1]; j++){	/* The indices of the nnz elements are retrieved by splicing COL */
			pv[i]+=v[m->row[j]];				/* The values corresponding in the vectors are summed together   */
		}
	}
	return;
}

void vectSparseMatAndMask (int *v, sparseMatrix *m, int lv, int *pv, char *maske1, char *maske2){
	int i, j; /* i=columns, j=rows */

	/* lv is the length of the vector, which is also the	*/
	/* number of rows and columns of the matrix 			*/

	for (i=0; i<lv; i++){						/* For each column/row 											 */
		pv[i]=0;
		for (j=m->col[i]; j<m->col[i+1]; j++){	/* The indices of the nnz elements are retrieved by splicing COL */

			#if 0
			printf("index = %d value = %d mask value = %d\n", m->row[j], v[m->row[j]], mask[m->row[j]]);
			#endif

			pv[i]+=v[m->row[j]] * (int)maske1[m->row[j]] * (int)maske2[m->row[j]];				/* The values corresponding in the vectors are summed together   */
		}
	}
	return;
}
