#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "tree.h"
#include "coevol.h"
#define EPSILON 0.001
#include <time.h>

double Gauss_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	return exp(- (X[0]*X[0] + X[1]*X[1] + X[2]*X[2]) );

}

double param1_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0],X[0],X[0]}, {X[0],X[0],X[0],X[0]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double indep_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0],X[0],X[0]}, {X[1],X[1],X[1],X[1]} };

	for (int t=0; t<ntree; t++){
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	}
	return res;
}

double indep_kappa_uniq_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0],X[1],X[1]}, {X[2],X[2],X[2],X[2]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double indep_kappa_uniq2_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0],X[0],X[0]}, {X[1],X[1],X[2],X[2]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}


double indep_kappa_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0],X[1],X[1]}, {X[2],X[2],X[3],X[3]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double param2_induc_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0]*X[1],X[0],X[0]*X[1]}, {X[0],X[0]*X[1],X[0],X[0]*X[1]} };
	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}


double asymetric_1on2_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0],X[0],X[0]}, {X[1],X[2],X[1],X[2]} };
	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double asymetric_2on1_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[1],X[0],X[1]}, {X[2],X[2],X[2],X[2]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}


double induc_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[1],X[0],X[1]}, {X[2],X[3],X[2],X[3]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double induc_2k_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[1],X[0]*X[4],X[1]*X[4]}, {X[2],X[3],X[2]*X[5],X[3]*X[5]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}


double symetric_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0]*X[2],X[0],X[0]*X[2]}, {X[1],X[1]*X[2],X[1],X[1]*X[2]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double symetric_1k_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){
	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0]*X[2],X[0]*X[3],X[0]*X[3]*X[2] }, {X[1],X[1]*X[2],X[1]*X[3],X[1]*X[3]*X[2]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double mu_lambda_kappa_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0]*X[1],X[0]*X[2],X[0]*X[1]*X[2]}, {X[0],X[0]*X[1],X[0]*X[2],X[0]*X[1]*X[2]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double asymetric_1on2_kappa1_kappa2_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[0],X[0]*X[3],X[0]*X[3]}, {X[1],X[2],X[1]*X[4],X[2]*X[4]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double asymetric_2on1_kappa1_kappa2_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree ){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[1],X[0]*X[3],X[1]*X[3]}, {X[2],X[2],X[2]*X[4],X[2]*X[4]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double allfree_multi( double *X, struct CoevolData *MyEpocsData, int IS[], int ntree){

	double res=0.0;
	double temp_rate[2][4]={ {X[0],X[1],X[2],X[3]}, {X[4],X[5],X[6],X[7]} };

	for (int t=0; t<ntree; t++)
		res+=ln_vraisemblance(MyEpocsData[t].root, *MyEpocsData[t].theTVector.current_vector, temp_rate, IS );
	return res;
}

double ( * pick_scenario_multi (char scenario, int *d ) )(double *, struct CoevolData *, int[], int){

		switch(scenario){

			case '1':
				*d=1;
				return &param1_multi;

			case 'i':
				*d=2;
				return &indep_multi;

			case 'x':
				*d=2;
				return &param2_induc_multi;

			case 'X':
				*d=3;
				return &mu_lambda_kappa_multi;

			case 'I':
				*d=4;
				return &indep_kappa_multi;

			case 'a':
				*d=3;
				return &asymetric_1on2_multi;

			case 'A':
				*d=5;
				return &asymetric_1on2_kappa1_kappa2_multi;

			case 'b':
				*d=3;
				return &asymetric_2on1_multi;

			case 'B':
				*d=5;
				return &asymetric_2on1_kappa1_kappa2_multi;

			case 'l':
				*d=4;
				return &induc_multi;

			case 'L':
				*d=6;
				return &induc_2k_multi;

			case 'r':
				*d=3;
				return &symetric_multi;

			case 'R':
				*d=4;
				return &symetric_1k_multi;

			case '8':
				*d=8;
				return &allfree_multi;

			case 'g':
				*d=3;
				return &Gauss_multi;

			case 'u':
				*d=3;
				return &indep_kappa_uniq_multi;

			case 'U':
				*d=3;
				return &indep_kappa_uniq2_multi;


			default:
				fprintf(stderr, "Sorry! Scenario '%c' is not implemented (yet?), bye\n", scenario), exit(4);
		}

	}

void set_rates_scenario( char scenario, double rate[2][4], double *X ){

	switch(scenario){

		case '1':
			rate[0][0]=rate[0][1]=rate[0][2]=rate[0][3]=X[0];
			rate[1][0]=rate[1][1]=rate[1][2]=rate[1][3]=X[0];
			break;

		case 'i':
			rate[0][0]=rate[0][1]=rate[0][2]=rate[0][3]=X[0];
			rate[1][0]=rate[1][1]=rate[1][2]=rate[1][3]=X[1];
			break;

		case 'x':
			rate[0][0]=rate[0][2]=rate[1][0]=rate[1][2]=X[0];
			rate[0][1]=rate[0][3]=rate[1][1]=rate[1][3]=X[0]*X[1];
			break;

		case 'X':
			rate[0][0]=rate[1][0]=X[0];
			rate[0][1]=rate[1][1]=X[0]*X[1];
			rate[0][2]=rate[1][2]=X[0]*X[2];
			rate[0][3]=rate[1][3]=X[0]*X[1]*X[2];
			break;

		case 'I':
			rate[0][0]=rate[0][1]=X[0];
			rate[0][2]=rate[0][3]=X[1];
			rate[1][0]=rate[1][1]=X[2];
			rate[1][2]=rate[1][3]=X[3];
			break;

		case 'l':
			rate[0][0]=rate[0][2]=X[0];
			rate[0][1]=rate[0][3]=X[1];
			rate[1][0]=rate[1][2]=X[2];
			rate[1][1]=rate[1][3]=X[3];
			break;

		case 'L':
			rate[0][0]=X[0];
			rate[0][1]=X[1];
			rate[0][2]=X[0]*X[4];
			rate[0][3]=X[1]*X[4];
			rate[1][0]=X[2];
			rate[1][1]=X[3];
			rate[1][2]=X[2]*X[5];
			rate[1][3]=X[3]*X[5];
			break;


		case 'r':
			rate[0][0] = X[0];
			rate[0][1] = X[0]*X[2];
			rate[0][2] = X[0];
			rate[0][3] = X[0]*X[2];
			rate[1][0] = X[1];
			rate[1][1] = X[1]*X[2];
			rate[1][2] = X[1];
			rate[1][3] = X[1]*X[2];
			break;

		case 'R':
			rate[0][0]=X[0];
			rate[0][1]=X[0]*X[2];
			rate[0][2]=X[0]*X[3];
			rate[0][3]=X[0]*X[2]*X[3];
			rate[1][0]=X[1];
			rate[1][1]=X[1]*X[2];
			rate[1][2]=X[1]*X[3];
			rate[1][3]=X[1]*X[2]*X[3];
			break;

		case 'a':
			rate[0][0]=X[0];
			rate[0][1]=X[0];
			rate[0][2]=X[0];
			rate[0][3]=X[0];
			rate[1][0]=X[1];
			rate[1][1]=X[2];
			rate[1][2]=X[1];
			rate[1][3]=X[2];
			break;

		case 'A':
			rate[0][0]=X[0];
			rate[0][1]=X[0];
			rate[0][2]=X[0]*X[3];
			rate[0][3]=X[0]*X[3];
			rate[1][0]=X[1];
			rate[1][1]=X[2];
			rate[1][2]=X[1]*X[4];
			rate[1][3]=X[2]*X[4];
			break;

		case 'b':
			rate[0][0]=X[0];
			rate[0][1]=X[1];
			rate[0][2]=X[0];
			rate[0][3]=X[1];
			rate[1][0]=X[2];
			rate[1][1]=X[2];
			rate[1][2]=X[2];
			rate[1][3]=X[2];
			break;

		case 'B':
			rate[0][0]=X[0];
			rate[0][1]=X[1];
			rate[0][2]=X[0]*X[3];
			rate[0][3]=X[1]*X[3];
			rate[1][0]=X[2];
			rate[1][1]=X[2];
			rate[1][2]=X[2]*X[4];
			rate[1][3]=X[2]*X[4];
			break;


		case '8':
			rate[0][0]=X[0];
			rate[0][1]=X[1];
			rate[0][2]=X[2];
			rate[0][3]=X[3];
			rate[1][0]=X[4];
			rate[1][1]=X[5];
			rate[1][2]=X[6];
			rate[1][3]=X[7];
			break;

		case 'u':
			rate[0][0]=X[0];
			rate[0][1]=X[0];
			rate[0][2]=X[1];
			rate[0][3]=X[1];
			rate[1][0]=X[2];
			rate[1][1]=X[2];
			rate[1][2]=X[2];
			rate[1][3]=X[2];
			break;

		case 'U':
			rate[0][0]=X[0];
			rate[0][1]=X[0];
			rate[0][2]=X[0];
			rate[0][3]=X[0];
			rate[1][0]=X[1];
			rate[1][1]=X[1];
			rate[1][2]=X[2];
			rate[1][3]=X[2];
			break;


		default:
			fprintf(stderr, "Sorry! Scenario '%c' is not implemented (yet?), bye\n", scenario), exit(4);
	}

}

/*
	Debugging
*/
void PrintVectMask( char *title, double *Vect, int d, char *mask){

	int i, skip=0;

//	printf("\n");
	printf("   %10s:",title);
//	printf("\t   ");
	for(i=0;i<d;i++)
		if(mask[i]){
			printf("              -");
			skip++;
		}
		else
			printf("%15.7g",Vect[i-skip]);
	printf("\n");


}
void PrintVect( char *title, double *Vect, int d){

	int i;

//	printf("\n");
	printf("   %10s:",title);
//	printf("\t   ");
	for(i=0;i<d;i++)
		printf("%15.7g",Vect[i]);
//		printf("%20.17g",Vect[i]);
	printf("\n");


}
void PrintStr( char *title, char *Str, int d){

	int i;

	printf("   %10s:",title);
	for(i=0;i<d;i++)
		printf("%15d",(int)Str[i]);
	printf("\n");


}


void PrintMat( char *title, double **Mat, int d ){

	int i,j;

//	printf("\n");
	printf("   %10s:\n",title);

	for(i=0;i<d;i++){
		printf("              ");
		for(j=0;j<d;j++)
			printf("%15.4g",Mat[i][j]);
			printf("\n");
	}

}





/*
   Recursive definition of determinate using expansion by minors.
*/
double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}

/*
   Find the cofactor matrix of a square matrix
*/
void CoFactor(double **a,int n,double **b)
{
   int i,j,ii,jj,i1,j1;
   double det;
   double **c;

   c = malloc((n-1)*sizeof(double *));
   for (i=0;i<n-1;i++)
     c[i] = malloc((n-1)*sizeof(double));

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[i1][j1] = a[ii][jj];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         det = Determinant(c,n-1);

         /* Fill in the elements of the cofactor */
         b[i][j] = pow(-1.0,i+j+2.0) * det;
      }
   }
   for (i=0;i<n-1;i++)
      free(c[i]);
   free(c);
}



void invert_H( double **H, double **Hinv, int d){

	int i,j;
	double det;

	if(d==1){
		Hinv[0][0]=1.0/H[0][0];
		return;
	}

	det=Determinant(H, d);
	CoFactor(H,d,Hinv);

	for(i=0;i<d;i++)
		for(j=0;j<d;j++)
			Hinv[i][j] /= det;

}

void prod_M_V( double **M, double *V, double *R, int d){

	int i,j;

	for(i=0;i<d;i++)
	{
		R[i]=0;
		for(j=0;j<d;j++)
			R[i]+=V[j]*M[i][j];

	}
}

char sign( double X ){
	return (X<0)?-1:( (X==0)?0:1);
}

/*
	First check if all dimensions are optimized in a direction coherent with g
int check_dX( double *dX, int d, double *g, char *mask ){

	int i;
	int skip=0;

	for(i=0;i<d;i++)
		if( mask[i] == 1 )
			skip++;
		else
			if( sign(g[i-skip]) != sign(-dX[i-skip]) && fabs(g[i-skip]) >= MIN_GRAD )
				return 0;

	return 1;

}
*/

/*
	This function
*/
void update_X( double *X, double *dX, int d, double *g, char *mask){

	int i;
	int skip=0;


	/*
		Update the X vector
	*/
	for(i=0;i<d;i++){

		if( mask[i] == 1 )
			skip++;
		else{

/*			if( 0 && sign(g[i-skip]) != sign(-dX[i-skip]) && fabs(dX[i-skip])>X[i]/2.0 ){
				fprintf(stderr, "update_X: movement in wrong direction for %d\n",i);
				fprintf(stderr, "[%d] X:%f g: %f, dX: %f\n",i, X[i], g[i-skip], -dX[i-skip]);
				exit(1);
			}
			else
			{*/
				X[i] -= dX[i-skip];
				if(X[i]<0){X[i]=MIN_MU;}
				if(X[i]>MAX_MU){X[i]=MAX_MU; dX[i-skip]=0;}
			//}
		}

	}

}


double dist_U_V( double *U, double *V, int d){

	int i;
	double dist=0;

	for(i=0;i<d;i++)
		dist += (U[i]-V[i])*(U[i]-V[i]);

	return dist;
}

double norm_vect( double *V, int d, char *mask){

	int i;
	double dist=0;

	if(mask == NULL)               // test only once for mask vector
		for(i=0;i<d;i++)
			dist += V[i]*V[i];
	else
		for(i=0;i<d;i++)
			if(mask[i] != 1)
				dist += V[i]*V[i];

	return sqrt(dist);
}

double factor( double X ){
	double factor=1.5;
	return (X<0)?(1.0/factor):( (X==0)?1:factor);
}

double grad_in_xi_multi( double *X, int i, struct CoevolData *MyEpocsData, int IS[], double (*scenar_multi)(double *, struct CoevolData *, int[], int), int ntree ){

	double Xi, e;
	double g;
	char ne=0;

	Xi=X[i];
	e=EPSILON*Xi;   /* e is always at the "good" scale */


	if(Xi != MAX_MU) { X[i] = Xi+e; ne++;}

	g =  (*scenar_multi)( X, MyEpocsData, IS, ntree );

	if(Xi != MIN_MU) { X[i] = Xi-e; ne++; }
	else             X[i] = Xi;

	g -= (*scenar_multi)( X, MyEpocsData, IS, ntree );

	X[i] = Xi;

	g /= ne*e;

	return g;
}


/*
	When mask is set to NULL, only compute g
	Otherwise, check mask and eventually update mask
               when g is too small or X is at max/min
*/
int compute_g_multi( double *g, int d, double *X, char *mask, char opt_mask_update, struct CoevolData *MyEpocsData, int IS[], double (*scenar_multi)(double *, struct CoevolData *, int[], int), int ntree ){

	int i;
	int skip=0;        /* the number of values that are masked */

	double gxi;

	skip=0;
	for(i=0; i<d; i++)
	{

		if( mask && mask[i]==1 ){
			skip++;
			continue;
		}

		gxi = grad_in_xi_multi( X, i, MyEpocsData, IS, scenar_multi, ntree);

		// max values for gradient
		if( gxi>100  )  gxi= 100;
		if( gxi<-100 )  gxi=-100;

		if( opt_mask_update ){

			if( (fabs(gxi) < MIN_GRAD ) || \
			    (X[i] <= MIN_MU && gxi<0) || \
			    (X[i] >= MAX_MU && gxi>0) )
			{
				if(verbose>1)
					printf("---> dim[%d]: gxi=%e Xi=%f min=%e max=%e --> freeze dim[%d]\n", \
					        i, gxi, X[i], MIN_MU, MAX_MU, i);

				mask[i] = 1;
				skip++;
				continue;
			}
		}

		g[i-skip] = gxi;


	}

	return d-skip;
}

//
// !! Here, H is not fully correct when we reach the MIN_MU or MAX_MU !!
//
int compute_H_multi(  double ** H, int d, double *X, char *mask, struct CoevolData *MyEpocsData, int IS[], double (*scenar_multi)(double *, struct CoevolData *, int[], int), int ntree ){

	int i,j;
	double e, ej;
	double f0 = (*scenar_multi)( X, MyEpocsData, IS, ntree );

	int skip=0;   /* the number of values that are masked */

	int skip_j=0;

	double Xi, Xj;

	skip=0;
	for(i=0; i<d; i++)
	{

		if(mask[i]==1)
		{
			skip++;
			continue;
		}

		Xi=X[i];
		e=EPSILON*Xi;   /* e is always at the "good" scale */


		/*  H(i,i) numerator */
		X[i] = Xi+e;
		H[i-skip][i-skip] =  (*scenar_multi)( X, MyEpocsData, IS, ntree );

		X[i] = Xi-e;
		H[i-skip][i-skip] += (*scenar_multi)( X, MyEpocsData, IS, ntree );

		X[i] = Xi;

		/* H(i,i) */
		H[i-skip][i-skip] = (H[i-skip][i-skip]-2.0*f0)/(e*e);


		/* H(i,j)=H(j,i) */
		skip_j=0;

		for(j=0;j<i;j++)
			if(mask[j] == 1)
				skip_j++;

		for(j=i+1;j<d;j++)
		{

			if(mask[j] == 1){
				skip_j++;
				continue;
			}

			Xi=X[i];
			Xj=X[j];
			ej=EPSILON*Xj;

			X[i] = Xi+e;
			X[j] = Xj+ej;
			H[i-skip][j-skip_j] = (*scenar_multi)( X, MyEpocsData, IS, ntree );

			X[i] = Xi-e;
			X[j] = Xj-ej;
			H[i-skip][j-skip_j] += (*scenar_multi)( X, MyEpocsData, IS, ntree );

			X[i] = Xi+e;
			X[j] = Xj-ej;
			H[i-skip][j-skip_j] -= (*scenar_multi)( X, MyEpocsData, IS, ntree );

			X[i] = Xi-e;
			X[j] = Xj+ej;
			H[i-skip][j-skip_j] -= (*scenar_multi)( X, MyEpocsData, IS, ntree );

			X[i] = Xi;
			X[j] = Xj;
			H[i-skip][j-skip_j] /= (4*e*ej);

			H[j-skip_j][i-skip] = H[i-skip][j-skip_j];

		}

	}

	return d-skip;

}




double **MemMat( int d ){

	int i;
	double **Mat;

	Mat = (double **)malloc( d*sizeof(double *) );
	if(!Mat)fprintf(stderr, "MemMat: cannot allocate Mat, sorry, bye"), exit(3);
	for(i=0; i<d; i++){
		Mat[i] = (double *)calloc( d, sizeof(double) );
		if(!Mat[i])fprintf(stderr, "MemMat: cannot allocate Mat[%d], sorry, bye",i), exit(3);
	}

	return Mat;

}

double *MemVect( int d ){

	double *V;

	V = (double *)malloc( d*sizeof(double) );
	if(!V)fprintf(stderr, "MemVect: cannot allocate V, sorry, bye"), exit(3);

	return V;
}

void MemSet( double **g, double **dX, double ***H, double ***Hinv, char **mask, int d ){

	*g  =   MemVect( d );
	*dX =   MemVect( d );
	*H    = MemMat( d );
	*Hinv = MemMat( d );

	*mask = (char*)calloc( d, sizeof(char) );
	if(!*mask)fprintf(stderr, "MemSet: cannot allocate mask, bye\n"),exit(3);

}
void FreeSet( double *g, double *dX, double **H, double **Hinv, char *mask, int d ){

	int i;

	free(g);
	free(dX);

	for(i=0;i<d;i++)free(H[i]);
	free(H);

	for(i=0;i<d;i++)free(Hinv[i]);
	free(Hinv);

	free(mask);
}



int build_mask_multi( double *g, int d,  double *X, char *mask, struct CoevolData *MyEpocsData, int IS[], double (*scenar_multi)(double *, struct CoevolData *, int[], int), int ntree ){

	int i;
	double Xi;
	int d2=d;

	for(i=0;i<d;i++)
	{
		mask[i]=0;

		if(g[i] == 0){    /* the value is not used in the calculation */
			X[i] = -1;
			mask[i]=1;
			d2--;
		}

		if( g[i] < 0){       /* check it ML is at 0 */

				double f0,fe;

				Xi = X[i];
				X[i]=0;
				f0 = (*scenar_multi)( X, MyEpocsData, IS, ntree );

				X[i]=MIN_MU;
				fe = (*scenar_multi)( X, MyEpocsData, IS, ntree );

				if( f0>=fe ){
					X[i]=0;
					mask[i]=1;
					d2--;
				}
				else
					X[i]=Xi;
		}

		if( g[i] > 0){       /* check if ML is at MAX_MU */

				double f0,fe;

				Xi = X[i];
				X[i]=MAX_MU;
				f0 = (*scenar_multi)( X, MyEpocsData, IS, ntree );

				X[i]=MAX_MU-MIN_MU;
				fe = (*scenar_multi)( X, MyEpocsData, IS, ntree );

				if( f0>=fe ){
					X[i]=MAX_MU;
					mask[i]=1;
					d2--;
				}
				else
					X[i]=Xi;
		}

	}

	return d2;
}


double check_limit_rate( double xi ){
	if( xi > MAX_MU ) return MAX_MU;
	if( xi < MIN_MU ) return MIN_MU;
	return xi;
}


double DistVect( double *U, double *V, int n){

	int i;
	double mydist=0;

	for(i=0;i<n;i++)
		mydist += pow(U[i]-V[i],2.0);

	return sqrt(mydist);

}
double NormVect( double *U, int n){

	int i;
	double mydist=0;

	for(i=0;i<n;i++)
		mydist += pow(U[i],2.0);

	return sqrt(mydist);

}





/*
	Although I haven't found any bug, this seems to not work efficiently for many cases.
*/
double Gradient_Climber_multi( double *g, int D, double *X, double *DX, char *mask, struct CoevolData *MyEpocsData, int IS[], double (*scenar_multi)(double *, struct CoevolData *, int[], int), int ntree ){

	extern int verbose;

	double Xb[4][8],
	       f[4]={0,},
	       dist=0,
	       *Xp[4];  /* pointer to Xb */

	int d,q,skip;
	int a=0;

	double *tmp;
	double c=0.1;

	double t=(sqrt(5)-1)/2.0;   // 1/\phi

	/*
		Copy X in all Xb
	*/
	for(d=0; d<D; d++)
		for(q=0;q<4;q++)
			Xb[q][d]=X[d];

	/*
		first ascent to the peak using the gradient
	*/
	f[0] = (*scenar_multi)( X, MyEpocsData, IS, ntree );

	/*
		use the 3 first bins (q=0,1,2) to hunt for the peak
	*/
	q = 2;
	a=0;
	do{

		q=(q+1)%3;                                     // q=0 then 1, 2, 0, etc.

		for(skip=0, d=0; d<D; d++)
		{
			if( mask[d] )	skip++;
			else			Xb[(q+1)%3][d] = check_limit_rate( Xb[q][d] + c*g[d-skip]);
		}

		// f[(q+1)%3] = (*scenar)( Xb[(q+1)%3], root, ts, IS );    // f[1]=   at first round
		f[(q+1)%3] = (*scenar_multi)( Xb[(q+1)%3], MyEpocsData, IS, ntree );

		a++;

//		printf("q = %d\n", q);
//		PrintVect( "Xb_low", Xb[q], D );
//		PrintVect( "Xb_high", Xb[(q+1)%3], D );

	}while( f[(q+1)%3] > f[q] );                                // Is f[1] > f[0] ? at first round


	/*
		At this stage
		if a == 1 :
			peak is between 0 and 1
			thus Xb[0] is fine
			     Xb[3] <- Xb[1]
		if a > 1:
			peak is in between (q+2)%3 and (q+1)%3
			thus Xb[3] <- Xb[ (q+1)%3 ]
			     Xb[0] <- Xb[ (q+2)%3 ]

	*/

	for(d=0; d<D; d++)
	{
		Xb[3][d] = Xb[(q+1)%3][d];
		if(a>1)
			Xb[0][d] = Xb[(q+2)%3][d];
	}



/*
	PrintVect( "X0", Xb[0], D);
	PrintVect( "X3", Xb[3], D);
	PrintVect( "f", f, 4);

	printf("There are %d steps in ascent\n", a);
*/



//	exit(1);

	/*
		At this point, the peak is in [ Xb[q] , Xb[(q+3)%6] ]
		thus sort by ascending order  Xb[0], Xb[1], Xb[2], Xb[3]
		Compute Xb[1] and Xb[2]
	*/
	for(skip=0, d=0; d<D; d++)
	{
		if( mask[d] )	skip++;
		else{
						Xb[1][d] = Xb[0][d] + (Xb[3][d] - Xb[0][d])*(1.0-t);
						Xb[2][d] = Xb[0][d] + (Xb[3][d] - Xb[0][d])*t;
			}
	}

	for(q=0;q<4;q++)
		Xp[q]=Xb[q];

	for(q=0;q<4;q++)
		// f[q] = (*scenar)( Xp[q], root, ts, IS );
		f[q] = (*scenar_multi)( Xp[q], MyEpocsData, IS, ntree );

/*	printf("-------------------------------\n");
	PrintVect( "X0", Xp[0], D);
	PrintVect( "X1", Xp[1], D);
	PrintVect( "X2", Xp[2], D);
	PrintVect( "X3", Xp[3], D);
	PrintVect( "f", f, 4);
	printf("-------------------------------\n");
*/
//	exit(1);

	/*
		Start splitting using 1 / golden number
	*/
	while( DistVect( Xp[0], Xp[3] ,D )/ NormVect(Xp[3], D)  > 1e-6 ){

		if( f[1]>f[2] ){

			tmp=Xp[3];

			Xp[3] = Xp[2];
			f[3] = f[2];

			Xp[2] = Xp[1];
			f[2] = f[1];

			Xp[1]=tmp;
			for(skip=0, d=0; d<D; d++)
			{
				if( mask[d] )	skip++;
				else			Xp[1][d] = Xp[0][d] + (Xp[3][d] - Xp[0][d])*(1.0-t);
			}

			// f[1] = (*scenar)( Xp[1], root, ts, IS );
			f[1] = (*scenar_multi)( Xp[1], MyEpocsData, IS, ntree );
		}else{

			tmp=Xp[0];

			Xp[0] = Xp[1];
			f[0] = f[1];

			Xp[1] = Xp[2];
			f[1] = f[2];

			Xp[2]=tmp;
			for(skip=0, d=0; d<D; d++)
			{
				if( mask[d] )	skip++;
				else			Xp[2][d] = Xp[0][d] + (Xp[3][d] - Xp[0][d])*t;
			}

			// f[2] = (*scenar)( Xp[2], root, ts, IS );
			f[2] = (*scenar_multi)( Xp[2], MyEpocsData, IS, ntree );
		}

/*		PrintVect( "> X0", Xp[0], D);
		PrintVect( "> X1", Xp[1], D);
		PrintVect( "> X2", Xp[2], D);
		PrintVect( "> X3", Xp[3], D);
		PrintVect( "> f", f, 4);
		printf("-------------------------------\n");
*/

	}

/*		printf(">>> FINAL\n");
		PrintVect( ">>> X0", Xp[0], D);
		PrintVect( ">>> X1", Xp[1], D);
		PrintVect( ">>> X2", Xp[2], D);
		PrintVect( ">>> X3", Xp[3], D);
		PrintVect( ">>> f", f, 4);
		printf("-------------------------------\n");
*/

	/*
		Compute dist
	*/
	for(dist=0, d=0; d<D; d++)
		dist += pow( X[d] - (Xp[0][d]+Xp[3][d])/2.0 , 2 );

	/*
		and update X
	*/
	for(d=0; d<D; d++)
		X[d] = (Xp[0][d]+Xp[3][d])/2.0;


	return sqrt(dist);
}

/*
	Update each dimension independently, ascending like a crab, using golden number optimisation
*/
double Stairs_Climber_multi( int D, double *X, double *DX, char *mask, struct CoevolData *MyEpocsData, int IS[], double (*scenar_multi)(double *, struct CoevolData *, int[], int), int ntree ){

	extern int verbose;

	double Xb[8],
		   xi[4],
	       f[4],
	       dist=0,
	       tmp=0,
	       dX;

	int d,q;
//	int a=0;

	double t=(sqrt(5)-1)/2.0;   // 1/\phi

	/*
		Copy X in Xb
	*/
	for(d=0; d<D; d++)
		Xb[d]=X[d];


	for(d=0;d<D;d++){

		if(mask[d]){continue;}

//		printf("Ascent dimension %d\n", d);

		/*
			first ascent to the peak in that dimension
		*/

		dX =  grad_in_xi_multi( Xb, d, MyEpocsData, IS, scenar_multi, ntree ) ;      // gradient in dimension d


		// f[0] = (*scenar)( Xb, root, ts, IS );

		f[0] = (*scenar_multi)( Xb, MyEpocsData, IS, ntree );

		xi[0] = Xb[d];

		/*
			Use xi[0], xi[1] and xi[2] to hunt for the peak
		*/
		q = 2;
//		a=0;
		do{
			q=(q+1)%3;                                     // q=0 then 1, 2, 0 etc.

			xi[(q+1)%3] = check_limit_rate( xi[q]+dX );    // xi[1]=xi[0]+dX at first round
			Xb[d]=xi[(q+1)%3];                             // replace Xb[ d ] before computing f[1]

			f[(q+1)%3] = (*scenar_multi)( Xb, MyEpocsData, IS, ntree );

//			a++;
		}while( (f[(q+1)%3] - f[q]) > 1e-10 );                       // Is f[1] > f[0] ? at first round
		// }while(f[(q+1)%3] > f[q]);


		/*
			At this point, the peak is in
				if a==1
					[ xi[0] , xi[1] ]
				else
					if( xi is at max/min and gradiant still in the direction)
						xi max/min is the peak
					else
						peak is in [ xi[(q+2)%3], xi[(q+1)%3 ]
		*/


		if( xi[(q+1)%3] == MIN_MU ){

			Xb[d] = (1+EPSILON) * MIN_MU;
			f[3] = (*scenar_multi)( Xb, MyEpocsData, IS, ntree );

			if( f[(q+1)%3 ] > f[3] )
				xi[0]=xi[3]=MIN_MU;

		}
		else if( xi[(q+1)%3] == MAX_MU ){

			Xb[d] = (1.0-EPSILON) * MAX_MU;
			f[3] = (*scenar_multi)( Xb, MyEpocsData, IS, ntree );

			if( f[(q+1)%3 ] > f[3] )
				xi[0]=xi[3]=MAX_MU;

		}
		else{

			xi[3] = xi[(q+1)%3];
//			if(a>1)
//				xi[0] = xi[(q+2)%3];                           // xi[0] is set to the ante-penultiem

		}



//		if(xi[0]>xi[3])          /* re-order if needed */
//		{
//			tmp = xi[0] ;
//			xi[0] = xi[3];
//			xi[3] = tmp;
//		}

		xi[1] = xi[0] + (xi[3] - xi[0])*(1.0-t);
		xi[2] = xi[0] + (xi[3] - xi[0]) * t;

		for(q=0;q<4;q++){
			Xb[d] = xi[q];
			f[q] = (*scenar_multi)( Xb, MyEpocsData, IS, ntree );
		}

//		PrintVect("xi", xi, 4);
//		PrintVect("f", f, 4);

		/*
			Start splitting using 1 / golden number
		*/
		while( fabs(xi[3] - xi[0])/(xi[3] + xi[0]) > 1e-7 ){

			if( f[1]>f[2] ){

//				printf("f[1] > f[2]\n");

				xi[3] = xi[2];
				f[3] = f[2];

				xi[2] = xi[1];
				f[2] = f[1];

				xi[1] = xi[0] + (xi[3] - xi[0])*(1.0-t);
				Xb[d] = xi[1];
				f[1] = (*scenar_multi)( Xb, MyEpocsData, IS, ntree );

			}else{

//				printf("f[2] >= f[1]\n");

				xi[0] = xi[1];
				f[0] = f[1];

				xi[1] = xi[2];
				f[1] = f[2];

				xi[2] = xi[0] + (xi[3] - xi[0])*t;
				Xb[d] = xi[2];
				f[2] = (*scenar_multi)( Xb, MyEpocsData, IS, ntree );

			}

//		printf("\n");
//		PrintVect("> xi", xi, 4);
//		PrintVect("> f", f, 4);

		}
		Xb[d] = (xi[0]+xi[3])/2.0;

	}

	/*
		Compute dist
	*/
	for(dist=0, d=0; d<D; d++)
		dist += pow( X[d]-Xb[d] , 2);

	/*
		and update X
	*/
	for(d=0; d<D; d++)
		X[d] = Xb[d];

//	printf("done\n");

	return sqrt(dist);
}

double ML_by_Grad_n_NR_multi( struct CoevolData *MyEpocsData, int ntree, double rate[2][4], int *IS, char scenario_multi ){

	double *g, *dX;        // the gradiant, dX: movement
	double **H, **Hinv;   // the Hessien, and H^-1

	int i;

	int d=-1, d2;              // dimensionnality
	double X[8],               // the vector to be optimized
	       Y[8];               // a backup used only for NR

	char *mask;

	double (*scenar_multi)(double *, struct CoevolData *, int[], int);

	double limit=1;
	double mydist=1;
	double f0=-1, f1=-1;

	extern int verbose;

//	verbose=3;

	/*
		Select the scenario
	*/
	scenar_multi = pick_scenario_multi ( scenario_multi, &d);

	/*
		Init X vector
	*/
	for(i=0;i<d;i++)
		X[i]=1;


	/*
		Memory Allocation
	*/
	MemSet( &g, &dX, &H, &Hinv, &mask, d );


	/*
		First pass and reduce dimension
	*/
	compute_g_multi(g, d, X, NULL, 0, MyEpocsData, IS, scenar_multi, ntree); // compute gradient for all dimensions


	d2=build_mask_multi(g, d, X, mask, MyEpocsData, IS, scenar_multi, ntree); // mask all dimensions with best ML value at 0 or with null gradient

	if(verbose>1){
		PrintVect("X0",X,d);
		PrintVect("g0",g,d);
		PrintStr("m0",mask,d);
	}

	int rounds=0;
	int NR=0;
	int GC=0;

	while( d2>0 && mydist > 1e-6 )   // while there are some dimensions to improve.
	{
		f0 = (*scenar_multi)( X, MyEpocsData, IS, ntree);

		if(verbose>1)
		{
			printf(">>>>>>> Round %d <<<<<<<<<\n",rounds);
			printf("        L= %e\n",f0);
//			if(rounds >15)exit(1);

		}

		d2=compute_g_multi( g, d, X, mask, 1, MyEpocsData, IS, scenar_multi, ntree );
		if(!d2)break;

		if(verbose>1){
			PrintVect("<<X>>",X, d);
			PrintVectMask("<<g>>",g, d, mask);
		}
		/*
			Oscillate between NR (default) or GC/SC (when NR fails)
		*/

		if( NR==0 || mydist>limit )
		{
			/*
				Oscillate between SC (first round) or GC (usually less good, but who knows)
			*/
			if(GC==1)
			{
				mydist = Gradient_Climber_multi( g, d, X, dX, mask, MyEpocsData, IS, scenar_multi, ntree );

				if(verbose>1){

					compute_g_multi( g, d, X, mask, 0, MyEpocsData, IS, scenar_multi, ntree );
					printf("<--------------------------\n");
					PrintVect("GC X",X,d);
					PrintVectMask("<g>",g,d, mask);
					PrintStr("m",mask,d);
					printf("** dist: %e dims: %d\n", mydist, d2);
					printf("-------------------------->\n");
				}

				GC=0;    // next round do SC

			}
			else
			{

				mydist = Stairs_Climber_multi( d, X, dX, mask, MyEpocsData, IS, scenar_multi, ntree );

				if(verbose>1){
					compute_g_multi( g, d, X, mask, 0, MyEpocsData, IS, scenar_multi, ntree );
					printf("<--------------------------\n");
					PrintVect("SC X",X,d);
					PrintVectMask("<g>",g,d, mask);
					PrintStr("m",mask,d);
					printf("** dist: %e dims: %d\n", mydist, d2);
					printf("-------------------------->\n");

				}

				if(rounds<200)GC=1;  // next round, do GC only in the early ones

			}

			NR=1;     // next round, do NR (if limit is ok)

//			PrintVect("X",X,d);
		}
		else
		{

			if(verbose>1){
				PrintVectMask("<g>",g,d, mask);
				PrintStr("m",mask,d);
			}
			compute_H_multi( H, d, X, mask, MyEpocsData, IS, scenar_multi, ntree );

			if(verbose)printf("try NR\n");

			invert_H( H, Hinv, d2);
			prod_M_V( Hinv, g, dX, d2);		// H, g and dX have dimension d2

			/*
				Try it on Y (a temporary copy of X)
			*/
			for(i=0;i<8;i++)Y[i]=X[i];
			update_X( Y, dX, d, g, mask );
			f1 = (*scenar_multi)( Y, MyEpocsData, IS, ntree);

			/*
				If it improves L, keep it
			*/
			if( f1>f0  )
//			if( f1>f0 && check_dX( dX, d, g, mask ) )   // same direction than g AND improve L
			{

				for(i=0;i<8;i++)X[i]=Y[i];
				mydist = norm_vect( dX, d2, NULL );

				if(verbose>1){
					PrintVect("NR X",X,d); //printf("Newton-Raphson\n");
					PrintVectMask("<g>",g,d, mask);
					printf("** dist: %e dims: %d\n", mydist, d2);
				}
				NR=1;
			}
			else
			{
				// if N-R goes in wrong direction, go back to Gradient_Climber
				limit/=10;
				if(verbose>1) printf("cannot use NR. go back to Gradient-Climber\n");
				NR=0;
			}


		}


		rounds++;

	}

	if(verbose)
		printf("#rounds: %d mydist: %e d2: %d\n", rounds, mydist, d2);


	FreeSet( g, dX, H, Hinv, mask, d );

	set_rates_scenario( scenario_multi, rate, X );

	return (*scenar_multi)( X, MyEpocsData, IS, ntree);
}

void ML_multi(struct CoevolData *MyEpocsData, int ntree, double *BestMaxlnL, double rate[2][4], double final_rate[2][4], int *IS, char scenario_multi, int *maxVector, int level, int *resVector, int *o, int *imax, int *jmax, int opt_IS, int round){
	int i, j, k;
	double maxlnL;

	for(i=0;i<2;i++){
		for(j=0;j<2;j++){

			if(!opt_IS && i+j>0)
				break;

			if( (i%2 + j%2)==2 )
				continue;

			if (verbose)
				printf("<< IS={%d,%d} >>\n",i,j);


			IS[0] = i; IS[1] = j;

			/*
			If I have one tree -> I'm moving only in the tvector corresponding to tree 0 (ntree-1)
			If I have >1 tree -> I'm moving in the tvectors with the recursion, and loop over the last tree
			*/
			MyEpocsData[ntree-1].theTVector.current_vector = MyEpocsData[ntree-1].theTVector.tvector;
			while(*MyEpocsData[ntree-1].theTVector.current_vector){

				// print_vector_types(MyEpocsData[ntree-1].theTVector.current_vector, MyEpocsData[ntree-1].nbranches, stdout);

				maxlnL = ML_by_Grad_n_NR_multi(MyEpocsData, ntree, rate, IS, scenario_multi);
				if (verbose)
					fprintf(stdout," IS[0]=%d, IS[1]=%d Maximum lnL = %f L=%g\n",IS[0],IS[1],maxlnL,exp(maxlnL));

				if (*BestMaxlnL==0 || maxlnL > *BestMaxlnL){
					*BestMaxlnL=maxlnL;
					*imax=i;
					*jmax=j;
					memcpy(o, resVector, ntree*sizeof(int));
					memcpy(final_rate, rate, 8*sizeof(double));
				}

				MyEpocsData[ntree-1].theTVector.current_vector++;
				resVector[ntree-1]+=1;
			}
			round++;

			if (ntree>1){
				for (k=level; k<ntree-1; k++){
					if ((resVector[k] - maxVector[k] + 1) < 0){

						level=k;

						resVector[k]+=1;

						MyEpocsData[k].theTVector.current_vector++;

						ML_multi(MyEpocsData, ntree, BestMaxlnL, rate, final_rate, IS, scenario_multi, maxVector, level, resVector, o, imax, jmax, opt_IS, round);

						resVector[k]-=1;
						level-=1;
						MyEpocsData[k].theTVector.current_vector--;
						print_vector_types(MyEpocsData[k].theTVector.current_vector, MyEpocsData[k].nbranches, stdout);
					}
				}
			}

		}
	}
	return;
}
