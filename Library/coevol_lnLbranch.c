#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "coevol.h"

#define EPSILON (1e-4)


void final_state0(int *IS, int *FS, int nstates);
void final_state1(int *IS, int *FS, int nstates);
void final_state2(int *IS, int *FS, int nstates);
void final_state3(int *IS, int *FS, int nstates);
void final_state4(int *IS, int *FS, int nstates);
void (*pfinal_state[5])(int *IS, int *FS, int nstates) = {final_state0,final_state1,final_state2,final_state3,final_state4};

double lnLbranch0( double l, int ev, int *IS, double tau[2][4], int *FS, int nstates );
double lnLbranch12( double l, int ev, int *IS, double tau[2][4], int *FS, int nstates );
double lnLbranch34( double l, int ev, int *IS, double tau[2][4], int *FS, int nstates );
double (*plnLbranch[5])(double l, int ev, int *IS, double tau[2][4], int *FS, int nstates) = {lnLbranch0,lnLbranch12,lnLbranch12,lnLbranch34,lnLbranch34};
/*
	We checked all cases with Pat
*/

/* ----------------------------------------------------------------- */
void final_state0(int *IS, int *FS, int nstates) {
	int i;
	for (i = 0; i < nstates; i++)
		FS[i] = IS[i];
	return;
}
/* ----------------------------------------------------------------- */
void final_state1(int *IS, int *FS, int nstates) {
	if ( IS[0]%2 ==0 ) {                  // if IS1 is 0 or 2, FS2 is enhanced
		if ( IS[1] <= 1 )
			FS[1] = 1;           // if IS2 is 0 or 1, FS2 is 1
		else
			FS[1] = 3;           // if IS2 is 2 or 3, FS2 is 3
	}
	else
		FS[1] = IS[1];
	// do it after IS[0] if FS==IS, same array ;-)
	FS[0] = (IS[0] <= 1 )?  2  :  0 ;     // if IS1 is 0 or 1, FS1 is 2; otherwise it is 0
	return;
}
/* ----------------------------------------------------------------- */
void final_state2(int *IS, int *FS, int nstates) {
	if ( IS[1] % 2 ==0 ) {
		if ( IS[0] <= 1 )
			FS[0] = 1;
		else
			FS[0] = 3;
	}
	else
		FS[0] = IS[0];
	FS[1] = ( IS[1] <= 1 )?  2  :  0 ;
	return;
}
/* ----------------------------------------------------------------- */
void final_state3(int *IS, int *FS, int nstates) {
	final_state1( IS,  FS, nstates );
	final_state2( FS, FS, nstates );
	return;
}
/* ----------------------------------------------------------------- */
void final_state4(int *IS, int *FS, int nstates) {
	final_state2( IS,  FS, nstates );
	final_state1( FS, FS, nstates );
	return;
}
/* ----------------------------------------------------------------- */
double lnLbranch0( double l, int ev, int *IS, double tau[2][4], int *FS, int nstates ){
	double lnL=0;

	final_state0( IS, FS, nstates );
	lnL = -l*( tau[0][IS[0]] + tau[1][IS[1]] );
	return lnL;
}

/* ----------------------------------------------------------------- */
double lnLbranch12( double l, int ev, int *IS, double tau[2][4], int *FS, int nstates ){
	double lnL=0;
	double A,B;

	pfinal_state[ev]( IS, FS, nstates );

	A =  tau[0][IS[0]]+tau[1][IS[1]];       // sum of first rates
	B =  tau[0][FS[0]]+tau[1][FS[1]];     // sum of second rates (here FS rates)

	//	printf("A: %f , B: %f\n", A,B);

	/* old jompo: replaced by ev-1 because ev-1==0 si ev==1 and ev-1==1 if ev==2 */
	/*
	   if(ev==1)
	   lnL = log( tau[0][IS[0]] );
	   else
	   lnL = log( tau[1][IS[1]] );
	 */
	lnL = log( tau[ev-1][IS[ev-1]] );

	if ( fabs(A-B) < EPSILON )
		lnL += log( l ) - l*(B+A)/2.0;
	else
		lnL += log( (exp(-l*A) - exp(-l*B)) /(B - A) );

	return lnL;
}

/* ----------------------------------------------------------------- */
double lnLbranch34( double l, int ev, int *IS, double tau[2][4], int *FS, int nstates ){
	double lnL=0;
	static int MS[2];
	double A,B,C;
	double ab,ac,bc,abc;

	if (ev == 3) {
		final_state1( IS,  MS, nstates );
		final_state2( MS, FS, nstates );
	}
	else {
		final_state2( IS,  MS, nstates );
		final_state1( MS, FS, nstates );
	}
	A =  tau[0][IS[0]]+tau[1][IS[1]];
	B =  tau[0][MS[0]]+tau[1][MS[1]];
	C =  tau[0][FS[0]]+tau[1][FS[1]];

	ab=(A+B)/2.0;
	ac=(A+C)/2.0;
	bc=(C+B)/2.0;
	abc=(A+B+C)/3.0;

	/* old jompo replaced because see below
	   if(ev==3)
	   lnL = log(tau[0][IS[0]] * tau[1][MS[1]]);
	   else
	   lnL = log(tau[1][IS[1]] * tau[0][MS[0]]);
	 */
	/* old jompo: replaced by ev-1 because ev-1==0 si ev==3 and ev-1==1 if ev==4 */
	/*            replaced by 4-ev because 4-ev==0 si ev==4 and 4-ev==1 if ev==3 */
	lnL = log(tau[ev-3][IS[ev-3]] * tau[4-ev][MS[4-ev]]);


	/* general case, when A != B != C */
	if ( fabs(A-B) > EPSILON && fabs(C-B) > EPSILON && fabs(C-A) > EPSILON) {
		lnL += log( ( (exp(-l*A) - exp(-l*B))/(B - A) - (exp(-l*A) - exp(-l*C))/(C - A) ) / (C - B) );
	}
	else {
		/* at least one pair is smaller than epsilon */
		if ( fabs(A-B) <= EPSILON && fabs(C-B) <= EPSILON && fabs(C-A) <= EPSILON ) {
			/* limit when A = B = C */
			lnL += 2 * log(l) - l*abc - log(2);
		}
		else {
			if ( fabs(A-B) <= EPSILON  ) {
				/* limit when A = B != C */
				lnL += -ab*l + log( (l*(C-ab) - 1 + exp( (ab-C)*l ))  / pow(C-ab,2) );
			}
			else {
				if ( fabs(A-C) <= EPSILON ) {
					/* limit when A = C != B */
					lnL += -ac*l + log( (exp( (ac-B)*l) -1 -l*(ac-B)) / pow(ac-B,2) );
				}
				else {
					/* limit when B = C != A */
					lnL += -bc*l + log( (exp( (bc-A)*l) -1 -l*(bc-A)) / pow(bc-A,2) );
				}
			}
		}
	}

	return lnL;
}


/* ----------------------------------------------------------------- */
void final_state( int ev, int *IS, int *FS, int nstates){

	int i;

	switch( ev ){                             // ev: 0: no event, 1: event_1, 2: event_2
	                                          // 12: event_1 followed by event_2 ; 21: vice-versa
		case 0:
			for (i = 0; i < nstates; i++)
			   FS[i] = IS[i];
			break;

		case 1:
			if( IS[0]%2 ==0 ){          // changed place of FS[0]= because if IS=FS, we test FS[0]
										// if IS1 is 0 or 2, FS2 is enhanced
				if( IS[1]<=1 )  FS[1]=1;           // if IS2 is 0 or 1, FS2 is 1
				else        FS[1]=3;           // if IS2 is 2 or 3, FS2 is 3
			}
			else
				FS[1]=IS[1];
							// do it after IS[0] if FS==IS, same array ;-)
			FS[0] = (IS[0]<=1 )?  2  :  0 ;     // if IS1 is 0 or 1, FS1 is 2; otherwise it is 0

			break;

		case 2:
			if( IS[1]%2 ==0 ){ // changed place of FS[1]= because if IS=FS, we test FS[1]
				if( IS[0]<=1 ) FS[0]=1;
				else        FS[0]=3;
			}
			else
				FS[0]=IS[0];

			FS[1] = ( IS[1]<=1 )?  2  :  0 ;

			break;


		case 3:	/* 12 */
			final_state( 1, IS,  FS, nstates );
			final_state( 2, FS, FS, nstates );
			break;

		case 4:	/* 21 */
			final_state( 2, IS,  FS, nstates );
			final_state( 1, FS, FS, nstates );
			break;


		default:
			fprintf(stderr, "I don't get it, sorry\n");

	}

}


/* ----------------------------------------------------------------- */
/* bug corrected feb 9 2017 when A=C */

double lnLbranch( double l, int ev, int *IS, double tau[2][4], int *FS, int nstates ){

	double lnL=0;

	static int MS[2];

	double A,B,C;

	double ab,ac,bc,abc;

	// printf("Lbranch: %f, %d, %d, %d\n", l, ev, IS[0], IS[1] );

	switch( ev ){

		case 0:
			final_state( ev, IS, FS, nstates );
			lnL = -l*( tau[0][IS[0]] + tau[1][IS[1]] );
			break;

		case 1:
		case 2:

			final_state( ev, IS, FS, nstates );

			A =  tau[0][IS[0]]+tau[1][IS[1]];       // sum of first rates
			B =  tau[0][FS[0]]+tau[1][FS[1]];     // sum of second rates (here FS rates)

		//	printf("A: %f , B: %f\n", A,B);

			/* old jompo: replaced by ev-1 because ev-1==0 si ev==1 and ev-1==1 if ev==2 */
			/*
			if(ev==1)
				lnL = log( tau[0][IS[0]] );
			else
				lnL = log( tau[1][IS[1]] );
			*/
			lnL = log( tau[ev-1][IS[ev-1]] );
			if ( fabs(A-B) < EPSILON )
				lnL += log( l ) - l*(B+A)/2.0;
			else
				lnL += log( (exp(-l*A) - exp(-l*B)) /(B - A) );


			break;

		case 3: /* 12 */
		case 4: /* 21 */
			// OLD jompo final_state( 1,  IS,  MS, nstates );
			final_state( !(4-ev)+1 ,  IS,  MS, nstates );
/*			{int i;
			for (i=0;i<2;i++) {
			fprintf(stderr,"IS[%d]=%d\n",i,IS[i]);
			}
			for (i=0;i<2;i++) {
			fprintf(stderr,"MS[%d]=%d\n",i,MS[i]);
			}
			}
*/
			// OLD jompo final_state( 2,  MS, FS, nstates );
			final_state( (4-ev)+1 ,  MS, FS, nstates );
/*
			{int i;
			for (i=0;i<2;i++) {
			fprintf(stderr,"MS[%d]=%d\n",i,MS[i]);
			}
			for (i=0;i<2;i++) {
			fprintf(stderr,"FS[%d]=%d\n",i,FS[i]);
			}
			}
*/
			A =  tau[0][IS[0]]+tau[1][IS[1]];
			B =  tau[0][MS[0]]+tau[1][MS[1]];
			C =  tau[0][FS[0]]+tau[1][FS[1]];

//			printf("A: %f B: %f C: %f\n",A,B,C);

			ab=(A+B)/2.0;
			ac=(A+C)/2.0;
			bc=(C+B)/2.0;
			abc=(A+B+C)/3.0;

			/* old jompo replaced because see below
			if(ev==3)
				lnL = log(tau[0][IS[0]] * tau[1][MS[1]]);
			else
				lnL = log(tau[1][IS[1]] * tau[0][MS[0]]);
			*/
/*			if(ev==3 )
				lnL = log(tau1[IS1] * tau2[MS2]);
			else
				lnL = log(tau2[IS2] * tau1[MS1]);
*/
			/* old jompo: replaced by ev-1 because ev-1==0 si ev==3 and ev-1==1 if ev==4 */
			/*            replaced by 4-ev because 4-ev==0 si ev==4 and 4-ev==1 if ev==3 */
			lnL = log(tau[ev-3][IS[ev-3]] * tau[4-ev][MS[4-ev]]);


			/* general case, when A != B != C */
			if ( fabs(A-B) > EPSILON && fabs(C-B) > EPSILON && fabs(C-A) > EPSILON) {
				lnL += log( ( (exp(-l*A) - exp(-l*B))/(B - A) - (exp(-l*A) - exp(-l*C))/(C - A) ) / (C - B) );
			}
			else {
				/* at least one pair is smaller than epsilon */
				if ( fabs(A-B) <= EPSILON && fabs(C-B) <= EPSILON && fabs(C-A) <= EPSILON ) {
					/* limit when A = B = C */
					lnL += 2 * log(l) - l*abc - log(2);
				}
				else {
					if ( fabs(A-B) <= EPSILON  ) {
						/* limit when A = B != C */
						lnL += -ab*l + log( (l*(C-ab) - 1 + exp( (ab-C)*l ))  / pow(C-ab,2) );

					}
					else {
						if ( fabs(A-C) <= EPSILON ) {
							/* limit when A = C != B */
							lnL += -ac*l + log( (exp( (ac-B)*l) -1 -l*(ac-B)) / pow(ac-B,2) );
						}
						else {
							/* limit when B = C != A */
							lnL += -bc*l + log( (exp( (bc-A)*l) -1 -l*(bc-A)) / pow(bc-A,2) );
						}
					}
				}
			}

	}

//	printf("Lbranch: %f, %d, %.10f\n", l, ev, exp(lnL) );

	return lnL;
}


/*
	We checked all cases with Pat
*/
void final_state_ori( int ev, int IS1, int IS2,  int *FS1, int *FS2  ){

	switch( ev ){                             // ev: 0: no event, 1: event_1, 2: event_2
	                                          // 12: event_1 followed by event_2 ; 21: vice-versa
		case 0:
			*FS1 = IS1;
			*FS2 = IS2;
			break;

		case 1:
			*FS1 = ( IS1<=1 )?  2  :  0 ;     // if IS1 is 0 or 1, FS1 is 2; otherwise it is 0

			if( IS1%2 ==0 ){                  // if IS1 is 0 or 2, FS2 is enhanced
				if( IS2<=1 )*FS2=1;           // if IS2 is 0 or 1, FS2 is 1
				else        *FS2=3;           // if IS2 is 2 or 3, FS2 is 3
			}
			else
				*FS2=IS2;

			break;

		case 2:
			*FS2 = ( IS2<=1 )?  2  :  0 ;

			if( IS2%2 ==0 ){
				if( IS1<=1 )*FS1=1;
				else        *FS1=3;
			}
			else
				*FS1=IS1;

			break;


		case 3 /* 12 */:
			final_state_ori( 1, IS1, IS2,  FS1, FS2 );
			final_state_ori( 2, *FS1, *FS2,  FS1, FS2 );
			break;

		case 4 /* 21 */:
			final_state_ori( 2, IS1, IS2,  FS1, FS2 );
			final_state_ori( 1, *FS1, *FS2,  FS1, FS2 );
			break;


		default:
			fprintf(stderr, "I don't get it, sorry\n");

	}

}


double lnLbranch_ori( double l, int ev, int IS1, int IS2, double *tau1, double *tau2, int *FS1, int *FS2  ){

	double lnL=0;

	int MS1,
	    MS2;

	double A,B,C;

	double ab,ac,bc,abc;

	switch( ev ){

		case 0:
			final_state_ori( ev,  IS1,  IS2,  FS1, FS2  );
			lnL = -l*( tau1[IS1] + tau2[IS2] );
			break;

		case 1:
		case 2:

			final_state_ori( ev,  IS1,  IS2,  FS1, FS2  );

			A =  tau1[IS1]+tau2[IS2];       // sum of first rates
			B =  tau1[*FS1]+tau2[*FS2];     // sum of second rates (here FS rates)

		//	printf("A: %f , B: %f\n", A,B);


			if(ev==1) lnL = log( tau1[IS1] );
			else      lnL = log( tau2[IS2] );

			if( fabs(A-B) < EPSILON )	lnL += log( l ) - l*(B+A)/2.0;
			else						lnL += log( (exp(-l*A) - exp(-l*B)) /(B - A) );


			break;

		case 3 /* 12 */:
		case 4 /* 21 */:

			final_state_ori( 1,  IS1,  IS2,  &MS1, &MS2  );
			final_state_ori( 2,  MS1,  MS2,  FS1, FS2  );

			A =  tau1[IS1]+tau2[IS2];
			B =  tau1[MS1]+tau2[MS2];
			C =  tau1[*FS1]+tau2[*FS2];

			ab=(A+B)/2.0;
			ac=(A+C)/2.0;
			bc=(C+B)/2.0;
			abc=(A+B+C)/3.0;

			if(ev==3 /* 12 */)
				lnL = log(tau1[IS1] * tau2[MS2]);
			else
				lnL = log(tau2[IS2] * tau1[MS1]);


			/* general case, when A != B != C */
			if( fabs(A-B) > EPSILON && fabs(C-B) > EPSILON && fabs(C-A) > EPSILON)
			{
				lnL += log( ( (exp(-l*A) - exp(-l*B))/(B - A) - (exp(-l*A) - exp(-l*C))/(C - A) ) / (C - B) );
			}
			else
			{
				/* at least one pair is smaller than epsilon */

				if( fabs(A-B) <= EPSILON && fabs(C-B) <= EPSILON && fabs(C-A) <= EPSILON )
				{
					/* limit when A = B = C */
					lnL += 2 * log(l) - l*abc - log(2);
				}
				else
				{

					if( fabs(A-B) <= EPSILON  )
					{
						/* limit when A = B != C */
						lnL += -ab*l + log( (l*(C-ab) - 1 + exp( (ab-C)*l ))  / pow(C-ab,2) );
					}
					else
					{
						if( fabs(A-C) <= EPSILON )
						{
							/* limit when A = C != B */
							lnL += -ac*l + log( (exp( (ac-B)*l) -1 -l*(ac-B)) / pow(ac-B,2) );

						}
						else
						{
							/* limit when B = C != A */
							lnL += -bc*l + log( (exp( (bc-A)*l) -1 -l*(bc-A)) / pow(bc-A,2) );
						}
					}
				}
			}

	}

	return lnL;
}

#if 0
/*
int main( int argc, char **argv ){

	int IS[2]={0,0},            /* states are 0: mu, 1: mu*, 2: nu, 3: nu* */
	    FS[2]={0,0};

	int ev=0, i;

	double mu[2]={0.2,0.1},
	       lambda[2]={1.0,1.0},
	       kappa[2]={1.0,1.0};


	double tau1[4],
	       tau2[4];

	double lnL=0.0;


	if( argc < 4 )
		fprintf(stderr, "usage is '%s ev IS1 IS2', bye\n", argv[0]),exit(1);

	ev= atoi(argv[1]);
	IS[0]= atoi(argv[2]);
	IS[1]= atoi(argv[3]);

	tau1[0] = mu[0];
	tau1[1] = mu[0]*lambda[0];
	tau1[2] = mu[0]*kappa[0];
	tau1[3] = mu[0]*kappa[0]*lambda[0];

	tau2[0] = mu[1];
	tau2[1] = mu[1]*lambda[1];
	tau2[2] = mu[1]*kappa[1];
	tau2[3] = mu[1]*kappa[1]*lambda[1];



	/*
		br 4 -  l=2, ev=0
	*/
	lnL = lnLbranch( 2, 0, IS, tau1, tau2, FS, 2 );


	/*
		br 1 -l=1, ev=1
	*/
	lnL += lnLbranch( 1, 1, IS, tau1, tau2, FS, 2 );

	for (i = 0; i < 2; i++)
	   IS[i]=FS[i];

	/*
		br 2 -  l=1, ev=2
	*/
	lnL += lnLbranch( 1, 2, IS, tau1, tau2, FS, 2 );

	/*
		br 3 -  l=1, ev=2
	*/
	lnL += lnLbranch( 1, 2, IS, tau1, tau2, FS, 2 );


	printf("L: %.10f\n",  exp(lnL));


	return 0;
}
*/
#endif



/* ----------------------------------------------------------------------------	*/
/* calcule les vraisemblances a partir du noeud n    				*/
/* vector_evt_type: vecteur des types d'evenement           			*/
/* tau: matrice des tau (mu, mu*, etc... 					*/
/* IS: intial states, FS: final states						*/
/* ----------------------------------------------------------------------------	*/
double ln_vraisemblance(Node *n, char *vector_evt_type, double tau[2][4], int *IS )
{
	// printf("Entering ln_vraisemblance\n");
	// printf("i'm at node %d\n", n->id);
	// printf("node length = %f\n", n->time);
	// print_vector_types(&vector_evt_type, 6, stderr);
	// printf("vector_evt_type: [");
	// for (int k=0; k<6; k++){
	// 	printf("%c, ", vector_evt_type[k]);
	// }
	// printf("]\n");
	int i, FS[2];
	/* int NIS[2], j */;
	double mylnL = 0.0L;

	/* at the root, nothing to do , just call child with lnL = 0 */
	if (n->id != 0) {
#if 1
		mylnL = lnLbranch(n->time, vector_evt_type[n->id-1], IS, tau, FS, 2 ); /* id-1 is the index of node id (from 1 to N) in vector (from 0 to N-1) */
#else
		mylnL = lnLbranch_ori( n->time , vector_evt_type[n->id-1], IS[0], IS[1], tau[0], tau[1], FS, FS+1  );
#endif
//		fprintf(stderr, "branch=%d l=%5.3f lnL=%f L=%f ev=%d ev propre=%d IS1=%d IS2=%d FS1=%d FS2=%d\n",n->id,n->time,mylnL,exp(mylnL),vector_evt_type[n->id-1],n->evtype,IS[0],IS[1],FS[0],FS[1]);
	}
	else {
		/* at root, take input state given */
		for (i = 0; i < 2; i++)
			FS[i]=IS[i];
	}

	/* calcule lnL des fils et les ajoute a son lnL */
	/* l'IS des fils est le FS de la branche */
	for ( i = 0; i < n->nbdesc; i++) {
		/* pour ne pas toucher FS: for (j=0;j<2;j++) NIS[j]=FS[j]; */
		mylnL += ln_vraisemblance(n->descs[i], vector_evt_type, tau, FS );
	}

	// printf("mylnL = %f\n", mylnL);
	return mylnL;

}
#if 0
		switch (n->evt[0]) {
			case 0:
				switch  (n->evt[1]) {
					case 0: ev = 0; break;
					case 1: ev = 1; break;
					case 2: ev = 2; break;
				}
				break;
			case 1:
				switch  (n->evt[1]) {
					case 0: ev = 1; break;
					case 1: ev = 1; break; /* should not occurs */
					case 2: ev = 12; break;
				}
				break;
			case 2:
				switch  (n->evt[1]) {
					case 0: ev = 2; break;
					case 1: ev = 21; break;
					case 2: ev = 2; break; /* should not occurs */
				}
				break;
		}

#endif
