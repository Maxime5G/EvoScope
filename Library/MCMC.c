#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "tree.h"
#include "coevol.h"
#define EPSILON 0.001
#include <time.h>

#ifndef NOMINMAX

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

#endif  /* NOMINMAX */

/* Priors */

double exponential_prior(double n, double l){
    /* generates the corresponding value from an exponential distribution with value n and lambda l */
    return (l * exp(-l*n));
}

//  f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))
double normal_prior(double sigma, double mu, double n){
    /* generates the corresponding value from a normal distribution with mean mu and sd sigma */
    return ( (1/(sqrt(2*M_PI) * sigma)) * (exp(-(pow(n-mu, 2) / (2* pow(sigma,2))))) );
}

/* Sampling random numbers between min and max */
double sample_number(double min, double max){
    double rand_number = rand() / (double)RAND_MAX;
    return min + rand_number * (max-min);
}

int acceptance(double current, double proposal){
    double flip_coin =(double)rand() / (double)RAND_MAX;
    double ratio = exp(proposal-current)>1?1:exp(proposal-current);

    if (ratio==1 || ratio>flip_coin)
        return 1;
    return 0;
}

void print_rates_vector(double rate_vector[2][4]){
    for (int i=0; i<2; i++){
        for (int j=0; j<4; j++){
            fprintf(stderr, "%.5f ", rate_vector[i][j]);
        }
        fprintf(stderr, "\n");
    }
    return;
}

void normalize_a_vector_bis(double (*rate_vector)[2][4]){
    double sum=0.0L;
    int i,j;
    for (i=0; i<2; i++){
        for (j=0; j<4; j++){
            sum += (*rate_vector)[i][j];
        }
    }
    for (i=0; i<2; i++){
        for (j=0; j<4; j++){
            (*rate_vector)[i][j] /= sum;
        }
    }
}

int MCMC2( struct CoevolData *MyEpocsData, int *IS, int Nrounds, int w, int *total_acc, int sampling, FILE *fptr ){

    int i,j,k;
    double proposal[2][4];  // parameter values of the proposal
    double current[2][4];   // parameter values of the current position in the MCMC chain
    double mylnL;           // lnML value
    double posterior;       // posterior value (log(prior) + log(ML))
    double current_value;   // current value (log(prior) + log(ML))

    int acc;                // binary value (1 = acceptance of the new parameter)

    int *indices=NULL;      // indices of the co-occurrences in the transition vector
    int which_cooc=0;       // which co-occurrence to update if we update the transition vector

    char tvectortoshow[MyEpocsData[0].theTVector.nbfork+1];             // INIT: array to store the transitions in the co-occurrences (3 or 4)- used in the output
    char tvector[MyEpocsData[0].nbranches+1];                           // INIT: array to store the transition vector naming it "tvector" for clarity

    memcpy(tvector, MyEpocsData[0].theTVector.tvector[0], MyEpocsData[0].nbranches);    // copying the content of the tvector to the newly created tvector

    indices = calloc(MyEpocsData[0].theTVector.nbfork, sizeof(int));    // storing the indices of the transitions in the tvector
    for (i=0; i<MyEpocsData[0].nbranches; i++){                         // useful to show all transitions in summary form
        if ((int)tvector[i] > 2){                                       // whenever there is a co-occurrence
            tvectortoshow[which_cooc]='3';                              // setting it (for starters) to 3
            tvector[i]=(char)3;                                         // (i.e., all co-occurrences are in the direction e1->e2)
            indices[which_cooc++]=i;
        }
    }
    tvectortoshow[MyEpocsData[0].theTVector.nbfork]='\0';

    /* Initializing the starting parameter values */
    /* Default is 5 (for now...) */
    for (i=0; i<2; i++){
        for (j=0; j<4; j++){
            current[i][j] = 5;
            proposal[i][j] = 5;
        }
    }

    /* First iteration of the MCMC chain */
    current_value = log(1) + ln_vraisemblance(MyEpocsData[0].root, tvector, current, IS); // If uniform distribution, prior is 1

    which_cooc = 0;
    fprintf(fptr, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "state","log-posterior","mu1","mu1star","mu2","mu2star","nu1","nu1star","nu2","nu2star","tvector");

    fprintf(stderr, "[");

    for (k=0; k<Nrounds; k++){

        if (k%(Nrounds/10) == 0){
            // fprintf(stderr,"e%d vs e%d\tmut. order : %ld / %ld  \r",evti,evtj, (ts-tvectors)+1, (long)pow(2,nbfork));
            // fprintf(stderr, "--> round %d / %d...\r", k, Nrounds);
            // fflush(stderr);
            fprintf(stderr, "%d%%...", (k/(Nrounds/10))*10);
        }

        if (k%2 == 0){
            /*
            Updating the transition vector every other round
            Basically, we flip the co-occurrence at the position which_cooc into the other direction
            So if we have '3' (i.e., the current direction is e1 -> e2)
            We set it to 4 (i.e., the proposed direction is e2 -> e1 )
            And vice-versa
            */

            if ((int)tvector[indices[which_cooc]]>3){
                tvector[indices[which_cooc]]=3;
                tvectortoshow[which_cooc]='3';
            }
            else{
                tvector[indices[which_cooc]]=4;
                tvectortoshow[which_cooc]='4';
            }

            mylnL = ln_vraisemblance(MyEpocsData[0].root, tvector, proposal, IS);

            posterior = log(1)+mylnL;   // If uniform distribution (which is the case now), the prior is 1

            acc = acceptance(current_value, posterior); // We check whether we accept the proposed value

            if (acc){
                // If accepted, we update the current value and add one to the total acceptance
                current_value = mylnL;
                (*total_acc)++;
            }else{
                // If not, we reset the values to their previous set
                if ((int)tvector[indices[which_cooc]]>3){
                    tvector[indices[which_cooc]]=3;
                    tvectortoshow[which_cooc]='3';
                }
                else{
                    tvector[indices[which_cooc]]=4;
                    tvectortoshow[which_cooc]='4';
                }
            }
            if (which_cooc < MyEpocsData[0].theTVector.nbfork-1)
                which_cooc++;
            else
                which_cooc=0;

        }
        else{
            // update the parameters
            for (i=0; i<2; i++){
                for (j=0; j<4; j++){

                    proposal[i][j] = sample_number(current[i][j]-w, current[i][j]+w);
                    if (proposal[i][j]<0 || proposal[i][j]>1000){
                        proposal[i][j] = current[i][j];
                        continue;
                    }

                    mylnL = ln_vraisemblance(MyEpocsData[0].root, tvector, proposal, IS);

                    posterior = log(1)+mylnL;   // If uniform distribution

                    acc = acceptance(current_value, posterior);

                    if (acc){
                        current[i][j] = proposal[i][j];
                        current_value = mylnL;
                        (*total_acc)++;
                    }else{
                        proposal[i][j] = current[i][j];
                    }
                }
            }
        }

        if ((k%sampling)==0){
            fprintf(fptr, "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n", k, posterior, proposal[0][0], proposal[0][1], proposal[1][0], proposal[1][1], proposal[0][2], proposal[0][3], proposal[1][2], proposal[1][3], tvectortoshow);
        }
    }
    fprintf(stderr, "100%%]\n\n");
    return 1;
}
