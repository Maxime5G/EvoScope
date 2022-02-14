#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "tree.h"
#include "coevol.h"
#define EPSILON 0.001
#include <time.h>

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

    // printf("proposal = %f current = %f flip_coin = %f, ratio = %f\n", proposal, current, flip_coin, ratio);

    if (ratio==1 || ratio>flip_coin)
        return 1;
    return 0;
}

void print_rates_vector(double rate_vector[2][4]){
    for (int i=0; i<2; i++){
        for (int j=0; j<4; j++){
            printf("%.5f ", rate_vector[i][j]);
        }
        printf("\n");
    }
    return;
}

int MCMC( struct CoevolData *MyEpocsData, double *results, int ntree, double rate[2][4], int *IS, int ML_tvector, int Nrounds, int w, int *total_acc ){

    int i,j,k;
    double mean[2][4];
    double proposal[2][4];
    double current[2][4];
    double sd=10.0;
    double mylnL;
    double posterior;
    double prior;
    int acc;
    double current_value;

    for (i=0; i<2; i++){
        for (j=0; j<4; j++){
            mean[i][j] = rate[i][j];
            current[i][j] = 1;
            proposal[i][j] = 1;
        }
    }

    // memcpy((void *)current, (void *)rate, (size_t)8*sizeof(double));
    // memcpy((void *)proposal, (void *)rate, (size_t)8*sizeof(double));

    // printf("MCMC_Chain\tmu1\tmu1*\tmu2\tmu2*\n");
    k=0;
    while (1){
        for (i=0; i<2; i++){
            // for (j=0; j<4; j++){
            for (j=0; j<2; j++){
                // current_value = log(normal_prior(sd, mean[i][j], current[i][j])) + ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], current, IS);
                current_value = log(1) + ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], current, IS); // If uniform distribution
                // proposal[i][j] = sample_number(mean[i][j]-w, mean[i][j]+w);
                proposal[i][j] = sample_number(current[i][j]-w, current[i][j]+w);

                if (proposal[i][j]<0)
                    proposal[i][j] = -proposal[i][j];

                proposal[i][j+2] = proposal[i][j];

                // print_rates_vector(current);
                // print_rates_vector(proposal);

                mylnL = ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], proposal, IS);
                // printf("mylnL = %f\n", mylnL);
                // posterior = log(normal_prior(sd, mean[i][j], proposal[i][j]))+mylnL;
                posterior = log(1)+mylnL;   // If uniform distribution

                acc = acceptance(current_value, posterior);
                // printf("acc = %d\n", acc);
                if (acc){
                    // printf("here\n");
                    current[i][j] = proposal[i][j];
                    current[i][j+2] = proposal[i][j];
                    (*total_acc)++;
                }else{
                    proposal[i][j]=current[i][j];
                    proposal[i][j+2]=current[i][j];
                }

                results[k] = posterior;
                // printf("%f\t%f\t%f\t%f\t%f\n", posterior, proposal[0][0], proposal[0][1], proposal[1][0], proposal[1][1]);
                if (k==Nrounds)
                    return 0;
                k++;
            }
        }
    }
    return 1;
}

int MCMC2( struct CoevolData *MyEpocsData, double *results, int ntree, double rate[2][4], int *IS, int ML_tvector, int Nrounds, int w, int *total_acc ){

    int i,j,k;
    double mean[2][4];
    double proposal[2][4];
    double current[2][4];
    double sd=10.0;
    double mylnL;
    double posterior;
    double prior;
    int acc;
    double current_value;

    for (i=0; i<2; i++){
        for (j=0; j<4; j++){
            mean[i][j] = rate[i][j];
            current[i][j] = 1;
            proposal[i][j] = 1;
        }
    }

    // memcpy((void *)current, (void *)rate, (size_t)8*sizeof(double));
    // memcpy((void *)proposal, (void *)rate, (size_t)8*sizeof(double));

    current_value = log(1) + ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], current, IS); // If uniform distribution

    printf("MCMC_Chain\tmu1\tmu1*\tmu2\tmu2*\n");
    // k=0;
    for (k=0; k<Nrounds; k++){
        for (i=0; i<2; i++){
            // for (j=0; j<4; j++){
            for (j=0; j<2; j++){

                proposal[i][j] = sample_number(current[i][j]-w, current[i][j]+w);

                if (proposal[i][j]<0)
                    proposal[i][j] = -proposal[i][j];

                proposal[i][j+2] = proposal[i][j];

            }
        }
        mylnL = ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], proposal, IS);
        // printf("mylnL = %f\n", mylnL);
        // posterior = log(normal_prior(sd, mean[i][j], proposal[i][j]))+mylnL;
        posterior = log(1)+mylnL;   // If uniform distribution
        acc = acceptance(current_value, posterior);
        // printf("acc = %d\n", acc);
        if (acc){
            for (i=0; i<2; i++){
                for (j=0; j<2; j++){
                    current[i][j] = proposal[i][j];
                    current[i][j+2] = proposal[i][j];
                }
            }
            current_value = mylnL;
            (*total_acc)++;
        }else{
            for (i=0; i<2; i++){
                for (j=0; j<2; j++){
                    proposal[i][j] = current[i][j];
                    proposal[i][j+2] = current[i][j];
                }
            }
        }

        printf("%f\t%f\t%f\t%f\t%f\n", posterior, proposal[0][0], proposal[0][1], proposal[1][0], proposal[1][1]);
        // results[k] = posterior;
    }
    return 1;
}

//

//
// int summary(){
//
// }

