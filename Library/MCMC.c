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

    // printf("proposal = %f current = %f flip_coin = %f, ratio = %f\n", proposal, current, flip_coin, ratio);

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
    // fprintf(stderr, "sum=%f\n", sum);
    for (i=0; i<2; i++){
        for (j=0; j<4; j++){
            (*rate_vector)[i][j] /= sum;
        }
    }
}

int MCMC2( struct CoevolData *MyEpocsData, int ntree, double rate[2][4], int *IS, int tvector, int Nrounds, int w, int *total_acc, int sampling ){

    int i,j,k;
    double mean[2][4];
    double proposal[2][4];
    double current[2][4];
    double mylnL;
    double posterior;
    int acc;
    double current_value;

    int *indices=NULL;
    double flipcoin;
    int rand_cooc;
    int count=0;

    char tvectortoshow[MyEpocsData[0].theTVector.nbfork+1];

    indices = calloc(MyEpocsData[0].theTVector.nbfork, sizeof(int));
    // tvectortoshow = malloc(sizeof(char)*MyEpocsData[0].theTVector.nbfork+1);
    // printf("nbfork=%d\n", MyEpocsData[0].theTVector.nbfork);
    for (i=0; i<MyEpocsData[0].nbranches; i++){
        if ((int)MyEpocsData[0].theTVector.tvector[tvector][i] > 2){
            // printf("count = %d\n", count);
            // printf("value = %c\n", atoi(&MyEpocsData[0].theTVector.tvector[tvector][i]));
            // sprintf(&tvectortoshow[count], "%d", MyEpocsData[0].theTVector.tvector[tvector][i]);
            // tvectortoshow[count]=MyEpocsData[0].theTVector.tvector[tvector][i];
            tvectortoshow[count]='3';
            MyEpocsData[0].theTVector.tvector[tvector][i]=(char)3;
            // printf("tv[count] = %c\n", tvectortoshow[count]);
            indices[count++]=i;
        }
    }
    tvectortoshow[MyEpocsData[0].theTVector.nbfork]='\0';
    //
    // printf("tvector=");
    // for (i=0; i<MyEpocsData[0].nbranches; i++)
    //     printf("%d", MyEpocsData[0].theTVector.tvector[tvector][i]);
    // printf("\nindices=");
    // for (i=0; i<MyEpocsData[0].theTVector.nbfork; i++){
    //     printf("%d ", indices[i]);
    // }
    // printf("\n");
    //
    // printf("test:");
    // for (i=0; i<MyEpocsData[0].theTVector.nbfork; i++){
    //     printf("%c", tvectortoshow[i]);
    // }
    // printf("\n");

    // printf("tvector to show=%s\n", tvectortoshow);

    // exit(1);

    for (i=0; i<2; i++){
        for (j=0; j<4; j++){
            mean[i][j] = rate[i][j];
            current[i][j] = 5;
            proposal[i][j] = 5;
        }
    }

    current_value = log(1) + ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[tvector], current, IS); // If uniform distribution
    rand_cooc = 0;
    // printf("state\tlog-posterior\tmu1\tmu1*\tmu2\tmu2*\tnu1\tnu1*\tnu2\tnu2*\n");
    printf("state\tlog-posterior\tmu1\tmu1*\tmu2\tmu2*\tnu1\tnu1*\tnu2\tnu2*\ttvector\n");
    for (k=0; k<Nrounds; k++){
        flipcoin = (double)rand() / (double)RAND_MAX;
        // printf("flipcoin = %f\n", flipcoin);
        if (flipcoin <= 0.5){
        // if (k%10 == 0){
            // update the tvectors
            // rand_cooc = (int) rand()% MyEpocsData[0].theTVector.nbfork;
            // printf("rand cooc = %d\n", rand_cooc);
            // printf("index = %d\n", indices[rand_cooc]);
            // printf("current = %d\n", MyEpocsData[0].theTVector.tvector[tvector][indices[rand_cooc]]);
            // printf("current vector = %s\n", tvectortoshow);
            // printf("tvector_current=");
            // for (i=0; i<MyEpocsData[0].nbranches; i++)
            //     printf("%d", MyEpocsData[0].theTVector.tvector[tvector][i]);
            // printf("\n");
            if ((int)MyEpocsData[0].theTVector.tvector[tvector][indices[rand_cooc]]>3){
                MyEpocsData[0].theTVector.tvector[tvector][indices[rand_cooc]]=3;
                tvectortoshow[rand_cooc]='3';
            }
            else{
                MyEpocsData[0].theTVector.tvector[tvector][indices[rand_cooc]]=4;
                tvectortoshow[rand_cooc]='4';
            }
            // (int)MyEpocsData[0].theTVector.tvector[tvector][rand_cooc]>3?MyEpocsData[0].theTVector.tvector[tvector][rand_cooc]=3:MyEpocsData[0].theTVector.tvector[tvector][rand_cooc]=4;
            // printf("updated = %d\n", MyEpocsData[0].theTVector.tvector[tvector][indices[rand_cooc]]);
            // printf("updated vector = %s\n", tvectortoshow);
            // printf("tvector_updated=");
            // for (i=0; i<MyEpocsData[0].nbranches; i++)
            //     printf("%d", MyEpocsData[0].theTVector.tvector[tvector][i]);
            // printf("\n");

            mylnL = ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[tvector], proposal, IS);

            posterior = log(1)+mylnL;   // If uniform distribution

            acc = acceptance(current_value, posterior);

            if (acc){
                // current[i][j] = proposal[i][j];
                current_value = mylnL;
                (*total_acc)++;
            }else{
                if ((int)MyEpocsData[0].theTVector.tvector[tvector][indices[rand_cooc]]>3){
                    MyEpocsData[0].theTVector.tvector[tvector][indices[rand_cooc]]=3;
                    tvectortoshow[rand_cooc]='3';
                }
                else{
                    MyEpocsData[0].theTVector.tvector[tvector][indices[rand_cooc]]=4;
                    tvectortoshow[rand_cooc]='4';
                }
            }
            if (rand_cooc < MyEpocsData[0].theTVector.nbfork-1)
                rand_cooc++;
            else
                rand_cooc=0;
            // rand_cooc < MyEpocsData[0].theTVector.nbfork? rand_cooc++ : rand_cooc = 0;
            // if ((k%sampling)==0)
                // printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n", k, posterior, proposal[0][0], proposal[0][1], proposal[1][0], proposal[1][1], proposal[0][2], proposal[0][3], proposal[1][2], proposal[1][3], tvectortoshow);
            // exit(1);
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

                    mylnL = ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[tvector], proposal, IS);

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

        if ((k%sampling)==0)
            printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\n", k, posterior, proposal[0][0], proposal[0][1], proposal[1][0], proposal[1][1], proposal[0][2], proposal[0][3], proposal[1][2], proposal[1][3], tvectortoshow);
    }
    return 1;
}


// OLD

// int MCMC2_NoTVECTOR( struct CoevolData *MyEpocsData, int ntree, double rate[2][4], int *IS, int ML_tvector, int Nrounds, int w, int *total_acc, int sampling ){
//
//     int i,j,k;
//     double mean[2][4];
//     double proposal[2][4];
//     double current[2][4];
//     double mylnL;
//     double posterior;
//     int acc;
//     double current_value;
//     double prior;
//
//     for (i=0; i<2; i++){
//         for (j=0; j<4; j++){
//             mean[i][j] = rate[i][j];
//             current[i][j] = 5;
//             proposal[i][j] = 5;
//         }
//     }
//     current_value = log(1) + ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], current, IS); // If uniform distribution
//
//     printf("state\tlog-posterior\tmu1\tmu1*\tmu2\tmu2*\tnu1\tnu1*\tnu2\tnu2*\n");
//     for (k=0; k<Nrounds; k++){
//         for (i=0; i<2; i++){
//             for (j=0; j<4; j++){
//
//                 proposal[i][j] = sample_number(current[i][j]-w, current[i][j]+w);
//                 if (proposal[i][j]<0 || proposal[i][j]>1000){
//                     proposal[i][j] = current[i][j];
//                     continue;
//                 }
//
//                 mylnL = ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], proposal, IS);
//                 posterior = log(1)+mylnL;   // If uniform distribution
//                 acc = acceptance(current_value, posterior);
//                 if (acc){
//                     current[i][j] = proposal[i][j];
//                     current_value = mylnL;
//                     (*total_acc)++;
//                 }else{
//                     proposal[i][j] = current[i][j];
//                 }
//             }
//         }
//
//         if ((k%sampling)==0)
//             printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", k, posterior, proposal[0][0], proposal[0][1], proposal[1][0], proposal[1][1], proposal[0][2], proposal[0][3], proposal[1][2], proposal[1][3]);
//     }
//     return 1;
// }
//
// int MCMC( struct CoevolData *MyEpocsData, double *results, int ntree, double rate[2][4], int *IS, int ML_tvector, int Nrounds, int w, int *total_acc ){
//
//     int i,j,k;
//     double mean[2][4];
//     double proposal[2][4];
//     double current[2][4];
//     // double sd=10.0;
//     double mylnL;
//     double posterior;
//     // double prior;
//     int acc;
//     double current_value;
//
//     for (i=0; i<2; i++){
//         for (j=0; j<4; j++){
//             mean[i][j] = rate[i][j];
//             current[i][j] = 1;
//             proposal[i][j] = 1;
//         }
//     }
//
//     // memcpy((void *)current, (void *)rate, (size_t)8*sizeof(double));
//     // memcpy((void *)proposal, (void *)rate, (size_t)8*sizeof(double));
//
//     // printf("MCMC_Chain\tmu1\tmu1*\tmu2\tmu2*\n");
//     k=0;
//     while (1){
//         for (i=0; i<2; i++){
//             // for (j=0; j<4; j++){
//             for (j=0; j<2; j++){
//                 // current_value = log(normal_prior(sd, mean[i][j], current[i][j])) + ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], current, IS);
//                 current_value = log(1) + ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], current, IS); // If uniform distribution
//                 // proposal[i][j] = sample_number(mean[i][j]-w, mean[i][j]+w);
//                 proposal[i][j] = sample_number(current[i][j]-w, current[i][j]+w);
//
//                 if (proposal[i][j]<0)
//                     proposal[i][j] = -proposal[i][j];
//
//                 proposal[i][j+2] = proposal[i][j];
//
//                 // print_rates_vector(current);
//                 // print_rates_vector(proposal);
//
//                 mylnL = ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], proposal, IS);
//                 // printf("mylnL = %f\n", mylnL);
//                 // posterior = log(normal_prior(sd, mean[i][j], proposal[i][j]))+mylnL;
//                 posterior = log(1)+mylnL;   // If uniform distribution
//
//                 acc = acceptance(current_value, posterior);
//                 // printf("acc = %d\n", acc);
//                 if (acc){
//                     // printf("here\n");
//                     current[i][j] = proposal[i][j];
//                     current[i][j+2] = proposal[i][j];
//                     (*total_acc)++;
//                 }else{
//                     proposal[i][j]=current[i][j];
//                     proposal[i][j+2]=current[i][j];
//                 }
//
//                 results[k] = posterior;
//                 // printf("%f\t%f\t%f\t%f\t%f\n", posterior, proposal[0][0], proposal[0][1], proposal[1][0], proposal[1][1]);
//                 if (k==Nrounds)
//                     return 0;
//                 k++;
//             }
//         }
//     }
//     return 1;
// }
//
//
// int MCMC3( struct CoevolData *MyEpocsData, double *results, int ntree, double rate[2][4], int *IS, int ML_tvector, int Nrounds, int w, int *total_acc ){
//
//     int i,j,k;
//     double mean[2][4];
//     double proposal[2][4];
//     double current[2][4];
//     // double sd=10.0;
//     double mylnL;
//     double posterior;
//     // double prior;
//     int acc;
//     double current_value;
//     double prior;
//     double current_w_vec[2][4];
//     double posterior_w_vec[2][4];
//     double weight_sugg;
//     double flipcoin;
//     double my_sum[8];
//     double my_val;
//
//     for (i=0; i<2; i++){
//         for (j=0; j<4; j++){
//             mean[i][j] = rate[i][j];
//             current[i][j] = 1;
//             proposal[i][j] = 1;
//             // w_vec[i][j] = 5;
//             current_w_vec[i][j] = .125;
//             posterior_w_vec[i][j] = .125;
//         }
//     }
//
//     current_value = log(1) + ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], current, IS); // If uniform distribution
//     // print_rates_vector(current_w_vec);
//     printf("state\tlog-posterior\tmu1\tmu1*\tmu2\tmu2*\tnu1\tnu1*\tnu2\tnu2*\n");
//     for (k=0; k<Nrounds; k++){
//         flipcoin = (double)rand() / (double)RAND_MAX;
//         // my_val=0.0L;
//         // for (i=0; i<2; i++){
//         //     for (j=0; j<4; j++){
//         //         posterior_w_vec[i][j] = min(1,fabs(sample_number(current_w_vec[i][j]-0.005, current_w_vec[i][j]+0.005)));
//         //     }
//         // }
//         // // print_rates_vector(posterior_w_vec);
//         // normalize_a_vector_bis(&posterior_w_vec);
//         // print_rates_vector(posterior_w_vec);
//
//         for (i=0; i<2; i++){
//             for (j=0; j<4; j++){
//                 if (i==0){
//                     if (j == 0){
//                         my_sum[j] = current_w_vec[i][j];
//                     }else{
//                         my_sum[j] = current_w_vec[i][j] + my_sum[j-1];
//                     }
//                 }else{
//                     if (j == 0){
//                         my_sum[j+4] = current_w_vec[i][j] + my_sum[j+3];
//                     }else{
//                         my_sum[j+4] = current_w_vec[i][j] + my_sum[j+3];
//                     }
//                 }
//             }
//         }
//         // print_rates_vector(posterior_w_vec);
//         // fprintf(stderr, "flipcoin = %f\n", flipcoin);
//         // fprintf(stderr, "myvec = ");
//         // for (i=0; i<8; i++){
//         //     fprintf(stderr, "%f, ", my_sum[i]);
//         // }
//         // fprintf(stderr,"\n");
//
//         for (i=0; i<2; i++){
//             for (j=0; j<4; j++){
//                 if (i==0) my_val=my_sum[j];
//                 if (i==1) my_val=my_sum[j+4];
//                 // fprintf(stderr, "myval %f flipcoin %f\n", my_val, flipcoin);
//                 if (flipcoin > my_val)
//                     continue;
//                 else{
//
//                     fprintf(stderr, "i %d j %d\n", i, j);
//                     proposal[i][j] = sample_number(current[i][j]-w, current[i][j]+w);
//                     // if (proposal[i][j]<0)
//                         // proposal[i][j] = -proposal[i][j];
//
//                     if (proposal[i][j]<0 || proposal[i][j]>1000){
//                         proposal[i][j] = current[i][j];
//                         // current_w_vec[i][j] = min(1,current_w_vec[i][j]+0.0001);
//                         // current_w_vec[i][j] = max(0, current_w_vec[i][j]-0.0001);
//                         // normalize_a_vector_bis(&current_w_vec);
//                         continue;
//                     }
//
//                     mylnL = ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], proposal, IS);
//                     posterior = log(1)+mylnL;   // If uniform distribution
//                     acc = acceptance(current_value, posterior);
//                     if (acc){
//                         current[i][j] = proposal[i][j];
//                         // current_w_vec[i][j] = max(0, current_w_vec[i][j]-0.0001);
//                         // current_w_vec[i][j] = min(1,current_w_vec[i][j]+0.0001);
//                         // normalize_a_vector_bis(&current_w_vec);
//                         // for (i=0; i<2; i++){
//                         //     for (j=0; j<4; j++){
//                         //         current_w_vec[i][j] = posterior_w_vec[i][j];
//                         //     }
//                         // }
//                         current_value = mylnL;
//                         (*total_acc)++;
//                     }else{
//                         proposal[i][j] = current[i][j];
//                         // current_w_vec[i][j] = min(1,current_w_vec[i][j]+0.0001);
//                         // current_w_vec[i][j] = max(0, current_w_vec[i][j]-0.0001);
//                         // normalize_a_vector_bis(&current_w_vec);
//                     }
//                     break;
//                 }
//                 break;
//
//                 // weight_sugg = fabs(sample_number(w_vec[i][j]-0.025, w_vec[i][j]+0.025));
//                 // // proposal[i][j] = sample_number(current[i][j]-w_vec[i][j], current[i][j]+w_vec[i][j]);
//                 // proposal[i][j] = sample_number(current[i][j]-w, current[i][j]+w);
//                 //
//                 // if (proposal[i][j]<0)
//                 //     proposal[i][j] = -proposal[i][j];
//                 //
//                 // if (proposal[i][j]>1000){
//                 //     proposal[i][j] = current[i][j];
//                 //     continue;
//                 // }
//                 //
//                 // // proposal[i][j+2] = proposal[i][j];
//                 // mylnL = ln_vraisemblance(MyEpocsData[0].root, MyEpocsData[0].theTVector.tvector[ML_tvector], proposal, IS);
//                 // posterior = log(1)+mylnL;   // If uniform distribution
//                 // acc = acceptance(current_value, posterior);
//                 // if (acc){
//                 //     current[i][j] = proposal[i][j];
//                 //     // w_vec[i][j] = weight_sugg;
//                 //     current_value = mylnL;
//                 //     (*total_acc)++;
//                 // }else{
//                 //     proposal[i][j] = current[i][j];
//                 // }
//             }
//         }
//
//         printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", k, posterior, proposal[0][0], proposal[0][1], proposal[1][0], proposal[1][1], proposal[0][2], proposal[0][3], proposal[1][2], proposal[1][3]);
//         print_rates_vector(current_w_vec);
//     }
//     // print_rates_vector(current_w_vec);
//     // fprintf(stderr, "myvec = ");
//     // for (i=0; i<8; i++){
//     //     fprintf(stderr, "%f, ", my_sum[i]);
//     // }
//     // fprintf(stderr,"\n");
//     print_rates_vector(current_w_vec);
//     return 1;
// }

//

//
// int summary(){
//
// }

