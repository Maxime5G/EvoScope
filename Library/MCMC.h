#include "tree.h"
#include <time.h>

// Priors
double exponential_prior(double n, double l);
double normal_prior(double sigma, double mu, double n);

// Sampling
double sample_number(double min, double max);

// MCMC
int MCMC2( struct CoevolData *MyEpocsData, int ntree, double rate[2][4], int *IS, int tvector, int Nrounds, int w, int *total_acc, int sampling );
int acceptance(double current, double proposal);


// OLD

int MCMC( struct CoevolData *MyEpocsData, double *results, int ntree, double rate[2][4], int *IS, int ML_tvector, int Nrounds, int w, int *total_acc );
int MCMC2_NoTVECTOR( struct CoevolData *MyEpocsData, int ntree, double rate[2][4], int *IS, int ML_tvector, int Nrounds, int w, int *total_acc, int sampling );
int MCMC3( struct CoevolData *MyEpocsData, double *results, int ntree, double rate[2][4], int *IS, int ML_tvector, int Nrounds, int w, int *total_acc );
