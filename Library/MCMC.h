#include "tree.h"
#include "coevol.h"
#include <time.h>

// Priors
double exponential_prior(double n, double l);
double normal_prior(double sigma, double mu, double n);

// Sampling
double sample_number(double min, double max);

// MCMC
int MCMC2( struct CoevolData *MyEpocsData, int *IS, int Nrounds, int w, int *total_acc, int sampling, FILE *fptr );
int acceptance(double current, double proposal);
