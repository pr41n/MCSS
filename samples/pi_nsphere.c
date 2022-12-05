/*
 * Estimation of pi from the volume of an n-sphere using the
 * most simple Monte Carlo method
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


double volume_factor (int n);

int main(int argc, char *argv[]){

    int dim, niter;             // niter up to 2e+9

    /* User input control */
    switch (argc) {                 ////////////
        case 2:
            sscanf(argv[1], "%d", &niter);
            dim = 2;
            break;
        case 3:
            sscanf(argv[1], "%d", &niter);
            sscanf(argv[2], "%d", &dim);
            break;
        default:
            printf("Usage:\t %s niter [OPTIONAL] dim\nDefault dim=2\n", argv[0]);
            exit(1);
    }

    if (niter <= 0) {printf("Error: niter must be positive\n"); exit(1);}
    if (dim <= 1) {printf("Error: dimension must be >=2\n"); exit(1);}

    /* BEGIN */
    double x[dim];
    double pi;
    long seed = time(0);
    int count = 0;

    /* Display title and seed for random numbers */
    printf("####### MONTE CARLO SIMULATION: PI #######\n");
    printf("niter: %g\tdimension: %d\tseed: %ld\n", (float) niter, dim, seed);

    /* Monte Carlo */
    srand(seed);

    double r2;
    for (int i=0; i<niter; i++) {
        r2 = 0;
        for (int k=0; k<dim; k++) {
            x[k] = (double) rand()/RAND_MAX * 2 - 1;        // uniform distribution in [-1,1]
            r2 += x[k]*x[k];                                // term of r² = x²
        }
        if (r2 <= 1.0) count++;       // add to count if point inside circle, including border
    }

    pi = (double) powf((double) 1/volume_factor(dim) * count/niter * pow(2,dim), (double)1 / ((int)dim/2));

    /* Error estimation */
    double  exact_prob = M_PI/4,
            empirical_prob = (double) count/niter,
            exact_err = pow(2,dim) * sqrt(exact_prob * (1-exact_prob) / niter),                 // only valid for dim=2 !!!!!
            empirical_err = pow(2, dim) * sqrt(empirical_prob * (1-empirical_prob) / (niter-1));        // only valid for dim=2 !!!!!

    /* END */

    printf("Number of points: %.2e\t Inside: %d (%.2f%%)\n", (float) niter, count, (float) count/niter*100);
    printf("Estimation of pi: %.5f +- %.5f (%.5f)\n", pi, exact_err, empirical_err);
    printf("Deviation: %.5g sigmas\n", (pi-M_PI)/exact_err);
    //printf("Volume of %d-sphere: %1.5f \n", dim, volume_factor(dim)*pi);

    return 0;
}

double volume_factor (int n) {
    /* Volume of a n-sphere with R=1 in units of pi
     * V_n = 2*pi/n * V_{n-2}; V_1=2; V_2=pi=1
     */

    if (n<0){printf("Error -> Volume factor: n must be positive\n"); exit(1);}

    switch (n) {
        case 1: return 2;
        case 2: return 1;       // V_2 = pi = 1
        default:
            return (float) 2/n * volume_factor(n-2);
    }
}
