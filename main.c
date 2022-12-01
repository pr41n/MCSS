#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// TODO error estimation

double volume_factor (int n);

int main(int argc, char *argv[]){
    int dim, niter;             //niter up to 2e+9
    switch (argc) {
        case 2:
            sscanf(argv[1], "%d", &niter);
            dim = 2;
            break;
        case 3:
            sscanf(argv[1], "%d", &niter);
            sscanf(argv[2], "%d", &dim);
            break;
        default:
            printf("Usage:\t %s niter [OPTIONAL] dim\n", argv[0]);
            exit(1);
    }

    if (niter <= 0) {printf("niter must be positive\n"); exit(1);}
    if (dim <= 1) {printf("dimension must be >=2\n"); exit(1);}

    double x[dim];
    double pi;
    long seed = time(0);
    int count = 0;

    // Display title and seed for random numbers
    printf("####### MONTE CARLO SIMULATION: PI #######\n");
    printf("niter: %g\tdimension: %d\tseed: %ld\n", (float) niter, dim, seed);

    // Use seed for random numbers
    srand(seed);

    // Monte Carlo
    double r2;
    for (int i=0; i<niter; i++) {
        r2 = 0;
        for (int k=0; k<dim; k++) {
            x[k] = (double) rand()/RAND_MAX * 2 - 1;        // uniform distribution in [-1,1]
            r2 += x[k]*x[k];                                // term of r² = x²
        }
        if (r2 <= 1.0) count++;       // add to count if point inside circle, including border
    }

    pi = (double) powf((double) volume_factor(dim) * count/niter * pow(2,dim), (double)1 / ((int)dim/2));
    printf("Number of points: %1.2e\t Inside: %d (%2.2f%%)\n", (float) niter, count, (float) count/niter*100);
    printf("Estimation of pi: %1.5f \n", pi);

    return 0;
}

double volume_factor (int n) {
    if (n<0){printf("Error -> Volume factor: n must be positive\n"); exit(1);}

    switch (n) {
        case 1: return 2;
        case 2: return 1;
        default:
            return (float) 2/n * volume_factor(n-2);
    }
}
