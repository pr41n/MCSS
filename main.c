#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NITER 1e6
#define DIM 2       // Dimension = 2 (circle)

int main(){
    double x[DIM];
    double pi;
    long seed = time(0);
    int count = 0;

    // Display title and seed for random numbers
    printf("####### MONTE CARLO SIMULATION: PI #######\n");
    printf("seed: %ld\n", seed);

    // Use seed for random numbers
    srand(seed);

    // Monte Carlo
    double r2;
    for (unsigned long int i=0; i<NITER; i++) {
        r2 = 0;
        for (int k=0; k<DIM; k++) {
            x[k] = (double) rand()/RAND_MAX * 2 - 1;        // uniform distribution in [-1,1]
            r2 += x[k]*x[k];                                // term of r² = x²
        }
        if (r2 <= 1.0) count++;       // add to count if point inside circle, including border
    }

    pi = (double) count/NITER;
    printf("Number of points: %1.2e\t Inside the circle: %d (%2.2f%%)\n", NITER, count, (float) count/NITER);
    printf("Estimation of pi: %1.5f \n", pi);

    return 0;
}

/*
 long random_at_most(long max) {
    // Code by Ryan Reich:
    // https://stackoverflow.com/questions/2509679/how-to-generate-a-random-integer-number-from-within-a-range
    //

    // Assumes 0 <= max <= RAND_MAX
    // Returns in the closed interval [0, max]
    unsigned long
        // max <= RAND_MAX < ULONG_MAX, so this is okay.
        num_bins = (unsigned long) max + 1,
        num_rand = (unsigned long) RAND_MAX + 1,
        bin_size = num_rand / num_bins,
        defect   = num_rand % num_bins;

    long x;
    do {
        x = random();        // random() has better distribution than rand(), not available in Windows
        printf("%ld ", x);
    }
    while (num_rand - defect <= (unsigned long) x);        // This is carefully written not to overflow

    return x/bin_size;        // Truncated division is intentional
}
 */
