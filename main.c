#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SEED time(0)

long random_at_most(long max) {
    /* Code by Ryan Reich: 
    * https://stackoverflow.com/questions/2509679/how-to-generate-a-random-integer-number-from-within-a-range
    */ 
    
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
int main(){
    printf("Monte Carlo simulation");
    return 0;
}
