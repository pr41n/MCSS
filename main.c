/* Estimation of pi from the volume of an n-sphere using the
 * most simple Monte Carlo method
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(int argc, char *argv[]){
    int dim, niter;             //niter up to 2e+9

    /* User input control */
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
            printf("Usage:\t %s ------\n", argv[0]);
            exit(1);
    }

    if (niter <= 0) {printf("Error: niter must be positive\n"); exit(1);}
    if (dim <= 1) {printf("Error: dimension must be >=2\n"); exit(1);}

    /* BEGIN */

    /* END */

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
