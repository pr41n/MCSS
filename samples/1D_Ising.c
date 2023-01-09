/* Estimated mean energy of a 1D system of spins
 * using Ising model
 * using Monte Carlo dynamic methods:
 *                      metropolis
 *                      heat bath
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int metropolis_update (float kappa, int L, int *lattice, int x);
int heat_bath (float *beta, float *J, int N);

double ener_sum (int L, int *lattice);
double magn_sum (int L, int *lattice);

int *stagger_lattice (int L);
int *up_lattice (int L);
int *down_lattice (int L);
int *random_lattice (int L);

////////// MAIN
int main (int argc, char *argv[]){
    // TODO optional latter as input (txt)

    // INPUT (TODO bof prevention)
    int L;
    float J, beta, niter;

    if (argc == 1 || argc > 5){
        printf("Usage:\t %s L [OPTIONAL] niter J beta \n", argv[0]);
        printf("\tdefault: niter=1E+03, J=1, beta=1\n");
        exit(1);
    }

    sscanf(argv[1], "%d", &L);
    sscanf((argc>=3) ? argv[2] : "1e3", "%f", &niter);        // set input if it exists else default value
    sscanf((argc>=4) ? argv[3] : "1", "%f", &J);
    sscanf((argc==5) ? argv[4] : "1", "%f", &beta);

    if (L%2 != 0){printf("ERROR: L must be even\n"); exit(1);}
    if (beta<=0){printf("ERROR: beta must be positive\n"); exit(1);}


    // BEGIN
    int N = L;
    float kappa = J*beta;
    srand(time(0));

    // Disp info
    printf("iterations: %3.2g\nparticles: %d\nbeta*J = %f ", niter, N, kappa);
    if (kappa - log(1+sqrt(2))/2 < 1e-6) printf("= ");
    else if (kappa < log(1+sqrt(2))/2) printf("< ");
    else printf("> ");
    printf("log(1+sqrt(2))/2 = 0.440687..\n");

    // Try an spsecific lattice
    int *lattice = random_lattice(L);          // May test with stagger_lattice, up_lattice, down_lattice
    double e = -J * ener_sum(L, lattice) / N;
    double m = magn_sum(L, lattice) / N;
    printf("Initial configuration:\te=%4.3f\tm=%4.3f\n", e, m);

    // Calculate energy and magnetization for beta
    // Sequential update
    for (int n=0; n<niter; n++){
        for (int i=0; i<L; i++){
            if (metropolis_update(kappa, L, lattice, i))
                lattice[i] = -lattice[i];

            if (lattice[i]==1)
                printf("+");
            else
                printf("-");
        }
        puts("");
        if (n%10 == 0) {
            e = -J * ener_sum(L, lattice) / N;
            m = magn_sum(L, lattice) / N;
            //printf("iteration %g:\te=%4.3f\tm=%4.3f\n", (float) n, e, m);
        }
    }
    e = -J * ener_sum(L, lattice) / N;
    m = magn_sum(L, lattice) / N;
    //printf("iteration %g:\te=%4.3f\tm=%4.3f\n", niter, e, m);

    return 0;
}
//////////

double ener_sum (int L, int *lattice) {
    /* E = J * sum(s_i*s_j)
     * sum through the closest neighbours of each element of lattice
     * considering only right neighbour (for not repeating)
     * considering periodical lattice
     */
    double E = 0;

    for (int k=0; k<L-1; k++)
        E += lattice[k] * lattice[k+1];
    E += lattice[L-1] * lattice[0];
    return E;
}

double magn_sum (int L, int *lattice) {
    /* M = sum(s_i)
     * sum of each element of lattice
     * considering periodical lattice
     */
    double M = 0;
    for (int k=0; k<L; k++)
        M += lattice[k];
    return M;
}

int metropolis_update (float kappa, int L, int *lattice, int x){
    /* Dynamic MC-method: Metropolis
     */

    /*-------Generally-------*/
    int *lat0 = lattice;
    int lat1[L];

    for (int k=0; k<L; k++)
        lat1[k] = (k!=x) ? lat0[k] : -lat0[k];


    double  P0 = exp(kappa * ener_sum(L, lat0)),
            P1 = exp(kappa * ener_sum(L, lat1));

    double Q = P1/P0;
    /*-------Specifically-------*/
    /*
    diff = lat1[x]*lat1[x+1] + lat1[x]*lat1[x-1] - lat0[x]*lat0[x+1] - lat1[x]*lat1[x-1]    // TODO problems in borders
    double Q = exp(kappa * diff)
    */

    if (Q >= 1) return 1;

    double r = (double) rand()/RAND_MAX;        // random number in [0, 1]      TODO should be in (0,1)
    if (r < Q) return 1;
    else return 0;
}

int *stagger_lattice (int L){
    /* stagger lattice (+ - + - + - ...)
     */
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * L);

    for (int i=0; i<L; i++)
        lat[i] = pow(-1, i%2);

    return lat;
}

int *up_lattice (int L){
    /* spins up lattice (+ + + + + + ...)
     */
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * L);

    for (int i=0; i<L; i++)
        lat[i] = +1;

    return lat;
}

int *down_lattice (int L){
    /* spins down lattice (- - - - - - ...)
     */
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * L);

    for (int i=0; i<L; i++)
        lat[i] = -1;

    return lat;
}

int *random_lattice (int L){
    /* random lattice (+ - + + + - ...)
     */
    int *lat = NULL;
    int r;          // random spin obtained each iteration (+1 or -1)

    lat = realloc(lat, sizeof(int) * L);

    printf("Randomly generated lattice:\n");

    for (int i=0; i<L; i++) {
        r = (rand()%2) * 2 - 1;     // +1 or -1 with equal probability
        lat[i] = r;
        if (r==1) printf("+");
        else printf("-");
    }
    puts("\n");
    return lat;
}
