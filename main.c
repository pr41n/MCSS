/* Estimated mean energy and magnetization of a 2D system of spins
 * using Ising model
 * using Monte Carlo dynamic methods:
 *                      metropolis
 *                      heat bath
*/

/*Differences with 1D_Ising: N, for_loop, lattice_def - in main
                             ener_sum, magn_sum,
                             metropolis_update (L -> L*L)
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int metropolis_update (float kappa, int L, int *lattice, int x);
int heat_bath (float *beta, float *J, int N);

double ener_sum (int L, int *lat);
double magn_sum (int L, int *lat);

int *stagger_lattice (int L);
int *up_lattice (int L);
int *down_lattice (int L);
int *random_lattice (int L);

void disp_lattice(int L, int *lat);

////////// MAIN
int main (int argc, char *argv[]){
    // TODO optional lattice as input (txt)
    // TODO input args using commands (-beta, -J, -o)

    /*INPUT*/
    int L;
    float J, beta, niter;

    if (argc == 1 || argc > 6){
        printf("Usage:\t %s L output_file [OPTIONAL] niter J beta \n", argv[0]);
        printf("\tdefault: niter=1E+03, J=1, beta=0.44\n");
        exit(1);
    }

    sscanf(argv[1], "%d", &L);
    FILE *file = fopen(argv[2], "w");
    sscanf((argc>=4) ? argv[3] : "1e3", "%f", &niter);        // set input if it exists else default value
    sscanf((argc>=5) ? argv[4] : "1", "%f", &J);
    sscanf((argc==6) ? argv[5] : "0.44", "%f", &beta);

    if (L%2 != 0){printf("ERROR: L must be even\n"); exit(1);}
    if (L<2){printf("ERROR: L must be >=2\n"); exit(1);}
    if (beta<=0){printf("ERROR: beta must be positive\n"); exit(1);}


    // BEGIN
    int N = L*L;
    float invN = (float) 1/N;
    float kappa = J*beta;

    int ndisp = (int) niter/10,
        nstore = (int) niter/1e7;

    if (ndisp==0)  ndisp=1;
    if (nstore==0) nstore=1;

    srand(time(0));

    // Disp info
    printf("iterations: %3.2g\nparticles: %d\nbeta*J = %f\t(", niter, N, kappa);
    if (kappa - log(1+sqrt(2))/2 < 1e-6) printf("= ");
    else if (kappa < log(1+sqrt(2))/2) printf("< ");
    else printf("> ");
    printf("log(1+sqrt(2))/2 = 0.440687..)\n\n");

    // Generating initial lattice
    int *lattice = random_lattice(L);          // May test with stagger_lattice, up_lattice, down_lattice
    double e = -J * ener_sum(L, lattice) * invN;
    double m = magn_sum(L, lattice) * invN;

    //printf("Initial configuration:\te=%4.3f\tm=%4.3f\n", e, m);
    //disp_lattice(L, lattice);
    //puts("\n");

    // Calculate energy and magnetization for beta
    // Sequential update
    for (int n=0; n<niter; n++){
        if (n%nstore == 0) {
            e = -J * ener_sum(L, lattice) * invN;
            m = magn_sum(L, lattice) * invN;
            fprintf(file, "\n%.8f\t%.8f", e, m);
        }
        if (n%ndisp == 0) {
            e = -J * ener_sum(L, lattice) * invN;
            m = magn_sum(L, lattice) * invN;
            printf("%3.0f%%:\te=%4.3f\tm=%4.3f\n", (float) n/niter*100, e, m);
            disp_lattice(L, lattice);
        }

        for (int k=0; k<N; k++){
            if (metropolis_update(kappa, L, lattice, k))
                lattice[k] = -lattice[k];
        }
    }
    e = -J * ener_sum(L, lattice) * invN;
    m = magn_sum(L, lattice) * invN;
    printf("100%%:\te=%4.3f\tm=%4.3f\n", e, m);
    disp_lattice(L, lattice);

    return 0;
}
//////////

double ener_sum (int L, int *lat) {
    /* E = J * sum(s_i*s_j)
     * sum through the closest neighbours of each element of lattice
     * considering only right and up neighbour (for not repeating)
     * assuming periodical lattice
     */
    double E = 0;
    /*Forward sum*/
    for (int i=0; i<L-1; i++){
    for (int j=0; j<L-1; j++)
        E += lat[i*L + j] * (lat[i*L + j+1] + lat[(i+1)*L + j]);
    E += lat[i*L + L-1] * (lat[i*L + 0] + lat[(i+1)*L + L-1]);
    }
    for (int j=0; j<L-1; j++)
        E += lat[(L-1)*L + j] * (lat[(L-1)*L + j+1] + lat[0 + j]);
    E += lat[(L-1)*L + L-1] * (lat[0 + L-1] + lat[(L-1)*L + L-2]);

    /* Backward sum

    */
    return E;
}

double magn_sum (int L, int *lat) {
    /* M = sum(s_i)
     * sum of each element of lattice
     * assuming periodical lattice
     */
    double M = 0;
    for (int i=0; i<L; i++){
    for (int j=0; j<L; j++)
        M += lat[i*L + j];
    }
    return M;
}

int metropolis_update (float kappa, int L, int *lattice, int x){
    /* Dynamic MC-method: Metropolis
     * returns boolean value whether to change or not
     */
    int *lat0 = lattice;

    /*-------Generally-------*//*
    int lat1[L*L];

    for (int k=0; k<L*L; k++)
        lat1[k] = (k!=x) ? lat0[k] : -lat0[k];

    double  P0 = exp(kappa * ener_sum(L, lat0)),
            P1 = exp(kappa * ener_sum(L, lat1));

    double Q = P1/P0;*/
    /*-------Specifically-------*/
    int s_up =    (x >= L*(L-1)) ? lat0[x%L]       : lat0[x+L],
        s_down =  (x <= L-1)     ? lat0[x+L*(L-1)] : lat0[x-L],
        s_left =  (x%L == 0)     ? lat0[x+L-1]     : lat0[x-1],
        s_right = (x%L == L-1)   ? lat0[x-L+1]     : lat0[x+1];

    double diff = -2*lat0[x] * (s_up + s_down + s_left + s_right);
    double Q = exp(kappa * diff);

    /*-------Always-------*/
    if (Q >= 1) return 1;

    double r = (double) rand()/RAND_MAX;        // random number in [0, 1]      TODO should be in (0,1)
    if (r < Q) return 1;
    else return 0;
}

int *stagger_lattice (int L){
    /* stagger lattice (+ - + - + - ...
     *                  - + - + - + ...
     *                  ...)
     */
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * L*L);

    for (int i=0; i<L; i++){
    for (int j=0; j<L; j++)
        lat[i*L + j] = pow(-1, i) * pow(-1, j%2);
    }

    return lat;
}

int *up_lattice (int L){
    /* spins up lattice (+ + + + + + ...
    *                    + + + + + + ...
    *                    ...)
     */
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * L*L);

    for (int i=0; i<L; i++){
    for (int j=0; j<L; j++)
        lat[i*L + j] = 1;
    }

    return lat;
}

int *down_lattice (int L){
    /* spins down lattice (- - - - - - ...
    *                      - - - - - - ...
    *                      ...)
     */
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * L*L);

    for (int i=0; i<L; i++){
    for (int j=0; j<L; j++)
        lat[i*L + j] = -1;
    }

    return lat;
}

int *random_lattice (int L){
    /* random lattice (+ - + + + - ...
     *                 - + - - + - ...
     *                 ...)
     */
    int r;          // random spin obtained each iteration (+1 or -1)
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * L*L);

    for (int i=0; i<L; i++) {
        for (int j=0; j<L; j++){
            r = (rand()%2) * 2 - 1;     // +1 or -1 with equal probability
            lat[i*L + j] = r;
        }
    }
    return lat;
}

void disp_lattice(int L, int *lat) {
    for (int i=0; i<L; i++){
        for (int j=0; j<L; j++){
            if (lat[i*L+j]==1) printf("+");
            else printf("-");
        }
    puts("");
    }
}
