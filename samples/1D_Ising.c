/* Estimated mean energy of a 1D system of spins
 * using Ising model
 * using Monte Carlo dynamic methods:
 *                      metropolis
 *                      heat bath
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double metropolis (float *beta, float *J, int N);

double heat_bath (float *beta, float *J, int N);

double energy (float *J, int *L, int *lattice);

double magnetization (int *L, int *lattice);

int *stagger_lattice (int *L);
int *up_lattice (int *L);
int *down_lattice (int *L);


int main (int argc, char *argv[]){

    int L = 2*2;  //Check if even, make them input
    float J = 1,
          beta = 1;

    int N = L;

    int *lattice = stagger_lattice(&L);

    double e = energy(&J, &L, lattice) / N;
    double m = magnetization(&L, lattice) / N;

    printf("e=%4.3f\tm=%4.3f", e, m);

    return 0;
}

double energy (float *J, int *L, int *lattice) {
    /* E = J * sum(s_i*s_j)
     * sum through the closest neighbours of each element of lattice
     * considering periodical lattice
     */

    double E = 0;

    for (int k=0; k<*L-1; k++)
        E += lattice[k] * lattice[k+1];

    E += lattice[*L-1] * lattice[0];
    E *= -(*J);

    return E;
}

double magnetization (int *L, int *lattice) {
    /* M = sum(s_i)
     * sum of each element of lattice
     * considering periodical lattice
     */

    double M = 0;

    for (int k=0; k<*L; k++)
        M += lattice[k];

    return M;
}

int *stagger_lattice (int *L){
    /* stagger lattice (+ - + - + - ...)
     */
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * (*L));

    for (int i=0; i<*L; i++)
        lat[i] = pow(-1, i%2);

    return lat;
}

int *up_lattice (int *L){
    /* spins up lattice (+ + + + + + ...)
     */
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * (*L));

    for (int i=0; i<*L; i++)
        lat[i] = +1;

    return lat;
}

int *down_lattice (int *L){
    /* spins down lattice (- - - - - - ...)
     */
    int *lat = NULL;
    lat = realloc(lat, sizeof(int) * (*L));

    for (int i=0; i<*L; i++)
        lat[i] = -1;

    return lat;
}
