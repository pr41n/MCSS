/* One Monte Carlo update in a single-cluster Wolff algorithm with recursion
 * Made for a 2D-Ising squared lattice, but easily expandable to others 2D Ising lattices
 */
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<time.h>
#include<math.h>

/*-----ALGORITHM-----*/
// Spin s in the lattice may be +1 or -1 if it is outside the cluster
/* If it's inside, we choose it to be -s_old and check possible
 * expansions as s_neighbour == s_old, storing s_old for efficiency
 */

/* For a random spin in the lattice:
 *      create cluster -> expand cluster (recursive function)
 */

int L = 1<<4; // 16, power of 2 (up to 256) is mandatory
float k = 0.4406868;

char *lat;               // lattice
int *rn, *ln, *un, *dn;  // relative positions of neighbours

uint64_t rnd,            // PRNG
         prob_bond;      // probability 1-exp(-2k)
uint32_t i,              // iterator and site in lat iteration
         shift=1,        // log_2(L), at least 1
         Nc;             // number of clusters
float frnd;              // normalised rnd: 0<frnd<=1
char spin;               // spin of the cluster

uint32_t Npn=0, Nn=0;    // number of possible neighbours and number of neighbours


void xorshift64(){           // generates a PRNG in (0, 2^64-1]
    rnd ^= rnd << 13;
    rnd ^= rnd >> 7;
    rnd ^= rnd << 17;
}

void expand_cluster(int site){
    int x = site & (L-1);                           // x,y position in lat
    int y = site >> shift;
    int n;
    int neigh[4] = {site+rn[x], site+ln[x],         // list neighbours
                    site+un[y], site+dn[y]};

    for (int k=0; k<4; k++){         // check all neighbours
        n = neigh[k];
        if (lat[n] == spin){         // could be inside the cluster
            Npn++;
            xorshift64();
            if (rnd > prob_bond){    // it is!
                Nn++;
                lat[n] = -spin;
                expand_cluster(n);   // now expand the cluster from it
            }
        }
    }
}

void disp_lattice(char *lattice) {
    for (int j=0; j<L*L; j++){
        if ((j & (L-1)) == 0) puts("");     // new line
        if (lattice[j]&1){
            if (lattice[j]==1) printf("+");
            else printf("-");
        }
        else {
            if (lattice[j]==0) printf("o");
            else printf(".");
        }
    }
    puts("");
}

int main(){
    /*-----INIT-----*/
    int N = L*L;
    while (L>>(shift+1)){
        shift++;
    }

    rnd = (uint64_t) time(0);
    prob_bond = (uint64_t) (exp(-2*k) * 0x1p+64);        // more precise than 1-exp(-2k)

    lat= (char*) malloc(sizeof(char) * N);       // realocating memory
    rn = (int *) malloc(sizeof(int) * L);
    ln = (int *) malloc(sizeof(int) * L);
    un = (int *) malloc(sizeof(int) * L);
    dn = (int *) malloc(sizeof(int) * L);

    for (i=0; i<N; i++){                           // set all spins to +1 or randomly
        xorshift64();
        lat[i] = +1;  //(rnd&1) * 2 - 1;
    }

    for (i=0; i<L; i++){                          // set relative neighbors
        rn[i] = +1;     // right neighbour
        ln[i] = -1;     // left neighbour
        un[i] = L;      // up neighbour
        dn[i] = -L;     // down neighbour

    }
    /*boundry conditions*/
    ln[0] = L-1;
    rn[L-1] = -ln[0];
    dn[0] = (L-1)*L;
    un[L-1] = -dn[0];

    /*-----SINGLE UPDATE-----*/
    disp_lattice(lat);

    xorshift64();               // choose randomly the spin for the new cluster
    //printf("%lu", rnd);
    frnd = (float) rnd * 0x1p-64;
    i = (uint32_t) (frnd * N);     // i in (0, N)
    printf("\ni=%u\tx=%u\ty=%u\n", i, i & (L-1), i >> shift);

    spin = lat[i];      // save spin value for expansion
    lat[i] = -spin;     // flip the spin
    expand_cluster(i);  // expand the cluster

    disp_lattice(lat);
    printf("\nRatio bonds and possible bonds: %7.2f%%\n", (float) Nn/Npn*100);
}
