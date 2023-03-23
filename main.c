/*Energy and magnetization for a 2D-Ising periodical
  lattice using Wolff algorithm*/

#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<math.h>
#include<time.h>

/*----VARIABLES----*/
/*-----------------*/
/*--Pre-definitions--*/
#define NBLCKDISP 10     // max number of blocks to display

/*--Global variables--*/
/*Lattice*/
int16_t L;                 // length of the lattice (up to 2^14)
uint32_t N;                 // number of particles (L*L)
uint8_t  shift=1;           // log_2(L), at least 1

int8_t *lat;                    // lattice (possible elements -1, +1, 0, -2)
int16_t *rn, *ln, *un, *dn,     // relative positions of neighbours (+1, -1, +L, -L)
        nn[4];                  // array for storing them
uint32_t n;                     // neighbour (0,1,2,...,N)

/*PRNG*/
uint64_t rnd,           // PRNG
         prob_bond;     // probability 1-exp(-2k)
float    frnd;          // normalised rnd 0<frnd<=1 (for Wolff)

/*Physical values*/
float   kappa;          // constant J/kT
double  e,              // observable energy density
        m;              // observable magnetization density

/*Simulation*/
int8_t spin;            // value of the selected spin from the cluster

uint8_t  nblock;        // number of blocks
uint16_t nmeas,         // number of measures per block
         nupdte,        // number of updates between 2 measures
         ntherm;        // number of thermalization updates
uint32_t Ncs,           // number of spins inside the cluster
         Nc;            // number of clusters (for SW)
uint64_t ntotal;        // total number of updates

/*External files*/
FILE    //*in_file,        // input file
        //*bk_file,        // back up file
        *out_file;       // output file

/*--Predefinitions for better performance--*/
/*iterators*/
int16_t x, y;        // coords in lattice
uint32_t i,          // current site in lattice (x+L*y)
         j, k;       // multi-purpose iterators
uint8_t  nit;        // iterator for neighbours (up to 4)
/*precalcs*/
double invN,         // 1/N
       JinvN;        // J/N


/*----FUNCTIONS----*/
/*-----------------*/
/*--PRNG--*/
void xorshift64(){           // generates a PRNG in (0, 2^64-1]
    rnd ^= rnd << 13;
    rnd ^= rnd >> 7;
    rnd ^= rnd << 17;
}
/*--init--*/
void get_data(int argc, char *argv[]){
    /*Get data from input (by now exec args, TODO entry file)*/
    if (argc == 1 || argc > 7){
        printf("Usage:\t %s L output_file [OPTIONAL] nblock nmeas nupdte ntherm kappa\n", argv[0]);
        printf("default: nblock=20, nmeas=1000, nupdte=5, ntherm=10, kappa=0.4406868\n");
        exit(1);
    }

    sscanf(argv[1], "%hd", &L);
    out_file = fopen(argv[2], "w");
    sscanf((argc>=4) ? argv[3] : "20", "%hhu", &nblock);        // set input if it exists elne default value
    sscanf((argc>=5) ? argv[4] : "1000", "%hu", &nmeas);
    sscanf((argc>=6) ? argv[5] : "5", "%hu", &nupdte);
    sscanf((argc>=7) ? argv[6] : "10", "%hu", &ntherm);
    sscanf((argc>=8) ? argv[7] : "0.4406868", "%f", &kappa);


    if (L<=1 || (L&(L-1))){printf("ERROR: L must be 2^n with n>0\n"); exit(1);}
    if (kappa<=0){printf("ERROR: kappa (=J*beta with J=1) must be positive\n"); exit(1);}
}

void setup(){
    /*--Remaining information--*/
    N = L*L;
    while (L>>(shift+1)){
        shift++;
    }
    ntotal = ntherm + nblock*nmeas*nupdte;

    /*--Predefinitions--*/
    int J=1;        //supposing J=1, else it should be input
    invN = (double) 1/N;
    JinvN = (double) J/N;

    /*--PRNG--*/
    prob_bond = (uint64_t) (exp(-2*kappa) * 0x1p+64);   // more precise than 1-exp(-2k)
    rnd = (uint64_t) time(0);

    /*--Lattice--*/
    lat= (int8_t *) malloc(sizeof(char) * N);       // realocating memory
    rn = (int16_t*) malloc(sizeof(short) * L);
    ln = (int16_t*) malloc(sizeof(short) * L);
    un = (int16_t*) malloc(sizeof(short) * L);
    dn = (int16_t*) malloc(sizeof(short) * L);

    for (i=0; i<N; i++){                            // set all spins to +1 or randomly
        xorshift64();                               // some updates to the random number
        lat[i] = +1;  //(rnd&1) * 2 - 1;
    }
    for (i=0; i<L; i++){                            // set relative neighbors
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
}
/*--MC update--*/
void expand_cluster(uint16_t x_cx, uint16_t y_cx, uint32_t i_cx){
    /* we need to create new cluster expansion versions of x,y,i,n,nn variables
     * which will be in different memory locations, one for each function call
     */
    uint8_t nit_cx;
    int16_t nn_cx[4] = {rn[x_cx], ln[x_cx],
                        un[y_cx], dn[y_cx]};
    uint32_t n_cx;

    for (nit_cx=0; nit_cx<4; nit_cx++){         // check all neighbours
        n_cx = i_cx + nn_cx[nit_cx];
        if (lat[n_cx] == spin){         // could be inside the cluster
            xorshift64();
            if (rnd > prob_bond){    // it is!
                Ncs++;
                lat[n_cx] = -spin;
                expand_cluster(n_cx&(L-1), n_cx>>shift ,n_cx);   // now expand the cluster from it
            }
        }
    }
}
void Wolff(){
    xorshift64();                   // choose randomly the spin for the new cluster in (0,N)
    frnd = (float) rnd * 0x1p-64;
    i = (uint32_t) (frnd * N);

    spin = lat[i];      // save spin value for expansion
    lat[i] = -spin;     // flip the spin
    expand_cluster(i&(L-1), i>>shift, i);  // expand the cluster
}
/*--measurements--*/
void measure() {
    static int64_t E, M;    // should be a 64-bit integer since N is 32-bit unsigned integer

    E = M = i = 0;      // reset measurements and site
    for (y=0; y<L; y++) {
        nn[1] = un[y];
        for (x=0; x<L; x++){
            spin = lat[i];
            M += spin;

            nn[0] = rn[x];
            for (nit=0; nit<2; nit++){
                n = i + nn[nit];
                E += spin * lat[n];
            }
            i++;
        }
    }
    e  = (double) -JinvN * E;
    m  = (double) invN * M;
}
/*--output--*/
void disp_lattice(int8_t *lattice) {
    for (j=0; j<N; j++){
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
void disp_init_info() {
    printf("total steps: %lu\nthermalization steps: %d\nmeasures: %d\nparticles: %d\nkappa: %.7f ", ntotal, ntherm, nblock*nmeas, N, kappa);
    if (kappa - log(1+sqrt(2))/2 < 1e-10)    printf("(near critical point ");
    else if (kappa < log(1+sqrt(2))/2)      printf("(below critical point ");
    else                                    printf("(above critical point ");
    printf("log(1+sqrt(2))/2 = 0.4406868..)\n\n");
}

/*----MAIN PROGRAM----*/
/*-----------------*/
int main(int argc, char *argv[]){
    get_data(argc, argv);       // Getting input data
    setup();                    // Setting up and generating the lattice

    uint8_t nbdisp = nblock/NBLCKDISP;
    if (nbdisp == 0) nbdisp = 1;        // nblock < NBLCKDISP, so display them all

    disp_init_info();               // disp some initial info

    /*Thermalization*/
    printf("Beginning thermalization\n");
    clock_t begin_timer = clock();
    for (int nt=0; nt<ntherm; nt++) {
        Wolff();
    }
    clock_t end_timer = clock();
    //disp_lattice(lat);
    printf("Thermalization finished!\n\n");

    /*Measurements*/
    printf("Beginning measures\n");
    float time_spent = (float) (end_timer-begin_timer)/CLOCKS_PER_SEC;
    printf("Estimated time: %5dmin\n", (int) (time_spent/ntherm * (ntotal-ntherm)/60));
    for (int nb=0; nb<nblock; nb++) {
        if (nb%nbdisp == 0 || nb == nblock-1){
            measure();
            printf("%3.0f%%:\te=%4.3f\tm^2=%4.3f\n", (float) 1.25*nb/nblock*100, e, m*m);
            //disp_lattice(lat);
        }
        for (j=0; j<nmeas; j++){
            for (k=0; k<nupdte; k++)
                Wolff();
            measure();
            fprintf(out_file, "\n%6.4f\t%6.4f", e, m);
        }
    }
    printf("Measures finished!\n");
    return 0;
}
