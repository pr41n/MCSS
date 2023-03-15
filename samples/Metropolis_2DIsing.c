/*Energy and magnetization for 2D periodical
  lattices using Metropolis method and Ising model*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/*----VARIABLES----*/
/*-----------------*/
/*--Pre-definitions--*/
#define NDISP 20        // number of points to display (every 5%)
#define NSTORE 1e7      // maximum number of points to store

/*--Global variables--*/        // iterators?
int Lx,                  // length of the lattice
    Ly,                  // height of the lattice
    N;                   // number of particles (Lx*Ly)

int *lat = NULL,         // lattice                                                 TODO as char
    //*ns[4] = {1, -1, Lx, -Lx},              // neighbour spins ({1, -1, L, -L} for 2D squared lattice)
    *rn, *ln, *un, *dn;  // relative positions of neighbours

float   J = 1,
        beta = 0.44,     // 1/kT
        kappa,           // J*beta
        nther=100,        // number of iterations for thermalization
        //nblks,           // number of blocks of measurements
        nmeas;           // number of iterations for measurements

double  e,               // energy density
        m,               // magnetization density
        prob[5];         // Metropolis probability table for exponents -8, -4, 0 (++--), +4 (+++-), +8 (++++)

FILE    //*in_file,        // input file
        //*bk_file,        // back up file
        *out_file;       // output file

/*--Pre-definitions for speeding up--*/
int i, j, n;        //iterators
double invN,         // 1/N
       JinvN;        // J/N


/*----FUNCTIONS----*/
/*-----------------*/
/*--init--*/
void get_data(int argc, char *argv[]){
    /*Get data from input (by now exec args, TODO entry file)*/
    if (argc == 1 || argc > 7){
        printf("Usage:\t %s Lx Ly output_file [OPTIONAL] niter J beta\n", argv[0]);
        printf("\tdefault: niter=1E+03, J=1, beta=0.44\n");
        exit(1);
    }

    sscanf(argv[1], "%d", &Lx);
    sscanf(argv[2], "%d", &Ly);
    out_file = fopen(argv[3], "w");
    sscanf((argc>=5) ? argv[4] : "1e3", "%f", &nmeas);        // set input if it exists elne default value
    sscanf((argc>=6) ? argv[5] : "1", "%f", &J);
    sscanf((argc==7) ? argv[6] : "0.44", "%f", &beta);

    if (Lx <= 1 || Ly <= 1){printf("ERROR: L must be >=2\n"); exit(1);}
    if (beta<=0){printf("ERROR: beta must be positive\n"); exit(1);}
}
void setup(){
    /*--Remaining global variables--*/
    N = Lx*Ly;
    kappa = J*beta;
    invN = (double) 1/N;
    JinvN = (double) J/N;

    for (i=0; i<5; i++)
        prob[i] = exp(kappa * (-8+4*i));

    /*--PRNG--*/
    srand(time(0));

    /*--Lattice--*/                            // TODO char, better PRNG
    lat = realloc(lat, sizeof(int) * N);       // realocating memory
    for (i=0; i<N; i++)
        lat[i] = (rand()&1) * 2 - 1;

    rn = realloc(rn, sizeof(int) * Lx);
    ln = realloc(ln, sizeof(int) * Lx);
    un = realloc(un, sizeof(int) * Ly);
    dn = realloc(dn, sizeof(int) * Ly);

    for (i=0; i<Lx && i<Ly; i++){
        if (i < Lx-1){
            rn[i] = +1;
            ln[i] = -1;
        }
        if (i < Ly-1){
            un[i] = Lx;
            dn[i] = -Lx;
        }
    }
    /*boundry conditions*/
    ln[0] = Lx-1;
    rn[Lx-1] = 0;
    dn[0] = (Ly-1)*Lx;
    un[Ly-1] = -dn[0];
}
/*--Monte Carlo update--*/
void metropolis_update(){
    static int  site_u, ss,   // site in lat, sum of neighbour spins
                py_u, my_u;   // TODO y-axis position of upper and lower neighbours

    static double Q;       // quotient Q=exp( -2*s_i * sum_neigh(s_j) )

    site_u = 0;
    for (j=0; j<Ly; j++) {
    for (i=0; i<Lx; i++) {
        ss = lat[site_u + rn[i]] +
             lat[site_u + ln[i]] +
             lat[site_u + un[j]] +
             lat[site_u + dn[j]];
        ss *= -lat[site_u];        // removing 2* for simplification

        Q = prob[(ss+4)>>1];       // {-8,-4,0,+4,+8} -> {0,1,2,3,4}   (simplification: 2x -> (2x+8)/4 = (x+4)/2)
        if (Q>=1)
            lat[site_u] = -lat[site_u];
        else if ((double) rand()/RAND_MAX < Q)
            lat[site_u] = -lat[site_u];

        site_u++;
    }
    }
}
/*--measurements--*/
void measure() {
    static int site_m, py_m;      // site in lat, TODO y-axis position of upper neighbour
    static double E, M;

    /*reset measurements*/
    E = 0;
    M = 0;
    site_m = 0;
    /*iterate lattice*/
    for (j=0; j<Ly; j++) {
        for (i=0; i<Lx; i++) {
            M += lat[site_m];
            E += lat[site_m] * ( lat[site_m + rn[i]] + lat[site_m + un[j]]);
            site_m++;
        }
    }

    e = (double) -JinvN * E;
    m = (double) invN * M;
}
/*--output--*/
void disp_lattice(int *la) {
    int y;
    for (j=0; j<Ly; j++){
        y = j*Lx;
        for (i=0; i<Lx; i++){
            if (la[y+i]==1) printf("+");
            else printf("-");
        }
        puts("");
    }
    puts("");       // new line
}
void disp_init_info() {
    printf("iterations: %3.2g\nparticles: %d\nbeta*J = %f ", nmeas, N, kappa);
    if (kappa - log(1+sqrt(2))/2 < 1e-6)    printf("= ");
    else if (kappa < log(1+sqrt(2))/2)      printf("< ");
    else                                    printf("> ");
    printf("log(1+sqrt(2))/2 = 0.440687..\n");
}

/*----PROGRAM----*/
/*-----------------*/
int main(int argc, char *argv[]){
    get_data(argc, argv);       // Getting input data
    setup();                    // Setting up and generating the lattice

    int ndisp = (int) nmeas/NDISP,      // When to display and store
        nstore = (int) nmeas/NSTORE;
    if (ndisp==0)  ndisp=1;
    if (nstore==0) nstore=1;

    disp_init_info();               // disp some initial info
    measure();
    printf("Initial lattice:\ne=%4.3f\tm=%4.3f\n", e, m);
    disp_lattice(lat);

    /*Thermalization*/
    for (n=0; n<nther; n++) {
        metropolis_update();
    }

    /*Measurements*/
    for (n=0; n<nmeas; n++) {
        if (n%nstore == 0) {
            measure();
            fprintf(out_file, "\n%6.4f\t%6.4f", e, m);
            if (n%ndisp == 0 || n == nmeas-1){
                printf("%3.0f%%:\te=%4.3f\tm=%4.3f\n", (float) n/nmeas*100, e, m);
                disp_lattice(lat);
            }
        }
        metropolis_update();
    }
    return 0;
}


