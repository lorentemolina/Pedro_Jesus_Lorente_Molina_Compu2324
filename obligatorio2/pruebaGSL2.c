#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gsl_rng.h"

#define SEED 1563196755

gsl_rng *tau;

int dim = 30, iter = 200;
float T = 0.01;

int rnd_int() {
    float num = gsl_rng_uniform_int(tau, 2);
    if (num == 0) return -1;
    else if (num == 1) return 1;    
    return 0;
}

int main() {
    extern gsl_rng *tau;
    tau = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, SEED);

    FILE *data;
    data = fopen("ising_data.dat", "w");

    int s[dim][dim];

    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            s[i][j] = rnd_int();
            s[dim-1][dim-1] = s[0][0];
            s[i][dim-1] = s[i][0];
            s[dim-1][j] = s[0][j];
            fprintf(data, "%d\t", s[i][j]);
        }
        if (i != dim-1) fprintf(data, "\n");
    }

    for (int k = 0; k < iter; k++) {
        for (int l = 0; l < dim*dim; l++) {
            int i = gsl_rng_uniform_int(tau, dim);
            int j = gsl_rng_uniform_int(tau, dim);
            int E = 2*s[i][j] * ( s[ ( dim + (i+1) % dim ) % dim ][j] + s[ ( dim + (i-1) % dim) % dim ][j] + s[i][ ( dim + (j+1) % dim) % dim ] + s[i][ ( dim + (j-1) % dim ) % dim ] );
            float p = fminf(1.0, expf(-E/T));
            float r = gsl_rng_uniform(tau);
            if (r < p) {
                s[i][j] *= -1;
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        s[dim-1][dim-1] = s[0][0];
                        s[i][dim-1] = s[i][0];
                        s[dim-1][j] = s[0][j];
                    }
                }
            }
        }
        for (int i = 0; i < dim; i++) {
            fprintf(data, "\n");
            for (int j = 0; j < dim; j++) {
                fprintf(data, "%d\t", s[i][j]);
            }
        }
    }

    fclose(data);
    return 0;
}
