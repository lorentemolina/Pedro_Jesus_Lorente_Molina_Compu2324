#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define dim 30 // Dimensiones del retículo
#define iter 200 // Número de iteraciones
#define T 0.01

int rnd_int() {
    return rand() % 2 == 0 ? 1 : -1;
}

int main() {
    srand(time(NULL)); // Semilla aleatoria

    FILE *data;
    data = fopen("ising_data.dat", "w");

    int s[dim][dim];

    // Generar configuración inicial del retículo
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            s[i][j] = rnd_int();
            fprintf(data, "%d", s[i][j]); // Imprimir el primer elemento sin coma
            if (j < dim - 1) { // Si no es el último elemento de la fila
                fprintf(data, ","); // Imprimir una coma
            }
        }
        fprintf(data, "\n");
    }
    fprintf(data, "\n");

    // Realizar las iteraciones
    for (int k = 0; k < iter; k++) {
        for (int l = 0; l < dim*dim; l++) {
            int i = rand() % dim;
            int j = rand() % dim;
            int E = 2 * s[i][j] * (s[(dim + (i+1) % dim) % dim][j] + s[(dim + (i-1) % dim) % dim][j] + s[i][(dim + (j+1) % dim) % dim] + s[i][(dim + (j-1) % dim) % dim]);
            float p = fminf(1.0, exp(-E/T));
            float r = (float)rand() / RAND_MAX;
            if (r < p) {
                s[i][j] *= -1;
            }
        }

        // Escribir los datos en el archivo después de cada iteración
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                fprintf(data, "%d", s[i][j]); // Imprimir el primer elemento sin coma
                if (j < dim - 1) { // Si no es el último elemento de la fila
                    fprintf(data, ","); // Imprimir una coma
                }
            }
            fprintf(data, "\n");
        }
        fprintf(data, "\n");
    }

    fclose(data);
    return 0;
}
