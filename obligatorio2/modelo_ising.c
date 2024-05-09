#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#define dim 80 // Dimensiones del retículo
#define iter 200 // Número de iteraciones
#define T 2.2 //Temperatura

int rnd_int() {
    return rand() % 2 == 0 ? 1 : -1;  // Para generar aleatoriamente 1 o -1
}

int main() {

    // Para observar los tiempos de compilación
    clock_t start_time, end_time;
    double total_time;

    // Marcamos el tiempo de inicio
    start_time = clock();

    // Semilla "aleatoria"
    srand(time(NULL)); 

    // Abrimos ficheros
    FILE *data = fopen("ising_data.dat", "w"); // para guardar los resultados
    FILE *time = fopen("tiempo_ejecucion.txt", "w"); // para obtener el tiempo que tarda en ejecutarse


    int s[dim][dim];

    // Genera la configuración inicial de la matriz
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            s[i][j] = rnd_int();
            fprintf(data, "%d", s[i][j]); // Imprimimos el primer elemento sin coma
            if (j < dim - 1) { // Si no es el último elemento de la fila
                fprintf(data, ","); // Imprime una coma
            }
        }
        fprintf(data, "\n");
    }
    fprintf(data, "\n");

    // Realizamos las iteraciones
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

        // Escribimos los datos en el archivo después de cada iteración
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                fprintf(data, "%d", s[i][j]); // Imprimimos el primer elemento sin coma
                if (j < dim - 1) { // Si no es el último elemento de la fila
                    fprintf(data, ","); // Imprime una coma
                }
            }
            fprintf(data, "\n");
        }
        fprintf(data, "\n");
    }

    fclose(data);

     // Marcamos el tiempo de finalización
    end_time = clock();

     // Calculamos el tiempo total de compilación
    total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Imprimimos el tiempo total de compilación
    fprintf(time, "dim: %-4d Tiempo de ejecución: %.4f\n", dim, total_time);

    fclose(time);

    return 0;
}
