#include <stdio.h>
#include <math.h>
#include <time.h>
#include "complex.h"

#define PI 3.14159265358979323846
#define H 6.62607004e-34
#define MAX 1000

int iter, N, n_ciclos;
double lambda, k, s, norm = 0;

double V[MAX]; // Potencial
fcomplex fonda[MAX]; // Función de onda
fcomplex A[MAX], b[MAX], alpha[MAX], beta[MAX], chi[MAX]; // Coeficientes

int main() {

    // Para observar los tiempos de compilación

    clock_t start_time, end_time;
    double total_time;

    // Marcamos el tiempo de inicio
    start_time = clock();

    // Abrimos los ficheros de salida

    FILE *poten;
    poten = fopen("potencial.dat", "w"); // potencial
    FILE *modfonda;
    modfonda = fopen("funcion_de_onda.dat", "w"); // módulo de la función de onda
    FILE *normfonda;
    normfonda = fopen("norma.dat", "w"); // norma de la función de onda
    FILE *time = fopen("tiempo_ejecucion.txt", "w"); // para obtener el tiempo que tarda en ejecutarse

    // Definimos los parámetros iniciales

    iter = 10000; // Nº de iteraciones temporales
    N = 500; // Divisiones del eje x
    n_ciclos = 100; // Ciclos que se realizan
    lambda = 0.8; // Parámetro lambda (Altura del potencial)
    k = 2 * PI * (double)n_ciclos / (double)N; // Factor k reescalado
    s = 1 / (4 * k * k); // Espaciado temporal reescalado
    fcomplex i = Complex(0, 1); // imaginario i

    // Cálculo del potencial

    for (int j = 0; j < N; j++) {
        if ((j < 3 * N / 5) && (j > 2 * N / 5)) V[j] = lambda * k * k; // Potencial de altura lambda*k*k
        else V[j] = 0.0; // Potencial nulo fuera del rango definido
        fprintf(poten, "%.15lf\n", V[j]); // Escritura en un fichero del potencial
    }

    // Función de onda inicial

    // Cálculo de la función de onda

    for (int j = 0; j < N; j++) {
        if (j == 0 || j == (N)) fonda[j] = Complex(0, 0); // Condiciones de contorno
        else {
            fonda[j] = Complex(cos(k * j) * exp((-8 * (4 * j - N) * (4 * j - N)) / (double)(N * N)),
                               sin(k * j) * exp((-8 * (4 * j - N) * (4 * j - N)) / (double)(N * N))); // se ha escrito el primer exponencial de la forma cos(x)+isen(x) y se multiplica
            norm += Cabs(Cmul(fonda[j], Conjg(fonda[j]))); // Cálculo de la norma de la función de onda
        }
    }

    // Normalización

    for (int j = 0; j < N; j++) {
        fonda[j] = Cdiv(fonda[j], Complex(sqrt((norm)), 0)); // dividiendo por la norma
        fprintf(modfonda, "%.15lf\n", Cabs(fonda[j])); // Escritura de la función de onda normalizada en cada punto j del eje x
    }
    fprintf(modfonda, "\n"); // para separar iteraciones

    // Cálculo de los vectores "A_j^0" y "alpha"

    for (int j = 0; j < N; j++) A[j] = Complex(-2 - V[j], 2 / s); // Vector A_j^0
    alpha[N - 1] = Complex(0.0, 0.0); // Condición de contorno alpha_N-1 = 0
    for (int j = N - 2; j > 0; j--) alpha[j - 1] = Cdiv(Complex(-1.0, 0.0), Cadd(A[j], alpha[j])); // Vector alpha_j

    // Aquí comienza el Algoritmo

    for (int t = 0; t < iter; t++) {

        norm = 0; // Inicializamos la norma en cada iteración

        // Paso 1: Cálculo de los vectores "b_j" y "beta"

        for (int j = 0; j < N; j++) b[j] = Cdiv(Cmul(Complex(0.0, 4.0), fonda[j]), Complex(s, 0.0)); // Vector b_j
        beta[N - 1] = Complex(0.0, 0.0); // Condición de contorno beta_N-1 = 0
        for (int j = N - 2; j > 0; j--) beta[j - 1] = Cdiv(Csub(b[j], beta[j]), Cadd(A[j], alpha[j])); // Vector beta_j

        // Paso 2: Cálculo del vector "chi"

        for (int j = 0; j < N - 1; j++) chi[j + 1] = Cadd(Cmul(alpha[j], chi[j]), beta[j]);

        // Paso 3: Cálculo de la función de onda

        for (int j = 0; j < N; j++) {
            fonda[j] = Csub(chi[j], fonda[j]);
            norm += Cabs(Cmul(fonda[j], Conjg(fonda[j])));
        }

        fprintf(normfonda, "%d\t%.15lf\n", t, sqrt(norm)); // Escritura en un fichero de la raíz de la norma al cuadrado para cada iteración
        for (int j = 0; j < N; j++) fprintf(modfonda, "%.15lf\n", Cabs(fonda[j])); // Escritura en un fichero del módulo de la función de onda
        fprintf(modfonda, "\n"); // Separación de las funciones de onda

    }

    fclose(modfonda);
    fclose(poten);
    fclose(normfonda);

    // Marcamos el tiempo de finalización
    end_time = clock();

     // Calculamos el tiempo total de compilación
    total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Imprimimos el tiempo total de compilación
    fprintf(time, "%.4f\n", total_time);

    fclose(time);

    return 0;
}
