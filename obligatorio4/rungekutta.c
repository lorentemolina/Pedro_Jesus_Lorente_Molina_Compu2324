#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// Definimos las constantes
#define PI 3.14159265358979323846
#define G 6.674e-11
#define RT 6.378160e6 // radio de la Tierra en metros
#define RL 1.7374e6 // radio de la Luna en metros
#define MT 5.9736e24 // Masa de la Tierra en kg
#define ML 0.07349e24 // Masa de la Luna en kg
#define DTL 3.844e8 // radio de la órbita Tierra-Luna en m
#define OMEGA 2.6617e-6 // periodo de la órbita de la Luna en s
#define MU 0.0123025 // valor mu
#define DELTA 7.01474e-12 // valor delta

// Definimos los parámetros de la simulación
#define ITER 1e6 // iteraciones
#define H 1 // paso
#define MAX 100 // tamaño máx

// Definimos las variables
double k1[MAX], k2[MAX], k3[MAX], k4[MAX]; // Vectores k_n^i donde se asigna un valor para cada ecuación diferencial
double y[MAX]; // Vector y(t) con las coordenadas y momentos generalizados del cohete
int paso = 0; // Para comprobar si el cohete ha pasado por el eje x a la Luna

// Función que devuelve la ecuación diferencial correspondiente
double f(int n, double r, double phi, double pr, double pphi, double t) {
    double rprima = sqrt(1 + r * r - 2 * r * cos(phi - OMEGA * t));
    if (n == 1) return pr;
    if (n == 2) return pphi / (r * r);
    if (n == 3) return (pphi * pphi) / (r * r * r) - DELTA * (1.0 / (r * r) + (MU / (rprima * rprima * rprima)) * (r - cos(phi - OMEGA * t)));
    if (n == 4) return -(DELTA * MU * r / (rprima * rprima * rprima)) * sin(phi - OMEGA * t);
    return 0;
}

// Función para calcular la distancia de la nave a la Luna
double calcular_rL(double r, double phi, double t) {
    return sqrt(r * r + DTL * DTL - 2 * r * DTL * cos(phi - OMEGA * t));
}

// Función auxiliar para calcular el Hamiltoniano H
double calcular_Ham(double r, double pr, double pphi, double rL) {
    double T = (pr * pr) / 2.0 + (pphi * pphi) / (2.0 * r * r);
    double V = -G * MT / r - G * ML / rL;
    return T + V;
}

int main() {

    // Para observar los tiempos de compilación
    
    clock_t start_time, end_time;
    double total_time;

    // Marcamos el tiempo de inicio

    start_time = clock();

    // Damos los valores iniciales
    double theta = PI / 3.387; // ángulo inicial de la velocidad respecto al eje x
    double v = 11200.0 / DTL; // velocidad inicial del lanzamiento del cohete en m/s
    double Ham = 0.0; // H
    double Ham_prima = 0.0; // H'

    // Abrimos los ficheros de salida
    FILE *out = fopen("rungekutta.dat", "w"); // para almacenar las trayectorias
    FILE *time = fopen("tiempo_ejecucion.txt", "w"); // para mostrar el tiempo que tarda en ejecutarse
    FILE *cons = fopen("cons_ham.dat", "w"); // para observar que se conserva H' con el tiempo

    // Asignamos los valores iniciales al vector y[n]
    y[1] = RT / DTL; // radio
    y[2] = PI / 2.0; // ángulo
    y[3] = v * cos(theta - y[2]); // momento radial
    y[4] = y[1] * v * sin(theta - y[2]); // momento angular

    fprintf(out, "%f\t%f\t%f\t%f\t%f\n", 1.0, 0.0, y[1] * cos(y[2]), y[1] * sin(y[2]), 0.0); // escribimos el instante inicial en el fichero de salida

    // Aquí comienza el algoritmo
    for (int k = 1; k < ITER; k++) {
        // Bucle para cada ecuación diferencial
        for (int n = 1; n <= 4; n++) k1[n] = H * f(n, y[1], y[2], y[3], y[4], k * H);
        for (int n = 1; n <= 4; n++) k2[n] = H * f(n, y[1] + k1[1] / 2.0, y[2] + k1[2] / 2.0, y[3] + k1[3] / 2.0, y[4] + k1[4] / 2.0, k * H + H / 2.0);
        for (int n = 1; n <= 4; n++) k3[n] = H * f(n, y[1] + k2[1] / 2.0, y[2] + k2[2] / 2.0, y[3] + k2[3] / 2.0, y[4] + k2[4] / 2.0, k * H + H / 2.0);
        for (int n = 1; n <= 4; n++) k4[n] = H * f(n, y[1] + k3[1], y[2] + k3[2], y[3] + k3[3], y[4] + k3[4], k * H + H);

        // Se actualiza el vector y[n] para t = t + h
        for (int n = 1; n <= 4; n++) y[n] += 1.0 / 6.0 * (k1[n] + 2 * k2[n] + 2 * k3[n] + k4[n]);

        // Escribimos las nuevas posiciones de la Luna, el cohete y también el tiempo
        fprintf(out, "%f\t%f\t%f\t%f\t%f\n", cos(OMEGA * k * H), sin(OMEGA * k * H), y[1] * cos(y[2]), y[1] * sin(y[2]), (double)(k * H));

        // Se hace una comprobación del ángulo del cohete y la Luna
        if (y[1] * cos(y[2]) >= cos(OMEGA * k * H) && !paso) {
            printf("Ángulo luna: %f\n", atan(sin(OMEGA * k * H) / cos(OMEGA * k * H)));
            printf("Ángulo cohete: %f\n", y[2]);
            printf("Distancia mínima cohete-luna: %f km\n", sqrt(1 + y[1] * y[1] - 2 * y[1] * cos(y[2] - OMEGA * H * k)) * DTL / 1000);
            paso = 1;
        }

       // Calculamos el Hamiltoniano
       double rL = calcular_rL(y[1], y[2], k * H + H);
       double Ham = calcular_Ham(y[1], y[3], y[4], rL);

       // Verificamos que H' sea constante
       Ham_prima = Ham - OMEGA * y[4]; // calculamos H'
       fprintf(cons, "%d\t%.10f\n", k, Ham_prima); // escribimos en el archivo el número de iteración y H'
    }

    fclose(out);
    fclose(cons);

    // Marcamos el tiempo de finalización
    end_time = clock();

     // Calculamos el tiempo total de compilación
    total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Imprimimos el tiempo total de compilación
    fprintf(time, "%.4f\n", total_time);

    fclose(time);

    return 0;
}
