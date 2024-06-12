#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// Definimos las constantes
#define PI 3.14159265358979323846
#define G 9.81
#define L1 1.0
#define L2 1.0
#define M1 1.0
#define M2 1.0

// Definimos los parámetros de la simulación
#define ITER 1e7 // iteraciones 
#define H 0.0001 // paso de tiempo 

// Estructura para almacenar el estado del sistema
typedef struct {
    double phi;
    double phi_dot;
    double psi;
    double psi_dot;
} State;

// Función que calcula las derivadas del sistema
void derivs(const State *s, State *ds) {
    double delta = s->psi - s->phi;
    double denom1 = (M1 + M2) * L1 - M2 * L1 * cos(delta) * cos(delta);
    double denom2 = (L2 / L1) * denom1;

    ds->phi = s->phi_dot;
    ds->psi = s->psi_dot;

    ds->phi_dot = (M2 * L1 * s->phi_dot * s->phi_dot * sin(delta) * cos(delta)
                  + M2 * G * sin(s->psi) * cos(delta)
                  + M2 * L2 * s->psi_dot * s->psi_dot * sin(delta)
                  - (M1 + M2) * G * sin(s->phi)) / denom1;

    ds->psi_dot = (-M2 * L2 * s->psi_dot * s->psi_dot * sin(delta) * cos(delta)
                  + (M1 + M2) * (G * sin(s->phi) * cos(delta)
                  - L1 * s->phi_dot * s->phi_dot * sin(delta)
                  - G * sin(s->psi))) / denom2;
}

// Función para ejecutar el método de Runge-Kutta de cuarto orden
void runge_kutta(State *s, double h) {
    State k1, k2, k3, k4, temp;
    derivs(s, &k1);

    temp.phi = s->phi + 0.5 * h * k1.phi;
    temp.phi_dot = s->phi_dot + 0.5 * h * k1.phi_dot;
    temp.psi = s->psi + 0.5 * h * k1.psi;
    temp.psi_dot = s->psi_dot + 0.5 * h * k1.psi_dot;
    derivs(&temp, &k2);

    temp.phi = s->phi + 0.5 * h * k2.phi;
    temp.phi_dot = s->phi_dot + 0.5 * h * k2.phi_dot;
    temp.psi = s->psi + 0.5 * h * k2.psi;
    temp.psi_dot = s->psi_dot + 0.5 * h * k2.psi_dot;
    derivs(&temp, &k3);

    temp.phi = s->phi + h * k3.phi;
    temp.phi_dot = s->phi_dot + h * k3.phi_dot;
    temp.psi = s->psi + h * k3.psi;
    temp.psi_dot = s->psi_dot + h * k3.psi_dot;
    derivs(&temp, &k4);

    s->phi += (h / 6.0) * (k1.phi + 2.0 * k2.phi + 2.0 * k3.phi + k4.phi);
    s->phi_dot += (h / 6.0) * (k1.phi_dot + 2.0 * k2.phi_dot + 2.0 * k3.phi_dot + k4.phi_dot);
    s->psi += (h / 6.0) * (k1.psi + 2.0 * k2.psi + 2.0 * k3.psi + k4.psi);
    s->psi_dot += (h / 6.0) * (k1.psi_dot + 2.0 * k2.psi_dot + 2.0 * k3.psi_dot + k4.psi_dot);
}

// Función para calcular la energía del sistema
double calcular_hamiltoniano(const State *s) {
    double delta = s->psi - s->phi;
    double T = 0.5 * (M1 + M2) * L1 * L1 * s->phi_dot * s->phi_dot
             + 0.5 * M2 * L2 * L2 * s->psi_dot * s->psi_dot
             + M2 * L1 * L2 * s->phi_dot * s->psi_dot * cos(delta);
    double V = (M1 + M2) * G * L1 * (1 - cos(s->phi)) + M2 * G * L2 * (1 - cos(s->psi));
    return T + V;
}

// Función para calcular la distancia entre dos estados
double calcular_distancia(const State *s1, const State *s2) {
    return sqrt((s1->phi - s2->phi) * (s1->phi - s2->phi) +
                (s1->phi_dot - s2->phi_dot) * (s1->phi_dot - s2->phi_dot) +
                (s1->psi - s2->psi) * (s1->psi - s2->psi) +
                (s1->psi_dot - s2->psi_dot) * (s1->psi_dot - s2->psi_dot));
}

// Función para calcular el coeficiente de Lyapunov
double calcular_coeficiente_lyapunov(State *s1, State *s2, double h, int iteraciones) {
    double d0 = calcular_distancia(s1, s2);
    double d = d0;
    double lyapunov_sum = 0.0;

    for (int i = 0; i < iteraciones; i++) {
        runge_kutta(s1, h);
        runge_kutta(s2, h);

        double d_new = calcular_distancia(s1, s2);
        lyapunov_sum += log(d_new / d);
        d = d_new;

        // Reescalar la segunda trayectoria para evitar el desbordamiento numérico
        // s2->phi_dot = s1->phi_dot + (s2->phi_dot - s1->phi_dot) * (d0 / d);
        // s2->psi = s1->psi + (s2->psi - s1->psi) * (d0 / d);
        // s2->psi_dot = s1->psi_dot + (s2->psi_dot - s1->psi_dot) * (d0 / d);
    }

    return lyapunov_sum / (iteraciones * h);
}

int main(int argc, char *argv[]) {

    // Para observar los tiempos de compilación
    clock_t start_time, end_time;
    double total_time;

    // Marcamos el tiempo de inicio
    start_time = clock();

    if (argc < 2) {
        printf("Usage: %s <energy>\n", argv[0]);
        return 1;
    }

    // Energía objetivo
    double E = atof(argv[1]);

    // Inicializamos el estado del sistema con algunas condiciones iniciales
    State s1 = {PI/18, 0.0, PI/18, 0.0};

    // Ajustamos la velocidad angular para alcanzar la energía deseada
    s1.phi_dot = sqrt((2 * (E - ((M1 + M2) * G * L1 * (1 - cos(s1.phi)) + M2 * G * L2 * (1 - cos(s1.psi))))) / ((M1 + M2) * L1 * L1));


    // Creamos una copia perturbada del estado inicial
    State s2 = s1;
    s2.phi += 1e-3;  // Perturbación inicial
    s2.phi_dot += 1e-3;
    s2.psi += 1e-3;
    s2.psi_dot += 1e-3;

    // Calculamos el coeficiente de Lyapunov
    double coeficiente_lyapunov = calcular_coeficiente_lyapunov(&s1, &s2, H, ITER);
    printf("Coeficiente de Lyapunov: %.6f\n", coeficiente_lyapunov);

    // Archivo para guardar los datos de Lyapunov
    FILE *lyapunov_file = fopen("lyapunov_data.txt", "a");
    if (lyapunov_file == NULL) {
        printf("No se pudo abrir el archivo para los datos de Lyapunov.\n");
        return 1;
    }

    // Guardamos las condiciones iniciales y el coeficiente de Lyapunov en el archivo
    fprintf(lyapunov_file, "%f\t%f\t%f\t%f\t%f\t%f\n", s1.phi, s1.phi_dot, s1.psi, s1.psi_dot, coeficiente_lyapunov, E);
    fclose(lyapunov_file);

    // Reinicializamos el estado del sistema para la simulación
    s1.phi = PI/18;
    s1.psi = PI/18;
    s1.psi_dot = 0.0;
    s1.phi_dot = sqrt((2 * (E - ((M1 + M2) * G * L1 * (1 - cos(s1.phi)) + M2 * G * L2 * (1 - cos(s1.psi))))) / ((M1 + M2) * L1 * L1));


    // Archivos de salida
    char filename[50];
    sprintf(filename, "poincare_phi_phi_dot_E%.0f.txt", E);
    FILE *out = fopen(filename, "w");
    if (out == NULL) {
        printf("No se pudo abrir el archivo de salida para el mapa de Poincaré.\n");
        return 1;
    }

    char filename2[50];
    sprintf(filename2, "poincare_phi_psi_E%.0f.txt", E);
    FILE *out2 = fopen(filename2, "w");
    if (out2 == NULL) {
        printf("No se pudo abrir el segundo archivo de salida.\n");
        fclose(out);
        return 1;
    }

    char filename3[50];
    sprintf(filename3, "poincare_psi_psi_dot_E%.0f.txt", E);
    FILE *out3 = fopen(filename3, "w");
    if (out3 == NULL) {
        printf("No se pudo abrir el tercer archivo de salida.\n");
        fclose(out);
        fclose(out2);
        return 1;
    }

    FILE *ham_file = fopen("hamiltoniano.txt", "w");
    if (ham_file == NULL) {
        printf("No se pudo abrir el archivo de salida para el Hamiltoniano.\n");
        fclose(out);
        fclose(out2);
        fclose(out3);
        return 1;
    }

    FILE *time = fopen("tiempo_ejecucion.txt", "w"); // para mostrar el tiempo que tarda en ejecutarse

    double psi_prev = s1.psi;
    double phi_prev = s1.phi;
    double phi_dot_prev = s1.phi_dot;

    // Ejecutamos el método de Runge-Kutta
    for (int i = 0; i < ITER; i++) {
        runge_kutta(&s1, H);
        double hamiltoniano = calcular_hamiltoniano(&s1);
        fprintf(ham_file, "%d\t%.10f\n", i, hamiltoniano);

       // if (i % 100 == 0) {  para hacer los vídeos más rápido
        //if (psi_prev < 0 && s1.psi >= 0) {
            fprintf(out, "%f\t%f\n", s1.phi, s1.phi_dot);
        //}
        //if (phi_dot_prev < 0 && s1.phi_dot >= 0) {  // comentar para graficar la trayectoria
            fprintf(out2, "%f\t%f\n", s1.phi, s1.psi);
        //}
        //if (phi_prev < 0 && s1.phi >= 0) {
            fprintf(out3, "%f\t%f\n", s1.psi, s1.psi_dot);
        //}
       // }

        // Actualizar el estado anterior
        psi_prev = s1.psi;
        phi_prev = s1.phi;
        phi_dot_prev = s1.phi_dot;
    }

     // Marcamos el tiempo de finalización
    end_time = clock();

     // Calculamos el tiempo total de compilación
    total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Imprimimos el tiempo total de compilación
    fprintf(time, "%.4f\n", total_time);

    fclose(time);

    // Cerramos los archivos cuando terminamos de escribir
    fclose(out);
    fclose(out2);
    fclose(out3);
    fclose(ham_file);

    return 0;
}
