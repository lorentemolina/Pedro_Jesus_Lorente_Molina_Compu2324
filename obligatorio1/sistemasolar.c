#include <stdio.h>
#include <math.h>

#define G 6.67430e-11    // Constante de gravitación universal en m^3/kg/s^2
#define Ms 1.989e30      // Masa del Sol en kg
#define distST 1.496e11  // Unidad de longitud: distancia entre la Tierra y el Sol en metros
#define Ft 5022004.955   // Factor de conversión tiempo

// Estructura para almacenar los datos de cada planeta
typedef struct {
    double mass;     // Masa del planeta en kg
    double distance; // Distancia del planeta al Sol en unidades de longitud
    double velocity; // Velocidad orbital del planeta en unidades de longitud/tiempo
} Planet;

// Función para calcular la aceleración gravitatoria
void acceleration(double positions[][2], double acc[][2], int n, Planet *planet_data) {
    int i, j;
    double rij[2], rij_mag;
    for (i = 0; i < n; i++) {
        acc[i][0] = 0.0;
        acc[i][1] = 0.0;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                rij[0] = positions[i][0] - positions[j][0];
                rij[1] = positions[i][1] - positions[j][1];
                rij_mag = sqrt(rij[0] * rij[0] + rij[1] * rij[1]);
                acc[i][0] += -planet_data[j].mass * rij[0] / pow(rij_mag, 3);
                acc[i][1] += -planet_data[j].mass * rij[1] / pow(rij_mag, 3);
            }
        }
    }
}

// Función para inicializar los valores de posiciones y velocidades
void initialize_positions_and_velocities(double positions[][2], double previous_positions[][2], double velocitiesV[][2], Planet *planet_data) {
    for (int i = 0; i < 5; i++) {
        positions[i][0] = planet_data[i].distance;
        positions[i][1] = 0.0;
        previous_positions[i][0] = positions[i][0];
        previous_positions[i][1] = positions[i][1];
        velocitiesV[i][0] = 0.0;
        velocitiesV[i][1] = planet_data[i].velocity;
    }
}

// Función que calcula el potencial gravitatorio y lo almacena en V
double calculate_gravitational_potential(double positions[][2], int n, Planet *planet_data) {
    double V = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double rij_x = positions[i][0] - positions[j][0];
                double rij_y = positions[i][1] - positions[j][1];
                double rij_mag = sqrt(rij_x * rij_x + rij_y * rij_y);
                V += -planet_data[i].mass * planet_data[j].mass / (2 * rij_mag);
            }
        }
    }
    return V;
}

// Función que calcula la energía cinética y la almacena en T
double calculate_kinetic_energy(double velocitiesV[][2], int n, Planet *planet_data) {
    double T = 0.0;
    for (int i = 0; i < n; i++) {
        double v_mag = sqrt(velocitiesV[i][0] * velocitiesV[i][0] + velocitiesV[i][1] * velocitiesV[i][1]);
        T += 0.5 * planet_data[i].mass * v_mag * v_mag;
    }
    return T;
}

// Función para calcular el momento angular de un planeta
double calculate_angular_momentum(double positions[][2], double velocitiesV[][2], int n, Planet *planet_data) {
    double L = 0.0;
    for (int i = 0; i < n; i++) {
        double x = positions[i][0];
        double y = positions[i][1];
        double vx = velocitiesV[i][0];
        double vy = velocitiesV[i][1];
        L += planet_data[i].mass * (x * vy - y * vx);
    }
    return L;
}

int main() {
     // Datos de los planetas (reescalados)
    Planet planet_data[] = {
        {1.0, 0.0, 0.0}, // Sol
        {3.301e23 / Ms, 5.791e10 / distST, 4.736e4 / distST * Ft}, // Mercurio
        {4.867e24 / Ms, 1.082e11 / distST, 3.502e4 / distST * Ft}, // Venus
        {5.972e24 / Ms, 1.496e11 / distST, 2.978e4 / distST * Ft}, // Tierra
        {6.39e23 / Ms, 2.279e11 / distST, 2.413e4 / distST * Ft} // Marte
    };

    // Parámetros de simulación
    int timesteps = 20000; // Pasos calculados en la simulación
    double h = 1.0 / 1000; // Intervalo de tiempo
    double t = 0.0; // Tiempo

    // Inicialización de posiciones y velocidades
    double positions[5][2];
    double previous_positions[5][2];
    double velocitiesV[5][2];
    double velocitiesW[5][2];
    initialize_positions_and_velocities(positions, previous_positions, velocitiesV, planet_data);

    // Variables para determinar y comprobar el periodo de los planetas
    int orbit[5] = {0}; // Comprobación de órbita completada (sin contar el Sol)
    double periods[5] = {0}; // Periodo

    // Simulación utilizando el algoritmo de Verlet en velocidad
    double acc[5][2];
    double V;      // Potencial gravitatorio
    double T;      // Energía cinética
    double E;      // Energía total
    double L;      // Momento angular

    // Abrimos los archivos de escritura
    FILE *f = fopen("planets_data.dat", "w");
    FILE *f2 = fopen("energiaplanetas.dat", "w");
    FILE *f3 = fopen("momangplanetas.dat", "w");
    FILE *f4 = fopen("geocent.dat", "w");
    FILE *f5 = fopen("periodos.dat", "w");

    // Calculamos la aceleración inicial (t=0.0)
    acceleration(positions, acc, 5, planet_data);

    for (int k = 0; k < timesteps; k++) {

         // Paso 0: Cálculo de la energía cinética y potencial para obtener la energía total en cada iteración
        T = calculate_kinetic_energy(velocitiesV, 5, planet_data);   // Inicializo a 0.0 dentro de la función
        V = calculate_gravitational_potential(positions, 5, planet_data);  // Inicializo a 0.0 dentro de la función
        E = T + V;

        // y cálculo del momento angular total en cada iteración
        L = calculate_angular_momentum(positions, velocitiesV, 5, planet_data);  // Inicializo a 0.0 dentro de la función

        // Periodo de los planetas
        for (int i = 0; i < 5; i++) {
            // Verifica si se cumple la condición para determinar una órbita completa del planeta
            if (positions[i][0] > 0.0 && positions[i][1] * previous_positions[i][1] < 0.0 && i != 0 && k > 100 && orbit[i] == 0) {
                orbit[i] = 1; // Marca la vuelta como completada
                periods[i] = k * h * Ft * 1 / 3600 * 1 / 24; // Calcula el periodo del planeta en días

                // Escribe el periodo calculado en el archivo
                fprintf(f5, "Periodo planeta %d %.2lf días.\n", i, periods[i]);
            }
        }

        // Paso 1: Calcular las nuevas posiciones y velocidadesw
        for (int j = 0; j < 5; j++) {

            previous_positions[j][1] = positions[j][1]; // Para usarlo en la comprobación del periodo

            positions[j][0] += velocitiesV[j][0] * h + 0.5 * acc[j][0] * h * h;
            positions[j][1] += velocitiesV[j][1] * h + 0.5 * acc[j][1] * h * h;

            velocitiesW[j][0] = velocitiesV[j][0] + 0.5 * acc[j][0] * h;
            velocitiesW[j][1] = velocitiesV[j][1] + 0.5 * acc[j][1] * h;
        }

        // Paso 2: Calcular las nuevas aceleraciones
        acceleration(positions, acc, 5, planet_data);

        // Paso 3: Calcular las nuevas velocidadesv
        for (int j = 0; j < 5; j++) {
            velocitiesV[j][0] = velocitiesW[j][0] + 0.5 * acc[j][0] * h;
            velocitiesV[j][1] = velocitiesW[j][1] + 0.5 * acc[j][1] * h;
        }

        // Paso 4: Escribir las posiciones de los planetas en el archivo
        if (k > 0) {
            fprintf(f, "\n");
            fprintf(f4, "\n");
        }
        for (int j = 0; j < 5; j++) {
            fprintf(f, "%.10f,%.10f\n", positions[j][0], positions[j][1]);
            fprintf(f4, "%.10f,%.10f\n", positions[j][0] - positions[3][0], positions[j][1] - positions[3][1]);
        }

        // Escribir la energía total para comprobar su conservación
        fprintf(f2, "%d %lf\n", k, E);

        // Escribir el momento angular total para comprobar su conservación
        fprintf(f3, "%d %lf\n", k, L);

        // Paso 5: Aumentar el tiempo
        t += h;
    }

    // Cerramos los archivos de escritura
    fclose(f);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);

    return 0;
}

