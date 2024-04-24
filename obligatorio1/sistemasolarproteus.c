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
    double rij[2], rij_mag;  // rij es un array de dos dimensiones para almacenar las distancias entre dos planetas en el espacio bidimensional
    for (i = 0; i < n; i++) {  // para inicializar a 0 las aceleraciones antes de volver a calcular
        acc[i][0] = 0.0;
        acc[i][1] = 0.0;
    }
    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++) {
            if (i != j) {
                rij[0] = positions[i][0] - positions[j][0]; // componente x
                rij[1] = positions[i][1] - positions[j][1]; // componente y
                rij_mag = sqrt(rij[0]*rij[0] + rij[1]*rij[1]);  // módulo
                acc[i][0] += - planet_data[j].mass * rij[0] / pow(rij_mag, 3);
                acc[i][1] += - planet_data[j].mass * rij[1] / pow(rij_mag, 3);
            }
        }
    }
}

// Función para inicializar los valores de posiciones y velocidades
void initialize_variables(double positions[][2], double previous_positions[][2], double velocitiesV[][2], int *orbit, double *periods, int n, Planet *planet_data) {
    for (int i = 0; i < n; i++) {
        positions[i][0] = planet_data[i].distance; // Posición inicial en x 
        positions[i][1] = 0.0; // Posición inicial en y (0 en este caso)

        previous_positions[i][0] = positions[i][0];  // Posición previa en x inicializada para la comprobación de las vueltas
        previous_positions[i][1] = positions[i][1]; // Posición previa en y inicializada para la comprobación de las vueltas

        velocitiesV[i][0] = 0.0; // Velocidad inicial en x (0 en este caso)
        velocitiesV[i][1] = planet_data[i].velocity; // Velocidad inicial en y

        orbit[i] = 0; // Órbita completa
        periods[i] = 0.0; // Periodos
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
                
                V += - planet_data[i].mass * planet_data[j].mass / (2*rij_mag);
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
    for (int i = 0; i < n; i++) {  // momento angular del Sol considerado nulo (no lo tomamos)
        double x = positions[i][0];
        double y = positions[i][1];
        double vx = velocitiesV[i][0];
        double vy = velocitiesV[i][1];
        L += planet_data[i].mass * (x * vy - y * vx);
    }
    return L;
}

int main() {

    // Para simular un número N de planetas
    int num_planets;

    // Solicitar al usuario el número de planetas
    printf("Ingrese el número de planetas que desea simular (incluyendo el Sol): ");
    scanf("%d", &num_planets);

    // Datos de los planetas (reescalados)
    Planet planet_data[] = {
        {1.0, 0.0, 0.0}, // Sol
        {3.301e23 / Ms, 5.791e10 / distST, 4.736e4 / distST * Ft}, // Mercurio
        {4.867e24 / Ms, 1.082e11 / distST, 3.502e4 / distST * Ft}, // Venus
        {5.972e24 / Ms, 1.496e11 / distST, 2.978e4 / distST * Ft}, // Tierra
        {6.39e23 / Ms, 2.279e11 / distST, 2.413e4 / distST * Ft}, // Marte
        {1.898e27 / Ms, 7.786e11 / distST, 1.307e4 / distST * Ft}, // Júpiter
        {5.683e26 / Ms, 1.434e12 / distST, 9.932e3 / distST * Ft}, // Saturno
        {8.681e25 / Ms, 2.871e12 / distST, 6.649e3 / distST * Ft}, // Urano
        {1.024e26 / Ms, 4.495e12 / distST, 5.447e3 / distST * Ft}  // Neptuno
    };

    // Calcular el número de Tierras necesarias
    int tierras_extra = num_planets - 9;
    if (tierras_extra < 0) {
        tierras_extra = 0;
    }

    // Agregar Tierras adicionales después de Neptuno
    for (int i = 0; i < tierras_extra; i++) {
        planet_data[9 + i] = planet_data[3]; // Copiar los datos de la Tierra
    }

    // Parámetros de simulación
    int timesteps = 10000; // Pasos calculados en la simulación
    double h = 1.0 / 100; // Intervalo de tiempo
    double t = 0.0; // Tiempo

    // Inicialización de posiciones y velocidades
    double positions[num_planets][2];
    double previous_positions[num_planets][2];
    double velocitiesV[num_planets][2];
    double periods[num_planets];
    int orbits[num_planets];
    initialize_variables(positions, previous_positions, velocitiesV, orbits, periods, num_planets, planet_data);

    // Simulación utilizando el algoritmo de Verlet en velocidad
    double acc[num_planets][2];
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
    acceleration(positions, acc, num_planets, planet_data);

    for (int k = 0; k < timesteps; k++) {

        // Paso 0: Cálculo de la energía cinética y potencial para obtener la energía total en cada iteración
        T = calculate_kinetic_energy(velocitiesV, num_planets, planet_data);   // Inicializo a 0.0 dentro de la función
        V = calculate_gravitational_potential(positions, num_planets, planet_data);  // Inicializo a 0.0 dentro de la función
        E = T + V;
        
        // y cálculo del momento angular total en cada iteración
        L = calculate_angular_momentum(positions, velocitiesV, num_planets, planet_data);  // Inicializo a 0.0 dentro de la función

        // Periodo de los planetas
        for (int i = 0; i < num_planets; i++) {
            // Verifica si se cumple la condición para determinar una órbita completa del planeta
            if (positions[i][0] > 0.0 && positions[i][1] * previous_positions[i][1] < 0.0 && i != 0 && k > 100 && orbits[i] == 0) {
                orbits[i] = 1; // Marca la vuelta como completada
                periods[i] = k * h * Ft * 1/3600 * 1/24; // Calcula el periodo del planeta en días

                // Escribe el periodo calculado en el archivo
                fprintf(f5, "Periodo planeta %d %.2lf días.\n", i, periods[i]);
            }
        }

        // Paso 1: Calcular las nuevas posiciones y velocidades
        for (int j = 0; j < num_planets; j++) {

            previous_positions[j][1] = positions[j][1]; // Para usarlo en la comprobación del periodo

            positions[j][0] += velocitiesV[j][0] * h + 0.5 * acc[j][0] * h*h;
            positions[j][1] += velocitiesV[j][1] * h + 0.5 * acc[j][1] * h*h;

            velocitiesV[j][0] += 0.5 * acc[j][0] * h;
            velocitiesV[j][1] += 0.5 * acc[j][1] * h;
        }

        // Paso 2: Calcular las nuevas aceleraciones
        acceleration(positions, acc, num_planets, planet_data);

        // Paso 3: Calcular las nuevas velocidades
        for (int j = 0; j < num_planets; j++) {
            velocitiesV[j][0] += 0.5 * acc[j][0] * h;
            velocitiesV[j][1] += 0.5 * acc[j][1] * h;
        }

        // Paso 4: Escribir las posiciones de los planetas en los archivos
        if (k > 0) {
            fprintf(f, "\n");
            fprintf(f4, "\n");
        }
        for (int j = 0; j < num_planets; j++) {
            fprintf(f, "%.10f,%.10f\n", positions[j][0], positions[j][1]);
            fprintf(f4, "%.10f,%.10f\n", positions[j][0] - positions[3][0], positions[j][1] - positions[3][1]); // Posición relativa a la Tierra
        }

        // Escribimos la energía total para comprobar su conservación
        fprintf(f2, "%d %lf\n", k, E);

        // Escribimos el momento angular total para comprobar su conservación
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


