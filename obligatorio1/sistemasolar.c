#include <stdio.h>
#include <math.h>

#define G 6.67430e-11    // Constante de gravitación universal en m^3/kg/s^2
#define Ms 1.989e30      // Masa del Sol en kg
#define distST 1.496e11  // Unidad de longitud: distancia entre la Tierra y el Sol en metros
#define Ft 5022004.955 // Factor de conversión tiempo

// Estructura para almacenar los datos de cada planeta
typedef struct {
    double mass;     // Masa del planeta en kg
    double distance; // Distancia del planeta al Sol en unidades de longitud
    double velocity; // Velocidad orbital del planeta en unidades de longitud/tiempo
} Planet;

// Función para calcular la aceleración gravitatoria
void acceleration(double *positions, double *acc, int n, Planet *planet_data) {
    int i, j;
    double rij[2], rij_mag;  //rij es un array de dos dimensiones para almacenar las distancias entre dos planetas en el espacio bidimensional
    for (i = 0; i < n; i++) {  //para inicializar a 0 las aceleraciones antes de volver a calcular
        acc[i*2] = 0.0;
        acc[i*2 + 1] = 0.0;
    }
    for (i = 0; i < n; i++) { 
        for (j = 0; j < n; j++) {
            if (i != j) {
                rij[0] = positions[i*2] - positions[j*2]; //componente x
                rij[1] = positions[i*2 + 1] - positions[j*2 + 1]; //componente y
                rij_mag = sqrt(rij[0]*rij[0] + rij[1]*rij[1]);  //módulo
                acc[i*2] += - planet_data[j].mass * rij[0] / pow(rij_mag, 3);
                acc[i*2+1] += - planet_data[j].mass * rij[1] / pow(rij_mag, 3);
            }
        }
    }
}

// Función para inicializar los valores de posiciones y velocidades
void initialize_positions_and_velocities(double *positions, double *previous_positions, double *velocitiesV, Planet *planet_data) {
    for (int i = 0; i < 5; i++) {
        positions[i*2] = planet_data[i].distance; // Posición inicial en x 
        positions[i*2 + 1] = 0.0; // Posición inicial en y (0 en este caso)
        previous_positions[i*2] = positions[i*2];  // Posición previa en x inicializada para la comprobación de las vueltas
        previous_positions[i*2 + 1] = positions[i*2 + 1]; // Posición previa en y inicializada para la comprobación de las vueltas
        velocitiesV[i*2] = 0.0; // Velocidad inicial en x (0 en este caso)
        velocitiesV[i*2 + 1] = planet_data[i].velocity; // Velocidad inicial en y
    }
}

// Función que calcula el potencial gravitatorio y lo almacena en V
double calculate_gravitational_potential(double *positions, int n, Planet *planet_data) {
    double V = 0.0;

     for (int i = 0; i < n; i++) { 
        for (int j = 0; j < n; j++) {
            if (i != j) {
                double rij_x = positions[i*2] - positions[j*2];
                double rij_y = positions[i*2 + 1] - positions[j*2 + 1];
                double rij_mag = sqrt(rij_x * rij_x + rij_y * rij_y);
                
                V += - planet_data[i].mass * planet_data[j].mass / (2*rij_mag);
            }
        }
     }
    return V;
}

// Función que calcula la energía cinética y la almacena en T
double calculate_kinetic_energy(double *velocitiesV, int n, Planet *planet_data) {
    double T = 0.0;
    for (int i = 0; i < n; i++) {
        double v_mag = sqrt(velocitiesV[i*2] * velocitiesV[i*2] + velocitiesV[i*2+1] * velocitiesV[i*2+1]);
        T += 0.5 * planet_data[i].mass * v_mag * v_mag;
    }
    return T;
}

// Función para calcular el momento angular de un planeta
double calculate_angular_momentum(double *positions, double *velocitiesV, int n, Planet *planet_data) {
    double L = 0.0; 
    for (int i = 0; i < n; i++) {  //momento angular del Sol considerado nulo (no lo tomamos)
        double x = positions[i*2];
        double y = positions[i*2 + 1];
        double vx = velocitiesV[i*2];
        double vy = velocitiesV[i*2 + 1];
        L += planet_data[i].mass * (x * vy - y * vx);
    }
    return L;
}

int main() {
     // Datos de los planetas (reescalados)
    Planet planet_data[] = {
        {1.0, 0.0, 0.0}, //Sol
        {3.301e23 / Ms, 5.791e10 / distST, 4.736e4 / distST * Ft}, //Mercurio
        {4.867e24 / Ms, 1.082e11 / distST, 3.502e4 / distST * Ft}, //Venus
        {5.972e24 / Ms, 1.496e11 / distST, 2.978e4 / distST * Ft}, //Tierra
        {6.39e23 / Ms, 2.279e11 / distST, 2.413e4 / distST * Ft} //Marte
    };

    //{1.898e27 / Ms, 7.786e11 / distST, 1.307e4 / distST * Ft}, // Júpiter
    //{5.683e26 / Ms, 1.434e12 / distST, 9.932e3 / distST * Ft}, // Saturno
    //{8.681e25 / Ms, 2.871e12 / distST, 6.649e3 / distST * Ft}, // Urano
    //{1.024e26 / Ms, 4.495e12 / distST, 5.447e3 / distST * Ft} // Neptuno

    // Parámetros de simulación
    int timesteps = 10000; //pasos calculados en la simulación
    double h = 1.0 / 1000; //intervalo de tiempo
    double t = 0.0; //tiempo

    // Inicialización de posiciones y velocidadesv
    double positions[10];
    double previous_positions[10];
    double velocitiesV[10];
    double velocitiesW[10];
    initialize_positions_and_velocities(positions, previous_positions, velocitiesV, planet_data);

    // Lista para almacenar las órbitas relativas de los planetas respecto al Sol y la Tierra
    double relative_orbitsH[timesteps][5][2]; 
    double relative_orbitsG[timesteps][5][2]; 

    // Variables para determinar y comprobar el periodo de los planetas
    int orbit[5] = {0}; //comprobación de órbita completada (sin contar el Sol)
    double periods[5] = {0}; //periodo

    // Simulación utilizando el algoritmo de Verlet en velocidad
    double acc[10];
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

         // Paso 0: Cálculo de la energía cinética y potencial para obtener la energía total en cada iteracción
        T = calculate_kinetic_energy(velocitiesV, 5, planet_data);   //inicializo a 0.0 dentro de la función
        V = calculate_gravitational_potential(positions, 5, planet_data);  //inicializo a 0.0 dentro de la función
        E = T + V;
        
        // y cálculo del momento angular total en cada iteracción
        L = calculate_angular_momentum(positions, velocitiesV, 5, planet_data);  //inicializo a 0.0 dentro de la función

        // Periodo de los planetas
        for (int i = 0; i < 5; i++) {
            // Verifica si se cumple la condición para determinar una órbita completa del planeta
            if (positions[i*2] > 0.0 && positions[i*2 + 1] * previous_positions[i*2 + 1] < 0.0 && i != 0 && k > 100 && orbit[i] == 0) {
                orbit[i] = 1; // Marca la vuelta como completada
                periods[i] = k * h * Ft * 1/3600 * 1/24; // Calcula el periodo del planeta en días

            // Escribe el periodo calculado en el archivo
            fprintf(f5, "Periodo planeta %d %.2lf días.\n", i, periods[i]);
            }
        }

       
        // Paso 1: Calcular las nuevas posiciones y velocidadesw
        for (int j = 0; j < 5; j++) {

            previous_positions[j*2 + 1] = positions[j*2 + 1]; //para usarlo en la comprobación del periodo

            positions[j*2] += velocitiesV[j*2] * h + 0.5 * acc[j*2] * h*h;
            positions[j*2 + 1] += velocitiesV[j*2 + 1] * h + 0.5 * acc[j*2 + 1] * h*h;

            velocitiesW[j*2] = velocitiesV[j*2] + 0.5 * acc[j*2] * h;
            velocitiesW[j*2 + 1] = velocitiesV[j*2 + 1] + 0.5 * acc[j*2 + 1] * h;
        }

        // Paso 2: Calcular las nuevas aceleraciones
        acceleration(positions, acc, 5, planet_data);

        // Paso 3: Calcular las nuevas velocidadesv
        for (int j = 0; j < 5; j++) {
            velocitiesV[j*2] = velocitiesW[j*2] + 0.5 * acc[j*2] * h;
            velocitiesV[j*2 + 1] = velocitiesW[j*2 +1] + 0.5 * acc[j*2 + 1] * h;
        }

        // Paso 4: Calcular la órbita relativa de los planetas respecto al Sol (modelo heliocéntrico)
        for (int j = 0; j < 5; j++) {
            relative_orbitsH[k][j][0] = positions[j*2];      
            relative_orbitsH[k][j][1] = positions[j*2 + 1];   
        }

        // y respecto a la Tierra (modelo geocéntrico)
        for (int j = 0; j < 5; j++) {
            relative_orbitsG[k][j][0] = (positions[j*2] - positions[3*2]);  
            relative_orbitsG[k][j][1] = (positions[j*2 + 1] - positions[3*2 +1]); 
        }

        // Paso 5: Escribir las órbitas relativas en los archivos
        if (k > 0) {
            fprintf(f, "\n");
            fprintf(f4, "\n");
        }
        for (int j = 0; j < 5; j++) {
            fprintf(f, "%.10f,%.10f\n", relative_orbitsH[k][j][0], relative_orbitsH[k][j][1]);
            fprintf(f4, "%.10f,%.10f\n", relative_orbitsG[k][j][0], relative_orbitsG[k][j][1]);
        }
        // escribimos la energía total para comprobar su conservación
        fprintf(f2, "%d %lf\n", k, E);

        // escribimos el momento angular total para comprobar su conservación
        fprintf(f3, "%d %lf\n", k, L);

        // Paso 6: Aumentar el tiempo
            t += h;
    }

    // cerramos los archivos de escritura
    fclose(f);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);

    return 0;
}

