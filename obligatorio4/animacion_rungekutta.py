import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np

# Leer el archivo de salida
with open("rungekutta.dat", "r") as file:  #fijo/adaptado
    data = file.readlines()

# Reducir la cantidad de datos seleccionando solo una muestra
data_sampled = [data[i] for i in range(0, len(data), 500)]  # Cambia el 500 por el factor de muestreo deseado

# Extraer los datos
time = []
x_rocket = []
y_rocket = []
x_moon = []
y_moon = []

for line in data_sampled:
    values = line.split("\t")
    time.append(float(values[4]))
    x_moon.append(float(values[0]))
    y_moon.append(float(values[1]))
    x_rocket.append(float(values[2]))
    y_rocket.append(float(values[3]))

# Duración máxima del video (en segundos)
duracion_maxima = 900  # 15 minutos

# Número máximo de frames basado en la duración máxima y el número de fps
max_frames = duracion_maxima * 10  # suponiendo 10 fps

# Limitar la cantidad de frames según la duración máxima
time = time[:max_frames]
x_moon = x_moon[:max_frames]
y_moon = y_moon[:max_frames]
x_rocket = x_rocket[:max_frames]
y_rocket = y_rocket[:max_frames]

# Preparar la figura
fig, ax = plt.subplots(figsize=(10, 6))

# Inicializar el contador de frames
frame_counter = 0

# Función para actualizar el gráfico en cada frame de la animación
def update(frame):
    ax.clear()

    # Dibujar la Tierra (más grande)
    ax.scatter(0, 0, color='blue', label='Tierra', s=200)

    # Dibujar la Luna como un círculo
    ax.scatter(x_moon[frame], y_moon[frame], color='gray', label='Luna', s=40)

    # Trayectoria del cohete
    ax.plot(x_rocket[:frame+1], y_rocket[:frame+1], label="Trayectoria del Cohete", color='orange', linewidth=1)

    # Dibujar el cohete en su posición actual (más pequeño)
    rocket_img = plt.imread('rocket.png')
    rocket_position = (x_rocket[frame], y_rocket[frame])
    rocket_box = OffsetImage(rocket_img, zoom=0.015)
    rocket_annotation = AnnotationBbox(rocket_box, rocket_position, frameon=False)
    ax.add_artist(rocket_annotation)

    # Dibujar la trayectoria de la Luna con líneas discontinuas
    ax.plot(x_moon[:frame+1], y_moon[:frame+1], color='gray', linestyle='dashed')

    ax.set_title("Trayectoria del Cohete hacia la Luna")
    ax.set_xlabel("Posición X (m)")
    ax.set_ylabel("Posición Y (m)")
    ax.legend()
    ax.grid(True)
    
    # Incrementar el contador de frames
    global frame_counter
    frame_counter += 1

# Crear la animación
ani = FuncAnimation(fig, update, frames=len(time), interval=100)

# Guardar la animación como un archivo de video
ani.save('trayectoria_cohete.mp4', fps=10)

# Agregar la leyenda con los parámetros
plt.text(0.02, 0.02, "1e6 iteraciones\n"
                      "theta = PI / 3.387\n"
                      "v = 11200.0 / DTL",
         transform=ax.transAxes, fontsize=8, color='black')

plt.show()
