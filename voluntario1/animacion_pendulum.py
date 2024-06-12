import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# Configuración de la opción de guardado
guardar_video = False  # Cambia a False si solo quieres mostrar la animación

# Cargar datos del archivo
data = np.loadtxt('poincare_phi_psi_E1.txt')
phi = data[::1000, 0]  # Tomar uno de cada mil datos
psi = data[::1000, 1]  # Tomar uno de cada mil datos

# Longitudes de los péndulos
L1 = 1.0
L2 = 1.0

# Definir posiciones del péndulo
x1 = L1 * np.sin(phi)
y1 = -L1 * np.cos(phi)
x2 = x1 + L2 * np.sin(psi)
y2 = y1 - L2 * np.cos(psi)

# Inicializar la figura y ejes
fig, ax = plt.subplots()
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)

# Inicializar líneas y puntos
line, = ax.plot([], [], 'o-', lw=2, color='blue')
line2, = ax.plot([], [], 'o-', lw=2, color='red')
trail1, = ax.plot([], [], '-', color='blue', alpha=0.5)
trail2, = ax.plot([], [], '-', color='red', alpha=0.5)

# Función de inicialización
def init():
    line.set_data([], [])
    line2.set_data([], [])
    trail1.set_data([], [])
    trail2.set_data([], [])
    return line, line2, trail1, trail2

# Función de actualización para la animación
def update(frame):
    line.set_data([0, x1[frame]], [0, y1[frame]])
    line2.set_data([x1[frame], x2[frame]], [y1[frame], y2[frame]])
    trail1.set_data(x1[:frame+1], y1[:frame+1])
    trail2.set_data(x2[:frame+1], y2[:frame+1])
    return line, line2, trail1, trail2

# Limitar a 5 minutos (300 segundos) con 10 fps -> 300 * 10 = 3000 frames
total_frames = min(len(phi), 300 * 10)

# Crear la animación
ani = FuncAnimation(fig, update, frames=range(total_frames), init_func=init, blit=True, interval=1000/10)

if guardar_video:
    # Configuración para guardar la animación
    writer = FFMpegWriter(fps=10, metadata=dict(artist='Me'), bitrate=1800)
    ani.save('doble_pendulo.mp4', writer=writer)
else:
    # Mostrar la animación
    plt.show()
