import numpy as np
import matplotlib.pyplot as plt

# Cargar los datos del archivo
data = np.loadtxt('hamiltoniano.txt')

# Separar los datos en dos arrays: paso y hamiltoniano
steps = data[:, 0]
hamiltoniano = data[:, 1]

# Recoger 1 de cada 1000 datos
steps_reduced = steps[::1000]
hamiltoniano_reduced = hamiltoniano[::1000]

# Crear la gráfica
plt.figure(figsize=(10, 6))
plt.plot(steps_reduced, hamiltoniano_reduced, label='Hamiltoniano')
plt.xlabel('Iteraciones')
plt.ylabel('Valor del Hamiltoniano')
plt.title('Conservación de la Energía del Sistema de Péndulo Doble')
plt.legend()
plt.grid(True)

# Guardar la gráfica en un archivo
plt.savefig('hamiltoniano.png')

# Mostrar la gráfica
plt.show()

