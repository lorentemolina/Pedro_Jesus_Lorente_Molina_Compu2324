import numpy as np
import matplotlib.pyplot as plt

def plot_poincare(data_file, x_label, y_label, title, output_file):
    # Leer los datos del archivo
    data = np.loadtxt(data_file)
    
    # Extraer los valores
    x_values = data[:, 0]
    y_values = data[:, 1]

    # Graficar el mapa de Poincaré
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, y_values, 'ro', markersize=1)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.grid(True)
    plt.savefig(output_file)  # Guardar la imagen como archivo .png
    plt.show()

# Generar los diferentes mapas de Poincaré
plot_poincare('poincare_phi_phi_dot_E15.txt', 'φ (phi)', 'φ̇ (phi_dot)', 'Mapa de Poincaré: φ vs φ̇', 'poincare_phi_phi_dot_E15.png')
plot_poincare('poincare_phi_psi_E15.txt', 'φ (phi)','ψ (psi)','Mapa de Poincaré: φ vs ψ', 'poincare_phi_psi_E15.png')
plot_poincare('poincare_psi_psi_dot_E15.txt', 'ψ (psi)', 'ψ̇ (psi_dot)', 'Mapa de Poincaré: ψ vs ψ̇', 'poincare_psi_psi_dot_E15.png')