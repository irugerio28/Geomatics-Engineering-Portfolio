"""
Autor: Ivan Rugerio
Título: Cálculo de Posicionamiento GNSS por Mínimos Cuadrados
Descripción: Algoritmo para determinar la posición del receptor (x, y, z) y el error del reloj 
             utilizando pseudodistancias y coordenadas satelitales (ECEF).
Tecnologías: Python, NumPy, Matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt

# --- Datos de entrada (Simulación) ---
# Coordenadas de los satélites en el sistema ECEF (metros)
# Formato: [X, Y, Z]
sat = np.array([
    [15600000,   7540000,  20140000],
    [18760000,   2750000,  18610000], 
    [17610000,  14640000,  13480000],
    [19170000,   6100000,  18390000]
])

# Pseudodistancias observadas (rho) en metros
rho = np.array([
    [21026883],
    [21549505],
    [20235671],
    [21741496]
])

# --- Configuración de la Matriz de Pesos ---
# Se asume igual peso para todas las observaciones (0.25)
v = np.array([0.25, 0.25, 0.25, 0.25])
Cl = np.diag(v)
ClInv = np.linalg.inv(Cl) 

# --- Inicialización del Receptor ---
# Punto de partida en el centro de la tierra (0,0,0) para la iteración
xr = 0.0
yr = 0.0
zr = 0.0
ErrorReloji = 0.0 

print('--- INICIO DEL PROCESO ITERATIVO ---')
print(f'Posición inicial: ({xr}, {yr}, {zr}) m')

# --- Parámetros de iteración ---
Tolerancia = 1e-3  # Convergencia requerida (1 mm)
maxIter = 100
dxtotal = 1.0 
c = 299792458      # Velocidad de la luz (m/s)
iter = 0
hist_dx = []       # Historial para graficar convergencia

# --- Bucle de Ajuste por Mínimos Cuadrados ---
while dxtotal > Tolerancia and iter < maxIter:
    iter += 1

    # Distancia geométrica actual (r)
    r_vec = np.sqrt((sat[:, 0] - xr)**2 + (sat[:, 1] - yr)**2 + (sat[:, 2] - zr)**2)
    r = r_vec.reshape(-1, 1) 

    # Construcción de la Matriz de Diseño (A)
    # Derivadas parciales respecto a x, y, z y el error de reloj
    A_col1 = -(sat[:, 0] - xr) / r_vec
    A_col2 = -(sat[:, 1] - yr) / r_vec
    A_col3 = -(sat[:, 2] - zr) / r_vec
    A_col4 = np.ones(4)
    A = np.column_stack((A_col1, A_col2, A_col3, A_col4))

    # Vector de residuos (Observado - Calculado)
    W = rho - (r + c * ErrorReloji)

    # Resolución del sistema normal: (A^T * P * A) * dx = A^T * P * W
    try:
        lhs = A.T @ ClInv @ A
        rhs = A.T @ ClInv @ W
        dx = np.linalg.solve(lhs, rhs)
    except np.linalg.LinAlgError:
        print(f"Error: Matriz singular en la iteración {iter}. Deteniendo.")
        break

    # Actualización de coordenadas
    xr = xr + dx[0, 0]
    yr = yr + dx[1, 0]
    zr = zr + dx[2, 0]
    ErrorReloji = ErrorReloji + dx[3, 0] / c

    dxtotal = np.linalg.norm(dx)
    hist_dx.append(dxtotal)

    print(f'Iteración {iter} | Ajuste aplicado: {dxtotal:.6f} m')

# --- Evaluación de Calidad (RMS) ---
r_final_vec = np.sqrt((sat[:, 0] - xr)**2 + (sat[:, 1] - yr)**2 + (sat[:, 2] - zr)**2)
r_final = r_final_vec.reshape(-1, 1)
V = rho - (r_final + c * ErrorReloji)
RMS = np.sqrt(np.mean(V**2))

# --- Reporte de Resultados ---
print('\n--- RESULTADOS FINALES ---')
if iter == maxIter:
    print(f'ADVERTENCIA: No se alcanzó la convergencia.')
else:
    print(f'Convergencia alcanzada en {iter} iteraciones.')

print(f'Posición X: {xr:.3f} m')
print(f'Posición Y: {yr:.3f} m')
print(f'Posición Z: {zr:.3f} m')
print(f'Error de Reloj: {ErrorReloji:.9f} s')
print(f'RMS del ajuste: {RMS:.6f} m')

# --- Visualización ---
plt.figure(figsize=(10, 6))
plt.plot(hist_dx, 'o-', color='navy', label='Magnitud del ajuste')
plt.yscale('log')
plt.xlabel('Número de Iteración')
plt.ylabel('Desplazamiento de corrección ||dx|| (m)')
plt.title('Convergencia del Ajuste de Posición GNSS')
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend()
plt.show()